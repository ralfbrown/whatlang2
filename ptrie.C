/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: ptrie.C - packed Word-frequency multi-trie			*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-09						*/
/*									*/
/*  (c) Copyright 2011,2012,2015,2019 Ralf Brown/CMU			*/
/*      This program is free software; you can redistribute it and/or   */
/*      modify it under the terms of the GNU General Public License as  */
/*      published by the Free Software Foundation, version 3.           */
/*                                                                      */
/*      This program is distributed in the hope that it will be         */
/*      useful, but WITHOUT ANY WARRANTY; without even the implied      */
/*      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         */
/*      PURPOSE.  See the GNU General Public License for more details.  */
/*                                                                      */
/*      You should have received a copy of the GNU General Public       */
/*      License (file COPYING) along with this program.  If not, see    */
/*      http://www.gnu.org/licenses/                                    */
/*                                                                      */
/************************************************************************/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include "mtrie.h"
#include "ptrie.h"
#include "framepac/file.h"
#include "framepac/utility.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

// since no node will ever point at the root, we can re-use the root
//   index as the null pointer
#define NOCHILD_INDEX 0

#define MULTITRIE_SIGNATURE "MulTrie\0"
#define MULTITRIE_FORMAT_MIN_VERSION 2 // earliest format we can read
#define MULTITRIE_FORMAT_VERSION 3

// reserve some space for future additions to the file format
#define MULTITRIE_PADBYTES_1  58

/************************************************************************/
/*	Types								*/
/************************************************************************/

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

double PackedTrieFreq::s_value_map[PACKED_TRIE_NUM_VALUES] ;
bool PackedTrieFreq::s_value_map_initialized = false ;

//----------------------------------------------------------------------

void write_escaped_key(Fr::CFile& f, const uint8_t* key, unsigned keylen) ;

/************************************************************************/
/*	Helper functions						*/
/************************************************************************/

#ifndef lengthof
#  define lengthof(x) (sizeof(x)/sizeof((x)[0]))
#endif /* lengthof */

/************************************************************************/
/*	Methods for class PackedTrieFreq				*/
/************************************************************************/

PackedTrieFreq::PackedTrieFreq(uint32_t freq, uint32_t langID, bool last,
			       bool is_stop)
{
   uint32_t data = (langID | (last * PACKED_TRIE_LASTENTRY) |
		    (is_stop * PACKED_TRIE_STOPGRAM)) ;
   uint32_t mant ;
   uint32_t expon ;
   quantize(freq,mant,expon) ;
   data |= mant ;
   data |= (expon << PACKED_TRIE_FREQ_EXP_SHIFT) ;
   m_freqinfo.store(data) ;
   return ;
}

//----------------------------------------------------------------------

void PackedTrieFreq::quantize(uint32_t freq, uint32_t &mantissa, uint32_t &exp)
{
   uint32_t e = 0 ;
   if (freq)
      {
      const uint32_t highbits = PACKED_TRIE_FREQ_HIBITS ;
      const uint32_t max_exponent
	 = PACKED_TRIE_FREQ_EXPONENT >> PACKED_TRIE_FREQ_EXP_SHIFT ;
      while ((freq & highbits) == 0 && e < max_exponent)
	 {
	 freq <<= PTRIE_EXPONENT_SCALE ;
	 e++ ;
	 }
      freq &= PACKED_TRIE_FREQ_MANTISSA ;
      if (freq == 0)
	 freq = PACKED_TRIE_MANTISSA_LSB ;
      }
   mantissa = freq ;
   exp = e ;
   return ;
}

//----------------------------------------------------------------------

PackedTrieFreq::~PackedTrieFreq()
{
   m_freqinfo.store(0U) ;
   return ;
}

//----------------------------------------------------------------------

double PackedTrieFreq::probability(uint32_t ID) const
{
   if (ID == languageID())
      return probability() ;
   else if (!isLast())
      return next()->probability(ID) ;
   return 0 ;
}

//----------------------------------------------------------------------

void PackedTrieFreq::isLast(bool last)
{
   uint32_t data = m_freqinfo.load() & ~PACKED_TRIE_LASTENTRY ;
   if (last)
      data |= PACKED_TRIE_LASTENTRY ;
   m_freqinfo.store(data) ;
   return ;
}

//----------------------------------------------------------------------

void PackedTrieFreq::initDataMapping(double (*mapfn)(uint32_t))
{
   for (size_t i = 0 ; i < PACKED_TRIE_NUM_VALUES ; i++)
      {
      uint32_t scaled = scaledScore(i << PACKED_TRIE_VALUE_SHIFT) ;
      double mapped_value ;
      if (mapfn)
	 {
	 scaled |= (i & 1) ; // add stop-gram bit as LSB
	 mapped_value = mapfn(scaled) ;
	 }
      else
	 {
	 mapped_value = (scaled / 100.0 * TRIE_SCALE_FACTOR) ;
	 if ((i & 1) != 0)
	    mapped_value = -mapped_value ;
	 }
      s_value_map[i] = mapped_value ;
      }
   s_value_map_initialized = true ;
   return ;
}

//----------------------------------------------------------------------

bool PackedTrieFreq::writeDataMapping(Fr::CFile& f)
{
   Fr::UInt32 count(PACKED_TRIE_NUM_VALUES) ;
   return f && f.writeValue(count) && f.writeValue(*s_value_map) ;
}

/************************************************************************/
/*	Methods for class PackedTrieNode				*/
/************************************************************************/

PackedTrieNode::PackedTrieNode()
{
   m_frequency_info.store(INVALID_FREQ) ;
   m_firstchild.store(0U) ;
   memset(m_children,'\0',sizeof(m_children)) ;
   return ;
}

//----------------------------------------------------------------------

bool PackedTrieNode::childPresent(unsigned int N) const 
{
   if (N >= (1<<PTRIE_BITS_PER_LEVEL))
      return false ;
   uint32_t children = m_children[N/32].load() ;
   uint32_t mask = (1U << (N % 32)) ;
   return (children & mask) != 0 ;
}

//----------------------------------------------------------------------

uint32_t PackedTrieNode::childIndex(unsigned int N) const
{
   if (N >= (1<<PTRIE_BITS_PER_LEVEL))
      return LangIDPackedMultiTrie::NULL_INDEX ;
   uint32_t children = m_children[N/32].load() ;
   uint32_t mask = (1U << (N % 32)) - 1 ;
   children &= mask ;
   return (firstChild() + m_popcounts[N/32] + Fr::popcount(children)) ;
}

//----------------------------------------------------------------------

uint32_t PackedTrieNode::childIndexIfPresent(uint8_t N) const
{
#if (1<<PTRIE_BITS_PER_LEVEL) < 256
   if (N >= (1<<PTRIE_BITS_PER_LEVEL))
      return LangIDPackedMultiTrie::NULL_INDEX ;
#endif
   uint32_t children = m_children[N/32].load() ;
   uint32_t mask = (1U << (N % 32)) ;
   if ((children & mask) == 0)
      return LangIDPackedMultiTrie::NULL_INDEX ;
   mask-- ;
   children &= mask ;
   return (firstChild() + m_popcounts[N/32] + Fr::popcount(children)) ;
}

//----------------------------------------------------------------------

uint32_t PackedTrieNode::childIndexIfPresent(unsigned int N) const
{
   if (N >= (1<<PTRIE_BITS_PER_LEVEL))
      return LangIDPackedMultiTrie::NULL_INDEX ;
   uint32_t children = m_children[N/32].load() ;
   uint32_t mask = (1U << (N % 32)) ;
   if ((children & mask) == 0)
      return LangIDPackedMultiTrie::NULL_INDEX ;
   mask-- ;
   children &= mask ;
   return (firstChild() + m_popcounts[N/32] + Fr::popcount(children)) ;
}

//----------------------------------------------------------------------

double PackedTrieNode::probability(const PackedTrieFreq *base,
				   uint32_t langID) const
{
   const PackedTrieFreq *freq = base + m_frequency_info.load() ;
   for ( ; ; freq++)
      {
      if (freq->languageID() == langID)
	 return freq->probability() ;
      if (freq->isLast())
	 break ;
      }
   return 0 ;
}

//----------------------------------------------------------------------

void PackedTrieNode::setChild(unsigned N)
{
   if (N < (1<<PTRIE_BITS_PER_LEVEL))
      {
      uint32_t mask = (1U << (N % 32)) ;
      m_children[N/32] |= mask ;
      }
   return ;
}

//----------------------------------------------------------------------

void PackedTrieNode::setPopCounts()
{
   // set up running population counts for faster lookup of children
   unsigned popcount = 0 ;
   for (size_t i = 0 ; i < lengthof(m_popcounts) ; i++)
      {
      m_popcounts[i] = (uint8_t)popcount ;
      uint32_t children = m_children[i].load() ;
      popcount += Fr::popcount(children) ;
      }
   return ;
}

/************************************************************************/
/*	Methods for class PackedMultiTrie				*/
/************************************************************************/

LangIDPackedMultiTrie::LangIDPackedMultiTrie(const LangIDMultiTrie *multrie)
{
   init() ;
   if (multrie)
      {
      m_numfreq = multrie->countFreqRecords() ;
      m_size = multrie->numFullByteNodes() ;
      m_numterminals = multrie->numTerminalNodes() ;
      m_size -= m_numterminals ;
      m_nodes = Fr::New<PackedTrieNode>(m_size) ;
      m_terminals = Fr::New<PackedTrieTerminalNode>(m_numterminals) ;
      m_freq = Fr::New<PackedTrieFreq>(m_numfreq) ;
      if (m_nodes && m_freq)
	 {
	 const auto mroot = multrie->rootNode() ;
	 auto proot = &m_nodes[ROOT_INDEX] ;
	 new (proot) PackedTrieNode ;
	 m_used = 1 ;
	 if (!insertChildren(proot,multrie,mroot,ROOT_INDEX))
	    {
	    m_size = 0 ;
	    m_numfreq = 0 ;
	    m_numterminals = 0 ;
	    }
	 cout << "   converted " << m_used << " full nodes, "
	      << m_termused << " terminals, and "
	      << m_freqused << " frequencies" << endl ;
	 }
      else
	 {
	 Fr::Free(m_nodes) ;
	 Fr::Free(m_freq) ;
	 m_nodes = nullptr ; 
	 m_freq = nullptr ;
	 m_size = 0 ;
	 m_numfreq = 0 ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

LangIDPackedMultiTrie::LangIDPackedMultiTrie(Fr::CFile& f, const char *filename)
{
   init() ;
   if (f && parseHeader(f))
      {
      size_t offset = f.tell() ;
      auto fmap = new Fr::MemMappedROFile(filename) ;
      if (fmap)
	 {
	 // we can memory-map the file, so just point our member variables
	 //   at the mapped data
	 m_fmap = fmap ;
	 m_nodes = (PackedTrieNode*)(**fmap + offset) ;
	 m_freq = (PackedTrieFreq*)(m_nodes + m_size) ;
	 m_terminals = (PackedTrieTerminalNode*)(m_freq + m_numfreq) ;
	 }
      else
	 {
	 // unable to memory-map the file, so read its contents into buffers
	 //   and point our variables at the buffers
	 char *nodebuffer = Fr::New<char>((m_size * sizeof(PackedTrieNode)) + (m_numterminals * sizeof(PackedTrieTerminalNode))) ;
	 m_nodes = (PackedTrieNode*)nodebuffer ;
	 m_terminals = (PackedTrieTerminalNode*)(m_nodes + m_size) ;
	 m_terminals_contiguous = true ;
	 m_freq = Fr::New<PackedTrieFreq>(m_numfreq) ;
	 if (!m_nodes || !m_freq ||
	    f.read(m_nodes,m_size,sizeof(PackedTrieNode)) != m_size ||
	    f.read(m_freq,m_numfreq,sizeof(PackedTrieFreq)) != m_numfreq ||
	    f.read(m_terminals,m_numterminals,sizeof(PackedTrieTerminalNode)) != m_numterminals)
	    {
	    Fr::Free(m_nodes) ;  m_nodes = nullptr ;
	    Fr::Free(m_freq) ;   m_freq = nullptr ;
	    m_terminals = nullptr ;
	    m_size = 0 ; 
	    m_numfreq = 0 ;
	    m_numterminals = 0 ;
	    }
	 }
      }
   return ;
}

//----------------------------------------------------------------------

LangIDPackedMultiTrie::~LangIDPackedMultiTrie()
{
   if (m_fmap)
      {
      delete m_fmap ;
      m_fmap = nullptr ;
      }
   else
      {
      Fr::Free(m_nodes) ;
      if (!m_terminals_contiguous)
	 Fr::Free(m_terminals) ;
      Fr::Free(m_freq) ;
      }
   init() ;				// clear all of the fields
   return ;
}

//----------------------------------------------------------------------

void LangIDPackedMultiTrie::init()
{
   m_fmap = nullptr ;
   m_nodes = nullptr ;
   m_terminals = nullptr ;
   m_freq = 0 ;
   m_size = 0 ;
   m_used = 0 ;
   m_numterminals = 0 ;
   m_termused = 0 ;
   m_numfreq = 0 ;
   m_freqused = 0 ;
   m_maxkeylen = 0 ;
   m_casesensitivity = CS_Full ;
   m_ignorewhitespace = false ;
   m_terminals_contiguous = false ;
   return ;
}

//----------------------------------------------------------------------

uint32_t LangIDPackedMultiTrie::allocateChildNodes(unsigned numchildren)
{
   uint32_t index = m_used ;
   m_used += numchildren ;
   if (m_used > m_size)
      {
      m_used = m_size ;
      return NOCHILD_INDEX ;		// error!  should never happen!
      }
   // initialize each of the new children
   for (size_t i = 0 ; i < numchildren ; i++)
      {
      new (m_nodes + index + i) PackedTrieNode ;
      }
   return index ;
}

//----------------------------------------------------------------------

uint32_t LangIDPackedMultiTrie::allocateTerminalNodes(unsigned numchildren)
{
   uint32_t index = m_termused ;
   m_termused += numchildren ;
   if (m_termused > m_numterminals)
      {
      m_termused = m_numterminals ;
      return NOCHILD_INDEX ;		// error!  should never happen!
      }
   // initialize each of the new children
   for (size_t i = 0 ; i < numchildren ; i++)
      {
      new (m_terminals + index + i) PackedTrieTerminalNode ;
      }
   return (index | PTRIE_TERMINAL_MASK) ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::insertTerminals(PackedTrieNode *parent,
				      const LangIDMultiTrie *mtrie,
				      const MultiTrieNode *mnode,
				      uint32_t mnode_index,
				      unsigned keylen)
{
   if (!parent || !mnode)
      return false ;
   unsigned numchildren = mnode->numExtensions(mtrie) ;
   if (numchildren == 0)
      return true ;
   keylen++ ;
   if (keylen > longestKey())
      m_maxkeylen = keylen ;
   uint32_t firstchild = allocateTerminalNodes(numchildren) ;
   parent->setFirstChild(firstchild) ;
   if (firstchild == NOCHILD_INDEX)
      {
      cerr << "insertTerminals: firstchild==NOCHILD_INDEX"<<endl;
      return false ;
      }
   unsigned index = 0 ;
   for (unsigned i = 0 ; i < (1<<PTRIE_BITS_PER_LEVEL) ; i++)
      {
      uint32_t nodeindex = mnode_index ;
      if (mtrie->extendKey(nodeindex,(uint8_t)i))
	 {
	 // set the appropriate bit in the child array
	 parent->setChild(i) ;
	 // add frequency info to the child node
	 auto pchild = node(firstchild + index) ;
	 index++ ;
	 const auto mchild = mtrie->node(nodeindex) ;
	 auto numfreq = mchild->numFrequencies() ;
	 if (numfreq > 0)
	    {
	    uint32_t freq_index = m_freqused ;
	    m_freqused += numfreq ;
	    pchild->setFrequencies(freq_index) ;
	    auto mfreq = mchild->frequencies() ;
	    while (mfreq && numfreq > 0)
	       {
	       (void)new (m_freq + freq_index) PackedTrieFreq
		  (mfreq->frequency(),mfreq->languageID(),numfreq <= 1) ;
	       freq_index++ ;
	       numfreq-- ;
	       mfreq = mfreq->next() ;
	       }
	    }
	 else
	    pchild->setFrequencies(INVALID_FREQ) ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::insertChildren(PackedTrieNode *parent,
				     const LangIDMultiTrie *mtrie,
				     const MultiTrieNode *mnode,
				     uint32_t mnode_index,
				     unsigned keylen)
{
   if (!parent || !mnode)
      return false ;
   // first pass: fill in all the children
   unsigned numchildren = mnode->numExtensions(mtrie) ;
   if (numchildren == 0)
      return true ;
   keylen++ ;
   if (keylen > longestKey())
      m_maxkeylen = keylen ;
   bool terminal = mnode->allChildrenAreTerminals(mtrie) ;
   uint32_t firstchild = (terminal
			  ? allocateTerminalNodes(numchildren)
			  : allocateChildNodes(numchildren)) ;
   parent->setFirstChild(firstchild) ;
   if (firstchild == NOCHILD_INDEX)
      {
      cerr << "insertChildren: firstchild==NOCHILD_INDEX" << endl ;
      return false ;
      }
   unsigned index = 0 ;
   for (unsigned i = 0 ; i < (1<<PTRIE_BITS_PER_LEVEL) ; i++)
      {
      uint32_t nodeindex = mnode_index ;
      if (mtrie->extendKey(nodeindex,(uint8_t)i))
	 {
	 // set the appropriate bit in the child array
	 parent->setChild(i) ;
	 // add frequency info to the child node
	 auto pchild = node(firstchild + index) ;
	 index++ ;
	 const auto mchild = mtrie->node(nodeindex) ;
	 auto numfreq = mchild->numFrequencies() ;
	 if (numfreq > 0)
	    {
	    uint32_t freq_index = m_freqused ;
	    m_freqused += numfreq ;
	    pchild->setFrequencies(freq_index) ;
	    auto mfreq = mchild->frequencies() ;
	    while (mfreq && numfreq > 0)
	       {
	       bool is_stop = (mfreq->isStopgram() || mfreq->frequency() == 0);
	       (void)new (m_freq + freq_index) PackedTrieFreq
		  (mfreq->frequency(),mfreq->languageID(),numfreq <= 1,is_stop) ;
	       freq_index++ ;
	       numfreq-- ;
	       mfreq = mfreq->next() ;
	       }
	    }
	 if (terminal)
	    {
	    if (!insertTerminals(pchild,mtrie,mchild,nodeindex,keylen))
	       return false ;
	    }
	 else if (!insertChildren(pchild,mtrie,mchild,nodeindex,keylen))
	    return false ;
	 }
      }
   parent->setPopCounts() ;
   return true ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::parseHeader(Fr::CFile& f)
{
   const size_t siglen = sizeof(MULTITRIE_SIGNATURE) ;
   char signature[siglen] ;
   if (f.read(signature,siglen,sizeof(char)) != siglen ||
       memcmp(signature,MULTITRIE_SIGNATURE,siglen) != 0)
      {
      // error: wrong file type
      return false ;
      }
   unsigned char version ;
   if (!f.readValue(&version)
       || version < MULTITRIE_FORMAT_MIN_VERSION
       || version > MULTITRIE_FORMAT_VERSION)
      {
      // error: wrong version of data file
      return false ;
      }
   unsigned char bits ;
   if (f.readValue(&bits) || bits != PTRIE_BITS_PER_LEVEL)
      {
      // error: wrong type of trie
      return false ;
      }
   char ignore_white ;
   char case_sens ;
   char padbuf[MULTITRIE_PADBYTES_1] ;
   Fr::UInt32 val_size, val_keylen, val_numfreq, val_numterm ;
   if (!f.readValue(&val_size) ||
      !f.readValue(&val_keylen) ||
      !f.readValue(&val_numfreq) ||
      !f.readValue(&val_numterm) ||
      !f.readValue(&ignore_white) ||
      !f.readValue(&case_sens) ||
      f.read(padbuf,sizeof(padbuf),1) != sizeof(padbuf))
      {
      // error reading header
      return false ;
      }
   m_maxkeylen = val_keylen.load() ;
   m_size = val_size.load() ;
   m_numterminals = val_numterm.load() ;
   m_numfreq = val_numfreq.load() ;
   return true ;
}

//----------------------------------------------------------------------

PackedTrieNode* LangIDPackedMultiTrie::findNode(const uint8_t *key, unsigned keylength) const
{
   uint32_t cur_index = ROOT_INDEX ;
   while (keylength > 0)
      {
      if (!extendKey(cur_index,*key))
	 return nullptr ;
      key++ ;
      keylength-- ;
      }
   return node(cur_index) ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::extendKey(uint32_t &nodeindex, uint8_t keybyte) const
{
   if ((nodeindex & PTRIE_TERMINAL_MASK) != 0)
      {
      nodeindex = LangIDPackedMultiTrie::NULL_INDEX ;
      return false ;
      }
#if 0 // currently not used, so don't waste time
   if (ignoringWhiteSpace() && keybyte == ' ')
      return true ;
   switch (caseSensitivity())
      {
      case CS_Full:
      default:
	 // do nothing
	 break ;
      case CS_ASCII:
	 if (isascii(keybyte))
	    keybyte = tolower(keybyte) ;
	 break ;
      case CS_Latin1:
	 keybyte = Fr::tolower(keybyte) ;
	 break ;
      }
#endif
   auto n = node(nodeindex) ;
   auto  index = n->childIndexIfPresent(keybyte) ;
   nodeindex = index ;
   return (index != LangIDPackedMultiTrie::NULL_INDEX) ;
}

//----------------------------------------------------------------------

uint32_t LangIDPackedMultiTrie::extendKey(uint8_t keybyte, uint32_t nodeindex) const
{
   if ((nodeindex & PTRIE_TERMINAL_MASK) != 0)
      {
      return LangIDPackedMultiTrie::NULL_INDEX ;
      }
#if 0 // currently not used, so don't waste time
   if (ignoringWhiteSpace() && keybyte == ' ')
      return true ;
   switch (caseSensitivity())
      {
      case CS_Full:
      default:
	 // do nothing
	 break ;
      case CS_ASCII:
	 if (isascii(keybyte))
	    keybyte = tolower(keybyte) ;
	 break ;
      case CS_Latin1:
	 keybyte = Fr::tolower(keybyte) ;
	 break ;
      }
#endif
   auto n = node(nodeindex) ;
   return n->childIndexIfPresent(keybyte) ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::enumerate(uint8_t *keybuf, unsigned maxkeylength,
				PackedTrieEnumFn *fn, void *user_data) const
{
   if (keybuf && fn && m_nodes && m_nodes[ROOT_INDEX].firstChild())
      {
      memset(keybuf,'\0',maxkeylength) ;
      return enumerateChildren(ROOT_INDEX,keybuf,maxkeylength*8,0,fn,user_data) ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::enumerateChildren(uint32_t nodeindex,
				       uint8_t *keybuf,
				       unsigned max_keylength_bits,
				       unsigned curr_keylength_bits,
				       PackedTrieEnumFn *fn,
				       void *user_data) const
{
   auto n = node(nodeindex) ;
   if (n->leaf() && !fn(n,keybuf,curr_keylength_bits/8,user_data))
      return false ;
   else if (terminalNode(n))
      return true ;
   if (curr_keylength_bits < max_keylength_bits)
      {
      unsigned curr_bits = curr_keylength_bits + PTRIE_BITS_PER_LEVEL ;
      for (unsigned i = 0 ; i < (1<<PTRIE_BITS_PER_LEVEL) ; i++)
	 {
	 uint32_t child = n->childIndexIfPresent(i) ;
	 if (child != LangIDPackedMultiTrie::NULL_INDEX)
	    {
	    unsigned byte = curr_keylength_bits / 8 ;
	    keybuf[byte] = i ;
	    if (!enumerateChildren(child,keybuf,max_keylength_bits,curr_bits,fn,user_data))
	       return false ;
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

LangIDPackedMultiTrie *LangIDPackedMultiTrie::load(Fr::CFile& f, const char *filename)
{
   if (f)
      {
      auto trie = new LangIDPackedMultiTrie(f,filename) ;
      if (!trie || !trie->good())
	 {
	 delete trie ;
	 return nullptr ;
	 }
      return trie ;
      }
   return nullptr ;
}

//----------------------------------------------------------------------

LangIDPackedMultiTrie *LangIDPackedMultiTrie::load(const char *filename)
{
   Fr::CInputFile fp(filename,Fr::CFile::binary) ;
   return fp ? load(fp,filename) : nullptr ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::writeHeader(Fr::CFile& f) const
{
   // write the signature string
   const size_t siglen = sizeof(MULTITRIE_SIGNATURE) ;
   if (f.write(MULTITRIE_SIGNATURE,siglen,sizeof(char)) != siglen)
      return false; 
   // follow with the format version number
   unsigned char version = MULTITRIE_FORMAT_VERSION ;
   unsigned char bits = PTRIE_BITS_PER_LEVEL ;
   if (!f.writeValue(version) || !f.writeValue(bits))
      return false ;
   // write out the size of the trie
   Fr::UInt32 val_used(size()), val_keylen(longestKey()), val_numfreq(m_numfreq), val_numterm(m_numterminals) ;
   char case_sens = caseSensitivity() ;
   if (!f.writeValue(val_used) ||
      !f.writeValue(val_keylen) ||
      !f.writeValue(val_numfreq) ||
      !f.writeValue(val_numterm) ||
      !f.writeValue(m_ignorewhitespace) ||
      !f.writeValue(case_sens))
      return false ;
   // pad the header with NULs for the unused reserved portion of the header
   if (f.putNulls(MULTITRIE_PADBYTES_1))
      f.writeComplete() ;
   return true ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::write(Fr::CFile& f) const
{
   if (!f || !writeHeader(f))
      return false ;
   // write the actual trie nodes
   if (f.write(m_nodes,m_size,sizeof(PackedTrieNode)) != m_size)
      return false ;
   // write the frequency information
   if (f.write(m_freq,m_numfreq,sizeof(PackedTrieFreq)) != m_numfreq)
      return false ;
   // write the terminals
   if (f.write(m_terminals,m_numterminals,sizeof(PackedTrieTerminalNode)) != m_numterminals)
      return false ;
   f.writeComplete() ;
   return true ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::write(const char *filename) const
{
   Fr::COutputFile fp(filename,Fr::CFile::safe_rewrite) ;
   return this->write(fp) ? fp.close() : false ;
}

//----------------------------------------------------------------------

static const PackedTrieFreq *base_frequency ;

static bool dump_ngram(const PackedTrieNode *node, const uint8_t *key,
		       unsigned keylen, void *user_data)
{
   Fr::CFile& f = *((Fr::CFile*)user_data) ;
   if (f && node)
      {
      f.printf("   ") ;
      write_escaped_key(f,key,keylen) ;
      f.printf("  ::") ;
      const PackedTrieFreq *freq = node->frequencies(base_frequency) ;
      for ( ; freq ; freq++)
	 {
	 f.printf(" %lu=%g",(unsigned long)freq->languageID(),freq->probability()) ;
	 if (freq->isLast())
	    break ;
	 }
      f.printf("\n") ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::dump(Fr::CFile& f) const
{
   Fr::LocalAlloc<uint8_t,10000> keybuf(longestKey()) ;
   base_frequency = m_freq ;
   return keybuf ? enumerate(keybuf,longestKey(),dump_ngram,&f) : false ;
}

/************************************************************************/
/*	Additional methods for class MultiTrie				*/
/************************************************************************/

// note: these global variables make add_ngram non-reentrant
static const PackedTrieFreq *frequency_base = nullptr ;
static const PackedTrieFreq *frequency_end = nullptr ;

static bool add_ngram(const PackedTrieNode *node, const uint8_t *key,
		      unsigned keylen, void *user_data)
{
   // not the most efficient method, since it does a separate insertion
   //   for each language ID, but we only use this during training so
   //   speed is not critical
   auto trie = (LangIDMultiTrie*)user_data ;
   const PackedTrieFreq *frequencies = node->frequencies(frequency_base) ;
   if (frequencies)
      {
      for ( ; frequencies < frequency_end ; frequencies++)
	 {
	 trie->insert(key,keylen,frequencies->languageID(),
		      (uint32_t)(frequencies->probability() * TRIE_SCALE_FACTOR + 0.5),
		      frequencies->isStopgram()) ;
	 if (frequencies->isLast())
	    break ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

LangIDMultiTrie::LangIDMultiTrie(const class LangIDPackedMultiTrie *ptrie)
{
   if (ptrie)
      {
      init(ptrie->size() * 3 / 2) ;
      Fr::LocalAlloc<uint8_t> keybuf(ptrie->longestKey()) ;
      if (keybuf)
	 {
	 frequency_base = ptrie->frequencyBaseAddress() ;
	 frequency_end = frequency_base + ptrie->numFrequencies() ;
	 ptrie->enumerate(keybuf,ptrie->longestKey(),add_ngram,this) ;
	 frequency_base = nullptr ;
	 frequency_end = nullptr ;
	 }
      }
   else
      init(1) ;
   return ;
}

// end of file ptrie.C //
