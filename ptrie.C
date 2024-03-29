/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: ptrie.C - packed Word-frequency multi-trie			*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-18						*/
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

#include <cstring>
#include "mtrie.h"
#include "ptrie.h"
#include "framepac/file.h"
#include "framepac/message.h"
#include "framepac/utility.h"

using namespace std ;
using namespace Fr ;

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
#define MULTITRIE_PADBYTES_1  59

/************************************************************************/
/*	Types								*/
/************************************************************************/

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

double PackedTrieFreq::s_value_map[PackedTrieFreq::TRIE_NUM_VALUES] ;
bool PackedTrieFreq::s_value_map_initialized = false ;

//----------------------------------------------------------------------

void write_escaped_key(CFile& f, const uint8_t* key, unsigned keylen) ;

/************************************************************************/
/*	Helper functions						*/
/************************************************************************/

/************************************************************************/
/*	Methods for class PackedTrieFreq				*/
/************************************************************************/

PackedTrieFreq::PackedTrieFreq(uint32_t freq, uint32_t langID, bool last, bool is_stop)
{
   uint32_t data = (langID | (last * TRIE_LASTENTRY) | (is_stop * TRIE_STOPGRAM)) ;
   uint32_t mant ;
   uint32_t expon ;
   quantize(freq,mant,expon) ;
   data |= mant ;
   data |= (expon << TRIE_FREQ_EXP_SHIFT) ;
   m_freqinfo.store(data) ;
   return ;
}

//----------------------------------------------------------------------

void PackedTrieFreq::quantize(uint32_t freq, uint32_t &mantissa, uint32_t &exp)
{
   uint32_t e = 0 ;
   if (freq)
      {
      constexpr uint32_t max_exponent = TRIE_FREQ_EXPONENT >> TRIE_FREQ_EXP_SHIFT ;
      while ((freq & TRIE_FREQ_HIBITS) == 0 && e < max_exponent)
	 {
	 freq <<= EXPONENT_SCALE ;
	 e++ ;
	 }
      freq &= TRIE_FREQ_MANTISSA ;
      if (freq == 0)
	 freq = TRIE_MANTISSA_LSB ;
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
   uint32_t data = m_freqinfo.load() & ~TRIE_LASTENTRY ;
   if (last)
      data |= TRIE_LASTENTRY ;
   m_freqinfo.store(data) ;
   return ;
}

//----------------------------------------------------------------------

void PackedTrieFreq::initDataMapping(double (*mapfn)(uint32_t))
{
   for (size_t i = 0 ; i < TRIE_NUM_VALUES ; i++)
      {
      uint32_t scaled = scaledScore(i << TRIE_VALUE_SHIFT) ;
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

bool PackedTrieFreq::writeDataMapping(CFile& f)
{
   UInt32 count(TRIE_NUM_VALUES) ;
   return f && f.writeValue(count) && f.writeValue(*s_value_map) ;
}

/************************************************************************/
/*	Methods for class PackedTrieNode				*/
/************************************************************************/

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
   return (firstChild() + m_popcounts[N/32] + popcount(children)) ;
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
   return (firstChild() + m_popcounts[N/32] + popcount(children)) ;
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
   unsigned pcount = 0 ;
   for (size_t i = 0 ; i < lengthof(m_popcounts) ; i++)
      {
      m_popcounts[i] = (uint8_t)pcount ;
      uint32_t children = m_children[i].load() ;
      pcount += popcount(children) ;
      }
   return ;
}

/************************************************************************/
/*	Methods for class PackedMultiTrie				*/
/************************************************************************/

LangIDPackedMultiTrie::LangIDPackedMultiTrie(const LangIDMultiTrie *multrie)
{
   if (multrie)
      {
      auto numterminals = multrie->numTerminalNodes() ;
      auto sz = multrie->numFullByteNodes() - numterminals ;
      m_terminals.reserve(numterminals) ;
      m_nodes.reserve(sz) ;
      m_freq.reserve(multrie->countFreqRecords()) ;
      if (m_nodes.capacity() && m_freq.capacity())
	 {
	 auto proot = m_nodes.alloc() ;
	 if (!insertChildren(node(proot),multrie,LangIDMultiTrie::ROOT_INDEX))
	    {
	    m_nodes.clear() ;
	    m_freq.clear() ;
	    m_terminals.clear() ;
	    }
	 SystemMessage::status("   converted %lu full nodes, %lu terminals, and %lu frequencies",
	    m_nodes.size(), m_terminals.size(), m_freq.size()) ;
	 }
      else
	 {
	 m_nodes.clear() ;
	 m_freq.clear() ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

LangIDPackedMultiTrie::LangIDPackedMultiTrie(CFile& f, const char *filename)
{
   size_t numfull ;
   size_t numfreq ;
   size_t numterminals ;
   if (f && parseHeader(f,numfull,numfreq,numterminals))
      {
      auto offset = f.tell() ;
      m_fmap.open(filename) ;
      if (m_fmap)
	 {
	 // we can memory-map the file, so just point our member variables
	 //   at the mapped data
	 const char* base = *m_fmap + offset ;
	 m_nodes.external_buffer(base,numfull) ;
	 base += numfull * sizeof(PackedTrieNode) ;
	 m_freq.external_buffer(base,numfreq) ;
	 base += numfreq * sizeof(PackedTrieFreq) ;
	 m_terminals.external_buffer(base,numterminals) ;
	 }
      else
	 {
	 // unable to memory-map the file, so read its contents into buffers
	 //   and point our variables at the buffers
	 if (!m_nodes.load(f,numfull) || !m_freq.load(f,numfreq) || !m_terminals.load(f,numterminals))
	    {
	    m_nodes.clear() ;
	    m_freq.clear() ;
	    m_terminals.clear() ;
	    }
	 }
      }
   return ;
}

//----------------------------------------------------------------------

uint32_t LangIDPackedMultiTrie::allocateChildNodes(unsigned numchildren)
{
   return m_nodes.allocBatch(numchildren) ;
}

//----------------------------------------------------------------------

uint32_t LangIDPackedMultiTrie::allocateTerminalNodes(unsigned numchildren)
{
   return m_terminals.allocBatch(numchildren) | TERMINAL_MASK ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::insertTerminals(PackedTrieNode *parent,
				      const LangIDMultiTrie *mtrie,
				      uint32_t mnode_index,
				      unsigned keylen)
{
   if (!parent)
      return false ;
   unsigned numchildren = mtrie->numExtensions(mnode_index) ;
   if (numchildren == 0)
      return true ;
   keylen++ ;
   if (keylen > longestKey())
      m_maxkeylen = keylen ;
   uint32_t firstchild = allocateTerminalNodes(numchildren) ;
   parent->setFirstChild(firstchild) ;
   if (firstchild == NOCHILD_INDEX)
      {
      SystemMessage::error("insertTerminals: firstchild==NOCHILD_INDEX") ;
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
	    uint32_t freq_index = m_freq.allocBatch(numfreq) ;
	    pchild->setFrequencies(freq_index) ;
	    auto mfreq = mchild->frequencies() ;
	    while (mfreq && numfreq > 0)
	       {
	       (void)new (&m_freq[freq_index]) PackedTrieFreq(mfreq->frequency(),mfreq->languageID(),numfreq <= 1) ;
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
				     uint32_t mnode_index,
				     unsigned keylen)
{
   if (!parent)
      return false ;
   // first pass: fill in all the children
   unsigned numchildren = mtrie->numExtensions(mnode_index) ;
   if (numchildren == 0)
      return true ;
   keylen++ ;
   if (keylen > longestKey())
      m_maxkeylen = keylen ;
   bool terminal = mtrie->allChildrenAreTerminals(mnode_index) ;
   uint32_t firstchild = (terminal
			  ? allocateTerminalNodes(numchildren)
			  : allocateChildNodes(numchildren)) ;
   parent->setFirstChild(firstchild) ;
   if (firstchild == NOCHILD_INDEX)
      {
      SystemMessage::error("insertChildren: firstchild==NOCHILD_INDEX") ;
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
	    uint32_t freq_index = m_freq.allocBatch(numfreq) ;
	    pchild->setFrequencies(freq_index) ;
	    auto mfreq = mchild->frequencies() ;
	    while (mfreq && numfreq > 0)
	       {
	       bool is_stop = (mfreq->isStopgram() || mfreq->frequency() == 0);
	       (void)new (&m_freq[freq_index]) PackedTrieFreq
		  (mfreq->frequency(),mfreq->languageID(),numfreq <= 1,is_stop) ;
	       freq_index++ ;
	       numfreq-- ;
	       mfreq = mfreq->next() ;
	       }
	    }
	 if (terminal)
	    {
	    if (!insertTerminals(pchild,mtrie,nodeindex,keylen))
	       return false ;
	    }
	 else if (!insertChildren(pchild,mtrie,nodeindex,keylen))
	    return false ;
	 }
      }
   parent->setPopCounts() ;
   return true ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::parseHeader(CFile& f, size_t& numfull, size_t& numfreq, size_t& numterminals)
{
   int version = f.verifySignature(MULTITRIE_SIGNATURE) ;
   if (version < 0)
      {
      if (version == -1) { /* read error */ }
      if (version == -2) { /* error: wrong file type */ }
      if (version == -3) { /* error: wrong byte order */ }
      return false ;
      }
   if (version < MULTITRIE_FORMAT_MIN_VERSION || version > MULTITRIE_FORMAT_VERSION)
      {
      // error: wrong version of data file
      return false ;
      }
   unsigned char bits ;
   if (!f.readValue(&bits) || bits != PTRIE_BITS_PER_LEVEL)
      {
      // error: wrong type of trie
      return false ;
      }
   char ignore_white ;
   char case_sens ;
   char padbuf[MULTITRIE_PADBYTES_1] ;
   UInt32 val_size, val_keylen, val_numfreq, val_numterm ;
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
   numfull = val_size.load() ;
   numterminals = val_numterm.load() ;
   numfreq = val_numfreq.load() ;
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
   if ((nodeindex & TERMINAL_MASK) != 0)
      {
      nodeindex = NULL_INDEX ;
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
   return (index != NULL_INDEX) ;
}

//----------------------------------------------------------------------

uint32_t LangIDPackedMultiTrie::extendKey(uint8_t keybyte, uint32_t nodeindex) const
{
   if ((nodeindex & TERMINAL_MASK) != 0)
      {
      return NULL_INDEX ;
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

bool LangIDPackedMultiTrie::enumerate(uint8_t *keybuf, unsigned maxkeylength, EnumFn *fn, void *user_data) const
{
   if (keybuf && fn && m_nodes && m_nodes[ROOT_INDEX].firstChild())
      {
      std::fill_n(keybuf,maxkeylength,'\0') ;
      return enumerateChildren(ROOT_INDEX,keybuf,maxkeylength*8,0,fn,user_data) ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::enumerateChildren(uint32_t nodeindex,
				       uint8_t *keybuf,
				       unsigned max_keylength_bits,
				       unsigned curr_keylength_bits,
				       EnumFn *fn, void *user_data) const
{
   auto n = node(nodeindex) ;
   if (n->leaf() && !fn(n,keybuf,curr_keylength_bits/8,user_data))
      return false ;
   else if (terminalNode(nodeindex))
      return true ;
   if (curr_keylength_bits < max_keylength_bits)
      {
      unsigned curr_bits = curr_keylength_bits + PTRIE_BITS_PER_LEVEL ;
      for (unsigned i = 0 ; i < (1<<PTRIE_BITS_PER_LEVEL) ; i++)
	 {
	 uint32_t child = n->childIndexIfPresent(i) ;
	 if (child != NULL_INDEX)
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

Owned<LangIDPackedMultiTrie> LangIDPackedMultiTrie::load(CFile& f, const char *filename)
{
   if (f)
      {
      //(if we use Owned trie(f,filename), template deduction tries to send 'f' by value instead of reference...)
      auto trie = new LangIDPackedMultiTrie(f,filename) ;
      if (trie && trie->good())
	 return trie ;
      }
   return nullptr ;
}

//----------------------------------------------------------------------

Owned<LangIDPackedMultiTrie> LangIDPackedMultiTrie::load(const char *filename)
{
   CInputFile fp(filename,CFile::binary) ;
   return fp ? load(fp,filename) : nullptr ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::writeHeader(CFile& f) const
{
   // write the signature string
   if (!f.writeSignature(MULTITRIE_SIGNATURE,MULTITRIE_FORMAT_VERSION))
      return false ;
   // follow with the number of bits per level of the trie
   unsigned char bits = PTRIE_BITS_PER_LEVEL ;
   if (!f.writeValue(bits))
      return false ;
   // write out the size of the trie
   UInt32 val_used(size()) ;
   UInt32 val_keylen(longestKey()) ;
   UInt32 val_numfreq(m_freq.size());
   UInt32 val_numterm(m_terminals.size()) ;
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

bool LangIDPackedMultiTrie::write(CFile& f) const
{
   if (!f || !writeHeader(f))
      return false ;
   // write the actual trie nodes
   if (!m_nodes.save(f))
      return false ;
   // write the frequency information
   if (!m_freq.save(f))
      return false ;
   // write the terminals
   if (!m_terminals.save(f))
      return false ;
   f.writeComplete() ;
   return true ;
}

//----------------------------------------------------------------------

bool LangIDPackedMultiTrie::write(const char *filename) const
{
   COutputFile fp(filename,CFile::safe_rewrite) ;
   return this->write(fp) ? fp.close() : false ;
}

//----------------------------------------------------------------------

static const PackedTrieFreq *base_frequency ;

static bool dump_ngram(const PackedTrieNode *node, const uint8_t *key,
		       unsigned keylen, void *user_data)
{
   CFile& f = *((CFile*)user_data) ;
   if (f && node)
      {
      f.printf("   ") ;
      write_escaped_key(f,key,keylen) ;
      f.printf("  ::") ;
      auto freq = node->frequencies(base_frequency) ;
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

bool LangIDPackedMultiTrie::dump(CFile& f) const
{
   LocalAlloc<uint8_t,10000> keybuf(longestKey()) ;
   base_frequency = m_freq.item(0) ;
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
   auto frequencies = node->frequencies(frequency_base) ;
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
      LocalAlloc<uint8_t> keybuf(ptrie->longestKey()) ;
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
