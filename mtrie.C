/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: mtrie.C - bit-slice-based Word-frequency multi-trie		*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-15						*/
/*									*/
/*  (c) Copyright 2011,2012,2015,2019 Carnegie Mellon University	*/
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

#include <cassert>
#include <cstring>
#include "mtrie.h"
#include "framepac/memory.h"
#include "framepac/message.h"
#include "framepac/texttransforms.h"

using namespace std ;
using namespace Fr ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define MULTITRIE_SIGNATURE "MulTrie\0"
#define MULTITRIE_FORMAT_VERSION 2

// reserve some space for future additions to the file format
#define MULTITRIE_PADBYTES_1  64

#if BITS_PER_LEVEL == 3
#  define LEVEL_SIZE 9
#else
#  define LEVEL_SIZE 8
#endif

/************************************************************************/
/*	Types								*/
/************************************************************************/

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

ItemPoolFlat<MultiTrieFrequency> MultiTrieFrequency::s_freq_records(100000) ;

/************************************************************************/
/*	Helper functions						*/
/************************************************************************/

#ifndef lengthof
#  define lengthof(x) (sizeof(x)/sizeof((x)[0]))
#endif /* lengthof */

//----------------------------------------------------------------------

// defined in trie.C
uint32_t scaled_frequency(uint32_t raw_freq, uint64_t total_count) ;

/************************************************************************/
/*	Methods for class MultiTrieFrequency				*/
/************************************************************************/

MultiTrieFrequency::~MultiTrieFrequency()
{
   m_next = 0 ;
   m_frequency = 0 ;
   return ;
}

//----------------------------------------------------------------------

MultiTrieFrequency* MultiTrieFrequency::allocate(uint32_t freq, uint32_t lang, bool stopgram)
{
   size_t recnum = s_freq_records.alloc() ;
   auto freq_record = &s_freq_records[recnum] ;
   new (freq_record) MultiTrieFrequency(freq,lang,stopgram) ;
   return freq_record ;
}

//----------------------------------------------------------------------

uint32_t MultiTrieFrequency::frequency(uint32_t ID) const
{
   if (ID == m_langID)
      return frequency() ;
   else if (next())
      return next()->frequency(ID) ;
   return 0 ;
}

//----------------------------------------------------------------------

void MultiTrieFrequency::newFrequency(uint32_t ID, uint32_t freq, bool stopgram)
{
   // add a record with the new language ID and frequency
   // because the allocation might reallocate the array containing
   //   'this', we need to get our index and use it to re-establish
   //   the correct object address after the allocation
   uint32_t index = (this - s_freq_records.item(0)) ;
   MultiTrieFrequency *f = allocate(freq,ID,stopgram) ;
   MultiTrieFrequency *prev = getAddress(index) ;
   prev->setNext(f) ;

}

//----------------------------------------------------------------------

void MultiTrieFrequency::setFrequency(uint32_t ID, uint32_t freq, bool stopgram)
{
   if (ID == languageID())
      {
      setFrequency(freq) ;
      return ;
      }
   MultiTrieFrequency *nxt = next() ;
   if (nxt)
      {
      nxt->setFrequency(ID,freq,stopgram) ;
      }
   else
      {
      newFrequency(ID,freq,stopgram) ;
      }
   return ;
}

//----------------------------------------------------------------------

void MultiTrieFrequency::incrFrequency(uint32_t ID, uint32_t incr)
{
   if (ID == languageID())
      {
      incrFrequency(incr) ;
      return ;
      }
   MultiTrieFrequency *nxt = next() ;
   if (nxt)
      {
      nxt->incrFrequency(ID,incr) ;
      }
   else
      {
      // add a record with the new ID
      newFrequency(ID,incr,false) ;
      }
   return ;
}

//----------------------------------------------------------------------

void MultiTrieFrequency::scaleFrequency(uint64_t total_count,
					uint32_t id)
{
   if (languageID() == id)
      m_frequency = scaled_frequency(m_frequency,total_count) ;
   else if (next())
      next()->scaleFrequency(total_count,id) ;
   return ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::readAll(CFile& f)
{
   UInt32 cnt ;
   if (f && f.readValue(&cnt))
      {
      if (s_freq_records.load(f,cnt.load()))
	 return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::write(CFile& f) const
{
   return f && f.writeValue(*this) ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::writeAll(CFile& f)
{
   UInt32 count(s_freq_records.size()) ;
   return f && f.writeValue(count) && s_freq_records.save(f) ;
}

/************************************************************************/
/*	Methods for class MultiTrieNode					*/
/************************************************************************/

MultiTrieNode::MultiTrieNode()
{
   m_frequency = LangIDMultiTrie::INVALID_FREQ ;
   return ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::isStopgram(unsigned langID) const
{
   auto freq = frequencies() ;
   for ( ; freq ; freq = freq->next())
      {
      if (freq->languageID() == langID)
	 return freq->isStopgram() ;
      }
   return false ;
}

//----------------------------------------------------------------------

unsigned MultiTrieNode::numFrequencies() const
{
   auto freq = frequencies() ;
   unsigned count = 0 ;
   for ( ; freq ; freq = freq->next())
      count++ ;
   return count ;
}

//----------------------------------------------------------------------

uint32_t MultiTrieNode::frequency(uint32_t langID) const
{
   auto freq = frequencies() ;
   for ( ; freq ; freq = freq->next())
      {
      if (freq->languageID() == langID)
	 return freq->frequency() ;
      }
   return 0 ;
}

//----------------------------------------------------------------------

void MultiTrieNode::setFrequency(uint32_t langID, uint32_t freq, bool stopgram)
{
   auto f = MultiTrieFrequency::getAddress(m_frequency) ;
   if (f)
      f->setFrequency(langID,freq,stopgram) ;
   else
      {
      // insert the initial frequency record
      f = MultiTrieFrequency::allocate(freq,langID,stopgram) ;
      m_frequency = MultiTrieFrequency::getIndex(f) ;
      }
   return ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::setFrequencies(MultiTrieFrequency *freqs)
{
   auto idx = MultiTrieFrequency::getIndex(freqs) ;
   if (idx == LangIDMultiTrie::INVALID_FREQ)
      return false ;
   m_frequency = idx ;
   return true ;
}

//----------------------------------------------------------------------

uint32_t MultiTrieNode::insertChild(unsigned int N, LangIDMultiTrie *trie)
{
   if (N < lengthof(m_children))
      {
      if (childPresent(N))
	 return childIndex(N) ;
      else
	 {
	 uint32_t new_index = trie->allocateNode() ;
	 if (new_index)
	    {
	    m_children[N] = new_index ;
	    return new_index ;
	    }
	 }
      }
   return NybbleTrie::INVALID_INDEX ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::load(CFile& f)
{
   return f && f.readValue(this) ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::write(CFile& f) const
{
   return f && f.writeValue(*this) ;
}

/************************************************************************/
/*	Methods for class LangIDMultiTrie				*/
/************************************************************************/

LangIDMultiTrie::LangIDMultiTrie(uint32_t cap)
{
   init(cap) ;
   return ;
}

//----------------------------------------------------------------------

LangIDMultiTrie::LangIDMultiTrie(const char *filename, uint32_t langID, bool verbose)
{
   init(1) ;
   LoadFn* insfn = [](NybbleTrie* trie, const uint8_t* key, unsigned keylen, uint32_t lang, uint32_t freq) -> bool
		      { return static_cast<LangIDMultiTrie*>(trie)->insert(key,keylen,lang,freq,false) ; } ;
   loadWords(filename,insfn,langID,verbose) ;
   return ;
}

//----------------------------------------------------------------------

void LangIDMultiTrie::init(uint32_t cap)
{
   MultiTrieFrequency::allocateDummy() ;
   m_maxkeylen = 0 ;
   m_totaltokens = 0 ;
   m_ignorewhitespace = false ;
   if (cap == 0)
      cap = 1 ;
   m_nodes.reserve(cap) ;
   auto root = m_nodes.alloc() ;
   // initialize the root node
   new (node(root)) MultiTrieNode ;
   return ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::insert(const uint8_t *key, unsigned keylength,
		       uint32_t langID, uint32_t frequency, bool stopgram)
{
   auto leaf = node(insertKey(key,keylength)) ;
   bool new_node = false ;
   if (leaf)
      {
      new_node = (leaf->frequency() == 0) ;
      leaf->setFrequency(langID,frequency,stopgram) ;
      leaf->markAsLeaf() ;
      if (frequency > 0)
	 {
	 addTokenCount() ;
	 }
      }
   return new_node ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::increment(const uint8_t *key, unsigned keylength,
			      uint32_t langID, uint32_t incr,
			      bool stopgram)
{
   uint32_t cur_index = ROOT_INDEX ;
   for (size_t i = 0 ; i < keylength ; i++)
      {
      if (!extendKey(cur_index,key[i]))
	 {
	 insert(key,keylength,langID,incr,stopgram) ;
	 return incr ;
	 }
      }
   auto n = node(cur_index) ;
   if (n)
      {
      uint32_t freq = n->frequency(langID) + incr ;
      n->setFrequency(langID,freq,n->isStopgram(langID)) ;
      return freq ;
      }
   else
      {
      insert(key,keylength,langID,incr,stopgram) ;
      return incr ;
      }
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::incrementExtensions(const uint8_t *key, unsigned prevlength,
				    unsigned keylength, uint32_t langID,
				    uint32_t incr)
{
   uint32_t cur_index = ROOT_INDEX ;
   // check whether the prevlength prefix is present in the trie
   for (size_t i = 0 ; i < prevlength ; i++)
      {
      if (!extendKey(cur_index,key[i]))
	 return false ;
      }
   // now add on one byte at a time, incrementing the count for each
   for (size_t i = prevlength ; i < keylength ; i++)
      {
      this->insertChild(cur_index,key[i]) ;
      auto n = node(cur_index) ;
      if (!n)
	 return false ;
      uint32_t freq = n->frequency(langID) + incr ;
      n->setFrequency(langID,freq,n->isStopgram(langID)) ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerate(uint8_t *keybuf, unsigned maxkeylength, EnumFn *fn, void *user_data) const
{
   if (keybuf && fn)
      {
      std::fill_n(keybuf,maxkeylength,'\0') ;
      return enumerateChildren(ROOT_INDEX,keybuf,maxkeylength*8,0,fn,user_data) ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerateChildren(NodeIndex nodeindex,
				      uint8_t *keybuf,
				      unsigned max_keylength_bits,
				      unsigned curr_keylength_bits,
   				      LangIDMultiTrie::EnumFn *fn,
				      void *user_data) const
{
   auto n = node(nodeindex) ;
//   assert(!(!n->leaf() && n->frequencies())) ;
   if (n->leaf() && !fn(this,nodeindex,keybuf,curr_keylength_bits/8,user_data))
      return false ;
   if (curr_keylength_bits < max_keylength_bits)
      {
      unsigned curr_bits = curr_keylength_bits + BITS_PER_LEVEL ;
      for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
	 {
	 uint32_t child = n->childIndex(i) ;
	 if (child != LangIDMultiTrie::NULL_INDEX)
	    {
	    auto childnode = node(child) ;
	    if (childnode)
	       {
	       unsigned byte = curr_keylength_bits / 8 ;
	       unsigned shift = LEVEL_SIZE - (curr_keylength_bits%8) - BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
	       if (shift == 0) curr_bits = curr_keylength_bits + BITS_PER_LEVEL - 1 ;
#endif
	       unsigned mask = (((1<<BITS_PER_LEVEL)-1) << shift) ;
	       keybuf[byte] &= ~mask ;
	       keybuf[byte] |= (i << shift) ;
	       if (!enumerateChildren(child,keybuf,max_keylength_bits,curr_bits,fn,user_data))
		  return false ;
	       }
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

size_t LangIDMultiTrie::countTerminalNodes(NodeIndex nodeindex, unsigned keylen_bits) const
{
   size_t count = 0 ;
   auto n = node(nodeindex)   ;
   if (!n->hasChildren())
      {
      return (keylen_bits % 8 == 0) ? 1 : 0 ;
      }
   keylen_bits += BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) --keylen_bits ;
#endif
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
      {
      uint32_t child = n->childIndex(i) ;
      if (child != NULL_INDEX)
	 count += countTerminalNodes(child,keylen_bits) ;
      }
   return count ;
}

//----------------------------------------------------------------------

size_t LangIDMultiTrie::countFullByteNodes(NodeIndex nodeindex, unsigned keylen_bits) const
{
   size_t count = 0 ;
   if (keylen_bits % 8 == 0)
      count++ ;
   keylen_bits = keylen_bits + BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   auto n = node(nodeindex) ;
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
      {
      uint32_t child = n->childIndex(i) ;
      if (child != NULL_INDEX)
	 count += countFullByteNodes(child,keylen_bits) ;
      }
   return count ;
}

//----------------------------------------------------------------------

unsigned LangIDMultiTrie::numExtensions(NodeIndex nodeindex, unsigned bits) const
{
   if (bits >= 8)
      {
      return 1 ;
      }
   auto n = node(nodeindex) ;
   unsigned count = 0 ;
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; ++i)
      {
      auto child = n->childIndex(i) ;
      if (child != NULL_INDEX)
	 count += numExtensions(child,bits+BITS_PER_LEVEL) ;
      }
   return count ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::allChildrenAreTerminals(NodeIndex nodeindex, unsigned bits) const
{
   auto n = node(nodeindex) ;
   if (bits >= 8)
      {
      return !n->hasChildren() ;
      }
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; ++i)
      {
      auto child = n->childIndex(i) ;
      if (child != NULL_INDEX && !allChildrenAreTerminals(child,bits+BITS_PER_LEVEL))
	 return false ;
      }
   return true ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::numFullByteNodes() const
{
   return countFullByteNodes(ROOT_INDEX) ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::numTerminalNodes() const
{
   return countTerminalNodes(ROOT_INDEX) ;
}

//----------------------------------------------------------------------

static bool count_freq_records(const LangIDMultiTrie* trie, NybbleTrie::NodeIndex nodeindex, const uint8_t *,
			       unsigned, void *user_data)
{
   auto node = trie->node(nodeindex) ;
   auto count = reinterpret_cast<uint32_t*>(user_data) ;
   (*count) += node->numFrequencies() ;
   return true ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::countFreqRecords() const
{
   uint32_t count = 0 ;
   uint8_t keybuf[longestKey()] ; 
   (void)this->enumerate(keybuf,longestKey(),count_freq_records,&count) ;
   return count ;
}

//----------------------------------------------------------------------

LangIDMultiTrie *LangIDMultiTrie::load(CFile& f)
{
   if (!f)
      return nullptr ;
   int version = f.verifySignature(MULTITRIE_SIGNATURE) ;
   if (version < 0)
      {
      if (version == -1) { /* read error */ }
      if (version == -2) { /* wrong file type */ }
      if (version == -3) { /* wrong byte order */ }
      return nullptr ;
      }
   if (version != MULTITRIE_FORMAT_VERSION)
      {
      // error: wrong version of data file
      return nullptr ;
      }
   unsigned char bits ;
   if (!f.readValue(&bits) || bits != BITS_PER_LEVEL)
      {
      // error: wrong type of trie
      return nullptr ;
      }
   UInt32 val_used, val_tokens, val_keylen ;
   if (!f.readValue(&val_used) ||
      !f.readValue(&val_tokens) ||
      !f.readValue(&val_keylen))
      {
      // error reading header
      return nullptr ;
      }
   uint32_t used = val_used.load() ;
   Owned<LangIDMultiTrie> trie(used) ;
   if (!trie)
      return nullptr ;
   trie->m_nodes.allocBatch(used) ;
   trie->addTokenCount(val_tokens.load()) ;
   trie->m_maxkeylen = (unsigned)val_keylen.load() ;
   // skip the padding (reserved bytes)
   f.seek(MULTITRIE_PADBYTES_1,SEEK_CUR) ;
   // read the actual trie nodes
   for (size_t i = 0 ; i < used ; i++)
      {
      auto node = trie->node(i) ;
      if (!node->load(f))
	 {
	 return nullptr ;
	 }
      }
   // finally, read the frequency information
   if (!MultiTrieFrequency::readAll(f))
      {
      return nullptr ;
      }
   return trie.move() ;
}

//----------------------------------------------------------------------

LangIDMultiTrie *LangIDMultiTrie::load(const char *filename)
{
   CInputFile fp(filename) ;
   return load(fp) ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::writeHeader(CFile& f) const
{
   // write the signature string
   if (!f.writeSignature(MULTITRIE_SIGNATURE,MULTITRIE_FORMAT_VERSION))
      return false ;
   // follow with the number of bits per level of the trie
   unsigned char bits = BITS_PER_LEVEL ;
   if (!f.writeValue(bits))
      return false ;
   // write out the size of the trie
   UInt32 val_used(size()), val_tokens(totalTokens()), val_keylen(longestKey()) ;
   if (!f.writeValue(val_used) ||
      !f.writeValue(val_tokens) ||
      !f.writeValue(val_keylen))
      return false ;
   // pad the header with NULs for the unused reserved portion of the header
   return f.putNulls(MULTITRIE_PADBYTES_1) ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::write(CFile& f) const
{
   if (!f || !writeHeader(f))
      return false ;
   // write the actual trie nodes
   for (size_t i = 0 ; i < size() ; i++)
      {
      const auto n = node(i) ;
      if (!n->write(f))
	 return false ;
      }
   // write the frequency information
   if (!MultiTrieFrequency::writeAll(f))
      return false ;
   f.writeComplete() ;
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::write(const char *filename) const
{
   COutputFile fp(filename,CFile::safe_rewrite) ;
   return this->write(fp) ? fp.close() : false ;
}

//----------------------------------------------------------------------

static const char hexdigit[] = "0123456789ABCDEF" ;

void write_escaped_char(uint8_t c, CFile& f)
{
   switch (c)
      {
      case '\\':
	 f.puts("\\\\") ;
	 break ;
      case ' ':
	 f.puts("\\ ") ;
	 break ;
      default:
	 if (c < ' ')
	    {
	    f.putc('\\') ;
	    f.putc(hexdigit[(c>>4)&0xF]) ;
	    f.putc(hexdigit[c&0xF]) ;
	    }
	 else
	    f.putc(c) ;
	 break ;
      }
   return ;
}

//----------------------------------------------------------------------

void write_escaped_key(CFile& f, const uint8_t* key, unsigned keylen)
{
   for (size_t i = 0 ; i < keylen ; i++)
      {
      write_escaped_char(key[i],f) ;
      }
   return ;
}

//----------------------------------------------------------------------

static bool dump_ngram(const LangIDMultiTrie* trie, NybbleTrie::NodeIndex nodeindex, const uint8_t *key,
		       unsigned keylen, void *user_data)
{
   auto node = trie->node(nodeindex) ;
   CFile& f = *((CFile*)user_data) ;
   if (f && node)
      {
      f.printf("   ") ;
      write_escaped_key(f,key,keylen) ;
      f.printf("  ::") ;
      auto freq = node->frequencies() ;
      for ( ; freq ; freq = freq->next())
	 {
	 f.printf(" %lu=%lu",(unsigned long)freq->languageID(), (unsigned long)freq->frequency()) ;
	 }
      f.printf("\n") ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::dump(CFile& f) const
{
   LocalAlloc<uint8_t,10000> keybuf(longestKey()) ;
   return keybuf ? enumerate(keybuf,longestKey(),dump_ngram,&f) : false ;
}

// end of file mtrie.C //
