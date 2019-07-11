/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: mtrie.C - bit-slice-based Word-frequency multi-trie		*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-07						*/
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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include "mtrie.h"
#include "framepac/byteorder.h"
#include "framepac/file.h"
#include "framepac/memory.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define BUCKET_SIZE 65536	// must be power of 2

// we want to store percentages for entries in the trie in 32 bits.  Since
//   it is very unlikely that any ngram in the trie will have a probability
//   greater than 4.2%, scale the percentage by a factor of one billion.
#define SCALE_FACTOR 1000000000L

#define MULTITRIE_SIGNATURE "MulTrie\0"
#define MULTITRIE_FORMAT_VERSION 1

// reserve some space for future additions to the file format
#define MULTITRIE_PADBYTES_1  64

#if MTRIE_BITS_PER_LEVEL == 3
#  define LEVEL_SIZE 9
#else
#  define LEVEL_SIZE 8
#endif

/************************************************************************/
/*	Types								*/
/************************************************************************/

typedef char LONGbuffer[4] ;

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

MultiTrieFrequency *MultiTrieFrequency::s_base_address = nullptr ;
uint32_t MultiTrieFrequency::s_max_alloc = 0 ;
uint32_t MultiTrieFrequency::s_curr_alloc = 0 ;

/************************************************************************/
/*	Helper functions						*/
/************************************************************************/

#ifndef lengthof
#  define lengthof(x) (sizeof(x)/sizeof((x)[0]))
#endif /* lengthof */

//----------------------------------------------------------------------

#ifndef round_up
inline uint32_t round_up(uint32_t value, uint32_t granularity)
{
   return granularity * ((value + granularity - 1) / granularity) ;
}
#endif

//----------------------------------------------------------------------

static uint32_t scaled_frequency(uint32_t raw_freq, uint64_t total_count)
{
   double percent = 100.0 * raw_freq / (double)total_count ;
   // avoid overflow by truncating excessively high percentages to the
   //   largest value representable in a uint32_t
   if (percent > ((uint32_t)~0) / (double)SCALE_FACTOR)
      return (uint32_t)~0 ;
   uint32_t scaled = (uint32_t)(SCALE_FACTOR * percent + 0.5) ;
   // avoid truncation to zero for very low percentages
   return (percent > 0.0 && scaled == 0) ? 1 : scaled ;
}

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
   if (s_curr_alloc >= s_max_alloc)
      {
      uint32_t new_alloc = s_max_alloc ? 2 * s_max_alloc : 1000000 ;
      auto new_base = Fr::New<MultiTrieFrequency>(new_alloc) ;
      if (new_base)
	 {
	 std::copy(s_base_address,s_base_address+s_max_alloc,new_base) ;
	 s_base_address = new_base ;
	 s_max_alloc = new_alloc ;
	 }
      else
	 return nullptr ;
      }
   MultiTrieFrequency *freq_record = baseAddress() + s_curr_alloc ;
   s_curr_alloc++ ;
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

void MultiTrieFrequency::setFrequency(uint32_t ID, uint32_t freq,
				      bool stopgram)
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
      // add a record with the new language ID and frequency
      // because the allocation might reallocate the array containing
      //   'this', we need to get our index and use it to re-establish
      //   the correct object address after the allocation
      uint32_t index = (this - baseAddress()) ;
      MultiTrieFrequency *f = allocate(freq,ID,stopgram) ;
      MultiTrieFrequency *prev = getAddress(index) ;
      prev->setNext(f) ;
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
      auto f = allocate(incr,ID) ;
      setNext(f) ;
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

MultiTrieFrequency *MultiTrieFrequency::read(Fr::CFile& f)
{
   if (f)
      {
      char buffer[sizeof(MultiTrieFrequency)] ;
      if (f.read(buffer,sizeof(buffer),sizeof(char)) == sizeof(buffer))
	 {
	 MultiTrieFrequency *freq = (MultiTrieFrequency*)buffer ;
	 return new MultiTrieFrequency(*freq) ;
	 }
      }
   return nullptr ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::readAll(Fr::CFile& f)
{
   Fr::UInt32 cnt ;
   if (f && f.readValue(&cnt))
      {
      uint32_t count = cnt.load() ;
      MultiTrieFrequency* base = nullptr ;
      if (f.readValues(&base,count))
	 {
	 Fr::Free(s_base_address) ;
	 s_base_address = base ;
	 s_max_alloc = s_curr_alloc = count ;
	 return true ;
	 }
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::write(Fr::CFile& f) const
{
   return f && f.writeValue(*this) ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::writeAll(Fr::CFile& f)
{
   Fr::UInt32 count(s_curr_alloc) ;
   return (f &&
           f.writeValue(count) &&
           f.write(baseAddress(),s_curr_alloc,sizeof(MultiTrieFrequency)) == s_curr_alloc) ;
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

void MultiTrieNode::setFrequency(uint32_t langID, uint32_t freq,
				 bool stopgram)
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
   if (freqs < MultiTrieFrequency::baseAddress())
      return false ;
   auto idx = MultiTrieFrequency::getIndex(freqs) ;
   if (idx == LangIDMultiTrie::INVALID_FREQ)
      return false ;
   m_frequency = idx ;
   return true ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::insertChild(unsigned int N, LangIDMultiTrie *trie)
{
   if (N < lengthof(m_children) && !childPresent(N))
      {
      uint32_t new_index = trie->allocateNode() ;
      if (new_index)
	 {
	 m_children[N] = new_index ;
	 return true ;
	 }
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::load(Fr::CFile& f)
{
   return f && f.readValue(this) ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::write(Fr::CFile& f) const
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

LangIDMultiTrie::LangIDMultiTrie(const char *filename, bool verbose)
{
   init(1) ;
   loadWords(filename,verbose) ;
   return ;
}

//----------------------------------------------------------------------

LangIDMultiTrie::~LangIDMultiTrie()
{
   unsigned num_buckets = m_capacity / BUCKET_SIZE ;
   for (size_t i = 0 ; i < num_buckets ; i++)
      {
      Fr::Free(m_nodes[i]) ;
      }
   Fr::Free(m_nodes) ;
   m_nodes = nullptr ;
   m_capacity = 0 ;
   m_used = 0 ;
   return ;
}

//----------------------------------------------------------------------

void LangIDMultiTrie::init(uint32_t cap)
{
   m_maxkeylen = 0 ;
   m_totaltokens = 0 ;
   m_ignorewhitespace = false ;
   if (cap == 0)
      cap = 1 ;
   cap = round_up(cap,BUCKET_SIZE) ;
   unsigned num_buckets = cap / BUCKET_SIZE ;
   m_nodes = Fr::New<MultiTrieNode*>(num_buckets) ;
   m_capacity = 0 ;
   if (m_nodes)
      {
      for (unsigned i = 0 ; i < num_buckets ; i++)
	 {
	 m_nodes[i] = Fr::New<MultiTrieNode>(BUCKET_SIZE) ;
	 if (m_nodes[i])
	    m_capacity += BUCKET_SIZE ;
	 else
	    break ;
	 }
      m_used = 1 ;
      // initialize the root node
      new (node(LangIDMultiTrie::ROOT_INDEX)) MultiTrieNode ;
      }
   return ;
}

//----------------------------------------------------------------------
// NOTE: currently doesn't work with encodings that include NUL bytes in
//   their representation of characters other than NUL.

bool LangIDMultiTrie::loadWords(const char *filename, uint32_t langID, bool verbose)
{
   if (!filename || !*filename)
      return false ;
   Fr::CInputFile fp(filename) ;
   bool warned = false;
   unsigned linenumber = 0 ;
   if (fp)
      {
      unsigned wc = 0 ;
      while (Fr::CharPtr line = fp.getTrimmedLine())
	 {
	 linenumber++ ;
	 char *lineptr = (char*)line ;
	 // check if blank or comment line
	 if (!*lineptr || *lineptr == ';' || *lineptr == '#')
	    continue ;
	 // extract the frequency
	 char *freq_end ;
	 uint32_t freq = (uint32_t)strtoul(lineptr,&freq_end,0) ;
	 if (freq == 0 || freq_end == lineptr)
	    {
	    if (!warned)
	       {
	       cerr << "Invalid text on line " << linenumber << " of file '"
		    << filename << "'" << endl ;
	       warned = true ;
	       }
	    continue ;
	    }
	 // trim leading and trailing whitespace from rest of line
	 lineptr = freq_end ;
	 while (*lineptr == ' ' || *lineptr == '\t')
	    lineptr++ ;
	 if (!*lineptr)
	    continue ;
	 char *lineend = strchr(lineptr,'\0') ;
	 while (lineend > lineptr && 
		(lineend[-1] <= ' ' || lineend[-1] == '\t'))
	    *--lineend = '\0' ;
	 unsigned len = strlen(lineptr) ;
	 insert((uint8_t*)lineptr,len,langID,freq,false) ;
	 wc++ ;
	 }
      if (verbose)
	 cerr << "Read " << wc << " words from '" << filename << "'" << endl ;
      return true ;
      }
   else
      {
      cerr << "Unable to read word list from '" << filename << "'" << endl ;
      return false ;
      }
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::allocateNode()
{
   if (m_used >= m_capacity)
      {
      // we've filled the node array, so add another bucket
      unsigned num_buckets = m_capacity / BUCKET_SIZE + 1 ;
      auto newbuckets = Fr::NewR<MultiTrieNode*>(m_nodes,num_buckets) ;
      if (newbuckets)
	 {
	 m_nodes = newbuckets ;
	 m_nodes[num_buckets-1] = Fr::New<MultiTrieNode>(BUCKET_SIZE) ;
	 if (m_nodes[num_buckets-1] == 0)
	    {
	    fprintf(stderr,"Out of memory!\n") ;
	    abort() ;
	    }
	 m_capacity += BUCKET_SIZE ;
	 }
      else
	 {
	 fprintf(stderr,"Out of memory!\n") ;
	 abort() ;
	 }
      }
   uint32_t node_index = m_used++ ;
   auto n = node(node_index) ;
   new (n) MultiTrieNode ;
   return node_index ;
}

//----------------------------------------------------------------------

MultiTrieNode* LangIDMultiTrie::node(uint32_t N) const
{
   if (N < m_used)
      {
      auto bucket = m_nodes[N / BUCKET_SIZE] ;
      return (bucket) ? &bucket[N % BUCKET_SIZE] : 0 ;
      }
   else
      return nullptr ;
}

//----------------------------------------------------------------------

MultiTrieNode *LangIDMultiTrie::rootNode() const
{
   return node(LangIDMultiTrie::ROOT_INDEX) ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::insertNybble(uint32_t nodeindex, uint8_t nybble)
{
   auto n = node(nodeindex) ;
   if (n->childPresent(nybble))
      return n->childIndex(nybble) ;
   if (n->insertChild(nybble,this))
      return n->childIndex(nybble) ;
   return (uint32_t)~0 ;
}

//----------------------------------------------------------------------

void LangIDMultiTrie::insertChild(uint32_t &nodeindex, uint8_t keybyte)
{
   if (ignoringWhiteSpace() && keybyte == ' ')
      return ;
#if MTRIE_BITS_PER_LEVEL == 8
   nodeindex = insertNybble(nodeindex,keybyte) ;
#elif MTRIE_BITS_PER_LEVEL == 4
   nodeindex = insertNybble(nodeindex,keybyte >> 4) ;
   nodeindex = insertNybble(nodeindex,keybyte & 0x0F) ;
#elif MTRIE_BITS_PER_LEVEL == 3
   nodeindex = insertNybble(nodeindex,(keybyte >> 6) & 0x03) ;
   nodeindex = insertNybble(nodeindex,(keybyte >> 3) & 0x07) ;
   nodeindex = insertNybble(nodeindex,keybyte & 0x07) ;
#elif MTRIE_BITS_PER_LEVEL == 2
   nodeindex = insertNybble(nodeindex,(keybyte >> 6) & 0x03) ;
   nodeindex = insertNybble(nodeindex,(keybyte >> 4) & 0x03) ;
   nodeindex = insertNybble(nodeindex,(keybyte >> 2) & 0x03) ;
   nodeindex = insertNybble(nodeindex,keybyte & 0x03) ;
#else
#  error No code for given MTRIE_BITS_PER_LEVEL
#endif      
   return ;
}


//----------------------------------------------------------------------

bool LangIDMultiTrie::insert(const uint8_t *key, unsigned keylength,
		       uint32_t langID, uint32_t frequency, bool stopgram)
{
   if (keylength > m_maxkeylen)
      m_maxkeylen = keylength ;
   uint32_t cur_index = ROOT_INDEX ;
   while (keylength > 0)
      {
      this->insertChild(cur_index,*key) ;
      key++ ;
      keylength-- ;
      }
   auto leaf = node(cur_index) ;
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

uint32_t LangIDMultiTrie::find(const uint8_t *key, unsigned keylength) const
{
   uint32_t cur_index = ROOT_INDEX ;
   while (keylength > 0)
      {
      if (!extendKey(cur_index,*key))
	 return 0 ;
      key++ ;
      keylength-- ;
      }
   auto n = node(cur_index) ;
   return n ? n->frequency() : 0 ;
}

//----------------------------------------------------------------------

MultiTrieNode* LangIDMultiTrie::findNode(const uint8_t *key, unsigned keylength) const
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

bool LangIDMultiTrie::extendNybble(uint32_t &nodeindex, uint8_t nybble) const
{
   auto n = node(nodeindex) ;
   if (n->childPresent(nybble))
      {
      nodeindex = n->childIndex(nybble) ;
      return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::extendKey(uint32_t &nodeindex, uint8_t keybyte) const
{
   if (ignoringWhiteSpace() && keybyte == ' ')
      return true ;
   uint32_t idx = nodeindex ;
#if MTRIE_BITS_PER_LEVEL == 8
   if (extendNybble(idx,keybyte))
#elif MTRIE_BITS_PER_LEVEL == 4
   if (extendNybble(idx,keybyte >> 4) &&
       extendNybble(idx,keybyte & 0x0F))
#elif MTRIE_BITS_PER_LEVEL == 3
   if (extendNybble(idx,(keybyte >> 6) & 0x03) &&
       extendNybble(idx,(keybyte >> 3) & 0x07) &&
       extendNybble(idx,keybyte & 0x07))
#elif MTRIE_BITS_PER_LEVEL == 2
   if (extendNybble(idx,(keybyte >> 6) & 0x03) &&
       extendNybble(idx,(keybyte >> 4) & 0x03) &&
       extendNybble(idx,(keybyte >> 2) & 0x03) &&
       extendNybble(idx,keybyte & 0x03))
#else
#  error No code for given MTRIE_BITS_PER_LEVEL
#endif
      {
      nodeindex = idx ;
      return true ;
      }
   nodeindex = 0 ;
   return false ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerate(uint8_t *keybuf, unsigned maxkeylength, EnumFn *fn, void *user_data) const
{
   if (keybuf && fn && m_nodes[ROOT_INDEX])
      {
      memset(keybuf,'\0',maxkeylength) ;
      return enumerateChildren(ROOT_INDEX,keybuf,maxkeylength*8,0,fn,user_data) ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerateChildren(uint32_t nodeindex,
				      uint8_t *keybuf,
				      unsigned max_keylength_bits,
				      unsigned curr_keylength_bits,
   				      LangIDMultiTrie::EnumFn *fn,
				      void *user_data) const
{
   auto n = node(nodeindex) ;
   assert(!(!n->leaf() && n->frequencies())) ;
   if (n->leaf() && !fn(n,keybuf,curr_keylength_bits/8,user_data))
      return false ;
   if (curr_keylength_bits < max_keylength_bits)
      {
      unsigned curr_bits = curr_keylength_bits + MTRIE_BITS_PER_LEVEL ;
      for (size_t i = 0 ; i < (1<<MTRIE_BITS_PER_LEVEL) ; i++)
	 {
	 uint32_t child = n->childIndex(i) ;
	 if (child != LangIDMultiTrie::NULL_INDEX)
	    {
	    auto childnode = node(child) ;
	    if (childnode)
	       {
	       unsigned byte = curr_keylength_bits / 8 ;
	       unsigned shift = LEVEL_SIZE - (curr_keylength_bits%8) - MTRIE_BITS_PER_LEVEL ;
#if MTRIE_BITS_PER_LEVEL == 3
	       if (shift == 0) curr_bits = curr_keylength_bits + MTRIE_BITS_PER_LEVEL - 1 ;
#endif
	       unsigned mask = (((1<<MTRIE_BITS_PER_LEVEL)-1) << shift) ;
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

bool LangIDMultiTrie::enumerateTerminalNodes(uint32_t nodeindex, uint32_t &count, unsigned keylen_bits) const
{
   auto n = node(nodeindex)   ;
   if (!n->hasChildren())
      return true ;
   else if ((keylen_bits % 8) == 0 && allChildrenAreTerminals(nodeindex))
      {
      count += numExtensions(nodeindex) ;
      return true ;
      }
   keylen_bits = keylen_bits + BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
      {
      uint32_t child = n->childIndex(i) ;
      if (child != LangIDMultiTrie::NULL_INDEX && !enumerateTerminalNodes(child,count,keylen_bits))
	 return false ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerateFullByteNodes(uint32_t nodeindex, uint32_t &count, unsigned keylen_bits) const
{
   if (keylen_bits % 8 == 0)
      count++ ;
   keylen_bits = keylen_bits + MTRIE_BITS_PER_LEVEL ;
#if MTRIE_BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   auto n = node(nodeindex) ;
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
      {
      uint32_t child = n->childIndex(i) ;
      if (child != LangIDMultiTrie::NULL_INDEX)
	 {
	 if (!enumerateFullByteNodes(child,count,keylen_bits))
	    return false ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

unsigned LangIDMultiTrie::numExtensions(uint32_t nodeindex, unsigned bits) const
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

bool LangIDMultiTrie::allChildrenAreTerminals(uint32_t nodeindex, unsigned bits) const
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

bool LangIDMultiTrie::singleChild(uint32_t nodeindex) const 
{
   auto node = this->node(nodeindex) ;
   for (size_t i = 0 ; i < 8 && node ; i += MTRIE_BITS_PER_LEVEL)
      {
      unsigned index = ~0 ;
      for (unsigned ch = 0 ; ch < (1<<MTRIE_BITS_PER_LEVEL) ; ch++)
	 {
	 if (node->childPresent(ch))
	    {
	    if (index != ~0U)
	       return false ; 		// multiple children
	    index = ch ;
	    }
	 }
      if (index == ~0U)
	 return false ; 		// no children at all
      node = this->node(node->childIndex(index)) ;
      }
   return node != nullptr ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::singleChildSameFreq(uint32_t nodeindex, bool allow_nonleaf, double ratio) const 
{
   auto node = this->node(nodeindex) ;
   auto parent_freq = node->frequency() ;
   for (size_t i = 0 ; i < 8 && node ; i += MTRIE_BITS_PER_LEVEL)
      {
      unsigned index = ~0 ;
      for (unsigned ch = 0 ; ch < (1<<MTRIE_BITS_PER_LEVEL) ; ch++)
	 {
	 if (node->childPresent(ch))
	    {
	    if (index != ~0U)
	       return false ; 		// multiple children
	    index = ch ;
	    }
	 }
      if (index == ~0U)
	 return false ; 		// no children at all
      node = this->node(node->childIndex(index)) ;
      }
   if (!node)
      return false ;
   auto freq = node->frequency() ;
   return ((freq <= parent_freq && freq >= ratio * parent_freq) || (allow_nonleaf && freq == 0)) ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::numFullByteNodes() const
{
   uint32_t count = 0 ;
   return enumerateFullByteNodes(ROOT_INDEX,count) ? count : 0 ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::numTerminalNodes() const
{
   uint32_t count = 0 ;
   return enumerateTerminalNodes(ROOT_INDEX,count) ? count : 0 ;
}

//----------------------------------------------------------------------

static bool count_freq_records(const MultiTrieNode *node, const uint8_t *,
			       unsigned, void *user_data)
{
   uint32_t *count = (uint32_t*)user_data ;
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

LangIDMultiTrie *LangIDMultiTrie::load(Fr::CFile& f)
{
   if (!f)
      return nullptr ;
   const size_t siglen = sizeof(MULTITRIE_SIGNATURE) ;
   char signature[siglen] ;
   if (f.read(signature,siglen,sizeof(char)) != siglen ||
      memcmp(signature,MULTITRIE_SIGNATURE,siglen) != 0)
      {
      // error: wrong file type
      return nullptr ;
      }
   unsigned char version ;
   if (!f.readValue(&version) || version != MULTITRIE_FORMAT_VERSION)
      {
      // error: wrong version of data file
      return nullptr ;
      }
   unsigned char bits ;
   if (!f.readValue(&bits) || bits != MTRIE_BITS_PER_LEVEL)
      {
      // error: wrong type of trie
      return nullptr ;
      }
   Fr::UInt32 val_used, val_tokens, val_keylen ;
   if (!f.readValue(&val_used) ||
      !f.readValue(&val_tokens) ||
      !f.readValue(&val_keylen))
      {
      // error reading header
      return nullptr ;
      }
   uint32_t used = val_used.load() ;
   auto trie = new LangIDMultiTrie(used) ;
   if (!trie)
      return nullptr ;
   trie->m_used = used ;
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
	 delete trie ;
	 return nullptr ;
	 }
      }
   // finally, read the frequency information
   if (!MultiTrieFrequency::readAll(f))
      {
      delete trie ;
      return nullptr ;
      }
   return trie ;
}

//----------------------------------------------------------------------

LangIDMultiTrie *LangIDMultiTrie::load(const char *filename)
{
   Fr::CInputFile fp(filename) ;
   return load(fp) ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::writeHeader(Fr::CFile& f) const
{
   // write the signature string
   const size_t siglen = sizeof(MULTITRIE_SIGNATURE) ;
   if (f.write(MULTITRIE_SIGNATURE,siglen,1) != siglen)
      return false; 
   // follow with the format version number
   unsigned char version = MULTITRIE_FORMAT_VERSION ;
   unsigned char bits = MTRIE_BITS_PER_LEVEL ;
   if (!f.writeValue(version) ||
      !f.writeValue(bits))
      return false ;
   // write out the size of the trie
   Fr::UInt32 val_used(size()), val_tokens(totalTokens()), val_keylen(longestKey()) ;
   if (!f.writeValue(val_used) ||
      !f.writeValue(val_tokens) ||
      !f.writeValue(val_keylen))
      return false ;
   // pad the header with NULs for the unused reserved portion of the header
   return f.putNulls(MULTITRIE_PADBYTES_1) ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::write(Fr::CFile& f) const
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
   Fr::COutputFile fp(filename,Fr::CFile::safe_rewrite) ;
   return this->write(fp) ? fp.close() : false ;
}

//----------------------------------------------------------------------

static const char hexdigit[] = "0123456789ABCDEF" ;

void write_escaped_char(uint8_t c, Fr::CFile& f)
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

void write_escaped_key(Fr::CFile& f, const uint8_t* key, unsigned keylen)
{
   for (size_t i = 0 ; i < keylen ; i++)
      {
      write_escaped_char(key[i],f) ;
      }
   return ;
}

//----------------------------------------------------------------------

static bool dump_ngram(const MultiTrieNode *node, const uint8_t *key,
		       unsigned keylen, void *user_data)
{
   Fr::CFile& f = *((Fr::CFile*)user_data) ;
   if (f && node)
      {
      f.printf("   ") ;
      write_escaped_key(f,key,keylen) ;
      f.printf("  ::") ;
      auto freq = node->frequencies() ;
      for ( ; freq ; freq = freq->next())
	 {
	 f.printf(" %lu=%lu",(unsigned long)freq->languageID(),
		 (unsigned long)freq->frequency()) ;
	 }
      f.printf("\n") ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::dump(Fr::CFile& f) const
{
   Fr::LocalAlloc<uint8_t,10000> keybuf(longestKey()) ;
   return keybuf ? enumerate(keybuf,longestKey(),dump_ngram,&f) : false ;
}

// end of file mtrie.C //
