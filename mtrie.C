/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: mtrie.C - bit-slice-based Word-frequency multi-trie		*/
/*  Version:  1.30				       			*/
/*  LastEdit: 27jun2019							*/
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
#define ROOT_INDEX 0

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

MultiTrieFrequency *MultiTrieFrequency::s_base_address = 0 ;
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

MultiTrieFrequency *MultiTrieFrequency::allocate(uint32_t freq, uint32_t lang,
						 bool stopgram)
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
	 return 0 ;
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
      MultiTrieFrequency *f = allocate(incr,ID) ;
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

MultiTrieFrequency *MultiTrieFrequency::read(FILE *fp)
{
   if (fp)
      {
      char buffer[sizeof(MultiTrieFrequency)] ;
      if (fread(buffer,sizeof(char),sizeof(buffer),fp) == sizeof(buffer))
	 {
	 MultiTrieFrequency *freq = (MultiTrieFrequency*)buffer ;
	 return new MultiTrieFrequency(*freq) ;
	 }
      }
   return 0 ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::readAll(FILE *fp)
{
   if (fp)
      {
      Fr::UInt32 cnt ;
      if (fread((char*)&cnt,sizeof(char),sizeof(cnt),fp) == sizeof(cnt))
	 {
	 uint32_t count = cnt.load() ;
	 MultiTrieFrequency *base = Fr::New<MultiTrieFrequency>(count) ;
	 if (base)
	    {
	    if (fread(base,sizeof(MultiTrieFrequency),count,fp) == count)
	       {
	       Fr::Free(s_base_address) ;
	       s_base_address = base ;
	       s_max_alloc = s_curr_alloc = count ;
	       return true ;
	       }
	    else
	       Fr::Free(base) ;
	    }
	 }
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::write(FILE *fp) const
{
   if (fp)
      {
      if (fwrite(this,sizeof(char),sizeof(*this),fp) != sizeof(*this))
	 return false ;
      return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::writeAll(FILE *fp)
{
   if (fp)
      {
      Fr::UInt32 count(s_curr_alloc) ;
      if (fwrite((char*)&count,sizeof(char),sizeof(count),fp) != sizeof(count))
	 return false ;
      return (fwrite(baseAddress(),sizeof(MultiTrieFrequency),s_curr_alloc,fp)
	      == s_curr_alloc) ;
      }
   return false ;
}

/************************************************************************/
/*	Methods for class MultiTrieNode					*/
/************************************************************************/

MultiTrieNode::MultiTrieNode()
{
   m_frequency_info = INVALID_FREQ ;
   m_isleaf = false ;
   for (unsigned i = 0 ; i < lengthof(m_children) ; i++)
      m_children[i] = LangIDMultiTrie::NULL_INDEX ;
   return ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::isStopgram(unsigned langID) const
{
   const MultiTrieFrequency *freq = frequencies() ;
   for ( ; freq ; freq = freq->next())
      {
      if (freq->languageID() == langID)
	 return freq->isStopgram() ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::hasChildren() const
{
   for (size_t i = 0 ; i < lengthof(m_children) ; i++)
      {
      if (m_children[i])
	 return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::childPresent(unsigned int N) const 
{
   return (N < lengthof(m_children)) ? (m_children[N] != LangIDMultiTrie::NULL_INDEX) : false ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::singleChild(const LangIDMultiTrie *trie) const 
{
   const MultiTrieNode *node = this ;
   for (size_t i = 0 ; i < 8 && node ; i += MTRIE_BITS_PER_LEVEL)
      {
      unsigned index = lengthof(m_children) ;
      for (unsigned ch = 0 ; ch < lengthof(m_children) ; ch++)
	 {
	 if (node->m_children[ch] != LangIDMultiTrie::NULL_INDEX)
	    {
	    if (index != lengthof(m_children))
	       return false ; 		// multiple children
	    index = ch ;
	    }
	 }
      if (index == lengthof(m_children))
	 return false ; 		// no children at all
      node = trie->node(node->m_children[index]) ;
      }
   return node != 0 ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::singleChildSameFreq(const LangIDMultiTrie *trie,
				        bool allow_nonleaf,
					double ratio) const 
{
   const MultiTrieNode *node = this ;
   for (size_t i = 0 ; i < 8 && node ; i += MTRIE_BITS_PER_LEVEL)
      {
      unsigned index = lengthof(m_children) ;
      for (unsigned ch = 0 ; ch < lengthof(m_children) ; ch++)
	 {
	 if (node->m_children[ch] != LangIDMultiTrie::NULL_INDEX)
	    {
	    if (index != lengthof(m_children))
	       return false ; 		// multiple children
	    index = ch ;
	    }
	 }
      if (index == lengthof(m_children))
	 return false ; 		// no children at all
      node = trie->node(node->m_children[index]) ;
      }
   if (!node)
      return false ;
   uint32_t freq = node->frequency() ;
   return ((freq <= frequency() && freq >= ratio * frequency()) ||
	   (allow_nonleaf && freq == 0)) ;
}

//----------------------------------------------------------------------

uint32_t MultiTrieNode::childIndex(unsigned int N) const
{
   return (N < lengthof(m_children)) ? m_children[N] : LangIDMultiTrie::NULL_INDEX ;
}

//----------------------------------------------------------------------

unsigned MultiTrieNode::numExtensions(const LangIDMultiTrie *trie) const
{
   unsigned count = 0 ;
   for (size_t i1 = 0 ; i1 < lengthof(m_children) ; i1++)
      {
      if (m_children[i1] == LangIDMultiTrie::NULL_INDEX) continue ;
      MultiTrieNode *child1 = trie->node(m_children[i1]) ;
#if MTRIE_BITS_PER_LEVEL < 8
      if (!child1)
	 continue ;
      for (size_t i2 = 0 ; i2 < lengthof(m_children) ; i2++)
	 {
	 if (child1->m_children[i2] == LangIDMultiTrie::NULL_INDEX) continue ;
	 MultiTrieNode *child2 = trie->node(child1->m_children[i2]) ;
#if MTRIE_BITS_PER_LEVEL < 4
	 if (!child2)
	    continue ;
	 for (size_t i3 = 0 ; i3 < lengthof(m_children) ; i3++)
	    {
	    if (child2->m_children[i3] == LangIDMultiTrie::NULL_INDEX) continue ;
	    MultiTrieNode *child3 = trie->node(child2->m_children[i3]) ;
#if MTRIE_BITS_PER_LEVEL == 2
	    if (!child3)
	       continue ;
	    for (size_t i4 = 0 ; i4 < lengthof(m_children) ; i4++)
	       {
	       if (child3->m_children[i4] == LangIDMultiTrie::NULL_INDEX) continue ;
	       MultiTrieNode *child4 = trie->node(child3->m_children[i4]) ;
	       if (child4)
		  count++ ;
	       }
#else // MTRIE_BITS_PER_LEVEL > 2
	    if (child3)
	       count++ ;
#endif
	    }
#else // MTRIE_BITS_PER_LEVEL >= 4
	 if (child2)
	    count++ ;
#endif
	 }
#else // MTRIE_BITS_PER_LEVEL == 8
      if (child1)
	 count++ ;
#endif
      }
   return count ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::allChildrenAreTerminals(const LangIDMultiTrie *trie) const
{
   for (size_t i1 = 0 ; i1 < lengthof(m_children) ; i1++)
      {
      if (m_children[i1] == LangIDMultiTrie::LangIDMultiTrie::NULL_INDEX) continue ;
      MultiTrieNode *child1 = trie->node(m_children[i1]) ;
#if MTRIE_BITS_PER_LEVEL < 8
      if (!child1)
	 continue ;
      for (size_t i2 = 0 ; i2 < lengthof(m_children) ; i2++)
	 {
	 if (child1->m_children[i2] == LangIDMultiTrie::NULL_INDEX) continue ;
	 MultiTrieNode *child2 = trie->node(child1->m_children[i2]) ;
#if MTRIE_BITS_PER_LEVEL < 4
	 if (!child2)
	    continue ;
	 for (size_t i3 = 0 ; i3 < lengthof(m_children) ; i3++)
	    {
	    if (child2->m_children[i3] == LangIDMultiTrie::NULL_INDEX) continue ;
	    MultiTrieNode *child3 = trie->node(child2->m_children[i3]) ;
#if MTRIE_BITS_PER_LEVEL == 2
	    if (!child3)
	       continue ;
	    for (size_t i4 = 0 ; i4 < lengthof(m_children) ; i4++)
	       {
	       if (child3->m_children[i4] == LangIDMultiTrie::NULL_INDEX) continue ;
	       MultiTrieNode *child4 = trie->node(child3->m_children[i4]) ;
	       if (child4 && child4->hasChildren())
		  return false ;
	       }
#else // MTRIE_BITS_PER_LEVEL > 2
	    if (child3 && child3->hasChildren())
	       return false ;
#endif
	    }
#else // MTRIE_BITS_PER_LEVEL >= 4
	 if (child2 && child2->hasChildren())
	    return false ;
#endif
	 }
#else // MTRIE_BITS_PER_LEVEL == 8
      if (child1 && child1->hasChildren())
	 return false ;
#endif
      }
   return true ;
}

//----------------------------------------------------------------------

unsigned MultiTrieNode::numFrequencies() const
{
   const MultiTrieFrequency *freq = frequencies() ;
   unsigned count = 0 ;
   for ( ; freq ; freq = freq->next())
      count++ ;
   return count ;
}

//----------------------------------------------------------------------

uint32_t MultiTrieNode::frequency(uint32_t langID) const
{
   const MultiTrieFrequency *freq = frequencies() ;
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
   MultiTrieFrequency *f = MultiTrieFrequency::getAddress(m_frequency_info) ;
   if (f)
      f->setFrequency(langID,freq,stopgram) ;
   else
      {
      // insert the initial frequency record
      f = MultiTrieFrequency::allocate(freq,langID,stopgram) ;
      m_frequency_info = MultiTrieFrequency::getIndex(f) ;
      }
   return ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::setFrequencies(MultiTrieFrequency *freqs)
{
   if (freqs < MultiTrieFrequency::baseAddress())
      return false ;
   uint32_t idx = MultiTrieFrequency::getIndex(freqs) ;
   if (idx == INVALID_FREQ)
      return false ;
   m_frequency_info = idx ;
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

bool MultiTrieNode::enumerateChildren(const LangIDMultiTrie *trie,
				      uint8_t *keybuf,
				      unsigned max_keylength_bits,
				      unsigned curr_keylength_bits,
				      MultiTrieEnumFn *fn,
				      void *user_data) const
{
  assert(!(!leaf() && frequencies())) ;
   if (leaf() && !fn(this,keybuf,curr_keylength_bits/8,user_data))
      return false ;
   if (curr_keylength_bits < max_keylength_bits)
      {
      unsigned curr_bits = curr_keylength_bits + MTRIE_BITS_PER_LEVEL ;
      for (size_t i = 0 ; i < lengthof(m_children) ; i++)
	 {
	 uint32_t child = childIndex(i) ;
	 if (child != LangIDMultiTrie::NULL_INDEX)
	    {
	    MultiTrieNode *childnode = trie->node(child) ;
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
	       if (!childnode->enumerateChildren(trie,keybuf,
						 max_keylength_bits,
						 curr_bits,fn,user_data))
		  return false ;
	       }
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::enumerateFullByteNodes(const LangIDMultiTrie *trie,
					   unsigned keylen_bits,
					   uint32_t &count) const
{
   if (keylen_bits % 8 == 0)
      count++ ;
   keylen_bits = keylen_bits + MTRIE_BITS_PER_LEVEL ;
#if MTRIE_BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   for (size_t i = 0 ; i < lengthof(m_children) ; i++)
      {
      uint32_t child = childIndex(i) ;
      if (child != LangIDMultiTrie::NULL_INDEX)
	 {
	 MultiTrieNode *childnode = trie->node(child) ;
	 if (childnode)
	    {
	    if (!childnode->enumerateFullByteNodes(trie,keylen_bits,count))
	       return false ;
	    }
	 else
	    return false ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::enumerateTerminalNodes(const LangIDMultiTrie *trie,
					   unsigned keylen_bits,
					   uint32_t &count) const
{
   if (!hasChildren())
      return true ;
   else if ((keylen_bits % 8) == 0 && this->allChildrenAreTerminals(trie))
      {
      count += this->numExtensions(trie) ;
      return true ;
      }
   keylen_bits = keylen_bits + MTRIE_BITS_PER_LEVEL ;
#if MTRIE_BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   for (size_t i = 0 ; i < lengthof(m_children) ; i++)
      {
      uint32_t child = childIndex(i) ;
      if (child != LangIDMultiTrie::NULL_INDEX)
	 {
	 MultiTrieNode *childnode = trie->node(child) ;
	 if (!childnode ||
	     !childnode->enumerateTerminalNodes(trie,keylen_bits,count))
	    return false ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::load(FILE *fp)
{
   if (fp)
      {
      if (fread(this,sizeof(*this),1,fp) == 1)
	 return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::write(FILE *fp) const
{
   if (fp && fwrite(this,sizeof(*this),1,fp) == 1)
      return true ;
   return false ;
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
   m_nodes = 0 ;
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
      new (node(ROOT_INDEX)) MultiTrieNode ;
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
   FILE *fp = fopen(filename,"rb") ;
   bool warned = false;
   unsigned linenumber = 0 ;
   if (fp)
      {
      unsigned wc = 0 ;
      char line[16384] ;
      memset(line,'\0',sizeof(line)) ;
      while (!feof(fp))
	 {
	 if (!fgets(line,sizeof(line),fp))
	    break ;
	 linenumber++ ;
	 char *lineptr = line ;
	 // start by trimming leading whitespace
	 while (*lineptr == ' ' || *lineptr == '\t')
	    lineptr++ ;
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
      fclose(fp) ;
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
      MultiTrieNode **newbuckets = Fr::NewR<MultiTrieNode*>(m_nodes,num_buckets) ;
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
   MultiTrieNode *n = node(node_index) ;
   new (n) MultiTrieNode ;
   return node_index ;
}

//----------------------------------------------------------------------

MultiTrieNode *LangIDMultiTrie::node(uint32_t N) const
{
   if (N < m_used)
      {
      MultiTrieNode *bucket = m_nodes[N / BUCKET_SIZE] ;
      return (bucket) ? &bucket[N % BUCKET_SIZE] : 0 ;
      }
   else
      return 0 ;
}

//----------------------------------------------------------------------

MultiTrieNode *LangIDMultiTrie::rootNode() const
{
   return node(ROOT_INDEX) ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::insertNybble(uint32_t nodeindex, uint8_t nybble)
{
   MultiTrieNode *n = node(nodeindex) ;
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
   MultiTrieNode *leaf = node(cur_index) ;
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
   MultiTrieNode *n = node(cur_index) ;
   return n ? n->frequency() : 0 ;
}

//----------------------------------------------------------------------

MultiTrieNode *LangIDMultiTrie::findNode(const uint8_t *key, unsigned keylength)
   const
{
   uint32_t cur_index = ROOT_INDEX ;
   while (keylength > 0)
      {
      if (!extendKey(cur_index,*key))
	 return 0 ;
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
   MultiTrieNode *n = node(cur_index) ;
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
      MultiTrieNode *n = node(cur_index) ;
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
   MultiTrieNode *n = node(nodeindex) ;
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

bool LangIDMultiTrie::enumerate(uint8_t *keybuf, unsigned maxkeylength,
			  MultiTrieEnumFn *fn, void *user_data) const
{
   if (keybuf && fn && m_nodes[0])
      {
      memset(keybuf,'\0',maxkeylength) ;
      return m_nodes[0]->enumerateChildren(this,keybuf,maxkeylength*8,0,fn,
					   user_data) ;
      }
   return false ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::numFullByteNodes() const
{
   uint32_t count = 0 ;
   return (m_nodes[0]->enumerateFullByteNodes(this,0,count) ? count : 0) ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::numTerminalNodes() const
{
   uint32_t count = 0 ;
   return (m_nodes[0]->enumerateTerminalNodes(this,0,count) ? count : 0) ;
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

LangIDMultiTrie *LangIDMultiTrie::load(FILE *fp)
{
   if (fp)
      {
      const size_t siglen = sizeof(MULTITRIE_SIGNATURE) ;
      char signature[siglen] ;
      if (fread(signature,sizeof(char),siglen,fp) != siglen ||
	  memcmp(signature,MULTITRIE_SIGNATURE,siglen) != 0)
	 {
	 // error: wrong file type
	 return 0 ;
	 }
      unsigned char version ;
      if (fread(&version,sizeof(char),sizeof(version),fp) != sizeof(version)
	  || version != MULTITRIE_FORMAT_VERSION)
	 {
	 // error: wrong version of data file
	 return 0 ;
	 }
      unsigned char bits ;
      if (fread(&bits,sizeof(char),sizeof(bits),fp) != sizeof(bits) ||
	  bits != MTRIE_BITS_PER_LEVEL)
	 {
	 // error: wrong type of trie
	 return 0 ;
	 }
      Fr::UInt32 val_used, val_tokens, val_keylen ;
      if (fread((char*)&val_used,sizeof(val_used),1,fp) != 1 ||
	 fread((char*)&val_tokens,sizeof(val_tokens),1,fp) != 1 ||
	 fread((char*)&val_keylen,sizeof(val_keylen),1,fp) != 1)
	 {
	 // error reading header
	 return 0 ;
	 }
      uint32_t used = val_used.load() ;
      auto trie = new LangIDMultiTrie(used) ;
      if (!trie)
	 return 0 ;
      trie->m_used = used ;
      trie->addTokenCount(val_tokens.load()) ;
      trie->m_maxkeylen = (unsigned)val_keylen.load() ;
      // skip the padding (reserved bytes)
      fseek(fp,MULTITRIE_PADBYTES_1,SEEK_CUR) ;
      // read the actual trie nodes
      for (size_t i = 0 ; i < used ; i++)
	 {
	 MultiTrieNode *node = trie->node(i) ;
	 if (!node->load(fp))
	    {
	    delete trie ;
	    return 0 ;
	    }
	 }
      // finally, read the frequency information
      if (!MultiTrieFrequency::readAll(fp))
	 {
	 delete trie ;
	 return 0 ;
	 }
//FIXME
      return trie ;
      }
   return 0 ;
}

//----------------------------------------------------------------------

LangIDMultiTrie *LangIDMultiTrie::load(const char *filename)
{
   Fr::CInputFile fp(filename) ;
   return fp ? load(fp.fp()) : nullptr ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::writeHeader(FILE *fp) const
{
   // write the signature string
   const size_t siglen = sizeof(MULTITRIE_SIGNATURE) ;
   if (fwrite(MULTITRIE_SIGNATURE,sizeof(char),siglen,fp) != siglen)
      return false; 
   // follow with the format version number
   unsigned char version = MULTITRIE_FORMAT_VERSION ;
   if (fwrite(&version,sizeof(char),sizeof(version),fp) != sizeof(version))
      return false ;
   unsigned char bits = MTRIE_BITS_PER_LEVEL ;
   if (fwrite(&bits,sizeof(char),sizeof(bits),fp) != sizeof(bits))
      return false ;
   // write out the size of the trie
   Fr::UInt32 val_used(size()), val_tokens(totalTokens()), val_keylen(longestKey()) ;
   if (fwrite((char*)&val_used,sizeof(val_used),1,fp) != 1 || 
      fwrite((char*)&val_tokens,sizeof(val_tokens),1,fp) != 1 ||
      fwrite((char*)&val_keylen,sizeof(val_keylen),1,fp) != 1)
      return false ;
   // pad the header with NULs for the unused reserved portion of the header
   for (size_t i = 0 ; i < MULTITRIE_PADBYTES_1 ; i++)
      {
      if (fputc('\0',fp) == EOF)
	 return false ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::write(FILE *fp) const
{
   if (fp)
      {
      if (writeHeader(fp))
	 {
	 // write the actual trie nodes
	 for (size_t i = 0 ; i < size() ; i++)
	    {
	    const MultiTrieNode *n = node(i) ;
	    if (!n->write(fp))
	       return false ;
	    }
	 // write the frequency information
	 if (!MultiTrieFrequency::writeAll(fp))
	    return false ;
//FIXME
	 return true ;
	 }
      }
   return false ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::write(const char *filename) const
{
   Fr::COutputFile fp(filename,Fr::CFile::safe_rewrite) ;
   //TODO: change LIMT::write to use CFile directly
   bool success = fp ? this->write(fp.fp()) : false ;
   return success ? fp.close() : false ;
}

//----------------------------------------------------------------------

static const char hexdigit[] = "0123456789ABCDEF" ;

static void write_escaped_char(uint8_t c, FILE *fp)
{
   switch (c)
      {
      case '\\':
	 fputs("\\\\",fp) ;
	 break ;
      case ' ':
	 fputs("\\ ",fp) ;
	 break ;
      default:
	 if (c < ' ')
	    {
	    fputc('\\',fp) ;
	    fputc(hexdigit[(c>>4)&0xF],fp) ;
	    fputc(hexdigit[c&0xF],fp) ;
	    }
	 else
	    fputc(c,fp) ;
	 break ;
      }
   return ;
}

//----------------------------------------------------------------------

static bool dump_ngram(const MultiTrieNode *node, const uint8_t *key,
		       unsigned keylen, void *user_data)
{
   FILE *fp = (FILE*)user_data ;
   if (fp && node)
      {
      fprintf(fp,"   ") ;
      for (size_t i = 0 ; i < keylen ; i++)
	 {
	 write_escaped_char(key[i],fp) ;
	 }
      fprintf(fp,"  ::") ;
      MultiTrieFrequency *freq = node->frequencies() ;
      for ( ; freq ; freq = freq->next())
	 {
	 fprintf(fp," %lu=%lu",(unsigned long)freq->languageID(),
		 (unsigned long)freq->frequency()) ;
	 }
      fprintf(fp,"\n") ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::dump(FILE *fp) const
{
   Fr::LocalAlloc<uint8_t,10000> keybuf(longestKey()) ;
   return keybuf ? enumerate(keybuf,longestKey(),dump_ngram,fp) : false ;
}

/************************************************************************/
/*	Methods for class MultiTriePointer				*/
/************************************************************************/

void MultiTriePointer::resetKey()
{
   m_nodeindex = ROOT_INDEX ;
   m_keylength = 0 ;
   m_failed = false ;
   return ;
}

//----------------------------------------------------------------------

bool MultiTriePointer::hasChildren(uint32_t node_index,
				   uint8_t nybble) const
{
   MultiTrieNode *n = m_trie->node(node_index) ;
   return n->childPresent(nybble) ;
}

//----------------------------------------------------------------------

bool MultiTriePointer::extendKey(uint8_t keybyte)
{
   if (m_failed)
      return false ;
   bool success = m_trie->extendKey(m_nodeindex,keybyte) ;
   if (success)
      m_keylength++ ;
   else
      m_failed = true ;
   return success ;
}

//----------------------------------------------------------------------

bool MultiTriePointer::lookupSuccessful() const
{
   MultiTrieNode *n = node() ;
   return n && n->leaf() ;
}

// end of file mtrie.C //
