/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: trie.C - Word-frequency trie					*/
/*  Version:  1.30				       			*/
/*  LastEdit: 27jun2019							*/
/*									*/
/*  (c) Copyright 2011,2012,2014,2015,2019 Ralf Brown/CMU		*/
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

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "trie.h"
#include "wildcard.h"
#include "framepac/config.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#if BITS_PER_LEVEL == 2
# define BUCKET_SIZE 65536	// must be power of 2
#else
# define BUCKET_SIZE 16384	// must be power of 2
#endif

#if BITS_PER_LEVEL == 3
#  define LEVEL_SIZE 9
#else
#  define LEVEL_SIZE 8
#endif

#define MAX_SCALE_RATIO 1E9
#define LOG_MAX_SCALE_RATIO 20.7232658369464d  /* log(MAX_SCALE_RATIO) */

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

/************************************************************************/
/*	Helper functions						*/
/************************************************************************/

#ifndef lengthof
#  define lengthof(x) (sizeof(x)/sizeof((x)[0]))
#endif /* lengthof */

//----------------------------------------------------------------------

#ifndef round_up
uint32_t round_up(uint32_t value, uint32_t granularity)
{
   return granularity * ((value + granularity - 1) / granularity) ;
}
#endif

//----------------------------------------------------------------------

double scaling_log_power(double power)
{
   return ::log(1.0 + fabs(power)) ;
}

//----------------------------------------------------------------------

uint32_t scaled_frequency(uint32_t raw_freq, uint64_t total_count)
{
   double percent = 100.0 * raw_freq / (double)total_count ;
   // avoid overflow by truncating excessively high percentages to the
   //   largest value representable in a uint32_t
   if (percent > ((uint32_t)~0) / (double)TRIE_SCALE_FACTOR)
      return (uint32_t)~0 ;
   uint32_t scaled = (uint32_t)(TRIE_SCALE_FACTOR * percent + 0.5) ;
   // avoid truncation to zero for very low percentages
   return (percent > 0.0 && scaled == 0) ? 1 : scaled ;
}

//----------------------------------------------------------------------

double scale_frequency(double freq, double power, double log_power)
{
   double scaled ;
   if (power < 0.0)
      {
      // avoid infinities and other problems by restricting the input
      //  to be a number large enough to produce a valid logarithm
      if (freq < DBL_MIN)
	 freq = DBL_MIN ;
      double shifted_ratio = 1.0 + (freq * (-power)) ; // avoid log(negative)
      scaled = ::log(shifted_ratio) / log_power ;
      // since we've mapped the value into [0,approx(1)] and we can only represent
      //   numbers up to slightly more than 4, scale by that much
      scaled *= 4.0 ;
      // the following should never trigger
      if (scaled > ((uint32_t)~0) / (double)TRIE_SCALE_FACTOR)
	 cerr << "DEBUG: truncating scaledpercent " << scaled << endl ;
      // ditto
      if (scaled <= 0.0)
	 {
	 cerr << "DEBUG: scaling underflow: " << freq << " -> " << scaled << endl  ;
	 scaled = DBL_MIN ;
	 }
      }
   else if (power > 0.0)
      {
      scaled = ::pow(100.0*freq,power) ;
      }
   else
      {
      scaled = 100.0 * freq ;
      }
   return scaled ;
}

//----------------------------------------------------------------------

uint32_t scaled_frequency(uint32_t raw_freq, uint64_t total_count,
			  double power, double log_power)
{
   double prop = raw_freq / (double)total_count ;
   double percent = scale_frequency(prop,power,log_power) ;
   // avoid overflow by truncating excessively high percentages to the
   //   largest value representable in a uint32_t
   if (percent > ((uint32_t)~0) / (double)TRIE_SCALE_FACTOR)
      return (uint32_t)~0 ;
   uint32_t scaled_value = (uint32_t)(TRIE_SCALE_FACTOR * percent + 0.5) ;
   // avoid truncation to zero for very low percentages
   return (percent > 0.0 && scaled_value == 0) ? 1 : scaled_value ;
}

//----------------------------------------------------------------------

double unscale_frequency(uint32_t freq, double power)
{
   double scaled = (freq / (double)TRIE_SCALE_FACTOR) ;
   double unscaled ;
   if (power < 0.0)
      {
      // divisor here must match multiplier in scale_frequency
      double prop = (scaled / 4.0) * scaling_log_power(power) ;
      unscaled = (::exp(prop) - 1.0) / (-power) ;
      return unscaled * 100.0 ;
      }
   else if (power > 0.0)
      {
      unscaled = ::pow(scaled,1.0/power) ;
      }
   else
      {
      unscaled = freq ;
      }
   return unscaled ;
}

/************************************************************************/
/*	Methods for class NybbleTrieNode				*/
/************************************************************************/

NybbleTrieNode::NybbleTrieNode()
{
   std::fill_n(m_children,lengthof(m_children),NybbleTrie::NULL_INDEX) ;
   m_frequency = 0 ;
   m_leaf = false ;
   m_stopgram = false ;
   return ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::hasChildren() const
{
   for (size_t i = 0 ; i < lengthof(m_children) ; i++)
      {
      if (m_children[i] != NybbleTrie::NULL_INDEX)
	 return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::hasChildren(const NybbleTrie *trie,
				 uint32_t min_freq) const
{
   for (size_t i = 0 ; i < lengthof(m_children) ; i++)
      {
      if (m_children[i] == NybbleTrie::NULL_INDEX)
	 continue ;
      NybbleTrieNode *child = trie->node(m_children[i]) ;
      if (child && child->frequency() >= min_freq)
	 return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::childPresent(unsigned int N) const 
{
   return (N < lengthof(m_children)) ? (m_children[N] != NybbleTrie::NULL_INDEX) : false ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::singleChild(const NybbleTrie *trie) const 
{
   const NybbleTrieNode *node = this ;
   for (size_t i = 0 ; i < 8 && node ; i += BITS_PER_LEVEL)
      {
      unsigned index = lengthof(m_children) ;
      for (unsigned ch = 0 ; ch < lengthof(m_children) ; ch++)
	 {
	 if (node->m_children[ch] != NybbleTrie::NULL_INDEX)
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

bool NybbleTrieNode::singleChildSameFreq(const NybbleTrie *trie,
					 bool allow_nonleaf, double ratio) const 
{
   const NybbleTrieNode *node = this ;
   for (size_t i = 0 ; i < 8 && node ; i += BITS_PER_LEVEL)
      {
      unsigned index = lengthof(m_children) ;
      for (unsigned ch = 0 ; ch < lengthof(m_children) ; ch++)
	 {
	 if (node->m_children[ch] != NybbleTrie::NULL_INDEX)
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

uint32_t NybbleTrieNode::childIndex(unsigned int N) const
{
   return (N < lengthof(m_children)) ? m_children[N] : NybbleTrie::NULL_INDEX ;
}

//----------------------------------------------------------------------

unsigned NybbleTrieNode::numExtensions(const NybbleTrie *trie,
				       uint32_t min_freq) const
{
   unsigned count = 0 ;
   for (size_t i1 = 0 ; i1 < lengthof(m_children) ; i1++)
      {
      if (m_children[i1] == NybbleTrie::NULL_INDEX) continue ;
      NybbleTrieNode *child1 = trie->node(m_children[i1]) ;
#if BITS_PER_LEVEL < 8
      if (!child1)
	 continue ;
      for (size_t i2 = 0 ; i2 < lengthof(m_children) ; i2++)
	 {
	 if (child1->m_children[i2] == NybbleTrie::NULL_INDEX) continue ;
	 NybbleTrieNode *child2 = trie->node(child1->m_children[i2]) ;
#if BITS_PER_LEVEL < 4
	 if (!child2)
	    continue ;
	 for (size_t i3 = 0 ; i3 < lengthof(m_children) ; i3++)
	    {
	    if (child2->m_children[i3] == NybbleTrie::NULL_INDEX) continue ;
	    NybbleTrieNode *child3 = trie->node(child2->m_children[i3]) ;
#if BITS_PER_LEVEL == 2
	    if (!child3)
	       continue ;
	    for (size_t i4 = 0 ; i4 < lengthof(m_children) ; i4++)
	       {
	       if (child3->m_children[i4] == NybbleTrie::NULL_INDEX) continue ;
	       NybbleTrieNode *child4 = trie->node(child3->m_children[i4]) ;
	       if (child4 && child4->frequency() >= min_freq)
		  count++ ;
	       }
#else // BITS_PER_LEVEL > 2
	    if (child3 && child3->frequency() >= min_freq)
	       count++ ;
#endif
	    }
#else // BITS_PER_LEVEL >= 4
	 if (child2 && child2->frequency() >= min_freq)
	    count++ ;
#endif
	 }
#else // BITS_PER_LEVEL == 8
      if (child1 && child1->frequency() >= min_freq)
	 count++ ;
#endif
      }
   return count ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::allChildrenAreTerminals(const NybbleTrie *trie) const
{
   for (size_t i1 = 0 ; i1 < lengthof(m_children) ; i1++)
      {
      if (m_children[i1] == NybbleTrie::NULL_INDEX) continue ;
      NybbleTrieNode *child1 = trie->node(m_children[i1]) ;
#if BITS_PER_LEVEL < 8
      if (!child1)
	 continue ;
      for (size_t i2 = 0 ; i2 < lengthof(m_children) ; i2++)
	 {
	 if (child1->m_children[i2] == NybbleTrie::NULL_INDEX) continue ;
	 NybbleTrieNode *child2 = trie->node(child1->m_children[i2]) ;
#if BITS_PER_LEVEL < 4
	 if (!child2)
	    continue ;
	 for (size_t i3 = 0 ; i3 < lengthof(m_children) ; i3++)
	    {
	    if (child2->m_children[i3] == NybbleTrie::NULL_INDEX) continue ;
	    NybbleTrieNode *child3 = trie->node(child2->m_children[i3]) ;
#if BITS_PER_LEVEL == 2
	    if (!child3)
	       continue ;
	    for (size_t i4 = 0 ; i4 < lengthof(m_children) ; i4++)
	       {
	       if (child3->m_children[i4] == NybbleTrie::NULL_INDEX) continue ;
	       NybbleTrieNode *child4 = trie->node(child3->m_children[i4]) ;
	       if (child4 && child4->hasChildren())
		  return false ;
	       }
#else // BITS_PER_LEVEL > 2
	    if (child3 && child3->hasChildren())
	       return false ;
#endif
	    }
#else // BITS_PER_LEVEL >= 4
	 if (child2 && child2->hasChildren())
	    return false ;
#endif
	 }
#else // BITS_PER_LEVEL == 8
      if (child1 && child1->hasChildren())
	 return false ;
#endif
      }
   return true ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::allChildrenAreTerminals(const NybbleTrie *trie,
					     uint32_t min_freq) const
{
   for (size_t i1 = 0 ; i1 < lengthof(m_children) ; i1++)
      {
      if (m_children[i1] == NybbleTrie::NULL_INDEX) continue ;
      NybbleTrieNode *child1 = trie->node(m_children[i1]) ;
#if BITS_PER_LEVEL < 8
      if (!child1)
	 continue ;
      for (size_t i2 = 0 ; i2 < lengthof(m_children) ; i2++)
	 {
	 if (child1->m_children[i2] == NybbleTrie::NULL_INDEX) continue ;
	 NybbleTrieNode *child2 = trie->node(child1->m_children[i2]) ;
#if BITS_PER_LEVEL < 4
	 if (!child2)
	    continue ;
	 for (size_t i3 = 0 ; i3 < lengthof(m_children) ; i3++)
	    {
	    if (child2->m_children[i3] == NybbleTrie::NULL_INDEX) continue ;
	    NybbleTrieNode *child3 = trie->node(child2->m_children[i3]) ;
#if BITS_PER_LEVEL == 2
	    if (!child3)
	       continue ;
	    for (size_t i4 = 0 ; i4 < lengthof(m_children) ; i4++)
	       {
	       if (child3->m_children[i4] == NybbleTrie::NULL_INDEX) continue ;
	       NybbleTrieNode *child4 = trie->node(child3->m_children[i4]) ;
	       if (child4 &&
		   (!child4->leaf() || child4->frequency() >= min_freq) &&
		   child4->numExtensions(trie,min_freq) > 0)
		  return false ;
	       }
#else // BITS_PER_LEVEL > 2
	    if (child3 && 
		(!child3->leaf() || child3->frequency() >= min_freq) &&
		child3->numExtensions(trie,min_freq) > 0)
	       return false ;
#endif
	    }
#else // BITS_PER_LEVEL >= 4
	 if (child2 && 
	     (!child2->leaf() || child2->frequency() >= min_freq) &&
	     child2->numExtensions(trie,min_freq) > 0)
	    return false ;
#endif
	 }
#else // BITS_PER_LEVEL == 8
      if (child1 && 
	  (!child1->leaf() || child1->frequency() >= min_freq) &&
	  child1->numExtensions(trie,min_freq) > 0)
	 return false ;
#endif
      }
   return true ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::insertChild(unsigned int N, NybbleTrie *trie)
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

void NybbleTrieNode::scaleFrequency(uint64_t total_count)
{
   m_frequency = scaled_frequency(m_frequency,total_count) ;
   return ;
}

//----------------------------------------------------------------------

void NybbleTrieNode::scaleFrequency(uint64_t total_count, double power,
				    double log_power)
{
   m_frequency = scaled_frequency(m_frequency,total_count,power,log_power) ;
   return ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::nextFrequencies(const NybbleTrie *trie,
				     uint32_t *frequencies,
				     uint8_t idx, unsigned bits) const
{
   if (bits >= 8)
      {
      frequencies[idx] = frequency() ;
      }
   else
      {
      unsigned curr_bits = bits + BITS_PER_LEVEL ;
      for (size_t i = 0 ; i < lengthof(m_children) ; i++)
	 {
	 uint32_t child = childIndex(i) ;
	 if (child == NybbleTrie::NULL_INDEX)
	    continue ;
	 NybbleTrieNode *childnode = trie->node(child) ;
	 if (!childnode)
	    continue ;
	 unsigned shift = LEVEL_SIZE - (bits % 8) - BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
	 if (shift == 0)
	    curr_bits = bits + BITS_PER_LEVEL - 1 ;
#endif
	 unsigned mask = (((1 << BITS_PER_LEVEL)-1) << shift) ;
	 idx &= ~mask ;
	 idx |= (i << shift) ;
	 if (!childnode->nextFrequencies(trie,frequencies,idx,curr_bits))
	    return false ;
	 }

      }
   return true ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::nextFrequencies(const NybbleTrie *trie,
				     uint32_t *frequencies) const
{
   if (!frequencies)
      return false ;
   return nextFrequencies(trie,frequencies,0,0) ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::enumerateChildren(const NybbleTrie *trie,
				       uint8_t *keybuf,
				       unsigned max_keylength_bits,
				       unsigned curr_keylength_bits,
				       NybbleTrieEnumFn *fn,
				       void *user_data) const
{
   if (leaf() && !fn(this,keybuf,curr_keylength_bits/8,user_data))
      return false ;
   if (curr_keylength_bits < max_keylength_bits)
      {
      unsigned curr_bits = curr_keylength_bits + BITS_PER_LEVEL ;
      for (size_t i = 0 ; i < lengthof(m_children) ; i++)
	 {
	 uint32_t child = childIndex(i) ;
	 if (child != NybbleTrie::NULL_INDEX)
	    {
	    NybbleTrieNode *childnode = trie->node(child) ;
	    if (childnode)
	       {
	       unsigned byte = curr_keylength_bits / 8 ;
	       unsigned shift = LEVEL_SIZE - (curr_keylength_bits%8) - BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
	       if (shift == 0)
		  {
		  curr_bits = curr_keylength_bits + BITS_PER_LEVEL - 1 ;
		  }
#endif
	       unsigned mask = (((1<<BITS_PER_LEVEL)-1) << shift) ;
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

bool NybbleTrieNode::enumerateTerminalNodes(const NybbleTrie *trie,
					    unsigned keylen_bits,
					    uint32_t &count,
					    uint32_t min_freq) const
{
   if (!hasChildren())
      return true ;
   else if ((keylen_bits % 8) == 0 &&
	    this->allChildrenAreTerminals(trie,min_freq))
      {
      count += this->numExtensions(trie,min_freq) ;
      return true ;
      }
   keylen_bits = keylen_bits + BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   for (size_t i = 0 ; i < lengthof(m_children) ; i++)
      {
      uint32_t child = childIndex(i) ;
      if (child != NybbleTrie::NULL_INDEX)
	 {
	 NybbleTrieNode *childnode = trie->node(child) ;
	 if (!childnode ||
	     !childnode->enumerateTerminalNodes(trie,keylen_bits,count,
						min_freq))
	    return false ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

bool NybbleTrieNode::enumerateFullNodes(const NybbleTrie *trie,
					uint32_t &count,
					uint32_t min_freq) const
{
   count++ ;
   if (this->allChildrenAreTerminals(trie,min_freq))
      return true ;
   for (size_t i1 = 0 ; i1 < lengthof(m_children) ; i1++)
      {
      if (m_children[i1] == NybbleTrie::NULL_INDEX) continue ;
      NybbleTrieNode *child1 = trie->node(m_children[i1]) ;
#if BITS_PER_LEVEL < 8
      if (!child1)
	 continue ;
      for (size_t i2 = 0 ; i2 < lengthof(m_children) ; i2++)
	 {
	 if (child1->m_children[i2] == NybbleTrie::NULL_INDEX) continue ;
	 NybbleTrieNode *child2 = trie->node(child1->m_children[i2]) ;
#if BITS_PER_LEVEL < 4
	 if (!child2)
	    continue ;
	 for (size_t i3 = 0 ; i3 < lengthof(m_children) ; i3++)
	    {
	    if (child2->m_children[i3] == NybbleTrie::NULL_INDEX) continue ;
	    NybbleTrieNode *child3 = trie->node(child2->m_children[i3]) ;
#if BITS_PER_LEVEL == 2
	    if (!child3)
	       continue ;
	    for (size_t i4 = 0 ; i4 < lengthof(m_children) ; i4++)
	       {
	       if (child3->m_children[i4] == NybbleTrie::NULL_INDEX) continue ;
	       NybbleTrieNode *child4 = trie->node(child3->m_children[i4]) ;
	       if (child4 && child4->frequency() >= min_freq)
		  child4->enumerateFullNodes(trie,count,min_freq) ;
	       }
#else // BITS_PER_LEVEL > 2
	    if (child3 && child3->frequency() >= min_freq)
	       child3->enumerateFullNodes(trie,count,min_freq) ;
#endif
	    }
#else // BITS_PER_LEVEL >= 4
	 if (child2 && child2->frequency() >= min_freq)
	    child2->enumerateFullNodes(trie,count,min_freq) ;
#endif
	 }
#else // BITS_PER_LEVEL == 8
      if (child1 && child1->frequency() >= min_freq)
	 child1->enumerateFullNodes(trie,count,min_freq) ;
#endif
      }
   return true ;
}

/************************************************************************/
/*	Methods for class NybbleTrie					*/
/************************************************************************/

NybbleTrie::NybbleTrie(uint32_t cap)
{
   init(cap) ;
   return ;
}

//----------------------------------------------------------------------

NybbleTrie::NybbleTrie(const char *filename, bool verbose)
{
   init(1) ;
   loadWords(filename,verbose) ;
   return ;
}

//----------------------------------------------------------------------

NybbleTrie::~NybbleTrie()
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

void NybbleTrie::init(uint32_t cap)
{
   m_userdata = nullptr ;
   m_maxkeylen = 0 ;
   m_totaltokens = 0 ;
   m_ignorewhitespace = false ;
   if (cap == 0)
      cap = 1 ;
   cap = round_up(cap,BUCKET_SIZE) ;
   unsigned num_buckets = cap / BUCKET_SIZE ;
   m_nodes = Fr::New<NybbleTrieNode*>(num_buckets) ;
   m_capacity = 0 ;
   if (m_nodes)
      {
      for (unsigned i = 0 ; i < num_buckets ; i++)
	 {
	 m_nodes[i] = Fr::New<NybbleTrieNode>(BUCKET_SIZE) ;
	 if (m_nodes[i])
	    m_capacity += BUCKET_SIZE ;
	 else
	    break ;
	 }
      m_used = 1 ;
      // initialize the root node
      new (node(NybbleTrie::ROOT_INDEX)) NybbleTrieNode ;
      }
   return ;
}

//----------------------------------------------------------------------
// NOTE: currently doesn't work with encodings that include NUL bytes in
//   their representation of characters other than NUL.

bool NybbleTrie::loadWords(const char *filename, bool verbose)
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
	 insert((uint8_t*)lineptr,len,freq,false) ;
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

uint32_t NybbleTrie::allocateNode()
{
   if (m_used >= m_capacity)
      {
      // we've filled the node array, so add another bucket
      unsigned num_buckets = m_capacity / BUCKET_SIZE + 1 ;
      NybbleTrieNode **newbuckets
	 = Fr::NewR<NybbleTrieNode*>(m_nodes,num_buckets) ;
      if (newbuckets)
	 {
	 m_nodes = newbuckets ;
	 m_nodes[num_buckets-1] = Fr::New<NybbleTrieNode>(BUCKET_SIZE) ;
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
   NybbleTrieNode *n = node(node_index) ;
   new (n) NybbleTrieNode ;
   return node_index ;
}

//----------------------------------------------------------------------

NybbleTrieNode* NybbleTrie::node(uint32_t N) const
{
   if (N < m_used)
      {
      NybbleTrieNode *bucket = m_nodes[N / BUCKET_SIZE] ;
      return (bucket) ? &bucket[N % BUCKET_SIZE] : 0 ;
      }
   else
      return nullptr ;
}

//----------------------------------------------------------------------

NybbleTrieNode *NybbleTrie::rootNode() const
{
   return node(NybbleTrie::ROOT_INDEX) ;
}

//----------------------------------------------------------------------

uint32_t NybbleTrie::insertNybble(uint32_t nodeindex, uint8_t nybble)
{
   NybbleTrieNode *n = node(nodeindex) ;
   uint32_t idx = n->childIndex(nybble) ;
   if (idx != NULL_INDEX)
      return idx ;
   if (n->insertChild(nybble,this))
      return n->childIndex(nybble) ;
   return (uint32_t)~0 ;
}

//----------------------------------------------------------------------

void NybbleTrie::insertChild(uint32_t &nodeindex, uint8_t keybyte)
{
   if (ignoringWhiteSpace() && keybyte == ' ')
      return ;
#if BITS_PER_LEVEL == 8
   nodeindex = insertNybble(nodeindex,keybyte) ;
#elif BITS_PER_LEVEL == 4
   nodeindex = insertNybble(nodeindex,(keybyte >> 4) & 0x0F) ;
   nodeindex = insertNybble(nodeindex,keybyte & 0x0F) ;
#elif BITS_PER_LEVEL == 3
   nodeindex = insertNybble(nodeindex,(keybyte >> 6) & 0x03) ;
   nodeindex = insertNybble(nodeindex,(keybyte >> 3) & 0x07) ;
   nodeindex = insertNybble(nodeindex,keybyte & 0x07) ;
#elif BITS_PER_LEVEL == 2
   nodeindex = insertNybble(nodeindex,(keybyte >> 6) & 0x03) ;
   nodeindex = insertNybble(nodeindex,(keybyte >> 4) & 0x03) ;
   nodeindex = insertNybble(nodeindex,(keybyte >> 2) & 0x03) ;
   nodeindex = insertNybble(nodeindex,keybyte & 0x03) ;
#else
#  error No code for given BITS_PER_LEVEL
#endif      
   return ;
}


//----------------------------------------------------------------------

bool NybbleTrie::insert(const uint8_t *key, unsigned keylength,
			uint32_t frequency, bool stopgram)
{
   if (keylength > m_maxkeylen)
      m_maxkeylen = keylength ;
   uint32_t cur_index = NybbleTrie::ROOT_INDEX ;
   while (keylength > 0)
      {
      this->insertChild(cur_index,*key) ;
      key++ ;
      keylength-- ;
      }
   NybbleTrieNode *leaf = node(cur_index) ;
   bool new_node = false ;
   if (leaf)
      {
      new_node = (leaf->frequency() == 0) ;
      leaf->setFrequency(frequency) ;
      leaf->markAsLeaf() ;
      leaf->markAsStopgram(stopgram) ;
      }
   return new_node ;
}

//----------------------------------------------------------------------

bool NybbleTrie::insertMax(const uint8_t *key, unsigned keylength,
			   uint32_t frequency, bool stopgram)
{
   if (keylength > m_maxkeylen)
      m_maxkeylen = keylength ;
   uint32_t cur_index = NybbleTrie::ROOT_INDEX ;
   while (keylength > 0)
      {
      this->insertChild(cur_index,*key) ;
      key++ ;
      keylength-- ;
      }
   NybbleTrieNode *leaf = node(cur_index) ;
   bool new_node = false ;
   if (leaf)
      {
      uint32_t oldfreq = leaf->frequency() ;
      new_node = (oldfreq == 0) ;
      if (frequency > oldfreq)
	 leaf->setFrequency(frequency) ;
      leaf->markAsLeaf() ;
      leaf->markAsStopgram(stopgram) ;
      }
   return new_node ;
}

//----------------------------------------------------------------------

uint32_t NybbleTrie::find(const uint8_t *key, unsigned keylength) const
{
   uint32_t cur_index = NybbleTrie::ROOT_INDEX ;
   while (keylength > 0)
      {
      if (!extendKey(cur_index,*key))
	 return 0 ;
      key++ ;
      keylength-- ;
      }
   NybbleTrieNode *n = node(cur_index) ;
   return n ? n->frequency() : 0 ;
}

//----------------------------------------------------------------------

NybbleTrieNode* NybbleTrie::findNode(const uint8_t *key, unsigned keylength)
   const
{
   uint32_t cur_index = NybbleTrie::ROOT_INDEX ;
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

uint32_t NybbleTrie::increment(const uint8_t *key, unsigned keylength,
			       uint32_t incr, bool stopgram)
{
   uint32_t cur_index = NybbleTrie::ROOT_INDEX ;
   for (size_t i = 0 ; i < keylength ; i++)
      {
      if (!extendKey(cur_index,key[i]))
	 {
	 insert(key,keylength,incr,stopgram) ;
	 return incr ;
	 }
      }
   NybbleTrieNode *n = node(cur_index) ;
   if (n)
      {
      uint32_t freq = n->frequency() + incr ;
      n->setFrequency(freq) ;
      n->markAsLeaf() ;
      return freq ;
      }
   else
      {
      insert(key,keylength,incr,stopgram) ;
      return incr ;
      }
}

//----------------------------------------------------------------------

bool NybbleTrie::incrementExtensions(const uint8_t *key,
				     unsigned prevlength,
				     unsigned keylength,
				     uint32_t incr)
{
   uint32_t cur_index = NybbleTrie::ROOT_INDEX ;
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
      NybbleTrieNode *n = node(cur_index) ;
      if (!n)
	 return false ;
      uint32_t freq = n->frequency() + incr ;
      n->setFrequency(freq) ;
      n->markAsLeaf() ;
      }
   if (keylength > m_maxkeylen)
      m_maxkeylen = keylength ;
   return true ;
}

//----------------------------------------------------------------------

bool NybbleTrie::extendNybble(uint32_t &nodeindex, uint8_t nybble) const
{
   NybbleTrieNode *n = node(nodeindex) ;
   if (n->childPresent(nybble))
      {
      nodeindex = n->childIndex(nybble) ;
      return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool NybbleTrie::extendKey(uint32_t &nodeindex, uint8_t keybyte) const
{
   if (ignoringWhiteSpace() && keybyte == ' ')
      return true ;
   uint32_t idx = nodeindex ;
#if BITS_PER_LEVEL == 8
   if (extendNybble(idx,keybyte))
#elif BITS_PER_LEVEL == 4
   if (extendNybble(idx,keybyte >> 4) &&
       extendNybble(idx,keybyte & 0x0F))
#elif BITS_PER_LEVEL == 3
   if (extendNybble(idx,(keybyte >> 6) & 0x03) &&
       extendNybble(idx,(keybyte >> 3) & 0x07) &&
       extendNybble(idx,keybyte & 0x07))
#elif BITS_PER_LEVEL == 2
   if (extendNybble(idx,(keybyte >> 6) & 0x03) &&
       extendNybble(idx,(keybyte >> 4) & 0x03) &&
       extendNybble(idx,(keybyte >> 2) & 0x03) &&
       extendNybble(idx,keybyte & 0x03))
#else
#  error No code for given BITS_PER_LEVEL
#endif
      {
      nodeindex = idx ;
      return true ;
      }
   nodeindex = 0 ;
   return false ;
}

//----------------------------------------------------------------------

bool NybbleTrie::enumerate(uint8_t *keybuf, unsigned maxkeylength,
			   NybbleTrieEnumFn *fn, void *user_data) const
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

static bool scale_frequency(const NybbleTrieNode *node,
			    const uint8_t * /*key*/, unsigned /*keylen*/,
			    void *user_data)
{
   const uint64_t *total_count = (uint64_t*)user_data ;
   ((NybbleTrieNode*)node)->scaleFrequency(*total_count) ;
   return true ;			// continue iterating
}

//----------------------------------------------------------------------

bool NybbleTrie::scaleFrequencies(uint64_t total_count)
{
   if (m_nodes[0])
      {
      uint8_t keybuf[10000] ;
      return m_nodes[0]->enumerateChildren(this,keybuf,8*sizeof(keybuf),0,
					   scale_frequency,&total_count) ;
      }
   else
      return false ;
}

//----------------------------------------------------------------------

class CountAndPower
   {
   public:
      uint64_t count ;
      double power ;
      double log_power ;
      CountAndPower(uint64_t c,double p, double l)
	 { count = c ; power = p ; log_power = l ; }
   } ;

static bool scale_frequency_smoothed(const NybbleTrieNode *node,
				     const uint8_t * /*key*/,
				     unsigned /*keylen*/,
				     void *user_data)
{
   const CountAndPower *c_p = (CountAndPower*)user_data ;
   const uint64_t total_count = c_p->count ;
   double power = c_p->power ;
   double log_power = c_p->log_power ;
   ((NybbleTrieNode*)node)->scaleFrequency(total_count,power,
					   log_power) ;
   return true ;			// continue iterating
}

//----------------------------------------------------------------------

bool NybbleTrie::scaleFrequencies(uint64_t total_count, double power,
				  double log_power)
{
   if (m_nodes[0])
      {
      uint8_t keybuf[10000] ;
      CountAndPower c_p(total_count,power,log_power) ;
      return m_nodes[0]->enumerateChildren(this,keybuf,8*sizeof(keybuf),0,
					   scale_frequency_smoothed,&c_p) ;
      }
   else
      return false ;
}

//----------------------------------------------------------------------

uint32_t NybbleTrie::numTerminalNodes(uint32_t min_freq) const
{
   uint32_t count = 0 ;
   return (m_nodes[0]->enumerateTerminalNodes(this,0,count,min_freq) ? count : 0) ;
}

//----------------------------------------------------------------------

uint32_t NybbleTrie::numFullNodes(uint32_t min_freq) const
{
   uint32_t count = 1 ; // root node is always a full node
   return (m_nodes[0]->enumerateFullNodes(this,count,min_freq) ? count : 1) ;
}

//----------------------------------------------------------------------

NybbleTrie *NybbleTrie::load(Fr::CFile& f)
{
   if (f)
      {
//TODO
      }
   return nullptr ;
}

//----------------------------------------------------------------------

NybbleTrie *NybbleTrie::load(const char *filename)
{
   Fr::CInputFile fp(filename) ;
   return load(fp) ;
}

//----------------------------------------------------------------------

bool NybbleTrie::write(Fr::CFile& f) const
{
   if (f)
      {
//TODO
      }
   return false ;
}


// end of file trie.C //
