/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: trie.C - Word-frequency trie					*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-09						*/
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
#include "framepac/texttransforms.h"

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
      auto child = trie->node(m_children[i]) ;
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

uint32_t NybbleTrieNode::childIndex(unsigned int N) const
{
   return (N < lengthof(m_children)) ? m_children[N] : NybbleTrie::NULL_INDEX ;
}

//----------------------------------------------------------------------

uint32_t NybbleTrieNode::insertChild(unsigned int N, NybbleTrie *trie)
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

void NybbleTrieNode::scaleFrequency(uint64_t total_count)
{
   m_frequency = scaled_frequency(m_frequency,total_count) ;
   return ;
}

//----------------------------------------------------------------------

void NybbleTrieNode::scaleFrequency(uint64_t total_count, double power, double log_power)
{
   m_frequency = scaled_frequency(m_frequency,total_count,power,log_power) ;
   return ;
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

void NybbleTrie::init(uint32_t cap)
{
   m_userdata = nullptr ;
   m_maxkeylen = 0 ;
   m_totaltokens = 0 ;
   m_ignorewhitespace = false ;
   if (cap == 0)
      cap = 1 ;
   m_nodes.reserve(cap) ;
   auto root = m_nodes.alloc() ;
   // initialize the root node
   new (node(root)) NybbleTrieNode ;
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
	 lineptr = Fr::trim_whitespace(freq_end) ;
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

NybbleTrieNode *NybbleTrie::rootNode() const
{
   return node(NybbleTrie::ROOT_INDEX) ;
}

//----------------------------------------------------------------------

NybbleTrie::NodeIndex NybbleTrie::insertNybble(NodeIndex nodeindex, uint8_t nybble)
{
   auto n = node(nodeindex) ;
   return n->insertChild(nybble,this) ;
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

NybbleTrie::NodeIndex NybbleTrie::insertKey(const uint8_t* key, unsigned keylength)
{
   if (keylength > m_maxkeylen)
      m_maxkeylen = keylength ;
   auto cur_index = ROOT_INDEX ;
   while (keylength > 0)
      {
      this->insertChild(cur_index,*key) ;
      key++ ;
      keylength-- ;
      }
   return cur_index ;
}

//----------------------------------------------------------------------

bool NybbleTrie::insert(const uint8_t *key, unsigned keylength,
			uint32_t frequency, bool stopgram)
{
   auto leaf = node(insertKey(key,keylength)) ;
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
   auto leaf = node(cur_index) ;
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

NybbleTrie::NodeIndex NybbleTrie::findKey(const uint8_t *key, unsigned keylength) const
{
   auto cur_index = ROOT_INDEX ;
   while (keylength > 0)
      {
      if (!extendKey(cur_index,*key))
	 return INVALID_INDEX ;
      key++ ;
      keylength-- ;
      }
   return cur_index ;
}

//----------------------------------------------------------------------

uint32_t NybbleTrie::increment(const uint8_t *key, unsigned keylength,
			       uint32_t incr, bool stopgram)
{
   auto cur_index = ROOT_INDEX ;
   for (size_t i = 0 ; i < keylength ; i++)
      {
      if (!extendKey(cur_index,key[i]))
	 {
	 insert(key,keylength,incr,stopgram) ;
	 return incr ;
	 }
      }
   auto n = node(cur_index) ;
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
   auto cur_index = ROOT_INDEX ;
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
      auto freq = n->frequency() + incr ;
      n->setFrequency(freq) ;
      n->markAsLeaf() ;
      }
   if (keylength > m_maxkeylen)
      m_maxkeylen = keylength ;
   return true ;
}

//----------------------------------------------------------------------

bool NybbleTrie::extendNybble(NodeIndex& nodeindex, uint8_t nybble) const
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

bool NybbleTrie::extendKey(NodeIndex& nodeindex, uint8_t keybyte) const
{
   if (ignoringWhiteSpace() && keybyte == ' ')
      return true ;
   auto idx = nodeindex ;
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
   nodeindex = NULL_INDEX ;
   return false ;
}

//----------------------------------------------------------------------

bool NybbleTrie::singleChild(NodeIndex nodeindex) const 
{
   auto node = this->node(nodeindex) ;
   for (size_t i = 0 ; i < 8 && node ; i += BITS_PER_LEVEL)
      {
      unsigned index = INVALID_INDEX ;
      for (unsigned ch = 0 ; ch < (1<<BITS_PER_LEVEL) ; ch++)
	 {
	 if (node->childPresent(ch))
	    {
	    if (index != INVALID_INDEX)
	       return false ; 		// multiple children
	    index = ch ;
	    }
	 }
      if (index == INVALID_INDEX)
	 return false ; 		// no children at all
      node = this->node(node->childIndex(index)) ;
      }
   return node != nullptr ;
}

//----------------------------------------------------------------------

bool NybbleTrie::singleChildSameFreq(NodeIndex nodeindex, bool allow_nonleaf, double ratio) const 
{
   auto node = this->node(nodeindex) ;
   auto parent_freq = node->frequency() ;
   for (size_t i = 0 ; i < 8 && node ; i += BITS_PER_LEVEL)
      {
      unsigned index = INVALID_INDEX ;
      for (unsigned ch = 0 ; ch < (1<<BITS_PER_LEVEL) ; ch++)
	 {
	 if (node->childPresent(ch))
	    {
	    if (index != INVALID_INDEX)
	       return false ; 		// multiple children
	    index = ch ;
	    }
	 }
      if (index == INVALID_INDEX)
	 return false ; 		// no children at all
      node = this->node(node->childIndex(index)) ;
      }
   if (!node)
      return false ;
   auto freq = node->frequency() ;
   return ((freq <= parent_freq && freq >= ratio * parent_freq) || (allow_nonleaf && freq == 0)) ;
}

//----------------------------------------------------------------------

bool NybbleTrie::countTerminalNodes(NodeIndex nodeindex, unsigned keylen_bits,
				    uint32_t& count, uint32_t min_freq) const
{
   auto node = this->node(nodeindex) ;
   if (!node->hasChildren())
      return true ;
   else if ((keylen_bits % 8) == 0 && allChildrenAreTerminals(nodeindex,min_freq))
      {
      count += numExtensions(nodeindex,min_freq) ;
      return true ;
      }
   keylen_bits = keylen_bits + BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; ++i)
      {
      uint32_t child = node->childIndex(i) ;
      if (child != NybbleTrie::NULL_INDEX && !countTerminalNodes(child,keylen_bits,count,min_freq))
	 return false ;
      }
   return true ;
}

//----------------------------------------------------------------------

unsigned NybbleTrie::numExtensions(NodeIndex nodeindex, uint32_t min_freq, unsigned bits) const
{
   unsigned count = 0 ;
   auto node = this->node(nodeindex) ;
   if (bits >= 8)
      {
      return (node && node->frequency() >= min_freq) ? 1 : 0 ;
      }
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; ++i)
      {
      auto childindex = node->childIndex(i) ;
      if (childindex != NULL_INDEX)
	 count += numExtensions(childindex,min_freq,bits+BITS_PER_LEVEL) ;
      }
   return count ;
}

//----------------------------------------------------------------------

bool NybbleTrie::allChildrenAreTerminals(NodeIndex nodeindex, uint32_t min_freq, unsigned bits) const
{
   auto node = this->node(nodeindex) ;
   if (bits >= 8)
      {
      return ! ((!node->leaf() || node->frequency() >= min_freq) && numExtensions(nodeindex,min_freq) > 0) ;
      }
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; ++i)
      {
      auto childindex = node->childIndex(i) ;
      if (childindex != NULL_INDEX && !allChildrenAreTerminals(childindex,min_freq,bits+BITS_PER_LEVEL))
	 return false ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool NybbleTrie::enumerate(uint8_t *keybuf, unsigned maxkeylength, EnumFn *fn, void *user_data) const
{
   if (keybuf && fn)
      {
      memset(keybuf,'\0',maxkeylength) ;
      return enumerateChildren(ROOT_INDEX,keybuf,maxkeylength*8,0,fn,user_data) ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool NybbleTrie::enumerateChildren(NodeIndex nodeindex,
				       uint8_t *keybuf,
				       unsigned max_keylength_bits,
				       unsigned curr_keylength_bits,
				       EnumFn *fn,
				       void *user_data) const
{
   auto n = node(nodeindex) ;
   if (n->leaf() && !fn(this,nodeindex,keybuf,curr_keylength_bits/8,user_data))
      return false ;
   if (curr_keylength_bits < max_keylength_bits)
      {
      unsigned curr_bits = curr_keylength_bits + BITS_PER_LEVEL ;
      for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
	 {
	 uint32_t child = n->childIndex(i) ;
	 if (child != NybbleTrie::NULL_INDEX)
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
	    if (!enumerateChildren(child,keybuf,max_keylength_bits,curr_bits,fn,user_data))
	       return false ;
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

static bool scale_frequency(const NybbleTrie* trie, NybbleTrie::NodeIndex nodeindex,
			    const uint8_t * /*key*/, unsigned /*keylen*/,
			    void *user_data)
{
   const uint64_t *total_count = (uint64_t*)user_data ;
   auto node = trie->node(nodeindex) ;
   node->scaleFrequency(*total_count) ;
   return true ;			// continue iterating
}

//----------------------------------------------------------------------

bool NybbleTrie::scaleFrequencies(uint64_t total_count)
{
   uint8_t keybuf[10000] ;
   return enumerateChildren(ROOT_INDEX,keybuf,8*sizeof(keybuf),0,scale_frequency,&total_count) ;
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

static bool scale_frequency_smoothed(const NybbleTrie* trie, NybbleTrie::NodeIndex nodeindex,
				     const uint8_t * /*key*/,
				     unsigned /*keylen*/,
				     void *user_data)
{
   const CountAndPower *c_p = (CountAndPower*)user_data ;
   const uint64_t total_count = c_p->count ;
   auto node = trie->node(nodeindex) ;
   node->scaleFrequency(total_count,c_p->power,c_p->log_power) ;
   return true ;			// continue iterating
}

//----------------------------------------------------------------------

bool NybbleTrie::scaleFrequencies(uint64_t total_count, double power,
				  double log_power)
{
   uint8_t keybuf[10000] ;
   CountAndPower c_p(total_count,power,log_power) ;
   return enumerateChildren(ROOT_INDEX,keybuf,8*sizeof(keybuf),0,scale_frequency_smoothed,&c_p) ;
}

//----------------------------------------------------------------------

uint32_t NybbleTrie::numTerminalNodes(uint32_t min_freq) const
{
   uint32_t count = 0 ;
   return countTerminalNodes(ROOT_INDEX,0,count,min_freq) ? count : 0 ;
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
