/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: trie.h - Word-frequency trie					*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-11						*/
/*									*/
/*  (c) Copyright 2011,2012,2014,2019 Carnegie Mellon University	*/
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

#ifndef __TRIE_H_INCLUDED
#define __TRIE_H_INCLUDED

#include <cstdio>
#include <limits.h>
#include <stdint.h>
#include "framepac/file.h"
#include "framepac/itempool.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

// we can trade off speed for memory by adjusting how many bits each
//   node in the trie represents.  Current supported values are 2, 3,
//   and 4; two bits per node uses about 60% as much total memory as 4
//   bits, but needs twice as many memory accesses for lookups; three bits
//   is in-between
#define BITS_PER_LEVEL 2

// we want to store percentages for entries in the trie in 32 bits.  Since
//   it is very unlikely that any ngram in the trie will have a probability
//   greater than 4.2%, scale the percentage by a factor of one billion.
#define TRIE_SCALE_FACTOR 1000000000L

/************************************************************************/
/************************************************************************/

class NybbleTrie ;

class WildcardSet ;
class WildcardCollection ;

//----------------------------------------------------------------------

class NybbleTrieNode
   {
   public:
      void *operator new(size_t, void *where) { return where ; }
      NybbleTrieNode() ;
      ~NybbleTrieNode() = default ;

      // accessors
      bool leaf() const { return m_leaf ; }
      bool isStopgram() const { return m_stopgram ; }
      bool hasChildren() const ;
      bool hasChildren(const NybbleTrie *trie, uint32_t min_freq) const ;
      bool childPresent(unsigned int N) const ;
      uint32_t childIndex(unsigned int N) const ;
      uint32_t frequency() const { return m_frequency ; }

      // modifiers
      void markAsLeaf() { m_leaf = true ; }
      void markAsStopgram(bool stop = true) { m_stopgram = stop ; }
      void setFrequency(uint32_t f) { m_frequency = f ; }
      void incrFrequency(uint32_t incr = 1) { m_frequency += incr ; }
      void scaleFrequency(uint64_t total_count) ;
      void scaleFrequency(uint64_t total_count, double power, double log_power) ;
      bool insertChild(unsigned int N, NybbleTrie *trie) ;

      // I/O
      static NybbleTrieNode *read(Fr::CFile& f) ;
      bool write(Fr::CFile& f) const ;
   protected:
      uint32_t	m_children[1<<BITS_PER_LEVEL] ;
      uint32_t  m_frequency ;
      bool	m_leaf ;
      bool	m_stopgram ;
   } ;

//----------------------------------------------------------------------

class NybbleTrie
   {
   public:
      static constexpr uint32_t ROOT_INDEX = 0U ;
      static constexpr uint32_t NULL_INDEX = 0U ;
      typedef NybbleTrieNode Node ;
      typedef bool EnumFn(const NybbleTrie *trie, uint32_t index, const uint8_t *key, unsigned keylen,
	                  void *user_data) ;
   public:
      NybbleTrie(uint32_t capacity = 0) ;
      NybbleTrie(const char *filename, bool verbose) ;
      ~NybbleTrie() = default ;

      bool loadWords(const char *filename, bool verbose) ;
      uint32_t allocateNode() { return (uint32_t)m_nodes.alloc() ; }

      // modifiers
      void setUserData(void *ud) { m_userdata = ud ; }
      void ignoreWhiteSpace(bool ignore = true) { m_ignorewhitespace = ignore ; }
      void addTokenCount(uint32_t incr = 1) { m_totaltokens += incr ; }
      bool insert(const uint8_t *key, unsigned keylength,
		  uint32_t frequency, bool stopgram) ;
      bool insertMax(const uint8_t *key, unsigned keylength,
		     uint32_t frequency, bool stopgram) ;
      void insertChild(uint32_t &nodeindex, uint8_t keybyte) ;
      uint32_t increment(const uint8_t *key, unsigned keylength,
			 uint32_t incr = 1, bool stopgram = false) ;
      bool incrementExtensions(const uint8_t *key, unsigned prevlength,
			       unsigned keylength, uint32_t incr = 1) ;

      // accessors
      void *userData() const { return m_userdata ; }
      uint32_t size() const { return m_nodes.size() ; }
      uint32_t capacity() const { return m_nodes.capacity() ; }
      uint32_t totalTokens() const { return m_totaltokens ; }
      unsigned longestKey() const { return m_maxkeylen ; }
      bool ignoringWhiteSpace() const { return m_ignorewhitespace ; }
      Node *node(uint32_t N) const { return m_nodes.item(N) ; }
      Node *rootNode() const ;
      Node *findNode(const uint8_t *key, unsigned keylength) const ;
      uint32_t find(const uint8_t *key, unsigned keylength) const ;
      bool extendKey(uint32_t &nodeindex, uint8_t keybyte) const ;
      bool enumerate(uint8_t *keybuf, unsigned maxkeylength, EnumFn *fn, void *user_data) const ;
      bool enumerateChildren(uint32_t nodeindex, uint8_t *keybuf, unsigned max_keylength_bits,
			     unsigned curr_keylength_bits, EnumFn *fn, void *user_data) const ;
      bool countTerminalNodes(uint32_t nodeindex, unsigned keylen_bits, uint32_t &count,
			      uint32_t min_freq = 0) const ;
      unsigned numExtensions(uint32_t nodeindex, uint32_t min_freq = 0, unsigned = 0) const ;
      bool allChildrenAreTerminals(uint32_t nodeindex, uint32_t min_freq = 0, unsigned = 0) const ;
      bool singleChild(uint32_t nodeindex) const ;
      bool singleChildSameFreq(uint32_t nodeindex, bool allow_nonleaf, double ratio) const ;

      bool scaleFrequencies(uint64_t total_count) ;
      bool scaleFrequencies(uint64_t total_count, double power, double log_power) ;
      uint32_t numTerminalNodes(uint32_t min_freq = 0) const ;

      // I/O
      static NybbleTrie *load(Fr::CFile& f) ;
      static NybbleTrie *load(const char *filename) ;
      bool write(Fr::CFile& f) const ;
   private:
      void init(uint32_t capacity) ;
      bool extendNybble(uint32_t &nodeindex, uint8_t nybble) const ;
      uint32_t insertNybble(uint32_t nodeindex, uint8_t nybble) ;
   private:
      Fr::ItemPool<Node> m_nodes ;
      void	       *m_userdata ;
      uint32_t	 	m_totaltokens ;
      unsigned	 	m_maxkeylen ;
      bool	 	m_ignorewhitespace ;
   } ;

//----------------------------------------------------------------------

template <class T>
class TriePointer
   {
   public:
      TriePointer() : m_trie(nullptr), m_valid(false) {}
      TriePointer(T* t) : m_trie(t) { this->resetKey() ; }
      TriePointer(const T* t) : m_trie(const_cast<T*>(t)) { this->resetKey() ; }
      ~TriePointer() { m_trie = nullptr ; m_index = 0 ; }
      TriePointer& operator= (const TriePointer& orig) = default ;

      // manipulators
      void setTrie(T* t) { m_trie = t ; }
      void resetKey() { m_index = T::ROOT_INDEX ; m_keylen = 0 ; m_failed = false ; m_valid = true ; }
      void invalidate() { m_valid = false ; }
      bool extendKey(uint8_t keybyte)
	 {
	    if (m_failed) return false ;
	    bool success = m_trie->extendKey(m_index,keybyte) ;
	    if (success)
	       ++m_keylen ;
	    else
	       m_failed = true ;
	    return success ;
	 }
	 
      // accessors
      bool OK() const
	 {
	    if (m_failed) return false ;
	    auto n = m_trie->node(m_index) ;
	    return n && n->leaf() ;
	 }
      bool valid() const { return m_valid ; }
      explicit operator bool () const { return m_valid ; }
      bool operator ! () const { return !m_valid ; }
      unsigned keyLength() const { return m_keylen ; }
      typename T::Node* node() const { return m_failed ? nullptr : m_trie->node(m_index) ; }

   private:
      T*  	 m_trie ;
      uint32_t	 m_index ;
      uint16_t	 m_keylen ;
      bool	 m_failed ;
      bool       m_valid ;
   } ;

typedef TriePointer<NybbleTrie> NybbleTriePointer ;

/************************************************************************/
/************************************************************************/

double scaling_log_power(double power) ;
double scale_frequency(double freq, double power, double log_power) ;
double unscale_frequency(uint32_t scaled, double power) ;
uint32_t scaled_frequency(uint32_t raw_freq, uint64_t total_count,
			  double power, double log_power) ;

#endif /* !__TRIE_H_INCLUDED */

/* end of file trie.h */
