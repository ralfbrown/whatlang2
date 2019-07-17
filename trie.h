/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: trie.h - Word-frequency trie					*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-16						*/
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

#include <cstdint>
#include <limits.h>
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
      uint32_t insertChild(unsigned int N, NybbleTrie *trie) ;

      // I/O
      static NybbleTrieNode *read(Fr::CFile& f) ;
      bool write(Fr::CFile& f) const ;
   protected:
      uint32_t	m_children[1<<BITS_PER_LEVEL] ;
      uint32_t  m_frequency { 0 } ;
      bool	m_leaf { false } ;
      bool	m_stopgram { false } ;
   } ;

//----------------------------------------------------------------------

class NybbleTrie
   {
   public:
      typedef uint32_t NodeIndex ;
      static constexpr NodeIndex ROOT_INDEX = 0U ;
      static constexpr NodeIndex NULL_INDEX = 0U ;
      static constexpr NodeIndex INVALID_INDEX = (uint32_t)~0 ;
      typedef NybbleTrieNode Node ;
      typedef bool LoadFn(NybbleTrie* trie, const uint8_t* key, unsigned keylen, uint32_t langID, uint32_t freq) ;
      typedef bool EnumFn(const NybbleTrie *trie, NodeIndex index, const uint8_t *key, unsigned keylen,
	                  void *user_data) ;
   public:
      NybbleTrie(NodeIndex capacity = 0) ;
      NybbleTrie(const char *filename, bool verbose) ;
      ~NybbleTrie() = default ;

      bool loadWords(const char *filename, LoadFn* insfn, uint32_t langID, bool verbose) ;
      NodeIndex allocateNode() { return (NodeIndex)m_nodes.alloc() ; }

      // modifiers
      void setUserData(void *ud) { m_userdata = ud ; }
      void ignoreWhiteSpace(bool ignore = true) { m_ignorewhitespace = ignore ; }
      void addTokenCount(uint32_t incr = 1) { m_totaltokens += incr ; }
      NodeIndex insertKey(const uint8_t* key, unsigned keylength) ;
      bool insert(const uint8_t *key, unsigned keylength, uint32_t frequency, bool stopgram) ;
      bool insertMax(const uint8_t *key, unsigned keylength, uint32_t frequency, bool stopgram) ;
      void insertChild(uint32_t &nodeindex, uint8_t keybyte) ;
      uint32_t increment(const uint8_t *key, unsigned keylength, uint32_t incr = 1, bool stopgram = false) ;
      bool incrementExtensions(const uint8_t *key, unsigned prevlength,
			       unsigned keylength, uint32_t incr = 1) ;

      // accessors
      void *userData() const { return m_userdata ; }
      uint32_t size() const { return m_nodes.size() ; }
      uint32_t capacity() const { return m_nodes.capacity() ; }
      uint32_t totalTokens() const { return m_totaltokens ; }
      unsigned longestKey() const { return m_maxkeylen ; }
      bool ignoringWhiteSpace() const { return m_ignorewhitespace ; }
      Node* node(NodeIndex N) const { return m_nodes.item(N) ; }
      Node* rootNode() const { return node(ROOT_INDEX) ; }
      NodeIndex findKey(const uint8_t* key, unsigned keylength) const ;
      Node* findNode(const uint8_t* key, unsigned keylength) const { return node(findKey(key,keylength)) ; }
      uint32_t find(const uint8_t* key, unsigned keylength) const ;
      bool extendKey(NodeIndex& nodeindex, uint8_t keybyte) const ;
      bool enumerate(uint8_t* keybuf, unsigned maxkeylength, EnumFn *fn, void *user_data) const ;
      bool enumerateChildren(NodeIndex nodeindex, uint8_t *keybuf, unsigned max_keylength_bits,
			     unsigned curr_keylength_bits, EnumFn *fn, void *user_data) const ;
      bool singleChild(NodeIndex nodeindex) const ;
      bool singleChildSameFreq(NodeIndex nodeindex, bool allow_nonleaf, double ratio) const ;
      bool allChildrenAreTerminals(NodeIndex nodeindex, uint32_t min_freq, unsigned = 0) const ;
      unsigned numExtensions(NodeIndex nodeindex, uint32_t min_freq, unsigned = 0) const ;
      uint32_t numFullByteNodes(uint32_t min_freq) const ;
      uint32_t numTerminalNodes(uint32_t min_freq) const ;

      bool scaleFrequencies(uint64_t total_count) ;
      bool scaleFrequencies(uint64_t total_count, double power, double log_power) ;

      // I/O
      static NybbleTrie *load(Fr::CFile& f) ;
      static NybbleTrie *load(const char *filename) ;
      bool write(Fr::CFile& f) const ;
   protected:
      NodeIndex insertNybble(NodeIndex nodeindex, uint8_t nybble) ;
   private:
      void init(uint32_t capacity) ;
      bool extendNybble(NodeIndex& nodeindex, uint8_t nybble) const ;
      size_t countTerminalNodes(NodeIndex nodeindex, uint32_t min_freq, unsigned keylen_bits = 0) const ;
      size_t countFullByteNodes(NodeIndex nodeindex, uint32_t min_freq, unsigned keylen_bits = 0) const ;
   protected:
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
      TriePointer() = default ;
      TriePointer(T* t) : m_trie(t) { this->resetKey() ; }
      TriePointer(const T* t) : m_trie(const_cast<T*>(t)) { this->resetKey() ; }
      ~TriePointer() { invalidate() ; }
      TriePointer& operator= (const TriePointer& orig) = default ;

      // manipulators
      void setTrie(T* t) { m_trie = t ; }
      void resetKey() { m_index = T::ROOT_INDEX ; m_keylen = 0 ; m_valid = true ; }
      void invalidate() { m_valid = false ; }
      bool extendKey(uint8_t keybyte)
	 {
	    if (!m_valid) return false ;
	    bool success = m_trie->extendKey(m_index,keybyte) ;
	    if (success)
	       ++m_keylen ;
	    else
	       m_valid = false ;
	    return success ;
	 }
	 
      // accessors
      bool OK() const
	 {
	    auto n = node() ;
	    return n && n->leaf() ;
	 }
      bool valid() const { return m_valid ; }
      explicit operator bool () const { return m_valid ; }
      bool operator ! () const { return !m_valid ; }
      unsigned keyLength() const { return m_keylen ; }
      typename T::Node* node() const { return m_valid ? m_trie->node(m_index) : nullptr ; }

   private:
      T*  	 m_trie { nullptr } ;
      uint32_t	 m_index ;
      uint16_t	 m_keylen ;
      bool       m_valid { false } ;
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
