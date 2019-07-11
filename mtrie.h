/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: mtrie.h - bit-slice-based Word-frequency multi-trie		*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-07						*/
/*									*/
/*  (c) Copyright 2011,2012,2019 Carnegie Mellon University		*/
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

#ifndef __MTRIE_H_INCLUDED
#define __MTRIE_H_INCLUDED

#include <cstdio>
#include <limits.h>
#include <stdint.h>
#include "trie.h"
#include "framepac/byteorder.h"
#include "framepac/file.h"
#include "framepac/trie.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

// we can trade off speed for memory by adjusting how many bits each
//   node in the trie represents.  Current supported values are 2 and 4;
//   two bits per node uses about 60% as much total memory as 4 bits,
//   but needs twice as many memory accesses for lookups
#define MTRIE_BITS_PER_LEVEL 2

#define LID_LANGID_MASK   0x0FFFFFFF
#define LID_STOPGRAM_MASK 0x8000000

/************************************************************************/
/************************************************************************/

class MultiTrieNode ;

class MultiTrie ;

//----------------------------------------------------------------------

class MultiTrieFrequency
   {
   public:
      static constexpr uint32_t INVALID_FREQ = (uint32_t)~0 ;
   public:
      MultiTrieFrequency() {} // only for internal use by readAll()
      MultiTrieFrequency(uint32_t freq, uint32_t langID,
			 bool stopgram = false,
			 MultiTrieFrequency *nxt = nullptr)
	 { m_frequency = freq ;
	   m_langID = (langID & LID_LANGID_MASK) | (stopgram * LID_STOPGRAM_MASK) ;
	   setNext(nxt) ; }
      ~MultiTrieFrequency() ;
      static MultiTrieFrequency *allocate(uint32_t freq = 0,
					  uint32_t langID = 0,
					  bool stopgram = false) ;

      // accessors
      static MultiTrieFrequency *baseAddress()
	 { return s_base_address ; }
      uint32_t frequency() const { return m_frequency ; }
      uint32_t frequency(uint32_t langID) const ;
      uint32_t languageID() const { return m_langID & LID_LANGID_MASK ; }
      bool isStopgram() const { return (m_langID & LID_STOPGRAM_MASK) != 0 ; }
      static MultiTrieFrequency *getAddress(uint32_t index)
         { return (index == INVALID_FREQ) ? 0 : baseAddress() + index ; }
      static uint32_t getIndex(MultiTrieFrequency *f)
         { return f ? (uint32_t)(f - baseAddress()) : INVALID_FREQ ; }
      MultiTrieFrequency *next() const { return getAddress(m_next) ; }

      // manipulators
      static void setBaseAddress(MultiTrieFrequency *base,
				 uint32_t max_alloc)
	 { s_base_address = base ; s_max_alloc = max_alloc ; }
      void setNext(MultiTrieFrequency *nxt)
	 { m_next = nxt ? (nxt - s_base_address) : INVALID_FREQ ; }
      void setNext(uint32_t nxt) { m_next = nxt ; }
      void setFrequency(uint32_t freq) { m_frequency = freq ; }
      void setFrequency(uint32_t ID, uint32_t freq, bool stopgram) ;
      void incrFrequency(uint32_t incr = 1) { m_frequency += incr ; }
      void incrFrequency(uint32_t ID, uint32_t incr) ;
      void scaleFrequency(uint64_t total_count, uint32_t ID) ;

      // I/O
      static MultiTrieFrequency *read(Fr::CFile& f) ;
      static bool readAll(Fr::CFile& f) ;
      bool write(Fr::CFile& f) const ;
      static bool writeAll(Fr::CFile& f) ;
   private:
      static MultiTrieFrequency *s_base_address ;
      static uint32_t s_max_alloc ;
      static uint32_t s_curr_alloc ;
      uint32_t m_next ;
      uint32_t m_frequency ;
      uint32_t m_langID ;
   } ;

//----------------------------------------------------------------------

class MultiTrieNode ;

//TODO: implement template in FramepaC-ng, make this a thin wrapper around it
class LangIDMultiTrie //: public Fr::MultiTrie<Fr::UInt32>
   {
   public:
      static constexpr uint32_t ROOT_INDEX = 0U ;
      static constexpr uint32_t NULL_INDEX = 0U ;
      static constexpr uint32_t INVALID_FREQ = MultiTrieFrequency::INVALID_FREQ ;

      typedef MultiTrieNode Node ;
      typedef bool EnumFn(const MultiTrieNode *node, const uint8_t *key,
	 		  unsigned keylen, void *user_data) ;
   public:
      LangIDMultiTrie(uint32_t capacity = 0) ;
      LangIDMultiTrie(const char *filename, bool verbose) ;
      LangIDMultiTrie(const class LangIDMultiTrie *) ;
      LangIDMultiTrie(const class LangIDPackedMultiTrie *) ;
      ~LangIDMultiTrie() ;

      bool loadWords(const char *filename, uint32_t langID,
		     bool verbose = false) ;
      uint32_t allocateNode() ;

      // modifiers
      void ignoreWhiteSpace(bool ignore = true) { m_ignorewhitespace = ignore ; }
      void setLanguage(uint32_t langID) { m_currentlangID = langID ; }
      void addTokenCount(uint32_t incr = 1) { m_totaltokens += incr ; }
      bool insert(const uint8_t *key, unsigned keylength,
		  uint32_t langID, uint32_t frequency, bool stopgram) ;
      void insertChild(uint32_t &nodeindex, uint8_t keybyte) ;
      uint32_t increment(const uint8_t *key, unsigned keylength,
			 uint32_t langID, uint32_t incr = 1,
			 bool stopgram = false) ;
      bool incrementExtensions(const uint8_t *key, unsigned prevlength,
			       unsigned keylength, uint32_t langID,
			       uint32_t incr = 1) ;

      // accessors
      uint32_t size() const { return m_used ; }
      uint32_t capacity() const { return m_capacity ; }
      uint32_t currentLanguage() const { return m_currentlangID ; }
      uint32_t totalTokens() const { return m_totaltokens ; }
      unsigned longestKey() const { return m_maxkeylen ; }
      bool ignoringWhiteSpace() const { return m_ignorewhitespace ; }
      MultiTrieNode *node(uint32_t N) const ;
      MultiTrieNode *rootNode() const ;
      MultiTrieNode *findNode(const uint8_t *key, unsigned keylength) const ;
      uint32_t find(const uint8_t *key, unsigned keylength) const ;
      bool extendKey(uint32_t &nodeindex, uint8_t keybyte) const ;
      bool enumerate(uint8_t *keybuf, unsigned maxkeylength, EnumFn *fn, void *user_data) const ;
      bool enumerateChildren(uint32_t nodeindex, uint8_t *keybuf, unsigned max_keylength_bits,
			     unsigned curr_keylength_bits, EnumFn *fn, void *user_data) const ;
      bool singleChild(uint32_t nodeindex) const ;
      bool singleChildSameFreq(uint32_t nodeindex, bool allow_nonleaf, double ratio) const ;
      uint32_t countFreqRecords() const ;
      uint32_t numFullByteNodes() const ;
      uint32_t numTerminalNodes() const ;

      // I/O
      static LangIDMultiTrie *load(Fr::CFile& f) ;
      static LangIDMultiTrie *load(const char *filename) ;
      bool write(Fr::CFile& f) const ;
      bool write(const char *filename) const ;
      bool dump(Fr::CFile& f) const ;
   private:
      void init(uint32_t capacity) ;
      bool writeHeader(Fr::CFile& f) const ;
      bool extendNybble(uint32_t &nodeindex, uint8_t nybble) const ;
      uint32_t insertNybble(uint32_t nodeindex, uint8_t nybble) ;
   private:
      MultiTrieNode     **m_nodes ;	 // list of buckets of nodes
      uint32_t	 	  m_capacity ;	 // number of nodes allocated
      uint32_t	 	  m_used ;	 // number of nodes in use
      uint32_t		  m_totaltokens ;
      uint32_t		  m_currentlangID ;
      unsigned		  m_maxkeylen ;
      bool		  m_ignorewhitespace ;
   } ;

//----------------------------------------------------------------------

class MultiTrieNode
   {
   public:
      void *operator new(size_t, void *where) { return where ; }
      MultiTrieNode() ;
      ~MultiTrieNode() {}

      // accessors
      bool leaf() const { return m_isleaf ; }
      bool isStopgram(unsigned langID) const ;
      bool hasChildren() const ;
      bool childPresent(unsigned int N) const ;
      uint32_t childIndex(unsigned int N) const ;
      unsigned numExtensions(const LangIDMultiTrie *trie) const ;
      bool allChildrenAreTerminals(const LangIDMultiTrie *trie) const ;
      uint32_t frequency(uint32_t ID = 0) const ;
      unsigned numFrequencies() const ;
      MultiTrieFrequency *frequencies() const
         { return MultiTrieFrequency::getAddress(m_frequency_info) ; }
      bool enumerateFullByteNodes(const LangIDMultiTrie *trie,
				  unsigned keylen_bits, uint32_t &count) const ;
      bool enumerateTerminalNodes(const LangIDMultiTrie *trie,
				  unsigned keylen_bits, uint32_t &count) const ;

      // modifiers
      void markAsLeaf() { m_isleaf = true ; }
      void setFrequency(uint32_t ID, uint32_t f, bool stopgram) ;
      bool setFrequencies(MultiTrieFrequency *freqs) ;
      bool insertChild(unsigned int N, LangIDMultiTrie *trie) ;

      // I/O
      bool load(Fr::CFile& f) ;
      bool write(Fr::CFile& f) const ;

   private:
      uint32_t	m_children[1<<MTRIE_BITS_PER_LEVEL] ;
      uint32_t  m_frequency_info ;
      bool	m_isleaf ;
   } ;

//----------------------------------------------------------------------

typedef TriePointer<LangIDMultiTrie> MultiTriePointer ;

#endif /* !__MTRIE_H_INCLUDED */

/* end of file trie.h */
