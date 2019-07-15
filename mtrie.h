/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: mtrie.h - bit-slice-based Word-frequency multi-trie		*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-11						*/
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
#include "framepac/itempool.h"
#include "framepac/file.h"
#include "framepac/trie.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define LID_LANGID_MASK   0x0FFFFFFF
#define LID_STOPGRAM_MASK 0x8000000

/************************************************************************/
/************************************************************************/

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
      uint32_t frequency() const { return m_frequency ; }
      uint32_t frequency(uint32_t langID) const ;
      uint32_t languageID() const { return m_langID & LID_LANGID_MASK ; }
      bool isStopgram() const { return (m_langID & LID_STOPGRAM_MASK) != 0 ; }
      static MultiTrieFrequency *getAddress(uint32_t index) { return s_freq_records.item(index) ; }
      static uint32_t getIndex(MultiTrieFrequency *f)
         { return f >= s_freq_records.begin() ? (uint32_t)(f - s_freq_records.begin()) : INVALID_FREQ ; }
      MultiTrieFrequency *next() const { return getAddress(m_next) ; }

      // manipulators
      void setNext(MultiTrieFrequency *nxt)
	 { m_next = nxt ? (nxt - s_freq_records.begin()) : INVALID_FREQ ; }
      void setNext(uint32_t nxt) { m_next = nxt ; }
      void setFrequency(uint32_t freq) { m_frequency = freq ; }
      void setFrequency(uint32_t ID, uint32_t freq, bool stopgram) ;
      void incrFrequency(uint32_t incr = 1) { m_frequency += incr ; }
      void incrFrequency(uint32_t ID, uint32_t incr) ;
      void scaleFrequency(uint64_t total_count, uint32_t ID) ;

      // I/O
      static bool readAll(Fr::CFile& f) ;
      bool write(Fr::CFile& f) const ;
      static bool writeAll(Fr::CFile& f) ;
   protected:
      void newFrequency(uint32_t ID, uint32_t freq, bool stopgram) ;

   private:
      static Fr::ItemPool<MultiTrieFrequency> s_freq_records ;
      uint32_t m_next ;
      uint32_t m_frequency ;
      uint32_t m_langID ;
   } ;

//----------------------------------------------------------------------

class LangIDMultiTrie ;

class MultiTrieNode : public NybbleTrieNode
   {
   public:
      void *operator new(size_t, void *where) { return where ; }
      MultiTrieNode() ;
      ~MultiTrieNode() = default ;

      // accessors
      bool isStopgram(unsigned langID) const ;
      uint32_t frequency(uint32_t ID = 0) const ;
      unsigned numFrequencies() const ;
      MultiTrieFrequency *frequencies() const { return MultiTrieFrequency::getAddress(m_frequency) ; }

      // modifiers
      void setFrequency(uint32_t ID, uint32_t f, bool stopgram) ;
      bool setFrequencies(MultiTrieFrequency *freqs) ;
      uint32_t insertChild(unsigned int N, LangIDMultiTrie *trie) ;

      // I/O
      bool load(Fr::CFile& f) ;
      bool write(Fr::CFile& f) const ;
   } ;

//----------------------------------------------------------------------

//TODO: implement template in FramepaC-ng, make this a thin wrapper around it
class LangIDMultiTrie : public NybbleTrie
   {
   public:
      static constexpr uint32_t ROOT_INDEX = 0U ;
      static constexpr uint32_t NULL_INDEX = 0U ;
      static constexpr uint32_t INVALID_FREQ = MultiTrieFrequency::INVALID_FREQ ;

      typedef MultiTrieNode Node ;
      typedef bool EnumFn(const LangIDMultiTrie* trie, NodeIndex nodeindex, const uint8_t *key,
	 		  unsigned keylen, void *user_data) ;
   public:
      LangIDMultiTrie(uint32_t capacity = 0) ;
      LangIDMultiTrie(const char *filename, bool verbose) ;
      LangIDMultiTrie(const class LangIDMultiTrie *) ;
      LangIDMultiTrie(const class LangIDPackedMultiTrie *) ;
      ~LangIDMultiTrie() = default ;

      bool loadWords(const char *filename, uint32_t langID, bool verbose = false) ;

      // modifiers
      void setLanguage(uint32_t langID) { m_currentlangID = langID ; }
      bool insert(const uint8_t *key, unsigned keylength,
		  uint32_t langID, uint32_t frequency, bool stopgram) ;
      uint32_t increment(const uint8_t *key, unsigned keylength,
			 uint32_t langID, uint32_t incr = 1,
			 bool stopgram = false) ;
      bool incrementExtensions(const uint8_t *key, unsigned prevlength,
			       unsigned keylength, uint32_t langID,
			       uint32_t incr = 1) ;

      // accessors
      Node* node(NodeIndex N) const { return static_cast<Node*>(m_nodes.item(N)) ; }
      uint32_t currentLanguage() const { return m_currentlangID ; }
      Node* findNode(const uint8_t* key, unsigned keylength) const { return node(findKey(key,keylength)) ; }
      bool enumerate(uint8_t* keybuf, unsigned maxkeylength, EnumFn* fn, void *user_data) const ;
      bool enumerateChildren(NodeIndex nodeindex, uint8_t* keybuf, unsigned max_keylength_bits,
			     unsigned curr_keylength_bits, EnumFn* fn, void *user_data) const ;
      unsigned numExtensions(NodeIndex nodeindex, unsigned = 0) const ;
      bool allChildrenAreTerminals(NodeIndex nodeindex, unsigned = 0) const ;
      uint32_t countFreqRecords() const ;
      uint32_t numFullByteNodes() const ;
      uint32_t numTerminalNodes() const ;

      // I/O
      static LangIDMultiTrie *load(Fr::CFile& f) ;
      static LangIDMultiTrie *load(const char *filename) ;
      bool write(Fr::CFile& f) const ;
      bool write(const char *filename) const ;
      bool dump(Fr::CFile& f) const ;

   protected:
      size_t countTerminalNodes(NodeIndex nodeindex, unsigned keylen_bits = 0) const ;
      size_t countFullByteNodes(NodeIndex nodeindex, unsigned keylen_bits = 0) const ;

   private:
      void init(uint32_t capacity) ;
      bool writeHeader(Fr::CFile& f) const ;
   private:
      uint32_t		  m_currentlangID ;
   } ;

//----------------------------------------------------------------------

typedef TriePointer<LangIDMultiTrie> MultiTriePointer ;

#endif /* !__MTRIE_H_INCLUDED */

/* end of file trie.h */
