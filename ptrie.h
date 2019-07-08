/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: ptrie.h - packed Word-frequency multi-trie			*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-07						*/
/*									*/
/*  (c) Copyright 2011,2012,2013,2015,2019 Carnegie Mellon University	*/
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

#ifndef __PTRIE_H_INCLUDED
#define __PTRIE_H_INCLUDED

#include <cstdio>
#include <limits.h>
#include <stdint.h>
#include "framepac/byteorder.h"
#include "framepac/file.h"
#include "framepac/memory.h"
#include "framepac/mmapfile.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define PTRIE_ROOT_INDEX 0
#define INVALID_FREQ ((uint32_t)~0)

#define PTRIE_BITS_PER_LEVEL 8
#define PTRIE_CHILDREN_PER_NODE 256

// how do we distinguish non-terminal from terminal nodes?
#define PTRIE_TERMINAL_MASK 0x80000000

// define the bitfields for PackedTrieFreq
#define PACKED_TRIE_LASTENTRY     0x00002000
#define PACKED_TRIE_STOPGRAM	  0x00004000
#define PACKED_TRIE_LANGID_MASK   0x00001FFF
#define PACKED_TRIE_FREQ_EXPONENT 0x00018000
#define PACKED_TRIE_FREQ_EXP_SHIFT 15
#define PACKED_TRIE_FREQ_MANTISSA 0xFFFE0000
#define PACKED_TRIE_FREQ_MAN_SHIFT 17
#define PACKED_TRIE_MANTISSA_LSB  0x00020000
#define PTRIE_EXPONENT_SCALE 2  // each count in exp is two bits of shift
#define PACKED_TRIE_FREQ_HIBITS   0xC0000000  // top two bits

#define PACKED_TRIE_VALUE (PACKED_TRIE_FREQ_EXPONENT | PACKED_TRIE_FREQ_MANTISSA | PACKED_TRIE_STOPGRAM)
#define PACKED_TRIE_VALUE_SHIFT (PACKED_TRIE_FREQ_EXP_SHIFT - 1) // incl sg-bit
#define PACKED_TRIE_NUM_VALUES (1UL << (32 - PACKED_TRIE_VALUE_SHIFT))

// we want to store percentages for entries in the trie in 32 bits.  Since
//   it is very unlikely that any ngram in the trie will have a probability
//   greater than 4.2%, scale the percentage by a factor of one billion.
#define TRIE_SCALE_FACTOR 1000000000L
#define SQRT_TRIE_SCALE_FACTOR (31622.7766016838)

/************************************************************************/
/************************************************************************/

class PackedTrieNode ;
typedef bool PackedTrieEnumFn(const PackedTrieNode *node, const uint8_t *key,
			      unsigned keylen, void *user_data) ;

class MultiTrieNode ;
class LangIDMultiTrie ;
class LangIDPackedMultiTrie ;

//----------------------------------------------------------------------

class PackedTrieFreq
   {
   private:
      Fr::UInt32 m_freqinfo ;
      static double s_value_map[PACKED_TRIE_NUM_VALUES] ;
      static bool s_value_map_initialized ;
   public:
      void *operator new(size_t, void *where) { return where ; }
      PackedTrieFreq()
	 { m_freqinfo.store(PACKED_TRIE_LASTENTRY) ; }
      PackedTrieFreq(uint32_t freq, uint32_t langID, bool last = true,
		     bool is_stop = false) ;
      ~PackedTrieFreq() ;

      // accessors
      static void quantize(uint32_t value, uint32_t &mantissa,
			   uint32_t &exponent) ;
      static uint32_t mantissa(uint32_t scaled)
	 {
	 return (scaled & PACKED_TRIE_FREQ_MANTISSA) >> PACKED_TRIE_FREQ_MAN_SHIFT ;
	 }
      uint32_t mantissa() const { return mantissa(m_freqinfo.load()) ; }
      static uint32_t exponent(uint32_t scaled)
	 {
	 return (scaled & PACKED_TRIE_FREQ_EXPONENT) >> PACKED_TRIE_FREQ_EXP_SHIFT ;
	 }
      uint32_t exponent() const { return exponent(m_freqinfo.load()) ; }
      static uint32_t scaledScore(uint32_t data)
	 {
	 uint32_t mant = data & PACKED_TRIE_FREQ_MANTISSA ;
	 uint32_t expon = data & PACKED_TRIE_FREQ_EXPONENT ;
	 expon >>= (PACKED_TRIE_FREQ_EXP_SHIFT - 1) ;
	 return (mant >> expon) ;
	 }
      uint32_t scaledScore() const
	 {
	 uint32_t data =m_freqinfo.load() ;
	 uint32_t mant = data & PACKED_TRIE_FREQ_MANTISSA ;
	 uint32_t expon = data & PACKED_TRIE_FREQ_EXPONENT ;
	 expon >>= (PACKED_TRIE_FREQ_EXP_SHIFT - 1) ;
	 return (mant >> expon) ;
	 }
      double mappedScore() const
	 {
	 uint32_t data = m_freqinfo.load() & PACKED_TRIE_VALUE ;
	 data >>= PACKED_TRIE_VALUE_SHIFT ;
	 return s_value_map[data] ;
	 }

      double probability(uint32_t langID) const ;
      double probability() const 
	 { return (scaledScore() / (100.0 * TRIE_SCALE_FACTOR)) ; }
      double percentage() const
	 { return (scaledScore() / (1.0 * TRIE_SCALE_FACTOR)) ; }
      uint32_t languageID() const
         { return m_freqinfo.load() & PACKED_TRIE_LANGID_MASK ; }
      bool isLast() const
	 { return (m_freqinfo.load() & PACKED_TRIE_LASTENTRY) != 0 ; }
      bool isStopgram() const
	 { return (m_freqinfo.load() & PACKED_TRIE_STOPGRAM) != 0 ; }
      const PackedTrieFreq *next() const { return isLast() ? 0 : (this + 1) ; }
      static bool dataMappingInitialized() { return s_value_map_initialized ; }

      // manipulators
      void isLast(bool last) ;
      static void initDataMapping(double (*mapfn)(uint32_t)) ;
      static bool writeDataMapping(Fr::CFile& f) ;
   } ;

//----------------------------------------------------------------------

class PackedTrieNode
   {
   private:
      Fr::UInt32 m_frequency_info ;
      Fr::UInt32 m_firstchild ;
#define LENGTHOF_M_CHILDREN (PTRIE_CHILDREN_PER_NODE / sizeof(Fr::UInt32) / 8)
      Fr::UInt32 m_children[LENGTHOF_M_CHILDREN] ;
      uint8_t	 m_popcounts[LENGTHOF_M_CHILDREN] ;
#undef LENGTHOF_M_CHILDREN
   public:
      void *operator new(size_t, void *where) { return where ; }
      PackedTrieNode() ;
      ~PackedTrieNode() {}

      // accessors
      bool leaf() const
         { return m_frequency_info.load() != INVALID_FREQ ; }
      bool childPresent(unsigned int N) const ;
      uint32_t firstChild() const { return m_firstchild.load() ; }
      uint32_t childIndex(unsigned int N) const ;
      uint32_t childIndexIfPresent(uint8_t N) const ;
      uint32_t childIndexIfPresent(unsigned int N) const ;
      double probability(const PackedTrieFreq *base, uint32_t ID = 0) const ;
      const PackedTrieFreq *frequencies(const PackedTrieFreq *base) const
         { return base + m_frequency_info.load() ; }
      bool enumerateChildren(const LangIDPackedMultiTrie *trie,
			     uint8_t *keybuf, unsigned max_keylength_bits,
			     unsigned curr_keylength_bits,
			     PackedTrieEnumFn *fn,
			     void *user_data) const ;

      // modifiers
      void setFirstChild(uint32_t index) { m_firstchild.store(index) ; }
      void setFrequencies(uint32_t index) { m_frequency_info.store(index) ; }
      void setChild(unsigned N) ;
      void setPopCounts() ;
   } ;

//----------------------------------------------------------------------

class PackedTrieTerminalNode
   {
   public:
      static constexpr uint32_t NULL_INDEX = 0U ;
   private:
      Fr::UInt32 m_frequency_info ;
   public:
      void *operator new(size_t, void *where) { return where ; }
      PackedTrieTerminalNode() { m_frequency_info.store(INVALID_FREQ); }
      ~PackedTrieTerminalNode() {}

      // accessors
      bool leaf() const { return true ; }
      bool childPresent(unsigned int /*N*/) const { return false ; }
      uint32_t firstChild() const { return 0 ; }
      uint32_t childIndex(unsigned int /*N*/) const { return NULL_INDEX ; }
      uint32_t childIndexIfPresent(unsigned int /*N*/) const { return NULL_INDEX ; }
      double probability(const PackedTrieFreq *base, uint32_t ID = 0) const ;
      const PackedTrieFreq *frequencies(const PackedTrieFreq *base) const
         { return base + m_frequency_info.load() ; }

      // modifiers
      void setFrequencies(uint32_t index)
	 { m_frequency_info.store(index) ; }
   } ;

//----------------------------------------------------------------------

enum PTrieCase
   {
      CS_Full = 0,
      CS_ASCII,
      CS_Latin1
   } ;

class LangIDPackedMultiTrie // : public Fr::PackedMultiTrie<...>
   {
   public:
      static constexpr uint32_t NULL_INDEX = 0U ;
   public:
      LangIDPackedMultiTrie() { init() ; }
      LangIDPackedMultiTrie(const LangIDMultiTrie *trie) ;
      LangIDPackedMultiTrie(Fr::CFile& f, const char *filename) ;
      ~LangIDPackedMultiTrie() ;

      bool parseHeader(Fr::CFile& f) ;

      // modifiers
      void ignoreWhiteSpace(bool ignore = true) { m_ignorewhitespace = ignore ; }
      void caseSensitivity(PTrieCase cs) { m_casesensitivity = cs ; }

      // accessors
      bool good() const
	 { return (m_nodes != 0) && (m_freq != 0) && m_size > 0 ; }
      bool terminalNode(const PackedTrieNode *n) const
	 { return (n < m_nodes) || (n >= m_nodes + m_size) ; }
      uint32_t size() const { return m_size ; }
      uint32_t numFrequencies() const { return m_numfreq ; }
      unsigned longestKey() const { return m_maxkeylen ; }
      bool ignoringWhiteSpace() const { return m_ignorewhitespace ; }
      PTrieCase caseSensitivity() const { return m_casesensitivity ; }
      const PackedTrieFreq *frequencyBaseAddress() const { return m_freq ; }
      PackedTrieNode *node(uint32_t N) const
	 { if (N < m_size) return &m_nodes[N] ;
	   if ((N & PTRIE_TERMINAL_MASK) != 0)
	      {
	      uint32_t termindex = (N & ~PTRIE_TERMINAL_MASK) ;
	      if (termindex < m_numterminals)
		 return (PackedTrieNode*)&m_terminals[termindex] ; 
	      }
	   return nullptr ;
	 }
      PackedTrieNode *findNode(const uint8_t *key, unsigned keylength) const ;
      bool extendKey(uint32_t &nodeindex, uint8_t keybyte) const ;
      uint32_t extendKey(uint8_t keybyte, uint32_t nodeindex) const ;
      bool enumerate(uint8_t *keybuf, unsigned maxkeylength,
		     PackedTrieEnumFn *fn, void *user_data) const ;

      // I/O
      static LangIDPackedMultiTrie *load(Fr::CFile& f, const char *filename) ;
      static LangIDPackedMultiTrie *load(const char *filename) ;
      bool write(Fr::CFile& f) const ;
      bool write(const char *filename) const ;
      bool dump(Fr::CFile& f) const ;
   private:
      void init() ;
      bool writeHeader(Fr::CFile& f) const ;
      uint32_t allocateChildNodes(unsigned numchildren) ;
      uint32_t allocateTerminalNodes(unsigned numchildren) ;
      bool insertChildren(PackedTrieNode *parent, const LangIDMultiTrie *mtrie,
			  const MultiTrieNode *mnode, uint32_t mnode_index,
			  unsigned keylen = 0) ;
      bool insertTerminals(PackedTrieNode *parent, const LangIDMultiTrie *mtrie,
			   const MultiTrieNode *mnode, uint32_t mnode_index,
			   unsigned keylen = 0) ;
   private:
      PackedTrieNode    *m_nodes ;	 // array of nodes
      PackedTrieTerminalNode *m_terminals ;
      PackedTrieFreq    *m_freq ;	 // array of frequency records
      Fr::MemMappedFile	*m_fmap ;	 // memory-map info
      uint32_t	 	 m_size ;	 // number of nodes in m_nodes
      uint32_t		 m_numterminals ;
      uint32_t		 m_numfreq ;	 // number of records in m_freq
      uint32_t		 m_used ;	 // #nodes in use (temp during ctor)
      uint32_t		 m_termused ;
      uint32_t		 m_freqused ;	 // #freq recs in use (ctor temp)
      unsigned		 m_maxkeylen ;
      enum PTrieCase	 m_casesensitivity ;
      bool		 m_ignorewhitespace ;
      bool		 m_terminals_contiguous ;
   } ;

//----------------------------------------------------------------------

class PackedMultiTriePointer
   {
   private:
//      static Fr::Allocator allocator ;
      LangIDPackedMultiTrie	*m_trie ;
      uint32_t    	 m_nodeindex ;
      int	   	 m_keylength ;
      bool		 m_failed ;
   protected:
      void initPointer(LangIDPackedMultiTrie *t)
	 { m_trie = t ; m_nodeindex = 0 ; m_failed = false ; }
   public:
//      void *operator new(size_t) { return allocator.allocate() ; }
//      void operator delete(void *blk) { allocator.release(blk) ; }
      PackedMultiTriePointer() { initPointer(0) ; } // for arrays
      PackedMultiTriePointer(LangIDPackedMultiTrie *t) { initPointer(t) ; }
      PackedMultiTriePointer(const LangIDPackedMultiTrie *t)
	 { initPointer(const_cast<LangIDPackedMultiTrie*>(t)) ; }
      ~PackedMultiTriePointer() { /*m_trie = 0 ; m_nodeindex = 0 ;*/ }

      void resetKey()
	 {
	 m_nodeindex = PTRIE_ROOT_INDEX ;
	 m_keylength = 0 ;
	 m_failed = false ;
	 return ;
	 }
      bool extendKey(uint8_t keybyte) ;

      // accessors
      bool lookupSuccessful() const ;
      bool hasChildren(uint32_t node_index, uint8_t nybble) const ;
      int keyLength() const { return m_keylength ; }
      PackedTrieNode *node() const
         { return m_failed ? 0 : m_trie->node(m_nodeindex) ; }
   } ;

#endif /* !__PTRIE_H_INCLUDED */

/* end of file ptrie.h */
