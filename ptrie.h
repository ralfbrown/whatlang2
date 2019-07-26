/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: ptrie.h - packed Word-frequency multi-trie			*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-26						*/
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
#include "trie.h"
#include "framepac/byteorder.h"
#include "framepac/file.h"
#include "framepac/memory.h"
#include "framepac/mmapfile.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define PTRIE_BITS_PER_LEVEL 8

/************************************************************************/
/************************************************************************/

class PackedTrieFreq
   {
   public:
      // define the bitfields used
      static constexpr uint32_t TRIE_LASTENTRY =      0x00002000 ;
      static constexpr uint32_t TRIE_STOPGRAM =       0x00004000 ;
      static constexpr uint32_t TRIE_LANGID_MASK =    0x00001FFF ;
      static constexpr uint32_t TRIE_FREQ_EXPONENT =  0x00018000 ;
      static constexpr uint32_t TRIE_FREQ_EXP_SHIFT = 15 ;
      static constexpr uint32_t TRIE_FREQ_MANTISSA =  0xFFFE0000 ;
      static constexpr uint32_t TRIE_FREQ_MAN_SHIFT = 17 ;
      static constexpr uint32_t TRIE_MANTISSA_LSB =   0x00020000 ;
      static constexpr uint32_t TRIE_FREQ_HIBITS =    0xC0000000 ; // top two bits
      static constexpr uint32_t EXPONENT_SCALE  = 2 ;  // each count in exp is two bits of shift

      // composite bitfields
      static constexpr uint32_t PACKED_TRIE_VALUE = (TRIE_FREQ_EXPONENT | TRIE_FREQ_MANTISSA | TRIE_STOPGRAM) ;
      static constexpr uint32_t TRIE_VALUE_SHIFT = (TRIE_FREQ_EXP_SHIFT - 1) ; // incl sign-bit
      static constexpr uint32_t TRIE_NUM_VALUES = (1UL << (32 - TRIE_VALUE_SHIFT)) ;

   public:
      void *operator new(size_t, void *where) { return where ; }
      PackedTrieFreq() = default ;
      PackedTrieFreq(uint32_t freq, uint32_t langID, bool last = true, bool is_stop = false) ;
      ~PackedTrieFreq() ;

      // accessors
      static constexpr unsigned maxExponent()
	 { return EXPONENT_SCALE * (TRIE_FREQ_EXPONENT >> TRIE_FREQ_EXP_SHIFT) ; }
      static constexpr double minWeight()
	 { return TRIE_FREQ_MANTISSA >> maxExponent() ; }
      static constexpr uint32_t maxLanguages() { return TRIE_LANGID_MASK+1 ; }
      static void quantize(uint32_t value, uint32_t &mantissa, uint32_t &exponent) ;
      static uint32_t mantissa(uint32_t scaled)
	 {
	 return (scaled & TRIE_FREQ_MANTISSA) >> TRIE_FREQ_MAN_SHIFT ;
	 }
      uint32_t mantissa() const { return mantissa(m_freqinfo.load()) ; }
      static uint32_t exponent(uint32_t scaled)
	 {
	 return (scaled & TRIE_FREQ_EXPONENT) >> TRIE_FREQ_EXP_SHIFT ;
	 }
      uint32_t exponent() const { return exponent(m_freqinfo.load()) ; }
      static uint32_t scaledScore(uint32_t data)
	 {
	 uint32_t mant = data & TRIE_FREQ_MANTISSA ;
	 uint32_t expon = (data & TRIE_FREQ_EXPONENT) >> (TRIE_FREQ_EXP_SHIFT - 1) ;
	 return (mant >> expon) ;
	 }
      uint32_t scaledScore() const { return scaledScore(m_freqinfo.load()) ; }
      double mappedScore() const
	 {
	 uint32_t data = (m_freqinfo.load() & PACKED_TRIE_VALUE) >> TRIE_VALUE_SHIFT ;
	 return s_value_map[data] ;
	 }

      double probability(uint32_t langID) const ;
      double probability() const  { return (scaledScore() / (100.0 * TRIE_SCALE_FACTOR)) ; }
      double percentage() const { return (scaledScore() / (1.0 * TRIE_SCALE_FACTOR)) ; }
      uint32_t languageID() const { return m_freqinfo.load() & TRIE_LANGID_MASK ; }
      bool isLast() const { return (m_freqinfo.load() & TRIE_LASTENTRY) != 0 ; }
      bool isStopgram() const { return (m_freqinfo.load() & TRIE_STOPGRAM) != 0 ; }
      const PackedTrieFreq *next() const { return isLast() ? nullptr : (this + 1) ; }
      static bool dataMappingInitialized() { return s_value_map_initialized ; }

      // manipulators
      void isLast(bool last) ;
      static void initDataMapping(double (*mapfn)(uint32_t)) ;
      static bool writeDataMapping(Fr::CFile& f) ;

   private:
      Fr::UInt32 m_freqinfo { TRIE_LASTENTRY } ;
      static double s_value_map[TRIE_NUM_VALUES] ;
      static bool s_value_map_initialized ;
   } ;

//----------------------------------------------------------------------

class PackedTrieTerminalNode
   {
   public:
      static constexpr uint32_t NULL_INDEX = 0U ;
      static constexpr uint32_t INVALID_FREQ = (uint32_t)~0 ;
   public:
      void *operator new(size_t, void *where) { return where ; }

      // accessors
      bool leaf() const { return m_frequency_info.load() != INVALID_FREQ ; }
      const PackedTrieFreq *frequencies(const PackedTrieFreq *base) const
         { return base + m_frequency_info.load() ; }

      // modifiers
      void reinit() { setFrequencies(INVALID_FREQ) ; }
      void setFrequencies(uint32_t index) { m_frequency_info.store(index) ; }
   protected:
      Fr::UInt32 m_frequency_info { INVALID_FREQ } ;
   } ;

//----------------------------------------------------------------------

class PackedTrieNode : public PackedTrieTerminalNode
   {
   public:
      static constexpr unsigned LENGTHOF_M_CHILDREN = (1<<PTRIE_BITS_PER_LEVEL) / (8*sizeof(Fr::UInt32)) ;
   public:
      void *operator new(size_t, void *where) { return where ; }

      // accessors
      bool childPresent(unsigned int N) const ;
      uint32_t firstChild() const { return m_firstchild.load() ; }
      uint32_t childIndex(unsigned int N) const ;
      uint32_t childIndexIfPresent(unsigned int N) const ;

      // modifiers
      void setFirstChild(uint32_t index) { m_firstchild.store(index) ; }
      void setChild(unsigned N) ;
      void setPopCounts() ;

   private:
      Fr::UInt32 m_firstchild { 0 } ;
      Fr::UInt32 m_children[LENGTHOF_M_CHILDREN] { 0 } ;
      uint8_t	 m_popcounts[LENGTHOF_M_CHILDREN] { 0 } ;
   } ;

//----------------------------------------------------------------------

class LangIDMultiTrie ;

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
      static constexpr uint32_t ROOT_INDEX = 0U ;
      static constexpr uint32_t INVALID_FREQ = (uint32_t)~0 ;
      typedef bool EnumFn(const PackedTrieNode *node, const uint8_t *key, unsigned keylen, void *user_data) ;

      // how do we distinguish non-terminal from terminal nodes?
      static constexpr uint32_t TERMINAL_MASK = 0x80000000 ;
   public:
      LangIDPackedMultiTrie() = default ;
      LangIDPackedMultiTrie(const LangIDMultiTrie *trie) ;
      LangIDPackedMultiTrie(Fr::CFile& f, const char *filename) ;
      ~LangIDPackedMultiTrie() = default ;

      bool parseHeader(Fr::CFile& f, size_t& numfull, size_t& numfreq, size_t& numterminals) ;

      // modifiers
      void ignoreWhiteSpace(bool ignore = true) { m_ignorewhitespace = ignore ; }
      void caseSensitivity(PTrieCase cs) { m_casesensitivity = cs ; }

      // accessors
      bool good() const { return size() > 0 && m_freq.size() ; }
      static bool terminalNode(uint32_t nodeindex) { return (nodeindex & TERMINAL_MASK) != 0 ; }
      uint32_t size() const { return m_nodes.size() ; }
      uint32_t numFrequencies() const { return m_freq.size(); }
      unsigned longestKey() const { return m_maxkeylen ; }
      bool ignoringWhiteSpace() const { return m_ignorewhitespace ; }
      PTrieCase caseSensitivity() const { return m_casesensitivity ; }
      const PackedTrieFreq* frequencyBaseAddress() const { return m_freq.item(0) ; }
      PackedTrieNode* node(uint32_t N) const
	 { if ((N & TERMINAL_MASK) != 0)
	      {
	      uint32_t termindex = (N & ~TERMINAL_MASK) ;
	      return (PackedTrieNode*)m_terminals.item(termindex) ; 
	      }
	   else
	      return m_nodes.item(N) ;
	 }
      PackedTrieNode *findNode(const uint8_t *key, unsigned keylength) const ;
      bool extendKey(uint32_t &nodeindex, uint8_t keybyte) const ;
      uint32_t extendKey(uint8_t keybyte, uint32_t nodeindex) const ;
      bool enumerate(uint8_t *keybuf, unsigned maxkeylength,
		     EnumFn *fn, void *user_data) const ;
      bool enumerateChildren(uint32_t nodeindex,
			     uint8_t *keybuf, unsigned max_keylength_bits,
			     unsigned curr_keylength_bits,
			     EnumFn *fn, void *user_data) const ;

      // I/O
      static Fr::Owned<LangIDPackedMultiTrie> load(Fr::CFile& f, const char *filename) ;
      static Fr::Owned<LangIDPackedMultiTrie> load(const char *filename) ;
      bool write(Fr::CFile& f) const ;
      bool write(const char *filename) const ;
      bool dump(Fr::CFile& f) const ;
   private:
      bool writeHeader(Fr::CFile& f) const ;
      uint32_t allocateChildNodes(unsigned numchildren) ;
      uint32_t allocateTerminalNodes(unsigned numchildren) ;
      bool insertChildren(PackedTrieNode *parent, const LangIDMultiTrie *mtrie,
			  uint32_t mnode_index, unsigned keylen = 0) ;
      bool insertTerminals(PackedTrieNode *parent, const LangIDMultiTrie *mtrie,
			   uint32_t mnode_index, unsigned keylen = 0) ;
   private:
      Fr::ItemPool<PackedTrieNode> m_nodes ;
      Fr::ItemPool<PackedTrieTerminalNode> m_terminals ;
      Fr::ItemPoolFlat<PackedTrieFreq> m_freq ;
      Fr::MemMappedFile	 m_fmap ;	 // memory-map info
      unsigned		 m_maxkeylen        { 0 } ;
      enum PTrieCase	 m_casesensitivity  { CS_Full } ;
      bool		 m_ignorewhitespace { false } ;
   } ;

//----------------------------------------------------------------------

typedef TriePointer<LangIDPackedMultiTrie> PackedMultiTriePointer ;

#endif /* !__PTRIE_H_INCLUDED */

/* end of file ptrie.h */
