/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: pstrie.h - packed simple word-frequency trie			*/
/*  Version:  1.30				       			*/
/*  LastEdit: 27jun2019							*/
/*									*/
/*  (c) Copyright 2012,2013,2015,2019 Ralf Brown/CMU			*/
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

#ifndef __PSTRIE_H_INCLUDED
#define __PSTRIE_H_INCLUDED

#include "ptrie.h"
#include "trie.h"
#include "framepac/byteorder.h"

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

/************************************************************************/
/************************************************************************/

class LangIDPackedTrie ;

typedef bool PackedSimpleTrieEnumFn(const uint8_t *key, unsigned keylen,
				    uint32_t frequency, void *user_data) ;

typedef bool PackedSimpleTrieMatchFn(const uint8_t *key, unsigned keylen,
				     const LangIDPackedTrie *trie,
				     const class PackedSimpleTrieNode *node,
				     void *user_data) ;

//----------------------------------------------------------------------

class PackedSimpleTrieNode ;

class PackedTrieMatch
   {
   private:
      const PackedSimpleTrieNode *m_node ;
      uint8_t			 *m_key ;
      unsigned			  m_keylen ;
   public:
      PackedTrieMatch() { m_node = 0 ; m_key = 0 ; m_keylen = 0 ; }
      ~PackedTrieMatch() {}

      // accessors
      const PackedSimpleTrieNode *node() const { return m_node ; }
      const uint8_t *key() const { return m_key ; }
      unsigned keyLength() const { return m_keylen ; }

      // modifiers
      void setNode(const PackedSimpleTrieNode *n) { m_node = n ; }
      void setKeyBuffer(uint8_t *buffer, unsigned len)
	 { m_key = buffer ; m_keylen = len ; } 
      void setKey(uint8_t *newkey, unsigned len)
	 { if (key())
	       { if (len < m_keylen) m_keylen = len ;
	       for (unsigned i = 0 ; i < keyLength() ; i++)
		  m_key[i] = newkey[i] ;
	       }
	 }
   } ;

//----------------------------------------------------------------------

class PackedSimpleTrieNode
   {
   private:
      Fr::UInt32 m_frequency ;
      Fr::UInt32 m_firstchild ;
#define LENGTHOF_M_CHILDREN (PTRIE_CHILDREN_PER_NODE / sizeof(Fr::UInt32) / 8)
      Fr::UInt32 m_children[LENGTHOF_M_CHILDREN] ;
      uint8_t	 m_popcounts[LENGTHOF_M_CHILDREN] ;
#undef LENGTHOF_M_CHILDREN
   public:
      void *operator new(size_t, void *where) { return where ; }
      PackedSimpleTrieNode() ;
      ~PackedSimpleTrieNode() {}

      // accessors
      bool leaf() const
         { return m_frequency.load() != INVALID_FREQ ; }
      unsigned numChildren() const ;
      bool childPresent(unsigned int N) const ;
      uint32_t firstChild() const { return m_firstchild.load() ; }
      uint32_t childIndex(unsigned int N) const ;
      uint32_t childIndexIfPresent(unsigned int N) const ;
      PackedSimpleTrieNode *childNodeIfPresent(unsigned int N,
					       const LangIDPackedTrie *trie) const ;
      uint32_t frequency() const { return m_frequency.load() ; }
      unsigned countMatches(const LangIDPackedTrie *trie,
			    const uint8_t *key, unsigned keylen,
			    const WildcardSet **alternatives,
			    unsigned max_matches,
			    bool nonterminals_only = false) const ;
      bool enumerateChildren(const LangIDPackedTrie *trie,
			     uint8_t *keybuf, unsigned max_keylength_bits,
			     unsigned curr_keylength_bits,
			     PackedSimpleTrieEnumFn *fn,
			     void *user_data) const ;
      bool enumerateMatches(const LangIDPackedTrie *trie,
			    uint8_t *keybuf, unsigned max_keylength_bits,
			    unsigned curr_keylength_bits,
			    const WildcardSet **alternatives,
			    PackedSimpleTrieMatchFn *fn,
			    void *user_data) const ;
      unsigned enumerateMatches(const class EnumerationInfo *info,
				unsigned keylen) const ;
      bool nextFrequencies(const LangIDPackedTrie *trie,
			   uint32_t *frequencies) const ;
      bool addToScores(const LangIDPackedTrie *trie, float *scores,
		       double weight) const ;
      bool addToScores(const LangIDPackedTrie *trie, double *scores,
		       double weight) const ;

      // modifiers
      void setFirstChild(uint32_t index) { m_firstchild.store(index) ; }
      void setFrequency(uint32_t freq) { m_frequency.store(freq) ; }
      void setChild(unsigned N) ;
      void setPopCounts() ;
   } ;

//----------------------------------------------------------------------

class LangIDPackedTrie
   {
   public:
      static constexpr uint32_t NULL_INDEX = 0U ;
   public:
      LangIDPackedTrie() { init() ; }
      LangIDPackedTrie(const NybbleTrie *trie, uint32_t min_freq = 1,
		 bool show_conversion = true) ;
      LangIDPackedTrie(FILE *fp, const char *filename) ;
      ~LangIDPackedTrie() ;

      bool parseHeader(FILE *fp) ;

      // modifiers

      // accessors
      bool good() const
	 { return (m_nodes != 0) && m_size > 0 ; }
      bool terminalNode(const PackedSimpleTrieNode *n) const
         { return (n < m_nodes) || (n >= m_nodes + m_size) ; }
      uint32_t size() const { return m_size ; }
      unsigned longestKey() const { return m_maxkeylen ; }
      PackedSimpleTrieNode *node(uint32_t N) const
	 { if (N < m_size) return &m_nodes[N] ;
	   if ((N & PTRIE_TERMINAL_MASK) != 0)
	      {
	      uint32_t termindex = (N & ~PTRIE_TERMINAL_MASK) ;
	      if (termindex < m_numterminals)
		 return (PackedSimpleTrieNode*)&m_terminals[termindex] ; 
	      }
	   return 0 ;
	 }
      static bool isTerminalNode(uint32_t N)
	 { return (N & PTRIE_TERMINAL_MASK) != 0 ; }
      PackedSimpleTrieNode *getFullNode(uint32_t N) const
	 { return &m_nodes[N] ; }
      PackedSimpleTrieNode *getTerminalNode(uint32_t N) const
	 { return (PackedSimpleTrieNode*)&m_terminals[N & ~PTRIE_TERMINAL_MASK] ; }

      PackedSimpleTrieNode *findNode(const uint8_t *key, unsigned keylength) const ;
      uint32_t find(const uint8_t *key, unsigned keylength) const ;
      bool extendKey(uint32_t &nodeindex, uint8_t keybyte) const ;
      unsigned countMatches(const uint8_t *key, unsigned keylen,
			    const WildcardSet **alternatives,
			    unsigned max_matches,
			    bool nonterminals_only = false) const ;
      bool enumerate(uint8_t *keybuf, unsigned maxkeylength,
		     PackedSimpleTrieEnumFn *fn, void *user_data) const ;
      bool enumerate(uint8_t *key, unsigned keylength,
		     const WildcardSet **alternatives,
		     PackedSimpleTrieMatchFn *fn, void *user_data) const ;
      unsigned enumerate(uint8_t *key, unsigned keylength,
			 const WildcardSet **alternatives,
			 PackedTrieMatch *matches,
			 unsigned max_matches,
			 bool require_extensible_match) const ;

      // I/O
      static LangIDPackedTrie *load(FILE *fp, const char *filename) ;
      static LangIDPackedTrie *load(const char *filename) ;
      bool write(FILE *fp) const ;
      bool write(const char *filename) const ;
      bool dump(FILE *fp) const ;
   private:
      void init() ;
      bool writeHeader(FILE *fp) const ;
      uint32_t allocateChildNodes(unsigned numchildren) ;
      uint32_t allocateTerminalNodes(unsigned numchildren) ;
      bool insertChildren(PackedSimpleTrieNode *parent, const NybbleTrie *trie,
			  const NybbleTrieNode *node, uint32_t node_index,
			  unsigned keylen, uint32_t min_freq = 1) ;
      bool insertTerminals(PackedSimpleTrieNode *parent, const NybbleTrie *trie,
			   const NybbleTrieNode *node, uint32_t node_index,
			   unsigned keylen, uint32_t min_freq = 1) ;
   private:
      PackedSimpleTrieNode    *m_nodes ; // array of nodes
      PackedTrieTerminalNode *m_terminals ;
      Fr::MemMappedFile	*m_fmap ;	 // memory-map info
      uint32_t	 	 m_size ;	 // number of nodes in m_nodes
      uint32_t		 m_used ;	 // #nodes in use (temp during ctor)
      uint32_t		 m_numterminals ;
      uint32_t		 m_termused ;
      unsigned		 m_maxkeylen ;
      bool		 m_terminals_contiguous ;
   } ;

//----------------------------------------------------------------------

class PackedTriePointer
   {
   private:
//      static Fr::Allocator allocator ;
      LangIDPackedTrie	*m_trie ;
      uint32_t    	 m_nodeindex ;
      int	   	 m_keylength ;
      bool		 m_failed ;
   protected:
      void initPointer(LangIDPackedTrie *t)
      { m_trie = t ; m_nodeindex = 0 ; m_keylength = 0 ; m_failed = false ; }
   public:
//      void *operator new(size_t) { return allocator.allocate() ; }
//      void operator delete(void *blk) { allocator.release(blk) ; }
      PackedTriePointer() { initPointer(0) ; } // for arrays
      PackedTriePointer(LangIDPackedTrie *t) { initPointer(t) ; }
      PackedTriePointer(const LangIDPackedTrie *t)
	 { initPointer((LangIDPackedTrie*)t) ; }
      ~PackedTriePointer() { /*m_trie = 0 ; m_nodeindex = 0 ;*/ }

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
      bool hasExtension(uint8_t keybyte) const ;
      bool hasChildren(uint32_t node_index, uint8_t nybble) const ;
      int keyLength() const { return m_keylength ; }
      PackedSimpleTrieNode *node() const
         { return m_failed ? 0 : m_trie->node(m_nodeindex) ; }
   } ;

#endif /* !__PSTRIE_H_INCLUDED */

/* end of file pstrie.h */
