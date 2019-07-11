/****************************** -*- C++ -*- *****************************/
/*									*/
/*	LangIdent: n-gram based language-identification			*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File: mtrie.C - bit-slice-based Word-frequency multi-trie		*/
/*  Version:  1.30				       			*/
/*  LastEdit: 2019-07-07						*/
/*									*/
/*  (c) Copyright 2011,2012,2015,2019 Carnegie Mellon University	*/
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
#include "framepac/texttransforms.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define BUCKET_SIZE 65536	// must be power of 2

#define MULTITRIE_SIGNATURE "MulTrie\0"
#define MULTITRIE_FORMAT_VERSION 1

// reserve some space for future additions to the file format
#define MULTITRIE_PADBYTES_1  64

#if BITS_PER_LEVEL == 3
#  define LEVEL_SIZE 9
#else
#  define LEVEL_SIZE 8
#endif

/************************************************************************/
/*	Types								*/
/************************************************************************/

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

MultiTrieFrequency *MultiTrieFrequency::s_base_address = nullptr ;
uint32_t MultiTrieFrequency::s_max_alloc = 0 ;
uint32_t MultiTrieFrequency::s_curr_alloc = 0 ;

/************************************************************************/
/*	Helper functions						*/
/************************************************************************/

#ifndef lengthof
#  define lengthof(x) (sizeof(x)/sizeof((x)[0]))
#endif /* lengthof */

//----------------------------------------------------------------------

// defined in trie.C
uint32_t scaled_frequency(uint32_t raw_freq, uint64_t total_count) ;

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

MultiTrieFrequency* MultiTrieFrequency::allocate(uint32_t freq, uint32_t lang, bool stopgram)
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
	 return nullptr ;
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
      auto f = allocate(incr,ID) ;
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

MultiTrieFrequency *MultiTrieFrequency::read(Fr::CFile& f)
{
   if (f)
      {
      char buffer[sizeof(MultiTrieFrequency)] ;
      if (f.read(buffer,sizeof(buffer),sizeof(char)) == sizeof(buffer))
	 {
	 MultiTrieFrequency *freq = (MultiTrieFrequency*)buffer ;
	 return new MultiTrieFrequency(*freq) ;
	 }
      }
   return nullptr ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::readAll(Fr::CFile& f)
{
   Fr::UInt32 cnt ;
   if (f && f.readValue(&cnt))
      {
      uint32_t count = cnt.load() ;
      MultiTrieFrequency* base = nullptr ;
      if (f.readValues(&base,count))
	 {
	 Fr::Free(s_base_address) ;
	 s_base_address = base ;
	 s_max_alloc = s_curr_alloc = count ;
	 return true ;
	 }
      }
   return false ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::write(Fr::CFile& f) const
{
   return f && f.writeValue(*this) ;
}

//----------------------------------------------------------------------

bool MultiTrieFrequency::writeAll(Fr::CFile& f)
{
   Fr::UInt32 count(s_curr_alloc) ;
   return (f &&
           f.writeValue(count) &&
           f.write(baseAddress(),s_curr_alloc,sizeof(MultiTrieFrequency)) == s_curr_alloc) ;
}

/************************************************************************/
/*	Methods for class MultiTrieNode					*/
/************************************************************************/

MultiTrieNode::MultiTrieNode()
{
   m_frequency = LangIDMultiTrie::INVALID_FREQ ;
   return ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::isStopgram(unsigned langID) const
{
   auto freq = frequencies() ;
   for ( ; freq ; freq = freq->next())
      {
      if (freq->languageID() == langID)
	 return freq->isStopgram() ;
      }
   return false ;
}

//----------------------------------------------------------------------

unsigned MultiTrieNode::numFrequencies() const
{
   auto freq = frequencies() ;
   unsigned count = 0 ;
   for ( ; freq ; freq = freq->next())
      count++ ;
   return count ;
}

//----------------------------------------------------------------------

uint32_t MultiTrieNode::frequency(uint32_t langID) const
{
   auto freq = frequencies() ;
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
   auto f = MultiTrieFrequency::getAddress(m_frequency) ;
   if (f)
      f->setFrequency(langID,freq,stopgram) ;
   else
      {
      // insert the initial frequency record
      f = MultiTrieFrequency::allocate(freq,langID,stopgram) ;
      m_frequency = MultiTrieFrequency::getIndex(f) ;
      }
   return ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::setFrequencies(MultiTrieFrequency *freqs)
{
   if (freqs < MultiTrieFrequency::baseAddress())
      return false ;
   auto idx = MultiTrieFrequency::getIndex(freqs) ;
   if (idx == LangIDMultiTrie::INVALID_FREQ)
      return false ;
   m_frequency = idx ;
   return true ;
}

//----------------------------------------------------------------------

uint32_t MultiTrieNode::insertChild(unsigned int N, LangIDMultiTrie *trie)
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

bool MultiTrieNode::load(Fr::CFile& f)
{
   return f && f.readValue(this) ;
}

//----------------------------------------------------------------------

bool MultiTrieNode::write(Fr::CFile& f) const
{
   return f && f.writeValue(*this) ;
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

void LangIDMultiTrie::init(uint32_t cap)
{
   m_maxkeylen = 0 ;
   m_totaltokens = 0 ;
   m_ignorewhitespace = false ;
   if (cap == 0)
      cap = 1 ;
   m_nodes.reserve(cap) ;
   auto root = m_nodes.alloc() ;
   // initialize the root node
   new (node(root)) MultiTrieNode ;
   return ;
}

//----------------------------------------------------------------------
// NOTE: currently doesn't work with encodings that include NUL bytes in
//   their representation of characters other than NUL.

bool LangIDMultiTrie::loadWords(const char *filename, uint32_t langID, bool verbose)
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
	 insert((uint8_t*)lineptr,len,langID,freq,false) ;
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

bool LangIDMultiTrie::insert(const uint8_t *key, unsigned keylength,
		       uint32_t langID, uint32_t frequency, bool stopgram)
{
   auto leaf = node(insertKey(key,keylength)) ;
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
   auto n = node(cur_index) ;
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
      auto n = node(cur_index) ;
      if (!n)
	 return false ;
      uint32_t freq = n->frequency(langID) + incr ;
      n->setFrequency(langID,freq,n->isStopgram(langID)) ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerate(uint8_t *keybuf, unsigned maxkeylength, EnumFn *fn, void *user_data) const
{
   if (keybuf && fn)
      {
      memset(keybuf,'\0',maxkeylength) ;
      return enumerateChildren(ROOT_INDEX,keybuf,maxkeylength*8,0,fn,user_data) ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerateChildren(NodeIndex nodeindex,
				      uint8_t *keybuf,
				      unsigned max_keylength_bits,
				      unsigned curr_keylength_bits,
   				      LangIDMultiTrie::EnumFn *fn,
				      void *user_data) const
{
   auto n = node(nodeindex) ;
   assert(!(!n->leaf() && n->frequencies())) ;
   if (n->leaf() && !fn(n,keybuf,curr_keylength_bits/8,user_data))
      return false ;
   if (curr_keylength_bits < max_keylength_bits)
      {
      unsigned curr_bits = curr_keylength_bits + BITS_PER_LEVEL ;
      for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
	 {
	 uint32_t child = n->childIndex(i) ;
	 if (child != LangIDMultiTrie::NULL_INDEX)
	    {
	    auto childnode = node(child) ;
	    if (childnode)
	       {
	       unsigned byte = curr_keylength_bits / 8 ;
	       unsigned shift = LEVEL_SIZE - (curr_keylength_bits%8) - BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
	       if (shift == 0) curr_bits = curr_keylength_bits + BITS_PER_LEVEL - 1 ;
#endif
	       unsigned mask = (((1<<BITS_PER_LEVEL)-1) << shift) ;
	       keybuf[byte] &= ~mask ;
	       keybuf[byte] |= (i << shift) ;
	       if (!enumerateChildren(child,keybuf,max_keylength_bits,curr_bits,fn,user_data))
		  return false ;
	       }
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerateTerminalNodes(NodeIndex nodeindex, uint32_t &count, unsigned keylen_bits) const
{
   auto n = node(nodeindex)   ;
   if (!n->hasChildren())
      return true ;
   else if ((keylen_bits % 8) == 0 && allChildrenAreTerminals(nodeindex))
      {
      count += numExtensions(nodeindex) ;
      return true ;
      }
   keylen_bits = keylen_bits + BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
      {
      uint32_t child = n->childIndex(i) ;
      if (child != LangIDMultiTrie::NULL_INDEX && !enumerateTerminalNodes(child,count,keylen_bits))
	 return false ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::enumerateFullByteNodes(NodeIndex nodeindex, uint32_t &count, unsigned keylen_bits) const
{
   if (keylen_bits % 8 == 0)
      count++ ;
   keylen_bits = keylen_bits + BITS_PER_LEVEL ;
#if BITS_PER_LEVEL == 3
   if (keylen_bits % 8 == 1) keylen_bits-- ;
#endif
   auto n = node(nodeindex) ;
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; i++)
      {
      uint32_t child = n->childIndex(i) ;
      if (child != LangIDMultiTrie::NULL_INDEX)
	 {
	 if (!enumerateFullByteNodes(child,count,keylen_bits))
	    return false ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

unsigned LangIDMultiTrie::numExtensions(NodeIndex nodeindex, unsigned bits) const
{
   if (bits >= 8)
      {
      return 1 ;
      }
   auto n = node(nodeindex) ;
   unsigned count = 0 ;
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; ++i)
      {
      auto child = n->childIndex(i) ;
      if (child != NULL_INDEX)
	 count += numExtensions(child,bits+BITS_PER_LEVEL) ;
      }
   return count ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::allChildrenAreTerminals(NodeIndex nodeindex, unsigned bits) const
{
   auto n = node(nodeindex) ;
   if (bits >= 8)
      {
      return !n->hasChildren() ;
      }
   for (size_t i = 0 ; i < (1<<BITS_PER_LEVEL) ; ++i)
      {
      auto child = n->childIndex(i) ;
      if (child != NULL_INDEX && !allChildrenAreTerminals(child,bits+BITS_PER_LEVEL))
	 return false ;
      }
   return true ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::numFullByteNodes() const
{
   uint32_t count = 0 ;
   return enumerateFullByteNodes(ROOT_INDEX,count) ? count : 0 ;
}

//----------------------------------------------------------------------

uint32_t LangIDMultiTrie::numTerminalNodes() const
{
   uint32_t count = 0 ;
   return enumerateTerminalNodes(ROOT_INDEX,count) ? count : 0 ;
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

LangIDMultiTrie *LangIDMultiTrie::load(Fr::CFile& f)
{
   if (!f)
      return nullptr ;
   const size_t siglen = sizeof(MULTITRIE_SIGNATURE) ;
   char signature[siglen] ;
   if (f.read(signature,siglen,sizeof(char)) != siglen ||
      memcmp(signature,MULTITRIE_SIGNATURE,siglen) != 0)
      {
      // error: wrong file type
      return nullptr ;
      }
   unsigned char version ;
   if (!f.readValue(&version) || version != MULTITRIE_FORMAT_VERSION)
      {
      // error: wrong version of data file
      return nullptr ;
      }
   unsigned char bits ;
   if (!f.readValue(&bits) || bits != BITS_PER_LEVEL)
      {
      // error: wrong type of trie
      return nullptr ;
      }
   Fr::UInt32 val_used, val_tokens, val_keylen ;
   if (!f.readValue(&val_used) ||
      !f.readValue(&val_tokens) ||
      !f.readValue(&val_keylen))
      {
      // error reading header
      return nullptr ;
      }
   uint32_t used = val_used.load() ;
   auto trie = new LangIDMultiTrie(used) ;
   if (!trie)
      return nullptr ;
   trie->m_nodes.allocBatch(used) ;
   trie->addTokenCount(val_tokens.load()) ;
   trie->m_maxkeylen = (unsigned)val_keylen.load() ;
   // skip the padding (reserved bytes)
   f.seek(MULTITRIE_PADBYTES_1,SEEK_CUR) ;
   // read the actual trie nodes
   for (size_t i = 0 ; i < used ; i++)
      {
      auto node = trie->node(i) ;
      if (!node->load(f))
	 {
	 delete trie ;
	 return nullptr ;
	 }
      }
   // finally, read the frequency information
   if (!MultiTrieFrequency::readAll(f))
      {
      delete trie ;
      return nullptr ;
      }
   return trie ;
}

//----------------------------------------------------------------------

LangIDMultiTrie *LangIDMultiTrie::load(const char *filename)
{
   Fr::CInputFile fp(filename) ;
   return load(fp) ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::writeHeader(Fr::CFile& f) const
{
   // write the signature string
   const size_t siglen = sizeof(MULTITRIE_SIGNATURE) ;
   if (f.write(MULTITRIE_SIGNATURE,siglen,1) != siglen)
      return false; 
   // follow with the format version number
   unsigned char version = MULTITRIE_FORMAT_VERSION ;
   unsigned char bits = BITS_PER_LEVEL ;
   if (!f.writeValue(version) ||
      !f.writeValue(bits))
      return false ;
   // write out the size of the trie
   Fr::UInt32 val_used(size()), val_tokens(totalTokens()), val_keylen(longestKey()) ;
   if (!f.writeValue(val_used) ||
      !f.writeValue(val_tokens) ||
      !f.writeValue(val_keylen))
      return false ;
   // pad the header with NULs for the unused reserved portion of the header
   return f.putNulls(MULTITRIE_PADBYTES_1) ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::write(Fr::CFile& f) const
{
   if (!f || !writeHeader(f))
      return false ;
   // write the actual trie nodes
   for (size_t i = 0 ; i < size() ; i++)
      {
      const auto n = node(i) ;
      if (!n->write(f))
	 return false ;
      }
   // write the frequency information
   if (!MultiTrieFrequency::writeAll(f))
      return false ;
   f.writeComplete() ;
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::write(const char *filename) const
{
   Fr::COutputFile fp(filename,Fr::CFile::safe_rewrite) ;
   return this->write(fp) ? fp.close() : false ;
}

//----------------------------------------------------------------------

static const char hexdigit[] = "0123456789ABCDEF" ;

void write_escaped_char(uint8_t c, Fr::CFile& f)
{
   switch (c)
      {
      case '\\':
	 f.puts("\\\\") ;
	 break ;
      case ' ':
	 f.puts("\\ ") ;
	 break ;
      default:
	 if (c < ' ')
	    {
	    f.putc('\\') ;
	    f.putc(hexdigit[(c>>4)&0xF]) ;
	    f.putc(hexdigit[c&0xF]) ;
	    }
	 else
	    f.putc(c) ;
	 break ;
      }
   return ;
}

//----------------------------------------------------------------------

void write_escaped_key(Fr::CFile& f, const uint8_t* key, unsigned keylen)
{
   for (size_t i = 0 ; i < keylen ; i++)
      {
      write_escaped_char(key[i],f) ;
      }
   return ;
}

//----------------------------------------------------------------------

static bool dump_ngram(const MultiTrieNode *node, const uint8_t *key,
		       unsigned keylen, void *user_data)
{
   Fr::CFile& f = *((Fr::CFile*)user_data) ;
   if (f && node)
      {
      f.printf("   ") ;
      write_escaped_key(f,key,keylen) ;
      f.printf("  ::") ;
      auto freq = node->frequencies() ;
      for ( ; freq ; freq = freq->next())
	 {
	 f.printf(" %lu=%lu",(unsigned long)freq->languageID(),
		 (unsigned long)freq->frequency()) ;
	 }
      f.printf("\n") ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LangIDMultiTrie::dump(Fr::CFile& f) const
{
   Fr::LocalAlloc<uint8_t,10000> keybuf(longestKey()) ;
   return keybuf ? enumerate(keybuf,longestKey(),dump_ngram,&f) : false ;
}

// end of file mtrie.C //
