/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     mklangid.C	build language-id model database	*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-10						*/
/*                                                                      */
/*  (c) Copyright 2010,2011,2012,2013,2014,2015,2019			*/
/*		 Ralf Brown/Carnegie Mellon University			*/
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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <iostream>
#include <iomanip>
#include "langid.h"
#include "prepfile.h"
#include "trie.h"
#include "mtrie.h"
#include "ptrie.h"
#include "framepac/bitvector.h"
#include "framepac/init.h"
#include "framepac/message.h"
#include "framepac/texttransforms.h"

using namespace std ;
using namespace Fr ;

/************************************************************************/
/*	Configuration and Manifest Constants				*/
/************************************************************************/

#define VERSION "1.30"

#define MAX_NGRAMS 5000

//#define DEBUG		// generate torrents of output in verbose mode

// since we always count at least trigrams, we don't support lengths less
//   than three
#define ABSOLUTE_MIN_LENGTH 3
#define DEFAULT_MAX_LENGTH  8
#define ABSOLUTE_MAX_LENGTH 500

#define MAX_OVERSAMPLE 2.5
#define MAX_INCREMENT 10

// if n-gram plus an affix appears at least this many times the count
//   of times the bare n-gram appears, ignore the bare n-gram
#define AFFIX_RATIO 0.90
#define MINLEN_AFFIX_RATIO 0.995        // be strict for minimum-length ngrams
// limit how low the user can set the affix ratio
#define MIN_AFFIX_RATIO 0.4

// assume that each character in a listed n-gram implies N characters
//   in the total training data when TotalCount: is not present; this discount
//   factor reduces false positive detections from mguesser ngram lists
#define ASSUMED_NGRAM_DENSITY 8

// how much to discount faked ngrams inserted via the UTF8: directive in
//   a frequency list file
#define FAKED_NGRAM_DISCOUNT 10

#define DEFAULT_SIMILARITY_THRESHOLD 0.50

// factor times minimum representable prob at which to cut off stopgrams
#define STOPGRAM_CUTOFF 2

// weighting of stopgrams based on amount of training data
#define TRAINSIZE_NO_WEIGHT 15000
#define TRAINSIZE_FULL_WEIGHT 2000000

// the amount by which to increase the weight of an n-gram which is unique
//   to a particular model among closely-related languages
#define UNIQUE_BOOST 1.0

// precompute smoothing before scaling to reduce quantization error when
//   packing the trie
#define SMOOTHING_POWER 0.14

/************************************************************************/
/************************************************************************/

#ifndef lengthof
#  define lengthof(x) (sizeof(x)/sizeof((x)[0]))
#endif /* lengthof */

#if BUFFER_SIZE < 2*FrMAX_LINE
#  undef BUFFER_SIZE
#  define BUFRER_SIZE (2*FrMAX_LINE)
#endif

/************************************************************************/
/*	Types for this module						*/
/************************************************************************/

struct NgramEnumerationData
   {
   public:
      NybbleTrie *m_ngrams ;
      uint32_t   *m_frequencies ;
      bool       &m_have_max_length ;
      bool	  m_inserted_ngram ;
      unsigned    m_min_length ;
      unsigned    m_max_length ;
      unsigned	  m_desired_length ;
      unsigned    m_topK ;
      unsigned    m_count ;
      unsigned	  m_alignment ;
      uint32_t    m_min_freq ;
   public:
      NgramEnumerationData(bool &have_max_length)
	 : m_have_max_length(have_max_length), m_inserted_ngram(false), m_alignment(1)
	 {}
      ~NgramEnumerationData() {}

   } ;

//----------------------------------------------------------------------

class StopGramInfo
   {
   private:
      NybbleTrie 	   *m_trie ;
      NybbleTrie	   *m_currngrams ;
      NybbleTrie	   *m_ngramweights ;
      const PackedTrieFreq *m_freqbase ;
      const LanguageScores *m_weights ;
      const BitVector  *m_selected ;
      unsigned		    m_activelang ;
   public:
      StopGramInfo(NybbleTrie *t, NybbleTrie *t2, NybbleTrie *t3,
		   const PackedTrieFreq *base,
	           const LanguageScores *wt, const BitVector *sel,
		   unsigned lang)
	 { m_trie = t ; m_currngrams = t2 ; m_ngramweights = t3 ;
	   m_freqbase = base ; m_weights = wt ; m_selected = sel ;
	   m_activelang = lang ; }
      ~StopGramInfo() {}

      // accessors
      NybbleTrie *trie() const { return m_trie ; }
      NybbleTrie *currLangTrie() const { return m_currngrams ; }
      NybbleTrie *weightTrie() const { return m_ngramweights ; }
      unsigned activeLanguage() const { return m_activelang ; }
      bool selected(size_t langid) const
	 { return m_selected ? m_selected->getBit(langid) : false ; }
      const PackedTrieFreq *freqBaseAddress() const { return m_freqbase ; }
      double weight(size_t langid) const
	 { return m_weights ? m_weights->score(langid) : 1.0 ; }
   } ;

//----------------------------------------------------------------------

class StopGramWeight
   {
   private:
      const NybbleTrie *m_weighttrie ;
      uint64_t          m_totalbytes ;
      bool	        m_scaled ;
   public:
      StopGramWeight(const NybbleTrie *wt, uint64_t tb, bool sc)
	 { m_weighttrie = wt ; m_totalbytes = tb ; m_scaled = sc ; }
      ~StopGramWeight() {}

      // accessors
      bool scaled() const { return m_scaled ; }
      uint64_t totalBytes() const { return m_totalbytes ; }
      const NybbleTrie *weights() const { return m_weighttrie ; }
      uint32_t weight(const uint8_t *key, unsigned keylen) const ;

      // modifiers

   } ;

//----------------------------------------------------------------------

typedef bool FileReaderFunc(PreprocessedInputFile *, va_list) ;

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

static bool verbose = false ;
static bool store_similarities = false ;
static bool do_dump_trie = false ;
static bool crubadan_format = false ;
static BigramExtension bigram_extension = BigramExt_None ;
static unsigned topK = MAX_NGRAMS ;
static unsigned minimum_length = ABSOLUTE_MIN_LENGTH ;
static unsigned maximum_length = DEFAULT_MAX_LENGTH ;
static unsigned alignment = 1 ;
static const char *vocabulary_file = nullptr ;
static double max_oversample = MAX_OVERSAMPLE ;
static double affix_ratio = AFFIX_RATIO ;
static double discount_factor = 1.0 ;
static LanguageIdentifier *language_identifier = nullptr ;
static bool skip_numbers = false ;
static bool subsample_input = false ;
static uint64_t byte_limit = ~0 ;
static double unique_boost = UNIQUE_BOOST ;
static double smoothing_power = SMOOTHING_POWER ;
static double log_smoothing_power = 1.0 ;

// multiple of min proportion in confusible models for an ngram to be
//   added to baseline model; 0 = disable the addition
static double confusibility_thresh = 0.0 ;

/************************************************************************/
/*	Helper functions						*/
/************************************************************************/

// in trigram.C
void insert_frequency(uint32_t newelt, uint32_t *heap, size_t heaplen) ;
uint32_t adjusted_threshold(const uint32_t *frequencies) ;

//----------------------------------------------------------------------

static void usage(const char *argv0, const char *bad_arg)
{
   if (bad_arg)
      cerr << "Unrecognized argument " << bad_arg << endl << endl ;
   cerr << "MKLANGID version " VERSION "  Copyright 2011,2012 Ralf Brown/CMU -- GNU GPLv3" << endl ;
   cerr << "Usage: " << argv0 << " [=DBFILE] {options} file ... [{options} file ...]"
	<< endl ;
   cerr << "  Specify =DBFILE to use language database DBFILE instead of the\n" ;
   cerr << "  default " DEFAULT_LANGID_DATABASE << "; with ==DBFILE, the database\n" ;
   cerr << "  will not be updated (use -w to store results)" << endl ; ;
   cerr << "Options:" << endl ;
   cerr << "   -h       show this usage summary" << endl ;
   cerr << "   -l LANG  specify language of following files (use ISO-639 two-letter code)" << endl ;
   cerr << "   -r REG   specify regional variant of the language (use ISO two-letter" << endl ;
   cerr << "            country codes, e.g. for locale 'en_US', use -l en -r US" << endl ;
   cerr << "            [optional])" << endl ;
   cerr << "   -e ENC   specify the character encoding, e.g. iso8859-1, utf-8, etc." << endl ;
   cerr << "   -s SRC   specify the source of the training data (optional)" << endl ;
   cerr << "   -W SCR   specify writing system (script) of the training data (optional)" << endl ;
   cerr << "   -kK      collect top K n-grams by frequency (default " << topK << ")\n" ;
   cerr << "   -mN      require n-grams to consist of at least N bytes (min 3)\n" ;
   cerr << "   -MN      limit n-grams to at most N bytes (default " << DEFAULT_MAX_LENGTH << ", max " << ABSOLUTE_MAX_LENGTH << ")\n" ;
   cerr << "   -i       ignore blanks when processing files\n" ;
   cerr << "   -n       skip ngrams containing newlines in following files\n" ;
   cerr << "   -nn      skip ngrams starting with digits as well\n" ;
   cerr << "   -LN      limit training to first N bytes of input\n" ;
   cerr << "   -L@N     limit training to N bytes uniformly sampled from input\n" ;
   cerr << "   -b       omit bigram table from model for following files\n" ;
   cerr << "   -ON      set maximum oversampling factor to N (default " << max_oversample << ")" << endl ;
   cerr << "   -aX      set affix ratio; remove 'ABC' if c(ABCD) >= X * c(ABC)" << endl ;
   cerr << "   -dX      set probability discount factor for X" << endl ;
   cerr << "   -R SPEC  compute stop-grams relative to related language(s) listed in SPEC" << endl ;
   cerr << "   -B BOOST increase smoothed scores of n-grams unique to model by BOOST*" << endl ;
   cerr << "   -S SMTH  set smoothing power to SMTH (negative for logarithmic)" << endl;
   cerr << "   -1       convert Latin-1 input to UTF-8" << endl ;
   cerr << "   -2b      pad input bytes to 16 bits (big-endian)" << endl ;
   cerr << "   -2l      pad input bytes to 16 bits (little-endian)" << endl ;
   cerr << "   -2-      don't pad input bytes to 16 bits" << endl ;
   cerr << "   -8b      convert UTF8 input to UTF-16 (big-endian)" << endl ;
   cerr << "   -8l      convert UTF8 input to UTF-16 (little-endian)" << endl ;
   cerr << "   -8-      don't convert UTF8" << endl ;
   cerr << "   -AN      alignment: only start ngram at multiple of N (1,2,4)" << endl ;
   cerr << "   -f       following files are frequency lists (count then string)" << endl ;
   cerr << "   -fc      following files are frequency lists (count/string, word delim)" << endl ;
   cerr << "   -ft      following files are frequency lists (string/tab/count)" << endl ;
   cerr << "   -v       run verbosely" << endl ;
   cerr << "   -wFILE   write resulting vocabulary list to FILE in plain text" << endl ;
   cerr << "   -D       dump computed multi-trie to standard output" << endl ;
   cerr << "Notes:" << endl ;
   cerr << "\tThe -1 -b -f -i -n -nn -R -w flags reset after each group of files." << endl;
   cerr << "\t-2 and -8 are mutually exclusive -- the last one specified is used." << endl ;
   exit(1) ;
}

//----------------------------------------------------------------------

static const char *get_arg(int &argc, const char **&argv)
{
   if (argv[1][2])
      return argv[1]+2 ;
   else
      {
      argc-- ;
      argv++ ;
      return argv[1] ;
      }
}

//----------------------------------------------------------------------

static void print_quoted_char(CFile& f, uint8_t ch)
{
   switch (ch)
      {
      case '\0': f.puts("\\0") ;		break ;
      case '\f': f.puts("\\f") ;		break ;
      case '\n': f.puts("\\n") ;		break ;
      case '\r': f.puts("\\r") ;		break ;
      case '\t': f.puts("\\t") ;		break ;
      case ' ':  f.puts("\\ ") ;		break ;
      case '\\': f.puts("\\\\") ;		break ;
      default:	 f.putc(ch) ;			break ;
      }
   return ;
}

//----------------------------------------------------------------------

static int UCS2_to_UTF8(unsigned long codepoint, char *buf)
{
   if (codepoint < 0x80)
      {
      // encode as single byte
      *buf = (char)codepoint ;
      return 1 ;
      }
   else if (codepoint > 0x10FFFF)
      return -1 ;
   if (codepoint < 0x800)
      {
      // encode in two bytes
      buf[0] = (unsigned char)(0xC0 | ((codepoint & 0x07C0) >> 6)) ;
      buf[1] = (unsigned char)(0x80 | (codepoint & 0x003F)) ;
      return 2 ;
      }
   else if (codepoint < 0x010000)
      {
      // encode as three bytes
      buf[0] = (unsigned char)(0xE0 | ((codepoint & 0xF000) >> 12)) ;
      buf[1] = (unsigned char)(0x80 | ((codepoint & 0x0FC0) >> 6)) ;
      buf[2] = (unsigned char)(0x80 | (codepoint & 0x003F)) ;
      return 3 ;
      }
   else
      {
      // encode as four bytes
      buf[0] = (unsigned char)(0xF0 | ((codepoint & 0x1C0000) >> 18)) ;
      buf[1] = (unsigned char)(0x80 | ((codepoint & 0x03F000) >> 12)) ;
      buf[2] = (unsigned char)(0x80 | ((codepoint & 0x000FC0) >> 6)) ;
      buf[3] = (unsigned char)(0x80 | (codepoint & 0x00003F)) ;
      return 4 ;
      }
}

//----------------------------------------------------------------------

static uint64_t read_files(const char **filelist, unsigned num_files,
			   bool show_error, FileReaderFunc *reader, ...)
{
   if (!filelist)
      return 0 ;
   va_list args ;
   va_start(args,reader) ;
   uint64_t total_bytes = 0 ;
   bool OK = true ;
   for (size_t i = 0 ; i < num_files && total_bytes < byte_limit && OK ; i++)
      {
      const char *filename = filelist[i] ;
      if (filename && *filename)
	 {
	 PreprocessedInputFile infile(filename,byte_limit - total_bytes, subsample_input) ;
	 if (infile.good())
	    {
	    cout << "  Processing " << filename << endl ;
	    va_list argcopy ;
	    va_copy(argcopy,args) ;
	    if (!reader(&infile,argcopy))
	       OK = false ;
	    }
	 else if (show_error)
	    {
	    cerr << "Error opening '" << filename << "' for reading" << endl ;
	    }
	 total_bytes += infile.bytesRead() ;
	 infile.close() ;
	 }
      }
   va_end(args) ;
   return total_bytes ;
}

/************************************************************************/
/*	Language-name manipulation functions				*/
/************************************************************************/

static int compare_langcode(const LanguageID *id1, const LanguageID *id2)
{
   if (!id1)
      return id2 ? +1 : 0 ;
   else if (!id2)
      return -1 ;
   const char *l1 = id1->language() ;
   const char *l2 = id2->language() ;
   if (!l1)
      return l2 ? +1 : 0 ;
   else if (!l2)
      return -1 ;
   return strcmp(l1,l2) ;
}

//----------------------------------------------------------------------

static int compare_codepair(const LanguageID *id1, const LanguageID *id2)
{
   if (!id1)
      return id2 ? +1 : 0 ;
   else if (!id2)
      return -1 ;
   const char *l1 = id1->language() ;
   const char *l2 = id2->language() ;
   if (!l1)
      return l2 ? +1 : 0 ;
   else if (!l2)
      return -1 ;
   int cmp = strcmp(l1,l2) ;
   if (cmp)
      return cmp ;
   const char *enc1 = id1->encoding() ;
   const char *enc2 = id2->encoding() ;
   if (!enc1)
      return enc2 ? +1 : 0 ;
   else if (!enc2)
      return -1 ;
   return strcmp(enc1,enc2) ;
}

/************************************************************************/
/*	Methods for class StopGramWeight				*/
/************************************************************************/

uint32_t StopGramWeight::weight(const uint8_t *key, unsigned keylen) const
{
   uint32_t raw = weights()->find(key,keylen) ;
   if (!scaled() && raw > 0)
      {
      double percent = raw / (double)TRIE_SCALE_FACTOR ;
      double proportion = percent / 100.0 ;
      uint32_t unscaled = (uint32_t)(proportion * totalBytes() + 0.5) ;
      return unscaled ;
      }
   else
      return raw ;
}

/************************************************************************/
/*	Stopgram Computation						*/
/************************************************************************/

static bool collect_language_ngrams_packed(const PackedTrieNode *node,
					   const uint8_t *key,
					   unsigned keylen, void *user_data)
{
   auto stop_gram_info = reinterpret_cast<StopGramInfo*>(user_data) ;
   assert(stop_gram_info != nullptr) ;
   assert(node != nullptr) ;
   if (node->leaf())
      {
      auto freq = node->frequencies(stop_gram_info->freqBaseAddress()) ;
      // the following is the smallest value which can be represented by
      //   the bitfields in PackedTrieFreq; anything less than that will
      //   round to zero unless smoothed
      unsigned max_exp = PTRIE_EXPONENT_SCALE * (PACKED_TRIE_FREQ_EXPONENT >> PACKED_TRIE_FREQ_EXP_SHIFT) ;
      double minweight = PACKED_TRIE_FREQ_MANTISSA >> max_exp ;
      // since extra stopgrams have a cost both in storage and runtime,
      //   we want to omit any that have little possible impact on the
      //   overall score, so bump up the above absolute minimum threshold
      minweight *= STOPGRAM_CUTOFF ;
      for ( ; ; freq++)
	 {
	 unsigned langid = freq->languageID() ;
	 if (stop_gram_info->activeLanguage() == langid)
	    {
	    uint32_t wt = freq->scaledScore() ;
	    stop_gram_info->currLangTrie()->insert(key,keylen,wt,false) ;
	    }
	 if (stop_gram_info->selected(langid) &&
	     !freq->isStopgram() && freq->percentage() > 0.0)
	    {
	    double wt = stop_gram_info->weight(langid) * freq->scaledScore() ;
	    uint32_t weight = (uint32_t)(wt + 0.5) ;
	    if (weight >= minweight)
	       {
	       // if the weight is large enough to be bothered with,
	       //   add the ngram to the presence trie and its computed
	       //   weight to the weight trie
	       stop_gram_info->trie()->insert(key,keylen,0,false) ;
	       stop_gram_info->weightTrie()->insertMax(key,keylen,weight,
						       false) ;
	       }
	    }
	 if (freq->isLast())
	    break ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

static bool select_models_by_name(const char *languages, BitVector &selected)
{
   auto descriptions = dup_string(languages) ;
   char *desc = *descriptions ;
   bool did_select = false ;
   while (desc && *desc)
      {
      char *desc_end = strchr(desc,',') ;
      if (desc_end)
	 *desc_end++ = '\0' ;
      else
	 desc_end = strchr(desc,'\0') ;
      unsigned langnum = language_identifier->languageNumber(desc) ;
      if (langnum != (unsigned)~0)
	 {
	 selected.setBit(langnum,true) ;
	 did_select = true ;
	 }
      else
	 {
	 cerr << "Warning: no match for language descriptor " << desc << endl ;
	 }
      desc = desc_end ;
      }
   return did_select ;
}

//----------------------------------------------------------------------

static bool select_models_by_similarity(size_t langid, BitVector &selected,
					LanguageScores *weights,
					const char *thresh)
{
   bool did_select = false ;
   if (!weights)
      {
      cout << "Unable to compute cross-language similarities, will not "
	      "compute stop-grams." << endl ;
      return did_select ;
      }
   char *endptr = nullptr ;
   double threshold = strtod(thresh,&endptr) ;
   if (threshold <= 0.0 || threshold > 1.0)
      threshold = DEFAULT_SIMILARITY_THRESHOLD ;
   const LanguageID *curr = language_identifier->languageInfo(langid) ;
   // figure out which, if any, language models are close enough to the
   //   current one to be the basis for stopgrams
   for (size_t langnum = 0 ; langnum < weights->numLanguages() ; langnum++)
      {
      if (langnum == langid || weights->score(langnum) < threshold)
	 continue ;
      const LanguageID *other = language_identifier->languageInfo(langnum) ;
      if (other)
	 {
	 // check that the other model isn't the same language,
	 //   region, AND encoding
	 if (curr->language() && other->language() &&
	     strcmp(curr->language(),other->language()) == 0 &&
	     curr->region() && other->region() &&
	     strcmp(curr->region(),other->region()) == 0 &&
	     curr->encoding() && other->encoding() &&
	     strcmp(curr->encoding(),other->encoding()) == 0)
	    continue ;
	 }
      did_select = true ;
      selected.setBit(langnum,true) ;
      if (/*verbose &&*/ other)
	 {
	 cout << "  similarity to "<< other->language() << "_"
	      << other->region() << "-" << other->encoding() << " is "
	      << weights->score(langnum) << endl ;
	 }
      }
   return did_select ;
}

//----------------------------------------------------------------------

static NybbleTrie *load_stop_grams_selected(unsigned langid,
					    LanguageScores *weights,
   					    LangIDPackedMultiTrie *ptrie,
					    const BitVector *selected,
					    NybbleTrie *&curr_ngrams,
					    NybbleTrie *&ngram_weights)
{
   // further discount the selected other languages by the absolute
   //   amount of training data in the primary language (since less
   //   data means a greater chance that the n-gram is not seen purely
   //   due to data sparsity)
   const LanguageID *curr = language_identifier->languageInfo(langid) ;
   if (!curr)
      {
      return new NybbleTrie ;
      }
   uint64_t train = curr->trainingBytes() ;
   if (train < TRAINSIZE_FULL_WEIGHT)
      {
      train = (train > TRAINSIZE_NO_WEIGHT) ? train - TRAINSIZE_NO_WEIGHT : 0 ;
      double scalefactor = (TRAINSIZE_FULL_WEIGHT - TRAINSIZE_NO_WEIGHT) ;
      double scale = ::pow(train / scalefactor,0.7) ;
      weights->scaleScores(scale) ;
      }
   // because the values stored in the tries have already been
   //   adjusted to smooth them, we need to apply the same
   //   adjustment to the inter-model weights to avoid overly
   //   reducing the weights of stopgrams
   for (size_t i = 0 ; i < weights->numLanguages() ; i++)
      {
      if (selected->getBit(i))
	 {
	 double sc = weights->score(i) ;
	 sc = ::pow(sc,smoothing_power) ;
#if 0 //!!!
	 if (i != langid)
	    {
	    // adjust the strength of the stop-grams for each model
	    //   pair by the coverage fractions of the two models
	    double cover = curr->coverageFactor() ;
	    if (cover > 0.0)
	       sc /= cover ;
	    // assume that the coverages are independent of each
	    //   other, which lets us simply multiple the two
	    //   coverage fractions
	    const LanguageID *other = language_identifier->languageInfo(i) ;
	    cover = other->coverageFactor() ;
	    cerr<<"adj="<<cover<<endl;
	    if (cover > 0.0)
	       sc /= cover ;
	    }
#endif
	 weights->setScore(i,sc) ;
	 }
      }
   uint8_t key[1000] ;
   unsigned maxkey = ptrie->longestKey() ;
   if (maxkey > sizeof(key))
      maxkey = sizeof(key) ;
   auto stop_grams = new NybbleTrie ;
   curr_ngrams = new NybbleTrie ;
   ngram_weights = new NybbleTrie ;
   auto freq_base = ptrie->frequencyBaseAddress() ;
   StopGramInfo stop_gram_info(stop_grams,curr_ngrams,ngram_weights,
			       freq_base,weights,selected,langid) ;
   ptrie->enumerate(key,maxkey,collect_language_ngrams_packed,
		    &stop_gram_info) ;
   return stop_grams ;
}

//----------------------------------------------------------------------

static NybbleTrie *load_stop_grams(const LanguageID *lang_info,
				   const char *languages,
				   NybbleTrie *&curr_ngrams,
				   NybbleTrie *&ngram_weights,
				   uint64_t &training_bytes)
{
   curr_ngrams = nullptr ;
   ngram_weights = nullptr ;
   training_bytes = 0 ;
   if (!languages)
      return nullptr ;
   auto ptrie = language_identifier->trie() ;
   if (!ptrie)
      {
      return nullptr ;
      }
   unsigned langid = language_identifier->languageNumber(lang_info) ;
   training_bytes = language_identifier->trainingBytes(langid) ;
   cout << "Computing similarities relative to "
	<< lang_info->language() << "_" << lang_info->region()
	<< "-" << lang_info->encoding() << endl ;
   Owned<LanguageScores> weights = language_identifier->similarity(langid) ;
   ScopedObject<BitVector> selected(language_identifier->numLanguages()) ;
   bool selected_models ;
   if (languages && *languages == '@')
      {
      selected_models = select_models_by_similarity(langid,*selected,weights,languages+1) ;
      }
   else
      {
      selected_models = select_models_by_name(languages,*selected) ;
      }
   NybbleTrie *stop_grams ;
   if (selected_models || 1) //!!! we need to run regardless, to create curr_ngrams
      {
      delete curr_ngrams ;
      delete ngram_weights ;
      stop_grams = load_stop_grams_selected(langid,weights,ptrie,&selected,curr_ngrams,ngram_weights) ;
      }
   else
      {
      stop_grams = new NybbleTrie ;
      }
   return stop_grams ;
}

//----------------------------------------------------------------------
// count the ngrams in the current training data which match ngrams in
//   the highly-similar models selected earlier

static void accumulate_confusible_ngrams(PreprocessedInputFile *infile,
					 NybbleTrie *confusible,
					 NybbleTriePointer *states)
{
   auto maxkey = confusible->longestKey() ;
   while (infile->moreData())
      {
      int keybyte = infile->getByte() ;
      if (keybyte == EOF)
	 break ;
      states[maxkey].invalidate() ;
      states[0].resetKey() ;
      for (auto i = maxkey ; i > 0 ; i--)
	 {
	 if (!states[i-1])
	    continue ;
	 if (states[i-1].extendKey((uint8_t)(keybyte&0xFF)))
	    {
	    // check whether we're at a leaf node; if so, increment its frequency
	    auto node = states[i-1].node() ;
	    if (node && node->leaf())
	       {
	       node->incrFrequency() ;
	       }
	    states[i] = states[i-1] ;
	    }
	 else
	    {
	    states[i].invalidate() ;
	    }
	 states[i-1].invalidate() ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static bool accumulate_confusible_ngrams(PreprocessedInputFile *infile, va_list args)
					 
{
   auto confusible = va_arg(args,NybbleTrie*) ;
   if (!infile || !infile->good() || !confusible)
      return false ;
   unsigned maxkey = confusible->longestKey() ;
   NewPtr<NybbleTriePointer> states(maxkey+2) ;
   if (states)
      {
      for (size_t i = 0 ; i < maxkey+2 ; ++i)
	 states[i].setTrie(confusible) ;
      accumulate_confusible_ngrams(infile,confusible,states.begin()) ;
      }
   return true ;
}

//----------------------------------------------------------------------

static bool add_stop_gram(const NybbleTrie* trie, uint32_t nodeindex,
			  const uint8_t *key, unsigned keylen,
			  void *user_data)
{
   auto node = trie->node(nodeindex) ;
   assert(node != nullptr) ;
   if (keylen <= 2 || !node->leaf())
      return true ;
   if (node->frequency() == 0 || node->isStopgram())
      {
      auto ngrams = reinterpret_cast<NybbleTrie*>(user_data) ;
      auto weights = reinterpret_cast<StopGramWeight*>(ngrams->userData()) ;
      auto weight = weights->weight(key,keylen) ;
      ngrams->insert(key,keylen,weight,true) ;
      }
   else if (confusibility_thresh > 0.0)
      {
      // optionally add the ngram as a regular ngram if the proportion
      //   in the current training data is higher than the given
      //   multiple of the lowest proportion in the confusible models
//FIXME
      // ngrams->insert(key,keylen,weight,false) ;
      }
   return true ;
}

//----------------------------------------------------------------------

static bool boost_unique_ngram(const NybbleTrie* trie, uint32_t nodeindex,
			       const uint8_t *key, unsigned keylen,
			       void *user_data)
{
   auto stop_grams = reinterpret_cast<NybbleTrie*>(user_data) ;
   auto n = trie->node(nodeindex) ;
   if (n && n->leaf() && n->frequency() > 0 && !n->isStopgram())
      {
      auto sgnode = stop_grams->findNode(key,keylen) ;
      if (!sgnode || !sgnode->leaf())
	 {
	 auto freq = n->frequency() ;
	 // boost the weight of this node
	 auto boost = scale_frequency(unique_boost,smoothing_power,
					log_smoothing_power) ;
	 auto boosted = (uint32_t)(freq * boost + 0.9) ;
	 if (boosted < n->frequency()) // did we roll over?
	    boosted = 0xFFFFFFFF ;
	 n->setFrequency(boosted) ;
	 }
      else if (0)
	 {
	 auto weights = reinterpret_cast<StopGramWeight*>(stop_grams->userData()) ;
	 auto freq = n->frequency() ;
	 auto sgfreq = scaled_frequency(sgnode->frequency(),weights->totalBytes(),
			       smoothing_power,log_smoothing_power) ;
	 if (freq > 2 * sgfreq)
	    {
	    cout << "alternate" << endl ;
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

static bool add_stop_grams(const char **filelist, unsigned num_files,
			   NybbleTrie *ngrams, NybbleTrie *stop_grams,
			   const NybbleTrie *ngram_weights, bool scaled)
{
   if (!stop_grams || stop_grams->size() <= 100)
      return true ;
   cout << "Computing Stop-Grams" << endl ;
   // accumulate counts for all the ngrams in the stop-gram list
   uint64_t total_bytes
      = read_files(filelist,num_files,false,&accumulate_confusible_ngrams,stop_grams) ;
   StopGramWeight stop_gram_weight(ngram_weights,total_bytes,scaled) ;
   ngrams->setUserData(&stop_gram_weight) ;
   stop_grams->setUserData(&stop_gram_weight) ;
   // 'stop_grams' contains the union of all n-grams in the models for
   //   closely-related languages.  Any n-grams in the current language's
   //   model which don't appear in any of the other models get a boost,
   //   because they are particularly strong evidence that the text is in
   //   the given language
   unsigned longkey = stop_grams->longestKey() ;
   if (ngrams->longestKey() > longkey)
      longkey = ngrams->longestKey() ;
   LocalAlloc<uint8_t> ngram_buf(longkey+1) ;
   if (unique_boost > 1.0)
      {
      ngrams->enumerate(&ngram_buf,ngrams->longestKey(),boost_unique_ngram,stop_grams) ;
      }
   // scan the stop-gram list and add any with a zero count to the main
   //   n-gram list, flagged as stop-grams; optionally, add any which
   //   have higher counts in the current language but are not in the
   //   baseline model to the model as well
   stop_grams->enumerate(ngram_buf,stop_grams->longestKey(),add_stop_gram,ngrams) ;
   return true ;
}

/************************************************************************/
/*	Model output							*/
/************************************************************************/

static unsigned count_languages(const LanguageIdentifier *id,
				int (*cmp)(const LanguageID*,const LanguageID*))
{
   unsigned count = 0 ;
   if (id)
      {
      size_t num_langs = id->numLanguages() ;
      LocalAlloc<LanguageID*> langcodes(num_langs) ;
      for (size_t i = 0 ; i < num_langs ; i++)
	 {
	 langcodes[i] = (LanguageID*)id->languageInfo(i) ;
	 }
      // sort the codes to get runs of equal codes
      std::sort(&langcodes,&langcodes+num_langs,cmp) ;
      count = 1 ;
      for (size_t i = 1 ; i < num_langs ; i++)
	 {
	 // increment the count every time the code differs from the
	 //   one before
	 if (cmp(langcodes[i-1],langcodes[i]) != 0)
	    count++ ;
	 }
      }
   return count ;
}

//----------------------------------------------------------------------

static bool save_database(const char *database_file)
{
   if (!database_file || !*database_file)
      database_file = DEFAULT_LANGID_DATABASE ;
   if (language_identifier->numLanguages() > 0)
      {
      if (do_dump_trie)
	 {
	 CFile f(stdout) ;
	 f.printf("=======================\n") ;
	 language_identifier->dump(f,verbose) ;
	 }
      unsigned num_languages = count_languages(language_identifier,compare_langcode) ;
      unsigned num_pairs = count_languages(language_identifier,compare_codepair) ;
      cout << "Database contains " << language_identifier->numLanguages()
	   << " models, " << num_languages
	   << " distinct language codes,\n\tand " << num_pairs
	   << " language/encoding pairs" << endl ;
      cout << "Saving database" << endl ;
      return language_identifier->write(database_file) ;
      }
   return false ;
}

//----------------------------------------------------------------------

static bool dump_ngrams(const NybbleTrie* trie, uint32_t nodeindex, const uint8_t *key,
			unsigned keylen, void *user_data)
{
   CFile& f = *((CFile*)user_data) ;
   auto node = trie->node(nodeindex) ;
   if (node->leaf())
      {
      uint32_t freq = node->frequency() ;
      if (node->isStopgram() && freq > 0)
	 f.putc('-') ;
      f.printf("%u\t",freq) ;
      for (size_t i = 0 ; i < keylen ; i++)
	 {
	 print_quoted_char(f,key[i]) ;
	 }
      f.printf("\n") ;
      }
   return true ;
}

//----------------------------------------------------------------------

// this global variable makes dump_vocabulary() non-reentrant
static uint64_t dump_total_bytes ;

static bool dump_ngrams_scaled(const NybbleTrie* trie, uint32_t nodeindex,
			       const uint8_t *key,
			       unsigned keylen, void *user_data)
{
   CFile& f = *((CFile*)user_data) ;
   auto node = trie->node(nodeindex) ;
   if (node->leaf())
      {
      uint32_t freq = node->frequency() ;
      if (node->isStopgram() && freq > 0)
	 f.putc('-') ;
      double unscaled = (unscale_frequency(freq,smoothing_power) * dump_total_bytes / 100.0) + 0.99 ;
      f.printf("%u\t",(unsigned)unscaled) ;
      for (size_t i = 0 ; i < keylen ; i++)
	 {
	 print_quoted_char(f,key[i]) ;
	 }
      f.printf("\n") ;
      }
   return true ;
}

//----------------------------------------------------------------------

static void dump_vocabulary(const NybbleTrie *ngrams, bool scaled,
			    const char *vocab_file, unsigned max_length,
			    uint64_t total_bytes, const LanguageID &opts)
{
   COutputFile f(vocab_file) ;
   if (!f)
      {
      cerr << "Unable to open '" << vocab_file << "' to write vocabulary"
	   << endl ;
      return ;
      }
   if (total_bytes)
      f.printf("TotalCount: %llu\n", (unsigned long long)total_bytes) ;
   f.printf("Lang: %s",opts.language()) ;
   if (opts.friendlyName() != opts.language())
      f.printf("=%s",opts.friendlyName()) ;
   f.printf("\nScript: %s\nRegion: %s\nEncoding: %s\nSource: %s\n",
      opts.script(),opts.region(),opts.encoding(),opts.source()) ;
   if (alignment > 1)
      f.printf("Alignment: %d\n",alignment) ;
   if (discount_factor > 1.0)
      f.printf("Discount: %g\n",discount_factor) ;
   if (ngrams->ignoringWhiteSpace())
      f.printf("IgnoreBlanks: yes\n") ;
   if (opts.coverageFactor() > 0.0 && opts.coverageFactor() != 1.0)
      f.printf("Coverage: %g\n",opts.coverageFactor()) ;
   if (opts.countedCoverage() > 0.0 && opts.countedCoverage() != 1.0)
      f.printf("WeightedCoverage: %g\n",opts.countedCoverage()) ;
   if (opts.freqCoverage() > 0.0)
      f.printf("FreqCoverage: %g\n",opts.freqCoverage()) ;
   if (opts.matchFactor() > 0.0)
      f.printf("MatchFactor: %g\n",opts.matchFactor()) ;
   LocalAlloc<uint8_t> keybuf(max_length+1) ;
   if (scaled)
      {
      dump_total_bytes = total_bytes ;
      (void)ngrams->enumerate(keybuf,max_length,dump_ngrams_scaled,&f) ;
      }
   else
      (void)ngrams->enumerate(keybuf,max_length,dump_ngrams,&f) ;
   return ;
}

/************************************************************************/
/************************************************************************/

static unsigned set_oversampling(unsigned top_K, unsigned abs_min_len,
				 unsigned min_len, bool aligned)
{
   if (abs_min_len < min_len)
      {
      double base = aligned ? 2.0 : 1.0 ;
      double oversample = ::pow(2,base+(min_len-abs_min_len)/5.0) ;
      if (oversample > max_oversample)
	 oversample = max_oversample ;
      return (unsigned)(top_K * oversample) ;
      }
   else
      return top_K ;
}

//----------------------------------------------------------------------

static bool count_raw_trigrams(PreprocessedInputFile *infile, va_list args)
{
   auto counts = va_arg(args,TrigramCounts*) ;
   // prime the trigrams with the first two bytes
   uint8_t c1 = infile->getByte() ;
   uint8_t c2 = infile->getByte() ;
   unsigned offset = 0 ;
   // now process the rest of the file, incrementing the count
   //   for every trigram encountered
   while (infile->moreData())
      {
      int c3 = infile->getByte() ;
      if (c3 == EOF)
	 break ;
      if (offset % alignment == 0)
	 counts->incr(c1,c2,(uint8_t)c3) ;
      c1 = c2 ;
      c2 = (uint8_t)c3 ;
      offset++ ;
      }
   return true ;
}

//----------------------------------------------------------------------

static uint64_t count_trigrams(const char **filelist, unsigned num_files,
			       TrigramCounts &counts, bool skip_newlines,
			       bool aligned, BigramCounts **bigrams)
{
   cout << "Counting trigrams" << endl ;
   uint64_t total_bytes = read_files(filelist,num_files,true,&count_raw_trigrams,&counts) ;
   if (bigrams)
      (*bigrams) = new BigramCounts(counts) ;
   if (bigram_extension == BigramExt_ASCIILittleEndian ||
       bigram_extension == BigramExt_UTF8LittleEndian)
      {
      // don't bother with trigrams starting with a NUL, as those will
      //  split a 16-bit character
      for (size_t i = 0 ; i < 256 ; i++)
	 {
	 for (size_t j = 0 ; j < 256 ; j++)
	    {
	    counts.clear('\0',i,j) ;
	    }
	 }
      if (skip_newlines)
	 {
	 counts.clear(' ','\0',' ') ;
	 for (size_t i = 0 ; i < 256 ; i++)
	    {
	    counts.clear(i,'\0','\r') ;
	    counts.clear(i,'\0','\n') ;
	    counts.clear('\r','\0',i) ;
	    counts.clear('\n','\0',i) ;
	    counts.clear('\t','\0',i) ;
	    }
	 }
      if (skip_numbers)
	 {
	 for (size_t c1 = '0' ; c1 <= '9' ; c1++)
	    {
	    counts.clear('.','\0',c1) ;
	    counts.clear(',','\0',c1) ;
	    counts.clear(c1,'\0','.') ;
	    counts.clear(c1,'\0',',') ;
	    for (size_t c3 = '0' ; c3 <= '9' ; c3++)
	       {
	       counts.clear(c1,'\0',c3) ;
	       }
	    }
	 // also skip certain repeated punctuation marks which are not informative
	 counts.clear('-','\0','-') ;
	 counts.clear('=','\0','=') ;
	 counts.clear('*','\0','*') ;
	 counts.clear('.','\0','.') ;
	 counts.clear('?','\0','?') ;
	 }
      }
   else if (bigram_extension == BigramExt_ASCIIBigEndian ||
	    bigram_extension == BigramExt_UTF8BigEndian)
      {
      if (bigram_extension == BigramExt_ASCIIBigEndian)
	 {
	 // force trigrams to start with a NUL, as otherwise we will
	 //  split a 16-bit character
	 for (size_t i = 1 ; i < 256 ; i++)
	    {
	    for (size_t j = 0 ; j < 256 ; j++)
	       {
	       counts.clear(i,'\0',j) ;
	       }
	    }
	 }
      if (skip_newlines)
	 {
	 for (size_t i = 0 ; i < 256 ; i++)
	    {
	    counts.clear(i,'\0','\r') ;
	    counts.clear(i,'\0','\n') ;
	    counts.clear('\0','\r',i) ;
	    counts.clear('\0','\n',i) ;
	    counts.clear('\0','\t',i) ;
	    }
	 }
      }
   else if (bigram_extension == BigramExt_None && !aligned)
      {
      if (skip_newlines)
	 {
	 for (size_t i = 0 ; i < 256 ; i++)
	    {
	    counts.clear(' ',' ',i) ;
	    for (size_t j = 0 ; j < 256 ; j++)
	       {
	       counts.clear(i,j,'\r') ;
	       counts.clear(i,j,'\n') ;
	       counts.clear(i,'\r',j) ;
	       counts.clear(i,'\n',j) ;
	       counts.clear('\r',i,j) ;
	       counts.clear('\n',i,j) ;
	       counts.clear('\t',i,j) ;
	       if (alignment == 1)
		  {
		  counts.clear(i,j,'\0') ;
		  counts.clear(i,'\0',j) ;
		  counts.clear('\0',i,j) ;
		  }
	       }
	    }
	 }
      if (skip_numbers)
	 {
	 for (size_t c1 = '0' ; c1 <= '9' ; c1++)
	    {
	    for (size_t c2 = '0' ; c2 <= '9' ; c2++)
	       {
	       counts.clear('.',c1,c2) ;
	       counts.clear(',',c1,c2) ;
	       for (size_t c3 = 0 ; c3 <= 0xFF ; c3++)
		  {
		  counts.clear(c1,c2,c3) ;
		  if (c3 >= '0' && c3 <= '9')
		     {
		     counts.clear(c1,'.',c3) ;
		     counts.clear(c1,',',c3) ;
		     }
		  }
	       }
	    }
	 // also skip certain repeated punctuation marks which are not informative
	 counts.clear('-','-','-') ;
	 counts.clear('=','=','=') ;
	 counts.clear('*','*','*') ;
	 counts.clear('.','.','.') ;
	 counts.clear('?','?','?') ;
	 }
      }
   cout << "  Processed " << total_bytes << " bytes" << endl ;
   return total_bytes ;
}

//----------------------------------------------------------------------

static bool count_ngrams(PreprocessedInputFile *infile, va_list args)
{
   auto ngrams = va_arg(args,NybbleTrie*) ;
   auto min_length = va_arg(args,unsigned) ;
   auto max_length = va_arg(args,unsigned) ;
   bool skip_newlines = (bool)va_arg(args,int) ;
   bool aligned = (bool)va_arg(args,int) ;
   if (max_length < min_length || max_length == 0)
      return false ;
   LocalAlloc<uint8_t> ngram(max_length) ;
   // fill the ngram buffer, except for the last byte
   for (size_t i = 0 ; i+1 < max_length ; i++)
      {
      int c = infile->getByte() ;
      if (c == EOF)
	 return false ;
      ngram[i] = (uint8_t)c ;
      }
   // now iterate through the file, counting the ngrams
   unsigned offset = 0 ;
   while (infile->moreData())
      {
      int c = infile->getByte() ;
      if (c == EOF)
	 break ;
      // fill the last byte of the buffer
      ngram[max_length-1] = (uint8_t)c ;
      // increment n-gram counts if they are an extension of a known n-gram,
      //   but don't include newlines if told not to do so
      size_t max_len = max_length ;
      if (skip_newlines)
         {
	 if (bigram_extension == BigramExt_ASCIIBigEndian ||
	     bigram_extension == BigramExt_UTF8BigEndian)
	    {
	    for (size_t i = (min_length - 1)/2 ; i < (max_length/2) ; i++)
	       {
	       if (ngram[2*i] == '\0' &&
		   (ngram[2*i+1] == '\n' || ngram[2*i+1] == '\r' || ngram[2*i+1] == '\0'))
		  {
		  max_len = 2*i ;
		  break ;
		  }
	       }
	    }
	 else if (bigram_extension == BigramExt_ASCIILittleEndian ||
		  bigram_extension == BigramExt_UTF8LittleEndian)
	    {
	    for (size_t i = (min_length - 1)/2 ; i < (max_length/2) ; i++)
	       {
	       if (ngram[2*i+1] == '\0' &&
		   (ngram[2*i] == '\n' || ngram[2*i] == '\r' || ngram[2*i] == '\0'))
		  {
		  max_len = 2*i ;
		  break ;
		  }
	       }
	    }
	 else
	    {
	    for (size_t i = min_length - 1 ; i < max_length ; i++)
	       {
	       if (ngram[i] == '\n' || ngram[i] == '\r' ||
		   (!aligned && bigram_extension == BigramExt_None && ngram[i] == '\0'))
		  {
		  max_len = i ;
		  break ;  
		  }
	       }
            }
         }
      if (alignment == 2 && min_length > 3 && ngram[0] == '\0' && ngram[2] == '\0')
	 {
	 // check for big-endian two-character sequences we weren't
	 //   able to filter at the trigram stage
	 if (skip_newlines)
	    {
	    if (ngram[1] == ' ' && ngram[3] == ' ')
	       max_len = 0 ;
	    }
	 if (skip_numbers)
	    {
	    if (isdigit(ngram[3]) &&
		(isdigit(ngram[1]) || ngram[1] == '.' || ngram[1] == ','))
	       max_len = 0 ;
	    else if (isdigit(ngram[1]) &&
		     (ngram[3] == '.' || ngram[3] == ','))
	       max_len = 0 ;
	    }
	 if ((ngram[1] == '-' && ngram[3] == '-') ||
	     (ngram[1] == '=' && ngram[3] == '=') ||
	     (ngram[1] == '*' && ngram[3] == '*') ||
	     (ngram[1] == '.' && ngram[3] == '.') ||
	     (ngram[1] == '?' && ngram[3] == '?'))
	    max_len = 0 ;
	 }
      if (max_len >= min_length && (offset % alignment) == 0)
	 ngrams->incrementExtensions(ngram,min_length-1,max_len) ;
      // shift the buffer by one byte
      memmove(ngram,ngram+1,max_length-1) ;
      }
   return true ;
}

//----------------------------------------------------------------------

static bool remove_suffix(const NybbleTrie* trie, uint32_t nodeindex, const uint8_t *key, unsigned keylen,
			  void *user_data)
{
   auto enum_data = reinterpret_cast<NgramEnumerationData*>(user_data) ;
   auto align = enum_data->m_alignment ;
   if (keylen == enum_data->m_desired_length && keylen >= align + enum_data->m_min_length)
      {
      auto node = trie->node(nodeindex) ;
      auto suffix = trie->findNode(key+align,keylen-align) ;
      if (suffix && node->frequency() >= affix_ratio * suffix->frequency())
	 suffix->setFrequency(0) ;
      }
   return true ;
}

//----------------------------------------------------------------------

static bool find_max_frequency(const NybbleTrie* trie, uint32_t nodeindex, const uint8_t * /*key*/,
			       unsigned /*keylen*/, void *user_data)
{
   auto freqptr = reinterpret_cast<uint32_t*>(user_data) ;
   auto node = trie->node(nodeindex) ;
   auto freq = node->frequency() ;
   if (*freqptr == (uint32_t)~0) // parent node?
      *freqptr = 0 ;
   else if (freq > *freqptr)
      *freqptr = freq ;
   return true ;
}

//----------------------------------------------------------------------

static bool find_ngram_cutoff(const NybbleTrie* trie, uint32_t nodeindex,
			      const uint8_t * key, unsigned keylen,
			      void *user_data)
{
   auto enum_data = reinterpret_cast<NgramEnumerationData*>(user_data) ;
   if (keylen >= minimum_length ||
       (enum_data->m_max_length < minimum_length &&
	keylen >= enum_data->m_max_length))
      {
      auto n = trie->node(nodeindex) ;
      auto freq = n->frequency() ;
      if (freq > 0 && freq > enum_data->m_frequencies[0])
	 {
	 // an optimization to eliminate extraneous n-grams: if the
	 //   node has only a single child, and that child has the
	 //   same frequency, this means that every occurrence of the
	 //   current n-gram is a prefix of that child n-gram, so
	 //   don't bother counting it; generalized to single-child
	 //   nodes where the child has almost the same frequency.
	 uint32_t max_freq = (uint32_t)~0 ;
	 if (trie->enumerateChildren(nodeindex,(uint8_t*)key,
				     8*(keylen+1), 8*keylen, 
				     find_max_frequency, &max_freq) &&
	      (max_freq < affix_ratio * freq ||
	       (keylen == minimum_length && affix_ratio < MINLEN_AFFIX_RATIO &&
		max_freq < MINLEN_AFFIX_RATIO * freq)))
	    {
	    insert_frequency(freq,enum_data->m_frequencies,
			     enum_data->m_topK) ;
	    enum_data->m_count++ ;
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

static bool filter_ngrams(const NybbleTrie* trie, uint32_t nodeindex, const uint8_t *key,
			  unsigned keylen, void *user_data)
{
   auto enum_data = reinterpret_cast<NgramEnumerationData*>(user_data) ;
   if (keylen >= minimum_length ||
       (enum_data->m_max_length < minimum_length &&
	keylen >= enum_data->m_max_length))
      {
      auto node = trie->node(nodeindex) ;
      auto freq = node->frequency() ;
      uint32_t max_freq = (uint32_t)~0 ;
      if (freq >= enum_data->m_min_freq &&
	  trie->enumerateChildren(nodeindex,(uint8_t*)key,
				  8*(keylen+1), 8*keylen,
				  find_max_frequency, &max_freq) &&
	  (max_freq < affix_ratio * freq ||
	   (keylen == minimum_length && affix_ratio < MINLEN_AFFIX_RATIO &&
	    max_freq < MINLEN_AFFIX_RATIO * freq)))
	 {
//for(size_t i = 0;i<keylen;i++)print_quoted_char(stdout,key[i]);
//cout <<" "<<keylen<<"@"<<freq<<endl;
	 enum_data->m_ngrams->insert(key,keylen,freq,node->isStopgram()) ;
	 enum_data->m_inserted_ngram = true ;
	 if (keylen == enum_data->m_max_length)
	    enum_data->m_have_max_length = true ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

static NybbleTrie *restrict_ngrams(NybbleTrie *ngrams, unsigned top_K,
				   unsigned min_length, unsigned max_length,
				   unsigned minlen, bool &have_max_length,
				   bool show_threshold = false)
{
   if (minlen > max_length)
      minlen = max_length ;
   LocalAlloc<uint32_t> top_frequencies(top_K,true) ;
   NewPtr<NybbleTrie> new_ngrams(1) ;
   NgramEnumerationData enum_data(have_max_length) ;
   enum_data.m_ngrams = new_ngrams.begin() ;
   enum_data.m_min_length = min_length ;
   enum_data.m_max_length = max_length ;
   enum_data.m_frequencies = top_frequencies ;
   enum_data.m_topK = top_K ;
   LocalAlloc<uint8_t> keybuf(max_length+1) ;
   // figure out the threshold we need to limit the total n-grams to the
   //   desired number
   enum_data.m_count = 0 ;
   enum_data.m_min_freq = 0 ;
   unsigned required = top_K / (maximum_length - max_length + 3) ;
   if (!ngrams->enumerate(keybuf,max_length,find_ngram_cutoff,&enum_data)
       || enum_data.m_count < required)
      {
      cout << "Only " << enum_data.m_count << " distinct ngrams at length "
	   << max_length << ": collect more data" << endl ;
      if (max_length < maximum_length)
	 {
	 return nullptr ;
	 }
      }
//   uint32_t threshold = adjusted_threshold(top_frequencies) ;
   uint32_t threshold = top_frequencies[0] ;
   if (enum_data.m_count < top_K && max_length == maximum_length)
      threshold = 1 ;
   enum_data.m_min_freq = threshold ;
   if (show_threshold)
      {
      cout << "  Enumerating ngrams of length " << minlen << " to "
	   << max_length << " occurring at least " << threshold << " times"
	   << endl ;
      }
   enum_data.m_have_max_length = false ;
   if (!ngrams->enumerate(keybuf,max_length,filter_ngrams,&enum_data) ||
       !enum_data.m_inserted_ngram)
      {
      new_ngrams = nullptr ;
      }
   gc() ;
   return new_ngrams.move() ;
}

//----------------------------------------------------------------------

static NybbleTrie *count_ngrams(const char **filelist, unsigned num_files,
				NybbleTrie *ngrams,
				unsigned min_length, unsigned max_length,
				bool &have_max_length,
				bool skip_newlines, bool aligned)
{
   cout << "Counting n-grams up to length " << max_length << endl ;
   (void)read_files(filelist,num_files,false,&count_ngrams,ngrams,min_length,max_length,skip_newlines,aligned) ;
   unsigned minlen = minimum_length ;
   if (minlen > max_length)
      minlen = max_length ;
   unsigned top_K = set_oversampling(topK,min_length,minimum_length,aligned) ;
   LocalAlloc<uint32_t> top_frequencies(top_K,true) ;
   NewPtr<NybbleTrie> new_ngrams(1) ;
   NgramEnumerationData enum_data(have_max_length) ;
   enum_data.m_ngrams = new_ngrams.begin() ;
   enum_data.m_min_length = min_length ;
   enum_data.m_max_length = max_length ;
   enum_data.m_frequencies = top_frequencies ;
   enum_data.m_topK = top_K ;
   enum_data.m_alignment = alignment ;
   LocalAlloc<uint8_t> keybuf(max_length+1) ;
   // find the threshold at which to cut off ngrams to limit to topK
   if (verbose)
      {
      cout << "  Determining threshold for ngrams of length " << minlen
	   << " to " << max_length << endl ;
      }
   // remove any suffixes of an n-gram which have nearly the same
   //   frequency as the n-gram containing them (but don't remove
   //   minimum-length n-grams)
   for (size_t len = min_length + 2 ; len <= max_length ; len++)
      {
      enum_data.m_desired_length = len ;
      (void)ngrams->enumerate(keybuf,len,remove_suffix,&enum_data) ;
      }
   // figure out the threshold we need to limit the total n-grams to the
   //   desired number
   enum_data.m_count = 0 ;
   enum_data.m_min_freq = 1 ;
   unsigned required = top_K / (maximum_length - max_length + 3) ;
   enum_data.m_frequencies[0] = 0 ;
   if (!ngrams->enumerate(keybuf,max_length,find_ngram_cutoff,&enum_data)
       || enum_data.m_count < required)
      {
      cout << "Only " << enum_data.m_count << " distinct ngrams at length "
	   << max_length << ": collect more data" << endl ;
      if (max_length < maximum_length)
	 {
	 return nullptr ;
	 }
      }
   uint32_t threshold = adjusted_threshold(top_frequencies) ;
   if (enum_data.m_count < top_K /*&& max_length == maximum_length*/)
      threshold = 1 ;
   enum_data.m_min_freq = threshold ;
   cout << "  Enumerating ngrams of length " << minlen << " to "
	<< max_length << " occurring at least " << threshold << " times"
	<< endl ;
   enum_data.m_have_max_length = false ;
   if (!ngrams->enumerate(keybuf,max_length,filter_ngrams,&enum_data) ||
       !enum_data.m_inserted_ngram)
      {
      new_ngrams = nullptr ;
      }
   gc() ;
   return new_ngrams.move() ;
}

//----------------------------------------------------------------------

static void merge_bigrams(NybbleTrie *ngrams, const BigramCounts *bigrams,
			  bool scaled, uint64_t total_bytes)
{
   if (!ngrams || !bigrams)
      return ;
   uint8_t keybuf[3] ;
   uint32_t min_count = 2 ;
   for (unsigned c1 = 0 ; c1 < 256 ; c1++)
      {
      keybuf[0] = (uint8_t)c1 ;
      for (unsigned c2 = 0 ; c2 < 256 ; c2++)
	 {
	 uint32_t count = bigrams->count(c1,c2) ;
	 if (count < min_count)
	    continue ;
	 keybuf[1] = (uint8_t)c2 ;
	 if (scaled)
	    count = scaled_frequency(count,total_bytes,smoothing_power,
				     log_smoothing_power) ;
	 ngrams->insert(keybuf,2,count,false) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static bool add_ngram(const NybbleTrie* trie, uint32_t nodeindex, const uint8_t *key,
   unsigned keylen, void* user_data)
{
   auto unpacked = reinterpret_cast<LangIDMultiTrie*>(user_data) ;
   if (!unpacked)
      return true ;
   auto node = trie->node(nodeindex) ;
   if (node)
      {
      uint32_t freq = node->frequency() ;
      unpacked->insert(key,keylen,unpacked->currentLanguage(),freq,
		   node->isStopgram()) ;
      }
   return true ;
}

//----------------------------------------------------------------------

static void add_ngrams(const NybbleTrie *ngrams, uint64_t total_bytes,
		       const LanguageID &opts, const char *filename)
{
   if (ngrams)
      {
      auto num_langs = language_identifier->numLanguages() ;
      // add the new language ID to the global database
      auto langID = language_identifier->addLanguage(opts,total_bytes) ;
      if (langID < num_langs)
	 {
	 CharPtr spec = language_identifier->languageDescriptor(langID) ;
	 cerr << "Duplicate language specification " << spec
	      << " encountered in " << filename
	      << ",\n  ignoring data to avoid database errors." << endl ;
	 }
      auto trie = language_identifier->unpackedTrie() ;
      if (trie)
	 {
	 trie->setLanguage(langID) ;
	 uint8_t keybuf[10000] ;
	 ngrams->enumerate(keybuf,sizeof(keybuf),add_ngram,trie) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static const char *add_UTF8_range(const char *range_spec,
				  NybbleTrie *ngrams,
				  uint64_t &total_bytes)
{
   bool bad_spec = false ;
   if (range_spec)
      {
      while (*range_spec && isspace(*range_spec))
	 range_spec++ ;
      if (!*range_spec)
	 return range_spec ;
      char *end ;
      unsigned long first = strtoul(range_spec,&end,0) ;
      if (end != range_spec)
	 {
	 range_spec = end ;
	 unsigned long last = first ;
	 while (*range_spec && isspace(*range_spec))
	    range_spec++ ;
	 if (*range_spec == '-')
	    {
	    range_spec++ ;
	    while (*range_spec && isspace(*range_spec))
	       range_spec++ ;
	    last = strtoul(range_spec,&end,0) ;
	    if (end == range_spec)
	       {
	       bad_spec = true ;
	       last = first ;
	       }
	    else
	       range_spec = end ;
	    }
	 for (unsigned long i = first ; i <= last ; i++)
	    {
	    char utf8[6] ;
	    int bytes = UCS2_to_UTF8(i,utf8) ;
	    if (bytes > 0)
	       {
               ngrams->insert((uint8_t*)utf8,bytes,1,false) ;
	       total_bytes += (bytes * FAKED_NGRAM_DISCOUNT) ;
	       }
	    }
	 while (*range_spec && isspace(*range_spec))
	    range_spec++ ;
	 if (*range_spec == ',')
	    range_spec++ ;
	 else if (*range_spec && !isdigit(*range_spec))
	    bad_spec = true ;
	 }
      }
   if (bad_spec)
      {
      // oops, invalid spec, skip to end
      cerr << "Error in language range specification near\n\t"
	   << range_spec << endl << endl ;
      range_spec = strchr(range_spec,'\0') ;
      }
   return range_spec ;
}

//----------------------------------------------------------------------

static void add_UTF8_codepoints(NybbleTrie *ngrams,
				const char *cp_list, uint64_t &total_bytes)
{
   while (*cp_list)
      {
      cp_list = add_UTF8_range(cp_list,ngrams,total_bytes) ;
      }
   return ;
}

//----------------------------------------------------------------------

static void coverage_matches(const uint8_t *buf, unsigned buflen, unsigned *cover,
			     double *freqtotal, double *matchcount,
			     const NybbleTrie *ngrams, bool scaled)
{
   NybbleTriePointer ptr(ngrams) ;
   for (size_t i = 0 ; i < buflen ; i++)
      {
      if (!ptr.extendKey(buf[i]))
	 return ;
      const NybbleTrieNode *n = ptr.node() ;
      if (n && n->leaf())
	 {
	 (*matchcount) += 1.0 ;
//	 (*matchcount) += ::pow(i,0.75) ;
	 for (size_t j = 0 ; j <= i ; j++)
	    {
	    cover[j]++ ;
	    double freq ;
	    if (scaled)
	       freq = unscale_frequency(n->frequency(),smoothing_power) ;
	    else
	       freq = (n->frequency() / (double)TRIE_SCALE_FACTOR) ;
	    freqtotal[j] += freq ;
	    }
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static bool compute_coverage(PreprocessedInputFile *infile, va_list args)
{
   auto ngrams = va_arg(args,const NybbleTrie*) ;
   if (!ngrams)
      return false ;
   auto overall_cover = va_arg(args,size_t*) ;
   auto counted_cover = va_arg(args,size_t*) ;
   auto freq_cover = va_arg(args,double*) ;
   auto match_count = va_arg(args,double*) ;
   bool scaled = va_arg(args,int) ;
   unsigned maxlen = ngrams->longestKey() ;
   if (maxlen > ABSOLUTE_MAX_LENGTH)
      maxlen = ABSOLUTE_MAX_LENGTH ;
   LocalAlloc<uint8_t> buf(maxlen+1) ;
   LocalAlloc<unsigned> cover(maxlen+1) ;
   LocalAlloc<double> freqtotal(maxlen+1) ;
   unsigned buflen = 0 ;
   // initialize the statistics and prime the buffer
   for (size_t i = 0 ; i < maxlen ; i++)
      {
      if (!infile->moreData())
	 {
	 break ;
	 }
      buflen++ ;
      buf[i] = infile->getByte() ;
      cover[i] = 0 ;
      freqtotal[i] = 0.0 ;
      }
   bool hit_eod = false ;
   (*match_count) = 0 ;
   while (buflen > 0)
      {
      // check matches against current buffer
      coverage_matches(buf,buflen,cover,freqtotal,match_count,ngrams,scaled) ;
      // update statistics
      (*overall_cover) += (cover[0] != 0) ;
      (*counted_cover) += cover[0] ;
      (*freq_cover) += freqtotal[0] ;
      // update buffer
      memmove(buf,buf+1,buflen-1) ;
      if (!hit_eod && infile->moreData())
	 {
	 int b = infile->getByte() ;
	 if (b == EOF)
	    {
	    hit_eod = true ;
	    break ;
	    }
	 buf[buflen-1] = (uint8_t)b ;
	 }
      else
	 {
	 buflen-- ;
	 if (buflen == 0)
	    break ;
	 }
      // update statistics buffers
      memmove(cover,cover+1,buflen*sizeof(cover[0])) ;
      memmove(freqtotal,freqtotal+1,buflen*sizeof(freqtotal[0])) ;
      cover[buflen-1] = 0 ;
      freqtotal[buflen-1] = 0.0 ;
      }
   return true ;
}

//----------------------------------------------------------------------

static void compute_coverage(LanguageID &lang_info,
			     const char **filelist, unsigned num_files,
			     const NybbleTrie *ngrams, bool scaled)
{
   size_t overall_coverage = 0 ;	// percentage of training bytes covered by ANY ngram
   size_t counted_coverage = 0 ;	// coverage weighted by number of ngrams covering a byte
   double match_count = 0.0 ;		// number of matching ngrams against training data
   double freq_coverage = 0.0 ;		// coverage weighted by freq of ngrams covering a bytes
   uint64_t training_bytes
      = read_files(filelist,num_files,false,&compute_coverage,ngrams,&overall_coverage,
		   &counted_coverage,&freq_coverage,&match_count,scaled) ;
   if (training_bytes > 0)
      {
      if (verbose)
	 {
	 cout << "    Coverage fraction " << (overall_coverage / (double)training_bytes) << endl ;
	 }
      lang_info.setCoverageFactor(overall_coverage / (double)training_bytes) ;
      lang_info.setCountedCoverage(counted_coverage / (double)training_bytes) ;
      freq_coverage = ::sqrt(freq_coverage) ; // adjust for multiple-counting of high-freq ngrams
      lang_info.setFreqCoverage(freq_coverage) ;
      lang_info.setMatchFactor(match_count / training_bytes) ;
      }
   else
      {
      lang_info.setCoverageFactor(0.0) ;
      }
   return ;
}

//----------------------------------------------------------------------

static bool load_frequencies(CFile& f, NybbleTrie *ngrams,
			     uint64_t &total_bytes, bool textcat_format,
			     LanguageID &opts,
			     BigramCounts *&bigrams, bool &scaled)
{
   scaled = false ;
   if (!f)
      return false ;
   bool have_total_bytes = false ;
   bool first_line = true ;
   bool have_bigram_counts = false ;
   double codepoint_discount = 1.0 ;
   BigramCounts *crubadan_bigrams = nullptr ;
   if (crubadan_format)
      crubadan_bigrams = new BigramCounts ;
   bool have_script = false ;
   bool try_guessing_script = false ;
   opts.setCoverageFactor(1.0) ;
   while (CharPtr buffer = f.getCLine())
      {
      char *endptr ;
      if (textcat_format)
	 {
	 char *tab = strchr((char*)buffer,'\t') ;
	 if (tab)
	    {
	    *tab = '\0' ;
	    size_t len = tab - buffer ;
	    tab++ ;
	    if (len > 0)
	       {
	       // we have a valid frequency-list entry, so add it to the
	       //   trie and update the global byte count
	       unsigned long count = strtoul(tab,&endptr,0) ;
	       ngrams->increment((uint8_t*)*buffer,len,count) ;
	       total_bytes += (len * count * ASSUMED_NGRAM_DENSITY) ;
	       }
	    }
	 }
      else if (buffer[0] == '#' || buffer[1] == ';')
	 {
	 // comment line
	 continue ;
	 }
      else if (first_line && strncasecmp(buffer,"TotalCount:",11) == 0)
	 {
	 total_bytes = strtoul(buffer+11,nullptr,0) ;
	 if (total_bytes > 0)
	    have_total_bytes = true ;
	 }
      else if (strncasecmp(buffer,"Lang:",5) == 0)
	 {
	 char *arg = trim_whitespace(*buffer + 5) ;
	 opts.setLanguage(arg) ;
	 }
      else if (strncasecmp(buffer,"Region:",7) == 0)
	 {
	 char *arg = trim_whitespace(*buffer + 7) ;
	 opts.setRegion(arg) ;
	 }
      else if (strncasecmp(buffer,"Encoding:",9) == 0)
	 {
	 char *arg = trim_whitespace(*buffer + 9) ;
	 opts.setEncoding(arg) ;
	 try_guessing_script = !have_script ;
	 }
      else if (strncasecmp(buffer,"Source:",7) == 0)
	 {
	 char *arg = trim_whitespace(*buffer + 7) ;
	 opts.setSource(arg) ;
	 }
      else if (strncasecmp(buffer,"Script:",7) == 0)
	 {
	 char *arg = trim_whitespace(*buffer + 7) ;
	 opts.setScript(arg) ;
	 have_script = true ;
	 try_guessing_script = false ;
	 }
      else if (strncasecmp(buffer,"Scaled:",7) == 0)
	 {
	 scaled = true ;
	 }
      else if (strncasecmp(buffer,"IgnoreBlanks:",13) == 0)
	 {
	 ngrams->ignoreWhiteSpace() ;
	 }
      else if (strncasecmp(buffer,"Alignment:",10) == 0)
	 {
	 char *arg = trim_whitespace(*buffer + 10) ;
	 opts.setAlignment(arg) ;
	 }
      else if (strncasecmp(buffer,"BigramCounts:",13) == 0)
	 {
	 have_bigram_counts = true ;
	 break ;
	 }
      else if (strncasecmp(buffer,"Discount:",9) == 0)
	 {
	 // permit discounting of probabilities to ensure that a very small
	 //   model doesn't overwhelm one built with adequate training data
	 codepoint_discount = atof(buffer+9) ;
	 if (codepoint_discount < 1.0)
	    codepoint_discount = 1.0 ;
	 }
      else if (strncasecmp(buffer,"Coverage:",9) == 0)
	 {
	 opts.setCoverageFactor(atof(buffer+9)) ;
	 }
      else if (strncasecmp(buffer,"WeightedCoverage:",17) == 0)
	 {
	 opts.setCountedCoverage(atof(buffer+17)) ;
	 }
      else if (strncasecmp(buffer,"FreqCoverage:",13) == 0)
	 {
	 opts.setFreqCoverage(atof(buffer+13)) ;
	 }
      else if (strncasecmp(buffer,"MatchFactor:",12) == 0)
	 {
	 opts.setMatchFactor(atof(buffer+12)) ;
	 }
      else if (strncasecmp(buffer,"UTF8:",5) == 0)
	 {
	 // add the UTF-8 representations of the listed codepoints to the
	 //  model, highly discounted
	 add_UTF8_codepoints(ngrams,buffer+5,total_bytes) ;
	 }
      else
	 {
	 if (try_guessing_script)
	    {
	    have_script = opts.guessScript() ;
	    try_guessing_script = false ;
	    }
	 endptr = nullptr ;
	 char *bufptr = skip_whitespace(*buffer) ;
	 bool stopgram = false ;
	 if (*bufptr == '-')
	    {
	    stopgram = true ;
	    bufptr++ ;
	    }
	 unsigned long count = strtoul(bufptr,&endptr,0) ;
	 if (count == 0)
	    {
	    count = 1 ;
	    stopgram = true ;
	    }
	 if (endptr && endptr > buffer)
	    {
	    unsigned shift_count = 0 ;
	    if (*endptr >= 'a' && *endptr <= 'z' &&
		(endptr[1] == '\t' || endptr[1] == ' '))
	       {
	       if (!scaled)
		  {
		  cerr << "Warning: found scaled count in counts file that does not contain" << endl ;
		  cerr << "the Scaled: directive; enabling scaled counts for remainder." << endl ;
		  scaled = true ;
		  }
	       shift_count = *endptr++ - 'a' + 1 ;
	       }
	    while (*endptr == '\t' || *endptr == ' ')
	       endptr++ ;
	    char *key = endptr ;
	    char *keyptr = key ;
	    while (*endptr && *endptr != '\n' && *endptr != '\r')
	       {
	       if (*endptr == '\\')
		  {
		  // un-quote the next character
		  endptr++ ;
		  if (!*endptr)
		     break ;
		  switch (*endptr)
		     {
		     case '0':
			*keyptr++ = '\0' ;
			break ;
		     case 'n':
			*keyptr++ = '\n' ;
			break ;
		     case 'r':
			*keyptr++ = '\r' ;
			break ;
		     case 't':
			*keyptr++ = '\t' ;
			break ;
		     case 'f':
			*keyptr++ = '\f' ;
			break ;
		     default:
			*keyptr++ = *endptr ;
		     }
		  endptr++ ;
		  }
	       else
		  *keyptr++ = *endptr++ ;
	       }
	    size_t len = keyptr - key ;
	    if (len > 0)
	       {
	       if (crubadan_format)
		  {
		  if (*key == '<')
		     *key = ' ' ;
		  if (keyptr[-1] == '>')
		     keyptr[-1] = ' ' ;
		  // add in bigram counts from the current ngram
		  for (size_t i = 0 ; i+1 < len ; i++)
		     {
		     crubadan_bigrams->incr(key[i],key[i+1],count) ;
		     }
		  }
	       // we have a valid frequency-list entry, so add it to the
	       //   trie and update the global byte count
	       if (count > 1 || !crubadan_format)
		  {
		  if (scaled)
		     {
		     count <<= shift_count ;
		     ngrams->insert((uint8_t*)key,len,count,stopgram) ;
		     }
		  else
		     ngrams->increment((uint8_t*)key,len,count,stopgram) ;
		  if (crubadan_format && len > 3)
		     {
		     if (keyptr[-1] == ' ')
			len-- ;
		     ngrams->increment((uint8_t*)key,len,count) ;
		     if (len > 3 && *key == ' ')
			ngrams->increment((uint8_t*)(key+1),len-1,count) ;
		     }
		  }
	       if (!have_total_bytes)
		  total_bytes += (len * count / 4) ;
	       }
	    }
	 }
      first_line = false ;
      }
   if (have_bigram_counts)
      {
      delete crubadan_bigrams ;
      bigrams = new BigramCounts ;
      if (!bigrams->read(f))
	 {
	 cerr << "Error reading bigram counts in vocabulary file" << endl ;
	 delete bigrams ;
	 bigrams = nullptr ;
	 }
      }
   else
      {
      if (crubadan_bigrams)
	 crubadan_bigrams->scaleTotal(100) ;
      bigrams = crubadan_bigrams ;
      }
   total_bytes = (uint64_t)(total_bytes * codepoint_discount) ;
   return true ;
}

//----------------------------------------------------------------------

static bool load_frequencies(const char **filelist, unsigned num_files,
			     LanguageID &opts, bool textcat_format, bool no_save)
{
   cout << "Loading frequency list " ;
   if (textcat_format)
      cout << "(TextCat format)" << endl ;
   else if (crubadan_format)
      cout << "(Crubadan format)" << endl ;
   else
      cout << "(MkLangID format)" << endl ;
   NewPtr<NybbleTrie> ngrams(1) ;
   if (!ngrams)
      {
      SystemMessage::no_memory("while loading frequency lists") ;
      return false ;
      }
   uint64_t total_bytes = 0 ;
   BigramCounts *bigrams = nullptr ;
   bool scaled = false ;
   for (size_t i = 0 ; i < num_files ; i++)
      {
      const char *filename = filelist[i] ;
      CInputFile fp(filename) ;
      if (fp)
	 {
	 cout << "  Reading " << filename << endl ;
	 load_frequencies(fp,*ngrams,total_bytes,textcat_format,opts,bigrams,scaled) ;
	 if (textcat_format || num_files > 1)
	    {
	    delete bigrams ;
	    bigrams = nullptr ;
	    }
	 }
      }
   merge_bigrams(*ngrams,bigrams,scaled,total_bytes) ;
   delete bigrams ;
   if (ngrams->size() > 0)
      {
      // output the merged vocabulary list as text if requested
      minimum_length = 1 ;
      if (vocabulary_file)
	 dump_vocabulary(ngrams,scaled,vocabulary_file,1000,total_bytes,
			 opts) ;
      // now that we have read in the n-grams, augment the database with that
      //   list for the indicated language and encoding
      if (no_save)
	 {
	 if (!vocabulary_file)
	    cerr << "*** N-grams WERE NOT SAVED (read-only database) ***" << endl ;
	 }
      else
	 {
	 cout << "Updating database" << endl ;
	 if (!scaled)
	    ngrams->scaleFrequencies(total_bytes,smoothing_power,log_smoothing_power) ;
	 add_ngrams(ngrams,total_bytes,opts,filelist[0]) ;
	 }
      return true ;
      }
   else
      {
      return false ;
      }
}

//----------------------------------------------------------------------

// this global makes merge_ngram non-reentrant
static const PackedTrieFreq *base_frequency = nullptr ;
static unsigned *model_sizes = nullptr ;

static bool merge_ngrams(const PackedTrieNode *node, const uint8_t *key,
			 unsigned keylen, void *user_data)
{
   NybbleTrie **merged = (NybbleTrie**)user_data ;
   auto freqlist = node->frequencies(base_frequency) ;
   for ( ; freqlist ; freqlist = freqlist->next())
      {
      if (!freqlist->isStopgram())
	 {
	 unsigned id = freqlist->languageID() ;
	 if (model_sizes)
	    model_sizes[id]++ ;
	 if (merged[id])
	    {
	    uint32_t freq = freqlist->scaledScore() ;
	    (void)merged[id]->insertMax(key,keylen,freq,false) ;
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

static NybbleTrie *find_encoding(const char *enc_name, NewPtr<NybbleTrie*>& encodings,
				 NewPtr<LanguageID*>& enc_info,
				 unsigned &num_encs, unsigned &encs_alloc)
{
   if (!enc_name || !*enc_name)
      return nullptr ;
   for (unsigned id = 0 ; id < num_encs ; id++)
      {
      if (enc_info[id] && enc_info[id]->encoding() &&
	  strcmp(enc_info[id]->encoding(),enc_name) == 0)
	 {
	 return encodings[id] ;
	 }
      }
   // if we get to this point, the specified encoding has not yet been seen
   //   so allocate a new trie and languageID record for the encoding
   if (num_encs >= encs_alloc)
      {
      unsigned new_alloc = 2 * encs_alloc ;
      if (encodings.reallocate(encs_alloc,new_alloc) &&
	  enc_info.reallocate(encs_alloc,new_alloc))
	 {
	 encs_alloc = new_alloc ;
	 }
      }
   if (num_encs < encs_alloc)
      {
      encodings[num_encs] = new NybbleTrie ;
      enc_info[num_encs] = new LanguageID("CLUS=Clustered","XX",enc_name,"merged") ;
      return encodings[num_encs++] ;
      }
   else
      return nullptr ;
}

//----------------------------------------------------------------------

static bool cluster_models_by_charset(LanguageIdentifier *clusterdb,
				      const char *cluster_dbfile)
{
   unsigned num_encs = 0 ;
   unsigned encs_alloc = 50 ;
   NewPtr<NybbleTrie*> encodings(encs_alloc) ;
   NewPtr< LanguageID*> enc_info(encs_alloc) ;
   // make a mapping from language ID to per-encoding merged models
   unsigned numlangs = language_identifier->numLanguages() ;
   NybbleTrie **merged = new NybbleTrie*[numlangs] ;
   for (unsigned langid = 0 ; langid < numlangs ; langid++)
      {
      // get the character encoding for the current model and find the
      //   merged trie for that encoding
      const char *enc_name = language_identifier->languageEncoding(langid) ;
      merged[langid] = find_encoding(enc_name,encodings,enc_info,num_encs,encs_alloc) ;
      if (!merged[langid])
	 {
	 SystemMessage::no_memory("while merging language models") ;
	 break ;
	 }
      }
   // iterate over all of the ngrams in the database, merging each frequency
   //   record into the appropriate per-encoding model
   auto ptrie = language_identifier->packedTrie() ;
   base_frequency = ptrie->frequencyBaseAddress() ;
   model_sizes = new unsigned[numlangs] ;
   std::fill_n(model_sizes,numlangs,0) ;
   unsigned maxkey = ptrie->longestKey() ;
   LocalAlloc<uint8_t> keybuf(maxkey) ;
   ptrie->enumerate(keybuf,maxkey,merge_ngrams,merged) ;
   base_frequency = nullptr ;
   // figure out the maximum size of an individual model for each encoding
   LocalAlloc<unsigned> max_sizes(num_encs) ;
   std::fill_n(&max_sizes,num_encs,0) ;
   for (size_t i = 0 ; i < num_encs ; i++)
      {
      for (unsigned langid = 0 ; langid < numlangs ; langid++)
	 {
	 if (merged[langid] == encodings[i] && model_sizes[langid] > max_sizes[i])
	    max_sizes[i] = model_sizes[langid] ;
	 }
      }
   // collect all of the merged models into the new language database
   LanguageIdentifier *ident = language_identifier ;
   (void)clusterdb->unpackedTrie() ; // ensure that we are unpacked
   language_identifier = clusterdb ;
   for (unsigned i = 0 ; i < num_encs ; i++)
      {
      bool have_max_length = true ;
      cerr << "adding encoding " << i << "  " << enc_info[i]->encoding() << endl;
      Owned<NybbleTrie> clustered = restrict_ngrams(encodings[i],2*max_sizes[i],1,maxkey,1,have_max_length) ;
      delete encodings[i] ;
      uint64_t total_bytes = 1 ; //FIXME!!!
      add_ngrams(clustered,total_bytes,*(enc_info[i]),"???") ;
      delete enc_info[i] ;
      }
   save_database(cluster_dbfile) ;
   language_identifier = ident ;
   delete[] model_sizes ;
   model_sizes = nullptr ;
   return true ;
}

//----------------------------------------------------------------------

static bool cluster_models(const char *cluster_db_name, double cluster_thresh)
{
   if (cluster_thresh < 0.0 || cluster_thresh > 1.0)
      return false ;
   auto clusterdb = new LanguageIdentifier(cluster_db_name) ;
   if (!clusterdb)
      {
      cerr << "Unable to create clustered language database " 
	   << cluster_db_name << endl ;
      return false ;
      }
   if (clusterdb->numLanguages() > 0)
      {
      cerr << "Non-empty language database specified: "
	   << cluster_db_name << endl ;
      return false ;
      }
   if (cluster_thresh == 0.0)
      {
      // cluster all models with the same character set together
      return cluster_models_by_charset(clusterdb,cluster_db_name) ;
      }
   else
      {
      //TODO: nonzero clustering thresholds
      cerr << "clustering thresholds other than 0.0 not implemented yet."
	   << endl ;
      return false ;
      }
}

//----------------------------------------------------------------------

static bool compute_ngrams(const char **filelist, unsigned num_files,
			   NybbleTrie *&ngrams,
			   bool skip_newlines, bool omit_bigrams,
			   bool ignore_whitespace, uint64_t &total_bytes,
			   bool aligned)
{
   // accumulate the most common n-grams (for n >= 3) by
   //   1. counting all trigrams
   //   2. filtering down to the K most frequent
   //   3. counting all n-grams (n > 3) which have the trigrams from step
   //         2 as a prefix
   //   4. filtering down the result of step 3 to the K most frequent overall
   // since we get bigrams counts essentially for free after step 1, also
   // add in the complete set of nonzero bigram counts while we're at it
   auto counts = new TrigramCounts ;
   ngrams = nullptr ;
   if (!counts)
      return false ;
   ngrams = new NybbleTrie ;
   if (!ngrams)
      {
      delete counts ;
      return false ;
      }
   BigramCounts *bi_counts = nullptr ;
   BigramCounts **bigram_ptr = omit_bigrams ? nullptr : &bi_counts ;
   total_bytes = count_trigrams(filelist,num_files,*counts,
				skip_newlines,aligned,bigram_ptr) ;
   unsigned top_K = set_oversampling(topK,ABSOLUTE_MIN_LENGTH,minimum_length,
				     aligned) ;
   counts->filter(top_K,maximum_length,verbose) ;
   ngrams->ignoreWhiteSpace(ignore_whitespace) ;
   if (counts->enumerate(*ngrams) && ngrams->longestKey() > 0)
      {
      delete counts ;
      counts = nullptr ;
      bool small_data = (total_bytes * top_K) < 1E11 ;
      // start with five-grams, and continue increasing length until we
      //   no longer have any maximal-length ngrams included in the top K
      bool have_max_length = false ;
      // if we're processing virtual 16-bit units, extend the length
      //   by multiples of 16 bits
      unsigned expansion = (aligned || bigram_extension != BigramExt_None) ? 2 : 1 ;
      unsigned min_length = 4 ; // we've already counted up to 3-grams
      unsigned max_length = 4 + (expansion * (small_data ? 2 : 1)) ;
      if (max_length > maximum_length)
	  max_length = maximum_length ;
      do {
         NybbleTrie *new_ngrams
	    = count_ngrams(filelist, num_files, ngrams, min_length,
			   max_length, have_max_length, skip_newlines,
			   aligned) ;
	 if (new_ngrams)
	    {
	    delete ngrams ;
	    ngrams = new_ngrams ;
	    }
	 else
	    have_max_length = false ;
	 // update lengths for next iteration
	 min_length = max_length + 1 ;
	 unsigned increment
	    = expansion + expansion*(max_length-ABSOLUTE_MIN_LENGTH + 1) / 2 ;
	 if (increment > expansion * MAX_INCREMENT)
	    increment = expansion * MAX_INCREMENT ;
	 if (small_data)
	    increment *= 2 ;
	 max_length += increment ;
	 if (max_length > maximum_length)
	    max_length = maximum_length ;
         } while (have_max_length && min_length <= maximum_length) ;
      }
   if (!omit_bigrams)
      {
      merge_bigrams(ngrams,bi_counts,false,total_bytes) ;
      }
   delete bi_counts ;
   delete counts ;
   return true ;
}

//----------------------------------------------------------------------

static bool process_files(const char **filelist, unsigned num_files,
			  const LanguageID &base_opts, NybbleTrie *curr_ngrams,
			  uint64_t training_bytes,
			  bool skip_newlines, bool omit_bigrams,
			  bool ignore_whitespace, NybbleTrie *stop_grams,
			  const NybbleTrie *ngram_weights, bool no_save,
			  bool /*check_script TODO*/)
{
   LanguageID opts(base_opts) ;
   NybbleTrie *ngrams ;
   uint64_t total_bytes = 0 ;
   bool scaled = false ;
   if (curr_ngrams && curr_ngrams->size() > 0)
      {
      cout << "Using baseline n-gram model from language database" << endl ;
      ngrams = curr_ngrams ;
      scaled = true ;
      total_bytes = training_bytes ;
      }
   else if (!compute_ngrams(filelist,num_files,ngrams,skip_newlines,
			    omit_bigrams,ignore_whitespace,total_bytes,
			    opts.alignment() > 1))
      return false ;
   compute_coverage(opts,filelist,num_files,ngrams,scaled) ;
   add_stop_grams(filelist,num_files,ngrams,stop_grams,ngram_weights,scaled) ;
   // output the vocabulary list as text if requested
   if (vocabulary_file)
      {
      unsigned max_length = ngrams->longestKey() ;
      dump_vocabulary(ngrams,scaled,vocabulary_file,max_length,
		      total_bytes,opts);
      }
   // now that we have the top K n-grams, augment the database with that
   //   list for the indicated language and encoding
   if (!no_save)
      {
      if (!scaled)
	 ngrams->scaleFrequencies(total_bytes,smoothing_power,log_smoothing_power) ;
      add_ngrams(ngrams,total_bytes,opts,filelist[0]) ;
      }
   else if (!vocabulary_file)
      {
      cerr << "*** N-grams WERE NOT SAVED (read-only database) ***" << endl ;
      }
   delete ngrams ;
   return true ;
}

/************************************************************************/
/*	Main Program							*/
/************************************************************************/

static void parse_bigram_extension(const char *arg)
{
   if (arg && *arg)
      {
      if (*arg == 'b' || *arg == 'B')
	 bigram_extension = BigramExt_ASCIIBigEndian ;
      else if (*arg == 'l' || *arg == 'L')
	 bigram_extension = BigramExt_ASCIILittleEndian ;
      else if (*arg == '-' || *arg == 'n' || *arg == 'N')
	 bigram_extension = BigramExt_None ;
      else
	 {
	 cerr << "Invalid value for -2 flag; bigram extension disabled."
	      << endl ;
	 bigram_extension = BigramExt_None ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static void parse_UTF8_extension(const char *arg)
{
   if (arg && *arg)
      {
      if (*arg == 'b' || *arg == 'B')
	 bigram_extension = BigramExt_UTF8BigEndian ;
      else if (*arg == 'l' || *arg == 'L')
	 bigram_extension = BigramExt_UTF8LittleEndian ;
      else if (*arg == '-' || *arg == 'n' || *arg == 'N')
	 bigram_extension = BigramExt_None ;
      else
	 {
	 cerr << "Invalid value for -8 flag; UTF8-to-UTF16 conversion disabled."
	      << endl ;
	 bigram_extension = BigramExt_None ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static void parse_clustering(const char *arg, double &threshold,
			     const char *&dbfile)
{
   if (arg && *arg)
      {
      char *endptr = nullptr ;
      double val = strtod(arg,&endptr) ;
      if (val < 0.0)
	 {
	 val = 0.0 ;
	 cerr << "-C threshold adjusted to 0.0" << endl ;
	 }
      else if (val > 1.0)
	 {
	 val = 1.0 ;
	 cerr << "-C threshold adjusted to 1.0" << endl ;
	 }
      if (endptr && *endptr == ',')
	 {
	 dbfile = endptr + 1 ;
	 threshold = val ;
	 }
      else
	 {
	 cerr << "-C flag missing filename" << endl ;
	 dbfile = nullptr ;
	 threshold = -1.0 ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static void parse_byte_limit(const char *spec)
{
   if (spec)
      {
      if (spec[0] == '@')
	 {
	 subsample_input = true ;
	 spec++ ;
	 }
      byte_limit = atol(spec) ;
      }
   return ;
}

//----------------------------------------------------------------------

static void parse_smoothing_power(const char *spec)
{
   if (spec && *spec)
      {
      double smooth = atof(spec) ;
      // negative smoothing is now logarithmic rather than exponentiation
      if (smooth < 0.0)
	 {
	 // since we're using logs, it makes sense for the user-provided
	 //   power to be a log as well.
	 smoothing_power = -::pow(10.0,-smooth) ;
	 log_smoothing_power = scaling_log_power(smoothing_power) ;
	 }
      else
	 {
	 if (smooth > 5.0)
	    smooth = 5.0 ;
	 smoothing_power = smooth ;
	 log_smoothing_power = 1.0 ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static void parse_translit(const char *spec, char *&from, char *&to)
{
   from = nullptr ;
   to = nullptr ;
   if (spec)
      {
      from = dup_string(spec).move() ;
      char *comma = strchr(from,',') ;
      if (comma)
	 {
	 if (comma[1])
	    to = dup_string(comma+1).move() ;
	 if (comma != from)
	    *comma = '\0' ;
	 else
	    {
	    cerr << "You may not omit the FROM encoding for -T" << endl ;
	    delete[] from ;
	    from = nullptr ;
	    }
	 }
      }
   return;
}

//----------------------------------------------------------------------

static bool process_argument_group(int &argc, const char **&argv,
				   LanguageID &lang_info, bool no_save,
				   const char *argv0)
{
   // reset any options which must be specified separately for each file group
   vocabulary_file = nullptr ;
   bool frequency_list = false ;
   bool frequency_textcat = false ;
   bool skip_newlines = false ;
   bool omit_bigrams = false ;
   bool end_of_args = false ;
   bool ignore_whitespace = false ;
   const char *related_langs = nullptr ;
   const char *cluster_db = nullptr ;
   double cluster_thresh = -1.0 ;  // never cluster
   char *from = nullptr ;
   char *to = nullptr ;
   crubadan_format = false ;
   PreprocessedInputFile::setDefaultConvertLatin1(false) ;
   byte_limit = ~0 ;
   // process any switches
   while (argc > 1 && argv[1][0] == '-')
      {
      if (argv[1][1] == '-')
	 {
	 // no more flag arguments
	 end_of_args = true ;
	 argv++ ;
	 argc-- ;
	 break ;
	 }
      switch (argv[1][1])
	 {
	 case '1': PreprocessedInputFile::setDefaultConvertLatin1(true) ; break ;
	 case '2': parse_bigram_extension(get_arg(argc,argv)) ; break ;
	 case '8': parse_UTF8_extension(get_arg(argc,argv)) ; 	break ;
	 case 'C': parse_clustering(get_arg(argc,argv),
				    cluster_thresh,cluster_db) ; break ;
	 case 'D': do_dump_trie = true ;			break ;
	 case 'l': lang_info.setLanguage(get_arg(argc,argv)) ;	break ;
	 case 'r': lang_info.setRegion(get_arg(argc,argv)) ;	break ;
	 case 'e': lang_info.setEncoding(get_arg(argc,argv)) ;	break ;
	 case 's': lang_info.setSource(get_arg(argc,argv)) ;	break ;
	 case 'W': lang_info.setScript(get_arg(argc,argv)) ;	break ;
	 case 'k': topK = atoi(get_arg(argc,argv)) ;		break ;
	 case 'm': minimum_length = atoi(get_arg(argc,argv)) ;	break ;
	 case 'M': maximum_length = atoi(get_arg(argc,argv)) ;	break ;
	 case 'i': ignore_whitespace = true ;			break ;
	 case 'n': skip_newlines = true ;
	           if (argv[1][2] == 'n') skip_numbers = true ;
		   break ;
	 case 'a': affix_ratio = atof(get_arg(argc,argv)) ;	break ;
	 case 'A': alignment = atoi(get_arg(argc,argv)) ;	break ;
	 case 'b': omit_bigrams = true ;			break ;
	 case 'B': unique_boost = atof(get_arg(argc,argv)) ;	break ;
	 case 'd': discount_factor = atof(get_arg(argc,argv)) ; break ;
	 case 'O': max_oversample = atof(get_arg(argc,argv)) ;  break ;
	 case 'f': frequency_list = true ;
		   frequency_textcat = argv[1][2] == 't' ;
		   crubadan_format = argv[1][2] == 'c' ;	break ;
         case 'L': parse_byte_limit(get_arg(argc,argv)) ;	break ;
	 case 'R': related_langs = get_arg(argc,argv) ;		break ;
	 case 'S': parse_smoothing_power(get_arg(argc,argv)) ;	break ;
	 case 'T': parse_translit(get_arg(argc,argv),from,to) ;	break ;
	 case 'v': verbose = true ;				break ;
	 case 'x': store_similarities = true ;			break ;
	 case 'w': vocabulary_file = argv[1]+2 ;		break ;
	 case 'h':
	 default: usage(argv0,argv[1]) ;			break ;
	 }
      argc-- ;
      argv++ ;
      }
   if (unique_boost < 1.0)
      unique_boost = 1.0 ;
   // enforce a valid alignment size
   if (alignment > 4)
      alignment = 4 ;
   else if (alignment == 3)
      alignment = 2 ;
   else if (alignment < 1)
      alignment = 1 ;
   lang_info.setAlignment(alignment) ;
   // setup defaults for reading in files
   PreprocessedInputFile::setSampling(byte_limit,subsample_input) ;
   PreprocessedInputFile::setDefaultBigramExt(bigram_extension) ;
   PreprocessedInputFile::setDefaultAlignment(alignment) ;
   PreprocessedInputFile::setIgnoreWhitespace(ignore_whitespace) ;
   if (byte_limit < (uint64_t)~0U && verbose)
      cout << "Limiting training to " << byte_limit << " bytes" << endl ;
   if (minimum_length < ABSOLUTE_MIN_LENGTH && !frequency_list)
      {
      minimum_length = ABSOLUTE_MIN_LENGTH ;
      cerr << "Minimum length adjusted to " << ABSOLUTE_MIN_LENGTH << endl ;
      }
   if (bigram_extension != BigramExt_None && minimum_length < 4)
      minimum_length = 4 ;
   if (maximum_length > ABSOLUTE_MAX_LENGTH)
      {
      maximum_length = ABSOLUTE_MAX_LENGTH ;
      cerr << "Maximum length adjusted to " << ABSOLUTE_MAX_LENGTH << endl ;
      }
   if (crubadan_format)
      {
      // set defaults for files from the Crubadan crawler
      lang_info.setEncoding("utf8") ;
      if (!lang_info.source() || !*lang_info.source())
	 lang_info.setSource("Crubadan-Project") ;
      }
   bool check_script = false ;
   if (!lang_info.guessScript())
      {
      if (strcasecmp(lang_info.encoding(),"utf8") == 0 ||
	  strcasecmp(lang_info.encoding(),"utf-8") == 0 ||
	  strncasecmp(lang_info.encoding(),"utf16",5) == 0 ||
	  strncasecmp(lang_info.encoding(),"utf-16",6) == 0)
	 {
	 check_script = true ;
	 }
      }
   if (maximum_length < minimum_length)
      maximum_length = minimum_length ;
   // check the affix ratio and adjust as needed
   if (affix_ratio > 1.0)
      affix_ratio = 2.0 ;  // disable affix filtering
   else if (affix_ratio < MIN_AFFIX_RATIO)
      affix_ratio = MIN_AFFIX_RATIO ;
   // now accumulate all the named files until we hit the end of the command
   //   line or another switch
   const char **filelist = argv + 1 ;
   while (argc > 1 && argv[1] && (end_of_args || argv[1][0] != '-'))
      {
      argc-- ;
      argv++ ;
      }
   bool success = false ;
   if (cluster_db && *cluster_db)
      {
      success = cluster_models(cluster_db,cluster_thresh) ;
      }
   else if (frequency_list)
      {
      while (filelist <= argv)
	 {
	 unsigned filecount = 1 ;
	 if (frequency_textcat)
	    filecount = (argv - filelist + 1) ;
	 LanguageID local_lang_info(&lang_info) ;
	 if (load_frequencies(filelist,filecount,local_lang_info,
			      frequency_textcat,no_save))
	    success = true ;
	 filelist += filecount ;
	 }
      }
   else
      {
      // check for a transliteration request
      CharPtr translit_to
	 = from ? aprintf("%s//TRANSLIT",to ? to : lang_info.encoding()) : nullptr ;
      if (!PreprocessedInputFile::setDefaultTransliteration(from,translit_to))
	 {
	 cerr << "Unable to perform conversion from " << from << " to " << translit_to
	      << endl ;
	 }
      NybbleTrie *curr_ngrams = nullptr ;
      NybbleTrie *ngram_weights = nullptr ;
      uint64_t training_bytes = 0 ;
      NybbleTrie *stop_grams = load_stop_grams(&lang_info,related_langs,
					       curr_ngrams,ngram_weights,
					       training_bytes) ;
      success = process_files(filelist,argv-filelist+1,lang_info,
			      curr_ngrams,training_bytes,
			      skip_newlines,omit_bigrams,ignore_whitespace,
			      stop_grams,ngram_weights,no_save,check_script) ;
      delete stop_grams ;
      delete ngram_weights ;
      }
   delete[] from ;
   delete[] to ;
   return success ;
}

//----------------------------------------------------------------------

static int real_main(int argc, const char **argv)
{
   const char *argv0 = argv[0] ;
   const char *database_file = DEFAULT_LANGID_DATABASE ;
   bool no_save = false ;
   if (argc > 1 && argv[1][0] == '=')
      {
      if (argv[1][1] == '=')
	 {
	 no_save = true ;
	 database_file = argv[1]+2 ;
	 }
      else
	 database_file = argv[1]+1 ;
      argv++ ;
      argc-- ;
      }
   if (argc < 2)
      {
      usage(argv0,nullptr) ;
      return 1 ;
      }
   language_identifier = load_language_database(database_file,"",true) ;
   bool success = false ;
   LanguageID lang_info("en","US","utf-8",nullptr) ;
   while (argc > 1)
      {
      if (process_argument_group(argc,argv,lang_info,no_save,argv0))
	 {
	 success = true ;
	 }
      }
   if (success && !no_save)
      save_database(database_file) ;
   delete language_identifier ;
   language_identifier = nullptr ;
   return 0 ;
}

//----------------------------------------------------------------------

int main(int argc, const char **argv)
{
   Initialize() ;
   int status = real_main(argc,argv) ;
   Shutdown() ;
   return status ;
}

// end of file mklangid.C //
