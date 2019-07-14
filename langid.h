/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     langid.h							*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-14						*/
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

#ifndef __LANGID_H_INCLUDED
#define __LANGID_H_INCLUDED

#include "ptrie.h"

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

// current binary file format version
#define LANGID_FILE_VERSION 6
#define LANGID_FILE_SIGNATURE "Language Identification Database\r\n\x1A\004\0"

// minimum file version still supported
#define LANGID_MIN_FILE_VERSION 6

// reserved space for future additions to the file header
#define LANGID_PADBYTES_1  63

#define LANGID_FILE_DMOFFSET  96

// version 1-4 file format uses fixed-length string fields for simplicity
#define LANGID_STRING_LENGTH 64

#ifndef DBDIR
#define DBDIR "/usr/share/langident"
#endif

#ifndef DEFAULT_LANGID_DATABASE
#define DEFAULT_LANGID_DATABASE DBDIR "/languages.db"
#endif

#ifndef ALTERNATE_LANGID_DATABASE
#define ALTERNATE_LANGID_DATABASE "~/.langident/languages.db"
#endif

#ifndef FALLBACK_LANGID_DATABASE
#define FALLBACK_LANGID_DATABASE "./languages.db"
#endif

#ifndef DEFAULT_CHARSET_DATABASE
#define DEFAULT_CHARSET_DATABASE DBDIR "/charsets.db"
#endif

#ifndef ALTERNATE_CHARSET_DATABASE
#define ALTERNATE_CHARSET_DATABASE "~/.langident/charsets.db"
#endif

#ifndef FALLBACK_CHARSET_DATABASE
#define FALLBACK_CHARSET_DATABASE "./charsets.db"
#endif

// since the bigram byte model is much weaker than the long-ngram
//   model, give it proportionally less weight so that it basically
//   acts as a tie-breaker when there are no long ngram hits
#ifndef DEFAULT_BIGRAM_WEIGHT
#define DEFAULT_BIGRAM_WEIGHT 0.15
#endif

// consider any language score up to this value to be the same as zero
//   to avoid random noise
#define LANGID_ZERO_SCORE 0.01

// how much above the minimal score must a language score be to be considered
//   even a guess
#define GUESS_CUTOFF (20 * LANGID_ZERO_SCORE)

// how much above the minimal score must a language score be to be considered
//   a reliable identification (and not get a question mark)?
#define UNSURE_CUTOFF (120 * LANGID_ZERO_SCORE)

// at what point are we so sure that we don't flag the identification even
//   if it is highly ambiguous?
#define SURE_THRESHOLD (800 * LANGID_ZERO_SCORE)

// the range of length- and frequency-weighted ngram coverages; this
//   determines the scaling of the 32-bit integer actually stored in
//   the binary model file
#define MAX_WEIGHTED_COVER  32.0
#define MAX_FREQ_COVER      100.0
#define MAX_MATCH_FACTOR    16.0

/************************************************************************/
/************************************************************************/

class NybbleTrie ;

class TrigramCounts
   {
   private:
      uint32_t m_counts[256 * 256 * 256] ;
   public:
      TrigramCounts() { std::fill_n(m_counts,lengthof(m_counts),0) ; }
      TrigramCounts(const TrigramCounts *orig) ;
      ~TrigramCounts() {}

      // accessors
      uint32_t count(uint8_t c1, uint8_t c2, uint8_t c3) const
	 { return m_counts[(c1 << 16) + (c2 << 8) + c3] ; }
      uint32_t totalCount(uint8_t c1, uint8_t c2) const ;
      bool enumerate(NybbleTrie &ngrams) const ;

      // modifiers
      void copy(const TrigramCounts *orig) ;
      void clear(uint8_t c1, uint8_t c2, uint8_t c3)
	 { m_counts[(c1 << 16) + (c2 << 8) + c3] = 0 ; }
      void incr(uint8_t c1, uint8_t c2, uint8_t c3, uint32_t cnt = 1)
	 { m_counts[(c1 << 16) + (c2 << 8) + c3] += cnt ; }
      void filter(unsigned K, unsigned max_len, bool verbose) ;
      void filter(int32_t threshold) ;

      // I/O
      static TrigramCounts *load(Fr::CFile& f) ;
      bool read(Fr::CFile& f) ;
      bool save(Fr::CFile& f) const ;
   } ;

//----------------------------------------------------------------------

class BigramCounts
   {
   private:
      uint64_t m_total ;
      uint32_t m_counts[256 * 256] ;
   public:
      BigramCounts() { std::fill_n(m_counts,lengthof(m_counts),0) ; m_total = 0 ; }
      BigramCounts(Fr::CFile& f) ;
      BigramCounts(const BigramCounts *) ;
      BigramCounts(const TrigramCounts &) ;
      BigramCounts(const TrigramCounts *) ;
      ~BigramCounts() {}

      // accessors
      uint32_t count(uint8_t c1, uint8_t c2) const { return m_counts[(c1 << 8) + c2] ; }
      uint64_t totalCount() const { return m_total ; }
      double probability(uint8_t c1, uint8_t c2) const
	 { return this->count(c1,c2) / (double)this->totalCount() ; }
      double averageProbability(const char *buffer, size_t buflen) const ;

      // modifiers
      void copy(const BigramCounts *orig) ;
      void clear(uint8_t c1, uint8_t c2)
	 { m_counts[(c1 << 8) + c2] = 0 ; }
      void set(uint8_t c1, uint8_t c2, uint32_t cnt)
	 { m_counts[(c1 << 8) + c2] = cnt ; }
      void incr(uint8_t c1, uint8_t c2, uint32_t cnt = 1)
	 { m_counts[(c1 << 8) + c2] += cnt ; }
      void scaleTotal(unsigned factor) { m_total *= factor ; }

      // I/O
      static BigramCounts *load(Fr::CFile& f) ;
      bool read(Fr::CFile& f) ;
      bool readBinary(Fr::CFile& f) ;
      bool dumpCounts(Fr::CFile& f) const ;
      bool save(Fr::CFile& f) const ;
   } ;

//----------------------------------------------------------------------

class LanguageID
   {
   public:
      LanguageID() ;
      LanguageID(const char *lang, const char *reg, const char *enc,
		 const char *source = nullptr, const char *script = "UNKNOWN") ;
      LanguageID(const LanguageID &orig) ;
      LanguageID(const LanguageID *orig) ;
      LanguageID& operator= (LanguageID&& orig) ;
      ~LanguageID() ;

      // accessors
      const char *language() const { return m_language ; }
      const char *friendlyName() const { return m_friendlyname ; }
      const char *region() const { return m_region ; }
      const char *encoding() const { return m_encoding ; }
      const char *source() const { return m_source ; }
      const char *script() const { return m_script ; }
      unsigned alignment() const { return m_alignment ; }
      double coverageFactor() const { return m_coverage > 0.0 ? m_coverage : 1.0 ; }
      double countedCoverage() const { return m_countcover ; }
      double freqCoverage() const { return m_freqcover ; }
      double matchFactor() const { return m_matchfactor ; }
      uint64_t trainingBytes() const { return m_trainbytes ; }

      // modifiers
      void setLanguage(const char *lang, const char *friendly = nullptr) ;
      void setRegion(const char *region) ;
      void setEncoding(const char *encoding) ;
      void setSource(const char *source) ;
      void setScript(const char *scr) ;
      void setAlignment(unsigned align) { m_alignment = align ; }
      void setAlignment(const char *align) ;
      void setCoverageFactor(double coverage) ;
      void setCountedCoverage(double coverage) ;
      void setFreqCoverage(double coverage) ;
      void setMatchFactor(double match) ;
      void setTraining(uint64_t train_bytes) { m_trainbytes = train_bytes ; }
      bool guessScript() ;
      // operators
      bool sameLanguage(const LanguageID &other, bool ignore_region) const ;
      bool matches(const LanguageID *lang_info) const ;
      bool matches(const char *language, const char *region,
		   const char *encoding, const char *source) const ;
      bool operator == (const LanguageID &) const ;

      // I/O
      static LanguageID* read(Fr::CFile& f, unsigned version) ;
      static bool read(Fr::CFile& f, LanguageID *langID, unsigned version) ;
      bool write(Fr::CFile& f) const ;

   protected:
      void clear() ;

   private:
      Fr::CharPtr m_language ;
      Fr::CharPtr m_region ;
      Fr::CharPtr m_encoding ;
      Fr::CharPtr m_source ;
      Fr::CharPtr m_script ;
      const char *m_friendlyname ;
      double m_coverage ;		// percent of training covered by ngrams
      double m_countcover ;		// coverage weighted by count of matches
      double m_freqcover ;		// coverage weighted by frequencies of matches
      double m_matchfactor ;
      uint64_t m_trainbytes ;
      unsigned m_alignment ;
   } ;

//----------------------------------------------------------------------

class LanguageScores
   {
   public: // types
      class Info
	 {
	 public:
//	    Info() {}
	    void init(double sc, unsigned short new_id)
	       { m_score = sc ; m_id = new_id ; }

	    // accessors
	    double score() const { return m_score ; }
	    unsigned short id() const { return m_id ; }

	    // manipulators
	    void setScore(double sc) { m_score = sc ; }
	    void incrScore(double inc) { m_score += inc ; }
	    void decrScore(double dec) { m_score -= dec ; }
	    void setLang(unsigned short id) { m_id = id ; }
	    Info& operator= (double sc) { m_score = sc ; return *this ; }
	    Info& operator= (int id) { m_id = id ; m_score = 0.0 ; return *this ; }

	    static void swap(Info&, Info&) ;

	    // comparison
	    static int compare(const Info&, const Info&) ;
	    bool operator< (const Info& other) const { return compare(*this,other) < 0 ; }
	    bool operator== (double sc) { return sc == m_score ; }

	 private:
	    double   	   m_score ;
	    unsigned short m_id ;
	 } ;
   public:
      LanguageScores(size_t num_languages) ;
      LanguageScores(const LanguageScores *orig) ;
      LanguageScores(const LanguageScores *orig, double scale) ;
      ~LanguageScores() = default ;

      // accessors
      void *userData() const { return m_userdata ; }
      bool sorted() const { return m_sorted ; }
      unsigned numLanguages() const { return m_num_languages ; }
      unsigned maxLanguages() const { return m_max_languages ; }
      unsigned activeLanguage() const { return m_active_language ; }
      unsigned topLanguage() const { return m_info[0].id() ; }
      unsigned languageNumber(size_t N) const
	 { return (N < numLanguages()) ? m_info[N].id() : ~0 ; }
      double score(size_t N) const
	 { return (N < numLanguages()) ? m_info[N].score() : -1.0 ; }
      double highestScore() const ;
      unsigned highestLangID() const ;
      unsigned nonzeroScores() const ;

      // manipulators
      void setUserData(void *u) { m_userdata = u ; }
      void clear() ;
      void reserve(size_t N) ;
      void setScore(size_t N, double val)
	 { if (N < numLanguages()) m_info[N].setScore(val) ; }
      void increment(size_t N, double incr = 1.0)
	 { if (N < numLanguages()) m_info[N].incrScore(incr) ; }
      void decrement(size_t N, double decr = 1.0)
	 { if (N < numLanguages()) m_info[N].incrScore(-decr) ; }
      void scaleScore(size_t N, double scale_factor)
	 { if (N < numLanguages()) m_info[N].setScore(m_info[N].score() * scale_factor) ; }
      void scaleScores(double scale_factor) ;
      void sqrtScores() ;
      void add(const LanguageScores *scores, double weight = 1.0) ;
      void addThresholded(const LanguageScores *scores, double threshold,
			  double weight = 1.0) ;
      void subtract(const LanguageScores *scores, double weight = 1.0) ;
      bool lambdaCombineWithPrior(LanguageScores *prior, double lambda,
				  double smoothing) ;
      void filter(double cutoff_ratio) ;
      void sort(double cutoff_ratio = 0.0) ;
      void sort(double cutoff_ratio, unsigned max_langs) ;
      void mergeDuplicateNamesAndSort(const LanguageID *langinfo) ;
      void filterDuplicates(const class LanguageIdentifier *,
			    bool ignore_region = false) ;
      void setLanguage(unsigned lang)
	 { m_active_language = lang ; }

      // iterator support
      Info* begin() const { return m_info.begin() ; }
      Info* end() const { return m_info.begin() + numLanguages() ; }
   protected: // methods
      void sortByName(const LanguageID *langinfo) ;
      void invalidate() { m_info = nullptr ; m_num_languages = m_max_languages = 0 ; }

   private: // members
      Fr::NewPtr<Info>   m_info ;
      void*		 m_userdata ;
   protected:
      unsigned	 	 m_num_languages { 0 } ;
      unsigned		 m_max_languages { 0 } ;
      unsigned		 m_active_language ;
      bool      	 m_sorted ;
   } ;

//----------------------------------------------------------------------

class WeightedLanguageScores : public LanguageScores
   {
   public:
      WeightedLanguageScores(size_t num_languages,
			     double def_weight = 1.0) ;
      ~WeightedLanguageScores() = default ;

      // accessors
      double weight(size_t N) const
	 { return (N < numLanguages()) ? m_weights[N] : 0.0 ; }

      // manipulators
      void setWeight(size_t N, double wt)
	 { if (N < numLanguages()) m_weights[N] = wt ; }
      void incrWeight(size_t N, double wt)
	 { if (N < numLanguages()) m_weights[N] += wt ; }
      void sqrtWeights() ;

   private:
      Fr::DoublePtr	m_weights ;
   } ;

//----------------------------------------------------------------------

class LanguageIdentifier
   {
   public:
      static const uint32_t unknown_lang = (uint32_t)~0 ;
   public:
      LanguageIdentifier(const char *language_data_file,
			 bool verbose = false) ;
      ~LanguageIdentifier() ;

      // accessors
      bool good() const { return m_langdata && m_langdata->good() ; }
      bool verbose() const { return m_verbose ; }
      bool applyCoverageFactor() const { return m_apply_cover_factor && m_adjustments ; }
      size_t allocLanguages() const { return m_alloc_languages ; }
      size_t numLanguages() const { return m_num_languages ; }
      double adjustmentFactor(size_t N) const { return m_adjustments[N] ; }
      LanguageIdentifier *charsetIdentifier() const { return m_charsetident ; }
      class LangIDPackedMultiTrie* trie() const { return m_langdata.get() ; }
      LangIDPackedMultiTrie *packedTrie() ;
      class LangIDMultiTrie *unpackedTrie() ;
      const char *databaseLocation() const { return m_directory ; }
      const char *languageName(size_t N) const ;
      const char *friendlyName(size_t N) const ;
      const char *languageScript(size_t N) const ;
      const uint8_t *alignments() const { return m_alignments ; }
      Fr::CharPtr languageDescriptor(size_t N) const ;
      const char *languageEncoding(size_t N) const ;
      const char *languageSource(size_t N) const ;
      const LanguageID *languageInfo(size_t N) const
	 { return N < numLanguages() ? &m_langinfo[N] : nullptr ; }
      uint64_t trainingBytes(size_t N) const
	 { return N < numLanguages() ? m_langinfo[N].trainingBytes() : 0 ; }
      unsigned languageNumber(const LanguageID *lang_info) const ;
      unsigned languageNumber(const char *langdescript) const ;
      bool identify(LanguageScores *scores, const char *buffer,
		    size_t buflen, const uint8_t *alignments,
		    bool ignore_whitespace = false,
		    bool apply_stop_grams = true,
		    size_t length_normalization = 0) const ;
      LanguageScores *identify(const char *buffer, size_t buflen,
			       bool ignore_whitespace = false,
			       bool apply_stop_grams = true,
			       bool enforce_alignments = true) const ;
      LanguageScores *identify(LanguageScores *scores, /* may be NULL */
			       const char *buffer, size_t buflen,
			       bool ignore_whitespace = false,
			       bool apply_stop_grams = true,
			       bool enforce_alignments = true) const ;
      bool finishIdentification(LanguageScores *scores, unsigned select_highestN = 0,
				double cutoff_ratio = 0.1) const ;
      LanguageScores *similarity(unsigned langid) const ;
      bool sameLanguage(size_t L1, size_t L2,
			bool ignore_region = false) const ;
      double bigramWeight() const { return m_bigram_weight ; }

      // modifiers
      uint32_t addLanguage(const LanguageID &info, uint64_t train_bytes) ;
      void charsetIdentifier(LanguageIdentifier *id) 
	 { m_charsetident = (id ? id : this) ; }
      void setBigramWeight(double weight) { m_bigram_weight = weight ; }
      void useFriendlyName(bool friendly = true) { m_friendly_name = friendly ; }
      void runVerbosely(bool v) { m_verbose = v ; }
      void applyCoverageFactor(bool apply) { m_apply_cover_factor = apply ; }
      void incrStringCount(size_t langnum) ;
      bool computeSimilarities() ;

      // I/O
      static bool checkSignature(Fr::CFile& f, unsigned *version = nullptr) ;
      bool writeStatistics(Fr::CFile& f) const ;
      bool writeHeader(Fr::CFile& f) const ;
      bool write(Fr::CFile& f) ;
      bool write(const char* filename) const ;
      bool dump(Fr::CFile& f, bool show_ngrams = false) const ;

   private:
      void setAlignments() ;
      bool setAdjustmentFactors() ;

   private:
      Fr::Owned<LangIDPackedMultiTrie> m_langdata ;
      Fr::Owned<LangIDMultiTrie> m_uncomplangdata ;
      Fr::NewPtr<LanguageID> m_langinfo ;
      Fr::DoublePtr          m_length_factors ;
      Fr::DoublePtr          m_adjustments ;
      Fr::UInt8Ptr           m_alignments ;
      Fr::UInt8Ptr           m_unaligned ;
      Fr::NewPtr<size_t>     m_string_counts ;
      Fr::CharPtr            m_directory ;
      LanguageIdentifier*    m_charsetident ;
      double 	             m_bigram_weight ;
      size_t	             m_alloc_languages ;
      size_t 	             m_num_languages ;
      bool   	             m_friendly_name ;
      bool	             m_apply_cover_factor ;
      bool                   m_verbose ;
   } ;

/************************************************************************/
/*	Procedural interface						*/
/************************************************************************/

LanguageIdentifier *load_language_database(const char *database_file,
					   const char *charset_file,
					   bool create = false,
					   bool verbose = false) ;
   // set charset_file to NULL for default search, "" to not use a separate
   //    database (use the main database for charset ID as well as lang ID)
void unload_language_database(LanguageIdentifier *id) ;
double set_stopgram_penalty(double wt) ;

bool smooth_language_scores(bool smooth) ;
bool smoothing_language_scores() ;
LanguageScores *smoothed_language_scores(LanguageScores *scores,
					 LanguageScores *&prior_scores,
					 size_t match_length) ;
// the following uses a global, so is not thread-safe
LanguageScores *smoothed_language_scores(LanguageScores *scores,
					 size_t match_length) ;


#endif /* !__LANGID_H_INCLUDED */

// end of file langid.h //
