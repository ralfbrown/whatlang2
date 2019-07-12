/************************************************************************/
/*                                                                      */
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     langid.C							*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-12						*/
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
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <float.h>
#include <numeric>
#include <stdint.h>
#include "langid.h"
#include "mtrie.h"
#include "framepac/byteorder.h"
#include "framepac/config.h"
#include "framepac/message.h"
#include "framepac/texttransforms.h"

using namespace Fr ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#ifdef __GNUC__
# define likely(x) __builtin_expect((x),1)
# ifndef unlikely
#   define unlikely(x) __builtin_expect((x),0)
# endif /* !unlikely */
#else
# define likely(x) (x)
# define unlikely(x) (x)
#endif

#ifndef UINT32_MAX
# define UINT32_MAX		0xFFFFFFFFU
#endif

/************************************************************************/
/*	Types								*/
/************************************************************************/

class ScoreAndID
   {
   private:
      double   m_score ;
      unsigned short m_id ;
   public:
      ScoreAndID() {}
      void init(double sc, unsigned short new_id)
	 { m_score = sc ; m_id = new_id ; }

      // accessors
      double score() const { return m_score ; }
      unsigned short id() const { return m_id ; }

      // manipulators
      static void swap(ScoreAndID &, ScoreAndID &) ;

      // comparison
      static int compare(const ScoreAndID &, const ScoreAndID &) ;
      bool operator< (const ScoreAndID& other) const { return compare(*this,other) < 0 ; }
   } ;

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

//Fr::Allocator LanguageID::allocator("LanguageID",sizeof(LanguageID)) ;

//Fr::Allocator LanguageScores::allocator("LanguageScores",sizeof(LanguageScores)) ;

//Fr::Allocator WeightedLanguageScores::allocator("WtLanguageScores", sizeof(WeightedLanguageScores)) ;

static double stop_gram_penalty = -9.0 ;

/************************************************************************/
/*	Helper functions						*/
/************************************************************************/

// not in standard library until C++17
inline double clamp(double v, double lo, double hi)
{
   return std::min(std::max(v,lo),hi) ;
}

//----------------------------------------------------------------------

static uint8_t read_byte(Fr::CFile& f, uint32_t default_value = (uint8_t)~0)
{
   uint8_t valbuf ;
   return f.readValue(&valbuf) ? valbuf : default_value ;
}

//----------------------------------------------------------------------

static bool write_uint8(Fr::CFile& f, uint8_t value)
{
   uint8_t valbuf = value ;
   return f.writeValue(valbuf) ;
}

//----------------------------------------------------------------------

static uint32_t read_uint32(Fr::CFile& f, uint32_t default_value = (uint32_t)~0)
{
   UInt32 val ;
   return f.readValue(&val) ? val.load() : default_value ;
}

//----------------------------------------------------------------------

static bool write_uint32(Fr::CFile& f, uint32_t value)
{
   UInt32 val(value) ;
   return f.writeValue(val) ;
}

//----------------------------------------------------------------------

static CharPtr read_fixed_field(Fr::CFile& f, size_t len)
{
   if (len == 0)
      return nullptr ;
   CharPtr buf(len) ;
   if (!buf || f.read(*buf,len) < len)
      {
      return nullptr ;
      }
   // ensure proper string termination if the file didn't have a NUL
   buf[len-1] = '\0' ;
   return buf ;
}

//----------------------------------------------------------------------

static bool read_uint64(Fr::CFile& f, uint64_t &value)
{
   UInt64 val ;
   if (f.readValue(&val))
      {
      value = 0 ;
      return false ;
      }
   else
      {
      value = val.load() ;
      return true ;
      }
}

//----------------------------------------------------------------------

static bool write_fixed_field(Fr::CFile& f, const char *s, size_t len)
{
   if (len == 0)
      return false ;
   size_t string_len = s ? strlen(s) : 0 ;
   size_t count = std::min(len-1,string_len) ;
   if (f.write(s,count) < count)
      return false ;
   for (size_t i = string_len ; i+1 < len ; i++)
      {
      if (!f.putc('\0'))
	 return false ;
      }
   return f.putc('\0') ;
}

//----------------------------------------------------------------------

static bool write_uint64(Fr::CFile& f, uint64_t value)
{
   UInt64 val(value) ;
   return f.writeValue(val) ;
}

//----------------------------------------------------------------------

static void parse_language_description(const char *descript,
				       char *&language,
				       char *&region,
				       char *&encoding,
				       char *&source)
{
   // the format of a desciptor is lang_REG-encoding/source
   language = nullptr ;
   region = nullptr ;
   encoding = nullptr ;
   source = nullptr ;
   if (descript && *descript)
      {
      const char *underscore = strchr(descript,'_') ;
      const char *dash = strchr(descript,'-') ;
      const char *slash = strchr(descript,'/') ;
      const char *nul = strchr(descript,'\0') ;
      const char *lang_end = underscore ;
      if (!lang_end || (dash && dash < lang_end)) lang_end = dash ;
      if (!lang_end || (slash && slash < lang_end)) lang_end = slash ;
      if (!lang_end) lang_end = nul ;
      language = New<char>(lang_end - descript + 1) ;
      if (language)
	 {
	 memcpy(language,descript,lang_end - descript) ;
	 language[lang_end - descript] = '\0' ;
	 }
      if (underscore)
	 {
	 // we have a region specified
	 const char *reg_end = dash ;
	 if (!reg_end || (slash && slash < reg_end)) reg_end = slash ;
	 if (!reg_end) reg_end = nul ;
	 region = New<char>(reg_end - underscore) ;
	 if (region)
	    {
	    memcpy(region,underscore+1,reg_end - underscore) ;
	    region[reg_end - (underscore + 1)] = '\0' ;
	    }
	 }
      if (dash)
	 {
	 // we have an encoding specified
	 const char *enc_end = slash ;
	 if (!enc_end) enc_end = nul ;
	 encoding = New<char>(enc_end - dash) ;
	 if (encoding)
	    {
	    memcpy(encoding,dash+1,enc_end - dash) ;
	    encoding[enc_end - (dash + 1)] = '\0' ;
	    }
	 }
      if (slash)
	 {
	 // we have a source specified
	 source = New<char>(nul - slash) ;
	 if (source)
	    {
	    memcpy(source,slash+1,nul - slash) ;
	    source[nul - (slash + 1)] = '\0' ;
	    }
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static double length_factor(unsigned len)
{
   return 270.0 * ::pow(len,0.75) ;
}

//----------------------------------------------------------------------

static NewPtr<double> make_length_factors(unsigned max_length, double bigram_weight)
{
   if (max_length < 3)
      max_length = 3 ;
   NewPtr<double> factors(max_length+1) ;
   if (factors)
      {
      factors[0] = 0.0 ;
      factors[1] = 1.0 ;
      factors[2] = bigram_weight * length_factor(2) ;
      for (size_t i = 3 ; i <= max_length ; i++)
	 {
	 factors[i] = length_factor(i) ;
	 }
      }
   return factors ;
}

//----------------------------------------------------------------------

static double scale_score(uint32_t score)
{
//   double scaled = ::sqrt(::sqrt(score / (100.0 * TRIE_SCALE_FACTOR))) ;
   // smoothing is now precomputed in the language model database
   double scaled = 1.0 ;
   if (score & 1)
      {
      scaled = stop_gram_penalty ;
      score &= ~1 ;
      }
   return scaled * score / (100.0 * TRIE_SCALE_FACTOR) ;
}

/************************************************************************/
/*	Methods for class ScoreAndID					*/
/************************************************************************/

int ScoreAndID::compare(const ScoreAndID &s1, const ScoreAndID &s2)
{
   if (s1.score() < s2.score())
      return +1 ;
   else if (s1.score() > s2.score())
      return -1 ;
   else
      return 0 ;
}

//----------------------------------------------------------------------

void ScoreAndID::swap(ScoreAndID &s1, ScoreAndID &s2)
{
   double tmpscore = s1.score() ;
   s1.m_score = s2.score() ;
   s2.m_score = tmpscore ;
   unsigned tmpid = s1.id() ;
   s1.m_id = s2.id() ;
   s2.m_id = tmpid ;
   return ;
}

/************************************************************************/
/*	Methods for class BigramCounts					*/
/************************************************************************/

BigramCounts::BigramCounts(const BigramCounts *orig)
{
   copy(orig) ;
   return ;
}

//----------------------------------------------------------------------

BigramCounts::BigramCounts(Fr::CFile& f)
{
   if (!f || f.read(m_counts,sizeof(m_counts),1) < sizeof(m_counts))
      {
      // we could not read the data, so clear all the counts, effectively
      //   removing this bigram model from consideration
      std::fill_n(m_counts,lengthof(m_counts),0) ;
      }
   m_total = std::accumulate(m_counts,m_counts+lengthof(m_counts),0) ;
   return ;
}

//----------------------------------------------------------------------

void BigramCounts::copy(const BigramCounts *orig)
{
   if (orig)
      {
      std::copy_n(orig->m_counts,lengthof(m_counts),m_counts) ;
      m_total = orig->m_total ;
      }
   else
      {
      std::fill_n(m_counts,lengthof(m_counts),0) ;
      m_total = 0 ;
      }
   return ;
}

//----------------------------------------------------------------------

double BigramCounts::averageProbability(const char *buffer, size_t buflen)
const
{
   double prob = 0.0 ;
   if (buffer && buflen > 1)
      {
      for (size_t i = 1 ; i < buflen ; i++)
	 {
	 prob += this->count(buffer[i-1],buffer[i]) ;
	 }
      prob /= m_total ;
      prob /= (buflen - 1) ;
      }
   return prob ;
}

//----------------------------------------------------------------------

BigramCounts* BigramCounts::load(Fr::CFile& f)
{
   if (f)
      {
      NewPtr<BigramCounts> model(1) ;
      if (model->read(f))
	 return model.move() ;
      }
   return nullptr ;
}

//----------------------------------------------------------------------

bool BigramCounts::read(Fr::CFile& f)
{
   memset(m_counts,'\0',sizeof(m_counts)) ;
   m_total = 0 ;
   for (size_t c1 = 0 ; c1 <= 0xFF ; c1++)
      {
      for (size_t c2 = 0 ; c2 <= 0xFF ; c2++)
	 {
	 char line[FrMAX_LINE] ;
	 if (!f.gets(line,sizeof(line)))
	    return false ;
	 char *end ;
	 char *lineptr = skip_whitespace(line) ;
	 unsigned long cnt = strtoul(lineptr,&end,0) ;
	 if (end != lineptr)
	    {
	    m_total += cnt ;
	    set(c1,c2,cnt) ;
	    }
	 else
	    return false ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

bool BigramCounts::readBinary(Fr::CFile& f)
{
   return f && f.read(m_counts,lengthof(m_counts),sizeof(m_counts[0])) == lengthof(m_counts) ;
}

//----------------------------------------------------------------------

bool BigramCounts::dumpCounts(Fr::CFile& f) const
{
   if (!f)
      return false ;
   for (size_t c1 = 0 ; c1 <= 0xFF ; c1++)
      {
      for (size_t c2 = 0 ; c2 <= 0xFF ; c2++)
	 {
	 f << (size_t)count(c1,c2) ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

bool BigramCounts::save(Fr::CFile& f) const
{
   return f && f.write(m_counts,sizeof(m_counts),1) == sizeof(m_counts) ;
}

/************************************************************************/
/*	Methods for class LanguageID					*/
/************************************************************************/

LanguageID::LanguageID()
{
   clear() ;
   return ;
}

//----------------------------------------------------------------------

LanguageID::LanguageID(const char *lang, const char *reg, const char *enc,
		       const char *src, const char *scrpt)
{
   clear() ;
   setLanguage(lang) ;
   setRegion(reg) ;
   setEncoding(enc) ;
   setSource(src) ;
   setScript(scrpt) ;
   return ;
}

//----------------------------------------------------------------------

LanguageID::LanguageID(const LanguageID &orig)
{
   clear() ;
   setLanguage(orig.language(),orig.friendlyName()) ;
   setRegion(orig.region()) ;
   setEncoding(orig.encoding()) ;
   setSource(orig.source()) ;
   setScript(orig.script()) ;
   setAlignment(orig.alignment()) ;
   setCoverageFactor(orig.coverageFactor()) ;
   setCountedCoverage(orig.countedCoverage()) ;
   setFreqCoverage(orig.freqCoverage()) ;
   setMatchFactor(orig.matchFactor()) ;
   return ;
}

//----------------------------------------------------------------------

LanguageID::LanguageID(const LanguageID *orig)
{
   clear() ;
   if (orig)
      {
      setLanguage(orig->language(),orig->friendlyName()) ;
      setRegion(orig->region()) ;
      setEncoding(orig->encoding()) ;
      setSource(orig->source()) ;
      setScript(orig->script()) ;
      setAlignment(orig->alignment()) ;
      setCoverageFactor(orig->coverageFactor()) ;
      setCountedCoverage(orig->countedCoverage()) ;
      setFreqCoverage(orig->freqCoverage()) ;
      setMatchFactor(orig->matchFactor()) ;
      }
   return ;
}

//----------------------------------------------------------------------

LanguageID& LanguageID::operator = (LanguageID&& orig)
{
   m_language = orig.m_language.move() ;
   m_region = orig.m_region.move() ;
   m_encoding = orig.m_encoding.move() ;
   m_source = orig.m_source.move() ;
   m_script = orig.m_script.move() ;
   m_friendlyname = orig.m_friendlyname ;
   m_coverage = orig.m_coverage ;
   m_countcover = orig.m_countcover ;
   m_trainbytes = orig.m_trainbytes ;
   m_alignment = orig.m_alignment ;
   m_freqcover = orig.m_freqcover ;
   m_matchfactor = orig.m_matchfactor ;
   m_alignment = orig.m_alignment ;
   return *this ;
}

//----------------------------------------------------------------------

LanguageID::~LanguageID()
{
   setLanguage(nullptr) ;
   setRegion(nullptr) ;
   setEncoding(nullptr) ;
   setSource(nullptr) ;
   setScript(nullptr) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::clear()
{
   m_language = nullptr ;
   m_friendlyname = nullptr ;
   m_region = nullptr ;
   m_encoding = nullptr ;
   m_source = nullptr ;
   m_script = nullptr ;
   m_trainbytes = 0 ;
   m_alignment = 1 ;
   m_coverage = 0.0 ;
   m_countcover = 0.0 ;
   m_freqcover = 0.0 ;
   m_matchfactor = 0.0 ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setLanguage(const char *lang, const char *friendly)
{
   if (lang && friendly && lang != friendly)
      m_language = aprintf("%s=%s",lang,friendly) ;
   else
      m_language = dup_string(lang) ;
   m_friendlyname = m_language ;
   if (m_language)
      {
      char *eq = strchr(*m_language,'=') ;
      if (eq)
	 {
	 m_friendlyname = eq + 1 ;
	 *eq = '\0' ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setRegion(const char *reg)
{
   m_region = dup_string(reg) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setEncoding(const char *enc)
{
   m_encoding = dup_string(enc) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setSource(const char *src)
{
   m_source = dup_string(src) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setScript(const char *scr)
{
   m_script = dup_string(scr) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setAlignment(const char *align)
{
   int a = align ? atoi(align) : 1 ;
   if (a < 1)
      a = 1 ;
   else if (a > 4)
      a = 4 ;
   else if (a == 3)
      a = 2 ;
   setAlignment(a) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setCoverageFactor(double coverage)
{
   if (coverage <= 0.0 || coverage > 1.0)
      coverage = 1.0 ;
   m_coverage = coverage ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setCountedCoverage(double coverage)
{
   m_countcover = clamp(coverage,0.0,MAX_WEIGHTED_COVER) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setFreqCoverage(double coverage)
{
   m_freqcover = clamp(coverage,0.0,MAX_FREQ_COVER) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageID::setMatchFactor(double match)
{
   m_matchfactor = clamp(match,0.0,MAX_MATCH_FACTOR) ;
   return ;
}

//----------------------------------------------------------------------

bool LanguageID::guessScript()
{
   if (!script()
       || strcasecmp(script(),"UNKNOWN") == 0
       || strcmp(script(),"") == 0)
      {
      // try to guess the script from the encoding
      const char *enc = encoding() ;
      if (strcasecmp(enc,"iso-8859-6") == 0)
	 {
	 setScript("Arabic") ;
	 }
      else if (strcasecmp(enc,"ArmSCII8") == 0)
	 {
	 setScript("Armenian") ;
	 }
      else if (strncasecmp(enc,"KOI8-",5) == 0 ||
	  strcasecmp(enc,"KOI7") == 0 ||
	  strcasecmp(enc,"CP866") == 0 ||
	  strcasecmp(enc,"RUSCII") == 0 ||
	  strcasecmp(enc,"Windows-1251") == 0 ||
	  strcasecmp(enc,"iso-8859-5") == 0 ||
	  strcasecmp(enc,"Latin-5") == 0 ||
	  strcasecmp(enc,"MacCyrillic") == 0)
	 {
	 setScript("Cyrillic") ;
	 }
      else if (strcasecmp(enc,"ISCII") == 0)
	 {
	 setScript("Devanagari") ;
	 }
      else if (strcasecmp(enc,"iso-8859-7") == 0 ||
	       strcasecmp(enc,"cp737") == 0)
	 {
	 setScript("Greek") ;
	 }
      else if (strcasecmp(enc,"GB2312") == 0 ||
	       strcasecmp(enc,"GB-2312") == 0 ||
	       strcasecmp(enc,"GB18030") == 0 ||
	       strcasecmp(enc,"GBK") == 0 ||
	       strcasecmp(enc,"Big5") == 0 ||
	       strcasecmp(enc,"EUC-CN") == 0 ||
	       strcasecmp(enc,"EUC-TW") == 0)
	 {
	 setScript("Han") ;
	 }
      else if (strcasecmp(enc,"EUC-KR") == 0)
	 {
	 setScript("Hangul") ;
	 }
      else if (strcasecmp(enc,"CP862") == 0 ||
	       strncasecmp(enc,"iso-8859-8",10) == 0)
	 {
	 setScript("Hebrew") ;
	 }
      else if (strcasecmp(enc,"ShiftJIS") == 0 ||
	       strcasecmp(enc,"Shift-JIS") == 0 ||
	       strcasecmp(enc,"ISO-2022") == 0 ||
	       strcasecmp(enc,"EUC-JP") == 0 ||
	       strncasecmp(enc,"EUC-JIS",7) == 0)
	 {
	 setScript("Kanji") ;
	 }
      else if (strcasecmp(enc,"TIS620") == 0 ||
	       strcasecmp(enc,"TSCII") == 0 ||
	       strcasecmp(enc,"iso-8859-11") == 0)
	 {
	 setScript("Thai") ;
	 }
      else if (strcasecmp(enc,"VISCII") == 0)
	 {
	 setScript("Vietnamese") ;
	 }
      else if (strcasecmp(enc,"ASCII") == 0 ||
	       strcasecmp(enc,"CP437") == 0 ||
	  strncasecmp(enc,"ASCII-16",8) == 0 ||
	  strncasecmp(enc,"iso-8859-",9) == 0 ||
	  strncasecmp(enc,"Latin",5) == 0)
	 {
	 setScript("Latin") ;
	 }
      else if (!script()
	       || strcmp(script(),"") == 0)
	 {
	 setScript("UNKNOWN") ;
	 return false ;
	 }
      else
	 return false ;
      return true ;	// successfully made a guess
      }
   return true ;	// already have a script assigned
}

//----------------------------------------------------------------------

LanguageID* LanguageID::read(Fr::CFile& f, unsigned file_version)
{
   if (!f)
      return nullptr ;
   NewPtr<LanguageID> langID(1) ;
   if (!read(f,*langID,file_version))
      {
      langID = nullptr ;
      }
   return langID.move() ;
}

//----------------------------------------------------------------------

bool LanguageID::sameLanguage(const LanguageID &info, bool ignore_region) const
{
   const char *info_1, *info_2 ;
   info_1 = language() ;
   info_2 = info.language() ;
   if (info_1 != info_2 && (!info_1 || !info_2 || strcmp(info_1,info_2) != 0))
      return false ;
   if (!ignore_region)
      {
      info_1 = region() ;
      info_2 = info.region() ;
      if (info_1 != info_2 &&
	  (!info_1 || !info_2 || strcmp(info_1,info_2) != 0))
	 return false ;
      }
   info_1 = encoding() ;
   info_2 = info.encoding() ;
   if (info_1 != info_2 && (!info_1 || !info_2 || strcmp(info_1,info_2) != 0))
      return false ;
   return true ;
}

//----------------------------------------------------------------------

bool LanguageID::matches(const LanguageID *lang_info) const
{
   if (!lang_info)
      return false ;			// MUST specify a language
   if (strcasecmp(language(),lang_info->language()) != 0)
      return false ;
   const char *reg = lang_info->region() ;
   if (reg && *reg && region() && strcasecmp(region(),reg) != 0)
      return false ;
   const char *enc = lang_info->encoding() ;
   if (enc && *enc && encoding() && strcasecmp(encoding(),enc) != 0)
      return false ;
   const char *src = lang_info->source() ;
   if (src && *src && source() && strcasecmp(source(),src) != 0)
      return false ;
   return true ;
}

//----------------------------------------------------------------------

bool LanguageID::matches(const char *lang, const char *reg,
			 const char *enc, const char *src) const
{
   if (!lang || !*lang)
      return false ;			// MUST specify a language
   if (strcasecmp(language(),lang) != 0)
      return false ;
   if (reg && *reg && region() && strcasecmp(region(),reg) != 0)
      return false ;
   if (enc && *enc && encoding() && strcasecmp(encoding(),enc) != 0)
      return false ;
   if (src && *src && source() && strcasecmp(source(),src) != 0)
      return false ;
   return true ;
}

//----------------------------------------------------------------------

bool LanguageID::operator == (const LanguageID &info) const
{
   const char *info_1, *info_2 ;
   info_1 = language() ;
   info_2 = info.language() ;
   if (info_1 != info_2 && (!info_1 || !info_2 || strcmp(info_1,info_2) != 0))
      return false ;
   info_1 = region() ;
   info_2 = info.region() ;
   if (info_1 != info_2 && (!info_1 || !info_2 || strcmp(info_1,info_2) != 0))
      return false ;
   info_1 = encoding() ;
   info_2 = info.encoding() ;
   if (info_1 != info_2 && (!info_1 || !info_2 || strcmp(info_1,info_2) != 0))
      return false ;
   info_1 = source() ;
   info_2 = info.source() ;
   if (info_1 != info_2 && (!info_1 || !info_2 || strcmp(info_1,info_2) != 0))
      return false ;
   return true ;
}

//----------------------------------------------------------------------

bool LanguageID::read(Fr::CFile& f, LanguageID *langID, unsigned version)
{
   if (!langID)
      return false ;
   langID->m_language = read_fixed_field(f,LANGID_STRING_LENGTH) ;
   langID->m_region = read_fixed_field(f,LANGID_STRING_LENGTH) ;
   langID->m_encoding = read_fixed_field(f,LANGID_STRING_LENGTH) ;
   langID->m_source = read_fixed_field(f,LANGID_STRING_LENGTH) ;
   langID->m_script = read_fixed_field(f,LANGID_STRING_LENGTH) ;
   (void)read_uint64(f,langID->m_trainbytes ) ;
   uint8_t align = read_byte(f, 1) ;
   if (align < 1) align = 1 ;
   (void)read_byte(f,0) ;
   (void)read_byte(f,0) ;
   (void)read_byte(f,0) ;
   uint32_t cover = 0 ;
   cover = read_uint32(f,0) ;
   double cover_factor = ((double)cover) / (double)UINT32_MAX ;
   langID->setCoverageFactor(cover_factor) ;
   if (version > 4)
      {
      cover = read_uint32(f,0) ;
      langID->setCountedCoverage(cover * MAX_WEIGHTED_COVER / UINT32_MAX ) ;
      cover = read_uint32(f,0) ;
      langID->setFreqCoverage(cover * MAX_FREQ_COVER / UINT32_MAX ) ;
      cover = read_uint32(f,0) ;
      langID->setMatchFactor(cover * MAX_MATCH_FACTOR / UINT32_MAX ) ;
      }
   langID->m_friendlyname = langID->m_language ;
   langID->setAlignment(align) ;
   if (langID->m_language)
      {
      char *eq = strchr(*langID->m_language,'=') ;
      if (eq)
	 {
	 langID->m_friendlyname = eq + 1 ;
	 *eq = '\0' ;
	 }
      }
   return langID->m_language != nullptr && langID->m_encoding != nullptr ;
}

//----------------------------------------------------------------------

bool LanguageID::write(Fr::CFile& f) const
{
   if (!f)
      return false ;
   bool success = false ;
   if (m_friendlyname && m_friendlyname != m_language && 
       m_friendlyname > m_language)
      {
      ((char*)m_friendlyname)[-1] = '=' ;
      }
   double count_cover = m_countcover / MAX_WEIGHTED_COVER ;
   double freq_cover = m_freqcover / MAX_FREQ_COVER ;
   double match_factor = matchFactor() / MAX_MATCH_FACTOR ;
   // to simplify the version 1-5 file formats, we'll use fixed-size fields
   if (write_fixed_field(f,m_language,LANGID_STRING_LENGTH) &&
       write_fixed_field(f,m_region,LANGID_STRING_LENGTH) &&
       write_fixed_field(f,m_encoding,LANGID_STRING_LENGTH) &&
       write_fixed_field(f,m_source,LANGID_STRING_LENGTH) &&
       write_fixed_field(f,m_script,LANGID_STRING_LENGTH) &&
       write_uint64(f,m_trainbytes) &&
       write_uint8(f,m_alignment) &&
       write_uint8(f,0) &&
       write_uint8(f,0) &&
       write_uint8(f,0) &&
       write_uint32(f,(uint32_t)(m_coverage * UINT32_MAX) ) &&
       write_uint32(f,(uint32_t)(count_cover * UINT32_MAX) ) &&
       write_uint32(f,(uint32_t)(freq_cover * UINT32_MAX) ) &&
       write_uint32(f,(uint32_t)(match_factor * UINT32_MAX) ))
      success = true ;
   if (m_friendlyname && m_friendlyname != m_language && 
       m_friendlyname > m_language)
      ((char*)m_friendlyname)[-1] = '\0' ;
   return success ;
}

/************************************************************************/
/*	Methods for class LanguageScores				*/
/************************************************************************/

LanguageScores::LanguageScores(size_t num_languages)
   : m_lang_ids(num_languages), m_scores(num_languages)
{
   m_sorted = false ;
   m_userdata = nullptr ;
   m_active_language = 0 ;
   if (m_lang_ids && m_scores)
      {
      m_num_languages = num_languages ;
      m_max_languages = num_languages ;
      std::fill_n(m_scores.begin(),num_languages,0.0) ;
      std::iota(m_lang_ids.begin(),m_lang_ids.begin()+num_languages,0) ;
      }
   else
      {
      invalidate() ;
      }
   return ;
}

//----------------------------------------------------------------------

LanguageScores::LanguageScores(const LanguageScores *orig)
{
   if (orig)
      {
      unsigned nlang = orig->numLanguages() ;
      m_num_languages = nlang ;
      m_lang_ids = NewPtr<unsigned short>(nlang) ;
      m_scores = NewPtr<double>(nlang) ;
      if (m_lang_ids && m_scores)
	 {
	 m_sorted = orig->m_sorted ;
	 std::copy_n(orig->m_scores.begin(),nlang,m_scores.begin()) ;
	 std::copy_n(orig->m_lang_ids.begin(),nlang,m_lang_ids.begin()) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

LanguageScores::LanguageScores(const LanguageScores *orig, double scale)
{
   if (orig)
      {
      unsigned nlang = orig->numLanguages() ;
      m_num_languages = nlang ;
      m_lang_ids = NewPtr<unsigned short>(nlang) ;
      m_scores = NewPtr<double>(nlang) ;
      if (m_lang_ids && m_scores)
	 {
	 m_sorted = orig->m_sorted ;
	 for (size_t i = 0 ; i < nlang ; i++)
	    {
	    m_scores[i] = orig->score(i) * scale ;
	    }
	 std::copy_n(orig->m_lang_ids.begin(),nlang,m_lang_ids.begin()) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

double LanguageScores::highestScore() const
{
   size_t nlang = numLanguages() ;
   if (sorted())
      {
      return m_scores[0] ;
      }
   else if (nlang == 0)
      {
      return 0.0 ;
      }
   else
      {
      return *std::max_element(m_scores.begin(),m_scores.begin()+nlang) ;
      }
}

//----------------------------------------------------------------------

unsigned LanguageScores::highestLangID() const
{
   auto nlang = numLanguages() ;
   if (sorted())
      {
      return m_lang_ids[0] ;
      }
   else if (nlang == 0)
      return (unsigned)~0 ;
   else
      {
      double* highest = std::max_element(m_scores.begin(),m_scores.begin()+nlang) ;
      return (highest - m_scores.begin()) ;
      }
}

//----------------------------------------------------------------------

unsigned LanguageScores::nonzeroScores() const
{
   if (sorted())
      {
      return (m_scores[0] > LANGID_ZERO_SCORE) ? numLanguages() : 0 ;
      }
   else
      {
      size_t count = 0 ;
      for (size_t i = 0 ; i < numLanguages() ; i++)
	 {
	 if (m_scores[i] > LANGID_ZERO_SCORE)
	    count++ ;
	 }
      return count ;
      }
}

//----------------------------------------------------------------------

void LanguageScores::clear()
{
   m_num_languages = maxLanguages() ;
   std::fill_n(m_scores.begin(),numLanguages(),0.0) ;
   m_sorted = false ;
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::reserve(size_t N)
{
   m_scores = NewPtr<double>(N) ;
   m_lang_ids = NewPtr<unsigned short>(N) ;
   m_num_languages = m_max_languages = N ;
   m_sorted = false ;
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::scaleScores(double scale_factor)
{
   for (size_t i = 0 ; i < numLanguages() ; i++)
      {
      m_scores[i] *= scale_factor ;
      }
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::sqrtScores()
{
   std::transform(m_scores.begin(),m_scores.begin()+numLanguages(),m_scores.begin(),::sqrt) ;
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::add(const LanguageScores *scores, double weight)
{
   if (scores && weight != 0)
      {
      size_t count = numLanguages() ;
      if (scores->numLanguages() < count)
	 count = scores->numLanguages() ;
      for (size_t i = 0 ; i < count ; i++)
	 {
	 m_scores[i] += (scores->m_scores[i] * weight) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::addThresholded(const LanguageScores *scores,
				    double threshold, double weight)
{
   if (scores && weight != 0)
      {
      size_t count = numLanguages() ;
      if (scores->numLanguages() < count)
	 count = scores->numLanguages() ;
      for (size_t i = 0 ; i < count ; i++)
	 {
	 double sc = scores->m_scores[i] ;
	 if (sc >= threshold)
	    m_scores[i] += (sc * weight) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::subtract(const LanguageScores *scores, double weight)
{
   if (scores && weight != 0)
      {
      size_t count = numLanguages() ;
      if (scores->numLanguages() < count)
	 count = scores->numLanguages() ;
      for (size_t i = 0 ; i < count ; i++)
	 {
	 m_scores[i] -= (scores->m_scores[i] * weight) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

bool LanguageScores::lambdaCombineWithPrior(LanguageScores *prior, double lambda,
					    double smoothing)
{
   size_t count = numLanguages() ;
   if (prior && prior->numLanguages())
      {
      for (size_t i = 0 ; i < count ; i++)
	 {
	 double priorscore = prior->m_scores[i] ;
	 double currscore = m_scores[i] ;
	 if (currscore >= LANGID_ZERO_SCORE)
	    prior->m_scores[i] += currscore * smoothing ;
	 m_scores[i] = lambda * currscore + (1.0 - lambda) * priorscore ;
	 }
      return true ;
      }
   return false ;
}

//----------------------------------------------------------------------

void LanguageScores::sort(double cutoff_ratio)
{
   if (!sorted() && numLanguages() > 0)
      {
      NewPtr<ScoreAndID> scores_and_ids(numLanguages()) ;
      if (!scores_and_ids)
	 {
	 //!!!
	 return ;
	 }
      double cutoff = LANGID_ZERO_SCORE ;
      if (cutoff_ratio > 0.0)
	 {
	 if (cutoff_ratio > 1.0)
	    cutoff_ratio = 1.0 ;
	 double threshold = highestScore() * cutoff_ratio ;
	 if (threshold > cutoff)
	    cutoff = threshold ;
	 }
      unsigned num_scores = 0 ;
      for (unsigned i = 0 ; i < numLanguages() ; i++)
	 {
	 if (m_scores[i] >= cutoff)
	    scores_and_ids[num_scores++].init(m_scores[i],m_lang_ids[i]) ;
	 }
      if (num_scores > 0)
	 {
	 if (num_scores > 1)
	    {
	    std::sort(*scores_and_ids,(*scores_and_ids)+num_scores) ;
	    }
	 for (unsigned i = 0 ; i < num_scores ; i++)
	    {
	    m_scores[i] = scores_and_ids[i].score() ;
	    m_lang_ids[i] = scores_and_ids[i].id() ;
	    }
	 m_num_languages = num_scores ;
	 }
      else
	 {
	 // nothing is above our cutoff, but we can't just discard everything,
	 //   so scan for the highest score and make that the sole score
	 for (unsigned i = 1 ; i < numLanguages() ; i++)
	    {
	    if (m_scores[i] > m_scores[0])
	       {
	       m_scores[0] = m_scores[i] ;
	       m_lang_ids[0] = m_lang_ids[i] ;
	       }
	    }
	 m_num_languages = 1 ;
	 }
      m_sorted = true ;
      }
   return ;
}

//----------------------------------------------------------------------

static void insert(double score, unsigned lang_id,
		   ScoreAndID *scores_and_ids,
		   unsigned &num_scores,
		   unsigned max_scores)
{
   for (size_t i = 0 ; i < num_scores ; i++)
      {
      if (score > scores_and_ids[i].score())
	 {
	 if (num_scores < max_scores)
	    num_scores++ ;
	 for (size_t j = num_scores - 1 ; j > i ; j--)
	    {
	    scores_and_ids[j] = scores_and_ids[j-1] ;
	    }
	 scores_and_ids[i].init(score,lang_id) ;
	 return ;
	 }
      }
   if (num_scores == 0)
      {
      scores_and_ids[0].init(score,lang_id) ;
      num_scores++ ;
      }
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::sort(double cutoff_ratio, unsigned max_langs)
{
   if (max_langs == 0 || max_langs > 10 || max_langs >= numLanguages())
      sort(cutoff_ratio) ;
   else if (!sorted() && numLanguages() > 0)
      {
      NewPtr<ScoreAndID> scores_and_ids(max_langs) ;
      if (!scores_and_ids)
	 {
	 //!!!
	 return ;
	 }
      double cutoff = LANGID_ZERO_SCORE ;
      unsigned num_scores = 0 ;
      for (unsigned i = 0 ; i < numLanguages() ; i++)
	 {
	 if (m_scores[i] >= cutoff)
	    {
	    insert(m_scores[i],m_lang_ids[i],*scores_and_ids,num_scores,max_langs) ;
	    double threshold = scores_and_ids[0].score() * cutoff_ratio ;
	    if (threshold > cutoff)
	       {
	       cutoff = threshold ;
	       }
	    }
	 }
      if (num_scores > 0)
	 {
	 for (unsigned i = 0 ; i < num_scores ; i++)
	    {
	    m_scores[i] = scores_and_ids[i].score() ;
	    m_lang_ids[i] = scores_and_ids[i].id() ;
	    }
	 m_num_languages = num_scores ;
	 }
      else
	 {
	 // nothing is above our cutoff, but we can't just discard everything,
	 //   so scan for the highest score and make that the sole score
	 for (unsigned i = 1 ; i < numLanguages() ; i++)
	    {
	    if (m_scores[i] > m_scores[0])
	       {
	       m_scores[0] = m_scores[i] ;
	       m_lang_ids[0] = m_lang_ids[i] ;
	       }
	    }
	 m_num_languages = 1 ;
	 }
      m_sorted = true ;
      }
   return ;
}

//----------------------------------------------------------------------

// the following variable makes sortByName() non-threadsafe
static const LanguageID *sort_langinfo = nullptr ;

static int compare_names(const ScoreAndID& s1, const ScoreAndID& s2)
{
   if (!sort_langinfo)
      return s1.id() - s2.id() ;
   const char *name1 = sort_langinfo[s1.id()].language() ;
   const char *name2 = sort_langinfo[s2.id()].language() ;
   if (name1 && name2)
      {
      return strcmp(name1,name2) ;
      }
   else if (!name2)
      return -1 ;
   else if (!name1)
      return +1 ;
   return 0 ; // items are identical by sort order
}

//----------------------------------------------------------------------

void LanguageScores::sortByName(const LanguageID *langinfo)
{
   if (numLanguages() > 0 && langinfo != nullptr)
      {
      NewPtr<ScoreAndID> scores_and_ids(numLanguages()) ;
      if (!scores_and_ids)
	 {
	 //!!!
	 return ;
	 }
      unsigned num_scores = 0 ;
      for (unsigned i = 0 ; i < numLanguages() ; i++)
	 {
	 if (m_scores[i] > 0.0)
	    scores_and_ids[num_scores++].init(m_scores[i],m_lang_ids[i]) ;
	 }
      sort_langinfo = langinfo ;
      std::sort(*scores_and_ids,(*scores_and_ids)+num_scores,compare_names) ;
      for (unsigned i = 0 ; i < num_scores ; i++)
	 {
	 m_scores[i] = scores_and_ids[i].score() ;
	 m_lang_ids[i] = scores_and_ids[i].id() ;
	 }
      m_num_languages = num_scores ;
      }
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::mergeDuplicateNamesAndSort(const LanguageID *langinfo)
{
   if (!langinfo)
      return ;
   sortByName(langinfo) ;
   for (size_t i = 0 ; i + 1 < numLanguages() ; i++)
      {
      const char *name1 = langinfo[languageNumber(i)].language() ;
      if (!name1 || !*name1)
	 continue ;
      for (size_t j = i + 1 ; j < numLanguages() ; j++)
	 {
	 const char *name2 = langinfo[languageNumber(j)].language() ;
	 if (!name2 || !*name2)
	    break ;
	 if (strcmp(name1,name2) == 0)
	    {
	    m_scores[i] += m_scores[j] ;
	    m_scores[j] = 0.0 ;
	    }
	 }
      }
   sort() ;
   return ;
}

//----------------------------------------------------------------------

void LanguageScores::filterDuplicates(const LanguageIdentifier *langid,
				      bool ignore_region)
{
   if (!langid)
      return ;
   unsigned dest = 1 ;
   for (size_t i = 1 ; i < numLanguages() ; i++)
      {
      bool is_dup = false ;
      for (size_t j = 0 ; j < i ; j++)
	 {
	 if (langid->sameLanguage(i,j,ignore_region))
	    {
	    is_dup = true ;
	    break ;
	    }
	 }
      if (!is_dup)
	 {
	 m_lang_ids[dest] = m_lang_ids[i] ;
	 m_scores[dest] = m_scores[i] ;
	 dest++ ;
	 }
      }
   m_num_languages = dest ;
   return ;
}

/************************************************************************/
/*	Methods for class WeightedLanguageScores			*/
/************************************************************************/

WeightedLanguageScores::WeightedLanguageScores(size_t num_languages,
					       double default_weight)
   : LanguageScores(num_languages), m_weights(num_languages)
{
   if (m_weights)
      {
      std::fill_n(*m_weights,numLanguages(),default_weight) ;
      }
   else
      {
      m_num_languages = 0 ;
      }
   return ;
}

//----------------------------------------------------------------------

void WeightedLanguageScores::sqrtWeights()
{
   std::transform(*m_weights,*m_weights+numLanguages(),*m_weights,::sqrt) ;
   return ;
}

/************************************************************************/
/*	Methods for class LanguageIdentifier				*/
/************************************************************************/

LanguageIdentifier::LanguageIdentifier(const char *language_data_file,
				       bool run_verbosely)
{
   m_alloc_languages = 0 ;
   m_num_languages = 0 ;
   m_apply_cover_factor = true ;
   useFriendlyName(false) ;
   charsetIdentifier(nullptr) ;
   setBigramWeight(DEFAULT_BIGRAM_WEIGHT) ;
   runVerbosely(run_verbosely) ;
   if (language_data_file && *language_data_file)
      {
      CInputFile fp(language_data_file,CFile::binary) ;
      if (fp)
	 {
	 unsigned version = 0 ;
	 if (checkSignature(fp,&version))
	    {
	    Fr::FilePath path(language_data_file) ;
	    m_directory = Fr::dup_string(path.directory()) ;
	    m_num_languages = read_uint32(fp,0) ;
	    if (numLanguages() > 0)
	       {
	       m_langinfo = NewPtr<LanguageID>(numLanguages()) ;
	       if (m_langinfo)
		  {
		  m_alloc_languages = numLanguages() ;
		  }
	       else
		  {
		  m_num_languages = 0 ;
		  }
	       uint8_t have_bigrams = false ;
	       (void)fp.readValue(&have_bigrams) ;
	       // skip the reserved padding
	       fp.seek(LANGID_PADBYTES_1,SEEK_CUR) ;
	       // read the language info records
	       for (size_t i = 0 ; i < numLanguages() ; i++)
		  {
		  if (!LanguageID::read(fp,&m_langinfo[i],version))
		     {
		     m_num_languages = 0 ;
		     break ;
		     }
		  }
	       // next, read the multi-trie
	       if (m_num_languages > 0)
		  {
		  m_langdata = LangIDPackedMultiTrie::load(fp,language_data_file) ;
		  if (m_langdata)
		     {
		     m_alloc_languages = numLanguages() ;
		     }
		  }
	       // finally, load the data mapping, if present
//FIXME
	       }
	    }
	 }
      }
   if (!PackedTrieFreq::dataMappingInitialized())
      {
      PackedTrieFreq::initDataMapping(scale_score) ;
      }
   setAlignments() ;
   setAdjustmentFactors() ;
   if (!m_langdata)
      m_langdata = new LangIDPackedMultiTrie ;
   if (!m_langinfo)
      m_langinfo = NewPtr<LanguageID>(1) ;
   m_string_counts = NewPtr<size_t>(numLanguages()) ;
   if (m_string_counts)
      std::fill_n(m_string_counts.begin(),numLanguages(),0) ;
   if (m_langdata)
      m_length_factors = make_length_factors(m_langdata->longestKey(),m_bigram_weight) ;
   return ;
}

//----------------------------------------------------------------------

LanguageIdentifier::~LanguageIdentifier()
{
   m_num_languages = 0 ;
   m_alloc_languages = 0 ;
   return ;
}

//----------------------------------------------------------------------

void LanguageIdentifier::setAlignments()
{
   m_alignments = NewPtr<uint8_t>(PACKED_TRIE_LANGID_MASK + 1) ;
   for (size_t i = 0 ; i < numLanguages() ; i++)
      {
      m_alignments[i] = languageInfo(i)->alignment() ;
      }
   for (size_t i = numLanguages() ; i <= PACKED_TRIE_LANGID_MASK ; i++)
      {
      m_alignments[i] = (uint8_t)~0 ;
      }
   if (!m_unaligned)
      {
      m_unaligned = NewPtr<uint8_t>(PACKED_TRIE_LANGID_MASK + 1) ;
      std::fill_n(m_unaligned.begin(),numLanguages(),1) ;
      for (size_t i = numLanguages() ; i <= PACKED_TRIE_LANGID_MASK ; i++)
	 {
	 m_unaligned[i] = (uint8_t)~0 ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::setAdjustmentFactors()
{
   m_adjustments = NewPtr<double>(numLanguages()) ;
   if (!m_adjustments)
      return false ;
   for (size_t i = 0 ; i < numLanguages() ; i++)
      {
      const LanguageID *lang_info = languageInfo(i) ;
      m_adjustments[i] = 1.0 ;
      if (lang_info)
	 {
//!!!	 double cover = lang_info->coverageFactor() ;
//	 double cover = lang_info->countedCoverage() ;
	 double cover = lang_info->matchFactor() ;
	 if (cover > 0.0)
	    {
	    cover = ::pow(cover,0.25) ;
	    double align = 1.0 ;
	    // adjust for the fact that an enforced alignment of
	    //   greater than 1 forces the match factor to be lower
	    //   since only 1/align bytes could possibly start a match
	    if (m_alignments && m_alignments[i] <= 8)
	       align = m_alignments[i] ;
	    m_adjustments[i] = align / cover ;
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

LangIDPackedMultiTrie* LanguageIdentifier::packedTrie()
{
   if (!m_langdata && m_uncomplangdata)
      {
      m_langdata = new LangIDPackedMultiTrie(m_uncomplangdata) ;
      m_uncomplangdata = nullptr ;
      }
   return m_langdata.get() ;
}

//----------------------------------------------------------------------

LangIDMultiTrie* LanguageIdentifier::unpackedTrie()
{
   if (!m_uncomplangdata && m_langdata)
      {
      m_uncomplangdata = new LangIDMultiTrie(m_langdata) ;
      m_langdata = nullptr ;
      }
   return m_uncomplangdata.get() ;
}

//----------------------------------------------------------------------

const char* LanguageIdentifier::languageName(size_t N) const
{
   if (N < numLanguages())
      {
      return (m_friendly_name
	      ? m_langinfo[N].friendlyName()
	      : m_langinfo[N].language()) ;
      }
   return nullptr ;
}

//----------------------------------------------------------------------

const char *LanguageIdentifier::friendlyName(size_t N) const
{
   return (N < numLanguages()) ? m_langinfo[N].friendlyName() : nullptr ;
}

//----------------------------------------------------------------------

const char *LanguageIdentifier::languageEncoding(size_t N) const
{
   return (N < numLanguages()) ? m_langinfo[N].encoding() : nullptr ;
}

//----------------------------------------------------------------------

const char *LanguageIdentifier::languageSource(size_t N) const
{
   return (N < numLanguages()) ? m_langinfo[N].source() : nullptr ;
}

//----------------------------------------------------------------------

const char *LanguageIdentifier::languageScript(size_t N) const
{
   if (N < numLanguages())
      {
      return m_langinfo[N].script() ? m_langinfo[N].script() : "UNKNOWN" ;
      }
   return nullptr ;
}

//----------------------------------------------------------------------

Fr::CharPtr LanguageIdentifier::languageDescriptor(size_t N) const
{
   if (N < numLanguages())
      {
      return aprintf("%s_%s-%s",
			m_langinfo[N].language(),
			m_langinfo[N].region(),
			m_langinfo[N].encoding()) ;
      }
   return nullptr ;
}

//----------------------------------------------------------------------

unsigned LanguageIdentifier::languageNumber(const LanguageID *lang_info)
   const
{
   unsigned modelnum = (unsigned)~0 ;
   if (lang_info)
      {
      // scan the list of models in the identifier for a uniquely-matching
      //   model
      for (size_t i = 0 ; i < numLanguages() ; i++)
	 {
	 const LanguageID *info = languageInfo(i) ;
	 if (info && info->matches(lang_info))
	    {
	    if (modelnum != (unsigned)~0)
	       {
	       modelnum = (unsigned)~0 ;
	       if (verbose())
		  {
		  cerr << "Multiple models match language specifier" << endl ;
		  }
	       break ;
	       }
	    modelnum = i ;
	    }
	 }
      }
   return modelnum ;
}

//----------------------------------------------------------------------

unsigned LanguageIdentifier::languageNumber(const char *langdescript) const
{
   unsigned modelnum = (unsigned)~0 ;
   if (langdescript)
      {
      // parse the description into language, region, encoding, and source
      char *language = nullptr ;
      char *region = nullptr ;
      char *encoding = nullptr ;
      char *source = nullptr ;
      parse_language_description(langdescript,language,region,encoding,source);
      // scan the list of models in the identifier for a uniquely-matching
      //   model
      for (size_t i = 0 ; i < numLanguages() ; i++)
	 {
	 const LanguageID *info = languageInfo(i) ;
	 if (info && info->matches(language,region,encoding,source))
	    {
	    if (modelnum != (unsigned)~0)
	       {
	       if (verbose())
		  {
		  cerr << "Multiple models match language specifier "
		       << langdescript << endl ;
		  }
	       modelnum = (unsigned)~0 ;
	       break ;
	       }
	    modelnum = i ;
	    }
	 }
      Free(language) ;
      Free(region) ;
      Free(encoding) ;
      Free(source) ;
      }
   return modelnum ;
}

//----------------------------------------------------------------------

static const unsigned max_alignments[4] = { 4, 1, 2, 1 } ;

static void identify_languages(const char *buffer, size_t buflen,
                               const LangIDPackedMultiTrie *langdata,
			       LanguageScores *scores,
			       const uint8_t *alignments,
			       const double *length_factors,
			       bool apply_stop_grams,
			       size_t length_normalizer)
{
   //assert(scores != nullptr) ;
   unsigned minhist = length_factors[2] ? 1 : 2 ;
   double *score_array = scores->scoreArray() ;
   double normalizer = (double)length_normalizer ;
   for (size_t index = 0 ; index + minhist < buflen ; index++)
      {
      uint32_t nodeindex = LangIDPackedMultiTrie::ROOT_INDEX ;
      if ((nodeindex = langdata->extendKey((uint8_t)buffer[index],nodeindex)) == LangIDPackedMultiTrie::NULL_INDEX)
	 continue ;
      if (minhist > 1 &&
	  (nodeindex = langdata->extendKey((uint8_t)buffer[index+1],nodeindex)) == LangIDPackedMultiTrie::NULL_INDEX)
	 continue ;
      // we have character sets with alignments of 1, 2, or 4 bytes; the
      //   low two bits of the offset from the start of the buffer tells
      //   us the maximum alignment which is valid at this point
      unsigned max_alignment = max_alignments[index%4] ;
      // since we'll almost always fail to extend the key before hitting
      //   the longest key in the trie, we can avoid conditional assignments
      //   and extra math by simply trying to extend the key all the way to
      //   the end of the buffer
      for (size_t i = index + minhist ; i < buflen ; i++)
	 {
	 uint8_t keybyte = (uint8_t)buffer[i] ;
	 if ((nodeindex = langdata->extendKey(keybyte,nodeindex)) == LangIDPackedMultiTrie::NULL_INDEX)
	    break ;
	 // check whether we're at a leaf node; if so, add all of the
	 //   frequencies to the scores
	 auto node = langdata->node(nodeindex) ;
	 if (node->leaf())
	    {
	    double len_factor = length_factors[i - index + 1] ;
	    const PackedTrieFreq *f = node->frequencies(langdata->frequencyBaseAddress()) ;
	    // normalize by text length so that scores are
	    //   comparable between different buffer sizes
	    len_factor /= normalizer ;
	    if (apply_stop_grams)
	       {
	       do {
		  unsigned id = f->languageID() ;
		  // ignore mis-aligned ngrams; we avoid a check that
		  //   'id' is in range by setting all possible IDs
		  //   above the number of models in the database such
		  //   that the alignment check never succeeds
		  if (likely(alignments[id] <= max_alignment))
		     {
		     double prob = f->mappedScore() ;
		     score_array[id] += (prob * len_factor) ;
		     }
		  f++ ;
	          } while (!f[-1].isLast()) ;
	       }
	    else
	       {
	       do {
		  unsigned id = f->languageID() ;
		  // ignore mis-aligned ngrams; we avoid a check that
		  //   'id' is in range by setting all possible IDs
		  //   above the number of models in the database such
		  //   that the alignment check never succeeds
		  if (likely(alignments[id] <= max_alignment))
		     {
		     double prob = f->mappedScore() ;
		     if (unlikely(prob <= 0.0))
			break ;		// only stopgrams from here on
		     score_array[id] += (prob * len_factor) ;
		     }
		  f++ ;
	          } while (!f[-1].isLast()) ;
	       }
	    }
	 }
      }
   return ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::identify(LanguageScores *scores,
				  const char *buffer, size_t buflen,
				  const uint8_t *alignments_,
				  bool ignore_whitespace,
				  bool apply_stop_grams,
				  size_t length_normalization) const
{
   if (!buffer || !scores || !m_langdata)
      return false ;
   if (scores && scores->maxLanguages() == numLanguages())
      {
      scores->clear() ;
      }
   else
      {
      delete[] scores ;
      scores = new LanguageScores(numLanguages()) ;
      }
   m_langdata.get()->ignoreWhiteSpace(ignore_whitespace) ;
   if (m_length_factors)
      m_length_factors[2] = m_bigram_weight * length_factor(2) ;
   if (!alignments_)
      alignments_ = m_unaligned ;
   if (length_normalization == 0)
      length_normalization = buflen ;
   identify_languages(buffer,buflen,m_langdata,scores,alignments_,
		      m_length_factors,apply_stop_grams,
		      length_normalization) ;
   m_langdata.get()->ignoreWhiteSpace(false) ;
   return true ;
}

//----------------------------------------------------------------------

LanguageScores* LanguageIdentifier::identify(const char* buffer,
					     size_t buflen,
					     bool ignore_whitespace,
					     bool apply_stop_grams,
					     bool enforce_alignment) const
{
   if (!buffer || !buflen || !m_langdata)
      return nullptr ;
   NewPtr<LanguageScores> scores(new LanguageScores(numLanguages())) ;
   const auto align = enforce_alignment ? m_alignments.begin() : nullptr ;
   if (!identify(scores.begin(),buffer,buflen,align,ignore_whitespace, apply_stop_grams,0))
      {
      scores = nullptr ;
      }
   return scores.move() ;
}

//----------------------------------------------------------------------

LanguageScores *LanguageIdentifier::identify(LanguageScores *scores,
					     const char *buffer,
					     size_t buflen,
					     bool ignore_whitespace,
					     bool apply_stop_grams,
					     bool enforce_alignment) const
{
   if (!buffer || !buflen || !m_langdata)
      return nullptr ;
   if (scores && scores->maxLanguages() == numLanguages())
      {
      scores->clear() ;
      }
   else
      {
      delete[] scores ;
      scores = new LanguageScores(numLanguages()) ;
      }
   const auto align = enforce_alignment ? m_alignments.get() : nullptr ;
   if (!identify(scores,buffer,buflen,align,ignore_whitespace,
		 apply_stop_grams,0))
      {
      delete[] scores ;
      scores = nullptr ;
      }
   return scores ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::finishIdentification(LanguageScores *scores, unsigned highestN,
					      double cutoff_ratio) const
{
   if (!scores)
      return false ;
   if (applyCoverageFactor())
      {
      for (size_t i = 0 ; i < scores->numLanguages() ; i++)
	 {
	 scores->setScore(i,scores->score(i) * adjustmentFactor(scores->languageNumber(i))) ;
	 }
      }
   if (highestN > 0)
      {
      if (highestN > scores->numLanguages())
	 highestN = scores->numLanguages() ;
      scores->sort(cutoff_ratio,highestN) ;
      }
   return true ;
}

//----------------------------------------------------------------------

static bool cosine_term(const PackedTrieNode *node, const uint8_t *,
			unsigned /*keylen*/, void *user_data)
{
   auto scores = (WeightedLanguageScores*)user_data ;
   double lang1prob = 0.0 ;
   auto frequencies = node->frequencies((PackedTrieFreq*)scores->userData()) ;
   unsigned langid = scores->activeLanguage() ;
   for (auto freq = frequencies ; freq ; freq = freq->next())
      {
      if (freq->languageID() == langid)
	 {
	 if (!freq->isStopgram())
	    lang1prob = freq->probability() ;
	 break ;
	 }
      }
   for (auto freq = frequencies ; freq ; freq = freq->next())
      {
      unsigned lang2 = freq->languageID() ;
      if (!freq->isStopgram())
	 {
	 double lang2prob = freq->probability() ;
	 scores->incrWeight(lang2,lang2prob * lang2prob) ;
	 scores->increment(lang2,lang1prob * lang2prob) ;
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

LanguageScores *LanguageIdentifier::similarity(unsigned langid) const
{
   if (langid >= numLanguages() || !trie())
      return nullptr ;
   auto scores = new WeightedLanguageScores(numLanguages(),0.0) ;
   if (scores && trie())
      {
      scores->setLanguage(langid) ;
      scores->setUserData((void*)trie()->frequencyBaseAddress()) ;
      auto maxkey = trie()->longestKey() ;
      LocalAlloc<uint8_t,512> keybuf(maxkey+1) ;
      trie()->enumerate(keybuf,maxkey,cosine_term,scores) ;
      scores->sqrtWeights() ;
      for (size_t i = 0 ; i < numLanguages() ; i++)
	 {
	 double wt = scores->weight(i) * scores->weight(langid) ;
	 if (wt > 0.0)
	    scores->setScore(i,(scores->score(i)/ wt)) ;
	 }
      }
   return scores ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::computeSimilarities()
{
//FIXME
   return true ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::sameLanguage(size_t L1, size_t L2,
				      bool ignore_region) const
{
   if (L1 < numLanguages() && L2 < numLanguages())
      {
      return m_langinfo[L1].sameLanguage(m_langinfo[L2],ignore_region) ;
      }
   return false ;
}

//----------------------------------------------------------------------

uint32_t LanguageIdentifier::addLanguage(const LanguageID &info,
					 uint64_t train_bytes)
{
   for (size_t i = 0 ; i < numLanguages() ; i++)
      {
      if (m_langinfo[i] == info)
	 return i ;
      }
   if (numLanguages() >= allocLanguages())
      {
      auto nlang = numLanguages() ;
      size_t new_alloc = nlang ? 2 * nlang : 120 ;
      NewPtr<LanguageID> new_info(new_alloc) ;
      if (!new_info)
	 return unknown_lang ;
      if (nlang)
	 std::move(m_langinfo.begin(),m_langinfo.begin()+nlang,new_info.begin()) ;
      m_langinfo = new_info ;
      m_alloc_languages = new_alloc ;
      }
   uint32_t langID = m_num_languages++ ;
   new (&m_langinfo[langID]) LanguageID(&info) ;
   m_langinfo[langID].setTraining(train_bytes) ;
   return langID ;
}

//----------------------------------------------------------------------

void LanguageIdentifier::incrStringCount(size_t langnum)
{
   if (m_string_counts && langnum < numLanguages())
      m_string_counts[langnum]++ ;
   return ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::checkSignature(Fr::CFile& f, unsigned *file_version)
{
   if (file_version)
      *file_version = 0 ;
   const unsigned siglen = sizeof(LANGID_FILE_SIGNATURE) ;
   char buffer[siglen] ;
   if (f.read(buffer,siglen,sizeof(char)) != siglen)
      {
      errno = EACCES ;
      return false ;
      }
   if (memcmp(buffer,LANGID_FILE_SIGNATURE,siglen) != 0)
      {
      errno = EINVAL ;
      return false ;
      }
   unsigned char version ;
   if (!f.readValue(&version) || version < LANGID_MIN_FILE_VERSION || version > LANGID_FILE_VERSION)
      {
      errno = EINVAL ;
      return false ;
      }
   if (file_version)
      *file_version = version ;
   return true ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::writeStatistics(Fr::CFile& f) const
{
   if (!f || !m_string_counts)
      return false ;
   f.printf("===================\n") ;
   f.printf("Number of strings extracted, by language:\n") ;
   LanguageScores counts(numLanguages()) ;
   for (size_t i = 0 ; i < numLanguages() ; i++)
      {
      counts.setScore(i,m_string_counts[i]) ;
      }
   counts.mergeDuplicateNamesAndSort(m_langinfo) ;
   for (size_t i = 0 ; i < numLanguages() ; i++)
      {
      double count = counts.score(i) ;
      if (count <= 0.0)
	 break ;
      unsigned langnum = counts.languageNumber(i) ;
      f.printf(" %7lu\t%s\n",(unsigned long)count,languageName(langnum)) ;
      }
   f.printf("===================\n") ;
   return true ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::writeHeader(Fr::CFile& f) const
{
   // write the signature string
   const unsigned siglen = sizeof(LANGID_FILE_SIGNATURE) ;
   if (f.write(LANGID_FILE_SIGNATURE,siglen) != siglen)
      return false ;
   // follow with the format version number
   unsigned char version = LANGID_FILE_VERSION ;
   if (f.write(&version,sizeof(version),1) != sizeof(version))
      return false ;
   // then the number of language info records
   uint32_t n_lang = numLanguages() ;
      if (!f.writeValue(n_lang))
      return false ;
   // set the flag for whether we have bigram models following the multi-trie
   uint8_t have_bigrams = 1 ;
   if (f.writeValue(have_bigrams))
      return false ;
   // pad the header with NULs for the unused reserved portion of the header
   return f.putNulls(LANGID_PADBYTES_1) ;
}

//----------------------------------------------------------------------

#if FIXME
static int compare_frequencies(const MultiTrieFrequency &f1,
			       const MultiTrieFrequency &f2)
{
   bool s1 = f1.isStopgram() ;
   bool s2 = f2.isStopgram() ;
   // all stopgrams go after non-stopgrams
   if (s1 && !s2)
      return +1 ;
   else if (!s1 && s2)
      return -1 ;
   // within stopgram/non-stopgram, sort by language ID for cache locality
   uint32_t id1 = f1.languageID() ;
   uint32_t id2 = f2.languageID() ;
   if (id1 < id2)
      return -1 ;
   else
      return +1 ;
}
#endif

//----------------------------------------------------------------------

static bool sort_frequencies(const LangIDMultiTrie* trie, NybbleTrie::NodeIndex nodeindex, const uint8_t * /*key*/,
			     unsigned /*keylen*/, void * /*user_data*/)
{
   auto node = trie->node(nodeindex) ;
   auto f = node->frequencies() ;
   if (f)
      {
      // sort the frequency records
//FIXME      f = std::stable_sort(f,compare_frequencies) ;
      if (!node->setFrequencies(f))
	 return false ;
      }
   return true ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::write(Fr::CFile& f)
{
   if (!f)
      return false ;
   bool success = writeHeader(f) ;
   if (success)
      {
      // sort the frequency records for each leaf node so that stop-grams
      //   come last
      LangIDMultiTrie *mtrie = unpackedTrie() ;
      uint8_t keybuf[500] ;
      if (mtrie &&
	  !mtrie->enumerate(keybuf,sizeof(keybuf),sort_frequencies,mtrie))
	 {
	 success = false ;
	 }
      // write out the languageID records
      for (size_t i = 0 ; i < numLanguages() ; i++)
	 {
	 m_langinfo[i].write(f) ;
	 }
      // now write out the trie
      auto ptrie = packedTrie() ;
      if (!ptrie || !ptrie->write(f))
	 {
	 success = false ;
	 }
      uint32_t all_ones = (uint32_t)~0 ;
      if (!f.writeValue(all_ones))
	 success = false ;
      // finally, write out the mapping from stored frequency value to
      //   actual weighted value
      uint64_t dm_offset = f.tell() ;
      if (PackedTrieFreq::writeDataMapping(f))
	 {
	 f.seek(LANGID_FILE_DMOFFSET) ;
	 if (!f.writeValue(dm_offset))
	    {
	    success = false ;
	    }
	 }
      if (success)
	 f.writeComplete() ;
      }
   return success ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::write(const char *filename) const
{
   COutputFile fp(filename,CFile::safe_rewrite) ;
   bool success = const_cast<LanguageIdentifier*>(this)->write(fp) ;
   return success ? fp.close() : false ;
}

//----------------------------------------------------------------------

bool LanguageIdentifier::dump(Fr::CFile& f, bool show_ngrams) const
{
   f.printf("LanguageIdentifier Begin\n") ;
   for (size_t i = 0 ; i < numLanguages() ; i++)
      {
      const LanguageID *info = &m_langinfo[i] ;
      f.printf("  Lang %2u: %s_%s-%s / %s\n",(unsigned)i,
	      info->language(),info->region(),info->encoding(),info->source());
      }
   bool success = true ;
   if (m_langdata && show_ngrams)
      {
      f.printf("LanguageIdentifier Trie\n") ;
      success = m_langdata->dump(f) ;
      }
   else if (m_uncomplangdata && show_ngrams)
      {
      f.printf("LanguageIdentifier Trie\n") ;
      success = m_uncomplangdata->dump(f) ;
      }
   f.printf("LanguageIdentifier End\n") ;
   return success ;
}

/************************************************************************/
/*	Procedural interface						*/
/************************************************************************/

static LanguageIdentifier* try_loading(const char* database_file, bool verbose)
{
   if (!database_file)
      return nullptr ;
   CharPtr db_filename ;
   if (database_file[0] == '~' && database_file[1] == '/')
      {
      const char *home = getenv("HOME") ;
      if (home)
	 {
	 db_filename = aprintf("%s%s",home,database_file+1) ;
	 }
      else
	 {
	 const char *user = getenv("USER") ;
	 if (user)
	    db_filename = aprintf("/home/%s%s",user,database_file+1) ;
	 }
      }
   if (!db_filename)
      db_filename = dup_string(database_file) ;
   auto id = new LanguageIdentifier(db_filename,verbose) ;
   if (!id)
      {
      SystemMessage::no_memory("loading language database") ;
      return nullptr ;
      }
   else if (id->numLanguages() == 0)
      {
      if (verbose)
	 cerr << "Unsuccessfully tried to open " << database_file << endl ;
      delete id ;
      return nullptr ;
      }
   else if (verbose)
      {
      cerr << "Opened language database " << database_file << endl ;
      }
   return id ;
}

//----------------------------------------------------------------------

LanguageIdentifier *load_language_database(const char *database_file,
					   const char *charset_file,
					   bool create, bool verbose)
{
   LanguageIdentifier *id = nullptr ;
   if (database_file && *database_file)
      id = try_loading(database_file, verbose) ;
   if (!id && !create)
      {
      id = try_loading(FALLBACK_LANGID_DATABASE, verbose) ;
      if (!id)
	 {
	 id = try_loading(ALTERNATE_LANGID_DATABASE, verbose) ;
	 if (!id)
	    id = try_loading(DEFAULT_LANGID_DATABASE, verbose) ;
	 }
      }
   if (!id && create)
      id = new LanguageIdentifier(database_file,verbose) ;
   if (!id)
      {
      if (database_file && *database_file)
	 {
	 cerr << "Warning: Unable to load database from " << database_file
	      << endl ;
	 }
      else
	 {
	 cerr << "Warning: unable to load database from standard locations"
	      << endl ;
	 }
      }
   else
      {
      LanguageIdentifier *cs = nullptr ;
      if (charset_file)
	 {
	 if (*charset_file)
	    cs = try_loading(charset_file, verbose) ;
	 else
	    cs = id ;
	 }
      if (!cs)
	 {
	 cs = try_loading(FALLBACK_CHARSET_DATABASE, verbose) ;
	 if (!cs)
	    {
	    cs = try_loading(ALTERNATE_CHARSET_DATABASE, verbose) ;
	    if (!cs)
	       cs = try_loading(DEFAULT_CHARSET_DATABASE, verbose) ;
	    }
	 }
      if (!cs)
	 cs = id ;
      id->charsetIdentifier(cs) ;
      }
   return id ;
}

//----------------------------------------------------------------------

void unload_language_database(LanguageIdentifier *id)
{
   if (!id)
      return ;
   if (id->charsetIdentifier() != id)
      {
      delete id->charsetIdentifier() ;
      id->charsetIdentifier(nullptr) ;
      }
   delete id ;
   return ;
}

//----------------------------------------------------------------------

double set_stopgram_penalty(double pen)
{
   double old_pen = stop_gram_penalty ;
   stop_gram_penalty = -10.0 * pen ;
   return old_pen ;
}

// end of file langid.C //
