/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LangIdent: long n-gram-based language identification		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     whatlang.C  main program/wrapper for simple identifier	*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-23 						*/
/*                                                                      */
/*  (c) Copyright 2011,2012,2013,2014,2019				*/
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

#include "langid.h"
#include "framepac/config.h"
#include "framepac/file.h"
#include "framepac/message.h"
#include "framepac/texttransforms.h"
#include "framepac/unicode.h"

using namespace Fr ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define DEFAULT_TOPN 3

#define DEFAULT_BLOCKSIZE 4096
#define FULL_FILE_BLOCKSIZE (256*1024)
#define MIN_BLOCKSIZE 80
#define BY_LINE_BLOCKSIZE (65536U)

#define CUTOFF_RATIO 0.8

#define VERSION "1.30"

/************************************************************************/
/*	Types for this Module						*/
/************************************************************************/

enum LineMode {
   LM_None,
   LM_8bit,
   LM_16bigendian,
   LM_16littleendian
   } ;

/************************************************************************/
/*	Global Variables						*/
/************************************************************************/

static bool terse_language = false ;
static bool verbose = false ;
static bool show_script = false ;

static double bigram_weight = DEFAULT_BIGRAM_WEIGHT ;

/************************************************************************/
/************************************************************************/

static void usage(const char *argv0)
{
   fprintf(stderr,
	   "WhatLang v" VERSION "  Copyright 2011,2012,2019 Ralf Brown/CMU -- GNU GPLv3\n"
	   "Usage: %s [flags] [file]\n"
	   "Flags:\n"
	   "  -h     show this usage summary\n"
	   "  -b0    make single identification for entire file\n"
	   "  -b1    identify languages line by line\n"
	   "  -bN    set block size to N bytes (default 4096)\n"
	   "  -f     use full (friendly) language name in terse mode\n"
	   "  -lF    use language identification database in file F\n"
	   "  -nN    output at most N guesses for the language of a block\n"
	   "  -rR    don't output languages scoring less than R times highest\n"
	   "  -s     show scores of multiple sources for a language (if present)\n"
	   "  -t     terse -- output only language name, not full description\n"
	   "  -v     verbose -- show all blocks, even if no language detected\n"
           "  -WSPEC set internal scoring weights according to SPEC:\n"
           "         b0.1,s1.5  would set bigram weights to 0.1 and stopgram weights\n"
           "                    to 1.5\n"
,
	   argv0) ;
   exit(1) ;
}

//----------------------------------------------------------------------

static bool same_language(const char *name1, const char *name2)
{
   if (!name1 || !name2)
      return false ;
   return strcmp(name1,name2) == 0 ;
}

//----------------------------------------------------------------------

static void write_as_UTF8(CFile& f, const char *buf, int buflen, LineMode line_mode)
{
   if (line_mode == LM_None || line_mode == LM_8bit)
      {
      f.write(buf,buflen) ;
      }
   else
      {
      auto line = (const unsigned char*)buf ;
      for (int i = 0 ; i < buflen ; i += 2)
	 {
	 char utf8[6] ;
	 wchar_t codepoint
	    = (line_mode == LM_16bigendian) 
	    ? ((line[i] << 8) | line[i+1])
	    : ((line[i+1] << 8) | line[i]) ;
	 bool byteswap = false ;
	 int bytes = Fr::Unicode_to_UTF8(codepoint,utf8,byteswap) ;
	 if (bytes < 0 && i + 2 < buflen)
	    {
	    i += 2 ;
	    wchar_t codepoint2
	       = (line_mode == LM_16bigendian) 
	       ? ((line[i] << 8) | line[i+1])
	       : ((line[i+1] << 8) | line[i]) ;
	    bytes = Fr::Unicode_surrogates_to_UTF8(codepoint,codepoint2,utf8,byteswap) ;
	    }
	 if (bytes > 0)
	    f.write(utf8,bytes) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static void identify(const char *buf, int buflen, 
		     const LanguageIdentifier &langid,
		     size_t offset, unsigned topN, double cutoff_ratio,
		     bool separate_sources, bool full_file,
		     LineMode line_mode)
{
   if (!buf || buflen == 0)
      return ;
   LanguageScores *rawscores = langid.identify(buf,buflen) ;
   langid.finishIdentification(rawscores) ;
   Owned<LanguageScores> scores = langid.smoothedScores(rawscores,buflen) ;
   if (!scores)
      return ;
   unsigned num_scores = langid.numLanguages() ;
   bool echo_text = (line_mode != LM_None) ;
   if (topN > num_scores)
      topN = num_scores ;
   scores->sort(cutoff_ratio,2*topN) ;
   if (!separate_sources && !(terse_language && echo_text))
      scores->filterDuplicates(&langid) ;
   double highest_score = scores->score(0) ;
   CFile out(stdout) ;
   if (highest_score > LANGID_ZERO_SCORE)
      {
      if (!full_file && !echo_text)
	 out.printf("@ %8.08lX-%8.08lX ",offset,offset+buflen-1) ;
      unsigned shown = 0 ;
      double threshold = highest_score * cutoff_ratio ;
      for (size_t i = 0 ; i < scores->numLanguages() && shown < topN ; i++)
	 {
	 double sc = scores->score(i) ;
	 if (sc <= LANGID_ZERO_SCORE || sc < threshold)
	    break ;
	 unsigned langnum = scores->languageNumber(i) ;
	 Fr::CharPtr langdesc ;
	 if (terse_language)
	    langdesc = Fr::dup_string(langid.languageName(langnum)) ;
	 else
	    langdesc = langid.languageDescriptor(langnum) ;
	 const char *sep = "" ;
	 const char *source = "" ;
	 if (separate_sources)
	    {
	    const char *src = langid.languageSource(langnum) ;
	    if (src && *src)
	       {
	       sep = "/" ;
	       source = src ;
	       }
	    }
	 if (terse_language && echo_text)
	    {
	    const char *langname
	       = langid.languageName(scores->languageNumber(i)) ;
	    for (size_t j = 0 ; j < i ; j++)
	       {
	       if (same_language(langid.languageName(scores->languageNumber(j)),
		     langname))
		  {
		  langname = nullptr ;
		  break ;
		  }
	       }
	    if (langname)
	       {
	       if (shown > 0)
		  {
		  out.putc(',') ;
		  }
	       out.puts(*langdesc) ;
	       if (show_script)
		  out.printf("@%s",langid.languageScript(langnum)) ;
	       shown++ ;
	       }
	    }
	 else
	    {
	    if (shown > 0)
	       {
	       out.putc(' ') ;
	       }
	    if (show_script)
	       out.printf("%s%s%s@%s:%f",*langdesc,sep,source,langid.languageScript(langnum),sc) ;
	    else
	       out.printf("%s%s%s:%f",*langdesc,sep,source,sc) ;
	    shown++ ;
	    }
	 }
      if (echo_text)
	 {
	 out.putc('\t') ;
	 write_as_UTF8(out,buf,buflen,line_mode) ;
	 }
      else
	 out.printf("\n") ;
      out.flush() ;
      }
   else if (echo_text)
      {
      out.puts("??\t") ;
      write_as_UTF8(out,buf,buflen,line_mode) ;
      }
   else if (verbose)
      {
      out.printf("@ %8.08lX-%8.08lX: no languages detected\n", offset,offset+buflen-1) ;
      }
   return ;
}

//----------------------------------------------------------------------

static const char* locate_newline(const char *buf, int buflen, LineMode line_mode)
{
   if (line_mode == LM_8bit)
      {
      const char *newline = (const char*)memchr(buf,'\n',buflen) ;
      return newline ? newline + 1 : nullptr ;
      }
   else if (line_mode == LM_16bigendian)
      {
      for (int i = 0 ; i < buflen ; i += 2)
	 {
	 if (buf[i] == '\0' && buf[i+1] == '\n')
	    return buf + i + 2 ;
	 }
      return nullptr ;
      }
   else // if (line_mode == LM_16littleendian)
      {
      for (int i = 0 ; i < buflen ; i += 2)
	 {
	 if (buf[i+1] == '\0' && buf[i] == '\n')
	    return buf + i + 2 ;
	 }
      return nullptr ;
      }
}

//----------------------------------------------------------------------

static void identify_languages(CFile& f,
			       const LanguageIdentifier &langid,
			       int blocksize, unsigned topN,
			       double cutoff_ratio, bool separate_sources,
			       LineMode line_mode)
{
   int overlap = blocksize / 4 ;
   int bufsize = blocksize < FULL_FILE_BLOCKSIZE ? 2*blocksize : blocksize ;
   CharPtr bufbase(bufsize) ;
   if (!bufbase)
      {
      fprintf(stderr,"Out of memory\n") ;
      return ;
      }
   char *highwater = (bufsize > blocksize
		      ? *bufbase + bufsize - blocksize
		      : *bufbase + bufsize) ;
   int buflen = f.read(*bufbase,bufsize) ;
   char *buf = *bufbase ;
   size_t offset = 0 ;
   while (buflen > 0)
      {
      int check_size = buflen > blocksize ? blocksize : buflen ;
      if (line_mode != LM_None)
	 {
	 const char *nextline = locate_newline(buf,buflen,line_mode) ;
	 if (nextline)
	    check_size = (nextline - buf) ;
	 }
      identify(buf,check_size,langid,offset,topN,cutoff_ratio,separate_sources,
	       blocksize >= FULL_FILE_BLOCKSIZE,line_mode) ;
      if (blocksize >= FULL_FILE_BLOCKSIZE)
	 {
	 break ;     // only do one block if "entire file" chosen as blocksize
	 }
      // slide the window forward by 3/4 the block size
      unsigned shift = (line_mode == LM_None) ? overlap : check_size ;
      buf += shift ;
      buflen -= shift ;
      if (buf >= highwater)
	 {
	 unsigned to_read = (buf - bufbase) ;
	 offset += to_read ;
	 std::copy_n(buf,buflen,bufbase.begin()) ;
	 buf = *bufbase ;
	 int additional = f.read(buf + buflen,to_read) ;
	 // stop if we've already identified up to the end of the file, to
	 //   prevent a small orphan block with inaccurate identification
	 if (additional <= 0 && line_mode == LM_None)
	    break ;
	 buflen += additional ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

static void identify_languages(const char *filename,
			       const LanguageIdentifier &langid,
			       int blocksize, unsigned topN,
			       double cutoff_ratio, bool separate_sources,
			       bool show_filename, LineMode line_mode)
{
   CInputFile fp(filename) ;
   if (fp)
      {
      if (show_filename)
	 printf("File %s\n",filename) ;
      identify_languages(fp,langid,blocksize,topN,cutoff_ratio,separate_sources,line_mode) ;
      }
   else
      {
      fprintf(stderr,"Unable to open '%s' for reading\n",filename) ;
      }
   return ;
}

//----------------------------------------------------------------------

static void set_weight(char type, double weight)
{
   switch (type)
      {
      case 'b':
	 bigram_weight = (weight >= 0.0 ? weight : DEFAULT_BIGRAM_WEIGHT) ;
	 break ;
      case 's':
	 set_stopgram_penalty(weight) ;
	 break ;
      default:
	 SystemMessage::error("Unknown weight type '%c' in -W argument",type) ;
	 break ;
      }
   return ;
}

//----------------------------------------------------------------------

static void parse_weights(const char *wtspec)
{
   if (!wtspec)
      return ;
   while (*wtspec)
      {
      while (*wtspec == ',')
	 wtspec++ ;		     // allow empty specs to simplify scripts
      char type = *wtspec++ ;
      char *spec_end ;
      double value = strtod(wtspec,&spec_end) ;
      if (spec_end && spec_end > wtspec)
	 {
	 set_weight(type,value) ;
	 wtspec = spec_end ;
	 if (*wtspec == ',')
	    wtspec++ ;
	 else
	    break ;
	 }
      else
	 break ;
      }
   return ;
}

//----------------------------------------------------------------------

int main(int argc, char **argv)
{
   unsigned topN = DEFAULT_TOPN ;
   int blocksize = DEFAULT_BLOCKSIZE ;
   double cutoff_ratio = CUTOFF_RATIO ;
   bool separate_sources = false ;
   bool apply_coverage = false ;
   bool use_friendly_name = false ;
   LineMode line_mode = LM_None ;
   LineMode line_type = LM_8bit ;
   const char *argv0 = argv[0] ;
   const char *language_db = nullptr ;

   while (argc > 1 && argv[1][0] == '-')
      {
      switch (argv[1][1])
	 {
	 case '1':
	    {
	    if (argv[1][2] == '6')
	       line_type = (toupper(argv[1][3]) == 'B') ? LM_16bigendian : LM_16littleendian ;
	    }
	    break ;
	 case '8':
	    line_type = LM_8bit ;
	    break ;
	 case 'A':
	    show_script = true ;
	    break ;
	 case 'b':
	    blocksize = atoi(argv[1]+2) ;
	    break ;
	 case 'C':
	    apply_coverage = !apply_coverage ;
	    break ;
	 case 'f':
	    use_friendly_name = true ;
	    break ;
	 case 'l':
	    language_db = argv[1]+2 ;
	    break ;
	 case 'n':
	    topN = atoi(argv[1]+2) ;
	    break ;
	 case 'r':
	    cutoff_ratio = strtod(argv[1]+2,nullptr) ;
	    break ;
	 case 's':
	    separate_sources = true ;
	    break ;
	 case 't':
	    terse_language = true ;
	    break ;
	 case 'v':
	    verbose = true ;
	    break ;
	 case 'W':
	    parse_weights(argv[1]+2) ;
	    break ;
	 default:
	    fprintf(stderr,"Unknown option '%s'\n",argv[1]) ;
	    /* FALLTHROUGH */
	 case 'h':
	    usage(argv0) ;
	    break ;
	 }
      argc-- ;
      argv++ ;
      }
   if (topN < 1)
      topN = 1 ;
   if (cutoff_ratio < 0.0001)
      cutoff_ratio = 0.0001 ;
   if (cutoff_ratio > 1.0)
      cutoff_ratio = 1.0 ;
   if (blocksize == 2)
      {
      line_mode = line_type ;
      blocksize = BY_LINE_BLOCKSIZE ;
      }
   else if (blocksize == 1)
      {
      line_mode = line_type ;
      blocksize = BY_LINE_BLOCKSIZE ;
      }
   else if (blocksize == 0 || blocksize > FULL_FILE_BLOCKSIZE)
      {
      blocksize = FULL_FILE_BLOCKSIZE ;
      }
   else if (blocksize < MIN_BLOCKSIZE)
      {
      blocksize = MIN_BLOCKSIZE ;
      fprintf(stderr,
	      "Specified block size is ridiculously small, adjusted to %d\n",
	      MIN_BLOCKSIZE) ;
      }
   auto langid = LanguageIdentifier::load(language_db, "", false, verbose) ;
   if (!langid)
      return 1 ;
   langid->setBigramWeight(bigram_weight) ;
   langid->applyCoverageFactor(apply_coverage) ;
   langid->useFriendlyName(use_friendly_name) ;
   langid->smoothScores(blocksize == 2) ;
   if (argc == 1)
      {
      // no filename specified on command line, so use stdin
      CFile in(stdin) ;
      identify_languages(in,*langid,blocksize,topN,cutoff_ratio,separate_sources,line_mode) ;
      }
   else
      {
      bool multiple_files = (argc > 2) ;
      for (int i = 1 ; i < argc ; i++)
	 {
	 identify_languages(argv[i],*langid,blocksize,topN,cutoff_ratio,separate_sources,multiple_files,line_mode) ;
	 }
      }
   return 0 ;
}

// end of file whatlang.C //
