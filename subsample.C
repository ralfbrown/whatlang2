/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LangIdent: long n-gram-based language identification		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     subsample.C  utility program to sample lines of text	*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-15 						*/
/*                                                                      */
/*  (c) Copyright 2012,2013,2014,2019 Carnegie Mellon University	*/
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
#include <cstring>
#include <iostream>
#include "framepac/file.h"
#include "framepac/random.h"

using namespace std ;
using namespace Fr ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define MAX_LINE 32768U

/************************************************************************/
/*	Types								*/
/************************************************************************/

class StringList
   {
   public:
      StringList(CharPtr strng) ;
      ~StringList() = default ;

      // accessors
      StringList* next() const { return m_next ; }
      const char* string() const { return m_string ; }
      unsigned length() const { return m_length ; }
      unsigned listlength() const ;

      // modifiers
      void next(StringList* nxt) { m_next = nxt ; }

      // factory
      static void append(CharPtr strng, StringList**& last_string) ;

   private:
      StringList* m_next { nullptr } ;
      CharPtr     m_string ;
      unsigned	  m_length ;
   } ;

/************************************************************************/
/*	Methods for class StringList					*/
/************************************************************************/

StringList::StringList(CharPtr strng) : m_string(strng)
{
   m_length = m_string ? strlen(m_string) : 0 ;
   return ;
}

//----------------------------------------------------------------------

unsigned StringList::listlength() const
{
   unsigned count = 0 ;
   for (const StringList *s = this ; s ; s = s->next())
      {
      count++ ;
      }
   return count ;
}

//----------------------------------------------------------------------

void StringList::append(CharPtr strng, StringList**& last_string)
{
   if (!last_string)
      return ;				// no list, so we can't do anything
   auto str = new StringList(strng) ;
   if (!str)
      return ;				// Out of memory!  Ignore line
   *last_string = str ;
   StringList** end_ptr = &str->m_next ;
   last_string = end_ptr ;
   return ;
}

/************************************************************************/
/************************************************************************/

static void usage(const char *argv0)
{
   cerr << "Usage: " << argv0 << " [options] count <inputfile" << endl ;
   cerr << "Extract 'count' lines from the input file and write them to "
           "standard output.\n"
           "Options:\n"
           "\t-b\tdeterministically sample approx. 'count' bytes in total\n"
           "\t-i\tsample every 'count'th line ('count' is interval, not total)\n"
           "\t-lX\tsample lines more than 'X' bytes in length ('count' ignored)\n"
           "\t-LX\tsample lines less than 'X' bytes in length\n"
           "\t-u\tsample uniformly-spaced lines\n"
           "\t-oF\twrite sampled lines to file F (default is standard output)\n"
           "\t-rF\twrite non-sampled (rejected) lines to file F\n"
           "\t-R\tgenerate the same 'random' sample every time\n"
	<< endl ;
   exit(1) ;
}

//----------------------------------------------------------------------

static void write_line(CFile& f, const char* line)
{
   if (f)
      {
      f.puts(line) ;
      f.putc('\n') ;
      }
   return ;
}
  
//----------------------------------------------------------------------

static bool select_line(bool selected, CFile& selectfp, CFile& rejectfp, const char* line)
{
   if (selected)
      write_line(selectfp,line) ;
   else
      write_line(rejectfp,line) ;
   return selected ;
}


//----------------------------------------------------------------------

static void take_uniform_bytes(StringList* lines, size_t sample_size, CFile& selectfp, CFile& rejectfp)
{
   if (!lines)
      return ;
   size_t total_bytes = 0 ;
   for (const StringList *l = lines ; l ; l = l->next())
      {
      total_bytes += l->length() ;
      }
   double sample_rate = (sample_size + 1.0) / (double)total_bytes ;
   size_t sampled_bytes = 0 ;
   total_bytes = 0 ;
   while (lines)
      {
      Owned<StringList> line = lines ;
      lines = lines->next() ;
      size_t len = line->length() ;
      if (select_line(sampled_bytes <= (total_bytes * sample_rate),selectfp,rejectfp,line->string()))
	 sampled_bytes += len ;
      total_bytes += len ;
      }
   return ;
}

//----------------------------------------------------------------------

static void take_uniform_sample(StringList* lines, size_t sample_size, CFile& selectfp, CFile& rejectfp)
{
   if (!lines)
      return ;
   size_t numlines = lines->listlength() ;
   double interval = sample_size / (double)numlines ;
   double count = interval/2.0 ;
   while (lines)
      {
      Owned<StringList> line = lines ;
      lines = lines->next() ;
      select_line((size_t)(count + interval) > (size_t)count,selectfp,rejectfp,line->string()) ;
      count += interval ;
      }
   return ;
}

//----------------------------------------------------------------------

static void take_random_sample(StringList* lines, size_t sample_size, CFile& selectfp, CFile& rejectfp)
{
   size_t numlines = lines->listlength() ;
   if (sample_size >= numlines)
      {
      while (lines)
	 {
         Owned<StringList> line = lines ;
	 lines = lines->next() ;
	 write_line(selectfp,line->string()) ;
         }
      }
   else
      {
      auto selected = RandomSample(numlines,sample_size) ;
      for (unsigned i = 0 ; i < numlines ; i++)
	 {
         Owned<StringList> line = lines ;
	 lines = lines->next() ;
	 select_line(selected[i],selectfp,rejectfp,line->string()) ;
	 }
      }
   return ;
}

//----------------------------------------------------------------------

int main(int argc, char **argv)
{
   bool uniform_sample = false ;
   bool use_bytes = false ;
   bool use_interval = false ;
   bool use_length = false ;
   bool randomized = true ;
   unsigned min_length = 0 ;
   unsigned max_length = (unsigned)~0 ;
   const char* output_file = "-" ;
   const char* reject_file = 0 ;
   const char* argv0 = argv[0] ;
   while (argc > 1 && argv[1][0] == '-')
      {
      switch (argv[1][1])
	 {
	 case 'b':
	    use_bytes = true ;
	    use_interval = false ;
	    use_length = false ;
	    uniform_sample = false ;
	    break ;
	 case 'i':
	    use_interval = true ;
	    use_bytes = false ;
	    use_length = false ;
	    uniform_sample = false ;
	    break ;
	 case 'l':
	    min_length=atoi(argv[1]+2) ;
	    use_length = true ;
	    use_bytes = false ;
	    use_interval = false ;
	    uniform_sample = false ;
	    break ;
	 case 'L':
	    max_length=atoi(argv[1]+2) ;
	    use_length = true ;
	    use_bytes = false ;
	    use_interval = false ;
	    uniform_sample = false ;
	    break ;
	 case 'o':
	    output_file = argv[1]+2 ;
	    break ;
	 case 'r':
	    reject_file = argv[1]+2 ;
	    break ;
	 case 'R':
	    randomized = false ;
	    break ;
	 case 'u':
	    uniform_sample = true ;
	    use_bytes = false ;
	    use_interval = false ;
	    use_length = false ;
	    break ;
	 default:
	    cerr << "Unrecognized option " << argv[1] << endl << endl ;
	    usage(argv0) ;
	    break ;
	 }
      argc-- ;
      argv++ ;
      }
   size_t sample_size = 0 ;
   if (use_length)
      {
      // ignore the rest of the commandline
      }
   else if (argc < 2)
      {
      usage(argv0) ;
      }
   else
      {
      sample_size = atoi(argv[1]) ;
      }
   StringList* lines = nullptr ;
   StringList** lastline = &lines ;
   COutputFile rejectfp(reject_file) ;
   size_t numlines = 0 ;
   CFile f(stdin) ;
   COutputFile selectfp(output_file) ;
   while (CharPtr line = f.getCLine())
      {
      if (use_length)
	 {
	 unsigned len = strlen(line) ;
	 select_line(len >= min_length && len <= max_length,selectfp,rejectfp,line) ;
	 }
      else if (use_interval)
	 {
	 select_line(numlines % sample_size == 0,selectfp,rejectfp,line) ;
	 }
      else
	 {
	 lines->append(line,lastline) ;
         }
      numlines++ ;
      }
   *lastline = nullptr ;		// terminate list of lines
   if (use_bytes)
      {
      take_uniform_bytes(lines,sample_size,selectfp,rejectfp) ;
      }
   else if (uniform_sample)
      {
      take_uniform_sample(lines,sample_size,selectfp,rejectfp) ;
      }
   else
      {
      if (randomized)
	 Randomize() ;
      take_random_sample(lines,sample_size,selectfp,rejectfp) ;
      }
   return 0 ;
}
