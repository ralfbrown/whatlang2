/************************************************************************/
/*                                                                      */
/*	LangIdent: long n-gram-based language identification		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     subsample.C  utility program to sample lines of text	*/
/*  Version:  1.24							*/
/*  LastEdit: 04aug2014 						*/
/*                                                                      */
/*  (c) Copyright 2012,2013,2014 Ralf Brown/Carnegie Mellon University	*/
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

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <memory.h>
#include <stdlib.h>
//#include <time.h>
using namespace std ;
#define FrHAVE_SRAND48

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define MAX_LINE 32768U

/************************************************************************/
/*	Types								*/
/************************************************************************/

class StringList
   {
   private:
      StringList *m_next ;
      char       *m_string ;
      unsigned	  m_length ;
   public:
      StringList(const char *strng) ;
      ~StringList() ;

      // accessors
      StringList *next() const { return m_next ; }
      const char *string() const { return m_string ; }
      unsigned length() const { return m_length ; }
      unsigned listlength() const ;

      // modifiers
      void next(StringList *nxt) { m_next = nxt ; }

      // factory
      static void append(const char *strng, StringList **&last_string) ;
      static void append(StringList *strng, StringList **&last_string) ;
   } ;

/************************************************************************/
/*	Methods for class StringList					*/
/************************************************************************/

StringList::StringList(const char *strng)
{
   m_next = nullptr ;
   if (strng)
      {
      m_length = strlen(strng) ;
      m_string = (char*)malloc(m_length+1) ;
      if (m_string)
	 memcpy(m_string,strng,m_length+1) ;
      else
	 m_length = 0 ;
      }
   else
      {
      m_length = 0 ;
      m_string = nullptr ;
      }
   return ;
}

//----------------------------------------------------------------------

StringList::~StringList()
{
   free(m_string) ;
   m_string = nullptr ;
   m_length = 0 ;
   m_next = nullptr ;
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

void StringList::append(const char *strng, StringList **&last_string)
{
   if (!last_string)
      return ;				// no list, so we can't do anything
   StringList *str = new StringList(strng) ;
   if (!str)
      return ;				// Out of memory!  Ignore line
   *last_string = str ;
   StringList **end_ptr = &str->m_next ;
   last_string = end_ptr ;
   return ;
}

//----------------------------------------------------------------------

void StringList::append(StringList *strng, StringList **&last_string)
{
   if (!last_string)
      return ;				// no list, so we can't do anything
   if (!strng)
      return ;				// Out of memory!  Ignore line
   strng->next(0) ;
   *last_string = strng ;
   StringList **end_ptr = &strng->m_next ;
   last_string = end_ptr ;
   return ;
}

/************************************************************************/
/*	Random sampling							*/
/************************************************************************/

size_t random_number(size_t range)
{
#ifdef FrHAVE_SRAND48
   long rn = lrand48() ;
#else
   if (range > RAND_MAX)
      {
      static bool warned = false ;
      if (!warned)
	 {
	 fprintf(stderr,"Warning: random number generator does not have a\n"
		 "\tsufficiently large range.") ;
	 warned = true ;
	 }
      }
   long rn = (long)rand() ;
#endif /* FrHAVE_SRAND48 */
   return rn % range ;
}

//----------------------------------------------------------------------

char* random_sample(size_t total_size, size_t sample_size)
{
   char *selected = (char*)calloc(total_size+1,sizeof(char)) ;
   if (!selected)
      {
      fprintf(stderr,"Out of memory while generating random sampling") ;
      return nullptr ;
      }
   // seed random number gen from time
#ifdef FrHAVE_SRAND48
   srand48(time(0)) ;
#else
   srand(time(0)) ;
#endif /* FrHAVE_SRAND48 */
   if (sample_size > total_size / 2)
      {
      // to avoid looping nearly forever on big samples, turn it into the
      //   equivalent problem of *deselecting* a small portion of the complete
      //   set of documents
      memset(selected,1,total_size) ;
      sample_size = total_size - sample_size ;
      for ( ; sample_size > 0 ; sample_size--)
	 {
	 size_t select ;
	 do {
	    select = random_number(total_size) ;
	    } while (!selected[select]) ;
	 selected[select] = (char)0 ;
	 }
      }
   else
      {
      for ( ; sample_size > 0 ; sample_size--)
	 {
	 size_t select ;
	 do {
	    select = random_number(total_size) ;
	    } while (selected[select]) ;
	 selected[select] = (char)1 ;
	 }
      }
   return selected ;
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
           "\t-rF\twrite non-sampled (rejected) lines to file F\n"
	<< endl ;
   exit(1) ;
}

//----------------------------------------------------------------------

static void take_uniform_bytes(StringList *lines, size_t sample_size,
			       FILE *rejectfp)
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
      StringList *line = lines ;
      lines = lines->next() ;
      size_t len = line->length() ;
      if (sampled_bytes <= (total_bytes * sample_rate))
	 {
	 fputs(line->string(),stdout) ;
	 sampled_bytes += len ;
	 }
      else if (rejectfp)
	 {
	 fputs(line->string(),rejectfp) ;
	 }
      total_bytes += len ;
      delete line ;
      }
   return ;
}

//----------------------------------------------------------------------

static void take_uniform_sample(StringList *lines, size_t sample_size,
				FILE *rejectfp)
{
   if (!lines)
      return ;
   size_t numlines = lines->listlength() ;
   double interval = sample_size / (double)numlines ;
   double count = interval/2.0 ;
   while (lines)
      {
      StringList *line = lines ;
      lines = lines->next() ;
      double increment = interval ;
      if ((size_t)(count + increment) > (size_t)count)
	 {
	 fputs(line->string(),stdout) ;
	 }
      else if (rejectfp)
	 {
	 fputs(line->string(),rejectfp) ;
	 }
      count += increment ;
      delete line ;
      }
   return ;
}

//----------------------------------------------------------------------

static void take_random_sample(StringList *lines, size_t sample_size, FILE *rejectfp)
{
   size_t numlines = lines->listlength() ;
   if (sample_size >= numlines)
      {
      while (lines)
	 {
         StringList *line = lines ;
	 lines = lines->next() ;
	 fputs(line->string(),stdout) ;
	 delete line ;
         }
      }
   else
      {
      char *selected = random_sample(numlines,sample_size) ;
      for (unsigned i = 0 ; i < numlines ; i++)
	 {
         StringList *line = lines ;
	 lines = lines->next() ;
	 if (selected[i])
	    fputs(line->string(),stdout) ;
	 else if (rejectfp)
	    fputs(line->string(),rejectfp) ;
	 delete line ;
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
   unsigned min_length = 0 ;
   unsigned max_length = (unsigned)~0 ;
   const char *reject_file = 0 ;
   const char *argv0 = argv[0] ;
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
	 case 'r':
	    reject_file = argv[1]+2 ;
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
   StringList *lines = nullptr ;
   StringList **lastline = &lines ;
   FILE *rejectfp = nullptr ;
   if (reject_file && *reject_file)
      {
      rejectfp = fopen(reject_file,"w") ;
      }
   size_t numlines = 0 ;
   while (!feof(stdin))
      {
      char line[MAX_LINE] ;
      line[0] = '\0' ;
      if (!fgets(line,sizeof(line),stdin))
	 break ;
      if (use_length)
	 {
	 unsigned len = strlen(line) ;
	 if (len >= min_length && len <= max_length)
	    fputs(line,stdout) ;
	 else if (rejectfp)
	    fputs(line,rejectfp) ;
	 }
      else if (use_interval)
	 {
	 if (numlines % sample_size == 0)
	    fputs(line,stdout) ;
	 else if (rejectfp)
	    fputs(line,rejectfp) ;
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
      take_uniform_bytes(lines,sample_size,rejectfp) ;
      }
   else if (uniform_sample)
      {
      take_uniform_sample(lines,sample_size,rejectfp) ;
      }
   else
      {
      take_random_sample(lines,sample_size,rejectfp) ;
      }
   if (rejectfp)
      fclose(rejectfp) ;
   return 0 ;
}
