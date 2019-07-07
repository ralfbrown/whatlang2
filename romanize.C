/************************************************************************/
/*                                                                      */
/*	LangIdent: long n-gram-based language identification		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     romanize.C						*/
/*  Version:  1.12							*/
/*  LastEdit: 08jan2012							*/
/*                                                                      */
/*  (c) Copyright 2011,2012 Ralf Brown/Carnegie Mellon University	*/
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

#include <cstdlib>
#include <cstdio>
#include <memory.h>
#include "roman.h"

/************************************************************************/
/************************************************************************/

static void usage(const char *argv0)
{
   fprintf(stderr,"Romanize UTF-8 text in Greek, Cyrillic, or Arabic\n") ;
   fprintf(stderr,"Usage: %s <utf8 >romanized\n",argv0) ;
   exit(1) ;
}

//----------------------------------------------------------------------

int main(int argc, const char **argv)
{
   if (argc >= 2)
      {
      usage(argv[0]) ;
      return 1 ;
      }
   char buffer[32768] ;
   unsigned bufpos = 0 ;
   unsigned buflen = 0 ;
   const unsigned highwater = 7 * sizeof(buffer) / 8 ;
   buflen = fread(buffer,sizeof(char),sizeof(buffer),stdin) ;
   while (buflen > bufpos)
      {
      if (bufpos >= highwater)
	 {
	 // refill the buffer
	 memcpy(buffer,buffer+bufpos,buflen - bufpos) ;
	 buflen -= bufpos ;
	 bufpos = 0 ;
	 buflen += fread(buffer+buflen,1,sizeof(buffer)-buflen,stdin) ;
	 }
      // get the next UTF-8 codepoint
      wchar_t codepoint ;
      unsigned len = get_UTF8_codepoint(buffer + bufpos,codepoint) ;
      if (len == 0)
	 {
	 fprintf(stderr,"Invalid UTF-8 in input, aborting\n") ;
	 break ;
	 }
      bufpos += len ;
      // romanize it
      char romanized[8] ;
      int bytes = romanize_codepoint(codepoint,romanized) ;
      // and write the result
      if (bytes > 0)
	 {
	 fwrite(romanized,sizeof(char),bytes,stdout) ;
	 }
      else
	 {
	 fprintf(stderr,"Error romanizing codepoint U%4.04X\n",codepoint) ;
	 break ;
	 }
      }
   return 0 ;
}

// end of file romanize.C //
