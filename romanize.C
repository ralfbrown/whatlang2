/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LangIdent: long n-gram-based language identification		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     romanize.C						*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-30						*/
/*                                                                      */
/*  (c) Copyright 2011,2012,2019 Ralf Brown/Carnegie Mellon University	*/
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
#include <cstdio>
#include <cstring>
#include "framepac/file.h"
#include "framepac/romanize.h"

using namespace Fr ;

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
   Initialize() ;
   CFile f(stdin) ;
   while (CharPtr line = f.getCLine())
      {
      CharPtr romanized = Romanizer::romanize(line) ;
      if (romanized)
	 fwrite(romanized,sizeof(char),strlen(romanized),stdout) ;
      fputc('\n',stdout) ;
      }
   Shutdown() ;
   return 0 ;
}

// end of file romanize.C //
