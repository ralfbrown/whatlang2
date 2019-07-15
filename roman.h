/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LangIdent: long n-gram-based language identification		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     roman.h							*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-15						*/
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

/************************************************************************/
/************************************************************************/

class Romanizer
   {
   public:
      static unsigned utf8codepoint(const char* buf, wchar_t& codepoint) ;
      static bool romanizable(wchar_t codepoint) ;
      static int romanize(wchar_t codepoint, char* buffer) ;
      static unsigned romanize(wchar_t codepoint, wchar_t& romanized1, wchar_t& romanized2) ;
   } ;

/************************************************************************/
/************************************************************************/

// end of file roman.h //
