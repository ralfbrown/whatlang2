/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     prepfile.h	preprocessed file input			*/
/*  Version:  1.30							*/
/*  LastEdit: 27jun2019							*/
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

//#include "FramepaC.h"
#include "framepac/file.h"
#ifndef NO_ICONV
# include <iconv.h>
#endif

using namespace std ;

/************************************************************************/
/*	Manifest Constants						*/
/************************************************************************/

#define BUFFER_SIZE 65536U


/************************************************************************/
/*	Types for this module						*/
/************************************************************************/

enum BigramExtension
   {
      BigramExt_None,
      BigramExt_ASCIILittleEndian,
      BigramExt_ASCIIBigEndian,
      BigramExt_UTF8LittleEndian,
      BigramExt_UTF8BigEndian
   } ;

//----------------------------------------------------------------------

class PreprocessedInputFile
   {
   private:
      Fr::CharPtr m_filename ;
      Fr::CFile*  m_fp ;
      size_t	  m_max_sample_bytes ;
      uint64_t	  m_bytes_read ;
      bool	  m_piped ;
      bool	  m_uniform_sample ;
      bool	  m_convert_Latin1 ;
      bool        m_ignore_whitespace ;
      BigramExtension m_bigram_ext ;
      unsigned    m_alignment ;
      Fr::List   *m_buffered_lines ;
      unsigned char buffered_char ;
      unsigned    original_buffer_len ;
      unsigned    translit_buffer_ptr ;
      unsigned    translit_buffer_len ;
#ifndef NO_ICONV
      iconv_t     m_conversion ;
      unsigned char original_buffer[BUFFER_SIZE] ;
#endif /* !NO_ICONV */
      unsigned char translit_buffer[2*BUFFER_SIZE] ;
      static Fr::CharPtr s_from_enc ;
      static Fr::CharPtr s_to_enc ;
      static uint64_t s_sample_bytes ;
      static bool   s_sample_uniformly ;
      static bool   s_convert_Latin1 ;
      static bool   s_ignore_whitespace ;
      static BigramExtension s_bigram_ext ;
      static unsigned s_alignment ;
   public://FIXME

   protected:
      Fr::CFile* open_sampled_input_file(const char *filename, size_t max_bytes) ;
      bool initializeTransliteration(const char *from, const char *to) ;
      bool shutdownTransliteration() ;
      int readInput(unsigned char *buf, size_t buflen) ;
      int fillBuffer() ;
      int peekAtBuffer() ;
      int getFromBuffer() ;
      unsigned getCodepoint() ;

   public:
      PreprocessedInputFile() ;
      PreprocessedInputFile(const char *filename, uint64_t sample_limit = s_sample_bytes,
			    bool uniform_sample = s_sample_uniformly,
			    const char *from_enc = s_from_enc, const char *to_enc = s_to_enc) ;
      ~PreprocessedInputFile() ;

      // accesse to state
      bool good() const { return m_fp != nullptr ; }
      bool ignoringWhitespace() const { return m_ignore_whitespace ; }
      BigramExtension bigramExt() const { return m_bigram_ext ; }
      uint64_t bytesRead() const { return m_bytes_read ; }

      bool open(const char *filename, uint64_t sample_limit = s_sample_bytes,
		bool uniform_sample = s_sample_uniformly,
		const char *from_enc = nullptr, const char *to_enc = nullptr) ;
      void close() ;

      // input
      bool moreData() const ;
      int peekByte() ;
      int getByte() ;

      // configuration
      static void setSampling(uint64_t sample_limit, bool uniform_sample = true)
         { s_sample_bytes = sample_limit ; s_sample_uniformly = uniform_sample ; }
      void setBigramExt(BigramExtension ext) { m_bigram_ext = ext ; }
      static void setDefaultBigramExt(BigramExtension ext) { s_bigram_ext = ext ; }
      void setConvertLatin1(bool cnv) { m_convert_Latin1 = cnv ; }
      static void setDefaultConvertLatin1(bool cnv) { s_convert_Latin1 = cnv ; }
      void setAlignment(unsigned align) { m_alignment = align ; }
      static void setDefaultAlignment(unsigned align) { s_alignment = align ; }
      void ignoreWhitespace(bool ignore) { m_ignore_whitespace = ignore ; }
      static void setIgnoreWhitespace(bool ignore) { s_ignore_whitespace = ignore ; }
      static bool setDefaultTransliteration(const char *from, const char *to) ;
   } ;

// end of file prepfile.h //
