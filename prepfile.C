/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     prepfile.C	preprocessed file input			*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-13						*/
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
#include <errno.h>
#include "prepfile.h"
#include "framepac/file.h"
#include "framepac/list.h"
#include "framepac/string.h"
#include "framepac/texttransforms.h"

using namespace Fr ;

/************************************************************************/
/*	Globals for this module						*/
/************************************************************************/

uint64_t PreprocessedInputFile::s_sample_bytes = ~0U ;
bool PreprocessedInputFile::s_sample_uniformly = true ;
bool PreprocessedInputFile::s_convert_Latin1 = false ;
bool PreprocessedInputFile::s_ignore_whitespace = false ;
BigramExtension PreprocessedInputFile::s_bigram_ext = BigramExt_None ;
unsigned PreprocessedInputFile::s_alignment = 1 ;
Fr::CharPtr PreprocessedInputFile::s_from_enc = nullptr ;
Fr::CharPtr PreprocessedInputFile::s_to_enc = nullptr ;

/************************************************************************/
/*	Methods for class PreprocessedInputFile				*/
/************************************************************************/

PreprocessedInputFile::PreprocessedInputFile()
{
   m_buffered_lines = Fr::List::emptyList() ;
   m_max_sample_bytes = (size_t)~0 ;
   m_uniform_sample = true ;
#ifndef NO_ICONV
   m_conversion = (iconv_t)-1 ;
#endif /* NO_ICONV */
   return ;
}

//----------------------------------------------------------------------

PreprocessedInputFile::PreprocessedInputFile(const char *filename,
					     size_t sample_limit,
					     bool uniform_sample,
					     const char *from_enc, const char *to_enc)
{
   m_buffered_lines = Fr::List::emptyList() ;
#ifndef NO_ICONV
   m_conversion = (iconv_t)-1 ;
#endif /* NO_ICONV */
   open(filename,sample_limit,uniform_sample,from_enc,to_enc) ;
   return ;
}

//----------------------------------------------------------------------

bool PreprocessedInputFile::initializeTransliteration(const char *from, const char *to)
{
#ifndef NO_ICONV
   m_conversion = (from && to) ? iconv_open(to,from) : (iconv_t)-1 ;
   if (m_conversion == (iconv_t)-1)
      return !(from && to) ;
   // reset conversion state
   iconv(m_conversion,0,0,0,0) ;
#endif /* !NO_ICONV */
   return true ;
}

//----------------------------------------------------------------------

bool PreprocessedInputFile::shutdownTransliteration()
{
#ifndef NO_ICONV
   if (m_conversion != (iconv_t)-1)
      iconv_close(m_conversion) ;
   m_conversion = (iconv_t)-1 ;
#endif /* !NO_ICONV */
   return true ;
}

//----------------------------------------------------------------------

Fr::CFile* PreprocessedInputFile::open_sampled_input_file(const char *filename,
						     size_t max_bytes)
{
   // load in the entire input file
   auto fp = new Fr::CInputFile(filename) ;
   if (!fp || !*fp)
      return fp ;
   ListBuilder lb ;
   size_t numlines = 0 ;
   size_t total_bytes = 0 ;
   while (Fr::String* line = fp->getline())
      {
      total_bytes += line->c_len() ;
      numlines++ ;
      lb += line ;
      }
   Fr::List* lines = lb.move() ;
   // subsample the lines we read to be just a little more than the
   //   desired number of bytes
   double interval = max_bytes / (double)total_bytes ;
   if (interval > 0.5)			// adjustment for high sampling rates
      interval += (interval-0.5)/6.0 ;
   size_t sampled = 0 ;
   if (interval >= 0.98)
      {
      m_buffered_lines = lines->reverse() ;
      lines = Fr::List::emptyList() ;
      sampled = total_bytes ;
      }
   else
      {
      double avgline = total_bytes / (double)numlines ;
      double count = interval / 2.0 ;
      m_buffered_lines = Fr::List::emptyList() ;
      while (lines)
	 {
	 StringPtr line = reinterpret_cast<String*>(poplist(lines)) ;
	 if (!line) continue ; //FIXME
	 size_t len = line->c_len() ;
	 double increment = interval * len / avgline ;
	 if (((size_t)(count + increment) > (size_t)count) ||
	     interval >= 1.0)
	    {
	    pushlist(line.move(),m_buffered_lines) ;
	    sampled += len ;
	    }
	 count += increment ;
	 }
      }
//   if (verbose)
      {
      cout << "  Sampled " << sampled << " bytes from input file" ;
      if (sampled < m_max_sample_bytes && max_bytes < total_bytes)
	 cout << " (requested " << max_bytes << ", filesize=" << total_bytes
	      << ")" ;
      cout << endl ;
      }
   // if we subsampled but didn't get enough bytes, try again with a higher
   //   limit
   if (sampled < m_max_sample_bytes && total_bytes >= max_bytes)
      {
      fp->close()  ;
      if (m_buffered_lines) m_buffered_lines->free() ;
      m_buffered_lines = Fr::List::emptyList() ;
      max_bytes *= (max_bytes / (double)sampled * 1.01) ;
      return open_sampled_input_file(filename,max_bytes) ;
      }
   else
      {
      // reset the EOF condition
      fp->seek(0) ;
      clearerr(fp->fp()) ;
      }
   return fp ;
}

//----------------------------------------------------------------------

bool PreprocessedInputFile::open(const char *filename, size_t sample_limit, bool uniform_sample,
				 const char *from_enc, const char *to_enc)
{
   close() ;
   m_filename = Fr::dup_string(filename) ;
   m_max_sample_bytes = sample_limit ;
   m_uniform_sample = uniform_sample ;
   m_bigram_ext = s_bigram_ext ;
   m_convert_Latin1 = s_convert_Latin1 ;
   m_ignore_whitespace = s_ignore_whitespace ;
   m_alignment = s_alignment ;
   m_bytes_read = 0 ;
   if (from_enc && to_enc)
      initializeTransliteration(from_enc,to_enc) ;
   if (sample_limit != ~0U
       && m_alignment == 1) // can't currently sample UTF16
      {
      m_fp = open_sampled_input_file(filename,sample_limit) ;
      }
   else
      m_fp = new Fr::CInputFile(filename) ;
   return good() ;
}

//----------------------------------------------------------------------

void PreprocessedInputFile::close()
{
   if (m_fp)
      {
      // close the input file
      m_fp = nullptr ;
      }
   original_buffer_len = 0 ;
   // discard any remnants of transliteration
   translit_buffer_ptr = 0 ;
   translit_buffer_len = 0 ;
   // discard any remnants of buffered data from subsampling
   m_buffered_lines->free() ;
   m_buffered_lines = Fr::List::emptyList() ;
   buffered_char = '\0' ;
   // free iconv() resources
   shutdownTransliteration() ;
   return ;
}

//----------------------------------------------------------------------

int PreprocessedInputFile::readInput(unsigned char *buf, size_t buflen)
{
   if (m_buffered_lines && m_buffered_lines != Fr::List::emptyList())
      {
      int count = 0 ;
      while ((unsigned long)count < buflen && m_buffered_lines != Fr::List::emptyList())
	 {
	 size_t len = ((String*)m_buffered_lines->front())->c_len() ;
	 if (count + len > buflen)
	    break ;
	 StringPtr line = (String*)poplist(m_buffered_lines) ;
	 if (!line) continue ;
	 std::copy_n(line->c_str(),len,buf) ;
	 buf += len ;
	 count += len ;
	 }
      return count ;
      }
   else if (m_fp)
      return m_fp.read(buf,buflen) ;
   else
      return -1 ;
}

//----------------------------------------------------------------------

int PreprocessedInputFile::fillBuffer()
{
#ifndef NO_ICONV
   if (m_conversion != (iconv_t)-1)
      {
      // try to fill the input buffer
      int remainder = sizeof(original_buffer) - original_buffer_len ;
      int read_count = readInput(original_buffer+original_buffer_len,
				 remainder) ;
      original_buffer_len += read_count ;
      // convert as much as possible
      char *origbuf = (char*)original_buffer ;
      size_t orig_len = original_buffer_len ;
      char *translitbuf = (char*)translit_buffer ;
      size_t translit_len = sizeof(translit_buffer) ;
      errno = 0 ;
      size_t count = iconv(m_conversion,&origbuf,&orig_len,&translitbuf,
			   &translit_len) ;
      if (count != (size_t)-1)
	 errno = 0 ;
      if (errno == EILSEQ)
	 {
	 // bad input, try to skip ahead
	 if (translit_len > 0)
	    {
	    *translitbuf++ = *origbuf++ ;
	    orig_len-- ;
	    translit_len-- ;
	    }
	 }
      else if (errno != EINVAL && errno != E2BIG)
	 {
	 // only need to deal with the error if it is not an
	 //   incomplete multibyte sequence or running out of room in
	 //   the output buffer; in this case, we just copy as much as
	 //   we can to the output buffer
	 while (orig_len > 0 && translit_len > 0)
	    {
	    *translitbuf++ = *origbuf++ ;
	    orig_len-- ;
	    translit_len-- ;
	    }
	 // reset the conversion state
	 iconv(m_conversion,0,0,0,0) ;
	 }
      else if (read_count == 0) // at end of file?
	 {
	 // anything still left in original_buffer is unconvertible
	 //  bytes, so just copy them over
	 translit_buffer_len = original_buffer_len ;
	 if (original_buffer_len)
	    {
	    std::copy_n(original_buffer,original_buffer_len,translit_buffer) ;
	    original_buffer_len = 0 ;
	    }
	 translit_buffer_ptr = 0 ;
	 return translit_buffer_len ;
	 }
      if (orig_len > 0)
	 {
	 // we have data remaining in the input buffer, so
	 //  copy it to the beginning for the next refill
	 std::copy_n(origbuf,orig_len,original_buffer) ;
	 }
      original_buffer_len = orig_len ;
      translit_buffer_ptr = 0 ;
      translit_buffer_len = (translitbuf - (char*)translit_buffer) ;
      }
   else
#endif /* NO_ICONV */
      {
      translit_buffer_ptr = 0 ;
      translit_buffer_len
	 = readInput(translit_buffer,sizeof(translit_buffer)) ;
      }
   return translit_buffer_len ;
}

//----------------------------------------------------------------------

bool PreprocessedInputFile::moreData() const
{
   if (m_bytes_read >= m_max_sample_bytes)
      return false ;
   return (translit_buffer_len > translit_buffer_ptr) || !m_fp.eof() ;
}

//----------------------------------------------------------------------

int PreprocessedInputFile::peekAtBuffer()
{
   if (translit_buffer_ptr >= translit_buffer_len)
      {
      if (fillBuffer() <= 0)
	 return EOF ;
      }
   return translit_buffer[translit_buffer_ptr] ; // don't advance pointer
}

//----------------------------------------------------------------------

int PreprocessedInputFile::peekByte()
{
   if (m_bytes_read >= m_max_sample_bytes)
      return EOF ;
   if (m_convert_Latin1)
      {
      if (buffered_char >= 0x80)
	 return buffered_char ;
      if (ignoringWhitespace())
	 {
	 while (peekAtBuffer() == ' ')
	    {
	    (void)getFromBuffer() ;
	    }
	 }
      return peekAtBuffer() ;
      }
   else if (m_bigram_ext == BigramExt_None)
      {
      if (ignoringWhitespace())
	 {
	 while (peekAtBuffer() == ' ')
	    {
	    (void)getFromBuffer() ;
	    }
	 }
      return peekAtBuffer() ;
      }
   else if ((m_bytes_read & 1) != 0)
      return buffered_char ;
   else
      return peekAtBuffer() ;
}

//----------------------------------------------------------------------

int PreprocessedInputFile::getFromBuffer()
{
   if (translit_buffer_ptr >= translit_buffer_len)
      {
      if (fillBuffer() <= 0)
	 return EOF ;
      }
   return translit_buffer[translit_buffer_ptr++] ; // advance pointer
}

//----------------------------------------------------------------------

unsigned PreprocessedInputFile::getCodepoint()
{
   int byte = getFromBuffer() ;
   if (byte == EOF)
      {
      return byte ;
      }
   // ensure no sign-extension
   unsigned codepoint = byte & 0xFF ;
   if ((byte & 0x80) != 0 &&
       !(m_bigram_ext == BigramExt_ASCIILittleEndian ||
	 m_bigram_ext == BigramExt_ASCIIBigEndian))
      {
      // figure out how many bytes make up the code point, and put the
      //   highest-order bits stored in the first byte into the
      //   codepoint variable
      unsigned extra = 0 ;
      if ((byte & 0xE0) == 0xC0)
	 {
	 codepoint = (byte & 0x1F) ;
	 extra = 1 ;
	 }
      else if ((byte & 0xF0) == 0xE0)
	 {
	 codepoint = (byte & 0x0F) ;
	 extra = 2 ;
	 }
      else if ((byte & 0xF8) == 0xF0)
	 {
	 codepoint = (byte & 0x07) ;
	 extra = 3 ;
	 }
      else if ((byte & 0xFC) == 0xF8)
	 {
	 codepoint = (byte & 0x03) ;
	 extra = 4 ;
	 }
      else if ((byte & 0xFE) == 0xFC)
	 {
	 codepoint = (byte & 0x01) ;
	 extra = 5 ;
	 }
      // read in the additional bytes specified by the high bits of the
      //   codepoint's first byte
      for (size_t i = 1 ; i <= extra ; i++)
	 {
	 byte = getFromBuffer() ;
	 if (byte == EOF)
	    return byte ;
	 else if ((byte & 0xC0) != 0x80)
	    break ;			// invalid UTF8
	 // each extra byte gives us six more bits of the codepoint
	 codepoint = (codepoint << 6) | (byte & 0x3F) ;
	 }
      }
   return codepoint ;
}

//----------------------------------------------------------------------

int PreprocessedInputFile::getByte()
{
   if (m_bytes_read >= m_max_sample_bytes)
      {
      translit_buffer_ptr = translit_buffer_len ;
      return EOF ;
      }
   if (m_convert_Latin1)
      {
      int cp ;
      if (buffered_char >= 0x80)
	 {
	 cp = buffered_char ;
	 buffered_char = 0 ;
	 }
      else
	 {
	 do {
	    cp = getFromBuffer() ;
	    } while (ignoringWhitespace() && cp == ' ') ;
	 if (cp != EOF)
	    {
	    if ((cp & 0xFF) >= 0x80)
	       {
	       buffered_char = (unsigned char)(0x80 | (cp & 0x3F)) ;
	       cp = 0xC0 | ((cp & 0xFF) >> 6) ;
	       }
	    }
	 }
      if (cp != EOF)
	 ++m_bytes_read ;
      return cp ;
      }
   if (m_bigram_ext == BigramExt_None)
      {
      int cp ;
      do {
         cp = getFromBuffer() ;
         } while (ignoringWhitespace() && cp == ' ') ;
      ++m_bytes_read ;
      return cp ;
      }
   if ((m_bytes_read & 1) != 0)
      {
      ++m_bytes_read ;
      return buffered_char ;
      }
   int cp = getCodepoint() ;
   if (cp == EOF)
      return cp ;
   ++m_bytes_read ;
   if (m_bigram_ext == BigramExt_ASCIILittleEndian ||
       m_bigram_ext == BigramExt_UTF8LittleEndian)
      {
      buffered_char = (unsigned char)((cp >> 8) & 0xFF) ;
      return (cp & 0xFF) ;
      }
   else
      {
      buffered_char = (unsigned char)(cp & 0xFF) ;
      return (cp >> 8) & 0xFF ;
      }
}

//----------------------------------------------------------------------

bool PreprocessedInputFile::setDefaultTransliteration(const char *from, const char *to)
{
   if (!from && !to)
      {
      s_from_enc = s_to_enc = nullptr ;
      return true ;
      }
   if (!from || !to)
      return false ;
#ifdef NO_ICONV
   (void)from ; (void)to ;
   return false ;
#else
   iconv_t conversion = iconv_open(to,from) ;
   bool can_translit = (conversion != (iconv_t)-1) ;
   if (can_translit)
      {
      iconv_close(conversion) ;
      s_from_enc = Fr::dup_string(from) ;
      s_to_enc = Fr::dup_string(to) ;
      }
   return can_translit ;
#endif
}


/************************************************************************/
/************************************************************************/

// end of file prepfile.C //

