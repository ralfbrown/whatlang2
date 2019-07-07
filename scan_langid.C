/************************************************************************/
/*									*/
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     scan_langid.C						*/
/*  Version:  1.00				       			*/
/*  LastEdit: 26aug2011							*/
/*									*/
/*  (c) Copyright 2011 Ralf Brown/CMU					*/
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

#ifdef BULK_EXTRACTOR

#include <cstdlib>
#include <cstdio>
#include "bulk_extractor.h"
#include "langid.h"
#include "FramepaC.h"

/************************************************************************/
/************************************************************************/

static unsigned sequence_number = 0 ;
static char *output_dir = "extract%" ;

/************************************************************************/
/************************************************************************/

static void process_buffer(const sbuf_t scanbuf, const char *output_directory)
{
   const uint8_t *buffer_start = scanbuf.buf ;
   const uint8_t *buffer_end = scanbuf.buf + scanbuf.size() ;
#if 0
   if (!process_file_data(buffer_start,buffer_end,"-",output_directory,false,
			  sequence_number,true))
#endif
      {

      }
   return ;
}

//----------------------------------------------------------------------

extern "C" void scan_langid(const class scanner_params &sp,
			    const class recursion_control_block &rcb)
{
   switch (sp.phase)
      {
      case 0:				// startup
	 //FIXME
	 break ;
      case 1:				// normal scan
	 process_buffer(sp.sbuf,output_dir) ;
	 break ;
      case 2:				// shutdown
	 //FIXME
	 break ;
      default:
	 fprintf(stderr,"Invalid 'phase' parameter to scan_ziprec\n") ;
	 break ;
      }
   return ;
}

#endif /* BULK_EXTRACTOR */

// end of file scan_strings.C //
