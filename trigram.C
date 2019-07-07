/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LA-Strings: language-aware text-strings extraction		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     mklangid.C						*/
/*  Version:  1.14							*/
/*  LastEdit: 14mar2012							*/
/*                                                                      */
/*  (c) Copyright 2010,2011,2012 Ralf Brown/Carnegie Mellon University	*/
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
#include "trie.h"

using namespace std ;

/************************************************************************/
/************************************************************************/

uint32_t adjusted_threshold(const uint32_t *frequencies)
{
   // if the threshold is zero, we have fewer than topK ngrams, so no
   //  pruning is necessary -- set the threshold to one; if greater
   //  than one, we probably won't go too far over topK upon pruning.
   //  But if the given threshold is one, the effect is no pruning,
   //  which definitely leaves us with too many n-grams, so increase
   //  it in that case
   uint32_t threshold = frequencies[0] ;
   return threshold < 2 ? 2 : threshold ;
}

//----------------------------------------------------------------------

void insert_frequency(uint32_t newelt, uint32_t *heap, size_t heaplen)
{
   // do a partial heap-sort -- we only need to know the lowest element
   //   in the fixed-size heap, so we can skip the final sort and merely
   //   maintain the heap as we insert values

   // discard the lowest-scoring element, place the new one in its place,
   //   and then bubble it down to the proper place in the heap
   heap[0] = newelt ;
   size_t N = 0 ;
   while (N < heaplen/2)
      {
      // left child = 2 * N + 1, right child = 2 * N + 2
      size_t child = (N << 1) + 1 ;
      // if either child is smaller than current node, swap with the
      //   smaller of the children
      if (child + 1 < heaplen && heap[child] > heap[child+1])
	 child++ ;
      // but if the node is smaller than both children, we have a valid heap
      //   and can stop
      if (heap[N] <= heap[child])
	 return ;
      uint32_t tmp = heap[child] ;
      heap[child] = heap[N] ;
      heap[N] = tmp ;
      N = child ;
      }
   return ;
}

/************************************************************************/
/*	Methods for class BigramCounts					*/
/************************************************************************/

BigramCounts::BigramCounts(const TrigramCounts &trigrams)
{
   m_total = 0 ;
   for (size_t c1 = 0 ; c1 <= 0xFF ; c1++)
      {
      for (size_t c2 = 0 ; c2 <= 0xFF ; c2++)
	 {
	 uint32_t cnt = trigrams.totalCount(c1,c2) ;
	 m_total += cnt ;
	 set(c1,c2,cnt) ;
	 }
      }
   return ;
}

/************************************************************************/
/*	Methods for class TrigramCounts					*/
/************************************************************************/

TrigramCounts::TrigramCounts(const TrigramCounts *orig)
{
   copy(orig) ;
   return ;
}

//----------------------------------------------------------------------

void TrigramCounts::copy(const TrigramCounts *orig)
{
   if (orig)
      {
      for (size_t i = 0 ; i < lengthof(m_counts) ; i++)
	 m_counts[i] = orig->m_counts[i] ;
      }
   else
      {
      for (size_t i = 0 ; i < lengthof(m_counts) ; i++)
	 m_counts[i] = 0 ;
      }
   return ;
}

//----------------------------------------------------------------------

uint32_t TrigramCounts::totalCount(uint8_t c1, uint8_t c2) const
{
   const uint32_t *values = &m_counts[(c1 << 16) + (c2 << 8)] ;
   uint32_t total = 0 ;
   for (size_t c3 = 0 ; c3 <= 0xFF ; c3++)
      {
      total += *(values++) ;
      }
   return total ;
}

//----------------------------------------------------------------------

bool TrigramCounts::enumerate(NybbleTrie &ngrams) const
{
   uint8_t c[3] ;
   for (unsigned c1 = 0 ; c1 < 256 ; c1++)
      {
      c[0] = c1 ;
      for (unsigned c2 = 0 ; c2 < 256 ; c2++)
	 {
	 c[1] = c2 ;
	 for (unsigned c3 = 0 ; c3 < 256 ; c3++)
	    {
	    uint32_t freq = count(c1,c2,c3) ;
	    if (freq > 0)
	       {
	       c[2] = c3 ;
	       ngrams.insert(c,3,freq,false) ;
	       }
	    }
	 }
      }
   return true ;
}

//----------------------------------------------------------------------

void TrigramCounts::filter(unsigned topK, unsigned max_len, bool verbose)
{
   if (topK == 0)
      return ;
   if (verbose)
      {
      cout << "Determining trigram cut-off-frequency" << endl ;
      }
   uint32_t top_frequencies[topK] ;
   memset(top_frequencies,'\0',topK*sizeof(uint32_t)) ;
   uint32_t min_freq = 1 ;
   // force skipping 00/00/00 and FF/FF/FF, since they are common fillers
   //   in binary files and will thus tend to clog the top-K list with
   //   long n-grams consisting of nothing but 00 or FF bytes.
   m_counts[0] = 0 ;
   m_counts[lengthof(m_counts)-1] = 0 ;
   for (size_t i = 0 ; i < lengthof(m_counts) ; i++)
      {
      if (m_counts[i] > min_freq)
	 {
	 insert_frequency(m_counts[i],top_frequencies,topK) ;
	 if (top_frequencies[0] > min_freq)
	    min_freq = top_frequencies[0] ;
	 }
      }
   // after processing the entire count table, the smallest value in
   //  'top_frequencies' is our cut-off; set all counts less than that
   //  value to zero
   uint32_t thresh = adjusted_threshold(top_frequencies) ;
   if (verbose)
      {
      cout << "Trigram cut-off frequency @ " << topK << " = " << thresh
	   << endl ;
      }
   unsigned distinct = 0 ;
   for (size_t i = 0 ; i < lengthof(m_counts) ; i++)
      {
      if (m_counts[i] < thresh)
	 m_counts[i] = 0 ;
      else
	 distinct++ ;
      }
   if (max_len < 3) max_len = 3 ;
   if (max_len > 100) max_len = 100 ;
   unsigned required = topK / (unsigned)::pow(1.5,max_len-3) ;
   if (distinct < required)
      {
      cout << "Fewer than " << required << " distinct trigrams -- you may need more training data" << endl ;
      }
   return ;
}

//----------------------------------------------------------------------

TrigramCounts *TrigramCounts::load(FILE *fp)
{
   if (fp)
      {
      TrigramCounts *model = new TrigramCounts ;
      if (model->read(fp))
	 return model ;
      delete model ;
      }
   return 0 ;
}

//----------------------------------------------------------------------

bool TrigramCounts::read(FILE *fp)
{
   if (fp)
      {
      return (fread(m_counts,sizeof(m_counts[0]),lengthof(m_counts),fp)
	      == lengthof(m_counts)) ;
      }
   return false ;
}

//----------------------------------------------------------------------

bool TrigramCounts::save(FILE *fp) const
{
   if (fp)
      {
      return (fwrite(m_counts,sizeof(m_counts[0]),lengthof(m_counts),fp)
	      == lengthof(m_counts)) ;
      }
   return false ;
}

// end of file trigram.C //
