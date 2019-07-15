/****************************** -*- C++ -*- *****************************/
/*                                                                      */
/*	LangIdent: long n-gram-based language identification		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     smooth.C  inter-string score smoothing			*/
/*  Version:  1.30							*/
/*  LastEdit: 2019-07-15 						*/
/*                                                                      */
/*  (c) Copyright 2011,2012,2013,2019 Carnegie Mellon University	*/
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

#include <cmath>
#include "langid.h"

using namespace Fr ;

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

/************************************************************************/
/*	Methods for class LanguageIdentifier				*/
/************************************************************************/

Owned<LanguageScores> LanguageIdentifier::smoothedScores(LanguageScores* scores, int match_length) const
{
   if (!scores || !smoothingScores())
      return scores ;
   // use exponential decay as the smoothing function, since it is far
   //   simpler and runs faster than the other functions tried, yet
   //   works at least as well
   if (!m_prior_scores)
      {
      m_prior_scores = new LanguageScores(scores->numLanguages()) ;
      m_prior_scores->addThresholded(scores,LANGID_ZERO_SCORE, ::log(match_length)) ;
      return scores ;
      }
   m_prior_scores->scaleScores(SMOOTHING_DECAY_FACTOR) ;
   // adaptively weight the current sentence relative to the smoothing
   //   scores
   double max_score = scores->highestScore() ;
   double scaled = max_score / UNSURE_CUTOFF ;
   double score_weight = 0.5 * ::cbrt(match_length) + 0.3 * ::pow(scaled,1.33) ;
   if (score_weight < 0.0)
      score_weight = 0.0 ;
   double lambda = score_weight / (1.0 + score_weight) ;
   // give a little more smoothing weight to longer strings, since their
   //   scores are more reliable
   double smoothwt = 2.0 + 0.25 * ::log(match_length) ;
   scores->lambdaCombineWithPrior(m_prior_scores,lambda,smoothwt) ;
   return scores ;
}

// end of file smooth.C //
