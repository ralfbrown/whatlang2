/************************************************************************/
/*                                                                      */
/*	LangIdent: long n-gram-based language identification		*/
/*	by Ralf Brown / Carnegie Mellon University			*/
/*									*/
/*  File:     smooth.C  inter-string score smoothing			*/
/*  Version:  1.23							*/
/*  LastEdit: 20aug2013 						*/
/*                                                                      */
/*  (c) Copyright 2011,2012,2013 Ralf Brown/Carnegie Mellon University	*/
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

/************************************************************************/
/************************************************************************/

// how much above the minimal score must a language score be to be considered
//   a reliable identification (and not get a question mark)?
#define UNSURE_CUTOFF (120 * LANGID_ZERO_SCORE)

// set the multiplicative factor by which to decay the prior scores for each
//   new string
#define SMOOTHING_DECAY_FACTOR 0.25

/************************************************************************/
/*	Global variables						*/
/************************************************************************/

static bool do_smoothing = true ;
static LanguageScores *prior_language_scores = 0 ;

/************************************************************************/
/************************************************************************/

bool smooth_language_scores(bool smooth)
{
   bool old_smooth = do_smoothing ;
   do_smoothing = smooth ;
   return old_smooth ;
}

//----------------------------------------------------------------------

bool smoothing_language_scores()
{
   return do_smoothing ;
}

//----------------------------------------------------------------------

LanguageScores *smoothed_language_scores(LanguageScores *scores,
					 LanguageScores *&prior_scores,
					 size_t match_length)
{
   if (!scores || !smoothing_language_scores())
      return scores ;
   // use exponential decay as the smoothing function, since it is far
   //   simpler and runs faster than the other functions tried, yet
   //   works at least as well
   if (!prior_scores)
      {
      prior_scores = new LanguageScores(scores->numLanguages()) ;
      prior_scores->addThresholded(scores,LANGID_ZERO_SCORE,
				   ::log(match_length)) ;
      return scores ;
      }
   prior_scores->scaleScores(SMOOTHING_DECAY_FACTOR) ;
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
   scores->lambdaCombineWithPrior(prior_scores,lambda,smoothwt) ;
   return scores ;
}

//----------------------------------------------------------------------

LanguageScores *smoothed_language_scores(LanguageScores *scores,
					 size_t match_length)
{
   return smoothed_language_scores(scores,prior_language_scores,match_length) ;
}

// end of file smooth.C //
