/*
This file is part of CASToR.

    CASToR is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    CASToR is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    CASToR (in file GNU_GPL.TXT). If not, see <http://www.gnu.org/licenses/>.

Copyright 2017-2018 all CASToR contributors listed below:

    --> current contributors: Thibaut MERLIN, Simon STUTE, Didier BENOIT, Claude COMTAT, Marina FILIPOVIC, Mael MILLARDET
    --> past contributors: Valentin VIELZEUF

This is CASToR version 2.0.1.
*/

/*!
  \file
  \ingroup management

  \brief Implementation of class sChronoManager
*/

#include "sChronoManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

sChronoManager* sChronoManager::mp_Instance = NULL;

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

sChronoManager::sChronoManager()
{
  // Simply default all members
  mp_Instance = NULL;
  m_nbThreadsForProjection = 0;
  m_nbThreadsForImageComputation = 0;
  m_nbThreadsMax = 0;
  // Permanent profiled steps
  mp_durationIterativeDataUpdateStep1 = NULL;
  mp_durationIterativeDataUpdateStep2 = NULL;
  mp_durationIterativeDataUpdateStep3 = NULL;
  mp_durationIterativeDataUpdateStep4 = NULL;
  mp_startIterativeDataUpdateStep1 = NULL;
  mp_startIterativeDataUpdateStep2 = NULL;
  mp_startIterativeDataUpdateStep3 = NULL;
  mp_startIterativeDataUpdateStep4 = NULL;
  // Custom profiled steps
  m_nbCustomSteps = 0;
  mpp_durationCustomSteps = NULL;
  mpp_startCustomSteps = NULL;
  // Verbose
  m_verbose = 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int sChronoManager::CheckParameters()
{
  // Check number of threads
  if (m_nbThreadsForProjection<=0)
  {
    Cerr("***** sChronoManager::CheckParameters() -> Number of threads for projections is incorrectly set !" << endl);
    return 1;
  }
  if (m_nbThreadsForImageComputation<=0)
  {
    Cerr("***** sChronoManager::CheckParameters() -> Number of threads for image computation is incorrectly set !" << endl);
    return 1;
  }
  // Check number of custom steps
  if (m_nbCustomSteps<0)
  {
    Cerr("***** sChronoManager::CheckParameters() -> Number of custom steps is incorrectly set !" << endl);
    return 1;
  }
  // Verbose
  if (m_verbose<0)
  {
    Cerr("***** sChronoManager::CheckParameters() -> Verbose cannot be negative !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int sChronoManager::Initialize()
{
  // Verbose
  if (m_verbose>=VERBOSE_LIGHT) Cout("sChronoManager::Initialize() -> Initialize all duration counters for profiling" << endl);
  // Initialize buffers for data update steps
  mp_durationIterativeDataUpdateStep1 = (DurationNano*)malloc(m_nbThreadsForProjection*sizeof(DurationNano));
  mp_durationIterativeDataUpdateStep2 = (DurationNano*)malloc(m_nbThreadsForProjection*sizeof(DurationNano));
  mp_durationIterativeDataUpdateStep3 = (DurationNano*)malloc(m_nbThreadsForProjection*sizeof(DurationNano));
  mp_durationIterativeDataUpdateStep4 = (DurationNano*)malloc(m_nbThreadsForProjection*sizeof(DurationNano));
  mp_startIterativeDataUpdateStep1 = (ChronoTime*)malloc(m_nbThreadsForProjection*sizeof(ChronoTime));
  mp_startIterativeDataUpdateStep2 = (ChronoTime*)malloc(m_nbThreadsForProjection*sizeof(ChronoTime));
  mp_startIterativeDataUpdateStep3 = (ChronoTime*)malloc(m_nbThreadsForProjection*sizeof(ChronoTime));
  mp_startIterativeDataUpdateStep4 = (ChronoTime*)malloc(m_nbThreadsForProjection*sizeof(ChronoTime));
  for (int th=0; th<m_nbThreadsForProjection; th++)
  {
    mp_durationIterativeDataUpdateStep1[th] = chrono::duration<int64_t,nano>::zero();
    mp_durationIterativeDataUpdateStep2[th] = chrono::duration<int64_t,nano>::zero();
    mp_durationIterativeDataUpdateStep3[th] = chrono::duration<int64_t,nano>::zero();
    mp_durationIterativeDataUpdateStep4[th] = chrono::duration<int64_t,nano>::zero();
  }
  // Initialize convolution
  m_durationConvolution = chrono::duration<int64_t,nano>::zero();
  // Initialize buffers for custom profiling
  if (m_nbCustomSteps>0)
  {
    mpp_durationCustomSteps = (DurationNano**)malloc(m_nbCustomSteps*sizeof(DurationNano*));
    mpp_startCustomSteps = (ChronoTime**)malloc(m_nbCustomSteps*sizeof(ChronoTime*));
    for (int nb=0; nb<m_nbCustomSteps; nb++)
    {
      mpp_durationCustomSteps[nb] = (DurationNano*)malloc(m_nbThreadsMax*sizeof(DurationNano));
      mpp_startCustomSteps[nb] = (ChronoTime*)malloc(m_nbThreadsMax*sizeof(ChronoTime));
      for (int th=0; th<m_nbThreadsMax; th++) mpp_durationCustomSteps[nb][th] = chrono::duration<int64_t,nano>::zero();
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sChronoManager::Display()
{
  if (m_verbose>=VERBOSE_LIGHT)
  {
    // Verbose
    Cout("sChronoManager::Display() -> Results from the profiling" << endl);
    // Sum over all threads
    for (int th=1; th<m_nbThreadsForProjection; th++)
    {
      mp_durationIterativeDataUpdateStep1[0] += mp_durationIterativeDataUpdateStep1[th];
      mp_durationIterativeDataUpdateStep2[0] += mp_durationIterativeDataUpdateStep2[th];
      mp_durationIterativeDataUpdateStep3[0] += mp_durationIterativeDataUpdateStep3[th];
      mp_durationIterativeDataUpdateStep4[0] += mp_durationIterativeDataUpdateStep4[th];
    }
    for (int nb=0; nb<m_nbCustomSteps; nb++) for (int th=1; th<m_nbThreadsMax; th++)
      mpp_durationCustomSteps[nb][0] += mpp_durationCustomSteps[nb][th];
    // Divide by number of threads (not for the custom steps as we cannot know if they are multi-threaded
    // using the number of threads for projection or image computation)
    mp_durationIterativeDataUpdateStep1[0] /= m_nbThreadsForProjection;
    mp_durationIterativeDataUpdateStep2[0] /= m_nbThreadsForProjection;
    mp_durationIterativeDataUpdateStep3[0] /= m_nbThreadsForProjection;
    mp_durationIterativeDataUpdateStep4[0] /= m_nbThreadsForProjection;
    // Display data update steps
    Cout("  --> Profiling of the data update step:" << endl);
    Cout("    | Datafile management: " << setfill( '0' )
         << setw( 2 ) << chrono::duration_cast<Hs>  ( mp_durationIterativeDataUpdateStep1[0] ).count() << " hours "
         << setw( 2 ) << chrono::duration_cast<Mins>( mp_durationIterativeDataUpdateStep1[0] % Hs( 1 ) ).count() << " mins "
         << setw( 2 ) << chrono::duration_cast<Secs>( mp_durationIterativeDataUpdateStep1[0] % Mins( 1 ) ).count() << " secs "
         << setw( 3 ) << chrono::duration_cast<Ms>  ( mp_durationIterativeDataUpdateStep1[0] % Secs( 1 ) ).count() << " ms" << endl);
    Cout("    | Dynamic management: " << setfill( '0' )
         << setw( 2 ) << chrono::duration_cast<Hs>  ( mp_durationIterativeDataUpdateStep2[0] ).count() << " hours "
         << setw( 2 ) << chrono::duration_cast<Mins>( mp_durationIterativeDataUpdateStep2[0] % Hs( 1 ) ).count() << " mins "
         << setw( 2 ) << chrono::duration_cast<Secs>( mp_durationIterativeDataUpdateStep2[0] % Mins( 1 ) ).count() << " secs "
         << setw( 3 ) << chrono::duration_cast<Ms>  ( mp_durationIterativeDataUpdateStep2[0] % Secs( 1 ) ).count() << " ms" << endl);
    Cout("    | Projection management: " << setfill( '0' )
         << setw( 2 ) << chrono::duration_cast<Hs>  ( mp_durationIterativeDataUpdateStep3[0] ).count() << " hours "
         << setw( 2 ) << chrono::duration_cast<Mins>( mp_durationIterativeDataUpdateStep3[0] % Hs( 1 ) ).count() << " mins "
         << setw( 2 ) << chrono::duration_cast<Secs>( mp_durationIterativeDataUpdateStep3[0] % Mins( 1 ) ).count() << " secs "
         << setw( 3 ) << chrono::duration_cast<Ms>  ( mp_durationIterativeDataUpdateStep3[0] % Secs( 1 ) ).count() << " ms" << endl);
    Cout("    | Optimizer management: " << setfill( '0' )
         << setw( 2 ) << chrono::duration_cast<Hs>  ( mp_durationIterativeDataUpdateStep4[0] ).count() << " hours "
         << setw( 2 ) << chrono::duration_cast<Mins>( mp_durationIterativeDataUpdateStep4[0] % Hs( 1 ) ).count() << " mins "
         << setw( 2 ) << chrono::duration_cast<Secs>( mp_durationIterativeDataUpdateStep4[0] % Mins( 1 ) ).count() << " secs "
         << setw( 3 ) << chrono::duration_cast<Ms>  ( mp_durationIterativeDataUpdateStep4[0] % Secs( 1 ) ).count() << " ms" << endl);
    // Display convolution performances
    Cout("  --> Profiling of the convolution steps: " << setfill( '0' )
         << setw( 2 ) << chrono::duration_cast<Hs>  ( m_durationConvolution ).count() << " hours "
         << setw( 2 ) << chrono::duration_cast<Mins>( m_durationConvolution % Hs( 1 ) ).count() << " mins "
         << setw( 2 ) << chrono::duration_cast<Secs>( m_durationConvolution % Mins( 1 ) ).count() << " secs "
         << setw( 3 ) << chrono::duration_cast<Ms>  ( m_durationConvolution % Secs( 1 ) ).count() << " ms" << endl);
    // Display custom update steps
    for (int nb=0; nb<m_nbCustomSteps; nb++)
    {
      Cout("  --> Custom update step " << nb+1 << ": " << setfill( '0' )
         << setw( 2 ) << chrono::duration_cast<Hs>  ( mpp_durationCustomSteps[nb][0] ).count() << " hours "
         << setw( 2 ) << chrono::duration_cast<Mins>( mpp_durationCustomSteps[nb][0] % Hs( 1 ) ).count() << " mins "
         << setw( 2 ) << chrono::duration_cast<Secs>( mpp_durationCustomSteps[nb][0] % Mins( 1 ) ).count() << " secs "
         << setw( 3 ) << chrono::duration_cast<Ms>  ( mpp_durationCustomSteps[nb][0] % Secs( 1 ) ).count() << " ms" << endl);
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
