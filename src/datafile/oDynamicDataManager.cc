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
  \ingroup datafile

  \brief Implementation of class oDynamicDataManager
*/

#include "oDynamicDataManager.hh"
#include "oImageDimensionsAndQuantification.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oDynamicDataManager::oDynamicDataManager()
{
  mp_ID = NULL;
  m_verbose = -1;
  m_nbTimeFrames = -1;
  mp_currentFrameIndex = NULL;
  m_respGatingFlag = false;
  m_rMotionCorrFlag = false;
  m_nbRespGates = -1;
  m2p_nbEventsPerRespGate = NULL;
  m2p_indexLastEventRespGate = NULL;
  mp_currentRespGateIndex = NULL;
  m_cardGatingFlag = false;
  m_cMotionCorrFlag = false;
  m_nbCardGates = -1;
  m2p_nbEventsPerCardGate = NULL;
  m2p_indexLastEventCardGate = NULL;
  mp_currentCardGateIndex = NULL;
  m_pMotionCorrFlag = false;
  m_nbPMotionTriggers = -1;
  mp_listPMotionTriggers = NULL;
  mp_frameNbPMotionTriggers = NULL;
  mp_framePMotionFirstIndex = NULL;
  mp_currentPMotionIndex = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oDynamicDataManager::~oDynamicDataManager()
{
  for (int fr=0; fr<m_nbTimeFrames; fr++)
  {
    if (m2p_nbEventsPerRespGate && m2p_nbEventsPerRespGate[fr]) 
      delete m2p_nbEventsPerRespGate[fr];
    if (m2p_nbEventsPerCardGate && m2p_nbEventsPerCardGate[fr]) 
      delete m2p_nbEventsPerCardGate[fr];
    if (m2p_indexLastEventRespGate && m2p_indexLastEventRespGate[fr]) 
      delete m2p_indexLastEventRespGate[fr];
    if (m2p_indexLastEventCardGate && m2p_indexLastEventCardGate[fr]) 
      delete m2p_indexLastEventCardGate[fr];
  }

  if (m2p_nbEventsPerRespGate!=NULL) delete m2p_nbEventsPerRespGate;
  if (m2p_nbEventsPerCardGate!=NULL) delete m2p_nbEventsPerCardGate;
  if (m2p_indexLastEventRespGate!= NULL) delete m2p_indexLastEventRespGate;
  if (m2p_indexLastEventCardGate!= NULL) delete m2p_indexLastEventCardGate;
  if (mp_currentFrameIndex!=NULL) delete mp_currentFrameIndex;
  if (mp_currentRespGateIndex!=NULL) delete mp_currentRespGateIndex;
  if (mp_currentCardGateIndex!=NULL) delete mp_currentCardGateIndex;
  if (mp_currentPMotionIndex!=NULL) delete mp_currentPMotionIndex;
  if (mp_frameNbPMotionTriggers!=NULL) delete mp_frameNbPMotionTriggers;
  if (mp_framePMotionFirstIndex!=NULL) delete mp_framePMotionFirstIndex;
  if (mp_listPMotionTriggers!=NULL) delete mp_listPMotionTriggers;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::InitDynamicData( int a_nbRespGates, int a_nbCardGates, const string& a_pathTo4DDataFile,
                                          int a_rmCorrFlag, int a_cmCorrFlag, int a_dmCorrFlag, int a_pmCorrFlag)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicData() -> Initialize dynamic data management" << endl);

  // Initialization
  m_nbTimeFrames = mp_ID->GetNbTimeFrames();
  m_nbRespGates = a_nbRespGates;
  m_nbCardGates = a_nbCardGates;
  m_nbPMotionTriggers = 0;
  
  if (m_nbRespGates > 1) m_respGatingFlag = true;
  if (m_nbCardGates > 1) m_cardGatingFlag = true;
  m_rMotionCorrFlag = a_rmCorrFlag;
  m_cMotionCorrFlag = a_cmCorrFlag;
  m_pMotionCorrFlag = a_pmCorrFlag;
  if(a_dmCorrFlag)
  {
    m_rMotionCorrFlag = true;
    m_cMotionCorrFlag = true;
  }

  // Instanciation & initialization for current indices (used during the loop on events to know which frame/gate the event belongs to)
  mp_currentFrameIndex     = new int[mp_ID->GetNbThreadsForProjection()];
  mp_currentRespGateIndex  = new int[mp_ID->GetNbThreadsForProjection()];
  mp_currentCardGateIndex  = new int[mp_ID->GetNbThreadsForProjection()];
  mp_currentPMotionIndex = new int[mp_ID->GetNbThreadsForProjection()];
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    mp_currentFrameIndex    [0] = 0;
    mp_currentRespGateIndex [0] = 0;
    mp_currentCardGateIndex [0] = 0;
    mp_currentPMotionIndex[0] = 0;
  }
  
  // Instanciation & initialization of number of events per frame for gated reconstructions
  m2p_nbEventsPerRespGate = new int64_t*[m_nbTimeFrames];
  m2p_nbEventsPerCardGate = new int64_t*[m_nbTimeFrames];
  m2p_indexLastEventRespGate = new int64_t*[m_nbTimeFrames];
  m2p_indexLastEventCardGate = new int64_t*[m_nbTimeFrames];
  
  for (int fr=0; fr<m_nbTimeFrames; fr++)
  {
    m2p_nbEventsPerRespGate[fr] = new int64_t[m_nbRespGates];
    m2p_nbEventsPerCardGate[fr] = new int64_t[m_nbRespGates*m_nbCardGates];
    m2p_indexLastEventRespGate[fr] = new int64_t[m_nbRespGates];
    m2p_indexLastEventCardGate[fr] = new int64_t[m_nbRespGates*m_nbCardGates];
  
    for (int rg=0; rg<m_nbRespGates; rg++)
    {
      m2p_nbEventsPerRespGate[fr][rg] = -1;
      m2p_indexLastEventRespGate[fr][rg] = -1;
      
      for (int cg=0; cg<m_nbCardGates; cg++)
      {
        m2p_nbEventsPerCardGate[fr][rg*m_nbCardGates+cg] = -1;
        m2p_indexLastEventCardGate[fr][rg*m_nbCardGates+cg] = -1;
      }
    }
  }

  // Check coherence of motion initialization
  // Throw error if cardiac gating is enabled alongside respiratory motion correction (not available)
  // TODO : warning in the documentation that Respiratory motion correction is not supported for cardiac-gated data in the current implementation
  if (m_rMotionCorrFlag && m_cardGatingFlag)
  {
    Cerr("***** oDynamicDataManager::InitDynamicData() -> Respiratory motion correction is enabled for cardiac-gated data. This is not supported in the current implementation  !" << endl);
    return 1;
  }

  // Check iPat motion vs physiological motion
  if (a_pmCorrFlag && (m_rMotionCorrFlag || m_cMotionCorrFlag) )
  {
    Cerr("***** oDynamicDataManager::InitDynamicData() -> Involuntary patient image-based motion correction cannot be used along with respiratory or cardiac motion correction  !" << endl);
    return 1;
  }

  // Incorrect motion initialization
  if ( (a_rmCorrFlag + a_cmCorrFlag + a_dmCorrFlag) >1 )
  {
    Cerr("***** oDynamicDataManager::InitDynamicData() -> Incorrect image-based motion initialization."<<endl);
    Cerr("*****                                          Please use -rm for respiratory motion, -cm for cardiac motion, or -dm for double respiratory and cardiac motion correction  !" << endl);
    return 1;
  }
  
  // Initialization for respiratory/cardiac gated reconstruction
  // Note : for Analytic Projection, gating flag could be enabled in order to project an image containing several gates, but no description file is required
  //        InitDynamicDataGating() is restricted to reconstruction
  //        Check on the existence of a_pathTo4DDataFile during reconstruction is performed onnly in Grecon
  if ((m_respGatingFlag || m_cardGatingFlag) && !a_pathTo4DDataFile.empty() )
  {
    if(m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicData() -> Initializing data for gated reconstruction... " << endl);
      
    if (InitDynamicDataGating(a_pathTo4DDataFile))
    {
      Cerr("***** oDynamicDataManager::InitDynamicData() -> A problem occured during the dynamic gating initialization !" << endl);
      return 1;
    }
  }

  // Initialization with involuntary patient motion corrected reconstruction
  if (m_pMotionCorrFlag && InitDynamicDataPatientMotion(a_pathTo4DDataFile) )
  {
    Cerr("***** oDynamicDataManager::InitDynamicData() -> A problem occured during the patient involuntary motion correction initialization !" << endl);
    return 1;
  }

  // Some feedback
  if (m_verbose>=VERBOSE_DETAIL)
  {
    if (m_respGatingFlag)  Cout("oDynamicDataManager::InitDynamicData() -> Respiratory gating enabled" << endl);
    if (m_rMotionCorrFlag) Cout("oDynamicDataManager::InitDynamicData() -> Respiratory image-based motion correction enabled" << endl);
    if (m_cardGatingFlag)  Cout("oDynamicDataManager::InitDynamicData() -> Cardiac gating enabled" << endl);
    if (m_cMotionCorrFlag) Cout("oDynamicDataManager::InitDynamicData() -> Cardiac image-based motion correction enabled" << endl);
    if (m_pMotionCorrFlag) Cout("oDynamicDataManager::InitDynamicData() -> Involuntary image-based patient motion correction enabled" << endl);
  }
    
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::InitDynamicDataGating(const string& a_pathToFile)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataGating() ... " << endl);
    
  // Respiratory gating enabled
  if (m_respGatingFlag)
  {
    // Read information about respiratory gating (number of events per respiratory gate per frame)
    if ( ReadDataASCIIFile(a_pathToFile, 
                          "respiratory_gates_data_splitting", 
                           m2p_nbEventsPerRespGate, 
                           m_nbRespGates, 
                           m_nbTimeFrames, 
                           KEYWORD_MANDATORY) )
    {
      Cerr("***** oDynamicDataManager::InitDynamicDataGating() -> Error while trying to read 'respiratory_gates_data_splitting' in file '" << a_pathToFile << "' !" << endl);
      return 1;
    }
    
    // Get the index of the last event of each respiratory gate using the number of events inside each gate
    uint64_t event_number_sum = 0;
    if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataGating() :" << endl);
            
    for (int fr=0; fr<m_nbTimeFrames; fr++)
      for (int rg=0; rg<m_nbRespGates; rg++)
      { 
        m2p_indexLastEventRespGate[fr][rg] += m2p_nbEventsPerRespGate[fr][rg] + event_number_sum;
        event_number_sum += m2p_nbEventsPerRespGate[fr][rg]; 
        if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Number of events for frame #" << fr << ", respiratory gate #" << rg << " = " << m2p_nbEventsPerRespGate[fr][rg] << endl);
      } 
  }

  // Cardiac gating enabled
  if (m_cardGatingFlag)
  {
    // Read information about cardiac gating (number of events per cardiac gate per respiratory gates times frames)
    if ( ReadDataASCIIFile(a_pathToFile, 
                          "cardiac_gates_data_splitting", 
                           m2p_nbEventsPerCardGate,
                           m_nbCardGates, 
                           m_nbTimeFrames*m_nbRespGates, 
                           KEYWORD_MANDATORY) )
    {
      Cerr("***** oDynamicDataManager::InitDynamicDataGating() -> Error while trying to read 'cardiac_gates_data_splitting' in file '" << a_pathToFile << "' !" << endl);
      return 1;
    }
    // Get the index of the last event of each cardiac gate using the number of events inside each gate
    uint64_t event_number_sum = 0;
    if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataGating() :" << endl);
        
    for (int fr=0; fr<m_nbTimeFrames; fr++)
      for (int g=0; g<m_nbRespGates*m_nbCardGates; g++)
      { 
        m2p_indexLastEventCardGate[fr][g] += m2p_nbEventsPerCardGate[fr][g] + event_number_sum;
        event_number_sum += m2p_nbEventsPerCardGate[fr][g]; 
        if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Number of events for frame #" << fr << ", cardiac gate #" << g << " = " << m2p_nbEventsPerCardGate[fr][g] << endl);
      } 
  }
  
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::InitDynamicDataPatientMotion(const string& a_pathToFile)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataPatientMotion() ... " << endl);

  // First get the number of time triggers into the dynamic file
  if (ReadDataASCIIFile(a_pathToFile, 
                       "nb_involuntary_motion_triggers", 
                        &m_nbPMotionTriggers, 
                        1, 
                        KEYWORD_MANDATORY))
    {
      Cerr("***** oDynamicDataManager::InitDynamicDataPatientMotion() -> Involuntary motion correction enabled but" 
        << " number of triggers could not be found in the dynamic file '" << a_pathToFile << "' !" << endl);
      return 1;
    }
  // Allocate time triggers and read them from the dynamic file
  mp_listPMotionTriggers = (uint32_t*)malloc(m_nbPMotionTriggers*sizeof(uint32_t));
  mp_frameNbPMotionTriggers = (uint16_t*)malloc(m_nbTimeFrames*sizeof(uint16_t));
  mp_framePMotionFirstIndex = (uint16_t*)malloc(m_nbTimeFrames*sizeof(uint16_t));
  for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
  {
    mp_frameNbPMotionTriggers[fr] = 1;
    mp_framePMotionFirstIndex[fr] = 0;
  }
  
  // Get timestamp of motion trigger
  if (ReadDataASCIIFile(a_pathToFile, 
                      "list_involuntary_motion_triggers", 
                        mp_listPMotionTriggers, 
                        m_nbPMotionTriggers, 
                        KEYWORD_MANDATORY))
  {
    Cerr("***** oDynamicDataManager::InitDynamicDataPatientMotion() -> Involuntary motion correction enabled but list of triggers"
      << " not found in the dynamic file '" << a_pathToFile << "' !" << endl);
    return 1;
  }

  // Recover the number of transformation required for each frame
  uint16_t fr=0;
  
  for(int t=0 ; t<m_nbPMotionTriggers ; t++)
  {
    while(mp_listPMotionTriggers[t] > mp_ID->GetFrameTimeStopInMs(0,fr) )
    {
      fr++;
      mp_framePMotionFirstIndex[fr]=t;
      mp_frameNbPMotionTriggers[fr]=t+1;
    }
    
    mp_frameNbPMotionTriggers[fr]++;
    
    // TODO check this (motion trigger after last frameTimeStop)
    if(fr >= m_nbTimeFrames)
      break;
  }
   
  if(fr<m_nbTimeFrames)
    for(int f=fr ; f<m_nbTimeFrames ; f++)
      mp_frameNbPMotionTriggers[f] = mp_frameNbPMotionTriggers[fr];

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::SetDynamicSpecificQuantificationFactors(FLTNB** a2p_quantificationFactors)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_NORMAL) Cout("oDynamicDataManager::SetDynamicSpecificQuantificationFactors() ... " << endl);
    
  // COMPUTE GATE-SPECIFIC QUANTITATIVE FACTORS
  if (m_nbRespGates>1 || m_nbCardGates>1)
  {
    for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
    {
      // We have cardiac gates and don't correct for motion -> quantification factors required for cardiac gates
      if (m_cardGatingFlag && !m_cMotionCorrFlag)
      {
        uint64_t total_events_in_frame = 0;
        for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++) total_events_in_frame += m2p_nbEventsPerCardGate[fr][g];

        if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Cardiac gating correction factors :" << endl);
        for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
        {
          a2p_quantificationFactors[fr][g] *= (FLTNB)total_events_in_frame/(FLTNB)m2p_nbEventsPerCardGate[fr][g];
          if (m_verbose>=VERBOSE_DETAIL) Cout("      Frame #" << fr << ", cardiac gate #" << g << " = " << a2p_quantificationFactors[fr][g] << endl);
        }
      }
      
      // We have resp gates and don't correct for resp motion (and no independent cardiac gate reconstruction) -> quantification factors required for cardiac gates
      else if(m_respGatingFlag && !m_rMotionCorrFlag)
      {
        uint64_t total_events_in_frame = 0;
        for(int rg=0 ; rg<m_nbRespGates ; rg++) total_events_in_frame += m2p_nbEventsPerRespGate[fr][rg];

        if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Respiratory gating correction factors :" << endl);
        for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
        {
          int rg = int(g/m_nbCardGates);
          a2p_quantificationFactors[fr][g] *= (FLTNB)total_events_in_frame/(FLTNB)m2p_nbEventsPerRespGate[fr][rg];
          if (m_verbose>=VERBOSE_DETAIL) Cout("      Frame #" << fr << ", gate #" << g << " = " << a2p_quantificationFactors[fr][g] << endl);
        }
      }
    }
  }
  
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::CheckParameters(int64_t a_nbEvents)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::CheckParameters() -> Check parameters for dynamic data settings" << endl);
  // Check image dimensions
  if (mp_ID==NULL)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  // Check time frames
  if (m_nbTimeFrames<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong number of time frames !" << endl);
    return 1;
  }
  // Check resp gates
  if (m_nbRespGates<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong number of respiratory gates !" << endl);
    return 1;
  }
  // Check card gates
  if (m_nbCardGates<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong number of respiratory gates !" << endl);
    return 1;
  }
  // Check involuntary motion triggers
  if (m_nbPMotionTriggers<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong number of involuntary motion subsets provided !" << endl);
    return 1;
  }
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong verbosity level provided !" << endl);
    return 1;
  }  
  
  // Feedback regarding the enabled/disabled options for the reconstruction
  if (m_verbose>=2)
  {
    if (m_respGatingFlag) Cout("oDynamicDataManager::CheckParameters() -> Respiratory gating is enabled" << endl);
    if (m_rMotionCorrFlag) Cout("oDynamicDataManager::CheckParameters() -> Respiratory motion correction enabled" << endl);
    if (m_cardGatingFlag) Cout("oDynamicDataManager::CheckParameters() -> Cardiac gating is enabled" << endl);
    if (m_cMotionCorrFlag) Cout("oDynamicDataManager::CheckParameters() -> Cardiac motion correction is enabled" << endl);
    if (m_pMotionCorrFlag) Cout("oDynamicDataManager::CheckParameters() -> Involuntary motion correction is enabled" << endl);
  }

  // Check consistency between number of events in the datafile and the potential splitting of the dynamic data into respiratory or cardiac gates (if any gating is enabled)
  if (m_respGatingFlag)
  {
    int last_fr = m_nbTimeFrames-1;
    int last_rg = m_nbRespGates-1;
    if (m2p_indexLastEventRespGate[last_fr][last_rg]+1 != a_nbEvents)
    {
      Cerr("***** oDynamicDataManager::CheckParameters() -> Problem while checking consistency of dynamic data !" << endl
        << "                                                The number of events in the datafile (" << a_nbEvents
        << ") is different from the total number of events in respiratory gates (" << m2p_indexLastEventRespGate[last_fr][last_rg] << ") as initialized in the gating file !" << endl);
      return 1;
    }
  }
  if (m_cardGatingFlag)
  {
    int last_fr = m_nbTimeFrames-1;
    int last_rg = m_nbCardGates-1;
    if (m2p_indexLastEventCardGate[last_fr][last_rg]+1 != a_nbEvents)
    {
      Cerr("***** oDynamicDataManager::CheckParameters() -> Problem while checking consistency of dynamic data !" << endl
        << "                                                The number of events in the datafile (" << a_nbEvents
        << ") is different to the total number of events in cardiac gates (" << m2p_indexLastEventRespGate[last_fr][last_rg] << ") as initialized in the gating file !" << endl);
      return 1;
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oDynamicDataManager::ResetCurrentDynamicIndices()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Reset indices for each thread (this is a thread-safe implementation)
  for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++)
  {
    // We initialize the frame index to -1 because at the beginning of the events loop, we don't know yet
    // if we are in the first frame (the user can reconstruct only a part of the whole datafile).
    mp_currentFrameIndex[th] = -1;
    mp_currentPMotionIndex[th] = 0;
    mp_currentRespGateIndex[th] = 0;
    mp_currentCardGateIndex[th] = 0;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::DynamicSwitch(int64_t a_currentEventIndex, uint32_t a_currentTime, int a_bed, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  // The value that will be returned
  int return_value = DYNAMIC_SWITCH_NOTHING;

  // -------------------------------------------------------------------------------------------------
  // Step 1: involuntary motion management
  // -------------------------------------------------------------------------------------------------

  // Only do this if the involuntary motion correction is enabled (meaning we must proceed to a deformation)
  if (m_pMotionCorrFlag)
  {
    // Search if we passed one of the next motion triggers (starting at current index)
    for (int mt=mp_currentPMotionIndex[a_th]; mt<m_nbPMotionTriggers; mt++)
    {

      // If we passed this trigger, set the return value to DEFORMATION. However, we continue the loop
      // in the case where we passed multiple triggers.
      if (a_currentTime>=mp_listPMotionTriggers[mt])
      {
        mp_currentPMotionIndex[a_th] = mt+1;
        #ifdef CASTOR_DEBUG
        if (m_verbose>=VERBOSE_DEBUG_EVENT) Cout("oDynamicDataManager::DynamicSwitch() -> Thread " << a_th << ", gate for patient motion correction switched to " << mp_currentPMotionIndex[a_th] << endl
                             << "                                        at event #" << a_currentEventIndex << ", timestamp = " << a_currentTime << endl);
        #endif
        return_value = DYNAMIC_SWITCH_DEFORMATION;
      }
    }
  }

  // -------------------------------------------------------------------------------------------------
  // Step 2: frame management
  // -------------------------------------------------------------------------------------------------

  // Special case if the frame number is still -1
  if (mp_currentFrameIndex[a_th]==-1)
  {
    // We check if we are before the first frame, in this case we directly return a CONTINUE signal to go to the next event
    if (a_currentTime<mp_ID->GetFrameTimeStartInMs(a_bed,0)) return DYNAMIC_SWITCH_CONTINUE;
    // Otherwise, we now are at least in the first frame (which will be updated right after this 'if' section)
    else mp_currentFrameIndex[a_th] = 0;
  }

  // A boolean to know later if the frame index has changed
  bool frame_has_changed = false;

  
  // Now we search for the first frame index in which the event belongs, starting from the current frame index. Note that we do that
  // this way (not simply incrementing the frame index) because we want to deal with the case where the next event managed by this
  // thread could possibly skip multiple frames at once.
  for (int fr=mp_currentFrameIndex[a_th]; fr<m_nbTimeFrames; fr++)
  {
    // If the current time is less than the time stop of the frame 'fr', then the event belongs to this frame
    if (a_currentTime<mp_ID->GetFrameTimeStopInMs(a_bed,fr))
    {
      // Test if the frame has actually changed
      if (mp_currentFrameIndex[a_th]!=fr)
      {
        // Set the boolean to true
        frame_has_changed = true;
        // Update the current frame index
        mp_currentFrameIndex[a_th] = fr;
      }
      break;
    }
  }

  #ifdef CASTOR_DEBUG
  if (frame_has_changed && m_verbose >=VERBOSE_DEBUG_EVENT)
    Cout("oDynamicDataManager::DynamicSwitch() -> Thread " << a_th << ", frame switched to " << mp_currentFrameIndex[a_th] << endl);
  #endif

  // -------------------------------------------------------------------------------------------------
  // Step 3: respiratory gating management
  // -------------------------------------------------------------------------------------------------

  // A boolean to know later if the respiratory index has changed
  bool resp_gate_has_changed = false;

  // Do this only if respiratory gating is enabled
  if (m_respGatingFlag)
  {
    // Test if the frame index has changed
    if (frame_has_changed)
    {
      // Reset the respiratory gate index
      mp_currentRespGateIndex[a_th] = 0;
      // No deformation signal need to be sent, as the forward/backward image matrices for the next frame will be in the reference position as required
    }
    // For this frame, search the first gate (from the current gate index) for which the provided event index is below the event-index-stop
    for (int rg=mp_currentRespGateIndex[a_th]; rg<m_nbRespGates; rg++)
    {
      // If the current event index is below the last event of this gate, then the event belongs to this gate
      // (We won't enter again in the if statement due to the flag setting to true)
      if (a_currentEventIndex<=m2p_indexLastEventRespGate[mp_currentFrameIndex[a_th]][rg] && resp_gate_has_changed == false)
      {
        // Test if the gate has actually changed
        if (mp_currentRespGateIndex[a_th]!=rg)
        {
          // Update the current gate index
          mp_currentRespGateIndex[a_th] = rg;
          // Verbose
          #ifdef CASTOR_DEBUG
          if (m_verbose>=VERBOSE_DEBUG_EVENT) Cout("oDynamicDataManager::DynamicSwitch() -> Thread " << a_th << ", respiratory gate switch to " << mp_currentRespGateIndex[a_th]
                                                                                   << " on event " << a_currentEventIndex << endl
                                                                                   << "current frame : " << mp_currentFrameIndex[a_th] << endl
                                                                                   << "current respiratory gate index " << mp_currentRespGateIndex[a_th] << endl);
          #endif
          // If motion correction is enabled, then we should return a DEFORMATION signal
          if (m_rMotionCorrFlag) return_value = DYNAMIC_SWITCH_DEFORMATION;
        }
        // Set the boolean to true
        resp_gate_has_changed = true;
      }
    }
  }
  
  // -------------------------------------------------------------------------------------------------
  // Step 4: cardiac gating management
  // -------------------------------------------------------------------------------------------------

  // A boolean to know later if the respiratory index has changed
  bool card_gate_has_changed = false;
  
  // Do this only if cardiac gating is enabled
  if (m_cardGatingFlag)
  {
    // Test if the frame or respiratory index have changed
    if (frame_has_changed || resp_gate_has_changed)
    {
      // Reset the cardiac gate index
      mp_currentCardGateIndex[a_th] = 0;
      // Thus if one apply cardiac motion correction, we should return a DEFORMATION signal
      if (m_cMotionCorrFlag) return_value = DYNAMIC_SWITCH_DEFORMATION;
    }
    // For this frame and respiratory gate, search the first gate (from the current gate index) for which the provided event index is below the event-index-stop
    for (int cg=mp_currentCardGateIndex[a_th]; cg<m_nbCardGates; cg++)
    {
      // If the current event index is below the event-index-stop of this gate, then the event belongs to this gate
      if (a_currentEventIndex<m2p_indexLastEventCardGate[mp_currentFrameIndex[a_th]][mp_currentRespGateIndex[a_th]*m_nbCardGates+cg]  && card_gate_has_changed == false)
      {
        // Test if the gate has actually changed
        if (mp_currentCardGateIndex[a_th]!=cg)
        {
          // Update the current gate index
          mp_currentCardGateIndex[a_th] = cg;
          // Verbose
          #ifdef CASTOR_DEBUG
          if (m_verbose>=VERBOSE_DEBUG_EVENT) Cout("oDynamicDataManager::DynamicSwitch() -> Thread " << a_th << ", cardiac gate switch to " <<  mp_currentCardGateIndex[a_th]
                                                                                   << " on event " << a_currentEventIndex << endl
                                                                                   << "current frame : " << mp_currentFrameIndex[a_th] << endl);
          #endif
          // If motion correction is enabled, then we should return a DEFORMATION signal
          if (m_cMotionCorrFlag) return_value = DYNAMIC_SWITCH_DEFORMATION;
        }
        // Set the boolean to true
        card_gate_has_changed = true;
      }
    }
  }
  
  // -------------------------------------------------------------------------------------------------
  // Step 5: just return the value !
  // -------------------------------------------------------------------------------------------------

  // Return the status of the dynamic switch
  return return_value;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
