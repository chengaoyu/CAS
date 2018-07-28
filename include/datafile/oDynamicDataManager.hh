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
  \ingroup  datafile
  \brief    Declaration of class oDynamicDataManager
*/

#ifndef ODYNAMICDATAMANAGER_HH
#define ODYNAMICDATAMANAGER_HH 1

#define DYNAMIC_SWITCH_NOTHING 0
#define DYNAMIC_SWITCH_CONTINUE 1
#define DYNAMIC_SWITCH_DEFORMATION 2

#include "gVariables.hh"
#include "gOptions.hh"
#include "sOutputManager.hh"

class oImageDimensionsAndQuantification;
class vDataFile;

/*!
  \class   oDynamicDataManager
  \brief   This class gathers the information about the dynamic splitting of the data
  \details It contains the functions dedicated to the reading of user-provided informations about the dynamic dataset. \n
           It is specific to each iDataFile object. \n
           It holds several arrays related to the 4D splitting of the data, the current gates indices, as well as the related functions.
*/
class oDynamicDataManager
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   oDynamicDataManager constructor. 
               Initialize the member variables to their default values.
    */
    oDynamicDataManager();
    /*!
      \brief   oDynamicDataManager destructor. 
    */
    ~oDynamicDataManager(); 

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      oDynamicDataManager::InitDynamicData
      \param   a_nbRespGates
      \param   a_nbCardGates
      \param   a_pathTo4DDataFile : path to an ASCII file containing dynamic metadata regarding the acquisition
      \param   a_rmCorrFlag : indicate whether respiratory motion correction is enabled (1) or disabled (0) 
      \param   a_cmCorrFlag : indicate whether cardiac motion correction is enabled (1) or disabled (0) 
      \param   a_dmCorrFlag : indicate whether simultaneous respiratory and cardiac motion corrections are enabled (1) or disabled (0) 
      \param   a_pmCorrFlag : indicate whether involuntary patient motion correction is enabled (1) or disabled (0) 
      \brief   Main function for instanciation and initialization of the member variables and arrays. Call the specific initialization function depending of the type of dataset.
      \todo     warning in the documentation that Respiratory motion correction is not supported for cardiac-gated data in the current implementation
      \return  0 if success, positive value otherwise.
    */
    int InitDynamicData( int a_nbRespGates, int a_nbCardGates, const string& a_pathTo4DDataSplittingFile,
                         int a_rmMCorrFlag, int a_cMmCorrFlag, int a_dmCorrFlag, int a_pMotionCorrFlag );
    /*!
      \fn      oDynamicDataManager::SetDynamicSpecificQuantificationFactors
      \param   FLTNB** a2p_quantificationFactors : 2 dimensional [timeframe][gate] set of quantitative factors to update
      \brief   Compute gate-specific quantificative factors using the number of events within each gate, and update the quantitative factors passed in argument
      \todo    should be updated if we support multi-bed gated dataset
      \return  0 if success, positive value otherwise.
    */
    int SetDynamicSpecificQuantificationFactors(FLTNB** a2p_quantificationFactors);
    /*!
      \fn      oDynamicDataManager::CheckParameters
      \param   a_nbEvents : Number of events in the acquisition, used to check consistency with the gating metadata
      \brief   Check all mandatory parameters.
      \return  0 if success, positive value otherwise.
    */
    int CheckParameters(int64_t a_nbEvents);
    /*!
      \fn      oDynamicDataManager::ResetCurrentDynamicIndices
      \brief   Reset to 0 the multithreaded dynamic arrays gathering the indices of current frame, gates and involuntary motion
    */
    void ResetCurrentDynamicIndices();
    /*!
      \fn      oDynamicDataManager::DynamicSwitch
      \param   a_currentEventIndex : Index of the current event in the iterative reconstruction loop
      \param   a_currentTime : Timestamp of the current event in the iterative reconstruction loop
      \param   a_bed
      \param   a_th
      \brief   This function is called in the reconstruction event loop. It is used to check if the current event belongs to a new respiratory/cardiac/involuntary motion gate \n
               Increment the index reporting the current relative gate in this case
      \todo    Check implementation for cardiac gating and involuntary patient motion correction
      \return  0 if nothing to do,
               1 if the current timestamp is inferior to the start time of the first frame,
               2 if a deformation is required (gate has changed and motion correction is enabled)
    */
    int DynamicSwitch(int64_t a_index, uint32_t a_time, int a_bed, int a_th);
    
  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      oDynamicDataManager::SetVerbose
      \param   a_verboseLevel
      \brief   set verbosity
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;};
    /*!
      \fn      oDynamicDataManager::SetImageDimensionsAndQuantification
      \param   ap_ImageDimensionsAndQuantification
      \brief   set the pointer to the oImageDimensionsAndQuantification object
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ID = ap_ImageDimensionsAndQuantification;};
    /*!
      \fn      oDynamicDataManager::SetRespMotionFlagOn
      \brief   set the respiratory motion flag to 'true'
    */
    inline void SetRespMotionFlagOn()
           {m_rMotionCorrFlag = true;}
    /*!
      \fn      oDynamicDataManager::SetCardMotionFlagOn
      \brief   set the cardiac motion flag to 'true'
    */
    inline void SetCardMotionFlagOn()
           {m_cMotionCorrFlag = true;} 
    /*!
      \fn      oDynamicDataManager::SetPMotionFlagOn
      \brief   set the involuntary patient motion flag to 'true'
    */
    inline void SetPMotionFlagOn()
           {m_pMotionCorrFlag = true;} 
    /*!
      \fn      oDynamicDataManager::IsRespMotionEnabled()
      \return  true if respiratory motion is enabled, false otherwise
    */
    inline bool IsRespMotionEnabled()
           {return m_rMotionCorrFlag;}
    /*!
      \fn      oDynamicDataManager::IsCardMotionEnabled()
      \return  true if cardiac motion is enabled, false otherwise
    */
    inline bool IsCardMotionEnabled()
           {return m_cMotionCorrFlag;}
    /*!
      \fn      oDynamicDataManager::IsPMotionEnabled()
      \return  true if involuntary patient motion is enabled, false otherwise
    */
    inline bool IsPMotionEnabled()
           {return m_pMotionCorrFlag;}
    /*!
      \fn      oDynamicDataManager::GetCurrentPMotionIndex
      \param   a_th
      \brief   return the number of the current involuntary patient motion image to be used in the ImageSpace matrices, for this thread
      \return  the index of the current involuntary patient motion gated image if this methodology is enabled, 0 otherwise
    */
    inline int GetCurrentPMotionIndex(int a_th)
           {return mp_currentPMotionIndex[a_th] ;}
    /*!
      \fn      oDynamicDataManager::GetCurrentTimeFrame(int a_th)
      \param   a_th
      \brief   Return the index indicating the current time frame for this thread.
      \return  time frame index.
    */
    inline int GetCurrentTimeFrame(int a_th)
           {return mp_currentFrameIndex[a_th];}
    /*!
      \fn      oDynamicDataManager::GetCurrentRespGate(int a_th)
      \param   a_th
      \brief   Return the index indicating the current respiratory gate for this thread.
      \return  respiratory gate index.
    */
    inline int GetCurrentRespGate(int a_th)
           {return mp_currentRespGateIndex[a_th];}
    /*!
      \fn      oDynamicDataManager::GetCurrentRespImage
      \param   a_th
      \brief   return the number of the current respiratory gated image, to be used in the ImageSpace matrices for the reconstruction
      \details if respiratory motion correction is enabled in reconstruction, the function returns 0 as only one image is reconstructed for each gate
      \return  the index of the gated image related to the reconstruction for the specific thread
    */
    inline int GetCurrentRespImage(int a_th)
           {return m_rMotionCorrFlag ? 0 : mp_currentRespGateIndex[a_th];}
    /*!
      \fn      oDynamicDataManager::GetCurrentCardGate
      \param   a_th
      \brief   Return the index indicating the current cardiac gate for this thread.
      \return  cardiac gate index.
    */
    inline int GetCurrentCardGate(int a_th)
           {return mp_currentCardGateIndex[a_th];}
    /*!
      \fn      oDynamicDataManager::GetCurrentCardImage
      \param   a_th
      \brief   return the number of the current cardiac gated image, to be used in the ImageSpace matrices for the reconstruction
      \details if cardiac motion correction is enabled in reconstruction, the function returns 0 as only one image is reconstructed for each gate
      \return  the index of the current gated image related to the reconstruction for the specific thread
    */
    inline int GetCurrentCardImage(int a_th)
           {return m_cMotionCorrFlag ? 0 : mp_currentCardGateIndex[a_th];}
    /*!
      \fn      oDynamicDataManager::GetNb2ndMotImgsForLMS
      \brief   return the number of secundary motion (typically cardiac) images to be used in the ImageSpace matrices for the list-mode sensitivity image generation
      \return  the number of cardiac images required for the sensitivity image generation
    */
    inline int GetNb2ndMotImgsForLMS()
           {return m_nbCardGates;}
    /*!
      \fn      oDynamicDataManager::PMotionFirstIndex
      \param   a_fr = frame index
      \brief   if patient motion is enabled, return the first index of patient motion for the frame fr, 0 otherwise
               This is required in the ImageSpace matrices for the list-mode sensitivity image generation
      \return  the total number of respiratory gates
    */
    inline int GetPMotionFirstIndexForLMS(int a_fr)
           //{return m_pMotionCorrFlag ? GetNbIPatMotionSubsets() : m_nbRespGates ;}
           {return m_pMotionCorrFlag ? mp_framePMotionFirstIndex[a_fr] : 0 ;}
    /*!
      \fn      oDynamicDataManager::GetNb1stMotImgsForLMS
      \param   a_fr = frame index
      \brief   return the number of first motion (either respiratory or involuntary patient motion) images to be used in the ImageSpace matrices for the list-mode sensitivity image generation
      \return  the total number of respiratory gates
    */
    inline int GetNb1stMotImgsForLMS(int a_fr)
           //{return m_pMotionCorrFlag ? GetNbIPatMotionSubsets() : m_nbRespGates ;}
           {return m_pMotionCorrFlag ? mp_frameNbPMotionTriggers[a_fr] : m_nbRespGates ;}
    /*!
      \fn      oDynamicDataManager::GetNbIPatMotionSubsets
      \brief   return the number of involuntary patient motion transformations
      \return  the total number of involuntary patient motion triggers
    */
    inline int GetNbIPatMotionSubsets()
           {return m_nbPMotionTriggers+1;}

  // -------------------------------------------------------------------
  // Private member functions  
  private:
    /*!
      \fn      oDynamicDataManager::InitDynamicDataGating
      \param   a_pathToFile : path to an ASCII file containing dynamic metadata regarding the acquisition
      \brief   Initialisation of arrays containing informations about the data splitting and respiratory/cardiac gated reconstruction
      \todo    check the reconstruction with cardiac gated reconstruction. Some adjustements will be required for double gating
      \return  0 if success, positive value otherwise.
    */
    int InitDynamicDataGating(const string& a_pathToGateFile);
    /*!
      \fn      oDynamicDataManager::InitDynamicDataPatientMotion
      \param   a_pathToFile : path to an ASCII file containing dynamic metadata regarding the acquisition
      \brief   Initialisation of involuntary patient motion correction information, if any.
      \todo    Implementation currently IN PROGRESS. Will have to define how motion subset information is provided as it depends of each datafile, if the acquisition contain several bed steps.
      \return  0 if success, positive value otherwise.
    */
    int InitDynamicDataPatientMotion(const string& a_pathToFile);

  // -------------------------------------------------------------------
  // Data members
  private:
    oImageDimensionsAndQuantification*  mp_ID; /*!< Pointer to the oImageDimensionsAndQuantification object */
    int m_verbose;                             /*!< Verbosity */
    
    // Framing
    int* mp_currentFrameIndex;                 /*!< Multithreaded array [nb threads] containing the index of the current frame in use*/
    int m_nbTimeFrames;                        /*!< Number of time frames */
     
    // Respiratory gating
    bool m_respGatingFlag;                     /*!< Flag indicating if respiratory gating is enabled */
    bool m_rMotionCorrFlag;                    /*!< Flag indicating if respiratory motion correction is enabled */
    int m_nbRespGates;                         /*!< Number of gates for respiratory motion correction */
    int64_t** m2p_nbEventsPerRespGate;         /*!< 2 dimensional array [nb frames][nb resp gates] containing the number of events in each respiratory gates */
    int64_t** m2p_indexLastEventRespGate;      /*!< 2 dimensional array [nb frames][nb resp gates] containing the last event in the respiratory gates */
    int* mp_currentRespGateIndex;              /*!< Multithreaded array [nb threads] containing the index of the current respiratory gate in use*/
    
    // Cardiac gating
    bool m_cardGatingFlag;                     /*!< Flag indicating if cardiac gating is enabled */
    bool m_cMotionCorrFlag;                    /*!< Flag indicating if patient involuntary motion correction is enabled */
    int m_nbCardGates;                         /*!< Number of gates for cardiac motion correction */
    int64_t** m2p_nbEventsPerCardGate;         /*!< 2 dimensional array [nb frames][nb resp gates] containing the number of events in each cardiac gates. */
    int64_t** m2p_indexLastEventCardGate;      /*!< 2 dimensional array [nb frames][nb resp gates] containing the last event in the cardiac gates */
    int* mp_currentCardGateIndex;              /*!< Multithreaded array [nb threads] containing the index of the current cardiac gate in use */
    
    // Involuntary motion
    bool m_pMotionCorrFlag;                    /*!< Flag indicating if patient involuntary motion correction is enabled */
    int m_nbPMotionTriggers;                   /*!< Total number of triggers for patient involuntary motion correction */
    uint16_t* mp_framePMotionFirstIndex;       /*!< First index of patient involuntary motion triggers by frame */
    uint16_t* mp_frameNbPMotionTriggers;       /*!< Number of patient involuntary motion triggers by frame */
    uint32_t* mp_listPMotionTriggers;          /*!< Array containing the timestamp of each trigger of the patient involuntary motion correction*/
    int* mp_currentPMotionIndex;               /*!< Multithreaded array [nb threads] containing the index for patient involuntary motion correction*/
};
#endif
