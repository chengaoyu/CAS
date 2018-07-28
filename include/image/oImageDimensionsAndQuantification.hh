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
  \ingroup  image
  \brief    Declaration of class oImageDimensionsAndQuantification
*/

#ifndef OIMAGEDIMENSIONSANDQUANTIFICATION_HH
#define OIMAGEDIMENSIONSANDQUANTIFICATION_HH 1

#include "gVariables.hh"
#include "sOutputManager.hh"
#include "gOptions.hh"
#include "oDynamicDataManager.hh"

/*!
  \class   oImageDimensionsAndQuantification
  \brief   This class is designed to manage all dimensions and quantification related stuff
  \details This class gather all dimensions information as well as quantification. It also
           manages the oDynamicDataManager which itself manages the dynamic data. Most
           classes in the project have a pointer to this class.
*/
class oImageDimensionsAndQuantification
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public oImageDimensionsAndQuantification::oImageDimensionsAndQuantification()
      \brief   The constructor of oImageDimensionsAndQuantification
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    oImageDimensionsAndQuantification();
    /*!
      \fn      public oImageDimensionsAndQuantification::~oImageDimensionsAndQuantification()
      \brief   The destructor of oImageDimensionsAndQuantification
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were build by this class.
    */
    ~oImageDimensionsAndQuantification();


  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public void oImageDimensionsAndQuantification::SetDefault()
      \brief   A function used to set number of threads and MPI instances to 1 and bypass the CheckParameters() function
      \details The m_checked flag is set to true
    */
    void SetDefault();
    /*!
      \fn      public int oImageDimensionsAndQuantification::CheckParameters()
      \brief   A function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int oImageDimensionsAndQuantification::Initialize()
      \brief   A function used to initialize all that is needed
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. It initializes image dimensions, dynamic basis functions and
               quantification factors.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int Initialize();
    /*!
      \fn      public void oImageDimensionsAndQuantification::CheckNumberOfProjectionThreadsConsistencyWithDataFileSize()
      \param   vDataFile** a2p_DataFile
      \details Reduce the number of projection threads to the lowest number of events in all datafiles, if the latter is below the former.
    */
    void CheckNumberOfProjectionThreadsConsistencyWithDataFileSize(vDataFile** a2p_DataFile);
    /*!
      \fn      public int oImageDimensionsAndQuantification::DealWithBedPositions()
      \param   vDataFile** a2p_DataFile
      \brief   Deal with provided or default bed relative positions.
      \details If relative bed positions are provided from datafiles, then reposition them with respect to the CASToR referential.
               Otherwise, take the default bed displacement from the scanner and compute the bed relative positions.
               Finally, the bed positions used in the reconstruction are stored in the mp_bedPositions member.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int DealWithBedPositions(vDataFile** a2p_DataFile);

  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int oImageDimensionsAndQuantification::InitializeFramingAndQuantification()
      \brief   A function used to initialize the framing and quantification tables
      \details This function is called by the Initialize() function.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int InitializeFramingAndQuantification();
    /*!
      \fn      private int oImageDimensionsAndQuantification::InitializeTimeBasisFunctions()
      \brief   A function used to initialize the time basis functions
      \details This function is called by the Initialize() function. If provided, it parses the file
               containing the time basis functions coefficients.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int InitializeTimeBasisFunctions();
    /*!
      \fn      private int oImageDimensionsAndQuantification::InitializeRespBasisFunctions()
      \brief   A function used to initialize the respiratory basis functions
      \details This function is called by the Initialize() function. If provided, it parses the file
               containing the respiratory basis functions coefficients.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int InitializeRespBasisFunctions();
    /*!
      \fn      private int oImageDimensionsAndQuantification::InitializeCardBasisFunctions()
      \brief   A function used to initialize the cardiac basis functions
      \details This function is called by the Initialize() function. If provided, it parses the file
               containing the cardiac basis functions coefficients.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int InitializeCardBasisFunctions();
    /*!
      \fn      private int oImageDimensionsAndQuantification::InitializeIgnoredCorrections()
      \brief   A function used to initialize the ignored corrections
      \details This function is called by the Initialize() function. It parses the ignored correction
               string (correction keywords separated by commas) and set the associated boolean members
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int InitializeIgnoredCorrections();


  // -----------------------------------------------------------------------------------------
  // Public dynamic data management functions
  public:
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::ResetCurrentDynamicIndices()
      \brief   Call the eponym function from the oDynamicDataManager class using the member object
    */
    inline void ResetCurrentDynamicIndices()
           {mp_DynamicDataManager->ResetCurrentDynamicIndices();}
    /*!
      \fn      public int oImageDimensionsAndQuantification::InitDynamicData()
      \param   a_pathTo4DDataSplittingFile
      \param   a_respMotionCorrectionFlag
      \param   a_cardMotionCorrectionFlag
      \param   a_doubleMotionCorrectionFlag
      \param   a_invMotionCorrectionFlag
      \param   a_nbRespGates
      \param   a_nbCardGates
      #\todo    we may get the information about rebinned datafile or not from the datafile header 
      \brief   Call the eponym function from the oDynamicDataManager object in order to initialize its data.
      \return  0 is success, positive value otherwise
    */
    int InitDynamicData( string a_pathTo4DDataSplittingFile,
                         int a_respMotionCorrectionFlag, int a_cardMotionCorrectionFlag, int a_doubleMotionCorrectionFlag, int a_invMotionCorrectionFlag,
                         int a_nbRespGates, int a_nbCardGates );
    /*!
      \fn      public int oImageDimensionsAndQuantification::CheckDynamicParameters()
      \param   int64_t a_nbEvents
      \brief   Call the eponym function from the oDynamicDataManager object in order to check its parameters.
      \return  0 is success, positive value otherwise
    */
    int CheckDynamicParameters(int64_t a_nbEvents);
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::DynamicSwitch()
      \param   a_currentEventIndex
      \param   a_currentTime
      \param   a_bed
      \param   a_th
      \brief   Call the eponym function from the oDynamicDataManager class using the member object
      \return  An integer reflecting the dynamic evoluation (gate/frame change, involuntary motion)
    */
    inline int DynamicSwitch( int64_t a_currentEventIndex, uint32_t a_currentTime, int a_bed, int a_th )
           {return mp_DynamicDataManager->DynamicSwitch(a_currentEventIndex, a_currentTime, a_bed, a_th);}


  // -----------------------------------------------------------------------------------------
  // Public Get & Set functions related to the dynamic data management
  public:
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetRespMotionFlagOn()
      \brief   Call the eponym function from the oDynamicDataManager class using the member object
    */
    inline void SetRespMotionFlagOn()
           {mp_DynamicDataManager->SetRespMotionFlagOn();}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetCardMotionFlagOn()
      \brief   Call the eponym function from the oDynamicDataManager class using the member object
    */
    inline void SetCardMotionFlagOn()
           {mp_DynamicDataManager->SetCardMotionFlagOn();}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetPMotionFlagOn()
      \brief   Call the eponym function from the oDynamicDataManager class using the member object
    */
    inline void SetPMotionFlagOn()
           {mp_DynamicDataManager->SetPMotionFlagOn();}   
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetCurrentRespGate()
      \param   a_th
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the current respiratory gate index
    */
    inline int GetCurrentRespGate(int a_th)
           {return mp_DynamicDataManager->GetCurrentRespGate(a_th);}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetCurrentCardGate()
      \param   a_th
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the current cardiac gate index
    */
    inline int GetCurrentCardGate(int a_th)
           {return mp_DynamicDataManager->GetCurrentCardGate(a_th);}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::PMotionFirstIndex
      \param   a_fr = frame index
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the first index for LMS
    */
    inline int GetPMotionFirstIndexForLMS(int a_fr)
           {return mp_DynamicDataManager->GetPMotionFirstIndexForLMS(a_fr);}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNb1stMotImgsForLMS()
      \param   a_fr = frame index
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the number of respiratory images for sensitivity
    */
    inline int GetNb1stMotImgsForLMS(int a_fr)
           {return mp_DynamicDataManager->GetNb1stMotImgsForLMS(a_fr);}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNb2ndMotImgsForLMS()
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the number of cardiac images for sensitivity
    */
    inline int GetNb2ndMotImgsForLMS()
           {return mp_DynamicDataManager->GetNb2ndMotImgsForLMS();}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetCurrentRespImage()
      \param   a_th
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the index of the current respiratory image for reconstruction
    */
    inline int GetCurrentRespImage(int a_th)
           {return mp_DynamicDataManager->GetCurrentRespImage(a_th);}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetCurrentCardImage()
      \param   a_th
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the index of the current cardiac image for reconstruction
    */
    inline int GetCurrentCardImage(int a_th)
           {return mp_DynamicDataManager->GetCurrentCardImage(a_th);}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetCurrentTimeFrame()
      \param   a_th
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the index of the current time frame
    */
    inline int GetCurrentTimeFrame(int a_th)
           {return mp_DynamicDataManager->GetCurrentTimeFrame(a_th);}
    /*!
      \fn      inline int oImageDimensionsAndQuantification::GetCurrentPMotionIndex()
      \param   a_th
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the index of the current involuntary motion gate
    */
    inline int GetCurrentPMotionIndex(int a_th)
           {return mp_DynamicDataManager->GetCurrentPMotionIndex(a_th);}

    /*!
      \fn      inline int oImageDimensionsAndQuantification::GetNbIPatMotionSubsets()
      \brief   call the eponym function from the oDynamicDataManager object
      \return  the total number of involuntary patient motion triggers
    */
    inline int GetNbIPatMotionSubsets()
           {return mp_DynamicDataManager->GetNbIPatMotionSubsets();}
           
           
  // -----------------------------------------------------------------------------------------
  // Get & Set functions
  public:
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the member m_verboseLevel to the provided value
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbThreadsForProjection()
      \brief   Get the number of threads used for projections
      \return  m_nbThreadsForProjection
    */
    inline int GetNbThreadsForProjection()
           {return m_nbThreadsForProjection;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbThreadsForImageComputation()
      \brief   Get the number of threads used for image operations
      \return  m_nbThreadsForImageComputation
    */
    inline int GetNbThreadsForImageComputation()
           {return m_nbThreadsForImageComputation;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbThreadsMax()
      \brief   Get the maximum between the number of threads used for projections and image operations
      \return  max(m_nbThreadsForProjection,m_nbThreadsForImageComputation)
    */
    inline int GetNbThreadsMax()
           { if (m_nbThreadsForProjection>m_nbThreadsForImageComputation) return m_nbThreadsForProjection;
             else return m_nbThreadsForImageComputation; }
    /*!
      \fn      public int oImageDimensionsAndQuantification::SetNbThreads()
      \param   const string& a_nbThreads
      \brief   Set the number of threads
      \details The string parameter can be a set of two parameters separated by a comma, the first for the number of threads for projections and
               the second for the number of threads for image operations. If no comma, then the number of threads is the same for all operations.
      \return  0 if correctly set, another value otherwise
    */
    int SetNbThreads(const string& a_nbThreads);
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetMPIRankAndSize()
      \param   int a_mpiRank
      \param   int a_mpiSize
      \brief   Set the MPI rank of the MPI instance, and the MPI size (the number of instances)
    */
    inline void SetMPIRankAndSize(int a_mpiRank, int a_mpiSize)
           {m_mpiRank = a_mpiRank; m_mpiSize = a_mpiSize;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetMPISize()
      \brief   Get the MPI size (the number of MPI instances)
      \return  m_mpiSize
    */
    inline int GetMPISize()
           {return m_mpiSize;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetMPIRank()
      \brief   Get the MPI instance number (the rank)
      \return  m_mpiRank
    */
    inline int GetMPIRank()
           {return m_mpiRank;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbBeds()
      \brief   Get the number of bed positions
      \return  m_nbBeds
    */
    inline int GetNbBeds()
           {return m_nbBeds;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbBeds()
      \param   int a_nbBeds
      \brief   Set the number of bed positions and allocate the bed positions if not already done.
    */
    inline void SetNbBeds(int a_nbBeds)
           {m_nbBeds = a_nbBeds; if (mp_bedPositions!=NULL) free(mp_bedPositions); mp_bedPositions = (FLTNB*)calloc(m_nbBeds,sizeof(FLTNB));}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetBedPosition()
      \param   int a_bedIndex
      \brief   Get the bed position associated to a bed index
      \return  The associated mp_bedPositions of the provided bed index. No check for overflow.
    */
    inline FLTNB GetBedPosition(int a_bedIndex)
           {return mp_bedPositions[a_bedIndex];}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetProvidedBedPositionFlag()
      \brief   Say if the bed relative position was provided from the datafile or not
      \return  m_providedBedPosition
    */
    inline bool GetProvidedBedPositionFlag() {return m_providedBedPosition;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbTimeFrames()
      \brief   Get the number of time frames
      \return  m_nbTimeFrames
    */
    inline int GetNbTimeFrames()
           {return m_nbTimeFrames;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbTimeBasisFunctions()
      \brief   Get the number of time basis functions
      \return  m_nbTimeBasisFunctions
    */
    inline int GetNbTimeBasisFunctions()
           {return m_nbTimeBasisFunctions;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetTimeBasisCoefficient()
      \param   int a_timeBasisFunction
      \param   int a_timeFrame
      \brief   Get the time basis coefficients for the provided frame and basis function
      \return  m2p_timeBasisFunctions[a_timeBasisFunction][a_timeFrame]
    */
    inline FLTNB GetTimeBasisCoefficient(int a_timeBasisFunction, int a_timeFrame)
           {return m2p_timeBasisFunctions[a_timeBasisFunction][a_timeFrame];}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetFrameDurationInSec()
      \param   int a_bed
      \param   int a_frame
      \brief   Get the frame duration for the given bed, in seconds as a FLTNB
      \return  ((FLTNB)(m2p_frameDurationsInMs[a_bed][a_frame]))/1000.
    */
    inline FLTNB GetFrameDurationInSec(int a_bed, int a_frame)
           {return ((FLTNB)(m2p_frameDurationsInMs[a_bed][a_frame]))/1000.;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetFrameTimeStartInSec()
      \param   int a_bed
      \param   int a_frame
      \brief   Get the frame time start for the given bed, in seconds as a FLTNB
      \return  ((FLTNB)(m2p_frameTimeStartInMs[a_bed][a_frame]))/1000.
    */
    inline FLTNB GetFrameTimeStartInSec(int a_bed, int a_frame)
           {return ((FLTNB)(m2p_frameTimeStartInMs[a_bed][a_frame]))/1000.;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetFrameTimeStopInSec()
      \param   int a_bed
      \param   int a_frame
      \brief   Get the frame time stop for the given bed, in seconds as a FLTNB
      \return  ((FLTNB)(m2p_frameTimeStopInMs[a_bed][a_frame]))/1000.
    */
    inline FLTNB GetFrameTimeStopInSec(int a_bed, int a_frame)
           {return ((FLTNB)(m2p_frameTimeStopInMs[a_bed][a_frame]))/1000.;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetFinalTimeStopInSec()
      \param   int a_bed
      \brief   Get the last frame time stop for the given bed, in seconds as a FLTNB
      \return  ((FLTNB)(m2p_frameTimeStopInMs[a_bed][m_nbTimeFrames-1]))/1000.
    */
    inline FLTNB GetFinalTimeStopInSec(int a_bed)
           {return ((FLTNB)(m2p_frameTimeStopInMs[a_bed][m_nbTimeFrames-1]))/1000.;}
    /*!
      \fn      public inline uint32_t oImageDimensionsAndQuantification::GetFrameDurationInMs()
      \param   int a_bed
      \param   int a_frame
      \brief   Get the frame duration for the given bed, in milliseconds as a uint32_t
      \return  m2p_frameDurationsInMs[a_bed][a_frame]
    */
    inline uint32_t GetFrameDurationInMs(int a_bed, int a_frame)
           {return m2p_frameDurationsInMs[a_bed][a_frame];}
    /*!
      \fn      public inline uint32_t oImageDimensionsAndQuantification::GetFrameTimeStartInMs()
      \param   int a_bed
      \param   int a_frame
      \brief   Get the frame time start for the given bed, in milliseconds as a uint32_t
      \return  m2p_frameTimeStartInMs[a_bed][a_frame]
    */
    inline uint32_t GetFrameTimeStartInMs(int a_bed, int a_frame)
           {return m2p_frameTimeStartInMs[a_bed][a_frame];}
    /*!
      \fn      public inline uint32_t oImageDimensionsAndQuantification::GetFrameTimeStopInMs()
      \param   int a_bed
      \param   int a_frame
      \brief   Get the frame time stop for the given bed, in milliseconds as a uint32_t
      \return  m2p_frameTimeStopInMs[a_bed][a_frame]
    */
    inline uint32_t GetFrameTimeStopInMs(int a_bed, int a_frame)
           {return m2p_frameTimeStopInMs[a_bed][a_frame];}
    /*!
      \fn      public inline uint32_t oImageDimensionsAndQuantification::GetFinalTimeStopInMs()
      \param   int a_bed
      \brief   Get the last frame time stop for the given bed, in milliseconds as a uint32_t
      \return  m2p_frameTimeStopInMs[a_bed][m_nbTimeFrames-1]
    */
    inline uint32_t GetFinalTimeStopInMs(int a_bed)
           {return m2p_frameTimeStopInMs[a_bed][m_nbTimeFrames-1];}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetFrames()
      \param   const string& a_frameList
      \brief   Set the frame list (a string that will be parsed by the InitializeFramingAndQuantification function)
    */
    inline void SetFrames(const string& a_frameList)
           {m_frameList = a_frameList;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbTimeBasisFunctions()
      \param   int a_nbTimeBasisFunctions
      \brief   Set the number of time basis functions
    */
    inline void SetNbTimeBasisFunctions(int a_nbTimeBasisFunctions)
           {m_nbTimeBasisFunctions = a_nbTimeBasisFunctions;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetTimeBasisFunctionsFile()
      \param   const string& a_timeBasisFunctionsFile
      \brief   Set the file name containing the time basis functions coefficients
    */
    inline void SetTimeBasisFunctionsFile(const string& a_timeBasisFunctionsFile)
           {m_timeBasisFunctionsFile = a_timeBasisFunctionsFile;}
    /*!
      \fn      public int oImageDimensionsAndQuantification::SetAcquisitionTime()
      \param   int a_bed
      \param   FLTNB a_timeStartInSec
      \param   FLTNB a_durationInSec
      \brief   Set the acquisition time if not already set by the SetTimeFrames()
      \details This function is called from the vDataFile once the acquisition start and duration are read in the header
      \return  0 if no problem, another value otherwise
    */
    int SetAcquisitionTime(int a_bed, FLTNB a_timeStartInSec, FLTNB a_durationInSec);
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetTimeStaticFlag()
      \brief   Get the time static flag that says if the reconstruction has only one frame or not
      \return  m_timeStaticFlag
    */
    inline bool GetTimeStaticFlag()
           {return m_timeStaticFlag;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbRespGates()
      \brief   Get the number of respiratory gates
      \return  m_nbRespGates
    */
    inline int GetNbRespGates()
           {return m_nbRespGates;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbRespBasisFunctions()
      \brief   Get the number of respiratory basis functions
      \return  m_nbRespBasisFunctions
    */
    inline int GetNbRespBasisFunctions()
           {return m_nbRespBasisFunctions;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetRespBasisCoefficient()
      \param   int a_respBasisFunction
      \param   int a_respGate
      \brief   Get the respiratory basis coefficients for the provided respiratory gate and basis function
      \return  m2p_respBasisFunctions[a_respBasisFunction][a_respGate]
    */
    inline FLTNB GetRespBasisCoefficient(int a_respBasisFunction, int a_respGate)
           {return m2p_respBasisFunctions[a_respBasisFunction][a_respGate];}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbRespGates()
      \param   int a_nbRespGates
      \brief   Set the number of respiratory gates
    */
    inline void SetNbRespGates(int a_nbRespGates)
           {m_nbRespGates = a_nbRespGates;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbRespBasisFunctions()
      \param   int a_nbRespBasisFunctions
      \brief   Set the number of respiratory basis functions
    */
    inline void SetNbRespBasisFunctions(int a_nbRespBasisFunctions)
           {m_nbRespBasisFunctions = a_nbRespBasisFunctions;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetRespBasisFunctionsFile()
      \param   const string& a_respBasisFunctionsFile
      \brief   Set the file name containing the respiratory basis functions coefficients
    */
    inline void SetRespBasisFunctionsFile(const string& a_respBasisFunctionsFile)
           {m_respBasisFunctionsFile = a_respBasisFunctionsFile;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetRespStaticFlag()
      \brief   Get the respiratory static flag that says if the reconstruction has only one respiratory gate or not
      \return  m_respStaticFlag
    */
    inline bool GetRespStaticFlag()
           {return m_respStaticFlag;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbCardGates()
      \brief   Get the number of cardiac gates
      \return  m_nbCardGates
    */
    inline int GetNbCardGates()
           {return m_nbCardGates;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbCardBasisFunctions()
      \brief   Get the number of cardiac basis functions
      \return  m_nbCardBasisFunctions
    */
    inline int GetNbCardBasisFunctions()
           {return m_nbCardBasisFunctions;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetCardBasisCoefficient()
      \param   int a_cardBasisFunction
      \param   int a_cardGate
      \brief   Get the cardiac basis coefficients for the provided cardiac gate and basis function
      \return  m2p_cardBasisFunctions[a_cardBasisFunction][a_cardGate]
    */
    inline FLTNB GetCardBasisCoefficient(int a_cardBasisFunction, int a_cardGate)
           {return m2p_cardBasisFunctions[a_cardBasisFunction][a_cardGate];}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbCardGates()
      \param   int a_nbCardGates
      \brief   Set the number of cardiac gates
    */
    inline void SetNbCardGates(int a_nbCardGates)
           {m_nbCardGates = a_nbCardGates;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbCardBasisFunctions()
      \param   int a_nbCardBasisFunctions
      \brief   Set the number of cardiac basis functions
    */
    inline void SetNbCardBasisFunctions(int a_nbCardBasisFunctions)
           {m_nbCardBasisFunctions = a_nbCardBasisFunctions;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetCardBasisFunctionsFile()
      \param   const string& a_cardBasisFunctionsFile
      \brief   Set the file name containing the cardiac basis functions coefficients
    */
    inline void SetCardBasisFunctionsFile(const string& a_cardBasisFunctionsFile)
           {m_cardBasisFunctionsFile = a_cardBasisFunctionsFile;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetCardStaticFlag()
      \brief   Get the cardiac static flag that says if the reconstruction has only one cardiac gate or not
      \return  m_cardStaticFlag
    */
    inline bool GetCardStaticFlag()
           {return m_cardStaticFlag;}
    /*!
      \fn      public inline INTNB oImageDimensionsAndQuantification::GetNbVoxX()
      \brief   Get the number of voxels along the X axis
      \return  m_nbVoxX
    */
    inline INTNB GetNbVoxX()
           {return m_nbVoxX;}
    /*!
      \fn      public inline INTNB oImageDimensionsAndQuantification::GetNbVoxY()
      \brief   Get the number of voxels along the Y axis
      \return  m_nbVoxY
    */
    inline INTNB GetNbVoxY()
           {return m_nbVoxY;}
    /*!
      \fn      public inline INTNB oImageDimensionsAndQuantification::GetNbVoxZ()
      \brief   Get the number of voxels along the Z axis
      \return  m_nbVoxZ
    */
    inline INTNB GetNbVoxZ()
           {return m_nbVoxZ;}
    /*!
      \fn      public inline INTNB oImageDimensionsAndQuantification::GetNbVoxXY()
      \brief   Get the number of voxels in a slice
      \return  m_nbVoxXY
    */
    inline INTNB GetNbVoxXY()
           {return m_nbVoxXY;}
    /*!
      \fn      public inline INTNB oImageDimensionsAndQuantification::GetNbVoxXYZ()
      \brief   Get the total number of voxels
      \return  m_nbVoxXYZ
    */
    inline INTNB GetNbVoxXYZ()
           {return m_nbVoxXYZ;}
    /*!
      \fn      public inline INTNB oImageDimensionsAndQuantification::GetNbVoxDiagonal()
      \brief   Get an estimation of the number of voxels along the image diagonal
      \return  An estimation of the number of voxels along the image diagonal
    */
    inline INTNB GetNbVoxDiagonal()
           { return ((INTNB)(sqrt( ((FLTNB)m_nbVoxX) * ((FLTNB)m_nbVoxX) +
                                   ((FLTNB)m_nbVoxY) * ((FLTNB)m_nbVoxY) +
                                   ((FLTNB)m_nbVoxZ) * ((FLTNB)m_nbVoxZ) ))); }
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbVoxX()
      \param   INTNB a_nbVoxX
      \brief   Set the number of voxels along the X axis
    */
    inline void SetNbVoxX(INTNB a_nbVoxX)
           {m_nbVoxX = a_nbVoxX;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbVoxY()
      \param   INTNB a_nbVoxY
      \brief   Set the number of voxels along the Y axis
    */
    inline void SetNbVoxY(INTNB a_nbVoxY)
           {m_nbVoxY = a_nbVoxY;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbVoxZ()
      \param   INTNB a_nbVoxZ
      \brief   Set the number of voxels along the Z axis
    */
    inline void SetNbVoxZ(INTNB a_nbVoxZ)
           {m_nbVoxZ = a_nbVoxZ;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetVoxSizeX()
      \brief   Get the voxel's size along the X axis, in mm
      \return  m_voxSizeX
    */
    inline FLTNB GetVoxSizeX()
           {return m_voxSizeX;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetVoxSizeY()
      \brief   Get the voxel's size along the Y axis, in mm
      \return  m_voxSizeY
    */
    inline FLTNB GetVoxSizeY()
           {return m_voxSizeY;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetVoxSizeZ()
      \brief   Get the voxel's size along the Z axis, in mm
      \return  m_voxSizeZ
    */
    inline FLTNB GetVoxSizeZ()
           {return m_voxSizeZ;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetFOVSizeX()
      \brief   Get the size of the field-of-view along the X axis, in mm
      \return  m_fovSizeX
    */
    inline FLTNB GetFOVSizeX()
           {return m_fovSizeX;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetFOVSizeY()
      \brief   Get the size of the field-of-view along the Y axis, in mm
      \return  m_fovSizeY
    */
    inline FLTNB GetFOVSizeY()
           {return m_fovSizeY;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetFOVSizeZ()
      \brief   Get the size of the field-of-view along the Z axis, in mm
      \return  m_fovSizeZ
    */
    inline FLTNB GetFOVSizeZ()
           {return m_fovSizeZ;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetVoxSizeX()
      \param   FLTNB a_voxSizeX
      \brief   Set the voxel's size along the X axis, in mm
    */
    inline void SetVoxSizeX(FLTNB a_voxSizeX)
           {m_voxSizeX = a_voxSizeX;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetVoxSizeY()
      \param   FLTNB a_voxSizeY
      \brief   Set the voxel's size along the Y axis, in mm
    */
    inline void SetVoxSizeY(FLTNB a_voxSizeY)
           {m_voxSizeY = a_voxSizeY;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetVoxSizeZ()
      \param   FLTNB a_voxSizeZ
      \brief   Set the voxel's size along the Z axis, in mm
    */
    inline void SetVoxSizeZ(FLTNB a_voxSizeZ)
           {m_voxSizeZ = a_voxSizeZ;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetFOVSizeX()
      \param   FLTNB a_fovSizeX
      \brief   Set the FOV's size along the X axis, in mm
    */
    inline void SetFOVSizeX(FLTNB a_fovSizeX)
           {m_fovSizeX = a_fovSizeX;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetFOVSizeY()
      \param   FLTNB a_fovSizeY
      \brief   Set the FOV's size along the Y axis, in mm
    */
    inline void SetFOVSizeY(FLTNB a_fovSizeY)
           {m_fovSizeY = a_fovSizeY;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetFOVSizeZ()
      \param   FLTNB a_fovSizeZ
      \brief   Set the FOV's size along the Z axis, in mm
    */
    inline void SetFOVSizeZ(FLTNB a_fovSizeZ)
           {m_fovSizeZ = a_fovSizeZ;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetFOVOutMasking()
      \param   FLTNB a_fovOutPercent
      \param   INTNB a_nbSliceOutMask
      \brief   Set the output FOV masking settings: transaxial unmasked FOV percent and number of extrem slices to removed
    */
    inline void SetFOVOutMasking(FLTNB a_fovOutPercent, INTNB a_nbSliceOutMask)
           {m_fovOutPercent = a_fovOutPercent; m_nbSliceOutMask = a_nbSliceOutMask;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetFOVOutPercent()
      \brief   Get the percentage of transaxial unmasked FOV
      \return  m_fovOutPercent
    */
    inline FLTNB GetFOVOutPercent()
           {return m_fovOutPercent;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetNbSliceOutMask()
      \brief   Get the number of extrem slices that will be masked at each side
      \return  m_nbSliceOutMask
    */
    inline INTNB GetNbSliceOutMask()
           {return m_nbSliceOutMask;}
    /*!
      \fn      public int oImageDimensionsAndQuantification::SetFlipOut()
      \param   const string& a_flipOut
      \brief   Set the output flip options, the parameter being a string potentially containing the letters x, y and z
      \return  0 if no problem, another value if the string is not correctly set
    */
    int SetFlipOut(const string& a_flipOut);
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetFlipOutX()
      \brief   Get the boolean saying if the output image must be flipped along the X axis
      \return  m_flipOutX
    */
    inline bool GetFlipOutX()
           {return m_flipOutX;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetFlipOutY()
      \brief   Get the boolean saying if the output image must be flipped along the Y axis
      \return  m_flipOutY
    */
    inline bool GetFlipOutY()
           {return m_flipOutY;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetFlipOutZ()
      \brief   Get the boolean saying if the output image must be flipped along the Z axis
      \return  m_flipOutZ
    */
    inline bool GetFlipOutZ()
           {return m_flipOutZ;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetOffsetX()
      \brief   Get the image offset along the X axis, in mm
      \return  m_offsetX
    */
    inline FLTNB GetOffsetX()
           {return m_offsetX;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetOffsetY()
      \brief   Get the image offset along the Y axis, in mm
      \return  m_offsetY
    */
    inline FLTNB GetOffsetY()
           {return m_offsetY;}
    /*!
      \fn      public inline FLTNB oImageDimensionsAndQuantification::GetOffsetZ()
      \brief   Get the image offset along the Z axis, in mm
      \return  m_offsetZ
    */
    inline FLTNB GetOffsetZ()
           {return m_offsetZ;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetOffsetX()
      \param   FLTNB a_offsetX
      \brief   Set the image offset along the X axis, in mm
    */
    inline void SetOffsetX(FLTNB a_offsetX)
           {m_offsetX = a_offsetX;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetOffsetY()
      \param   FLTNB a_offsetY
      \brief   Set the image offset along the Y axis, in mm
    */
    inline void SetOffsetY(FLTNB a_offsetY)
           {m_offsetY = a_offsetY;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetOffsetZ()
      \param   FLTNB a_offsetZ
      \brief   Set the image offset along the Z axis, in mm
    */
    inline void SetOffsetZ(FLTNB a_offsetZ)
           {m_offsetZ = a_offsetZ;}
    /*!
      \fn      public FLTNB oImageDimensionsAndQuantification::GetQuantificationFactor()
      \param   int a_bed
      \param   int a_frame
      \param   int a_respGate
      \param   int a_cardGate
      \brief   Get the quantification factor corresponding to the provided bed, frame, respiratory and cardiac gates
    */
    FLTNB GetQuantificationFactor(int a_bed, int a_frame, int a_respGate, int a_cardGate);
    /*!
      \fn      public int oImageDimensionsAndQuantification::SetCalibrationFactor()
      \param   int a_bed
      \param   FLTNB a_calibrationFactor
      \brief   Set the calibration factor for the provided bed
      \details Apply it to all frames and gates
      \return  0 if no problem, another value otherwise
    */
    int SetCalibrationFactor(int a_bed, FLTNB a_calibrationFactor);
    /*!
      \fn      public int oImageDimensionsAndQuantification::SetPETIsotope()
      \param   int a_bed
      \param   const string& a_isotope
      \brief   Set the PET isotope for the provided bed
      \details Apply decay and branching ratio corrections if the PET isotope is found in the database of the config/ directory
      \return  0 if no problem, another value otherwise
    */
    int SetPETIsotope(int a_bed, const string& a_isotope);
    /*!
      \fn      public int oImageDimensionsAndQuantification::SetSPECTIsotope()
      \param   int a_bed
      \param   const string& a_isotope
      \brief   Set the SPECT isotope for the provided bed
      \details Apply decay and branching ratio corrections if the SPECT isotope is found in the database of the config/ directory
      \return  0 if no problem, another value otherwise
    */
    int SetSPECTIsotope(int a_bed, const string& a_isotope);
    /*!
      \fn      public int oImageDimensionsAndQuantification::SetDynamicSpecificQuantificationFactors()
      \param   const string& a_quantificationFile
      \brief   Apply specific quantification factors manually provided as an option
      \details It does not set the quantification factors, it applies them (it is an update, not an affectation)
      \return  0 if no problem, another value otherwise
    */
    int SetDynamicSpecificQuantificationFactors(const string& a_quantificationFile);
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetIgnoredCorrections()
      \param   const string& a_ignoredCorrectionsList
      \brief   Set the string specifying the corrections that will be ignored
    */
    inline void SetIgnoredCorrections(const string& a_ignoredCorrectionsList)
           {m_ignoredCorrectionsList = a_ignoredCorrectionsList;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetIgnoreAttnCorrectionFlag()
      \brief   Get the boolean that says if the attenuation correction is ignored or not
      \return  m_ignoreAttnCorrectionFlag
    */
    inline bool GetIgnoreAttnCorrectionFlag()
           {return m_ignoreAttnCorrectionFlag;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetIgnoreNormCorrectionFlag()
      \brief   Get the boolean that says if the normalization correction is ignored or not
      \return  m_ignoreNormCorrectionFlag
    */
    inline bool GetIgnoreNormCorrectionFlag()
           {return m_ignoreNormCorrectionFlag;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetIgnoreScatCorrectionFlag()
      \brief   Get the boolean that says if the scatter correction is ignored or not
      \return  m_ignoreScatCorrectionFlag
    */
    inline bool GetIgnoreScatCorrectionFlag()
           {return m_ignoreScatCorrectionFlag;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetIgnoreRandCorrectionFlag()
      \brief   Get the boolean that says if the random correction is ignored or not
      \return  m_ignoreRandCorrectionFlag
    */
    inline bool GetIgnoreRandCorrectionFlag()
           {return m_ignoreRandCorrectionFlag;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetIgnoreDecaCorrectionFlag()
      \brief   Get the boolean that says if the decay correction is ignored or not
      \return  m_ignoreDecaCorrectionFlag
    */
    inline bool GetIgnoreDecaCorrectionFlag()
           {return m_ignoreDecaCorrectionFlag;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetIgnoreBratCorrectionFlag()
      \brief   Get the boolean that says if the branching ratio correction is ignored or not
      \return  m_ignoreBratCorrectionFlag
    */
    inline bool GetIgnoreBratCorrectionFlag()
           {return m_ignoreBratCorrectionFlag;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetIgnoreFdurCorrectionFlag()
      \brief   Get the boolean that says if the frame duration correction is ignored or not
      \return  m_ignoreFdurCorrectionFlag
    */
    inline bool GetIgnoreFdurCorrectionFlag()
           {return m_ignoreFdurCorrectionFlag;}
    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::GetIgnoreCaliCorrectionFlag()
      \brief   Get the boolean that says if the calibration correction is ignored or not
      \return  m_ignoreCaliCorrectionFlag
    */
    inline bool GetIgnoreCaliCorrectionFlag()
           {return m_ignoreCaliCorrectionFlag;}
    /*!
      \fn      public inline int oImageDimensionsAndQuantification::GetNbMultiModalImages()
      \brief   Get the number of additional multimodal images
      \return  m_nbMultiModalImages
    */
    inline int GetNbMultiModalImages()
           {return m_nbMultiModalImages;}
    /*!
      \fn      public inline void oImageDimensionsAndQuantification::SetNbMultiModalImages()
      \param   int a_nbMultiModalImages
      \brief   Set the number of additional multimodal images.
    */
    inline void SetNbMultiModalImages(int a_nbMultiModalImages)
           {m_nbMultiModalImages = a_nbMultiModalImages;}

    /*!
      \fn      public inline bool oImageDimensionsAndQuantification::IsInitialized()
      \brief   Returns true if the object has been initialized
    */
    inline bool IsInitialized() {return m_initialized;}
           

  // -----------------------------------------------------------------------------------------
  // Data members
  private:
    oDynamicDataManager* mp_DynamicDataManager; /*!< oDynamicDataManager object related to dynamic raw data management*/
    // Number of threads
    int m_nbThreadsForProjection;       /*!< The number of threads for projections */
    int m_nbThreadsForImageComputation; /*!< The number of threads for image operations */
    // MPI size and rank
    int m_mpiRank;                      /*!< The rank of the MPI instance */
    int m_mpiSize;                      /*!< The size of the MPI process (i.e. the number of instances) */
    // Number of beds
    int m_nbBeds;                       /*!< The number of bed positions simultaneously reconstructed */
    FLTNB* mp_bedPositions;             /*!< The bed positions used during the reconstruction */
    bool m_providedBedPosition;         /*!< A flag saying if the bed position has been provided from the datafile */
    // Framing of the acquisition
    string m_frameList;                 /*!< A string containing the list of reconstructed frames */
    int m_nbTimeBasisFunctions;         /*!< The number of time basis functions */
    int m_nbTimeFrames;                 /*!< The number of time frames */
    FLTNB** m2p_timeBasisFunctions;     /*!< The table of time basis functions coefficients */
    uint32_t** m2p_frameDurationsInMs;  /*!< The table of frame durations, per bed position */
    uint32_t** m2p_frameTimeStartInMs;  /*!< The table of frame time start, per bed position */
    uint32_t** m2p_frameTimeStopInMs;   /*!< The table of frame time stop, per bed position */
    string m_timeBasisFunctionsFile;    /*!< The file containing the time basis coefficients */
    bool m_timeStaticFlag;              /*!< A boolean saying if only one frame is being reconstructed */
    // Quantification
    FLTNB*** m3p_quantificationFactors; /*!< The table of quantification factors, per bed, per frame, per gate */
    string m_ignoredCorrectionsList;    /*!< The string containing the list of corrections to be ignored */
    bool m_ignoreAttnCorrectionFlag;    /*!< A boolean saying if the attenuation correction is ignored */
    bool m_ignoreNormCorrectionFlag;    /*!< A boolean saying if the normalization correction is ignored */
    bool m_ignoreRandCorrectionFlag;    /*!< A boolean saying if the random correction is ignored */
    bool m_ignoreScatCorrectionFlag;    /*!< A boolean saying if the scatter correction is ignored */
    bool m_ignoreDecaCorrectionFlag;    /*!< A boolean saying if the decay correction is ignored */
    bool m_ignoreBratCorrectionFlag;    /*!< A boolean saying if the branching ratio correction is ignored */
    bool m_ignoreFdurCorrectionFlag;    /*!< A boolean saying if the frame duration correction is ignored */
    bool m_ignoreCaliCorrectionFlag;    /*!< A boolean saying if the calibration correction is ignored */
    // Respiratory gating
    int m_nbRespBasisFunctions;         /*!< The number of respiratory basis functions */
    int m_nbRespGates;                  /*!< The number of respiratory gates */
    FLTNB** m2p_respBasisFunctions;     /*!< The table of respiratory basis functions coefficients */
    string m_respBasisFunctionsFile;    /*!< The file containing the respiratory basis coefficients */
    bool m_respStaticFlag;              /*!< A boolean saying if only one respiratory gate is being reconstructed */
    // Cardiac gating
    int m_nbCardBasisFunctions;         /*!< The number of cardiac basis functions */
    int m_nbCardGates;                  /*!< The number of cardiac gates */
    FLTNB** m2p_cardBasisFunctions;     /*!< The table of cardiac basis functions coefficients */
    string m_cardBasisFunctionsFile;    /*!< The file containing the cardiac basis coefficients */
    bool m_cardStaticFlag;              /*!< A boolean saying if only one cardiac gate is being reconstructed */
    // Image discretization
    INTNB m_nbVoxX;                    /*!< The number of voxels along the X axis */
    INTNB m_nbVoxY;                    /*!< The number of voxels along the Y axis */
    INTNB m_nbVoxZ;                    /*!< The number of voxels along the Z axis */
    INTNB m_nbVoxXY;                   /*!< The number of voxels in a slice */
    INTNB m_nbVoxXYZ;                  /*!< The total number of voxels */
    FLTNB m_voxSizeX;                  /*!< The voxel's size along the X axis, in mm */
    FLTNB m_voxSizeY;                  /*!< The voxel's size along the Y axis, in mm */
    FLTNB m_voxSizeZ;                  /*!< The voxel's size along the Z axis, in mm */
    FLTNB m_fovSizeX;                  /*!< The size of the field-of-view along the X axis, in mm */
    FLTNB m_fovSizeY;                  /*!< The size of the field-of-view along the Y axis, in mm */
    FLTNB m_fovSizeZ;                  /*!< The size of the field-of-view along the Z axis, in mm */
    FLTNB m_offsetX;                   /*!< The image offset along the X axis, in mm (default: 0.) */
    FLTNB m_offsetY;                   /*!< The image offset along the Y axis, in mm (default: 0.) */
    FLTNB m_offsetZ;                   /*!< The image offset along the Z axis, in mm (default: 0.) */
    // Percentage of transaxial FOV unmasked before saving
    FLTNB m_fovOutPercent;             /*!< The percentage of the transaxial field-of-view not masked before saving reconstructed images (default: 0. = none) */
    INTNB m_nbSliceOutMask;            /*!< The number of extrem slices to be masked from both sides of the axial FOV before saving reconstruction images (default: 0) */
    // Output flip
    bool m_flipOutX;                   /*!< A boolean saying if the images must be flipped along the X axis before being saved (default: false) */
    bool m_flipOutY;                   /*!< A boolean saying if the images must be flipped along the Y axis before being saved (default: false) */
    bool m_flipOutZ;                   /*!< A boolean saying if the images must be flipped along the Z axis before being saved (default: false) */
    // Parameters checked
    bool m_checked;                    /*!< A boolean that says if the function CheckParameters() has been called */
    // Object initialized
    bool m_initialized;                /*!< A boolean that says if the function Initialize() has been called */
    // Verbose
    int m_verbose;                     /*!< The verbose level */
    int m_nbMultiModalImages;          /*!< The number of additional multimodal images */
};

#endif
