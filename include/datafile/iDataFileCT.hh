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
  \brief    Declaration of class iDataFileCT
*/

#ifndef IDATAFILECT_HH
#define IDATAFILECT_HH 1

#include "gVariables.hh"
#include "vDataFile.hh"
#include "vScanner.hh"
#include "iEventHistoCT.hh"
#include "iEventListCT.hh"


/*!
  \class   iDataFileCT
  \brief   Inherit from vDataFile. Class that manages the reading of a CT input file (header + data).
  \details It contains several arrays corresponding to the different kind of informations the data file could contain. \n
           As many booleans as arrays say if the data are here or not. The data file can be either completely loaded, or read event by event during reconstruction. \n
           MPI is coming here to cut the data file into peaces (also either can be loaded or read on-the-fly).
*/
class iDataFileCT : public vDataFile
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iDataFileCT constructor. 
               Initialize the member variables to their default values.
    */
    iDataFileCT();
    /*!
      \brief   iDataFileCT destructor. 
    */
    ~iDataFileCT();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      iDataFileCT::ReadSpecificInfoInHeader()
      \param   bool a_affectQuantificationFlag
      \brief   Read through the header file and recover specific CT information.
      \details If the parameter flag is on, then affect the quantification factors from the oImageDimensionsAndQuantification after
               reading relevant information
      \return  0 is success, positive value otherwise
    */
    int ReadSpecificInfoInHeader(bool a_affectQuantificationFlag);
    /*!
      \fn      iDataFileCT::WriteHeader()
      \brief   Generate a header file according to the data output information.
      \return  0 if success, and positive value otherwise.
    */
    int WriteHeader();
    /*!
      \fn      iDataFileCT::ComputeSizeEvent()
      \brief   Computation of the size of each event according to the mandatory/optional correction fields
      \return  0 is success, positive value otherwise
    */
    int ComputeSizeEvent();
    /*!
      \fn      iDataFileCT::PrepareDataFile()
      \brief   Store different kind of information inside arrays (data relative to specific correction as well as basic raw data for the case data is loaded in RAM) \n
               Use the flag provided by the user to determine how the data has to be sorted (preloaded or read on the fly)
      \return  0 is success, positive value otherwise
    */
    int PrepareDataFile();
    /*!
      \fn      iDataFileCT::WriteEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write event according to the chosen type of data
      \return  0 if success, and positive value otherwise.
    */
    int WriteEvent(vEvent* ap_Event, int a_th=0);
    /*!
      \fn      iDataFileCT::GetEventSpecific()
      \param   ap_buffer : address pointing to the event to recover
      \param   a_th : index of the thread from which the function was called
      \brief   Read an event from the position pointed by 'ap_buffer', parse the generic or modality-specific information, and store them in the (multithreaded) 'm2p_BufferEvent' object
      \return  the thread-specific  'm2p_BufferEvent' object containing the modality-specific information for the event
    */
    vEvent* GetEventSpecific(char* ap_buffer, int a_th);
    /*!
      \fn      iDataFileCT::InitAngles()
      \param   ap_angles
      \brief   allocate memory for the mp_angles variable using m_nbProjections
               and initialize the projection angles with the provided list of values
      \return  0 if success, positive value otherwise
    */
    int InitAngles(FLTNB* ap_angles);
    /*!
      \fn      iDataFileCT::DescribeSpecific()
      \brief   Implementation of the pure virtual eponym function that simply prints info about the datafile
    */
    void DescribeSpecific();

  // -------------------------------------------------------------------
  // Public functions dedicated to the projection script
  public:
    /*!
      \fn      iDataFileCT::PROJ_InitFile()
      \brief   Initialize the fstream objets for output writing as well as some other variables specific to the Projection script (Event-based correction flags, Estimated size of data file)
      \return  0 if success, and positive value otherwise.
    */
    int PROJ_InitFile();
    /*!
      \fn      iDataFileCT::PROJ_GetScannerSpecificParameters()
      \brief   Get SPECT specific parameters for projections from the scanner object, through the scannerManager.
      \return  0 if success, positive value otherwise
    */
    int PROJ_GetScannerSpecificParameters();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      iDataFileCT::GetNbProjections() 
      \return  total number of projections in the CT acquisition
    */
    inline uint16_t GetNbProjections()
           {return m_nbOfProjections;};
    /*!
      \fn      iDataFileCT::GetAngles() 
      \return  Angles for each projection of the acquisition
    */
    inline FLTNB* GetAngles()
           {return mp_angles;};
    /*!
      \fn      iDataFileCT::SetEventKindFlagOn()
      \brief   set to true the flag indicating the presence of the kind of a list-mode event in the datafile
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetEventKindFlagOn()
           {m_eventKindFlag = true;}
    /*!
      \fn      iDataFileCT::SetNbProjections()
      \param   a_nbProjections
      \brief   initialize the number of projections
    */
    inline void SetNbProjections(uint16_t a_nbProjections)
           {m_nbOfProjections = a_nbProjections;} 
    /*!
      \fn      iDataFileCT::GetDetectorRotDirection()
      \brief   Simply return m_detectorRotDirection
      \return  m_detectorRotDirection
    */
    inline int GetDetectorRotDirection()
           {return m_detectorRotDirection;}
    /*!
      \fn      iDataFileCT::SetDetectorRotDirection()
      \param   a_direction
      \brief   initialize the rotation direction of the CT detector
    */
    inline void SetDetectorRotDirection(int a_direction)
           {m_detectorRotDirection = a_direction;}
    /*!
      \fn      iDataFileCT::GetEventKindFlag()
      \brief   Simply return m_eventKindFlag
      \return  m_eventKindFlag
    */
    inline bool GetEventKindFlag()
           {return m_eventKindFlag;}
    /*!
      \fn      iDataFileCT::GetScatCorrectionFlag()
      \brief   Simply return m_scatCorrectionFlag
      \return  m_scatCorrectionFlag
    */
    inline bool GetScatCorrectionFlag()
           {return m_scatCorrectionFlag;}
    /*!
      \fn      iDataFileCT::GetBlankCorrectionFlag()
      \brief   Simply return m_blankCorrectionFlag
      \return  m_blankCorrectionFlag
    */
    inline bool GetBlankCorrectionFlag()
           {return m_blankCorrectionFlag;}

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      int iDataFileCT::SetSpecificParametersFrom()
      \brief   Initialize all parameters specific to CT from the provided datafile.
      \return  0 if success, and positive value otherwise
    */
    int SetSpecificParametersFrom(vDataFile* ap_DataFile);
    /*!
      \fn      iDataFileCT::CheckSpecificParameters()
      \brief   Check parameters specific to CT data
      \return  0 if success, and positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iDataFileCT::CheckFileSizeConsistency()
      \brief   This function is implemented in child classes \n
               Check if file size is consistent.
      \return  0 if success, and positive value otherwise.
    */
    int CheckFileSizeConsistency();
    /*!
      \fn      iDataFileCT::WriteHistoEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write a CT histogram event
      \return  0 if success, and positive value otherwise.
    */
    int WriteHistoEvent(iEventHistoCT* ap_Event, int a_th);
    /*!
      \fn      iDataFileCT::WriteListEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write a CT list-mode event
      \return  0 if success, and positive value otherwise.
    */
    int WriteListEvent(iEventListCT* ap_Event, int a_th);
    /*!
      \fn      iDataFileCT::CheckSpecificConsistencyWithAnotherDataFile()
      \param   vDataFile* ap_DataFile
      \brief   Check consistency between 'this' and the provided datafile, for specific characteristics.
      \details Implementation of the pure virtual function from vDataFile. It checks correction flags, etc.
      \return  0 if the provided datafile is consistent with 'this', another value otherwise
    */
    int CheckSpecificConsistencyWithAnotherDataFile(vDataFile* ap_DataFile);

  // -------------------------------------------------------------------
  // Data members
  private:
    bool m_eventKindFlag;             /*!< Flag for informations about the event nature (true, scatter) in the data. Default value = false */
    bool m_blankCorrectionFlag;       /*!< Flag that says if normalization correction terms are included in the data. Default = false */
    bool m_ignoreBlankCorrectionFlag; /*!< Flag to say if we ignore the normalization correction even if present. Default = false */
    bool m_scatCorrectionFlag;        /*!< Flag that says if scatter correction terms are included in the data. Default = false */
    bool m_ignoreScatCorrectionFlag;  /*!< Flag to say if we ignore the scatter correction even if present. Default = false */
    uint16_t m_nbOfProjections;       /*!< Total number of projections during the acquisition(for all the heads). No Default*/
    FLTNB* mp_angles;                 /*!< Angle [for each projection] in degrees. No Default */
    int m_detectorRotDirection;       /*!< Head rotation direction (0=clockwise, 1=counterclockwise)*/
};

#endif
