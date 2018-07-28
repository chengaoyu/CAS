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
  \brief    Declaration of class iDataFilePET
*/

#ifndef IDATAFILEPET_HH
#define IDATAFILEPET_HH 1

#include "gVariables.hh"
#include "vDataFile.hh"

/*!
  \class   iDataFilePET
  \brief   Inherit from vDataFile. Class that manages the reading of a PET input file (header + data).
  \details It contains several arrays corresponding to the different kind of informations the data file could contain. \n
           As many booleans as arrays say if the data are here or not. The data file can be either completely loaded, or read event by event during reconstruction. \n
           MPI is coming here to cut the data file into peaces (also either can be loaded or read on-the-fly).
*/
class iDataFilePET : public vDataFile
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iDataFilePET constructor. 
               Initialize the member variables to their default values.
    */
    iDataFilePET();
    /*!
      \brief   iDataFilePET destructor. 
    */
    ~iDataFilePET();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      iDataFilePET::ReadSpecificInfoInHeader()
      \param   bool a_affectQuantificationFlag
      \brief   Read through the header file and gather specific PET information.
      \details If the parameter flag is on, then affect the quantification factors from the oImageDimensionsAndQuantification after
               reading relevant information
      \return  0 is success, positive value otherwise
    */
    int ReadSpecificInfoInHeader(bool a_affectQuantificationFlag);
    /*!
      \fn      iDataFilePET::WriteHeader()
      \brief   Generate a header file according to the data output information.
      \return  0 if success, and positive value otherwise.
    */
    int WriteHeader();
    /*!
      \fn      iDataFilePET::ComputeSizeEvent()
      \brief   Computation of the size of each event according to the mandatory/optional correction fields
      \return  0 is success, positive value otherwise
    */
    int ComputeSizeEvent();
    /*!
      \fn      iDataFilePET::PrepareDataFile()
      \brief   Store different kind of information inside arrays (data relative to specific correction as well as basic raw data for the case data is loaded in RAM) \n
               Use the flag provided by the user to determine how the data has to be sorted (preloaded or read on the fly)
      \return  0 is success, positive value otherwise
    */
    int PrepareDataFile();
    /*!
      \fn      iDataFilePET::WriteEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write event according to the chosen type of data
      \return  0 if success, and positive value otherwise.
    */
    int WriteEvent(vEvent* ap_Event, int a_th);
    /*!
      \fn      iDataFilePET::GetEventSpecific()
      \param   ap_buffer : address pointing to the event to recover
      \param   a_th : index of the thread from which the function was called
      \brief   Read an event from the position pointed by 'ap_buffer', parse the generic or modality-specific information, and store them in the (multithreaded) 'm2p_BufferEvent' object
      \return  the thread-specific  'm2p_BufferEvent' object containing the modality-specific information for the event
    */
    vEvent* GetEventSpecific(char* ap_buffer, int a_th);
    /*!
      \fn      iDataFilePET::DescribeSpecific()
      \brief   Implementation of the pure virtual eponym function that simply prints info about the datafile
    */
    void DescribeSpecific();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      iDataFilePET::GetTOFInfoFlag()
      \return  m_TOFInfoFlag
    */
    inline bool GetTOFInfoFlag()
           {return m_TOFInfoFlag;}
    /*!
      \fn      iDataFilePET::GetIgnoreTOFFlag()
      \return  m_ignoreTOFFlag
    */
    inline bool GetIgnoreTOFFlag()
           {return m_ignoreTOFFlag;}
    /*!
      \fn      iDataFilePET::GetTOFResolution() 
      \return  TOF resolution for the acquisition
    */
    inline FLTNB GetTOFResolution()
           {return m_resolutionTOF;}
    /*!
      \fn      iDataFilePET::GetNbTOFBins() 
      \return  number of TOF bins in the acquisition
    */
    inline int GetNbTOFBins()
           {return m_nbBinsTOF;}
    /*!
      \fn      iDataFilePET::GetTOFBinSize() 
      \return  size of TOF bins in the acquisition
    */
    inline FLTNB GetTOFBinSize()
           {return m_binSizeTOF;}
    /*!
      \fn      iDataFilePET::GetTOFMeasurementRange() 
      \return  Maximum range of values for TOF delta measurements
    */
    inline FLTNB GetTOFMeasurementRange()
           {return m_TOFMeasurementRange;}
    /*!
      \fn      iDataFilePET::GetMaxAxialDiffmm()
      \return  max ring difference in the acquisition
    */
    inline FLTNB GetMaxAxialDiffmm()
           {return m_maxAxialDiffmm;}
    /*!
      \fn      iDataFilePET::SetMaxNumberOfLinesPerEvent()
      \brief   set the max number of line per event in the datafile
    */
    inline void SetMaxNumberOfLinesPerEvent(uint16_t a_value)
           {m_maxNumberOfLinesPerEvent = a_value;}
    /*!
      \fn      iDataFilePET::GetMaxNumberOfLinesPerEvent()
      \return  the max number of line per event in the datafile
    */
    inline uint16_t GetMaxNumberOfLinesPerEvent()
           {return m_maxNumberOfLinesPerEvent;}
    /*!
      \fn      iDataFilePET::SetIsotope()
      \param   a_value
      \brief   initialize the isotope string value
      \details The name should corresponds to one corresponding tag in the isotope configuration file in config/.
               This function is dedicated to datafile conversion scripts
    */
    inline void SetIsotope(string a_value)
           {m_isotope = a_value;}
    /*!
      \fn      iDataFilePET::GetIsotope()
      \return  the isotope string value
    */
    inline string GetIsotope()
           {return m_isotope;}
    /*!
      \fn      iDataFilePET::SetIgnoreTOFFlag()
      \param   a_flag
      \brief   Set a boolean that that if we ignore TOF information or not
    */
    inline void SetIgnoreTOFFlag(bool a_ignoreTOFFlag)
           {m_ignoreTOFFlag = a_ignoreTOFFlag;}
    /*!
      \fn      iDataFilePET::SetEventKindFlagOn()
      \brief   set to true the flag indicating the presence of the kind of a list-mode event in the datafile
      //TODO check if consistent with datafile type
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetEventKindFlagOn()
           {m_eventKindFlag = true;} 
    /*!
      \fn      void iDataFilePET::SetAtnCorrectionFlagOn()
      \brief   set to true the flag indicating the presence of attenuation correction factors in the datafile
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetAtnCorrectionFlagOn()
           {m_atnCorrectionFlag = true;} 
    /*!
      \fn      iDataFilePET::SetNormCorrectionFlagOn()
      \brief   set to true the flag indicating the presence of normalization correction factors in the datafile
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetNormCorrectionFlagOn()
           {m_normCorrectionFlag = true;} 
    /*!
      \fn      iDataFilePET::SetScatterCorrectionFlagOn()
      \brief   set to true the flag indicating the presence of scatter correction factors in the datafile
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetScatterCorrectionFlagOn()
           {m_scatCorrectionFlag = true;} 
    /*!
      \fn      iDataFilePET::SetRandomCorrectionFlagOn()
      \brief   set to true the flag indicating the presence of random correction factors in the datafile
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetRandomCorrectionFlagOn()
           {m_randCorrectionFlag = true;} 
    /*!
      \fn      iDataFilePET::GetNormCorrectionFlag()
      \brief   Simply return m_normCorrectionFlag
      \return  m_normCorrectionFlag
    */
    inline bool GetNormCorrectionFlag()
           {return m_normCorrectionFlag;}
    /*!
      \fn      iDataFilePET::GetAtnCorrectionFlag()
      \brief   Simply return m_atnCorrectionFlag
      \return  m_atnCorrectionFlag
    */
    inline bool GetAtnCorrectionFlag()
           {return m_atnCorrectionFlag;}
    /*!
      \fn      iDataFilePET::GetEventKindFlag()
      \brief   Simply return m_eventKindFlag
      \return  m_eventKindFlag
    */
    inline bool GetEventKindFlag()
           {return m_eventKindFlag;}
    /*!
      \fn      iDataFilePET::GetScatCorrectionFlag()
      \brief   Simply return m_scatCorrectionFlag
      \return  m_scatCorrectionFlag
    */
    inline bool GetScatCorrectionFlag()
           {return m_scatCorrectionFlag;}
    /*!
      \fn      iDataFilePET::GetRandCorrectionFlag()
      \brief   Simply return m_randCorrectionFlag
      \return  m_randCorrectionFlag
    */
    inline bool GetRandCorrectionFlag()
           {return m_randCorrectionFlag;}

  // -------------------------------------------------------------------
  // Public functions dedicated to the projection script
  public:
    /*!
      \fn      iDataFilePET::PROJ_InitFile()
      \brief   Initialize the fstream objets for output writing as well as some other variables specific to the Projection script (Event-based correction flags, Estimated size of data file)
      \return  0 if success, and positive value otherwise.
    */
    int PROJ_InitFile();
    /*!
      \fn      iDataFilePET::PROJ_GetScannerSpecificParameters()
      \brief   Get PET specific parameters for projections from the scanner object, through the scannerManager.
      \return  0 if success, positive value otherwise
    */
    int PROJ_GetScannerSpecificParameters();

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      int iDataFilePET::SetSpecificParametersFrom()
      \brief   Initialize all parameters specific to PET from the provided datafile.
      \return  0 if success, and positive value otherwise
    */
    int SetSpecificParametersFrom(vDataFile* ap_DataFile);
    /*!
      \fn      iDataFilePET::CheckSpecificParameters()
      \brief   Check parameters specific to PET data
      \return  0 if success, and positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iDataFilePET::WriteListEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write a PET list-mode event
      \return  0 if success, and positive value otherwise.
    */
    int WriteListEvent(iEventListPET* ap_Event, int a_th=0);
    /*!
      \fn      iDataFilePET::WriteHistoEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write a PET histogram event
      \return  0 if success, and positive value otherwise.
    */
    int WriteHistoEvent(iEventHistoPET* ap_Event, int a_th);
    /*!
      \fn      iDataFilePET::CheckFileSizeConsistency()
      \brief   This function is implemented in child classes \n
               Check if file size is consistent.
      \return  0 if success, and positive value otherwise.
    */
    int CheckFileSizeConsistency();
    /*!
      \fn      iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile()
      \param   vDataFile* ap_DataFile
      \brief   Check consistency between 'this' and the provided datafile, for specific characteristics.
      \details Implementation of the pure virtual function from vDataFile. It checks correction flags, etc.
      \return  0 if the provided datafile is consistent with 'this', another value otherwise
    */
    int CheckSpecificConsistencyWithAnotherDataFile(vDataFile* ap_DataFile);

  // -------------------------------------------------------------------
  // Data members
  private:
    uint16_t m_maxNumberOfLinesPerEvent; /*!< Number of lines in each event in the datafile. Default = 1 */
    FLTNB  m_maxAxialDiffmm;             /*!< Max axial difference in mm between 2 crystals in a LOR. Default value calculated from the scanner files */
    string m_isotope;                    /*!< Isotope. Default = unknown */
    bool m_eventKindFlag;                /*!< Flag for informations about the event nature (true, scatter, random) in the data. Default = false */
    bool m_atnCorrectionFlag;            /*!< Flag that says if attenuation correction terms are included in the data. Default = false */
    bool m_ignoreAttnCorrectionFlag;     /*!< Flag to say if we ignore the attenuation correction even if present. Default = false */
    bool m_normCorrectionFlag;           /*!< Flag that says if normalization correction terms are included in the data. Default = false */
    bool m_ignoreNormCorrectionFlag;     /*!< Flag to say if we ignore the normalization correction even if present. Default = false */
    bool m_scatCorrectionFlag;           /*!< Flag that says if scatter correction terms are included in the data. Default = false */
    bool m_ignoreScatCorrectionFlag;     /*!< Flag to say if we ignore the scatter correction even if present. Default = false */
    bool m_randCorrectionFlag;           /*!< Flag that says if random correction terms are included in the data. Default = false */
    bool m_ignoreRandCorrectionFlag;     /*!< Flag to say if we ignore the random correction even if present. Default = false */
    bool m_TOFInfoFlag;                  /*!< Flag that says if TOF information is included in the data. Default = false */
    bool m_ignoreTOFFlag;                /*!< Flag to say if we ignore the TOF data even if present, or not. Default = false */
    FLTNB m_resolutionTOF;               /*!< TOF resolution in ps. Default = -1.0 */
    int m_nbBinsTOF;                     /*!< Number of TOF bins for histogram mode. Default = 1 */
    FLTNB m_binSizeTOF;                  /*!< Size of TOF bins for histogram mode. Default = -1.0 */
    FLTNB m_TOFMeasurementRange;         /*!< Maximum range of values for TOF delta time measurements (delta t max - delta t min) in ps*/
};

#endif
