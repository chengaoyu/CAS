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
  \brief    Declaration of class iDataFileSPECT
*/

#ifndef IDATAFILESPECT_HH
#define IDATAFILESPECT_HH 1

#include "gVariables.hh"
#include "vDataFile.hh"
#include "vScanner.hh"
#include "iEventHistoSPECT.hh"
#include "iEventListSPECT.hh"

/*!
  \class   iDataFileSPECT
  \brief   Inherit from vDataFile. Class that manages the reading of a SPECT input file (header + data).
  \details It contains several arrays corresponding to the different kind of informations the data file could contain. \n
           As many booleans as arrays say if the data are here or not. The data file can be either completely loaded, or read event by event during reconstruction. \n
           MPI is coming here to cut the data file into peaces (also either can be loaded or read on-the-fly).
*/
class iDataFileSPECT : public vDataFile
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iDataFileSPECT constructor. 
               Initialize the member variables to their default values.
    */
    iDataFileSPECT();
    /*!
      \brief   iDataFileSPECT destructor. 
    */
    ~iDataFileSPECT();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      iDataFileSPECT::ReadSpecificInfoInHeader()
      \param   bool a_affectQuantificationFlag
      \brief   Read through the header file and recover specific SPECT information.
      \details If the parameter flag is on, then affect the quantification factors from the oImageDimensionsAndQuantification after
               reading relevant information
      \return  0 is success, positive value otherwise
    */
    int ReadSpecificInfoInHeader(bool a_affectQuantificationFlag);
    /*!
      \fn      iDataFileSPECT::WriteHeader()
      \brief   Generate a header file according to the data output information.
      \return  0 if success, and positive value otherwise.
    */
    int WriteHeader();
    /*!
      \fn      iDataFileSPECT::ComputeSizeEvent()
      \brief   Computation of the size of each event according to the mandatory/optional correction fields
      \return  0 is success, positive value otherwise
    */
    int ComputeSizeEvent();
    /*!
      \fn      iDataFileSPECT::PrepareDataFile()
      \brief   Store different kind of information inside arrays (data relative to specific correction as well as basic raw data for the case data is loaded in RAM) \n
               Use the flag provided by the user to determine how the data has to be sorted (preloaded or read on the fly)
      \return  0 is success, positive value otherwise
    */
    int PrepareDataFile();
    /*!
      \fn      iDataFileSPECT::WriteEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write event according to the chosen type of data
      \return  0 if success, and positive value otherwise.
    */
    int WriteEvent(vEvent* ap_Event, int a_th=0);
    /*!
      \fn      iDataFileSPECT::GetEventSpecific()
      \param   ap_buffer : address pointing to the event to recover
      \param   a_th : index of the thread from which the function was called
      \brief   Read an event from the position pointed by 'ap_buffer', parse the generic or modality-specific information, and store them in the (multithreaded) 'm2p_BufferEvent' object
      \return  the thread-specific  'm2p_BufferEvent' object containing the modality-specific information for the event
    */
    vEvent* GetEventSpecific(char* ap_buffer, int a_th);
    /*!
      \fn      iDataFileSPECT::InitAngles()
      \param   ap_angles
      \brief   allocate memory for the mp_angles variable using m_nbProjections
               and initialize the projection angles with the provided list of values
      \return  0 if success, positive value otherwise
    */
    int InitAngles(FLTNB* ap_angles);
    /*!
      \fn      iDataFileSPECT::InitCorToDetectorDistance()
      \param   ap_CORtoDetectorDistance
      \brief   allocate memory for the ap_CORtoDetectorDistance variable
               using m_nbProjections, and initialize the projection angles
               with the provided list of values
      \return  0 if success, positive value otherwise
    */
    int InitCorToDetectorDistance(FLTNB* ap_CORtoDetectorDistance);
    /*!
      \fn      iDataFileSPECT::DescribeSpecific()
      \brief   Implementation of the pure virtual eponym function that simply prints info about the datafile
    */
    void DescribeSpecific();

  // -------------------------------------------------------------------
  // Public functions dedicated to the projection script
  public:
    /*!
      \fn      iDataFileSPECT::PROJ_InitFile()
      \brief   Initialize the fstream objets for output writing as well as some other variables specific to the Projection script (Event-based correction flags, Estimated size of data file)
      \return  0 if success, and positive value otherwise.
    */
    int PROJ_InitFile();
    /*!
      \fn      iDataFileSPECT::PROJ_GetScannerSpecificParameters()
      \brief   Get SPECT specific parameters for projections from the scanner object, through the scannerManager.
      \return  0 if success, positive value otherwise
    */
    int PROJ_GetScannerSpecificParameters();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      iDataFileSPECT::GetNbHeads() 
      \return  the number of heads
    */
    inline uint16_t GetNbHeads()
           {return m_nbHeads;};
    /*!
      \fn      iDataFileSPECT::GetNbProjections() 
      \return  total number of projections in the SPECT acquisition
    */
    inline uint16_t GetNbProjections()
           {return m_nbOfProjections;};
    /*!
      \fn      iDataFileSPECT::GetAngles() 
      \return  pointers to angles of projections
    */
    inline FLTNB* GetAngles()
           {return mp_angles;};
    /*!
      \fn      iDataFileSPECT::GetCORtoDetectorDistance() 
      \return  pointers to the COR to detector distances for each projection
    */
    inline FLTNB* GetCORtoDetectorDistance()
           {return mp_CORtoDetectorDistance;};
    /*!
      \fn      iDataFileSPECT::GetNbBins()
      \param   axis : axis corresponding to transaxial (0) or axial (1) bins 
      \return  number of bins corresponding to the axis passed in parameter
    */
    inline uint16_t GetNbBins(int axis)
           {return mp_nbOfBins[axis];};
    /*!
      \fn      iDataFileSPECT::SetEventKindFlagOn()
      \brief   set to true the flag indicating the presence of the kind of a list-mode event in the datafile
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetEventKindFlagOn()
           {m_eventKindFlag = true;}
    /*!
      \fn      iDataFileSPECT::SetIsotope()
      \param   a_value
      \brief   initialize the isotope string value
      \details The name should corresponds to one corresponding tag in the isotope configuration file in config/.
               This function is dedicated to datafile conversion scripts
    */
    inline void SetIsotope(string a_value)
           {m_isotope = a_value;} 
    /*!
      \fn      iDataFileSPECT::GetIsotope()
      \return  the isotope string value
    */
    inline string GetIsotope()
           {return m_isotope;}
    /*!
      \fn      iDataFileSPECT::SetNbBins()
      \param   a_binTrs
      \param   a_binAxl
      \brief   initialize the bin values
    */
    inline void SetNbBins(uint16_t a_binTrs, uint16_t a_binAxl)
           {mp_nbOfBins[0] = a_binTrs; mp_nbOfBins[1]=a_binAxl;} 
    /*!
      \fn      iDataFileSPECT::SetNbProjections()
      \param   a_nbProjections
      \brief   initialize the number of projections
    */
    inline void SetNbProjections(uint16_t a_nbProjections)
           {m_nbOfProjections = a_nbProjections;} 
    /*!
      \fn      iDataFileSPECT::GetHeadRotDirection()
      \brief   Simply return m_headRotDirection
      \return  m_headRotDirection
    */
    inline int GetHeadRotDirection()
           {return m_headRotDirection;}
    /*!
      \fn      iDataFileSPECT::SetHeadRotDirection()
      \param   a_direction
      \brief   initialize the rotation direction of the gamma camera(s)
    */
    inline void SetHeadRotDirection(int a_direction)
           {m_headRotDirection = a_direction;}
    /*!
      \fn      iDataFileSPECT::GetEventKindFlag()
      \brief   Simply return m_eventKindFlag
      \return  m_eventKindFlag
    */
    inline bool GetEventKindFlag()
           {return m_eventKindFlag;}
    /*!
      \fn      iDataFileSPECT::GetScatCorrectionFlag()
      \brief   Simply return m_scatCorrectionFlag
      \return  m_scatCorrectionFlag
    */
    inline bool GetScatCorrectionFlag()
           {return m_scatCorrectionFlag;}
    /*!
      \fn      iDataFileSPECT::GetNormCorrectionFlag()
      \brief   Simply return m_normCorrectionFlag
      \return  m_normCorrectionFlag
    */
    inline bool GetNormCorrectionFlag()
           {return m_normCorrectionFlag;}

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      int iDataFileSPECT::SetSpecificParametersFrom()
      \brief   Initialize all parameters specific to SPECT from the provided datafile.
      \return  0 if success, and positive value otherwise
    */
    int SetSpecificParametersFrom(vDataFile* ap_DataFile);
    /*!
      \fn      iDataFileSPECT::CheckSpecificParameters()
      \brief   Check parameters specific to SPECT data
      \return  0 if success, and positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iDataFileSPECT::CheckFileSizeConsistency()
      \brief   This function is implemented in child classes \n
               Check if file size is consistent.
      \return  0 if success, and positive value otherwise.
    */
    int CheckFileSizeConsistency();
    /*!
      \fn      iDataFileSPECT::WriteHistoEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write a SPECT histogram event
      \return  0 if success, and positive value otherwise.
    */
    int WriteHistoEvent(iEventHistoSPECT* ap_Event, int a_th);
    /*!
      \fn      iDataFileSPECT::WriteListEvent()
      \param   ap_Event : event containing the data to write
      \param   a_th : index of the thread from which the function was called
      \brief   Write a SPECT list-mode event
      \return  0 if success, and positive value otherwise.
    */
    int WriteListEvent(iEventListSPECT* ap_Event, int a_th);
    /*!
      \fn      iDataFileSPECT::CheckSpecificConsistencyWithAnotherDataFile()
      \param   vDataFile* ap_DataFile
      \brief   Check consistency between 'this' and the provided datafile, for specific characteristics.
      \details Implementation of the pure virtual function from vDataFile. It checks correction flags, etc.
      \return  0 if the provided datafile is consistent with 'this', another value otherwise
    */
    int CheckSpecificConsistencyWithAnotherDataFile(vDataFile* ap_DataFile);

  // -------------------------------------------------------------------
  // Data members
  private:
    string m_isotope;                /*!< Isotope. Default value =unknown */
    bool m_eventKindFlag;            /*!< Flag for informations about the event nature (true, scatter) in the data. Default value = false */
    bool m_normCorrectionFlag;       /*!< Flag that says if normalization correction terms are included in the data. Default = false */
    bool m_ignoreNormCorrectionFlag; /*!< Flag to say if we ignore the normalization correction even if present. Default = false */
    bool m_scatCorrectionFlag;       /*!< Flag that says if scatter correction terms are included in the data. Default = false */
    bool m_ignoreScatCorrectionFlag; /*!< Flag to say if we ignore the scatter correction even if present. Default = false */
    uint16_t mp_nbOfBins[2];         /*!< Transaxial/Axial number of bins. Default value =1,1 */
    FLTNB m_acquisitionZoom;         /*!< Zoom used during the acquisition to limit the area of detection for monolithic detectors */
    uint16_t m_nbOfProjections;      /*!< Total number of projections during the acquisition(for all the heads). No Default*/
    FLTNB* mp_angles;                /*!< Angle [for each projection] in degrees. If SPECT system contains several heads, first head angles should be entered first, 
                                          followed by 2nd head angles, etc.. No Default */
    uint16_t m_nbHeads;              /*!< Number of heads in the SPECT systems. Default =1*/
    FLTNB* mp_CORtoDetectorDistance; /*!< Distance camera surface to COR (mm) [for each projection]. \n
                                        if not provided, the distance given for each heads in the camera description file is taken and considered constant for each projections related to each head \n
                                        if provided then: if positive value (either a constant value, or a value specific to each projection) then it overwrites the one given in the camera file \n
                                                          if negative, the distance given for each heads in the camera description file is taken and considered constant for each projections related to each head \n
                                        Default value = Recovered from the camera description file */
    int m_headRotDirection;          /*!< Head rotation direction (0=clockwise, 1=counterclockwise)*/
};

#endif
