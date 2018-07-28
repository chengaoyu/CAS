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
  \brief    Declaration of class vDataFile
*/

#ifndef VDATAFILE_HH
#define VDATAFILE_HH 1

#include "gVariables.hh"
#include "vEvent.hh"
#include "iEventPET.hh"
#include "iEventListPET.hh"
#include "iEventHistoPET.hh"
#include "iEventSPECT.hh"
#include "iEventCT.hh"
#include "iEventNorm.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oMemoryMapped.hh"

/**
 * @defgroup DATA_MODE Data mode
 *
 *    \brief DataFile modes (histogram, list, ...) \n
 *           Defined in vDataFile.hh
 * @{
 */
/** Constant corresponding to an event with undefined mode (=-1) */
#define MODE_UNKNOWN -1
/** Constant corresponding to an event of List mode (=0) */
#define MODE_LIST 0
/** Constant corresponding to an event of Histogram mode (=1) */
#define MODE_HISTOGRAM 1
/** Constant corresponding to a normalization value (=2) */
#define MODE_NORMALIZATION 2
/** @} */

/**
 * @defgroup DATA_TYPE Data type
 *
 *    \brief Data modality (PET, SPECT, CT, PETTrans ...) \n
 *           Defined in vDataFile.hh
 * @{
 */
/** Constant corresponding to an event with undefined type (=-1) */
#define TYPE_UNKNOWN -1
/** Constant corresponding to an event of PET type (=0) */
#define TYPE_PET 0
/** Constant corresponding to an event of SPECT type (=1) */
#define TYPE_SPECT 1
/** Constant corresponding to an event of CT type (=2) */
#define TYPE_CT 2
/** @} */

/**
 * @defgroup DATA_SPEC Data physical specificity
 *
 *    \brief Data physical specificity (EMISSION or TRANSMISSION) \n
 *           Defined in vDataFile.hh
 * @{
 */
/** Constant corresponding to an unknown event (=-1) */
#define SPEC_UNKNOWN -1
/** Constant corresponding to an emission event (=0) */
#define SPEC_EMISSION 0
/** Constant corresponding to a transmission event (=1) */
#define SPEC_TRANSMISSION 1
/** @} */


/*!
  \class   vDataFile
  \brief   This class is designed to be a mother virtual class for DataFile
  \details This class manages the reading of the generic input file (header + data). \n
           It uses some events as buffers to get data informations during run-time.
*/
class vDataFile
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   vDataFile constructor. 
      \details Initialize the member variables to their default values.
    */
    vDataFile ();
    /*!
      \brief   vDataFile destructor. 
    */
    virtual ~vDataFile();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      vDataFile::ReadInfoInHeader()
      \param   bool a_affectQuantificationFlag = true
      \brief   Read and check general information from the header datafile \n
               Call the ReadSpecificInformationInHeader() function implemented in child classes
      \details If the parameter flag is on, then affect the quantification factors from the oImageDimensionsAndQuantification after
               reading relevant information
      \return  0 if success, and positive value otherwise.
    */
    int ReadInfoInHeader(bool a_affectQuantificationFlag = true);
    /*!
      \fn      vDataFile::SetParametersFrom()
      \brief   Initialize all parameters from the provided datafile.
      \details This function is intended to be used right after the constructor
      \return  0 if success, and positive value otherwise
    */
    int SetParametersFrom(vDataFile* ap_DataFile);
    /*!
      \fn      vDataFile::CheckParameters()
      \brief   Check the initialization of member variables \n
               Call the CheckSpecificParameters() function implemented in child classes
      \return  0 if success, and positive value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      vDataFile::InitializeMappedFile()
      \brief   Check the datafile existency, map it to memory and get the raw char* pointer. \n
      \return  0 if success, and positive value otherwise.
    */
    int InitializeMappedFile();
    /*!
      \fn      vDataFile::ComputeSizeEvent() = 0
      \brief   This function is implemented in child classes \n
               Computation of the size of each event according to the mandatory/optional correction fields
      \return  0 is success, positive value otherwise
    */
    virtual int ComputeSizeEvent() = 0;
    /*!
      \fn      vDataFile::PrepareDataFile() = 0
      \brief   This function is implemented in child classes \n
               Store different kind of information inside arrays (data relative to specific correction as well as basic raw data for the case data is loaded in RAM) \n
               Use the flag provided by the user to determine how the data has to be sorted (preloaded or read on the fly)
      \return  0 is success, positive value otherwise
    */
    virtual int PrepareDataFile() = 0;
    /*!
      \fn      vDataFile::OpenFileForWriting()
      \param   string a_suffix
      \brief   Open a binary file stream for writing, with eventually the suffix appended to the file name
      \details Use the m2p_dataFile for that purpose. This function along with the CloseFile function use only one thread
      \return  0 upon success, another value otherwise
    */
    int OpenFileForWriting(string a_suffix="");
    /*!
      \fn      vDataFile::CloseFile()
      \brief   Close as many binary file stream for writing
      \details Use the m2p_dataFile for that purpose, see OpenFileForWriting() function
      \return  0 upon success, another value otherwise
    */
    int CloseFile();
    /*!
      \fn      vDataFile::WriteHeader() = 0
      \brief   This function is implemented in child classes. \n
               Generate a header file according to the data output information.
      \return  0 if success, and positive value otherwise.
    */
    virtual int WriteHeader() = 0;
    /*!
      \fn      vDataFile::WriteEvent() = 0
      \param   ap_Event : Event to write
      \param   a_th : index of the thread from which the function was called
      \brief   This function is implemented in child classes \n
               Write event according to the chosen type of data
      \return  0 if success, and positive value otherwise.
    */
    virtual int WriteEvent(vEvent* ap_Event, int a_th=0) = 0;
    /*!
      \fn      vDataFile::GetEvent()
      \param   int64_t a_eventIndex
      \param   int a_th
      \brief   
      \details x
    */
    vEvent* GetEvent(int64_t a_eventIndex, int a_th = 0);
    /*!
      \fn      vDataFile::GetEventSpecific() = 0
      \param   ap_buffer : address pointing to the event to recover in the event buffer array
      \param   a_th : index of the thread from which the function was called
      \brief   This function is implemented in child classes \n
               Read an event from the position pointed by 'ap_buffer', parse the generic or modality-specific information, and store them in the (multithreaded) 'm2p_BufferEvent' object
      \return  the thread-specific  'm2p_BufferEvent' object containing the modality-specific information for the event
    */
    virtual vEvent* GetEventSpecific(char* ap_buffer, int a_th) = 0;
    /*!
      \fn      vDataFile::GetEventIndexStartAndStop()
      \param   ap_indexStart : pointer to recover the index of the first event
      \param   ap_indexStop : pointer to recover the index of the last event 
      \param   a_subsetNum : actual subset index of the iteration (0 if no subsets)
      \param   a_nbSubsets : max number of subsets in this iteration (1 if no subsets)
      \brief   Compute the index start and stop of the events loop with respect to the current subset and MPI size and rank
    */
    void GetEventIndexStartAndStop(int64_t* ap_indexStart, int64_t* ap_indexStop, int a_subsetNum = 0, int a_NbSubsets = 1);
    /*!
      \fn      vDataFile::CheckConsistencyWithAnotherBedDataFile()
      \param   vDataFile* ap_DataFile
      \brief   Check consistency between 'this' and the provided datafile as two bed positions
      \details It checks data type, mode, number of events if histogram, event size, calibration factor, and scanner.
               For characteristics specific to the modality, it finally calls the pure virtual
               CheckSpecificConsistyWithAnotherBedDataFile() function implemented by children.
      \return  0 if the provided datafile is consistent with 'this', another value otherwise
    */
    int CheckConsistencyWithAnotherBedDataFile(vDataFile* ap_DataFile);
    /*!
      \fn      vDataFile::Describe()
      \brief   A function used to describe the generic parts of the datafile
    */
    void Describe();
    /*!
      \fn      vDataFile::DescribeSpecific() = 0
      \brief   A pure virtual function used to describe the specific parts of the datafile
    */
    virtual void DescribeSpecific() = 0;

  // -------------------------------------------------------------------
  // Functions dedicated to analytical projections
  public:
    /*!
      \fn      vDataFile::PROJ_InitFile() = 0
      \brief   This function is implemented in child classes \n
               Initialize the fstream objets for output writing as well as some other variables specific to the Projection script
      \return  0 if success, and positive value otherwise.
    */
    virtual int PROJ_InitFile() = 0;
    /*!
      \fn      vDataFile::PROJ_GetScannerSpecificParameters() = 0
      \brief   This function is implemented in child classes \n
               It is used to set several variables of the datafile when using the projection script. \n
               Get modality specific parameters from the scanner object, through the scannerManager.
      \return  0 if success, positive value otherwise
    */
    virtual int PROJ_GetScannerSpecificParameters() = 0;
    /*!
      \fn      vDataFile::PROJ_WriteData()
      \brief   Write/Merge chunk of data in a general data file. 
      \todo    adapt to the data loading/writing in RAM
      \return  0 if success, and positive value otherwise.
    */
    int PROJ_WriteData();
    /*!
      \fn      vDataFile::PROJ_DeleteTmpDataFile()
      \brief   Delete temporary datafile used for multithreaded output writing if needed
      \todo    More checks (in this functions and in the calls to this function)
      \return  0 if success, and positive value otherwise.
    */
    int PROJ_DeleteTmpDataFile();
    /*!
      \fn      vDataFile::PROJ_GenerateEvent()
      \param   idx_elt1 : first ID of the event
      \param   idx_elt2 : second ID of the event
      \param   a_th : index of the thread from which the function was called
      \brief   Generate a standard event and set up its ID \n
               Used by the projection, list-mode sensitivity generation, and datafile converter scripts
      \return  the thread specific m2p_BufferEvent array containing the event
    */
    vEvent* PROJ_GenerateEvent(int idx_elt1, int idx_elt2, int a_th);

  // -------------------------------------------------------------------
  // Public Get & Set FUNCTIONS
  public:
    /*!
      \fn      inline int vDataFile::GetBedIndex()
      \return  its bed index
    */
    inline int GetBedIndex()
           {return m_bedIndex;}
    /*!
      \fn      inline int vDataFile::GetDataMode()
      \return  data mode
    */
    inline int GetDataMode()
           {return m_dataMode;}
    /*!
      \fn      string vDataFile::GetDataModeToString()
      \return  The data mode as a human-readable string
    */
    string GetDataModeToString();
    /*!
      \fn      inline int vDataFile::GetDataType()
      \return  data type
    */
    inline int GetDataType()
           {return m_dataType;}
    /*!
      \fn      string vDataFile::GetDataTypeToString()
      \return  The data type as a human-readable string
    */
    string GetDataTypeToString();
    /*!
      \fn      inline int vDataFile::GetDataSpec()
      \return  data spec
    */
    inline int GetDataSpec()
           {return m_dataSpec;}
    /*!
      \fn      string vDataFile::GetDataSpecToString()
      \return  The data spec as a human-readable string
    */
    string GetDataSpecToString();
    /*!
      \fn      inline int64_t vDataFile::GetSize()
      \return  number of events in the datafile
    */
    int64_t GetSize()
         {return m_nbEvents;}
    /*!
      \fn      inline int64_t vDataFile::GetSizeEvent()
      \return  the event size
    */
    int64_t GetEventSize() {return m_sizeEvent;}
    /*!
      \fn      inline string vDataFile::GetHeaderDataFileName()
      \return  headerdatafile name
    */
    inline string GetHeaderDataFileName()
           {return m_headerFileName;};
    /*!
      \fn      inline string vDataFile::GetDataFileName()
      \return  datafile name
    */
    inline string GetDataFileName()
           {return m_dataFileName;};
    /*!
      \fn      inline int64_t vDataFile::GetStartTime()
      \return  int64_t corresponding to the acquisition start time (s)
    */
    inline int64_t GetStartTime()
           {return m_startTimeInSec;};
    /*!
      \fn      inline int64_t vDataFile::GetDuration()
      \return  int64_t corresponding to the acquisition duration (s)
    */
    inline int64_t GetDuration()
           {return m_durationInSec;};
    /*!
      \fn      inline FLTNB vDataFile::GetCalibrationFactor()
      \return  calibration factor
    */
    inline FLTNB GetCalibrationFactor()
           {return m_calibrationFactor;}
    /*!
      \fn      inline FLTNB* vDataFile::GetPOIResolution()
      \return  pointer to the array containing POI (position of interaction) resolution
    */
    inline FLTNB* GetPOIResolution()
           {return mp_POIResolution;}
    /*!
      \fn      inline bool* vDataFile::GetPOIDirectionFlag()
      \return  pointer to the array containing POI (position of interaction) direction flags
    */
    inline bool* GetPOIDirectionFlag()
           {return mp_POIDirectionFlag;}
    /*!
      \fn      inline bool vDataFile::GetPOIInfoFlag()
      \return  m_POIInfoFlag
    */
    inline bool GetPOIInfoFlag()
           {return m_POIInfoFlag;}
    /*!
      \fn      inline bool vDataFile::GetIgnorePOIFlag()
      \return  m_ignorePOIFlag
    */
    inline bool GetIgnorePOIFlag()
           {return m_ignorePOIFlag;}
    /*!
      \fn      virtual int vDataFile::GetMaxRingDiff()
      \brief   Return an error by default. \n
               This function is surcharged by the PET (and CT) scanner daughter class
      \return  -1 (error) if not surcharged by a daughter class
    */
    virtual int GetMaxRingDiff();
    /*!
      \fn      inline void vDataFile::SetDataMode()
      \param   a_dataMode
      \brief   set the data mode
    */
    inline void SetDataMode(int a_dataMode)
           {m_dataMode = a_dataMode;}
    /*!
      \fn      inline void vDataFile::SetDataType()
      \param   a_dataType
      \brief   set the data type
    */
    inline void SetDataType(int a_dataType)
           {m_dataType = a_dataType;}
    /*!
      \fn      inline void vDataFile::SetBedIndex()
      \param   a_bedIndex
      \brief   set the bed index corresponding to this data file
    */
    inline void SetBedIndex(int a_bedIndex)
           {m_bedIndex = a_bedIndex;}
    /*!
      \fn      inline bool vDataFile::GetBedPositionFlag()
      \return  The m_bedPositionFlag
    */
    inline bool GetBedPositionFlag()
           {return m_bedPositionFlag;}
    /*!
      \fn      inline FLTNB vDataFile::GetRelativeBedPosition()
      \return  The m_relativeBedPositionInMm
    */
    inline FLTNB GetRelativeBedPosition()
           {return m_relativeBedPosition;}
    /*!
      \fn      inline void vDataFile::SetVerbose()
      \param   a_verboseLevel
      \brief   set verbosity
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      inline int vDataFile::GetVerbose()
      \brief   Get the verbose level
      \return  m_verbose
    */
    inline int GetVerbose()
           {return m_verbose;}
    /*!
      \fn      inline void vDataFile::SetImageDimensionsAndQuantification()
      \param   ap_ImageDimensionsAndQuantification
      \brief   set the pointer to the oImageDimensionsAndQuantification object
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ID = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      inline void vDataFile::SetPOIResolution()
      \param   ap_value : vector of 3 elements (x,y,z)
      \brief   initialize the POI resolution (for list-mode)
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetPOIResolution(FLTNB ap_value[3])
           {for (int i=0 ; i<3 ; i++) mp_POIResolution[i] = ap_value[i];} 
    /*!
      \fn      inline void vDataFile::SetIgnorePOIFlag()
      \param   a_flag
      \brief   Set a boolean that that if we ignore POI information or not
    */
    inline void SetIgnorePOIFlag(bool a_ignorePOIFlag)
           {m_ignorePOIFlag = a_ignorePOIFlag;}
    /*!
      \fn      inline void vDataFile::SetHeaderDataFileName()
      \param   const string& a_headerFileName
      \brief   set the data header file name
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetHeaderDataFileName(const string& a_headerFileName)
           {m_headerFileName = a_headerFileName;}
    /*!
      \fn      inline void vDataFile::SetBinaryDataFileName()
      \param   const string& a_dataFileName
      \brief   set the data binary file name
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetBinaryDataFileName(const string& a_dataFileName)
           {m_dataFileName = a_dataFileName;}
    /*!
      \fn      inline void vDataFile::SetCalibrationFactor()
      \param   a_value
      \brief   initialize the global calibration factor with a FLTNB value
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetCalibrationFactor(FLTNB a_value)
           {m_calibrationFactor = a_value;} 
    /*!
      \fn      inline void vDataFile::SetNbEvents()
      \param   a_value
      \brief   initialize the number of events with a int64_t value
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetNbEvents(int64_t a_value)
           {m_nbEvents = a_value;} 
    /*!
      \fn      inline void vDataFile::SetStartTimeInSec()
      \param   a_value
      \brief   initialize the acquisition start time (s) with a FLTNB value
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetStartTime(FLTNB a_value)
           {m_startTimeInSec = a_value;} 
    /*!
      \fn      inline void vDataFile::SetDurationInSec()
      \param   a_value
      \brief   initialize the acquisition duration (s) with a FLTNB value
      \details This function is dedicated to datafile conversion scripts
    */
    inline void SetDuration(FLTNB a_value)
           {m_durationInSec = a_value;}
    /*!
      \fn      inline string vDataFile::GetScannerName()
      \return  the scanner name
    */
    inline string GetScannerName()
           {return m_scannerName;}

    virtual int Shuffle( int64_t );

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      virtual int vDataFile::SetSpecificParametersFrom()
      \brief   This function is implemented in child classes \n
               Initialize all specific parameters from the provided datafile.
      \return  0 if success, and positive value otherwise
    */
    virtual int SetSpecificParametersFrom(vDataFile* ap_DataFile) = 0;
    /*!
      \fn      virtual int vDataFile::CheckSpecificParameters() = 0
      \brief   This function is implemented in child classes \n
               Check specific parameters of child classes
      \return  0 if success, and positive value otherwise.
    */
    virtual int CheckSpecificParameters() = 0;
    /*!
      \fn      virtual int vDataFile::CheckFileConsistency() = 0
      \brief   This function is implemented in child classes \n
               Check if file size is consistent.
      \return  0 if success, and positive value otherwise.
    */
    virtual int CheckFileSizeConsistency() = 0;
    /*!
      \fn      virtual int vDataFile::ReadSpecificInfoInHeader() = 0
      \param   bool a_affectQuantificationFlag = true
      \brief   This function is implemented in child classes \n
               Read and check modality-specific information from the header datafile
      \details If the parameter flag is on, then affect the quantification factors from the oImageDimensionsAndQuantification after
               reading relevant information
      \return  0 if success, and positive value otherwise.
    */
    virtual int ReadSpecificInfoInHeader(bool a_affectQuantificationFlag = true) = 0;
    /*!
      \fn      vDataFile::CheckSpecificConsistencyWithAnotherDataFile()
      \param   vDataFile* ap_DataFile
      \brief   Check consistency between 'this' and the provided datafile, for specific characteristics.
      \details Pure virtual function implemented by children. It checks correction flags, etc.
      \return  0 if the provided datafile is consistent with 'this', another value otherwise
    */
    virtual int CheckSpecificConsistencyWithAnotherDataFile(vDataFile* ap_DataFile) = 0;

  // -------------------------------------------------------------------
  // Data members
  protected:
    oImageDimensionsAndQuantification* mp_ID; /*!< Pointer to the oImageDimensionsAndQuantification object */
    int m_verbose;                            /*!< Verbosity */
    
    // Variables related to the acquisition
    string m_headerFileName;                  /*!< String containing the path to the header file */
    string m_dataFileName;                    /*!< String containing the path to the raw datafile */
    int64_t m_nbEvents;                       /*!< Total number of events in the raw data */
    int m_dataMode;                           /*!< Flag indicating if the data is List (=0) or Histogram (=1) mode */
    int m_dataType;                           /*!< Flag indicating if the data is PET (=0),SPECT (=1) or TRANSMISSION type (=2) */
    int m_dataSpec;                           /*!< Flag indicating the physical specificity of the data: SPEC_EMISSION or SPEC_TRANSMISSION */
    FLTNB m_startTimeInSec;                   /*!< Start time of the acquisition (s) */
    FLTNB m_durationInSec;                    /*!< Duration of the acquisition (s) */
    FLTNB m_calibrationFactor;                /*!< Calibration factor for the data. Default value =1.0 */
    int m_bedIndex;                           /*!< Bed position index corresponding to this data file */
    FLTNB m_relativeBedPosition;              /*!< Bed relative position in mm */
    bool m_bedPositionFlag;                   /*!< Flag indicating that a relative bed position has been provided */
    string m_scannerName;                     /*!< Scanner name */
    
    // POI: Position Of Interaction
    bool m_POIInfoFlag;                       /*!< Flag to say if POI information is included in the datafile for each event */
    bool m_ignorePOIFlag;                     /*!< Flag to say if we ignore the POI data if present, or not. Default = false */
    bool mp_POIDirectionFlag[3];              /*!< Flag to say which direction is included in the POI for each event; radial, tangential and axial */
    FLTNB mp_POIResolution[3];                /*!< POI resolution (position of interaction) for each direction: radial, tangential and axial */
    
    // Members for I/O actions
    int64_t m_sizeEvent;                      /*!< Size of an event in the datafile (calculated from mandatory and optional fields) */
    fstream** m2p_dataFile;                   /*!< File associated to the raw data file (multithreaded) */
    vEvent** m2p_BufferEvent;                 /*!< vEvent structure, used to read and transfer the raw data to each part of the algorithm (multithreaded) */
    int64_t m_mpi1stEvent;                    /*!< First index managed by this MPI instance */
    int64_t m_mpiLastEvent;                   /*!< Last index (included) managed by this MPI instance */
    int64_t m_mpiNbEvents;                    /*!< Number of events managed by this MPI instance */
    oMemoryMapped* mp_MappedFile;             /*!< The object managing the mapping of the datafile */
    char* mp_mappedMemory;                    /*!< The raw pointer directly mapped to the datafile */
};

#endif
