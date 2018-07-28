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

  \brief Implementation of class vDataFile
*/

#include "vDataFile.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vDataFile::vDataFile() 
{
  mp_ID = NULL;
  m_verbose = -1;

  // Variables related to the acquisition
  m2p_dataFile = NULL;
  m_headerFileName = "";
  m_dataFileName = "";
  m_scannerName = "";
  m_dataMode = MODE_UNKNOWN;
  m_dataType = TYPE_UNKNOWN;
  m_dataSpec = SPEC_UNKNOWN;
  m_nbEvents = -1;
  m_bedIndex = -1;
  m_relativeBedPosition = 0.;
  m_bedPositionFlag = false;
  m_startTimeInSec = 0.;
  m_durationInSec = -1.;
  m_calibrationFactor = 1.;
  m_sizeEvent = -1;

  // Default POI (meaning we do not have any)
  mp_POIResolution[0] = -1.;
  mp_POIResolution[1] = -1.;
  mp_POIResolution[2] = -1.;
  mp_POIDirectionFlag[0] = false;
  mp_POIDirectionFlag[1] = false;
  mp_POIDirectionFlag[2] = false;
  m_POIInfoFlag = false;
  m_ignorePOIFlag = false;
  
  // Variable related to Buffer/Container arrays
  m2p_BufferEvent = NULL;
  m_mpi1stEvent = -1;
  m_mpiLastEvent = -1;
  m_mpiNbEvents = -1; 
  mp_MappedFile = NULL;
  mp_mappedMemory = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vDataFile::~vDataFile() 
{
  // Free the (vEvent) buffer event
  if (m2p_BufferEvent)
  {
    for(int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
      if (m2p_BufferEvent[th]) delete m2p_BufferEvent[th];
    delete[] m2p_BufferEvent;
  }
  // Close datafiles and delete them
  if (m2p_dataFile)
  {
    for(int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
      if (m2p_dataFile[th]) m2p_dataFile[th]->close();
    delete m2p_dataFile;
  }
  // Delete the mapped memory object
  if (mp_MappedFile!=NULL) delete mp_MappedFile;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::ReadInfoInHeader(bool a_affectQuantificationFlag)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("vDataFile::ReadInfoInHeader() -> Read datafile header from '" << m_headerFileName << " ...'" << endl);

  // Check that the image dimensions and quantification object has been set if the a_affectQuantificationFlag is true
  if (a_affectQuantificationFlag && mp_ID==NULL)
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Image dimensions and quantification object has not been set while it is asked to affect" << endl);
    Cerr("                                       its parameters from information read into the header !" << endl);
    return 1;
  }

  // Read mandatory general fields in the header, check if errors (either mandatory tag not found or issue during data reading/conversion (==1) )
  if (ReadDataASCIIFile(m_headerFileName, "Number of events", &m_nbEvents, 1, KEYWORD_MANDATORY) )
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while trying to read number of events in the header data file !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(m_headerFileName, "Scanner name", &m_scannerName, 1, KEYWORD_MANDATORY) ) 
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while trying to read scanner name in the header data file !" << endl);
    return 1;
  }

  // Get data file name
  if (ReadDataASCIIFile(m_headerFileName, "Data filename", &m_dataFileName, 1, KEYWORD_MANDATORY) ) 
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while trying to read data filename in the header data file !" << endl);
    return 1;
  }
  // If it is an absolute path (start with a OS_SEP) then let it as is, otherwise past the path of the header file name (because we suppose the two
  // files are side by side).
  // This is currently only unix compatible...
  if (m_dataFileName.substr(0,1)!=OS_SEP && m_headerFileName.find(OS_SEP)!=string::npos)
  {
    // Extract the path from the header file name
    size_t last_slash = m_headerFileName.find_last_of(OS_SEP);
    // Paste the path to the data file name
    m_dataFileName = m_headerFileName.substr(0,last_slash+1) + m_dataFileName;
  }

  // For data mode, multiple declarations are allowed, so we get the mode as a string and do some checks
  string data_mode = "";
  if (ReadDataASCIIFile(m_headerFileName, "Data mode", &data_mode, 1, KEYWORD_MANDATORY) ) 
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while trying to read data mode in the header data file !" << endl);
    return 1;
  }
  if ( data_mode=="LIST" || data_mode=="LISTMODE" || data_mode=="LIST-MODE" ||
       data_mode=="list" || data_mode=="listmode" || data_mode=="list-mode" ||  data_mode=="0" ) m_dataMode = MODE_LIST;
  else if ( data_mode=="HISTOGRAM" || data_mode=="histogram" || data_mode=="Histogram" ||
            data_mode=="HISTO" || data_mode=="histo" || data_mode=="Histo" || data_mode=="1" ) m_dataMode = MODE_HISTOGRAM;
  else if ( data_mode=="NORMALIZATION" || data_mode=="normalization" || data_mode=="Normalization" ||
            data_mode=="NORM" || data_mode=="norm" || data_mode=="Norm" || data_mode=="2" ) m_dataMode = MODE_NORMALIZATION;
  else
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Unknown data mode '" << data_mode << "' found in header file !" << endl);
    return 1;
  }

  // For data type, multiple declarations are allowed, so we get the mode as a string and do some checks
  string data_type = "";
  if (ReadDataASCIIFile(m_headerFileName, "Data type", &data_type, 1, KEYWORD_MANDATORY) ) 
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while trying to read data type in the header data file !" << endl);
    return 1;
  }
  if ( data_type=="PET" || data_type=="pet" || data_type=="0" )
  {
    m_dataType = TYPE_PET;
    m_dataSpec = SPEC_EMISSION;
  }
  else if ( data_type=="SPECT" || data_type=="spect" || data_type=="1" )
  {
    m_dataType = TYPE_SPECT;
    m_dataSpec = SPEC_EMISSION;
  }
  else if ( data_type=="CT" || data_type=="ct" || data_type=="2" )
  {
    m_dataType = TYPE_CT;
    m_dataSpec = SPEC_TRANSMISSION;
  }
  else
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Unknown data type '" << data_type << "' found in header file !" << endl);
    return 1;
  }

  // Get start time and duration (optional for the normalization mode
  if (m_dataMode!=MODE_NORMALIZATION)
  {
    if (ReadDataASCIIFile(m_headerFileName, "Start time (s)", &m_startTimeInSec, 1, KEYWORD_MANDATORY) )
    {
      Cerr("***** vDataFile::ReadInfoInHeader() -> Error while trying to read acquisition start time in the header data file !" << endl);
      return 1;
    }
    if (ReadDataASCIIFile(m_headerFileName, "Duration (s)", &m_durationInSec, 1, KEYWORD_MANDATORY) )
    {
      Cerr("***** vDataFile::ReadInfoInHeader() -> Error while trying to read acquisition stop time in the header data file !" << endl);
      return 1;
    }
  }
  else
  {
    if (ReadDataASCIIFile(m_headerFileName, "Start time (s)", &m_startTimeInSec, 1, KEYWORD_OPTIONAL) == 1 )
    {
      Cerr("***** vDataFile::ReadInfoInHeader() -> Error while reading acquisition start time in the header data file !" << endl);
      return 1;
    }
    if (ReadDataASCIIFile(m_headerFileName, "Duration (s)", &m_durationInSec, 1, KEYWORD_OPTIONAL) == 1 )
    {
      Cerr("***** vDataFile::ReadInfoInHeader() -> Error while reading acquisition stop time in the header data file !" << endl);
      return 1;
    }
  }

  // Set the acquisition timing to the quantification factors
  if (a_affectQuantificationFlag && mp_ID->SetAcquisitionTime(m_bedIndex, m_startTimeInSec, m_durationInSec))
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while passing acquisition start and stop time to the ImageDimensionsAndQuantification !" << endl);
    return 1;
  }

  // Read optional bed relative position
  int read_bed_position_status = ReadDataASCIIFile(m_headerFileName, "Horizontal bed relative position (mm)", &m_relativeBedPosition, 1, KEYWORD_OPTIONAL);
  if (read_bed_position_status==KEYWORD_OPTIONAL_SUCCESS) m_bedPositionFlag = true;
  else if (read_bed_position_status==KEYWORD_OPTIONAL_SUCCESS)
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while reading the 'Bed relative position (mm)' field in the data file header !" << endl);
    return 1;
  }

  // Read optional fields in the header, check if errors (issue during data reading/conversion (==1) )
  if (ReadDataASCIIFile(m_headerFileName, "Calibration factor", &m_calibrationFactor, 1, KEYWORD_OPTIONAL) == 1 ||
      ReadDataASCIIFile(m_headerFileName, "POI capability", mp_POIResolution, 3, KEYWORD_OPTIONAL) == 1 ||
      ReadDataASCIIFile(m_headerFileName, "POI correction flag", &m_POIInfoFlag, 3, KEYWORD_OPTIONAL) == 1 )
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while reading optional field in the header data file !" << endl);
    return 1;
  }
  
  // Read fields specific to the modality (call to the related function implemented in child classes)
  if (ReadSpecificInfoInHeader(a_affectQuantificationFlag) )
  {
    Cerr("***** vDataFile::ReadInfoInHeader() -> Error while trying to read modality-specific informations from the header data file !" << endl);
    return 1;
  }

  // Give the calibration factor to the oImageDimensionsAndQuantification that manages the quantification factors
  if (a_affectQuantificationFlag && mp_ID->SetCalibrationFactor(m_bedIndex, m_calibrationFactor))
  {
    Cerr("***** vDataFile::ReadSpecificInfoInHeader() -> A problem occured while setting the calibration factor to oImageDimensionsAndQuantification !" << endl);
    return 1;
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::SetParametersFrom(vDataFile* ap_DataFile)
{
  // Verbose
  m_verbose = ap_DataFile->GetVerbose();
  // Variables related to the acquisition
  m2p_dataFile = NULL;
  m_scannerName = ap_DataFile->GetScannerName();
  m_dataMode = ap_DataFile->GetDataMode();
  m_dataType = ap_DataFile->GetDataType();
  m_dataSpec = ap_DataFile->GetDataSpec();
  m_bedIndex = ap_DataFile->GetBedIndex();
  m_relativeBedPosition = ap_DataFile->GetRelativeBedPosition();
  m_bedPositionFlag = ap_DataFile->GetBedPositionFlag();
  m_startTimeInSec = ap_DataFile->GetStartTime();
  m_durationInSec = ap_DataFile->GetDuration();
  m_calibrationFactor = ap_DataFile->GetCalibrationFactor();
  m_sizeEvent = ap_DataFile->GetEventSize();
  // Default POI (meaning we do not have any)
  m_POIInfoFlag = ap_DataFile->GetPOIInfoFlag();
  mp_POIResolution[0] = ap_DataFile->GetPOIResolution()[0];
  mp_POIResolution[1] = ap_DataFile->GetPOIResolution()[1];
  mp_POIResolution[2] = ap_DataFile->GetPOIResolution()[2];
  mp_POIDirectionFlag[0] = ap_DataFile->GetPOIDirectionFlag()[0];
  mp_POIDirectionFlag[1] = ap_DataFile->GetPOIDirectionFlag()[1];
  mp_POIDirectionFlag[2] = ap_DataFile->GetPOIDirectionFlag()[2];
  // Call the specific function
  if (SetSpecificParametersFrom(ap_DataFile))
  {
    Cerr("***** vDataFile::SetParametersFrom() -> An error occured while setting specific parameters fro mthe provided datafile !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::CheckParameters()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("vDataFile::CheckParameters() -> Check mandatory parameters" << endl);
  // Check mandatory parameters
  if (mp_ID == NULL)
  {
    Cerr("***** vDataFile::CheckParameters() -> ImageDimensionsAndQuantification object not initialized !" << endl);
    return 1;
  }
  if (m_headerFileName == "")
  {
    Cerr("***** vDataFile::CheckParameters() -> String containing path to header file not initialized !" << endl);
    return 1;
  }  
  if (m_dataFileName == "")
  {
    Cerr("***** vDataFile::CheckParameters() -> String containing path to raw data file not initialized !" << endl);
    return 1;
  }
  if (m_nbEvents<0)
  {
    Cerr("***** vDataFile::CheckParameters() -> Number of events incorrectly initialized !" << endl);
    return 1;
  }
  if (m_dataMode!=MODE_LIST && m_dataMode!=MODE_HISTOGRAM && m_dataMode!=MODE_NORMALIZATION)
  {
    Cerr("***** vDataFile::CheckParameters() -> Data mode incorrectly initialized !" << endl);
    return 1;
  }
  if (m_dataMode!= MODE_NORMALIZATION && m_durationInSec<0)
  {
    Cerr("***** vDataFile::CheckParameters() -> Acquisition duration (s) incorrectly initialized !" << endl);
    return 1;
  }
  if (m_dataType!=TYPE_PET && m_dataType!=TYPE_SPECT && m_dataType!=TYPE_CT)
  {
    Cerr("***** vDataFile::CheckParameters() -> Data type incorrectly initialized !" << endl);
    return 1;
  }
  if (m_dataSpec!=SPEC_EMISSION && m_dataSpec!=SPEC_TRANSMISSION)
  {
    Cerr("***** vDataFile::CheckParameters() -> Data physical property incorrectly initialized !" << endl);
    return 1;
  }
  if (m_bedIndex<0)
  {
    Cerr("***** vDataFile::CheckParameters() -> Bed position index incorrectly initialized !" << endl);
    return 1;
  }
  if (m_verbose<0)
  {
    Cerr("***** vDataFile::CheckParameters() -> Verbosity incorrectly initialized !" << endl);
    return 1;
  }
  // Call to the related function implemented in child classes
  if (CheckSpecificParameters())
  {
    Cerr("***** vDataFile::CheckParameters() -> Error while checking specific parameters !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::CheckConsistencyWithAnotherBedDataFile(vDataFile* ap_DataFile)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Check consistency between this and the provided datafiles" << endl);
  // Check data type
  if (m_dataType!=ap_DataFile->GetDataType())
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Data types are inconsistent !" << endl);
    return 1;
  }
  // Check data mode
  if (m_dataMode!=ap_DataFile->GetDataMode())
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Data modes are inconsistent !" << endl);
    return 1;
  }
  // For histogram and normalization modes, check the number of events (i.e. the number of histogram bins)
  if ( (m_dataMode==MODE_HISTOGRAM || m_dataMode==MODE_NORMALIZATION) && m_nbEvents!=ap_DataFile->GetSize())
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Total number of events are inconsistent !" << endl);
    return 1;
  }
  // Check the scanner name
  if (m_scannerName!=ap_DataFile->GetScannerName())
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Scanner names are inconsistent !" << endl);
    return 1;
  }
  // Check the calibration factor
  if (m_calibrationFactor!=ap_DataFile->GetCalibrationFactor())
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Calibration factors are inconsistent !" << endl);
    return 1;
  }
  // Check event size
  if (m_sizeEvent!=ap_DataFile->GetEventSize())
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Events sizes are inconsistent !" << endl);
    return 1;
  }
  // Check POI info flag
  if (m_POIInfoFlag!=ap_DataFile->GetPOIInfoFlag())
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> POI flags are inconsistent !" << endl);
    return 1;
  }
  // Check POI resolutions
  if ( (mp_POIResolution[0]!=ap_DataFile->GetPOIResolution()[0]) ||
       (mp_POIResolution[1]!=ap_DataFile->GetPOIResolution()[1]) ||
       (mp_POIResolution[2]!=ap_DataFile->GetPOIResolution()[2]) )
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> POI resolutions are inconsistent !" << endl);
    return 1;
  }
  // Check the bed position flag
  if (m_bedPositionFlag!=ap_DataFile->GetBedPositionFlag())
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Bed relative positions must be provided for all beds or not at all !" << endl);
    return 1;
  }
  // Call specific function
  if (CheckSpecificConsistencyWithAnotherDataFile(ap_DataFile))
  {
    Cerr("***** vDataFile::CheckConsistencyWithAnotherBedDataFile() -> Inconsistency detected for specific characteristics !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::InitializeMappedFile()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_NORMAL) Cout("vDataFile::InitializeMappedFile() -> Map datafile to memory" << endl);

  // Check file size consistency here (this is a pure virtual function implemented by children)
  if (CheckFileSizeConsistency())
  {
    Cerr("***** vDataFile::InitializeMappedFile() -> A problem occured while checking file size consistency !" << endl);
    return 1;
  }

  // Create mapped file
  mp_MappedFile = new oMemoryMapped();
  // Open the file
  if (mp_MappedFile->Open(m_dataFileName.c_str() , oMemoryMapped::WholeFile, oMemoryMapped::Normal))
  {
    Cerr("***** vDataFile::InitializeMappedFile() -> Failed to open data file '" << m_dataFileName << "' !" << endl);
    Cerr("                                           Provided in data file header '" << m_headerFileName << "'." << endl);
    return 1;
  }

  // Get raw pointer to mapped memory
  mp_mappedMemory = (char*)mp_MappedFile->GetData();

  // ------------------ Compute here the piece of file that each MPI instance manages --------------------

  // 1. Compute the number of events that each instance has to manage
  int64_t instance_size = m_nbEvents / mp_ID->GetMPISize();
  // All instances manage an equal part of the file but the last one which also manages the few last events
  if (mp_ID->GetMPIRank()!=mp_ID->GetMPISize()-1) m_mpiNbEvents = instance_size;
  else m_mpiNbEvents = instance_size + (m_nbEvents - instance_size*mp_ID->GetMPISize());
  // 2. Compute the first event managed by the instance
  m_mpi1stEvent = mp_ID->GetMPIRank() * instance_size;
  // 3. Compute the last event managed by the instance
  m_mpiLastEvent = m_mpi1stEvent + m_mpiNbEvents - 1;

  // End
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vDataFile::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose==0) return;
  // Describe the datafile
  Cout("vDataFile::Describe() -> Here is some generic content of the datafile" << endl);
  Cout("  --> Header file name: " << m_headerFileName << endl);
  Cout("  --> Data file name: " << m_dataFileName << endl);
  Cout("  --> Data type: " << GetDataTypeToString() << endl);
  Cout("  --> Data mode: " << GetDataModeToString() << endl);
  Cout("  --> Data spec: " << GetDataSpecToString() << endl);
  Cout("  --> Scanner name: " << m_scannerName << endl);
  Cout("  --> Number of specific events: " << m_nbEvents << endl);
  Cout("  --> Size of an event: " << m_sizeEvent << " bytes" << endl);
  Cout("  --> Time start: " << m_startTimeInSec << " sec" << endl);
  Cout("  --> Duration: " << m_durationInSec << " sec" << endl);
  Cout("  --> Calibration factor: " << m_calibrationFactor << endl);
  if (m_bedPositionFlag) Cout("  --> Relative axial bed position: " << m_relativeBedPosition << " mm" << endl);
  // Call the specific function
  DescribeSpecific();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::OpenFileForWriting(string a_suffix)
{
  // Check that file is not already open
  if (m2p_dataFile && m2p_dataFile[0]->is_open())
  {
    Cerr("***** vDataFile::OpenFileForWriting() -> Output data file is already open !" << endl);
    return 1;
  }

  // Allocate the m2p_dataFile for each thread
  m2p_dataFile = new fstream*[1];

  // Set the name of output files
  sOutputManager* p_OutMgr = sOutputManager::GetInstance();
  string path_name = p_OutMgr->GetPathName();
  string base_name = p_OutMgr->GetBaseName();
  m_headerFileName = path_name + base_name + a_suffix + ".cdh";
  m_dataFileName = path_name + base_name + a_suffix + ".cdf";

  // Verbose
  if (m_verbose>=2) Cout("vDataFile::OpenFileForWriting() -> Output data file header is '" << m_headerFileName << "'" << endl);

  // Open file   
  m2p_dataFile[0] = new fstream( m_dataFileName.c_str(), ios::binary | ios::out );

  // Check
  if (!m2p_dataFile[0]->is_open())
  {
    Cerr("***** vDataFile::OpenFileForWriting() -> Failed to create output file '" << m_dataFileName << "' !" << endl);
    return 1;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::CloseFile()
{
  // Close binary data file
  if (m2p_dataFile)
  {
    // Verbose
    if (m_verbose>=2) Cout("vDataFile::CloseFile() -> Closing output binary data file" << endl);
    // Close file
    if (m2p_dataFile[0]) m2p_dataFile[0]->close();
    delete m2p_dataFile;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vEvent* vDataFile::GetEvent(int64_t a_eventIndex, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  // The event pointer that will be returned
  vEvent* p_event;
  // Compute the adress into the buffer
  char* event_address_in_array = mp_mappedMemory + a_eventIndex*m_sizeEvent;
  // Get the event
  p_event = GetEventSpecific(event_address_in_array, a_th);
  // Return the event
  return p_event;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

  //!!!!\ This function has been modified to be used specifically with a
  //      dev tool, but it can be used currently by the datafile
  //      Ignored for now in windows compilation

int vDataFile::Shuffle( int64_t nb_events_to_load )
{ 
  #ifndef _WIN32
  if (m_verbose >=2) Cout("vDataFile::Shuffle() : Everyday I shuffle in... " << endl);
  // Buffer storing ordered indices from ( 0 ) to ( nb_events_to_load - 1 )
  int64_t *rndmIdx = new int64_t[ nb_events_to_load ];
  std::iota( rndmIdx, rndmIdx + nb_events_to_load, 0 );

  // Initializing random
  std::random_device rd;
  std::mt19937 rndm( rd() );
  rndm.seed( 1100001001 );

  // Shuffling the buffer of indices
  std::shuffle( rndmIdx, rndmIdx + nb_events_to_load, rndm );

  // Creating a tmp buffer
  char *mp_arrayEvents_tmp = new char[ nb_events_to_load*m_sizeEvent ];

  // Copy sorted buffer to shuffled buffer
  for( int64_t i = 0; i < nb_events_to_load; ++i )
  {
    for( int64_t j = 0; j < m_sizeEvent; ++j )
    {
      //mp_arrayEvents_tmp[ i * m_sizeEvent + j ] = mp_arrayEvents[ rndmIdx[ i ] * m_sizeEvent + j ];
      mp_arrayEvents_tmp[ i * m_sizeEvent + j ] = mp_mappedMemory[ rndmIdx[ i ] * m_sizeEvent + j ];
    }
  }

  // Freeing memory
  delete[] rndmIdx;
  //delete[] mp_mappedMemory;

  //mp_arrayEvents = mp_arrayEvents_tmp;
  mp_mappedMemory = mp_arrayEvents_tmp;

  if (m_verbose >=2) Cout("vDataFile::Shuffle OK" << endl);
  #endif
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vDataFile::GetEventIndexStartAndStop( int64_t* ap_indexStart, int64_t* ap_indexStop, int a_subsetNum, int a_nbSubsets )
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_NORMAL)
  // Basically here, the index start for a single MPI instance will be the current subset number.
  // If multiple instances are used, the whole datafile is split into equal-size concatenated pieces.
  // So for each instance, we just have to found the first index falling in its range (assuming that
  // the event index step is always equal to the number of subsets).

  // For the first machine, index start is the current subset number
  if (mp_ID->GetMPIRank()==0) *ap_indexStart = a_subsetNum;
  // For the other machines, we must search for the first index falling in their range belonging to this
  // subset (a subset being starting at a_subsetNum with a step equal to the number of subsets, a_nbSubsets)
  else
  {
    // Compute the modulo of the first index of this machine minus the subset number with respect to the number of subsets
    int64_t modulo = (m_mpi1stEvent-a_subsetNum) % a_nbSubsets;
    // If this modulo is null, then the index start is the first index
    if (modulo==0) *ap_indexStart = m_mpi1stEvent;
    // Otherwise, the index start is equal to the first index plus the number of subsets minus the modulo
    else *ap_indexStart = m_mpi1stEvent + (a_nbSubsets - modulo);
  }

  // For index stop, we simply get the last event of the MPI instance (+1 because the for loop is exclusive)
  *ap_indexStop = m_mpiLastEvent + 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::GetMaxRingDiff()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_NORMAL)
  Cerr("*****vDataFile::GetMaxRingDiff() -> This function is not implemented for the used system" << endl);
  Cerr("                                   (this error may be prompted if the present function is erroneously called for a SPECT system)" << endl);
  return -1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::PROJ_WriteData()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  
  // Close all multithreaded datafiles before merging
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
    if (m2p_dataFile[th]) m2p_dataFile[th]->close();

  // If the output file doesn't exist yet, simply rename the first temporary file name using the ouput file name
  if (!ifstream(m_dataFileName))
  {
    // If only one projection is required (i.e one threads && no projection of frame or gates)...
    if(mp_ID->GetNbThreadsForProjection()==1  // No multithreading so no multiple tmp datafiles 
       && mp_ID->GetNbTimeFrames()
       *  mp_ID->GetNb1stMotImgsForLMS(mp_ID->GetNbTimeFrames()-1)
       *  mp_ID->GetNb2ndMotImgsForLMS()  == 1)
       {
         // rename the first temporary file name to the global name
         string tmp_file_name = m_dataFileName + "_0";
    
         // ... then just rename the first temporary file name using the ouput file name
         rename(tmp_file_name.c_str(),m_dataFileName.c_str());
         // no need to concatenate, so we leave here.
         return 0 ;
       }
  }

  // Create the final output file which will concatenate the events inside the thread-specific data files
  ofstream merged_file(m_dataFileName.c_str(), ios::out | ios::binary | ios::app);

  // Concatenation : generate input file ("data_file") to read the buffer of the thread-specific data files and store the information in the final output file
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    // Build thread file name
    stringstream ss; ss << th;
    string file_name = m_dataFileName;
    file_name.append("_").append(ss.str());
    // Open it
    // Note SS: Maybe we can use the m2p_dataFile[th] here by just rewarding to the beginning of the file ?
    // Note TM: There were some issues involving the use of rdbuf and concatenation of the file in this case (ifstream were needed instead of fstream for some reasons)
    //          But we should have another look on how the projection data writing works with the implementation of the new analytical simulator.
    ifstream data_file(file_name.c_str(), ios::binary | ios::in);

    if (!data_file)
    {
      Cerr(endl);
      Cerr("***** vDataFile::PROJ_WriteData() -> Input temporary thread file '" << file_name << "' is missing or corrupted !" << endl);
      return 1;
    }

    // Concatenate it to the merged file
    merged_file << data_file.rdbuf();
    // Close file
    data_file.close();
    
    // Re-open datafiles (needed if projecting frame/gates, as the contents of the temporary datafile are copied to the main datafile after each frame/gate)
    m2p_dataFile[th]->open( file_name.c_str(), ios::binary | ios::out | ios::trunc);
  }

  // Close merged file
  merged_file.close();

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vDataFile::PROJ_DeleteTmpDataFile()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("vDataFile::PROJ_DeleteTmpDataFile() ..." << endl);
  // Generate input file ("data_file") to read the buffer of the thread-specific data files and store the information in the final output file
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    // Build thread file name
    stringstream ss; ss << th;

    if (m2p_dataFile[th]) m2p_dataFile[th]->close();

    string file_name = m_dataFileName;
    file_name.append("_").append(ss.str());

    // Remove temporary file for data output writing (Projection script only)
    ifstream fcheck(file_name.c_str());
    if(fcheck.good())
    {
      fcheck.close();
      #ifdef _WIN32
      string dos_instruction = "del " + file_name;
      system(dos_instruction.c_str());
      #else
      remove(file_name.c_str());
      #endif
    }
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vEvent* vDataFile::PROJ_GenerateEvent(int a_idxElt1, int a_idxElt2, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  // Only 1 line required for projection/sensitivity generation
  m2p_BufferEvent[a_th]->SetNbLines(1);  
  m2p_BufferEvent[a_th]->SetID1(0, a_idxElt1);
  m2p_BufferEvent[a_th]->SetID2(0, a_idxElt2);
  // Return the event
  return m2p_BufferEvent[a_th];
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

string vDataFile::GetDataTypeToString()
{
  // Simply switch between the different data types
  string data_type = "";
  if (m_dataType==TYPE_CT) data_type = "CT";
  else if (m_dataType==TYPE_PET) data_type = "PET";
  else if (m_dataType==TYPE_SPECT) data_type = "SPECT";
  else data_type = "unknown";
  return data_type;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

string vDataFile::GetDataModeToString()
{
  // Simply switch between the different data modes
  string data_mode = "";
  if (m_dataMode==MODE_HISTOGRAM) data_mode = "histogram";
  else if (m_dataMode==MODE_LIST) data_mode = "list-mode";
  else if (m_dataMode==MODE_NORMALIZATION) data_mode = "normalization";
  else data_mode = "unknown";
  return data_mode;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

string vDataFile::GetDataSpecToString()
{
  // Simply switch between the different data specs
  string data_spec = "";
  if (m_dataSpec==SPEC_EMISSION) data_spec = "emission";
  else if (m_dataSpec==SPEC_TRANSMISSION) data_spec = "transmission";
  else data_spec = "unknown";
  return data_spec;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
