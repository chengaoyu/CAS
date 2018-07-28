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

  \brief Implementation of class iDataFileSPECT
*/

#include "iDataFileSPECT.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iDataFileSPECT::iDataFileSPECT() : vDataFile() 
{
  // Set all members to default values
  m_dataType = TYPE_SPECT;
  m_dataSpec = SPEC_EMISSION;
  m_isotope = "unknown";
  mp_nbOfBins[0] = 1;
  mp_nbOfBins[1] = 1;
  m_acquisitionZoom = 1.;
  m_nbOfProjections = 0;
  mp_angles = NULL;
  m_nbHeads = 1;
  mp_CORtoDetectorDistance = NULL;
  m_eventKindFlag = false;
  m_normCorrectionFlag = false;
  m_ignoreNormCorrectionFlag = false;
  m_scatCorrectionFlag = false;
  m_ignoreScatCorrectionFlag = false;
  m_headRotDirection = GEO_ROT_CW;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iDataFileSPECT::~iDataFileSPECT() 
{
  if (mp_angles) delete[] mp_angles;
  if (mp_CORtoDetectorDistance) delete[] mp_CORtoDetectorDistance;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::ReadSpecificInfoInHeader(bool a_affectQuantificationFlag)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("iDataFileSPECT::ReadSpecificInfoInHeader() -> Read information specific to SPECT" << endl);

  // Create pointers dedicated to recover the addresses of the member variables of the scanner object
  FLTNB* p_angles = NULL;
  FLTNB* p_CORtoDetectorDistance = NULL; 
  FLTNB p_pixSizeXY[2]; // not used here
  
  // Get Geometric parameters recovered from the scanner object
  sScannerManager* p_scannermanager;
  p_scannermanager = sScannerManager::GetInstance(); 
  p_scannermanager->GetSPECTSpecificParameters(&m_nbOfProjections, 
                                               &m_nbHeads,
                                               &m_acquisitionZoom,
                                                mp_nbOfBins,
                                                p_pixSizeXY, 
                                                p_angles, 
                                                p_CORtoDetectorDistance,
                                               &m_headRotDirection);

  // Check m_nbOfProjections first before allocating projection angles and radius using this variable
  if (m_nbOfProjections==0)
  {
    Cerr("***** iDataFileSPECT::ReadSpecificInfoInHeader() -> Number of projections should be strictly positive !" << endl);
    return 1;
  }

  // Allocation and initialization of Projection angles and Center of rotation to SPECT detector surface distance (radius) :
  mp_angles = new FLTNB[m_nbOfProjections];
  mp_CORtoDetectorDistance = new FLTNB[m_nbOfProjections];

  // Recover values
  for (int a=0 ; a<m_nbOfProjections ; a++)
  {
    mp_angles[a] = p_angles[a];
    mp_CORtoDetectorDistance[a] = p_CORtoDetectorDistance[a];
  }

  // Feedback to user
  if (m_verbose==VERBOSE_DETAIL) 
  {
    Cout("  --> Provided projection angles | distance Center of Rotation (COR) to detector: " << endl);
    for (int a=0 ; a<m_nbOfProjections ; a++) Cout("      " << mp_angles[a] << " | "  << mp_CORtoDetectorDistance[a] << endl);
  }

  // Read optional fields in the header, check if errors (issue during data reading/conversion (==1) )
  if (ReadDataASCIIFile(m_headerFileName, "Isotope", &m_isotope, 1, 0) == 1 ||
      ReadDataASCIIFile(m_headerFileName, "Event kind flag", &m_eventKindFlag, 1, 0) == 1 ||
      ReadDataASCIIFile(m_headerFileName, "Normalization correction flag", &m_normCorrectionFlag, 1, 0) == 1 ||
      ReadDataASCIIFile(m_headerFileName, "Scatter correction flag", &m_scatCorrectionFlag, 1, 0) == 1 )
      {
        Cerr("***** iDataFileSPECT::ReadSpecificInfoInHeader() -> Error while reading optional fields in the header data file !" << endl);
        return 1;
      }

  // Give the SPECT isotope to the oImageDimensionsAndQuantification that manages the quantification factors
  if (a_affectQuantificationFlag && mp_ID->SetSPECTIsotope(m_bedIndex, m_isotope))
  {
    Cerr("***** iDataFileSPECT::ReadSpecificInfoInHeader() -> A problem occured while setting the isotope to oImageDimensionsAndQuantification !" << endl);
    return 1;
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::SetSpecificParametersFrom(vDataFile* ap_DataFile)
{
  iDataFileSPECT* p_DataFileSPECT = (dynamic_cast<iDataFileSPECT*>(ap_DataFile));
  m_isotope = p_DataFileSPECT->GetIsotope();
  mp_nbOfBins[0] = p_DataFileSPECT->GetNbBins(0);
  mp_nbOfBins[1] = p_DataFileSPECT->GetNbBins(1);
  m_nbOfProjections = p_DataFileSPECT->GetNbProjections();
  mp_angles = p_DataFileSPECT->GetAngles();
  m_nbHeads = p_DataFileSPECT->GetNbHeads();
  mp_CORtoDetectorDistance = p_DataFileSPECT->GetCORtoDetectorDistance();
  m_eventKindFlag = p_DataFileSPECT->GetEventKindFlag();
  m_normCorrectionFlag = p_DataFileSPECT->GetNormCorrectionFlag();
  m_scatCorrectionFlag = p_DataFileSPECT->GetScatCorrectionFlag();
  m_headRotDirection = p_DataFileSPECT->GetHeadRotDirection();
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::ComputeSizeEvent()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("iDataFileSPECT::ComputeSizeEvent() -> In bytes" << endl);

  // For MODE_LIST events
  if (m_dataMode == MODE_LIST) 
  {
    // Size of the mandatory element in a list-mode event: Time + 2*crystalID
    m_sizeEvent = sizeof(uint32_t) + 2*sizeof(uint32_t);
    // Optional flags
    if (m_eventKindFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_scatCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_normCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    for (int i=0 ; i<3; i++)
    {
      // POI available for this direction
      if (mp_POIDirectionFlag[i]) m_sizeEvent += sizeof(FLTNBDATA);
    }
  }
  // For MODE_HISTOGRAM events
  if(m_dataMode == MODE_HISTOGRAM)
  {
    // Size of the mandatory element in a histo event: Time + event_value + 2*crystalID
    m_sizeEvent = sizeof(uint32_t) + sizeof(FLTNBDATA) + 2*sizeof(uint32_t);
    // Optional flags
    if (m_scatCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_normCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
  }
  // Unknown event type
  else
  {
    Cerr("***** iDataFileSPECT::ComputeSizeEvent() -> Unknown event mode !" << endl);
    return 1;
  }

  // Check
  if (m_sizeEvent<=0) 
  {
    Cerr("***** iDataFileSPECT::ComputeSizeEvent() -> Error, the Event size in bytes should be >= 0 !" << endl;);
    return 1;
  }

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Event size: " << m_sizeEvent << " bytes" << endl);
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::PrepareDataFile()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_NORMAL)
  {
    if (m_dataMode==MODE_HISTOGRAM) Cout("iDataFileSPECT::PrepareDataFile() -> Build histogram events" << endl);
    else if (m_dataMode==MODE_LIST) Cout("iDataFileSPECT::PrepareDataFile() -> Build listmode events" << endl);
    else if (m_dataMode==MODE_NORMALIZATION) Cout("iDataFileSPECT::PrepareDataFile() -> Build normalization events" << endl);
  }

  // ==============================================================================
  // Allocate event buffers (one for each thread)
  // ==============================================================================

  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Allocate an event buffer for each thread" << endl);
  // Instanciation of the event buffer according to the data type
  m2p_BufferEvent = new vEvent*[mp_ID->GetNbThreadsForProjection()];
  
  // Allocate the events per each thread
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    // For MODE_LIST events
    if (m_dataMode == MODE_LIST) 
    {
      m2p_BufferEvent[th] = new iEventListSPECT();
    }
    // For MODE_HISTOGRAM events
    if (m_dataMode == MODE_HISTOGRAM)
    {
      m2p_BufferEvent[th] = new iEventHistoSPECT();
    }
    // Allocate crystal/view IDs
    if (m2p_BufferEvent[th]->AllocateID())
    {
      Cerr("*****iDataFileSPECT::PrepareDataFile() -> Error while trying to allocate memory for the Event object!" << endl);
      return 1;
    }
  }

  // ==============================================================================
  // Deal with specific corrections
  // ==============================================================================

  // In case of normalization correction flag, see if we ignore this correction
  if (m_normCorrectionFlag)
  {
    // Affect the ignored flag from the ignored corrections list processed by the oImageDimensionsAndQuantification
    m_ignoreNormCorrectionFlag = mp_ID->GetIgnoreNormCorrectionFlag();
    // Verbose
    if (m_verbose>=VERBOSE_DETAIL)
    {
      if (m_ignoreNormCorrectionFlag) Cout("  --> Ignore normalization correction" << endl);
      else Cout("  --> Correct for normalization" << endl);
    }
  }
  // In case of scatter correction flag, see if we ignore this correction
  if (m_scatCorrectionFlag)
  {
    // Affect the ignored flag from the ignored corrections list processed by the oImageDimensionsAndQuantification
    m_ignoreScatCorrectionFlag = mp_ID->GetIgnoreScatCorrectionFlag();
    // Verbose
    if (m_verbose>=VERBOSE_DETAIL)
    {
      if (m_ignoreScatCorrectionFlag) Cout("  --> Ignore scatter correction" << endl);
      else Cout("  --> Correct for scatter events" << endl);
    }
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vEvent* iDataFileSPECT::GetEventSpecific(char* ap_buffer, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // Work on a copy of the input pointer
  char* file_position = ap_buffer;

  // For MODE_LIST SPECT data
  if (m_dataMode == MODE_LIST)
  {
    // Cast the event pointer
    iEventListSPECT* event = (dynamic_cast<iEventListSPECT*>(m2p_BufferEvent[a_th]));
    // Mandatory time field: [uint32_t (time)]
    event->SetTimeInMs(*reinterpret_cast<uint32_t*>(file_position)); 
    file_position += sizeof(uint32_t);
    // Optional kind: [uint8_t kind]
    if (m_eventKindFlag)
    {
      event->SetKind(*reinterpret_cast<uint8_t*>(file_position));
      file_position += sizeof(uint8_t);
    }
    // Optional scatter correction field: [FLTNBDATA (scatter)]
    if (m_scatCorrectionFlag)
    {
      if (!m_ignoreScatCorrectionFlag) event->SetScatterRate(*reinterpret_cast<FLTNBDATA*>(file_position)); 
      file_position += sizeof(FLTNBDATA);
    }
    // Optional normalization correction field: [FLTNBDATA (norm)]
    if (m_normCorrectionFlag)
    {
      if (!m_ignoreNormCorrectionFlag) event->SetNormalizationFactor(*reinterpret_cast<FLTNBDATA*>(file_position)); 
      file_position += sizeof(FLTNBDATA);
    }
    // Optional POI correction field: [FLTNBDATA (POI1[3])]
    for (int i=0; i<3; i++)
    {
      if (mp_POIDirectionFlag[i])
      {
        if (!m_ignorePOIFlag) event->SetPOI(i,*reinterpret_cast<FLTNBDATA*>(file_position)); 
        file_position += sizeof(FLTNBDATA);
      }
    }
    // Mandatory angular projection ID: [uint32_t (ID1)]
    event->SetID1(0, *reinterpret_cast<uint32_t*>(file_position)); 
    file_position += sizeof(uint32_t);
    // Mandatory crystal ID: [uint32_t (ID2)]
    event->SetID2(0, *reinterpret_cast<uint32_t*>(file_position)); 
    file_position += sizeof(uint32_t);
  }

  // For MODE_HISTOGRAM SPECT DATA
  if (m_dataMode == MODE_HISTOGRAM) 
  {
    // Cast the event pointer
    iEventHistoSPECT* event = (dynamic_cast<iEventHistoSPECT*>(m2p_BufferEvent[a_th]));
    // Mandatory time field: [uint32_t (time)]
    event->SetTimeInMs(*reinterpret_cast<uint32_t*>(file_position)); 
    file_position += sizeof(uint32_t);
    // Mandatory bin value: [FLTNBDATA bin value]
    event->SetEventValue(0, *reinterpret_cast<FLTNBDATA*>(file_position));
    file_position += sizeof(FLTNBDATA);
    // Optional scatter correction field: [FLTNBDATA (scatter)]
    if (m_scatCorrectionFlag)
    {
      if (!m_ignoreScatCorrectionFlag) event->SetScatterRate(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }
    // Optional normalization correction field: [FLTNBDATA (norm)]
    if (m_normCorrectionFlag)
    {
      if (!m_ignoreNormCorrectionFlag) event->SetNormalizationFactor(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }
    // Mandatory angular projection ID: [uint32_t (ID1)]
    event->SetID1(0, *reinterpret_cast<uint32_t*>(file_position)); 
    file_position += sizeof(uint32_t);
    // Mandatory crystal ID: [uint32_t (ID2)]
    event->SetID2(0, *reinterpret_cast<uint32_t*>(file_position)); 
    file_position += sizeof(uint32_t);
  }

  // Return the updated event
  return m2p_BufferEvent[a_th];
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iDataFileSPECT::DescribeSpecific()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose==0) return;
  // Describe the datafile
  Cout("iDataFileSPECT::DescribeSpecific() -> Here is some specific content of the SPECT datafile" << endl);
  Cout("  --> Isotope: " << m_isotope << endl);
  if (m_dataMode==MODE_LIST && m_eventKindFlag) Cout("  --> Event kind is present" << endl);
  if (m_scatCorrectionFlag) Cout("  --> Scatter correction is present" << endl);
  if (m_normCorrectionFlag) Cout("  --> Normalization correction is present" << endl);
  Cout("  --> Number of heads: " << m_nbHeads << endl);
  if (mp_nbOfBins[0]!=1) Cout("  --> Transaxial number of bins: " << mp_nbOfBins[0] << endl);
  if (mp_nbOfBins[1]!=1) Cout("  --> Axial number of bins: " << mp_nbOfBins[1] << endl);
  if (m_acquisitionZoom!=1.) Cout("  --> Acquisition zoom: " << m_acquisitionZoom << endl);
  if (m_headRotDirection==GEO_ROT_CW) Cout("  --> Head rotation is clockwise" << endl);
  else if (m_headRotDirection==GEO_ROT_CCW) Cout("  --> Head rotation is counter-clockwise" << endl);
  else Cout("  --> Head rotation is undefined !!!" << endl);
  Cout("  --> Number of acquisition projections: " << m_nbOfProjections << endl);
  for (uint16_t p=0; p<m_nbOfProjections; p++) Cout("    | Projection " << p << " at " << mp_angles[p] << " deg and "
                                               << mp_CORtoDetectorDistance[p] << " mm from center of rotation" << endl);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::CheckSpecificParameters()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Error if m_dataType != SPECT
  if (m_dataType != TYPE_SPECT)
  {
    Cerr("***** iDataFileSPECT::CheckSpecificParameters() -> Data type should be SPECT !'" << endl);
    return 1;
  }
  // Check number of projections
  if (m_nbOfProjections == 0)
  {
    Cerr("***** iDataFileSPECT::CheckSpecificParameters() -> Number of projections not initialized (should be >0) !" << endl);
    return 1;
  }
  // Check projection angles
  if (mp_angles == NULL)
  {
    Cerr("***** iDataFileSPECT::CheckSpecificParameters() -> Projection angles not initialized !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::CheckFileSizeConsistency()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Create the stream
  fstream* p_file = new fstream( m_dataFileName.c_str(), ios::binary| ios::in );
  // Check that datafile exists
  if (!p_file->is_open()) 
  {
    Cerr("***** iDataFilePET::CheckFileSizeConsistency() -> Failed to open input file '" << m_dataFileName.c_str() << "' !" << endl);
    Cerr("                                                 (Provided in the data file header: " << m_headerFileName << ")" << endl);
    return 1;
  }
  // Get file size in bytes
  p_file->seekg(0, ios::end);
  int64_t sizeInBytes = p_file->tellg();
  // Close stream and delete it
  p_file->close();
  delete p_file;
  // Check datafile self-consistency
  if (m_nbEvents*m_sizeEvent != sizeInBytes)
  {
    Cerr("----------------------------------------------------------------------------------------------------------------------------------------" << endl);
    Cerr("***** iDataFileSPECT::CheckFileSizeConsistency() -> DataFile size is not consistent with the information provided by the user/datafile !" << endl);
    Cerr("  --> Expected size: "<< m_nbEvents*m_sizeEvent << endl);
    Cerr("  --> Actual size: "<< sizeInBytes << endl << endl);
    Cerr("      ADDITIONAL INFORMATION ABOUT THE DATAFILE INITIALIZATION" << endl);
    if (m_eventKindFlag) Cerr("  --> Event kind term is enabled" << endl);
    else Cerr("  --> No information about the kind of events in the data" << endl);
    if (m_normCorrectionFlag) Cerr("  --> Normalization correction term is enabled" << endl);
    else Cerr("  --> No normalization correction term in the data" << endl);
    if (m_scatCorrectionFlag) Cerr("  --> Scatter correction term is enabled" << endl);
    else  Cerr("  --> No scatter correction term in the data" << endl);
    Cerr("  --> Calibration factor value is: " << m_calibrationFactor << endl);
    Cerr("  --> Number of bins are equal to " << mp_nbOfBins[0] << "," << mp_nbOfBins[1] << " for transaxial,axial axis respectively" << endl);
    if (mp_POIResolution[0]<0.) Cerr("  --> No POI enabled on the radial axis" << endl);
    else Cerr("  --> POI resolution on the radial axis is: " << mp_POIResolution[0] << endl);
    if (mp_POIResolution[1]<0.) Cerr("  --> No POI enabled on the tangential axis" << endl);
    else Cerr("  --> POI resolution on the tangential axis is: " << mp_POIResolution[1] << endl);
    if (mp_POIResolution[2]<0.) Cerr("  --> No POI enabled on the axial axis" << endl);
    else Cerr("  --> POI resolution on the axial axis is: " << mp_POIResolution[2] << endl);
    Cerr("----------------------------------------------------------------------------------------------------------------------------------------" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::CheckSpecificConsistencyWithAnotherDataFile(vDataFile* ap_DataFile)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Dynamic cast the vDataFile to a iDataFilePET
  iDataFileSPECT* p_data_file = (dynamic_cast<iDataFileSPECT*>(ap_DataFile));
  // Check isotope
  if (m_isotope!=p_data_file->GetIsotope())
  {
    Cerr("***** iDataFileSPECT::CheckSpecificConsistencyWithAnotherDataFile() -> Isotopes are inconsistent !" << endl);
    return 1;
  }
  // Check event kind flag
  if (m_eventKindFlag!=p_data_file->GetEventKindFlag())
  {
    Cerr("***** iDataFileSPECT::CheckSpecificConsistencyWithAnotherDataFile() -> Event kind flags are inconsistent !" << endl);
    return 1;
  }
  // Check normalization correction flag
  if (m_normCorrectionFlag!=p_data_file->GetNormCorrectionFlag())
  {
    Cerr("***** iDataFileSPECT::CheckSpecificConsistencyWithAnotherDataFile() -> Normalization correction flags are inconsistent !" << endl);
    return 1;
  }
  // Check scatter correction flag
  if (m_scatCorrectionFlag!=p_data_file->GetScatCorrectionFlag())
  {
    Cerr("***** iDataFileSPECT::CheckSpecificConsistencyWithAnotherDataFile() -> Scatter correction flags are inconsistent !" << endl);
    return 1;
  }
  // Check data mode
  if (m_dataMode!=p_data_file->GetDataMode())
  {
    Cerr("***** iDataFileSPECT::CheckSpecificConsistencyWithAnotherDataFile() -> Data modes are inconsistent (list-mode or histogram) !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::InitAngles(FLTNB* ap_angles)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Check
  if (m_nbOfProjections == 0)
  {
    Cerr("***** iDataFileSPECT::InitAngles() -> Number of projection angles not initialized !'" << endl);
    return 1;
  }
  // Allocate
  mp_angles = new FLTNB[m_nbOfProjections];
  // Affect
  for (uint16_t a=0 ; a<m_nbOfProjections ; a++) mp_angles[a] = ap_angles[a];
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::InitCorToDetectorDistance(FLTNB* ap_CORtoDetectorDistance)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Check
  if (m_nbOfProjections == 0)
  {
    Cerr("***** iDataFileSPECT::InitAngles() -> Number of projection angles not initialized !'" << endl);
    return 1;
  }
  // Allocate
  mp_CORtoDetectorDistance = new FLTNB[m_nbOfProjections];
  // Affect
  for (uint16_t a=0 ; a<m_nbOfProjections ; a++) mp_CORtoDetectorDistance[a] = ap_CORtoDetectorDistance[a];
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::PROJ_InitFile()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("iDataFileSPECT::PROJ_InitFile() -> Initialize datafile for writing" << endl);
  
  m_startTimeInSec = 0.; //Std initialization for projection
  m_durationInSec = 1.; //Std initialization for projection
  m_nbEvents = 0; //Std initialization for projection
  m_calibrationFactor =  1.;
  m_isotope = "unknown";
  m_scatCorrectionFlag = false;
  m_normCorrectionFlag = false;
  mp_POIResolution[0] = -1.;
  mp_POIResolution[1] = -1.;
  mp_POIResolution[2] = -1.;

  // Instanciate a fstream datafile for each thread
  m2p_dataFile = new fstream*[mp_ID->GetNbThreadsForProjection()];

  // Set name of the projection data header
  sOutputManager* p_OutMgr;
  p_OutMgr = sOutputManager::GetInstance();
  string path_name = p_OutMgr->GetPathName();
  string img_name = p_OutMgr->GetBaseName();
  m_headerFileName = path_name.append(img_name).append("_CstrProj").append(".Cdh");
  
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    m_dataFileName = m_headerFileName.substr(0, m_headerFileName.find_last_of(".")).append(".Cdf"); // Generate datafile name from header file

    // Projeted data will be written in several files (corresponding to the number of thread) to be concatenated at the end of the projection process.   
    stringstream ss;
    ss << th;

    string datafile_name = m_dataFileName;
    datafile_name.append("_").append(ss.str());
    
    m2p_dataFile[th] = new fstream( datafile_name.c_str(), ios::binary | ios::out | ios::trunc);
  }

  //remove content from the output data file, in case it already exists
  //todo warn the user a datafile with the same name will be erased (eventually add an option to disable the warning)
  ofstream output_file(m_dataFileName.c_str(), ios::out | ios::trunc);
  output_file.close();
  
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Output datafile name :" << m_dataFileName << endl); 
    
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::WriteEvent(vEvent* ap_Event, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  if (m_dataMode == MODE_LIST) 
  {
    // TODO should create as many vEvent as (int)fp.
    if (WriteListEvent((iEventListSPECT*)ap_Event, a_th))
    {
      Cerr("*****iDataFileSPECT::WriteEvent() -> Error while trying to write projection datafile (list-mode)" << endl);
      return 1;
    }
  }

  if (m_dataMode == MODE_HISTOGRAM)
  {
    if (WriteHistoEvent((iEventHistoSPECT*)ap_Event, a_th))
    {
      Cerr("*****iDataFileSPECT::WriteEvent() -> Error while trying to write projection datafile (histogram)" << endl);
      return 1;
    }
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::WriteHistoEvent(iEventHistoSPECT* ap_Event, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  // Write sequentially each field of the event according to the type of the event.
  m2p_dataFile[a_th]->clear();
  
  uint32_t time = ap_Event->GetTimeInMs();
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&time), sizeof(uint32_t));
  
  FLTNB event_value = ap_Event->GetEventValue(0);
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&event_value), sizeof(FLTNB));

  if(m_scatCorrectionFlag)   
  {
    FLTNB scat_rate = ap_Event->GetEventScatRate();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&scat_rate), sizeof(FLTNB));
  }
  
  if(m_normCorrectionFlag)   
  {
    FLTNB norm_corr_factor = ap_Event->GetNormFactor();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&norm_corr_factor), sizeof(FLTNB));
  }

  uint32_t id1 = ap_Event->GetID1(0);
  uint32_t id2 = ap_Event->GetID2(0);
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&id1), sizeof(uint32_t));
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&id2), sizeof(uint32_t));

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::WriteListEvent(iEventListSPECT* ap_Event, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  // Write sequentially each field of the event according to the type of the event.
  m2p_dataFile[a_th]->clear();
  
  uint32_t time = ap_Event->GetTimeInMs();
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&time), sizeof(uint32_t));

  if(m_eventKindFlag) 
  {
    uint8_t event_kind = ap_Event->GetKind();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&event_kind), sizeof(uint8_t));
  }
  
  if(m_scatCorrectionFlag)   
  {
    FLTNB scat_rate = ap_Event->GetEventScatRate();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&scat_rate), sizeof(FLTNB));
  }
  
  if(m_normCorrectionFlag)   
  {
    FLTNB norm_corr_factor = ap_Event->GetNormFactor();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&norm_corr_factor), sizeof(FLTNB));
  }
  
  if(mp_POIResolution[0]>0)
  {
    FLTNB POI = ap_Event->GetPOI(0);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI), sizeof(FLTNB));
  }

  if(mp_POIResolution[1]>0)
  {
    FLTNB POI = ap_Event->GetPOI(1);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI), sizeof(FLTNB));
  }

  if(mp_POIResolution[2]>0)
  {
    FLTNB POI = ap_Event->GetPOI(2);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI), sizeof(FLTNB));
  }

  uint32_t id1 = ap_Event->GetID1(0);
  uint32_t id2 = ap_Event->GetID2(0);
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&id1), sizeof(uint32_t));
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&id2), sizeof(uint32_t));

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::WriteHeader()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_NORMAL) Cout("iDataFileSPECT::WriteHeader() -> Write header file '" << m_headerFileName << "'" << endl); 
  // Open file
  fstream headerFile;
  headerFile.open(m_headerFileName.c_str(), ios::out);
  if (!headerFile.is_open())
  {
    Cerr("***** iDataFileSPECT::WriteHeader() -> Failed to open output header file '" << m_headerFileName << "' !" << endl);
    return 1;
  }
  // Data file name
  headerFile << "Data filename: " << GetFileFromPath(m_dataFileName) << endl;
  // Number of events
  headerFile << "Number of events: " << m_nbEvents << endl; 
  // Data mode
  if (m_dataMode==MODE_HISTOGRAM) headerFile << "Data mode: histogram" << endl;
  else if (m_dataMode==MODE_LIST) headerFile << "Data mode: list-mode" << endl;
  else if (m_dataMode==MODE_NORMALIZATION) headerFile << "Data mode: normalization" << endl;
  // SPECT data type
  headerFile << "Data type: SPECT" << endl;
  // Acquisition start time in seconds
  headerFile << "Start time (s): " << m_startTimeInSec << endl; 
  // Acquisition duration in seconds
  headerFile << "Duration (s): " << m_durationInSec << endl; 
  // Scanner name
  headerFile << "Scanner name: " << sScannerManager::GetInstance()->GetScannerName() << endl; 
  // Number of bins
  if (mp_nbOfBins[0] != 0 && mp_nbOfBins[1] != 0)
    headerFile << "Number of bins: " << mp_nbOfBins[0] << ", " << mp_nbOfBins[1] << endl; 
  // Number of projections
  headerFile << "Number of projections: " << m_nbOfProjections << endl; 
  // Projection angles
  headerFile << "Projection angles: " << mp_angles[0];
  for (int a=1 ; a<m_nbOfProjections ; a++) headerFile << ", " << mp_angles[a]; 
  headerFile << endl;
  // Distance camera surface to COR
  headerFile << "Distance camera surface to COR: " << mp_CORtoDetectorDistance[0];
  for (int h=1 ; h<m_nbHeads ; h++) headerFile << ", " << mp_CORtoDetectorDistance[h] << endl;
  headerFile << endl;
  // Calibration factor
  headerFile << "Calibration factor: " << m_calibrationFactor << endl;
  // Isotope
  headerFile << "Isotope: " << m_isotope << endl;
  // Normalization correction flag
  headerFile << "Normalization correction flag: " << m_normCorrectionFlag << endl;
  // Scatter correction flag
  headerFile << "Scatter correction flag: " << m_scatCorrectionFlag << endl;
  // Head rotation direction
  string rot_direction = (m_headRotDirection == GEO_ROT_CCW) ? "CCW" : "CW";
  headerFile << "Head rotation direction: " << rot_direction << endl; 
  // Close file
  headerFile.close();
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFileSPECT::PROJ_GetScannerSpecificParameters()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("iDataFileSPECT::PROJ_GetScannerSpecificParameters() ..."<< endl);
  // Create pointers dedicated to recover the addresses of the member variables of the scanner object
  FLTNB* p_angles = NULL;
  FLTNB* p_CORtoDetectorDistance = NULL; 
  FLTNB p_pixSizeXY[2] ; // not used here
  // Get pointers to SPECT specific parameters in scanner
  if (sScannerManager::GetInstance()->GetSPECTSpecificParameters(&m_nbOfProjections, &m_nbHeads, &m_acquisitionZoom, mp_nbOfBins, p_pixSizeXY, p_angles, p_CORtoDetectorDistance, &m_headRotDirection) )
  {
    Cerr("*****iDataFileSPECT::PROJ_GetScannerSpecificParameters() -> An error occurred while trying to get SPECT geometric parameters from the scanner object !" << endl);
    return 1;
  }
  // In projection, mp_CORtoDetectorDistance is allocated according to the number of heads, and not the number of projections
  // (Projection does not currently allow to set a radius specific to the projection angles)
  mp_CORtoDetectorDistance = new FLTNB[m_nbHeads];
  // Retrieve SPECT distance between the detector and center of rotation
  for (int h=0 ; h<m_nbHeads ; h++)
    mp_CORtoDetectorDistance[h] = p_CORtoDetectorDistance[0];
  // Retrieve projection angles
  mp_angles = new FLTNB[m_nbOfProjections];
  for (int a=0 ; a<m_nbOfProjections ; a++)
    mp_angles[a] = p_angles[a];
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
