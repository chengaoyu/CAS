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

  \brief Implementation of class iDataFilePET
*/

#include "iDataFilePET.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iDataFilePET::iDataFilePET() : vDataFile() 
{
  // Set all members to default values
  m_dataType = TYPE_PET;
  m_dataSpec = SPEC_EMISSION;
  m_maxNumberOfLinesPerEvent = 1;
  m_maxAxialDiffmm = -1.;
  m_isotope = "unknown";
  m_resolutionTOF = -1.;
  m_nbBinsTOF = -1;
  m_binSizeTOF = -1.;
  m_eventKindFlag = false;
  m_atnCorrectionFlag = false;
  m_ignoreAttnCorrectionFlag = false;
  m_normCorrectionFlag = false;
  m_ignoreNormCorrectionFlag = false;
  m_scatCorrectionFlag = false;
  m_ignoreScatCorrectionFlag = false;
  m_randCorrectionFlag = false;
  m_ignoreRandCorrectionFlag = false;
  m_TOFInfoFlag = false;
  m_ignoreTOFFlag = false;
  m_TOFMeasurementRange = 0.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iDataFilePET::~iDataFilePET() {;}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::ReadSpecificInfoInHeader(bool a_affectQuantificationFlag)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("iDataFilePET::ReadSpecificInfoInHeader() -> Read information specific to PET" << endl);
  // Read optional fields in the header, check if errors (issue during data reading/conversion (==1) )
  if (ReadDataASCIIFile(m_headerFileName, "Maximum number of lines per event", &m_maxNumberOfLinesPerEvent, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Maximum axial difference mm", &m_maxAxialDiffmm, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "TOF resolution (ps)", &m_resolutionTOF, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "TOF information flag", &m_TOFInfoFlag, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "List TOF measurement range (ps)", &m_TOFMeasurementRange, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Histo TOF number of bins", &m_nbBinsTOF, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Histo TOF bin size (ps)", &m_binSizeTOF, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Coincidence kind flag", &m_eventKindFlag, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Attenuation correction flag", &m_atnCorrectionFlag, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Normalization correction flag", &m_normCorrectionFlag, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Scatter correction flag", &m_scatCorrectionFlag, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Random correction flag", &m_randCorrectionFlag, 1, KEYWORD_OPTIONAL)==1 ||
      ReadDataASCIIFile(m_headerFileName, "Isotope", &m_isotope, 1, KEYWORD_OPTIONAL)==1 )
  {
    Cerr("***** iDataFilePET::ReadSpecificInfoInHeader() -> Error while reading optional fields in the header data file !" << endl);
    return 1;
  }
  // Give the PET isotope to the oImageDimensionsAndQuantification that manages the quantification factors
  if (a_affectQuantificationFlag && mp_ID->SetPETIsotope(m_bedIndex, m_isotope))
  {
    Cerr("***** iDataFilePET::ReadSpecificInfoInHeader() -> A problem occured while setting the isotope to oImageDimensionsAndQuantification !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::SetSpecificParametersFrom(vDataFile* ap_DataFile)
{
  iDataFilePET* p_DataFilePET = (dynamic_cast<iDataFilePET*>(ap_DataFile));
  m_maxNumberOfLinesPerEvent = p_DataFilePET->GetMaxNumberOfLinesPerEvent();
  m_maxAxialDiffmm = p_DataFilePET->GetMaxAxialDiffmm();
  m_isotope = p_DataFilePET->GetIsotope();
  m_resolutionTOF = p_DataFilePET->GetTOFResolution();
  m_nbBinsTOF = p_DataFilePET->GetNbTOFBins();
  m_binSizeTOF = p_DataFilePET->GetTOFBinSize();
  m_TOFInfoFlag = p_DataFilePET->GetTOFInfoFlag();
  m_eventKindFlag = p_DataFilePET->GetEventKindFlag();
  m_atnCorrectionFlag = p_DataFilePET->GetAtnCorrectionFlag();
  m_normCorrectionFlag = p_DataFilePET->GetNormCorrectionFlag();
  m_scatCorrectionFlag = p_DataFilePET->GetScatCorrectionFlag();
  m_randCorrectionFlag = p_DataFilePET->GetRandCorrectionFlag();
  m_TOFMeasurementRange = p_DataFilePET->GetTOFMeasurementRange();
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::ComputeSizeEvent()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("iDataFilePET::ComputeSizeEvent() -> In bytes" << endl);

  // For MODE_LIST events
  if (m_dataMode == MODE_LIST) 
  {
    // Size of the mandatory element in a list-mode event: Time + max_nb_lines*2*crystalID
    m_sizeEvent = sizeof(uint32_t) + m_maxNumberOfLinesPerEvent*2*sizeof(uint32_t);
    // Number of lines
    if (m_maxNumberOfLinesPerEvent>1) m_sizeEvent += sizeof(uint16_t);
    // Optional flags
    if (m_eventKindFlag) m_sizeEvent += sizeof(uint8_t);
    if (m_atnCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_scatCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_randCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_normCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_TOFInfoFlag) m_sizeEvent += sizeof(FLTNBDATA);
    // POI availability in each direction
    for (int i=0 ; i<3; i++) if (mp_POIDirectionFlag[i]) m_sizeEvent += 2*sizeof(FLTNBDATA);
  }
  // For MODE_HISTOGRAM events
  else if (m_dataMode == MODE_HISTOGRAM)
  {
    // Size of the mandatory element in a histo event: Time + nbBinsTOF*event_value + max_nb_line*2*crystalID
    m_sizeEvent = sizeof(uint32_t) + m_nbBinsTOF*sizeof(FLTNBDATA) + m_maxNumberOfLinesPerEvent*2*sizeof(uint32_t);
    // Number of lines
    if (m_maxNumberOfLinesPerEvent>1) m_sizeEvent += sizeof(uint16_t);
    // Optional flags
    if (m_atnCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_randCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_normCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    if (m_scatCorrectionFlag) m_sizeEvent += m_nbBinsTOF*sizeof(FLTNBDATA);
  }
  // For MODE_NORMALIZATION events
  else if (m_dataMode == MODE_NORMALIZATION)
  {
    // max_nb_lines*2*crystalID
    m_sizeEvent = m_maxNumberOfLinesPerEvent*2*sizeof(uint32_t);
    // Number of lines
    if (m_maxNumberOfLinesPerEvent>1) m_sizeEvent += sizeof(uint16_t);
    // Optional flag for normalization
    if (m_normCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
    // Optional flag for attenuation
    if (m_atnCorrectionFlag) m_sizeEvent += sizeof(FLTNBDATA);
  }
  // Unknown event type
  else
  {
    Cerr("***** iDataFilePET::ComputeSizeEvent() -> Unknown event mode !" << endl);
    return 1;
  }

  // Check
  if (m_sizeEvent<=0) 
  {
    Cerr("***** iDataFilePET::ComputeSizeEvent() -> Error, the Event size in bytes should be >= 0 !" << endl;);
    return 1;
  }

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Event size = " << m_sizeEvent << " bytes" << endl);
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::PrepareDataFile()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_NORMAL)
  {
    if (m_dataMode==MODE_HISTOGRAM) Cout("iDataFilePET::PrepareDataFile() -> Build histogram events" << endl);
    else if (m_dataMode==MODE_LIST) Cout("iDataFilePET::PrepareDataFile() -> Build listmode events" << endl);
    else if (m_dataMode==MODE_NORMALIZATION) Cout("iDataFilePET::PrepareDataFile() -> Build normalization events" << endl);
  }

  // ==============================================================================
  // Allocate event buffers (one for each thread)
  // ==============================================================================
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Allocating an event buffer for each thread" << endl);
  // Instanciation of the event buffer according to the data type
  m2p_BufferEvent = new vEvent*[mp_ID->GetNbThreadsForProjection()];
  
  // Allocate the events per each thread
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    
    // For MODE_LIST events
    if (m_dataMode == MODE_LIST) 
    {
      m2p_BufferEvent[th] = new iEventListPET();
      // If TOF, convert and transmit information about the total range of TOF measurements ( conversion from delta time into delta length )
      if (m_TOFInfoFlag && !m_ignoreTOFFlag) 
      {
        ((iEventListPET*)m2p_BufferEvent[th])->SetHasTOFinfo(true);
        ((iEventListPET*)m2p_BufferEvent[th])->SetSpatialTOFMeasurementRange(m_TOFMeasurementRange*SPEED_OF_LIGHT*0.5);
      }
    }
    // For MODE_HISTOGRAM events
    else if (m_dataMode == MODE_HISTOGRAM)
    {
      m2p_BufferEvent[th] = new iEventHistoPET();
      // If we ignore the TOF information, then the event buffer will have only one TOF bin;
      // the TOF contributions will be sum up when reading the event.
      if (m_ignoreTOFFlag) ((iEventHistoPET*)m2p_BufferEvent[th])->SetEventNbTOFBins(1);
      else ((iEventHistoPET*)m2p_BufferEvent[th])->SetEventNbTOFBins(m_nbBinsTOF);
    }
    // For MODE_NORMALIZATION events
    else if (m_dataMode == MODE_NORMALIZATION)
    {
      m2p_BufferEvent[th] = new iEventNorm();
    }

    // Set the maximum number of lines per event
    m2p_BufferEvent[th]->SetNbLines(m_maxNumberOfLinesPerEvent);
    
    // Allocate crystal IDs
    if (m2p_BufferEvent[th]->AllocateID())
    {
      Cerr("*****iDataFilePET::PrepareDataFile() -> Error while trying to allocate memory for the Event object for thread " << th << " !" << endl;);
      return 1;
    }
  }

  // ==============================================================================
  // Deal with specific corrections
  // ==============================================================================

  // In case of attenuation correction flag, see if we ignore this correction
  if (m_atnCorrectionFlag)
  {
    // Affect the ignored flag from the ignored corrections list processed by the oImageDimensionsAndQuantification
    m_ignoreAttnCorrectionFlag = mp_ID->GetIgnoreAttnCorrectionFlag();
    // Verbose
    if (m_verbose>=VERBOSE_DETAIL)
    {
      if (m_ignoreAttnCorrectionFlag) Cout("  --> Ignore attenuation correction" << endl);
      else Cout("  --> Correct for attenuation" << endl);
    }
  }
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
  // In case of random correction flag, see if we ignore this correction
  if (m_randCorrectionFlag)
  {
    // Affect the ignored flag from the ignored corrections list processed by the oImageDimensionsAndQuantification
    m_ignoreRandCorrectionFlag = mp_ID->GetIgnoreRandCorrectionFlag();
    // Verbose
    if (m_verbose>=VERBOSE_DETAIL)
    {
      if (m_ignoreRandCorrectionFlag) Cout("  --> Ignore random correction" << endl);
      else Cout("  --> Correct for random events" << endl);
    }
  }
  
  // ==============================================================================
  // Deal with TOF
  // ==============================================================================
  
  // TOF information availability and use
  if (m_TOFInfoFlag)
  {
    // Special case when TOF is ignored but we are in list-mode, we have scatter correction factors and we do not ignore them.
    // In this case, the scatter correction will be completely wrong, so we make the code crash here
    if (m_ignoreTOFFlag && m_dataMode==MODE_LIST && m_scatCorrectionFlag && !m_ignoreScatCorrectionFlag)
    {
      Cerr("***** iDataFilePET::PrepareDataFile() -> Scatter correction for list-mode TOF data will be erroneous if TOF information is ignored in the projections !" << endl);
      return 1;
    }
    // Verbose
    if (m_verbose>=VERBOSE_DETAIL)
    {
      if (m_ignoreTOFFlag) Cout("  --> TOF information available but ignored " << endl);
      else Cout("  --> Use TOF information" << endl);
    }
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vEvent* iDataFilePET::GetEventSpecific(char* ap_buffer, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // Work on a copy of the input pointer
  char* file_position = ap_buffer;

  // For MODE_LIST PET data
  if (m_dataMode == MODE_LIST)
  {
    // Cast the event pointer
    iEventListPET* event = (dynamic_cast<iEventListPET*>(m2p_BufferEvent[a_th]));
    // Mandatory time field: [uint32_t (time)]
    event->SetTimeInMs(*reinterpret_cast<uint32_t*>(file_position)); 
    file_position += sizeof(uint32_t);
    // Optional attenuation correction field: [FLTNBDATA (attnCorrectionFactor)]
    if (m_atnCorrectionFlag)
    {
      if (!m_ignoreAttnCorrectionFlag) event->SetAttenuationCorrectionFactor(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }
    // Optional kind: [uint8_t kind]
    if (m_eventKindFlag)
    {
      event->SetKind(*reinterpret_cast<uint8_t*>(file_position));
      file_position += sizeof(uint8_t);
    }
    // Optional scatter correction field: [FLTNBDATA (scatter)]
    if (m_scatCorrectionFlag)
    {
      if (!m_ignoreScatCorrectionFlag) event->SetScatterRate(0, *reinterpret_cast<FLTNBDATA*>(file_position)); 
      file_position += sizeof(FLTNBDATA);
    }
    // Optional random correction field: [FLTNBDATA (random)]
    if (m_randCorrectionFlag)
    {
      if (!m_ignoreRandCorrectionFlag) event->SetRandomRate(*reinterpret_cast<FLTNBDATA*>(file_position)); 
      file_position += sizeof(FLTNBDATA);
    }
    // Optional normalization correction field: [FLTNBDATA (norm)]
    if (m_normCorrectionFlag)
    {
      if (!m_ignoreNormCorrectionFlag) event->SetNormalizationFactor(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }
    // Optional TOF information field: [FLTNBDATA (TOF)]
    if (m_TOFInfoFlag)
    {
      if (!m_ignoreTOFFlag) event->SetTOFMeasurement(*reinterpret_cast<FLTNBDATA*>(file_position)); 
      file_position += sizeof(FLTNBDATA);
    }
    // Optional POI correction fields: [FLTNBDATA (POI1[3])] [FLTNBDATA (POI2[3])]
    for (int i=0; i<3; i++)
    {
      if (mp_POIDirectionFlag[i])
      {
        if (!m_ignorePOIFlag) event->SetPOI1(i,*reinterpret_cast<FLTNBDATA*>(file_position)); 
        file_position += sizeof(FLTNBDATA);
        if (!m_ignorePOIFlag) event->SetPOI2(i,*reinterpret_cast<FLTNBDATA*>(file_position)); 
        file_position += sizeof(FLTNBDATA);
      }
    }
    // Mandatory nb contributing LORs if m_maxNumberOfLinesPerEvent>1: [uint16_t nbLines]
    if (m_maxNumberOfLinesPerEvent>1)
    {
      event->SetNbLines(*reinterpret_cast<uint16_t*>(file_position));
      file_position += sizeof(uint16_t);
    }
    // Mandatory crystal IDs: [uint32_t (ID1)] [uint32_t (ID2)]
    for (int i=0 ; i<event->GetNbLines() ; i++)
    {
      event->SetID1(i, *reinterpret_cast<uint32_t*>(file_position)); 
      file_position += sizeof(uint32_t);
      event->SetID2(i, *reinterpret_cast<uint32_t*>(file_position)); 
      file_position += sizeof(uint32_t);
    }
  }

  // For MODE_HISTOGRAM PET data
  else if (m_dataMode == MODE_HISTOGRAM) 
  {
    // Cast the event pointer
    iEventHistoPET* event = (dynamic_cast<iEventHistoPET*>(m2p_BufferEvent[a_th]));
    // Mandatory time field: [uint32_t (time)]
    event->SetTimeInMs(*reinterpret_cast<uint32_t*>(file_position)); 
    file_position += sizeof(uint32_t);
    // Optional attenuation correction field: [FLTNBDATA (attnCorrectionFactor)]
    if (m_atnCorrectionFlag)
    {
      if (!m_ignoreAttnCorrectionFlag) event->SetAttenuationCorrectionFactor(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }

    // Optional random correction field: [FLTNBDATA (random)]
    if (m_randCorrectionFlag)
    {
      if (!m_ignoreRandCorrectionFlag) event->SetRandomRate(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }

    // Optional normalization correction field: [FLTNBDATA (norm)]
    if (m_normCorrectionFlag)
    {
      if (!m_ignoreNormCorrectionFlag) event->SetNormalizationFactor(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }
    // If we ignore the TOF information, then we have to sum up the bins' contributions
    if (m_ignoreTOFFlag)
    {
      // Compute the total of event values and scatter rates
      FLTNBDATA total_event_value = 0.;
      FLTNBDATA total_scatter_rate = 0.;
      for (int tb=0 ; tb<m_nbBinsTOF ; tb++)
      {
        // Mandatory bin value: [FLTNBDATA bin value]
        total_event_value += *reinterpret_cast<FLTNBDATA*>(file_position);
        file_position += sizeof(FLTNBDATA);
        // Optional scatter correction field: [FLTNBDATA (scatter)]
        if (m_scatCorrectionFlag)
        {
          total_scatter_rate += *reinterpret_cast<FLTNBDATA*>(file_position);
          file_position += sizeof(FLTNBDATA);
        }
      }
      // Affect the total to the event
      event->SetEventValue(0,total_event_value);
      if (!m_ignoreScatCorrectionFlag) event->SetScatterRate(0,total_scatter_rate);
    }
    // Otherwise we set all different TOF bin contributions
    else
    {
      for (int tb=0 ; tb<m_nbBinsTOF ; tb++)
      {
        // Mandatory bin value: [FLTNBDATA bin value]
        event->SetEventValue(tb, *reinterpret_cast<FLTNBDATA*>(file_position));
        file_position += sizeof(FLTNBDATA);
        // Optional scatter correction field: [FLTNBDATA (scatter)]
        if (m_scatCorrectionFlag)
        {
          if (!m_ignoreScatCorrectionFlag) event->SetScatterRate(tb, *reinterpret_cast<FLTNBDATA*>(file_position));
          file_position += sizeof(FLTNBDATA);
        }
      }
    }
    // Mandatory nb contributing LORs if m_maxNumberOfLinesPerEvent>1: [uint16_t nbLines]
    if (m_maxNumberOfLinesPerEvent>1)
    {
      event->SetNbLines(*reinterpret_cast<uint16_t*>(file_position));
      file_position += sizeof(uint16_t);
    }
    // Mandatory crystal IDs: [uint32_t (c1)] [uint32_t (c2)]
    for (int i=0 ; i<event->GetNbLines() ; i++)
    {
      event->SetID1(i, *reinterpret_cast<uint32_t*>(file_position));
      file_position += sizeof(uint32_t);
      event->SetID2(i, *reinterpret_cast<uint32_t*>(file_position));
      file_position += sizeof(uint32_t); 
    }
  }

  // For MODE_NORMALIZATION PET data
  else if (m_dataMode == MODE_NORMALIZATION)
  {
    // Cast the event pointer
    iEventNorm* event = (dynamic_cast<iEventNorm*>(m2p_BufferEvent[a_th]));
    // Optional attenuation correction field: [FLTNBDATA (norm)]
    if (m_atnCorrectionFlag)
    {
      if (!m_ignoreAttnCorrectionFlag) event->SetAttenuationCorrectionFactor(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }
    // Optional normalization correction field: [FLTNBDATA (norm)]
    if (m_normCorrectionFlag)
    {
      if (!m_ignoreNormCorrectionFlag) event->SetNormalizationFactor(*reinterpret_cast<FLTNBDATA*>(file_position));
      file_position += sizeof(FLTNBDATA);
    }
    // Mandatory nb contributing LORs if m_maxNumberOfLinesPerEvent>1: [uint16_t nbLines]
    if (m_maxNumberOfLinesPerEvent>1)
    {
      event->SetNbLines(*reinterpret_cast<uint16_t*>(file_position));
      file_position += sizeof(uint16_t);
    }
    // Mandatory crystal IDs: [uint32_t (c1)] [uint32_t (c2)]
    for (int i=0 ; i<event->GetNbLines() ; i++)
    {
      event->SetID1(i, *reinterpret_cast<uint32_t*>(file_position));
      file_position += sizeof(uint32_t);
      event->SetID2(i, *reinterpret_cast<uint32_t*>(file_position));
      file_position += sizeof(uint32_t); 
    }
  }

  // Return the updated event
  return m2p_BufferEvent[a_th];
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iDataFilePET::DescribeSpecific()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose==0) return;
  // Describe the datafile
  Cout("iDataFilePET::DescribeSpecific() -> Here is some specific content of the PET datafile" << endl);
  Cout("  --> Isotope: " << m_isotope << endl);
  if (m_TOFInfoFlag)
  {
    Cout("  --> TOF resolution: " << m_resolutionTOF << " ps" << endl);
    if (m_dataMode==MODE_LIST)
    {
      Cout("  --> TOF measurement range (for list-mode data): " << m_TOFMeasurementRange << " ps" << endl);
    }
    else if (m_dataMode==MODE_HISTOGRAM)
    {
      Cout("  --> Number of TOF bins: " << m_nbBinsTOF << endl);
      Cout("  --> TOF bin size: " << m_binSizeTOF << " ps" << endl);
    }
  }
  if (m_dataMode==MODE_LIST && m_eventKindFlag) Cout("  --> Event kind is present" << endl);
  if (m_scatCorrectionFlag) Cout("  --> Scatter correction is present" << endl);
  if (m_randCorrectionFlag) Cout("  --> Random correction is present" << endl);
  if (m_atnCorrectionFlag) Cout("  --> Attenuation correction is present" << endl);
  if (m_normCorrectionFlag) Cout("  --> Normalization correction is present" << endl);
  if (m_maxNumberOfLinesPerEvent>1) Cout("  --> Maximum number of lines per event: " << m_maxNumberOfLinesPerEvent << endl);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::CheckSpecificParameters()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Error if m_dataType != PET
  if (m_dataType != TYPE_PET)
  {
    Cerr("***** iDataFilePET::CheckSpecificParameters() -> Data type should be PET !'" << endl);
    return 1;
  }
  // Check if Maximum ring difference has been initialized, use default value (negative) otherwise.
  if (m_maxAxialDiffmm < 0.)
  {
    m_maxAxialDiffmm = -1.;
  }
  // Some checks for POI use
  if (m_POIInfoFlag)
  {
    // Not compatible with histogram data
    if (m_dataMode==MODE_HISTOGRAM)
    {
      Cerr("***** iDataFilePET::CheckSpecificParameters() -> POI correction flag is enabled while the data are histogrammed, no sense !" << endl);
      return 1;
    }
    // For each direction (radial, tangential and axial), look at the resolution. If not negative, the info should be the datafile
    for (int i=0; i<3; i++) if (mp_POIResolution[i]>=0.) mp_POIDirectionFlag[i] = true;
  }
  // Some checks for TOF use
  if (m_TOFInfoFlag)
  {
    // Check if the resolution was provided
    if (m_resolutionTOF<0.)
    {
      Cerr("***** iDataFilePET::CheckSpecificParameters() -> TOF information is used while there is no TOF resolution specified in the datafile header !" << endl);
      return 1;
    }
    // For histogram data
    if (m_dataMode==MODE_HISTOGRAM)
    {
      // Check if the number of bins has been provided
      if (m_nbBinsTOF<=1)
      {
        Cerr("***** iDataFilePET::CheckSpecificParameters() -> TOF information is used while there is only one TOF bin (specified in the datafile header) !" << endl);
        return 1;
      }
      // Check if the bin size has been provided
      if (m_binSizeTOF<=0.)
      {
        Cerr("***** iDataFilePET::CheckSpecificParameters() -> TOF information is used while there is no bin size specified in the datafile header !" << endl);
        return 1;
      }
    }
    // For list mode data
    else if (m_dataMode==MODE_LIST)
    {
      // set the number of TOF bins to 1
      m_nbBinsTOF = 1;
      return 0;
    }
  }
  else
  {
    // Set the number of TOF bins to 1
    m_nbBinsTOF = 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::CheckFileSizeConsistency()
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
    Cerr("--------------------------------------------------------------------------------------------------------------------------------------" << endl);
    Cerr("***** iDataFilePET::CheckFileSizeConsistency() -> DataFile size is not consistent with the information provided by the user/datafile !" << endl);
    Cerr("  --> Expected size : "<< m_nbEvents*m_sizeEvent << endl);
    Cerr("  --> Actual size : "<< sizeInBytes << endl << endl);
    Cerr("      ADDITIONAL INFORMATION ABOUT THE DATAFILE INITIALIZATION : " << endl);
    if (m_maxNumberOfLinesPerEvent > 1) Cerr("  --> Compression enabled. Each event should contain " << m_maxNumberOfLinesPerEvent <<
                                             "LORs, or an equivalent number of LOR + garbage LORs" << endl);
    else Cerr("  --> No compression in the data" << endl);
    if (m_eventKindFlag) Cerr("  --> Coincidence kind term is enabled" << endl);
    else Cerr("  --> No information about the kind of coincidences in the data" << endl);
    if (m_normCorrectionFlag) Cerr("  --> Normalization correction term is enabled" << endl);
    else Cerr("  --> No normalization correction term in the data" << endl);
    if (m_atnCorrectionFlag) Cerr("  --> Attenuation correction term is enabled" << endl);
    else Cerr("  --> No attenuation correction term in the data" << endl);
    if (m_scatCorrectionFlag) Cerr("  --> Scatter correction term is enabled" << endl);
    else Cerr("  --> No scatter correction term in the data" << endl);
    if (m_randCorrectionFlag) Cerr("  --> Random correction term is enabled" << endl);
    else Cerr("  --> No random correction term in the data" << endl);
    if (m_TOFInfoFlag) Cerr("  --> TOF information is enabled. Resolution is: " << m_resolutionTOF << " ps"<< endl);
    else Cerr("  --> No TOF information in the data" << endl);
    Cerr("  --> Calibration factor value is: " << m_calibrationFactor << endl);
    if (mp_POIResolution[0]<0.) Cerr("  --> No POI enabled on the radial axis" << endl);
    else Cerr("  --> POI resolution on the radial axis is: " << mp_POIResolution[0] << endl);
    if (mp_POIResolution[1]<0.) Cerr("  --> No POI enabled on the tangential axis" << endl);
    else Cerr("  --> POI resolution on the tangential axis is: " << mp_POIResolution[1] << endl);
    if (mp_POIResolution[2]<0.) Cerr("  --> No POI enabled on the axial axis" << endl);
    else Cerr("  --> POI resolution on the axial axis is: " << mp_POIResolution[2] << endl);
    Cerr("--------------------------------------------------------------------------------------------------------------------------------------" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile(vDataFile* ap_DataFile)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Dynamic cast the vDataFile to a iDataFilePET
  iDataFilePET* p_data_file = (dynamic_cast<iDataFilePET*>(ap_DataFile));
  // Check maximum number of lines per event
  if (m_maxNumberOfLinesPerEvent!=p_data_file->GetMaxNumberOfLinesPerEvent())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Maximum numbers of lines by event are inconsistent !" << endl);
    return 1;
  }
  // Check maximum ring difference
  if (m_maxAxialDiffmm!=p_data_file->GetMaxAxialDiffmm())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Maximum ring differences are inconsistent !" << endl);
    return 1;
  }
  // Check isotope
  if (m_isotope!=p_data_file->GetIsotope())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Isotopes are inconsistent !" << endl);
    return 1;
  }
  // Check event kind flag
  if (m_eventKindFlag!=p_data_file->GetEventKindFlag())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Event kind flags are inconsistent !" << endl);
    return 1;
  }
  // Check attenuation correction flag
  if (m_atnCorrectionFlag!=p_data_file->GetAtnCorrectionFlag())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Attenuation correction flags are inconsistent !" << endl);
    return 1;
  }
  // Check normalization correction flag
  if (m_normCorrectionFlag!=p_data_file->GetNormCorrectionFlag())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Normalization correction flags are inconsistent !" << endl);
    return 1;
  }
  // Check scatter correction flag
  if (m_scatCorrectionFlag!=p_data_file->GetScatCorrectionFlag())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Scatter correction flags are inconsistent !" << endl);
    return 1;
  }
  // Check random correction flag
  if (m_randCorrectionFlag!=p_data_file->GetRandCorrectionFlag())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Random correction flags are inconsistent !" << endl);
    return 1;
  }
  // Check TOF info flag
  if (m_TOFInfoFlag!=p_data_file->GetTOFInfoFlag())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> TOF information flags are inconsistent !" << endl);
    return 1;
  }
  // Check TOF resolution
  if (m_resolutionTOF!=p_data_file->GetTOFResolution())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> TOF resolutions are inconsistent !" << endl);
    return 1;
  }
  // Check number of TOF bins
  if (m_nbBinsTOF!=p_data_file->GetNbTOFBins())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Numbers of TOF bins are inconsistent !" << endl);
    return 1;
  }
  // Check TOF bin size
  if (m_binSizeTOF!=p_data_file->GetTOFBinSize())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> TOF bin sizes are inconsistent !" << endl);
    return 1;
  }
  // Check data mode
  if (m_dataMode!=p_data_file->GetDataMode())
  {
    Cerr("***** iDataFilePET::CheckSpecificConsistencyWithAnotherDataFile() -> Data modes are inconsistent (list-mode or histogram) !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::PROJ_InitFile()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("iDataFilePET::PROJ_InitFile() -> Initialize datafile for writing"<< endl); 
    
  m_nbEvents = 0; //Std initialization for projection
  m_nbBinsTOF = 1; // Initialization of this variable required for PROJ_Write functions
  
  m_isotope = "unknown";
  
  // Default time initialization if image dimensions object has not been initialized
  if(mp_ID->IsInitialized())
  {
    m_startTimeInSec = mp_ID->GetFrameTimeStartInSec(0, 0);
    
    for(uint16_t fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      m_durationInSec += mp_ID->GetFrameDurationInSec(0, fr);
  
    // Instanciate a fstream datafile for each thread
    m2p_dataFile = new fstream*[mp_ID->GetNbThreadsForProjection()];
  }
  else
  {
    m_startTimeInSec = 0.;
    m_durationInSec = 1.;
    m2p_dataFile = new fstream*[1];
  }



  // Set name of the projection data header
  sOutputManager* p_OutMgr;
  p_OutMgr = sOutputManager::GetInstance();
  string path_name = p_OutMgr->GetPathName();
  string img_name = p_OutMgr->GetBaseName();
  m_headerFileName = path_name.append(img_name).append("_CstrProj").append(".Cdh");

  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    m_dataFileName = m_headerFileName.substr(0, m_headerFileName.find_last_of(".")).append(".Cdf"); // Generate datafile name from header file
    
    // We'll write projeted data in several files (corresponding to the number of thread) to be concatenated at the end of the projection process.   
    stringstream ss;
    ss << th;
    string datafile_name = m_dataFileName;
    datafile_name.append("_").append(ss.str());
    
    m2p_dataFile[th] = new fstream( datafile_name.c_str(), ios::binary | ios::out | ios::trunc);
  }


  // Check there is no existing datafile with such name. Remove it otherwise
  ifstream fcheck(m_dataFileName.c_str());
  if(fcheck.good())
  {
    fcheck.close();
    #ifdef _WIN32
    string dos_instruction = "del " + m_dataFileName;
    system(dos_instruction.c_str());
    #else
    remove(m_dataFileName.c_str());
    #endif
  }
  
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Output datafile name :" << m_dataFileName << endl); 
    
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::WriteEvent(vEvent* ap_Event, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  if (m_dataMode == MODE_LIST) 
  {
    if (WriteListEvent((iEventListPET*)ap_Event, a_th))
    {
      Cerr("*****iDataFilePET::WriteEvent() -> Error while trying to write projection datafile (list-mode)" << endl;);
      return 1;
    }
  }

  if (m_dataMode == MODE_HISTOGRAM)
  {
    if (WriteHistoEvent((iEventHistoPET*)ap_Event, a_th))
    {
      Cerr("*****iDataFilePET::WriteEvent() -> Error while trying to write projection datafile (histogram)" << endl;);
      return 1;
    }
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::WriteHistoEvent(iEventHistoPET* ap_Event, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // Write sequentially each field of the event according to the type of the event.
  m2p_dataFile[a_th]->clear();
  
  uint32_t time = ap_Event->GetTimeInMs();
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&time), sizeof(uint32_t));
  
  if(m_atnCorrectionFlag)
  {
    FLTNB atn_corr_factor = ap_Event->GetAtnCorrFactor(); 
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&atn_corr_factor), sizeof(FLTNB));
  }
  
  if(m_randCorrectionFlag) 
  {
    FLTNB rdm_rate = ap_Event->GetEventRdmRate();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&rdm_rate), sizeof(FLTNB));
  }
  
  if(m_normCorrectionFlag)   
  {
    FLTNB norm_corr_factor = ap_Event->GetNormFactor();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&norm_corr_factor), sizeof(FLTNB));
  }
  
  for(int b=0 ; b<m_nbBinsTOF ; b++)
  {
    FLTNB event_value = ap_Event->GetEventValue(b);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&event_value), sizeof(FLTNB));
    if(m_scatCorrectionFlag) 
    {
      FLTNB scat_rate = ap_Event->GetEventScatRate(b);
      m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&scat_rate), sizeof(FLTNB));
    }
  }
  
  uint16_t nb_lines = ap_Event->GetNbLines();
  // Write the number of lines only if the m_maxNumberOfLinesPerEvent is above 1
  if (m_maxNumberOfLinesPerEvent>1) m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&nb_lines), sizeof(uint16_t));

  for(int i=0 ; i<nb_lines ; i++)
  {
    uint32_t id1 = ap_Event->GetID1(i);
    uint32_t id2 = ap_Event->GetID2(i);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&id1), sizeof(uint32_t));
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&id2), sizeof(uint32_t));
  }

  // 0-filling if needed (nb lines inferior to the max number for PET data with compression)
  if(nb_lines<m_maxNumberOfLinesPerEvent)
    for(int i=0 ; i<m_maxNumberOfLinesPerEvent-nb_lines ; i++)
    {
      uint32_t gbg = 0;
      m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&gbg), sizeof(uint32_t));
      m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&gbg), sizeof(uint32_t));
    }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::WriteListEvent(iEventListPET* ap_Event, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // Write sequentially each field of the event according to the type of the event.
  m2p_dataFile[a_th]->clear();
  
  uint32_t time = ap_Event->GetTimeInMs();
  m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&time), sizeof(uint32_t));
  
  
  if(m_atnCorrectionFlag)
  {
    FLTNB atn_corr_factor = ap_Event->GetAtnCorrFactor(); 
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&atn_corr_factor), sizeof(FLTNB));
  }
  
  if(m_eventKindFlag) 
  {
    uint8_t event_kind = ap_Event->GetKind();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&event_kind), sizeof(uint8_t));
  }
  
  if(m_scatCorrectionFlag)   
  {
    FLTNB scat_rate = ap_Event->GetEventScatRate(0);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&scat_rate), sizeof(FLTNB));
  }
  
  if(m_randCorrectionFlag) 
  {
    FLTNB rdm_rate = ap_Event->GetEventRdmRate();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&rdm_rate), sizeof(FLTNB));
  }
  
  if(m_normCorrectionFlag)   
  {
    FLTNB norm_corr_factor = ap_Event->GetNormFactor();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&norm_corr_factor), sizeof(FLTNB));
  }
    
  if(m_TOFInfoFlag)
  {
    FLTNBDATA TOF_measurement = ap_Event->GetTOFMeasurement();
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&TOF_measurement), sizeof(FLTNBDATA));
  }
  
  
  if(mp_POIResolution[0]>0)
  {
    FLTNB POI_1 = ap_Event->GetPOI1(0);
    FLTNB POI_2 = ap_Event->GetPOI2(0);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI_1), sizeof(FLTNB));
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI_2), sizeof(FLTNB));
  }
  if(mp_POIResolution[1]>0)
  {
    FLTNB POI_1 = ap_Event->GetPOI1(1);
    FLTNB POI_2 = ap_Event->GetPOI2(1);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI_1), sizeof(FLTNB));
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI_2), sizeof(FLTNB));
  }
  if(mp_POIResolution[2]>0)
  {
    FLTNB POI_1 = ap_Event->GetPOI1(2);
    FLTNB POI_2 = ap_Event->GetPOI2(2);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI_1), sizeof(FLTNB));
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&POI_2), sizeof(FLTNB));
  }
  
  uint16_t nb_lines = ap_Event->GetNbLines();
  // Write the number of lines only if the m_maxNumberOfLinesPerEvent is above 1
  if (m_maxNumberOfLinesPerEvent>1) m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&nb_lines), sizeof(uint16_t));

  for(int i=0 ; i<nb_lines ; i++)
  {
    uint32_t id1 = ap_Event->GetID1(i);
    uint32_t id2 = ap_Event->GetID2(i);
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&id1), sizeof(uint32_t));
    m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&id2), sizeof(uint32_t));
  }
  
  // 0-filling if needed (nb lines inferior to the max number for PET data with compression)
  if(nb_lines<m_maxNumberOfLinesPerEvent)
    for(int i=0 ; i<m_maxNumberOfLinesPerEvent-nb_lines ; i++)
    {
      uint32_t gbg = 0;
      m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&gbg), sizeof(uint32_t));
      m2p_dataFile[a_th]->write(reinterpret_cast<char*>(&gbg), sizeof(uint32_t));
    }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::WriteHeader()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_NORMAL) Cout("iDataFilePET::WriteHeader() -> Write header file '" << m_headerFileName << "'" << endl);
  // Open file
  fstream headerFile;
  headerFile.open(m_headerFileName.c_str(), ios::out);
  
  if (!headerFile.is_open())
  {
    Cerr("***** iDataFilePET::WriteHeader() -> Failed to open output header file '" << m_headerFileName << "' !" << endl);
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
  // PET data type
  headerFile << "Data type: PET" << endl;
  // Acquisition start time in seconds
  headerFile << "Start time (s): " << m_startTimeInSec << endl; 
  // Acquisition duration in seconds
  headerFile << "Duration (s): " << m_durationInSec << endl; 
  // Scanner name
  headerFile << "Scanner name: " << sScannerManager::GetInstance()->GetScannerName() << endl; 
  // Maximum axial difference
  if(m_maxAxialDiffmm > 0.)
    headerFile << "Maximum axial difference mm: " << m_maxAxialDiffmm << endl; 
  // Maximum number of lines per event (write it only if higher than 1
  if (m_maxNumberOfLinesPerEvent >1) headerFile << "Maximum number of lines per event: " << m_maxNumberOfLinesPerEvent << endl; 

  // Calibration factor
  headerFile << "Calibration factor: " << m_calibrationFactor << endl;
  // Isotope
  headerFile << "Isotope: " << m_isotope << endl; 
  // TOF info
  if (m_TOFInfoFlag)
  {
    // Flag
    headerFile << "TOF information flag: 1" << endl;
    // TOF resolution
    headerFile << "TOF resolution (ps): " << m_resolutionTOF << endl; 
    // List-mode TOF measurement range
    if (m_dataMode == MODE_LIST) headerFile << "List TOF measurement range (ps): " << m_TOFMeasurementRange << endl;
    // Histogram
    else if (m_dataMode == MODE_HISTOGRAM) 
    {
      // Number of TOF bins
      headerFile << "Histo TOF number of bins: " << m_nbBinsTOF << endl;
      // TOF bin size
      headerFile << "Histo TOF bin size (ps): " << m_binSizeTOF << endl; 
    }
  }
  // Event kind
  if (m_eventKindFlag)
  {
    // Error if histogram data are used
    if (m_dataMode==MODE_HISTOGRAM)
    {
      Cerr("***** iDataFilePET::WriteHeader -> Request for writing event type in histogram mode (this field is specific to list-mode) !" << endl);
      return 1;
    }
    else headerFile << "Coincidence kind flag: " << m_eventKindFlag << endl;
  }
  // Correction flags
  if (m_atnCorrectionFlag)
    headerFile << "Attenuation correction flag: " << m_atnCorrectionFlag << endl;
  if (m_normCorrectionFlag)
    headerFile << "Normalization correction flag: " << m_normCorrectionFlag << endl; 
  if (m_scatCorrectionFlag)
    headerFile << "Scatter correction flag: " << m_scatCorrectionFlag << endl; 
  if (m_randCorrectionFlag)
    headerFile << "Random correction flag: " << m_randCorrectionFlag << endl; 
  // Close file
  headerFile.close();
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iDataFilePET::PROJ_GetScannerSpecificParameters()
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("PROJ_GetScannerSpecificParameters() ..."<< endl);
  // Get pointers to SPECT specific parameters in scanner
  if (sScannerManager::GetInstance()->PROJ_GetPETSpecificParameters(&m_maxAxialDiffmm) )
  {
    Cerr("***** iDataFilePET::PROJ_GetScannerSpecificParameters() -> An error occurred while trying to get PET geometric parameters from the scanner object !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
