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
  \ingroup  projector
  \brief    Implementation of class oProjectorManager
*/

#include "oProjectorManager.hh"
#include "sOutputManager.hh"
#include "sAddonManager.hh"
#include "iDataFilePET.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oProjectorManager::oProjectorManager()
{
  // Scanner and image dimensions
  mp_Scanner = NULL;
  mp_ImageDimensionsAndQuantification = NULL;
  // Data file (used to get some info about TOF bins and POI)
  mp_DataFile = NULL;
  // TOF and POI options
  m_applyTOF = false;
  m_applyPOI = false;
  m_nbTOFBins = -1;
  // Computation strategy for projection lines
  m_computationStrategy = -1;
  // Forward and backward options for projectors
  m_optionsForward = "";
  m_optionsBackward = "";
  // Common options
  m_optionsCommon = "";
  // Forward and backward projectors
  m_forwardProjectorName = "";
  m_backwardProjectorName = "";
  mp_SystemMatrixForward = NULL;
  mp_SystemMatrixBackward = NULL;
  mp_ProjectorForward = NULL;
  mp_ProjectorBackward = NULL;
  m_useSystemMatrixForward = false;
  m_useSystemMatrixBackward = false;
  m_useProjectorForward = false;
  m_useProjectorBackward = false;
  m_useMatchedProjectors = false;
  // Forward and backward projection lines (as many as threads)
  m2p_ProjectionLines = NULL;
  // Verbosity
  m_verbose = 0;
  // Not checked yet
  m_checked = false;
  // Not initialized yet
  m_initialized = false;
  m_applyMask = false;
  mp_mask = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oProjectorManager::~oProjectorManager() 
{
  // Go through the destructor only if the object was initialized
  if (m_initialized)
  {
    // Delete projection lines
    if (m2p_ProjectionLines)
    {
      for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++)
        if (m2p_ProjectionLines[th]) delete m2p_ProjectionLines[th];
      delete[] m2p_ProjectionLines;
    }
    // Delete projectors
    if (m_useProjectorForward) delete mp_ProjectorForward;
    if (m_useSystemMatrixForward) delete mp_SystemMatrixForward;
    if (m_applyMask) delete[] mp_mask;
    if (!m_useMatchedProjectors)
    {
      if (m_useProjectorBackward) delete mp_ProjectorBackward;
      if (m_useSystemMatrixBackward) delete mp_SystemMatrixBackward;
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

bool oProjectorManager::IsForwardOperatorCompatibleWithSPECTAttenuationCorrection()
{
  // Case for a vProjector
  if (m_useProjectorForward) return mp_ProjectorForward->GetCompatibilityWithSPECTAttenuationCorrection();
  // Case for a system matrix
  else return mp_SystemMatrixForward->GetCompatibilityWithSPECTAttenuationCorrection();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

bool oProjectorManager::IsBackwardOperatorCompatibleWithSPECTAttenuationCorrection()
{
  // Case for a vProjector
  if (m_useProjectorBackward) return mp_ProjectorBackward->GetCompatibilityWithSPECTAttenuationCorrection();
  // Case for a system matrix
  else return mp_SystemMatrixBackward->GetCompatibilityWithSPECTAttenuationCorrection();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oProjectorManager::CheckParameters()
{
  // Check scanner
  if (mp_Scanner==NULL)
  {
    Cerr("***** oProjectorManager::CheckParameters() -> No scanner provided !" << endl);
    return 1;
  }
  // Check image dimensions
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** oProjectorManager::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  // Check data file
  if (mp_DataFile==NULL)
  {
    Cerr("***** oProjectorManager::CheckParameters() -> No data file provided !" << endl);
    return 1;
  }
  // Check computation strategy
  if ( m_computationStrategy != IMAGE_COMPUTATION_STRATEGY &&
       m_computationStrategy != FIXED_LIST_COMPUTATION_STRATEGY &&
       m_computationStrategy != ADAPTATIVE_LIST_COMPUTATION_STRATEGY )
  {
    Cerr("***** oProjectorManager::CheckParameters() -> Unknown computation strategy provided !" << endl);
    return 1;
  }
  // Check forward projector options
  if (m_optionsForward=="")
  {
    Cerr("***** oProjectorManager::CheckParameters() -> No forward projector options provided !" << endl);
    return 1;
  }
  // Check backward projector options
  if (m_optionsBackward=="")
  {
    Cerr("***** oProjectorManager::CheckParameters() -> No backward projector options provided !" << endl);
    return 1;
  }
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** oProjectorManager::CheckParameters() -> Wrong verbosity level provided !" << endl);
    return 1;
  }
  // Normal end
  m_checked = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oProjectorManager::CheckSPECTAttenuationCompatibility(const string& a_pathToAttenuationImage)
{
  // In SPECT with attenuation correction, there are some requirements with the projection method
  if (a_pathToAttenuationImage!="" && mp_DataFile->GetDataType()==TYPE_SPECT)
  {
    // Check that the projection line strategy is not IMAGE_COMPUTATION_STRATEGY (this is not compatible)
    // Note that we cannot do this check in the oProjectionLine directly, because we cannot now in advance
    // that the projection method including attenuation will be used...
    if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
    {
      Cerr("***** oProjectorManager::CheckSPECTAttenuationCompatibility() -> The image-computation strategy of the oProjectionLine is not compatible with SPECT attenuation correction !");
      return 1;
    }
    // Check that the forward projector is compatible with SPECT with attenuation correction
    if (!IsForwardOperatorCompatibleWithSPECTAttenuationCorrection())
    {
      Cerr("***** oProjectorManager::CheckSPECTAttenuationCompatibility() -> The forward projector is not compatible with SPECT attenuation correction !" << endl);
      return 1;
    }
    // Check that the backward projector is compatible with SPECT with attenuation correction
    if (!IsBackwardOperatorCompatibleWithSPECTAttenuationCorrection())
    {
      Cerr("***** oProjectorManager::CheckSPECTAttenuationCompatibility() -> The backward projector is not compatible with SPECT attenuation correction !" << endl);
      return 1;
    }
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oProjectorManager::Initialize()
{
  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** oProjectorManager::Initialize() -> Must call CheckParameters() before Initialize() !" << endl);
    return 1;
  }

  // Verbose
  if (m_verbose>=1) Cout("oProjectorManager::Initialize() -> Initialize projectors and projection lines" << endl);

  // -------------------------------------------------------------------
  // Manage TOF
  // -------------------------------------------------------------------

  FLTNB tof_resolution = -1.;
  FLTNB tof_bin_size = -1.;
  // TOF is only for PET data
  if (mp_DataFile->GetDataType()==TYPE_PET)
  {
    // Cast the datafile pointer
    iDataFilePET* p_pet_datafile = (dynamic_cast<iDataFilePET*>(mp_DataFile));
    // Case of TOF information in the datafile
    if (p_pet_datafile->GetTOFInfoFlag())
    {
      // Case where it is asked to ignore TOF info
      if (p_pet_datafile->GetIgnoreTOFFlag())
      {
        m_applyTOF = USE_NOTOF;
        m_nbTOFBins = 1;
      }
      // Otherwise, we use TOF
      else
      {
        // Get the TOF resolution from the data file, if we apply TOF
        tof_resolution = p_pet_datafile->GetTOFResolution();
        // Get the TOF bin size in ps, if we apply TOF
        tof_bin_size = p_pet_datafile->GetTOFBinSize();
        // For list-mode data
        if (mp_DataFile->GetDataMode()==MODE_LIST)
        {
          m_nbTOFBins = 1;
          m_applyTOF = USE_TOFPOS;
        }
        // For histogram data
        else
        {
          m_applyTOF = USE_TOFBIN;
          // Get the number of TOF bins from the data file
          m_nbTOFBins = p_pet_datafile->GetNbTOFBins();
        }
      }
      // Verbose
      if (m_verbose>=VERBOSE_NORMAL)
      {
        if (m_applyTOF==USE_NOTOF) Cout("  --> Do not use the TOF projections" << endl);
        else if (m_applyTOF==USE_TOFBIN) Cout("  --> Use TOF projector for histogram data" << endl);
        else if (m_applyTOF==USE_TOFPOS) Cout("  --> Use TOF projector for listmode data" << endl);
      }
    }
    // Case no TOF
    else
    {
      m_applyTOF = USE_NOTOF;
      m_nbTOFBins = 1;
    }
  }
  // Not PET data
  else
  {
    m_nbTOFBins = 1;
    m_applyTOF = USE_NOTOF;
  }

  // -------------------------------------------------------------------
  // Manage POI
  // -------------------------------------------------------------------

  // Case we have POI and it is not asked to ignore the information
  if (mp_DataFile->GetPOIInfoFlag() && !mp_DataFile->GetIgnorePOIFlag())
  {
    // Cannot use POI with histogram mode
    if (mp_DataFile->GetDataMode()==MODE_HISTOGRAM)
    {
      Cerr("***** oProjectorManager::Initialize() -> POI information has no sense with histogram data !" << endl);
      return 1;
    }
    // Apply POI
    m_applyPOI = true;
  }
  // Case we do not use POI
  else m_applyPOI = false;
  // Get POI resolution from the data file
  FLTNB* poi_resolution = mp_DataFile->GetPOIResolution();

  // -------------------------------------------------------------------
  // Initialize projectors and or system matrix
  // -------------------------------------------------------------------

  // Compare projector options to know if we use matched ones for forward and backward operations
  if (m_optionsForward==m_optionsBackward) m_useMatchedProjectors = true;
  else m_useMatchedProjectors = false;

  // Parse projector options and initialize them
  if (ParseOptionsAndInitializeProjectors())
  {
    Cerr("***** oProjectorManager::Initialize() -> A problem occured while parsing projector options and initializing it !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------
  // If compression, loaded SM and some projectors are not compatible
  // -------------------------------------------------------------------

  // Compression can currently occur only in PET
  if (mp_DataFile->GetDataType()==TYPE_PET)
  {
    // Check if there is some compression
    if ( (dynamic_cast<iDataFilePET*>(mp_DataFile))->GetMaxNumberOfLinesPerEvent() > 1 )
    {
      // It is not compatible with loaded forward system matrices
      if (m_useSystemMatrixForward)
      {
        Cerr("***** oProjectorManager::Initialize() -> Cannot use a loaded forward system matrix with compression in the datafile !" << endl);
        return 1;
      }
      // It is not compatible with forward projectors declared as not compatible !
      else if (!mp_ProjectorForward->GetCompatibilityWithCompression())
      {
        Cerr("***** oProjectorManager::Initialize() -> Selected forward projector '" << m_forwardProjectorName << "' is not compatible with compression in the datafile !" << endl);
        return 1;
      }
      // It is not compatible with loaded backward system matrices
      if (m_useSystemMatrixBackward)
      {
        Cerr("***** oProjectorManager::Initialize() -> Cannot use a loaded backward system matrix with compression in the datafile !" << endl);
        return 1;
      }
      // It is not compatible with backward projectors declared as not compatible !
      else if (!mp_ProjectorBackward->GetCompatibilityWithCompression())
      {
        Cerr("***** oProjectorManager::Initialize() -> Selected backward projector '" << m_backwardProjectorName << "' is not compatible with compression in the datafile !" << endl);
        return 1;
      }
    }
  }

  // -------------------------------------------------------------------
  // Transfer the mask to projectors
  // -------------------------------------------------------------------

  if (m_applyMask)
  {
    mp_ProjectorForward->SetMask(mp_mask);
    if (!m_useMatchedProjectors) mp_ProjectorBackward->SetMask(mp_mask);
  }


  // -------------------------------------------------------------------
  // Initialize projection lines
  // -------------------------------------------------------------------

  // Initialize as many projection lines as threads
  m2p_ProjectionLines = new oProjectionLine*[mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()];
  for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++)
  {
    m2p_ProjectionLines[th] = new oProjectionLine();
    m2p_ProjectionLines[th]->SetThreadNumber(th);
    m2p_ProjectionLines[th]->SetMatchedProjectors(m_useMatchedProjectors);
    m2p_ProjectionLines[th]->SetNbTOFBins(m_nbTOFBins);
    m2p_ProjectionLines[th]->SetTOFBinSize(tof_bin_size);
    m2p_ProjectionLines[th]->SetTOFResolution(tof_resolution);
    m2p_ProjectionLines[th]->SetPOIResolution(poi_resolution);
    m2p_ProjectionLines[th]->SetImageDimensionsAndQuantification(mp_ImageDimensionsAndQuantification);
    m2p_ProjectionLines[th]->SetComputationStrategy(m_computationStrategy);
    m2p_ProjectionLines[th]->SetForwardProjector(mp_ProjectorForward);
    m2p_ProjectionLines[th]->SetBackwardProjector(mp_ProjectorBackward);
    m2p_ProjectionLines[th]->SetVerbose(m_verbose);
    if (m2p_ProjectionLines[th]->CheckParameters())
    {
      Cerr("***** oProjectorManager::Initialize() -> An error occured while checking parameters of an oProjectionLine !" << endl);
      return 1;
    }
    if (m2p_ProjectionLines[th]->Initialize())
    {
      Cerr("***** oProjectorManager::Initialize() -> An error occured while initializing an oProjectionLine !" << endl);
      return 1;
    }
  }

  // Normal end
  m_initialized = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oProjectorManager::ParseOptionsAndInitializeProjectors()
{
  string list_options = "";
  string file_options = "";

  // This is for the automatic initialization of the projectors
  typedef vProjector *(*maker_projector) ();

  // ---------------------------------------------------------------------------------------------------
  // Manage forward projector
  // ---------------------------------------------------------------------------------------------------

  // ______________________________________________________________________________
  // Get the projector name in the options and isolate the real projector's options

  // Search for a colon ":", this indicates that a configuration file is provided after the projector name
  size_t colon = m_optionsForward.find_first_of(":");
  size_t comma = m_optionsForward.find_first_of(",");

  // Case 1: we have a colon
  if (colon!=string::npos)
  {
    // Get the projector name before the colon
    m_forwardProjectorName = m_optionsForward.substr(0,colon);
    // Get the configuration file after the colon
    file_options = m_optionsForward.substr(colon+1);
    // List of options is empty
    list_options = "";
  }
  // Case 2: we have a comma
  else if (comma!=string::npos)
  {
    // Get the projector name before the first comma
    m_forwardProjectorName = m_optionsForward.substr(0,comma);
    // Get the list of options after the first comma
    list_options = m_optionsForward.substr(comma+1);
    // Configuration file is empty
    file_options = "";
  }
  // Case 3: no colon and no comma (a single projector name)
  else
  {
    // Get the projector name
    m_forwardProjectorName = m_optionsForward;
    // List of options is empty
    list_options = "";
    // Build the default configuration file
    file_options = sOutputManager::GetInstance()->GetPathToConfigDir() + "/projector/" + m_forwardProjectorName + ".conf";
  }

  // ______________________________________________________________________________
  // Case 1: projector is equal to keyword 'matrix', then use a system matrix
  if (m_forwardProjectorName==SYSTEM_MATRIX_KEYWORD)
  {
    mp_ProjectorForward = NULL;
    m_useProjectorForward = false;
    m_useSystemMatrixForward = true;
    // TODO: put all these limitations into a dedicated function from oSystemMatrix
    // TODO: forbid TOF in PET with system matrix
    // TODO: forbid simultaneous bed reconstruction with system matrix
    // TODO: forbid image offset with system matrix
    Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> Loading of custom system matrices is not yet implemented !" << endl);
    return 1;
  }
  // ______________________________________________________________________________
  // Case 2: on-the-fly projector
  else
  {
    // Unset system matrix
    mp_SystemMatrixForward = NULL;
    m_useSystemMatrixForward = false;
    // Set projector on
    m_useProjectorForward = true;
    // Get projector's listfrom addon manager
    std::map <string,maker_projector> list = sAddonManager::GetInstance()->mp_listOfProjectors;
    // Create the projector
    if (list[m_forwardProjectorName]) mp_ProjectorForward = list[m_forwardProjectorName]();
    else
    {
      Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> Projector '" << m_forwardProjectorName << "' does not exist !" << endl);
      return 1;
    }
    // Set parameters
    mp_ProjectorForward->SetScanner(mp_Scanner);
    if (mp_ProjectorForward->SetImageDimensionsAndQuantification(mp_ImageDimensionsAndQuantification))
    {
      Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while setting the image dimensions of the forward projector !" << endl);
      return 1;
    }
    mp_ProjectorForward->SetApplyTOF(m_applyTOF);
    mp_ProjectorForward->SetApplyPOI(m_applyPOI);
    mp_ProjectorForward->SetVerbose(m_verbose);
    // Provide common options list
    if (mp_ProjectorForward->ReadCommonOptionsList(m_optionsCommon))
    {
      Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while parsing and reading forward projector's common options !" << endl);
      return 1;
    }
    // Provide configuration file if any
    if (file_options!="" && mp_ProjectorForward->ReadConfigurationFile(file_options))
    {
      Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while parsing and reading forward projector's configuration file !" << endl);
      return 1;
    }
    // Provide options if any
    if (list_options!="" && mp_ProjectorForward->ReadOptionsList(list_options))
    {
      Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while parsing and reading forward projector's options !" << endl);
      return 1;
    }
    // Check parameters
    if (mp_ProjectorForward->CheckParameters())
    {
      Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while checking forward projector parameters !" << endl);
      return 1;
    }
    // Initialize the projector
    if (mp_ProjectorForward->Initialize())
    {
      Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while initializing the forward projector !" << endl);
      return 1;
    }
  }

  // ---------------------------------------------------------------------------------------------------
  // Manage backward projector
  // ---------------------------------------------------------------------------------------------------

  // If options are the same, then forward and backward are the same
  if (m_useMatchedProjectors)
  {
    // In this case, matched projectors
    m_useSystemMatrixBackward = m_useSystemMatrixForward;
    m_useProjectorBackward = m_useProjectorForward;
    mp_SystemMatrixBackward = mp_SystemMatrixForward;
    mp_ProjectorBackward = mp_ProjectorForward;
  }
  // Else, unmatched projectors
  else
  {
    // ______________________________________________________________________________
    // Get the projector name in the options and isolate the real projector's options

    // Search for a colon ":", this indicates that a configuration file is provided after the projector name
    colon = m_optionsBackward.find_first_of(":");
    comma = m_optionsBackward.find_first_of(",");

    // Case 1: we have a colon
    if (colon!=string::npos)
    {
      // Get the projector name before the colon
      m_backwardProjectorName = m_optionsBackward.substr(0,colon);
      // Get the configuration file after the colon
      file_options = m_optionsBackward.substr(colon+1);
      // List of options is empty
      list_options = "";
    }
    // Case 2: we have a comma
    else if (comma!=string::npos)
    {
      // Get the projector name before the first comma
      m_backwardProjectorName = m_optionsBackward.substr(0,comma);
      // Get the list of options after the first comma
      list_options = m_optionsBackward.substr(comma+1);
      // Configuration file is empty
      file_options = "";
    }
    // Case 3: no colon and no comma (a single projector name)
    else
    {
      // Get the projector name
      m_backwardProjectorName = m_optionsBackward;
      // List of options is empty
      list_options = "";
      // Build the default configuration file
      file_options = sOutputManager::GetInstance()->GetPathToConfigDir() + "/projector/" + m_backwardProjectorName + ".conf";
    }

    // ______________________________________________________________________________
    // Case 1: projector is equal to keyword 'matrix', then use a system matrix
    if (m_backwardProjectorName==SYSTEM_MATRIX_KEYWORD)
    {
      mp_ProjectorBackward = NULL;
      m_useProjectorBackward = false;
      m_useSystemMatrixBackward = true;
      // TODO: put all these limitations into a dedicated function from oSystemMatrix
      // TODO: forbid TOF in PET with system matrix
      // TODO: forbid simultaneous bed reconstruction with system matrix
      // TODO: forbid image offset with system matrix
      Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> Loading of custom system matrices is not yet implemented !" << endl);
      return 1;
    }
    // ______________________________________________________________________________
    // Case 2: on-the-fly projector
    else
    {
      // Unset system matrix
      mp_SystemMatrixBackward = NULL;
      m_useSystemMatrixBackward = false;
      // Set projector on
      m_useProjectorBackward = true;
      // Get projector's listfrom addon manager
      std::map <string,maker_projector> list = sAddonManager::GetInstance()->mp_listOfProjectors;
      // Create the projector
      if (list[m_backwardProjectorName]) mp_ProjectorBackward = list[m_backwardProjectorName]();
      else
      {
        Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> Projector '" << m_backwardProjectorName << "' does not exist !" << endl);
        return 1;
      }
      // Set parameters
      mp_ProjectorBackward->SetScanner(mp_Scanner);
      if (mp_ProjectorBackward->SetImageDimensionsAndQuantification(mp_ImageDimensionsAndQuantification))
      {
        Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while setting the image dimensions of the backward projector !" << endl);
        return 1;
      }
      mp_ProjectorBackward->SetApplyTOF(m_applyTOF);
      mp_ProjectorBackward->SetApplyPOI(m_applyPOI);
      mp_ProjectorBackward->SetVerbose(m_verbose);
      // Provide common options list
      if (mp_ProjectorBackward->ReadCommonOptionsList(m_optionsCommon))
      {
        Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while parsing and reading backward projector's common options !" << endl);
        return 1;
      }
      // Provide configuration file if any
      if (file_options!="" && mp_ProjectorBackward->ReadConfigurationFile(file_options))
      {
        Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while parsing and reading backward projector's configuration file !" << endl);
        return 1;
      }
      // Provide options if any
      if (list_options!="" && mp_ProjectorBackward->ReadOptionsList(list_options))
      {
        Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while parsing and reading backward projector's options !" << endl);
        return 1;
      }
      // Check parameters
      if (mp_ProjectorBackward->CheckParameters())
      {
        Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while checking backward projector parameters !" << endl);
        return 1;
      }
      // Initialize the projector
      if (mp_ProjectorBackward->Initialize())
      {
        Cerr("***** oProjectorManager::ParseOptionsAndInitializeProjectors() -> A problem occured while initializing backward projector !" << endl);
        return 1;
      }
    }
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectorManager::SetSensitivityModeOn()
{
  // Set the sensitivity ON in the projectors, if used
  if (m_useProjectorForward) mp_ProjectorForward->SetSensitivityMode(true);
  if (m_useProjectorBackward) mp_ProjectorBackward->SetSensitivityMode(true);
  // If the data file is a list-mode, then deactivate the TOF usage in the projector
  if (mp_DataFile->GetDataMode()==MODE_LIST)
  {
    if (m_useProjectorForward) mp_ProjectorForward->SetApplyTOF(USE_NOTOF);
    if (m_useProjectorBackward) mp_ProjectorBackward->SetApplyTOF(USE_NOTOF);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectorManager::SetSensitivityModeOff()
{
  // Set the sensitivity OFF in the projectors, if used
  if (m_useProjectorForward) mp_ProjectorForward->SetSensitivityMode(false);
  if (m_useProjectorBackward) mp_ProjectorBackward->SetSensitivityMode(false);
  // Reset TOF and POI parameters to original values
  if (m_useProjectorForward)
  {
    mp_ProjectorForward->SetApplyTOF(m_applyTOF);
    mp_ProjectorForward->SetApplyPOI(m_applyPOI);
  }
  if (m_useProjectorBackward)
  {
    mp_ProjectorBackward->SetApplyTOF(m_applyTOF);
    mp_ProjectorBackward->SetApplyPOI(m_applyPOI);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectorManager::ApplyBedOffset(int a_bed)
{
  // Apply it to all projection lines (only if more than one bed position)
  if (mp_ImageDimensionsAndQuantification->GetNbBeds()>1)
    for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++)
      m2p_ProjectionLines[th]->SetBedOffset(mp_ImageDimensionsAndQuantification->GetBedPosition(a_bed));
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oProjectionLine* oProjectorManager::ComputeProjectionLine(vEvent* ap_Event, int a_th) 
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectorManager::ComputeProjectionLine() -> Called while not initialized !" << endl);
    return NULL;
  }
  #endif

  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // Get the list of indices
  uint32_t *index1 = ap_Event->GetEventID1();
  uint32_t *index2 = ap_Event->GetEventID2();
  int nb_indices = ap_Event->GetNbLines();

  // Clean the projection line
  m2p_ProjectionLines[a_th]->Reset();

  // With list-mode data, we may need POI and/or TOF measurements
  if (ap_Event->GetDataMode()==MODE_LIST)
  {
    // Set POI measurement
    if (m_applyPOI)
    {
      // For PET
      if (ap_Event->GetDataType()==TYPE_PET)
      {
        // Have to dynamic_cast the event into a iEventPETList to access the GetPOI functions
        m2p_ProjectionLines[a_th]->SetPOI1((dynamic_cast<iEventListPET*>(ap_Event))->GetPOI1());
        m2p_ProjectionLines[a_th]->SetPOI2((dynamic_cast<iEventListPET*>(ap_Event))->GetPOI2());
      }
      // For SPECT
      else if (ap_Event->GetDataType()==TYPE_SPECT)
      {
        // By convention in SPECT, the second end point is in the camera (the first one being outside)
        //m2p_ProjectionLines[a_th]->SetPOI2((dynamic_cast<iEventListModeSPECT*>(ap_Event))->GetPOI());
      }
      // For CT
      else if (ap_Event->GetDataType()==TYPE_CT)
      {
        ;
      }
    }
    // Set TOF measurement (only for PET obviously)
    if (m_applyTOF && ap_Event->GetDataType()==TYPE_PET)
    {
      // Have to dynamic_cast the event into a iEventPETList to access the GetTOF function
      m2p_ProjectionLines[a_th]->SetTOFMeasurement((dynamic_cast<iEventListPET*>(ap_Event))->GetTOFMeasurement());
    }
  }

  // Project forward (and also compute line length)
  int return_value = 0;
  if (m_useProjectorForward)
    return_value = mp_ProjectorForward->Project( FORWARD, m2p_ProjectionLines[a_th], index1, index2, nb_indices );
  else if (m_useSystemMatrixForward)
    return_value = mp_SystemMatrixForward->Project( FORWARD, m2p_ProjectionLines[a_th], index1, index2, nb_indices );
  if (return_value)
  {
    Cerr("***** oProjectorManager::ComputeProjectionLine() -> A problem occured while forward projecting an event !" << endl);
    return NULL;
  }

  // Project backward
  if (!m_useMatchedProjectors)
  {
    // Then project
    if (m_useProjectorBackward)
      return_value = mp_ProjectorBackward->Project( BACKWARD, m2p_ProjectionLines[a_th], index1, index2, nb_indices );
    else if (m_useSystemMatrixBackward)
      return_value = mp_SystemMatrixBackward->Project( BACKWARD, m2p_ProjectionLines[a_th], index1, index2, nb_indices);
    if (return_value)
    {
      Cerr("***** oProjectorManager::ComputeProjectionLine() -> A problem occured while backward projecting an event !" << endl);
      return NULL;
    }
  }

  // Return the line
  return m2p_ProjectionLines[a_th];
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
int oProjectorManager::ProcessAndSetMask(FLTNB* ap_maskImage)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  if (mp_mask!=NULL)
  {
    Cerr("***** oProjectorManager::ProcessAndSetMask() -> Mask already initialized !" << endl);
    return 1;
  }

  mp_mask = new bool[mp_ImageDimensionsAndQuantification->GetNbVoxXYZ()];
  for (INTNB v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ();v++)
  {
    mp_mask[v] = (ap_maskImage[v]>0.);
  }
  m_applyMask = true;
  
  return 0;
}
