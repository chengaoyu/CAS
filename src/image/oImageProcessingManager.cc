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
  \brief    Implementation of class oImageProcessingManager
*/

#include "oImageProcessingManager.hh"
#include "sAddonManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oImageProcessingManager::oImageProcessingManager()
{
  // Image dimensions
  mp_ImageDimensionsAndQuantification = NULL;
  // Options
  m_options = {};
  // Image processing objects and associated bool
  m_nbImageProcessingModules = 0;
  m2p_ImageProcessingModules = NULL;
  mp_applyForward = NULL;
  mp_applyIntra = NULL;
  mp_applyPost = NULL;
  // Booleans
  m_checked = false;
  m_initialized = false;
  // Verbosity
  m_verbose = -1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oImageProcessingManager::~oImageProcessingManager() 
{
  // Delete object
  if (m2p_ImageProcessingModules)
  {
    for (int c=0; c<m_nbImageProcessingModules; c++) if (m2p_ImageProcessingModules[c]) delete m2p_ImageProcessingModules[c];
    free(m2p_ImageProcessingModules);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageProcessingManager::CheckParameters()
{
  // Check image dimensions
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** oImageProcessingManager::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** oImageProcessingManager::CheckParameters() -> Wrong verbosity level provided !" << endl);
    return 1;
  }
  // All set
  m_checked = true;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageProcessingManager::ShowCommonHelp()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "------------------------------------------------------------------" << endl;
  cout << "-----  How to use an image processing module" << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout << endl;
  cout << "An image processing module is called through the -proc option. The provided argument describes the processing module to be used," << endl;
  cout << "its options, and when to include it within the algorithm. The syntax of the argument must obey one of the three following options:" << endl;
  cout << "  proc::when (in this case, the default configuration file of the processing module is used to set the options values)" << endl;
  cout << "  proc:file.conf::when (in this case, the provided configuration is used)" << endl;
  cout << "  proc,param1,param2,...::when (in this case, the options values are directly provided in the argument)" << endl;
  cout << "In any case, the description of the options specific to each processing module, their order in the list and their configuration" << endl;
  cout << "files syntax are provided in the specific help of each module." << endl;
  cout << "The 'when' parameter is an argument describing when to include the processing module within the algorithm. It is a list of keywords" << endl;
  cout << "separating by commas. The following keywords can be used:" << endl;
  cout << "  forward  (include module into forward model; the processed current estimate is forward-projected)" << endl;
  cout << "  post     (apply module before saving the image; the processed image is not put back as the estimate for the next update)" << endl;
  cout << "  intra    (apply module to the updated image use it as the current estimate for the next update)" << endl;
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageProcessingManager::Initialize()
{
  // Check if parameters have been checked
  if (!m_checked)
  {
    Cerr("***** oImageProcessingManager::Initialize() -> Parameters have not been checked ! Please call CheckParameters() before." << endl);
    return 1;
  }
  // Case with no options (no image processing module)
  if (m_options.size()==0)
  {
    m_initialized = true;
    return 0;
  }
  // Verbose
  if (m_verbose>=1) Cout("oImageProcessingManager::Initialize() -> Initialize image processing modules" << endl);
  // Parse image processing modules options and initialize them
  if (ParseOptionsAndInitializeImageProcessingModules())
  {
    Cerr("***** oImageProcessingManager::Initialize() -> A problem occured while parsing image processing modules options and initializing them !" << endl);
    return 1;
  }
  // All set
  m_initialized = true;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules()
{
  // ===================================================================
  // First get the number of processing modules from the list of options
  // ===================================================================

  m_nbImageProcessingModules = m_options.size();

  // Allocate the tables
  m2p_ImageProcessingModules = (vImageProcessingModule**)malloc(m_nbImageProcessingModules*sizeof(vImageProcessingModule*));
  mp_applyForward = (bool*)malloc(m_nbImageProcessingModules*sizeof(bool));
  mp_applyIntra = (bool*)malloc(m_nbImageProcessingModules*sizeof(bool));
  mp_applyPost = (bool*)malloc(m_nbImageProcessingModules*sizeof(bool));

  // ===================================================================
  // Then we loop over all modules, read options and initialize them
  // ===================================================================

  // This is for the automatic initialization of the processing modules
  typedef vImageProcessingModule *(*maker_image_processing_module) ();
  // Get image processing modules list from addon manager
  std::map <string,maker_image_processing_module> list = sAddonManager::GetInstance()->mp_listOfImageProcessingModules;

  // Start the loop
  for (int c=0; c<m_nbImageProcessingModules; c++)
  {
    // Default initializations
    m2p_ImageProcessingModules[c] = NULL;
    mp_applyForward[c] = false;
    mp_applyIntra[c] = false;
    mp_applyPost[c] = false;

    // ___________________________________________________________________________________
    // Search for a double-colon and isolate the module's options from the 'when' actions

    size_t double_colon = m_options[c].find_first_of("::");

    // Send an error if no double-colon
    if (double_colon==string::npos)
    {
      Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> Wrong syntax in the " << c+1 << "th image processing module !" << endl);
      Cerr("                                                                                    No double-colon \"::\" found." << endl);
      ShowCommonHelp();
      return 1;
    }

    // Separate the two arguments
    string proc_part_options = m_options[c].substr(0,double_colon);
    string when_part_options = m_options[c].substr(double_colon+2);

    // ___________________________________________________________________________________
    // Get the module name in the options and isolate the actual module's options

    // Useful strings
    string module = "";
    string list_options = "";
    string file_options = "";

    // Search for a colon ":", this indicates that a configuration file is provided after the module's name
    size_t colon = proc_part_options.find_first_of(":");
    size_t comma = proc_part_options.find_first_of(",");

    // Case 1: we have a colon
    if (colon!=string::npos)
    {
      // Get the image processing module name before the colon
      module = proc_part_options.substr(0,colon);
      // Get the configuration file after the colon
      file_options = proc_part_options.substr(colon+1);
      // List of options is empty
      list_options = "";
    }
    // Case 2: we have a comma
    else if (comma!=string::npos)
    {
      // Get the image processing module name before the first comma
      module = proc_part_options.substr(0,comma);
      // Get the list of options after the first comma
      list_options = proc_part_options.substr(comma+1);
      // Configuration file is empty
      file_options = "";
    }
    // Case 3: no colon and no comma (a single image processing module name)
    else
    {
      // Get the image processing module name
      module = proc_part_options;
      // List of options is empty
      list_options = "";
      // Build the default configuration file
      file_options = sOutputManager::GetInstance()->GetPathToConfigDir() + "/processing/" + module + ".conf";
    }

    // ___________________________________________________________________________________
    // Read the 'when' actions

    // Loop while commas are found
    while ((comma=when_part_options.find_first_of(",")) != string::npos)
    {
      // Extract the first option
      string option = when_part_options.substr(0,comma);
      // Extract the rest
      when_part_options = when_part_options.substr(comma+1);
      // Check the meaning of the option
      if (option=="forward") {mp_applyForward[c] = true;}
      else if (option=="post") {mp_applyPost[c] = true;}
      else if (option=="intra") {mp_applyIntra[c] = true;}
      else
      {
        Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> Unknown keyword '" << option << "' provided in options list !" << endl);
        ShowCommonHelp();
        return 1;
      }
    }
    // Last option
    if (when_part_options=="forward") {mp_applyForward[c] = true;}
    else if (when_part_options=="post") {mp_applyPost[c] = true;}
    else if (when_part_options=="intra") {mp_applyIntra[c] = true;}
    else
    {
      Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> Unknown keyword '" << when_part_options << "' provided in options list !" << endl);
      ShowCommonHelp();
      return 1;
    }

    // ______________________________________________________________________________
    // Create processing module and call associated functions

    // Create the image processing module
    if (list[module]) m2p_ImageProcessingModules[c] = list[module]();
    else
    {
      Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> Image processing module '" << module << "' does not exist !" << endl);
      sAddonManager::GetInstance()->ShowHelpImageProcessingModule();
      return 1;
    }
    // Set parameters
    m2p_ImageProcessingModules[c]->SetImageDimensionsAndQuantification(mp_ImageDimensionsAndQuantification);
    m2p_ImageProcessingModules[c]->SetVerbose(m_verbose);
    // Provide configuration file if any
    if (file_options!="" && m2p_ImageProcessingModules[c]->ReadConfigurationFile(file_options))
    {
      Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> A problem occured while reading and checking configuration file for image processing module '" << module << "' !" << endl);
      return 1;
    }
    // Provide options if any
    if (list_options!="" && m2p_ImageProcessingModules[c]->ReadOptionsList(list_options))
    {
      Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> A problem occured while parsing and reading options list for image processing module '" << module << "' !" << endl);
      return 1;
    }
    // Check parameters
    if (m2p_ImageProcessingModules[c]->CheckParameters())
    {
      Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> A problem occured while checking parameters for image processing module '" << module << "' !" << endl);
      return 1;
    }
    // Initialize the image processing module
    if (m2p_ImageProcessingModules[c]->Initialize())
    {
      Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> A problem occured while initializing image processing module '" << module << "' !" << endl);
      return 1;
    }
    // Check if processing module is dynamic and if dynamic basis functions are used, then it is not compatible if used inside the reconstruction (all but 'post')
    bool intra_reconstruction = mp_applyForward[c] || mp_applyIntra[c];
    bool condition1 = m2p_ImageProcessingModules[c]->GetAffectTimeDimensionFlag() && !mp_ImageDimensionsAndQuantification->GetTimeStaticFlag();
    bool condition2 = m2p_ImageProcessingModules[c]->GetAffectRespDimensionFlag() && !mp_ImageDimensionsAndQuantification->GetRespStaticFlag();
    bool condition3 = m2p_ImageProcessingModules[c]->GetAffectCardDimensionFlag() && !mp_ImageDimensionsAndQuantification->GetCardStaticFlag();
    if (intra_reconstruction && (condition1 || condition2 || condition3))
    {
      Cerr("***** oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules() -> Cannot use dynamic image processing module '" << module << "' along with dynamic basis functions inside the reconstruction !" << endl);
      return 1;
    }
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageProcessingManager::ApplyProcessingForward(oImageSpace* ap_ImageSpace)
{
  #ifdef CASTOR_DEBUG
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageProcessingManager::ApplyProcessingForward() -> Called while not initialized !" << endl);
    return 1;
  }
  #endif
  // Loop on processing modules
  for (int c=0; c<m_nbImageProcessingModules; c++)
  {
    // Apply it only if asked for
    if (mp_applyForward[c])
    {
      // Verbose
      if (m_verbose>=2) Cout("oImageProcessingrManager::ApplyProcessingForward() -> Apply image processing module " << c+1 << " to forward image" << endl);
      // Get the pointer to the image
      FLTNB**** image = ap_ImageSpace->m4p_forwardImage;
      // Apply convolution
      m2p_ImageProcessingModules[c]->Process(image);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageProcessingManager::ApplyProcessingIntra(oImageSpace* ap_ImageSpace)
{
  #ifdef CASTOR_DEBUG
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageProcessingManager::ApplyProcessingIntra() -> Called while not initialized !" << endl);
    return 1;
  }
  #endif
  // Loop on processing modules
  for (int c=0; c<m_nbImageProcessingModules; c++)
  {
    // Apply it only if asked for
    if (mp_applyIntra[c])
    {
      // Verbose
      if (m_verbose>=2) Cout("oImageProcessingrManager::ApplyProcessingIntra() -> Apply image processing module " << c+1 << " to current image" << endl);
      // Get the pointer to the image
      FLTNB**** image = ap_ImageSpace->m4p_image;
      // Apply processing module
      m2p_ImageProcessingModules[c]->Process(image);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageProcessingManager::ApplyProcessingPost(oImageSpace* ap_ImageSpace)
{
  #ifdef CASTOR_DEBUG
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageProcessingManager::ApplyProcessingPost() -> Called while not initialized !" << endl);
    return 1;
  }
  #endif
  // Loop on processing modules
  for (int c=0; c<m_nbImageProcessingModules; c++)
  {
    // Apply it only if asked for
    if (mp_applyPost[c])
    {
      // Verbose
      if (m_verbose>=2) Cout("oImageProcessingManager::ApplyProcessingPost() -> Apply image processing module " << c+1 << " to output image" << endl);
      // Get the pointer to the output image
      FLTNB**** image = ap_ImageSpace->m4p_outputImage;
      // Apply processing module
      m2p_ImageProcessingModules[c]->Process(image);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
