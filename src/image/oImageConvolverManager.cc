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
  \brief    Implementation of class oImageConvolverManager
*/

#include "oImageConvolverManager.hh"
#include "sAddonManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oImageConvolverManager::oImageConvolverManager()
{
  // Image dimensions
  mp_ImageDimensionsAndQuantification = NULL;
  // Options
  m_options = {};
  // Image convolver objects and associated bool
  m_nbImageConvolvers = 0;
  m2p_ImageConvolvers = NULL;
  mp_applyForward = NULL;
  mp_applyBackward = NULL;
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

oImageConvolverManager::~oImageConvolverManager() 
{
  // Delete object
  if (m2p_ImageConvolvers)
  {
    for (int c=0; c<m_nbImageConvolvers; c++) if (m2p_ImageConvolvers[c]) delete m2p_ImageConvolvers[c];
    free(m2p_ImageConvolvers);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageConvolverManager::CheckParameters()
{
  // Check image dimensions
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** oImageConvolverManager::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** oImageConvolverManager::CheckParameters() -> Wrong verbosity level provided !" << endl);
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

void oImageConvolverManager::ShowCommonHelp()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "------------------------------------------------------------------" << endl;
  cout << "-----  How to use an image convolver" << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout << endl;
  cout << "An image convolver is called through the -conv option. The provided argument describes the convolver to be used, its options," << endl;
  cout << "and when to include it within the algorithm. The syntax of the argument must obey one of the three following options:" << endl;
  cout << "  conv::when (in this case, the default configuration file of the convolver is used to set the options values)" << endl;
  cout << "  conv:file.conf::when (in this case, the provided configuration is used)" << endl;
  cout << "  conv,param1,param2,...::when (in this case, the options values are directly provided in the argument)" << endl;
  cout << "In any case, the description of the options specific to each convolver, their order in the list and their configuration" << endl;
  cout << "files syntax are provided in the specific help of each convolver." << endl;
  cout << "The 'when' parameter is an argument describing when to include the convolver within the algorithm. It is a list of keywords" << endl;
  cout << "separating by commas. The following keywords can be used:" << endl;
  cout << "  forward  (include convolver into forward model; a convolution of the current estimate is forward-projected)" << endl;
  cout << "  backward (include convolver into backward model; a convolution of the correction terms is used for the update)" << endl;
  cout << "  post     (apply convolver before saving the image; the convolved image is not put back as the estimate for the next update)" << endl;
  cout << "  psf      (include both 'forward' and 'backward'; the standard image-based PSF modelling)" << endl;
  cout << "  sieve    (include both 'psf' and 'post'; the standard method of sieve)" << endl;
  cout << "  intra    (apply convolver to the updated image and use it as the current estimate for the next update)" << endl;
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageConvolverManager::Initialize()
{
  // Check if parameters have been checked
  if (!m_checked)
  {
    Cerr("***** oImageConvolverManager::Initialize() -> Parameters have not been checked ! Please call CheckParameters() before." << endl);
    return 1;
  }
  // Case with no options (no convolver)
  if (m_options.size()==0)
  {
    m_initialized = true;
    return 0;
  }
  // Verbose
  if (m_verbose>=1) Cout("oImageConvolverManager::Initialize() -> Initialize image convolvers" << endl);
  // Parse image convolver options and initialize them
  if (ParseOptionsAndInitializeImageConvolvers())
  {
    Cerr("***** oImageConvolverManager::Initialize() -> A problem occured while parsing image convolvers options and initializing them !" << endl);
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

int oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers()
{
  // ==============================================================
  // First get the number of convolvers from the list of options
  // ==============================================================

  m_nbImageConvolvers = m_options.size();

  // Allocate the tables
  m2p_ImageConvolvers = (vImageConvolver**)malloc(m_nbImageConvolvers*sizeof(vImageConvolver*));
  mp_applyForward = (bool*)malloc(m_nbImageConvolvers*sizeof(bool));
  mp_applyBackward = (bool*)malloc(m_nbImageConvolvers*sizeof(bool));
  mp_applyIntra = (bool*)malloc(m_nbImageConvolvers*sizeof(bool));
  mp_applyPost = (bool*)malloc(m_nbImageConvolvers*sizeof(bool));

  // ==============================================================
  // Then we loop over all convolvers, read options and initialize
  // ==============================================================

  // This is for the automatic initialization of the convolvers
  typedef vImageConvolver *(*maker_image_convolver) ();
  // Get image convolver's list from addon manager
  std::map <string,maker_image_convolver> list = sAddonManager::GetInstance()->mp_listOfImageConvolvers;

  // Start the loop
  for (int c=0; c<m_nbImageConvolvers; c++)
  {
    // Default initializations
    m2p_ImageConvolvers[c] = NULL;
    mp_applyForward[c] = false;
    mp_applyBackward[c] = false;
    mp_applyIntra[c] = false;
    mp_applyPost[c] = false;

    // ___________________________________________________________________________________
    // Search for a double-colon and isolate the convolver's options from the 'when' actions

    size_t double_colon = m_options[c].find_first_of("::");

    // Send an error if no double-colon
    if (double_colon==string::npos)
    {
      Cerr("***** oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers() -> Wrong syntax in the " << c+1 << "th image convolver !" << endl);
      Cerr("                                                                            No double-colon \"::\" found." << endl);
      ShowCommonHelp();
      return 1;
    }

    // Separate the two arguments
    string conv_part_options = m_options[c].substr(0,double_colon);
    string when_part_options = m_options[c].substr(double_colon+2);

    // ___________________________________________________________________________________
    // Get the image convolver name in the options and isolate the actual image convolver's options

    // Useful strings
    string convolver = "";
    string list_options = "";
    string file_options = "";

    // Search for a colon ":", this indicates that a configuration file is provided after the image convolver name
    size_t colon = conv_part_options.find_first_of(":");
    size_t comma = conv_part_options.find_first_of(",");

    // Case 1: we have a colon
    if (colon!=string::npos)
    {
      // Get the image convolver name before the colon
      convolver = conv_part_options.substr(0,colon);
      // Get the configuration file after the colon
      file_options = conv_part_options.substr(colon+1);
      // List of options is empty
      list_options = "";
    }
    // Case 2: we have a comma
    else if (comma!=string::npos)
    {
      // Get the image convolver name before the first comma
      convolver = conv_part_options.substr(0,comma);
      // Get the list of options after the first comma
      list_options = conv_part_options.substr(comma+1);
      // Configuration file is empty
      file_options = "";
    }
    // Case 3: no colon and no comma (a single image convolver name)
    else
    {
      // Get the image convolver name
      convolver = conv_part_options;
      // List of options is empty
      list_options = "";
      // Build the default configuration file
      file_options = sOutputManager::GetInstance()->GetPathToConfigDir() + "/convolver/" + convolver + ".conf";
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
      else if (option=="backward") {mp_applyBackward[c] = true;}
      else if (option=="post") {mp_applyPost[c] = true;}
      else if (option=="psf") {mp_applyForward[c] = true; mp_applyBackward[c] = true;}
      else if (option=="sieve") {mp_applyForward[c] = true; mp_applyBackward[c] = true; mp_applyPost[c] = true;}
      else if (option=="intra") {mp_applyIntra[c] = true;}
      else
      {
        Cerr("***** oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers() -> Unknown keyword '" << option << "' provided in options list !" << endl);
        ShowCommonHelp();
        return 1;
      }
    }
    // Last option
    if (when_part_options=="forward") {mp_applyForward[c] = true;}
    else if (when_part_options=="backward") {mp_applyBackward[c] = true;}
    else if (when_part_options=="post") {mp_applyPost[c] = true;}
    else if (when_part_options=="psf") {mp_applyForward[c] = true; mp_applyBackward[c] = true;}
    else if (when_part_options=="sieve") {mp_applyForward[c] = true; mp_applyBackward[c] = true; mp_applyPost[c] = true;}
    else if (when_part_options=="intra") {mp_applyIntra[c] = true;}
    else
    {
      Cerr("***** oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers() -> Unknown keyword '" << when_part_options << "' provided in options list !" << endl);
      ShowCommonHelp();
      return 1;
    }

    // ______________________________________________________________________________
    // Create convolver and call associated functions

    // Create the image convolver
    if (list[convolver]) m2p_ImageConvolvers[c] = list[convolver]();
    else
    {
      Cerr("***** oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers() -> Image convolver '" << convolver << "' does not exist !" << endl);
      sAddonManager::GetInstance()->ShowHelpImageConvolver();
      return 1;
    }
    // Set parameters
    m2p_ImageConvolvers[c]->SetImageDimensionsAndQuantification(mp_ImageDimensionsAndQuantification);
    m2p_ImageConvolvers[c]->SetVerbose(m_verbose);
    // Provide configuration file if any
    if (file_options!="" && m2p_ImageConvolvers[c]->ReadConfigurationFile(file_options))
    {
      Cerr("***** oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers() -> A problem occured while reading and checking configuration file for image convolver '" << convolver << "' !" << endl);
      return 1;
    }
    // Provide options if any
    if (list_options!="" && m2p_ImageConvolvers[c]->ReadOptionsList(list_options))
    {
      Cerr("***** oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers() -> A problem occured while parsing and reading options list for image convolver '" << convolver << "' !" << endl);
      return 1;
    }
    // Check parameters
    if (m2p_ImageConvolvers[c]->CheckParameters())
    {
      Cerr("***** oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers() -> A problem occured while checking parameters for image convolver '" << convolver << "' !" << endl);
      return 1;
    }
    // Initialize the image convolver
    if (m2p_ImageConvolvers[c]->Initialize())
    {
      Cerr("***** oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers() -> A problem occured while initializing image convolver '" << convolver << "' !" << endl);
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

int oImageConvolverManager::ConvolveForward(oImageSpace* ap_ImageSpace)
{
  #ifdef CASTOR_DEBUG
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageConvolverManager::ConvolveForward() -> Called while not initialized !" << endl);
    return 1;
  }
  #endif
  // Loop on convolvers
  for (int c=0; c<m_nbImageConvolvers; c++)
  {
    // Apply it only if asked for
    if (mp_applyForward[c])
    {
      // Verbose
      if (m_verbose>=2)
      {
        if (m_nbImageConvolvers>1) Cout("oImageConvolverManager::ConvolveForward() -> Apply convolution " << c+1 << " to forward image" << endl);
        else Cout("oImageConvolverManager::ConvolveForward() -> Apply convolution to forward image" << endl);
      }
      // Loop on basis functions
      for (int tb=0; tb<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tb++)
      {
        for (int rb=0; rb<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rb++)
        {
          for (int cb=0; cb<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cb++)
          {
            // Get the pointer to the image
            FLTNB* image = ap_ImageSpace->m4p_forwardImage[tb][rb][cb];
            // Apply convolution
            m2p_ImageConvolvers[c]->ApplyConvolution(image);
          }
        }
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

int oImageConvolverManager::ConvolveBackward(oImageSpace* ap_ImageSpace)
{
  #ifdef CASTOR_DEBUG
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageConvolverManager::ConvolveBackward() -> Called while not initialized !" << endl);
    return 1;
  }
  #endif
  // Loop on convolvers
  for (int c=0; c<m_nbImageConvolvers; c++)
  {
    // Apply it only if asked for
    if (mp_applyBackward[c])
    {
      // Verbose
      if (m_verbose>=2)
      {
        if (m_nbImageConvolvers>1)
        {
          if (ap_ImageSpace->IsLoadedSensitivity()) Cout("oImageConvolverManager::ConvolveBackward() -> Apply convolution " << c+1 << " to backward image" << endl);
          else Cout("oImageConvolverManager::ConvolveBackward() -> Apply convolution " << c+1 << " to backward and sensitivity images" << endl);
        }
        else
        {
          if (ap_ImageSpace->IsLoadedSensitivity()) Cout("oImageConvolverManager::ConvolveBackward() -> Apply convolution to backward image" << endl);
          else Cout("oImageConvolverManager::ConvolveBackward() -> Apply convolution to backward and sensitivity images" << endl);
        }
      }
      // Deal with backward images, loop on number of images
      for (int img=0; img<ap_ImageSpace->GetNbBackwardImages(); img++)
      {
        // Loop on basis functions
        for (int tb=0; tb<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tb++)
        {
          for (int rb=0; rb<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rb++)
          {
            for (int cb=0; cb<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cb++)
            {
              // Get the pointer to the image
              int thread_0 = 0;
              FLTNB* image = ap_ImageSpace->m6p_backwardImage[img][thread_0][tb][rb][cb];
              // Apply convolution
              m2p_ImageConvolvers[c]->ApplyConvolutionTranspose(image);
            }
          }
        }
      }
      // Deal with sensitivity images for histogram reconstructions
      if (!ap_ImageSpace->IsLoadedSensitivity())
      {
        // Loop on frames and gates
        for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
        {
          for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr); rg++)
          {
            for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS(); cg++)
            {
              // Get the pointer to the image
              int thread_0 = 0;
              FLTNB* image = ap_ImageSpace->m5p_sensitivity[thread_0][fr][rg][cg];
              // Apply convolution
              m2p_ImageConvolvers[c]->ApplyConvolutionTranspose(image);
            }
          }
        }
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

int oImageConvolverManager::ConvolveIntra(oImageSpace* ap_ImageSpace)
{
  #ifdef CASTOR_DEBUG
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageConvolverManager::ConvolveIntra() -> Called while not initialized !" << endl);
    return 1;
  }
  #endif
  // Loop on convolvers
  for (int c=0; c<m_nbImageConvolvers; c++)
  {
    // Apply it only if asked for
    if (mp_applyIntra[c])
    {
      // Verbose
      if (m_verbose>=2)
      {
        if (m_nbImageConvolvers>1) Cout("oImageConvolverManager::ConvolveIntra() -> Apply convolution " << c+1 << " to current image" << endl);
        else Cout("oImageConvolverManager::ConvolveIntra() -> Apply convolution to current image" << endl);
      }
      // Loop on basis functions
      for (int tb=0; tb<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tb++)
      {
        for (int rb=0; rb<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rb++)
        {
          for (int cb=0; cb<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cb++)
          {
            // Get the pointer to the image
            FLTNB* image = ap_ImageSpace->m4p_image[tb][rb][cb];
            // Apply convolution
            m2p_ImageConvolvers[c]->ApplyConvolution(image);
          }
        }
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

int oImageConvolverManager::ConvolveSensitivity(oImageSpace* ap_ImageSpace)
{
  #ifdef CASTOR_DEBUG
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageConvolverManager::ConvolveSensitivity() -> Called while not initialized !" << endl);
    return 1;
  }
  #endif
  // Loop on convolvers
  for (int c=0; c<m_nbImageConvolvers; c++)
  {
    // Apply it only if asked for
    if (mp_applyBackward[c])
    {
      // Verbose
      if (m_verbose>=2)
      {
        if (m_nbImageConvolvers>1) Cout("oImageConvolverManager::ConvolveSensitivity() -> Apply convolution " << c+1 << " to sensitivity image" << endl);
        else Cout("oImageConvolverManager::ConvolveSensitivity() -> Apply convolution to sensitivity image" << endl);
      }
      // Loop on frames and gates
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
      {
        for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr); rg++)
        {
          for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS(); cg++)
          {
            // Get the pointer to the image
            int thread_0 = 0;
            FLTNB* image = ap_ImageSpace->m5p_sensitivity[thread_0][fr][rg][cg];
            // Apply convolution
            m2p_ImageConvolvers[c]->ApplyConvolutionTranspose(image);
          }
        }
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

int oImageConvolverManager::ConvolvePost(oImageSpace* ap_ImageSpace)
{
  #ifdef CASTOR_DEBUG
  // Check if initialized
  if (!m_initialized)
  {
    Cerr("***** oImageConvolverManager::ConvolvePost() -> Called while not initialized !" << endl);
    return 1;
  }
  #endif
  // Loop on convolvers
  for (int c=0; c<m_nbImageConvolvers; c++)
  {
    // Apply it only if asked for
    if (mp_applyPost[c])
    {
      // Verbose
      if (m_verbose>=2)
      {
        if (m_nbImageConvolvers>1) Cout("oImageConvolverManager::ConvolvePost() -> Apply convolution " << c+1 << " to output image" << endl);
        else Cout("oImageConvolverManager::ConvolvePost() -> Apply convolution to output image" << endl);
      }
          
      // Loop on frames and respiratory/cardiac gates
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
      {
        for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
        {
          for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
          {
            // Get the pointer to the output image
            FLTNB* image = ap_ImageSpace->m4p_outputImage[fr][rg][cg];
            // Apply convolution
            m2p_ImageConvolvers[c]->ApplyConvolution(image);
          }
        }
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
