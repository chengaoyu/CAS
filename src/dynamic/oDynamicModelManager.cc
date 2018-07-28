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
  \ingroup  dynamic
  \brief    Implementation of class oDynamicModelManager
*/

#include "oDynamicModelManager.hh"
#include "sAddonManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn oDynamicModelManager
  \brief Constructor of oDynamicModelManager. Simply set all data members to default values.
*/
oDynamicModelManager::oDynamicModelManager()
{
  // Image dimensions
  mp_ID = NULL;
    // Options for each model type
  m_options = "";

  // Model object and associated bool
  mp_DynamicModel = NULL;
  m_UseModel = false;

  // Verbosity
  m_verbose = -1;
  m_checked = false;
  m_initialized = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~oDynamicModelManager
  \brief Destructor of oDynamicModelManager. Free memory from all allocated tabs.
*/
oDynamicModelManager::~oDynamicModelManager() 
{
  // Delete model objects
  if (mp_DynamicModel) delete mp_DynamicModel;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckParameters
  \brief This function is used to check parameters after the latter
         have been all set using Set functions.
  \return 0 if success, positive value otherwise.
*/
int oDynamicModelManager::CheckParameters()
{
  #ifdef CASTOR_VERBOSE
  if (m_verbose>=2) Cout("oDynamicModelManager::CheckParameters() ..."<< endl);
  #endif
  
  // Check image dimensions
  if (mp_ID==NULL)
  {
    Cerr("***** oDynamicModelManager::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** oDynamicModelManager::CheckParameters() -> Wrong verbosity level provided !" << endl);
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
/*
  \fn Initialize
  \brief Set the dynamic model flag and instanciate/initialize model objects
         through the ParseOptionsAndInitializeDeformations() private function.
  \return 0 if success, positive value otherwise.
*/
int oDynamicModelManager::Initialize()
{
  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** oDynamicModelManager::Initialize() -> Must call CheckParameters() before Initialize() !" << endl);
    return 1;
  }
  
  // Check options
  if (m_options=="")
  {
    m_initialized = true;
    m_UseModel = false;
    return 0;
  }

  // We have a model
  m_UseModel = true;

  // Verbose
  if (m_verbose>=1) Cout("oDynamicModelManager::Initialize() -> Initialize models" << endl);
  
  // Parse model options and initialize them
  if (ParseOptionsAndInitializeModel())
  {
    Cerr("***** oDynamicModelManager::Initialize() -> A problem occured while parsing model options and initializing them !" << endl);
    return 1;
  }

  // Normal end
  m_initialized = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ParseOptionsAndInitializeProjectors
  \brief Parse dynamic model options contained in the previously provided
         strings. This function is called inside the Initialize() function.
  \details Manage the options reading and initialize specific vDynamicModel
           Options are a string containing first the name of the model,
           then either a ':' and a configuration file specific to the model
           - or - as many ',' as needed parameters for this model.
           Specific pure virtual functions of the vDynamicModel are used to read parameters and initialize them.
  \return 0 if success, positive value otherwise
*/
int oDynamicModelManager::ParseOptionsAndInitializeModel()
{
  #ifdef CASTOR_VERBOSE
  if (m_verbose>=3) Cout("oDynamicModelManager::ParseOptionsAndInitializeModel ..."<< endl); 
  #endif

  string dynamic_model = "";
  string list_options = "";
  string file_options = "";

  // This is for the automatic initialization of the models
  typedef vDynamicModel *(*maker_dynamic_model) ();

  // Get model's list from addon manager
  std::map <string,maker_dynamic_model> list = sAddonManager::GetInstance()->mp_listOfDynamicModels;

  size_t colon, comma;

  // ---------------------------------------------------------------------------------------------------
  // Manage model for respiratory motion
  // ---------------------------------------------------------------------------------------------------
    
  // First, check if we have dynamic data
  if (mp_ID->GetNbTimeFrames() <= 1 &&
      mp_ID->GetNbRespGates() <= 1 &&
      mp_ID->GetNbCardGates() <= 1)
  {
    Cerr("***** oDynamicModelManager::CheckParameters() -> Dynamic model should be used with more than one time frame/dynamic gate !" << endl);
    return 1;
  }

  // ______________________________________________________________________________
  // Get the model name in the options and isolate the real model's options

  // Search for a colon ":", this indicates that a configuration file is provided after the model name
  colon = m_options.find_first_of(":");
  comma = m_options.find_first_of(",");

  // Case 1: we have a colon
  if (colon!=string::npos)
  {
    // Get the model name before the colon
    dynamic_model = m_options.substr(0,colon);
    // Get the configuration file after the colon
    file_options = m_options.substr(colon+1);
    // List of options is empty
    list_options = "";
  }
  // Case 2: we have a comma
  else if (comma!=string::npos)
  {
    // Get the model name before the first comma
    dynamic_model = m_options.substr(0,comma);
    // Get the list of options after the first comma
    list_options = m_options.substr(comma+1);
    // Configuration file is empty
    file_options = "";
  }
  // Case 3: no colon and no comma (a single model name)
  else
  {
    // Get the model name
    dynamic_model = m_options;
    // Configuration file is empty
    file_options = "";
    // List of options is empty
    list_options = "";
  }

  // Create the model
  if (list[dynamic_model]) mp_DynamicModel = list[dynamic_model]();
  else
  {
    Cerr("***** oDynamicModelManager::ParseOptionsAndInitializeModel() -> Model '" << dynamic_model << "' does not exist !" << endl);
    sAddonManager::GetInstance()->ShowHelpDynamicModel();
    return 1;
  }
  mp_DynamicModel->SetImageDimensionsAndQuantification(mp_ID);
  mp_DynamicModel->SetVerbose(m_verbose);
  // Provide configuration file if any
  if (file_options!="" && mp_DynamicModel->ReadAndCheckConfigurationFile(file_options))
  {
    Cerr("***** oDynamicModelManager::ParseOptionsAndInitializeModel() -> A problem occured while reading and checking frame dynamic model's configuration file !" << endl);
    return 1;
  }
  // Provide options if any
  if (list_options!="" && mp_DynamicModel->ReadAndCheckOptionsList(list_options))
  {
    Cerr("***** oDynamicModelManager::ParseOptionsAndInitializeModel() -> A problem occured while parsing and reading frame dynamic model's options !" << endl);
    return 1;
  }  
  // Check parameters
  if (mp_DynamicModel->CheckParameters())
  {
    Cerr("***** oDynamicModelManager::ParseOptionsAndInitializeModel() -> A problem occured while checking frame dynamic model parameters !" << endl);
    return 1;
  }
  // Initialize the model
  if (mp_DynamicModel->Initialize())
  {
    Cerr("***** oDynamicModelManager::ParseOptionsAndInitializeModel() -> A problem occured while initializing frame dynamic model !" << endl);
    return 1;
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ApplyDynamicModel
  \param ap_ImageS : pointer to the ImageSpace
  \param a_iteration : index of the actual iteration
  \param a_subset : index of the actual subset
  \brief Call successively EstimateModelParameters() ans EstimateImageWithModel() 
         functions of the dynamic model object is 'm_UseModel' is on.
  \return 0 if success, positive value otherwise
*/
int oDynamicModelManager::ApplyDynamicModel(oImageSpace* ap_ImageS, int a_iteration, int a_subset)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oDynamicModelManager::ApplyDynamicModel() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  if (m_UseModel)
  {
    // Verbose
    if(m_verbose>=2) Cout("oDynamicModelManager::ApplyDynamicModel ..."<< endl);

    // Estimate model parameters
    if( mp_DynamicModel->EstimateModelParameters(ap_ImageS, a_iteration, a_subset) )
      {
      Cerr("***** oDynamicModelManager::StepPostProcessInsideSubsetLoop() -> A problem occured while applying dynamic model to current estimate images !" << endl);
      return 1;
    }
    
    // Generate the serie of dynamic images using the model parameters
    if( mp_DynamicModel->EstimateImageWithModel(ap_ImageS, a_iteration, a_subset) )
    {
      Cerr("***** oDynamicModelManager::StepPostProcessInsideSubsetLoop() -> A problem occured while applying dynamic model to current estimate images !" << endl);
      return 1;
    }
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SaveParametricImages
  \param   a_iteration : current iteration index
  \param   a_subset : current number of subsets (or -1 by default)
  \brief Call SaveCoeffImages() function of the dynamic model object is 
        'm_UseModel' is on, in order to save any parameter image
  \return 0 if success, positive value otherwise
*/
int oDynamicModelManager::SaveParametricImages(int a_iteration, int a_subset)
{
  if (m_UseModel)
  {
    // Verbose
    if(m_verbose>=2) Cout("oDynamicModelManager::SaveParametricImages ..."<< endl);
    
    // Write output parametric images if required
    mp_DynamicModel->ComputeOutputParImage();
    
    // Apply Masking if required
    if (mp_DynamicModel->ApplyOutputFOVMaskingOnParametricImages())
    {
      Cerr("***** oDynamicModelManager::SaveParametricImages() -> A problem occured while trying to apply FOV masking on parametric images !" << endl);
      return 1;
    }
    
    // Save images
    if (mp_DynamicModel->SaveParametricImages(a_iteration, a_subset))
    {
      Cerr("***** oDynamicModelManager::SaveParametricImages() -> A problem occured while trying to save image coefficients !" << endl);
      return 1;
    }
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
