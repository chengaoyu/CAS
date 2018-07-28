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
  \ingroup  optimizer
  \brief    Implementation of class iOptimizerMLTR
*/

#include "iOptimizerMLTR.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerMLTR::iOptimizerMLTR() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 0
  m_initialValue = 0.;
  // Only one backward image for MLTR
  m_nbBackwardImages = 1;
  // MLTR does not accept penalties
  //m_penaltyEnergyFunctionDerivativesOrder = 0;
  // MLTR is compatible with histogram data but not listmode
  m_listmodeCompatibility = false;
  m_histogramCompatibility = true;
  // MLTR is only compatible with transmission data
  m_emissionCompatibility = false;
  m_transmissionCompatibility = true;

  // --------------------------
  // Specific member parameters
  // --------------------------

  // The alpha image
  mp_alpha = NULL;
  // The alpha ratio between exterior and interior of the eliptical FOV
  m_alphaRatio = 1.;
  // Relaxation parameters
  m_currentRelaxationFactor = 1.;
  m_initialRelaxationFactor = 1.;
  m_finalRelaxationFactor = 1.;
  // Non-negativity constraint
  m_nonNegativityConstraint = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerMLTR::~iOptimizerMLTR()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerMLTR::ShowHelpSpecific()
{
  cout << "This optimizer is a version of the MLTR algorithm implemented from equation 16 of the paper from K. Van Slambrouck and J. Nuyts:" << endl;
  cout << "\"Reconstruction scheme for accelerated maximum lihelihood reconstruction: the patchwork structure\"," << endl;
  cout << "IEEE Trans. Nucl. Sci., vol. 61, pp. 173-81, 2014." << endl;
  cout << "An additional empiric relaxation factor has been added onto the additive update. Its value onto the first and last updates" << endl;
  cout << "can be parameterized. Its value for all updates in between is computed linearly from these first and last provided values." << endl;
  cout << "Subsets can be used." << endl;
  cout << "The following options can be used (in this particular order when provided as a list):" << endl;
  cout << "  initial image value: to set the uniform voxel value for the initial image" << endl;
  cout << "  alpha ratio: to set the ratio between exterior and interior of the cylindrical FOV alpha values (0 value means 0 inside exterior)" << endl;
  cout << "  initial relaxation factor: to set the empiric multiplicative factor on the additive update used at the first update" << endl;
  cout << "  final relaxation factor: to set the empiric multiplicative factor on the additive update used at the last update" << endl;
  cout << "  non-negativity constraint: 0 if no constraint or 1 in order to apply the constraint during the image update" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLTR::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLTR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the alpha ratio option
  key_word = "alpha ratio";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_alphaRatio, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLTR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the initial relaxation factor option
  key_word = "initial relaxation factor";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialRelaxationFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLTR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the final relaxation factor option
  key_word = "final relaxation factor";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_finalRelaxationFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLTR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the non-negativity constraint option
  key_word = "non-negativity constraint";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_nonNegativityConstraint, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLTR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLTR::ReadOptionsList(const string& a_optionsList)
{
  // There are 5 floating point variables as options
  const int nb_options = 5;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "MLTR configuration"))
  {
    Cerr("***** iOptimizerMLTR::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_initialValue = options[0];
  m_alphaRatio = options[1];
  m_initialRelaxationFactor = options[2];
  m_finalRelaxationFactor = options[3];
  if (options[4]==0.) m_nonNegativityConstraint = false;
  else m_nonNegativityConstraint = true;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLTR::CheckSpecificParameters()
{
  // Check that the alpha ratio is positive
  if (m_alphaRatio<0.)
  {
    Cerr("***** iOptimizerMLTR->CheckSpecificParameters() -> Provided alpha ratio (" << m_alphaRatio << ") must be positive !" << endl);
    return 1;
  }
  // Check that the initial and final relaxation factors are strictly positive
  if (m_initialRelaxationFactor<0.)
  {
    Cerr("***** iOptimizerMLTR->CheckSpecificParameters() -> Provided initial relaxation factor (" << m_initialRelaxationFactor << ") must be positive !" << endl);
    return 1;
  }
  if (m_finalRelaxationFactor<0.)
  {
    Cerr("***** iOptimizerMLTR->CheckSpecificParameters() -> Provided final relaxation factor (" << m_finalRelaxationFactor << ") must be positive !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLTR::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerMLTR::InitializeSpecific() -> Use the MLTR optimizer" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      if (m_initialRelaxationFactor!=1. || m_finalRelaxationFactor!=1.)
      {
        Cout("  --> Initial relaxation factor: " << m_initialRelaxationFactor << endl);
        Cout("  --> Final relaxation factor: " << m_finalRelaxationFactor << endl);
      }
      if (m_nonNegativityConstraint) Cout("  --> Apply a non-negativity constraint during image update" << endl);
    }
  }

  // ===================================================================
  // Step 1: Allocate the alpha image from a miscellaneous image
  mp_alpha = mp_ImageSpace->AllocateMiscellaneousImage();

  // ===================================================================
  // Step 2: Compute interior and exterior alpha values

  // We assume the interior value to be always 1, then we change the exterior value according to the alpha ratio provided
  FLTNB alpha_interior = 1.;
  FLTNB alpha_exterior = m_alphaRatio;
  if (m_alphaRatio!=1. && m_verbose>=VERBOSE_DETAIL) Cout("  --> Use alpha values in interior/exterior of eliptical FOV of " << alpha_interior << "/" << alpha_exterior << endl);

  // ===================================================================
  // Step 3: Set the alpha values inside and outside the cylindrical FOV

  // Precast half the number of voxels over X and Y minus 1 (for efficiency)
  FLTNB flt_base_x = 0.5*((FLTNB)(mp_ImageDimensionsAndQuantification->GetNbVoxX()-1));
  FLTNB flt_base_y = 0.5*((FLTNB)(mp_ImageDimensionsAndQuantification->GetNbVoxY()-1));
  // Compute FOV elipse radius over X and Y, then squared
  FLTNB squared_radius_x = 0.5 * ((FLTNB)(mp_ImageDimensionsAndQuantification->GetNbVoxX())) * mp_ImageDimensionsAndQuantification->GetVoxSizeX();
  squared_radius_x *= squared_radius_x;
  FLTNB squared_radius_y = 0.5 * ((FLTNB)(mp_ImageDimensionsAndQuantification->GetNbVoxY())) * mp_ImageDimensionsAndQuantification->GetVoxSizeY();
  squared_radius_y *= squared_radius_y;
  // We assume that the computation of the distance from the center for a given
  // voxel and comparing it with the output FOV percentage costs more than performing
  // the loops in an inverse order compared to how the image is stored in memory.
  // Thus we begin the loops over X, then Y, then we test and if test passes, we
  // do the remaining loop over Z and over all dynamic dimensions.
  int x;
  #pragma omp parallel for private(x) schedule(guided)
  for (x=0; x<mp_ImageDimensionsAndQuantification->GetNbVoxX(); x++)
  {
    // Compute X distance from image center, then squared
    FLTNB squared_distance_x = (((FLTNB)x)-flt_base_x) * mp_ImageDimensionsAndQuantification->GetVoxSizeX();
    squared_distance_x *= squared_distance_x;
    // Loop over Y
    for (int y=0; y<mp_ImageDimensionsAndQuantification->GetNbVoxY(); y++)
    {
      // Compute Y distance from image center, then squared
      FLTNB squared_distance_y = (((FLTNB)y)-flt_base_y) * mp_ImageDimensionsAndQuantification->GetVoxSizeY();
      squared_distance_y *= squared_distance_y;
      // The alpha value
      FLTNB alpha_value = 0.;
      // Set the alpha value to the interior if INSIDE the elipse
      if ( squared_distance_x/squared_radius_x + squared_distance_y/squared_radius_y <= 1. ) alpha_value = alpha_interior;
      // Else OUTSIDE the elipse
      else alpha_value = alpha_exterior;
      // Loop over Z
      for (int z=0; z<mp_ImageDimensionsAndQuantification->GetNbVoxZ(); z++)
      {
        // Compute global voxel index
        INTNB index = z*mp_ImageDimensionsAndQuantification->GetNbVoxXY() + y*mp_ImageDimensionsAndQuantification->GetNbVoxX() + x;
        // Set the alpha value
        mp_alpha[index] = alpha_value;
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

int iOptimizerMLTR::PostDataUpdateSpecificStep(int a_iteration, int a_nbIterations, int a_subset, int* ap_nbSubsets)
{
  // ===================================================================
  // In this function, we update the value of the relaxation factor to
  // be used in the next image update.
  // ===================================================================

  // Special case if only one update
  if (a_nbIterations==1 && ap_nbSubsets[0]==1)
  {
    // Test if the initial and final relaxation factors differ, then throw an error
    if (m_initialRelaxationFactor != m_finalRelaxationFactor)
    {
      Cerr("***** iOptimizerMLTR::PostDataUpdateSpecificStep() -> Initial and final relaxation differ while there is only one update to do !" << endl);
      return 1;
    }
    // Set the value
    m_currentRelaxationFactor = m_initialRelaxationFactor;
  }

  // Compute the total number of updates minus 1
  int total_number_of_updates_minus_one = -1;
  for (int it=0; it<a_nbIterations; it++) total_number_of_updates_minus_one += ap_nbSubsets[it];
  // Compute the current update index
  int current_update_index = 0;
  for (int it=0; it<a_iteration; it++) current_update_index += ap_nbSubsets[it];
  current_update_index += a_subset;
  // Compute the position on the slope of updates
  FLTNB ratio_of_updates = ((FLTNB)current_update_index) / ((FLTNB)total_number_of_updates_minus_one);

  // Finally update the current relaxation factor to be used at this update, from the initial and final relaxation factors
  m_currentRelaxationFactor = m_initialRelaxationFactor * (1.-ratio_of_updates) + m_finalRelaxationFactor * ratio_of_updates;

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("iOptimizerMLTR::PostDataUpdateSpecificStep() -> Current relaxation factor value: " << m_currentRelaxationFactor << endl);

  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLTR::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                   FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                   FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Line weight here is a projection of a uniform 1 image multiplied by the model
  *ap_weight = ForwardProject(ap_Line, mp_alpha) * a_forwardModel;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLTR::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                 FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                 FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Truncate data to 0 if negative
  if (a_data<0.) a_data = 0.;
  // The correction value is simply the subtraction of the data from the model
  *ap_backwardValues = a_forwardModel - a_data;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLTR::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                  FLTNB a_sensitivity, FLTNB* ap_correctionValues, INTNB a_voxel )
{
  // Compute image update factor
  FLTNB image_update_factor = *ap_correctionValues * mp_alpha[a_voxel] * m_currentRelaxationFactor / a_sensitivity;
  // Update image
  *ap_newImageValue = a_currentImageValue + image_update_factor;
  // Apply non-negativity constraint if asked for
  if (m_nonNegativityConstraint && *ap_newImageValue<0.) *ap_newImageValue = 0.;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
