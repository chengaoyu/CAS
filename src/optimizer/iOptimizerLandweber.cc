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
  \brief    Implementation of class iOptimizerLandweber
*/

#include "iOptimizerLandweber.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerLandweber::iOptimizerLandweber() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 0
  m_initialValue = 0.;
  // Only one backward image for Landweber
  m_nbBackwardImages = 1;
  // Landweber does not accept penalties
  //m_penaltyEnergyFunctionDerivativesOrder = 0;
  // Landweber is not compatible with listmode data
  m_listmodeCompatibility = false;
  // Landweber is compatible with histogram data only
  m_histogramCompatibility = true;
  // Landweber is compatible with both emission and transmission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = true;

  // --------------------------
  // Specific member parameters
  // --------------------------

  // The custom relaxation parameter
  m_relaxationFactor = -1.;
  // Non-negativity constraint
  m_nonNegativityConstraint = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerLandweber::~iOptimizerLandweber()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerLandweber::ShowHelpSpecific()
{
  cout << "This optimizer implements the standard Landweber algorithm for least-squares optimization." << endl;
  cout << "With transmission data, it uses the log-converted model to derive the update." << endl;
  cout << "Be aware that the relaxation parameter is not automatically set, so it often requires some" << endl;
  cout << "trials and errors to find an optimal setting. Also, remember that this algorithm is particularly" << endl;
  cout << "slow to converge." << endl;
  cout << "The following options can be used (in this particular order when provided as a list):" << endl;
  cout << "  initial image value: to set the uniform voxel value for the initial image" << endl;
  cout << "  relaxation factor: to set the relaxation factor applied to the update" << endl;
  cout << "  non-negativity constraint: 0 if no constraint or 1 in order to apply the constraint during the image update" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerLandweber::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerLandweber::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the relaxation factor option
  key_word = "relaxation factor";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_relaxationFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerLandweber::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
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

int iOptimizerLandweber::ReadOptionsList(const string& a_optionsList)
{
  // There are 3 floating point variables as options
  const int nb_options = 3;
  FLTNB options[nb_options];
  
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "Landweber configuration"))
  {
    Cerr("***** iOptimizerLandweber::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_initialValue = options[0];
  m_relaxationFactor = options[1];
  if (options[2]==0.) m_nonNegativityConstraint = false;
  else m_nonNegativityConstraint = true;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerLandweber::CheckSpecificParameters()
{
  // Check that relaxation factor value is strictly positive
  if (m_relaxationFactor<=0.)
  {
    Cerr("***** iOptimizerLandweber->CheckSpecificParameters() -> Provided relaxation factor (" << m_relaxationFactor << ") must be strictly positive !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerLandweber::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerLandweber::InitializeSpecific() -> Use the Landweber algorithm" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      Cout("  --> Relaxation factor: " << m_relaxationFactor << endl);
      if (m_nonNegativityConstraint) Cout("  --> Apply a non-negativity constraint during image update" << endl);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerLandweber::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                        FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                        FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Line weight here is simply 1
  *ap_weight = 1.;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerLandweber::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                      FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                      FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Case for emission tomography
  if (m_dataSpec==SPEC_EMISSION)
  {
    // Simply subtract the model from the data
    *ap_backwardValues = (a_data - a_forwardModel);
  }
  // Case for transmission tomography
  else if (m_dataSpec==SPEC_TRANSMISSION)
  {
    // Subtract scatters
    a_data -= a_additiveCorrections;
    a_forwardModel -= a_additiveCorrections;
    // Truncate both if less than 1
    if (a_data<1.) a_data = 1.;
    if (a_forwardModel<1.) a_forwardModel = 1.;
    *ap_backwardValues = log( a_forwardModel / a_data );
  }
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerLandweber::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                       FLTNB a_sensitivity, FLTNB* ap_correctionValues, INTNB a_voxel )
{
  // Compute image update factor
  FLTNB image_update_factor = *ap_correctionValues * m_relaxationFactor / a_sensitivity;
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
