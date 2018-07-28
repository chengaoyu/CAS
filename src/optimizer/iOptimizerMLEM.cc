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
  \brief    Implementation of class iOptimizerMLEM
*/

#include "iOptimizerMLEM.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerMLEM::iOptimizerMLEM() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image for MLEM
  m_nbBackwardImages = 1;
  // MLEM does not accept penalties
  //m_penaltyEnergyFunctionDerivativesOrder = 0;
  // MLEM is compatible with listmode and histogram data
  m_listmodeCompatibility = true;
  m_histogramCompatibility = true;
  // MLEM is compatible with both emission and log-converted transmission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = true;

  // --------------------------
  // Specific member parameters
  // --------------------------

  m_dataSpaceDenominatorThreshold = -1.;
  m_minimumImageUpdateFactor = -1.;
  m_maximumImageUpdateFactor = -1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerMLEM::~iOptimizerMLEM()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerMLEM::ShowHelpSpecific()
{
  cout << "This optimizer is the standard MLEM (for Maximum Likelihood Expectation Maximization)." << endl;
  cout << "If subsets are used, it naturally becomes the OSEM optimizer." << endl;
  cout << "With transmission data, the log-converted pre-corrected data are used as in J. Nuyts et al: \"Iterative reconstruction" << endl;
  cout << "for helical CT: a simulation study\", Phys. Med. Biol., vol. 43, pp. 729-737, 1998." << endl;
  cout << "The following options can be used (in this particular order when provided as a list):" << endl;
  cout << "  initial image value: to set the uniform voxel value for the initial image" << endl;
  cout << "  denominator threshold: to set the threshold of the data space denominator under which the ratio is set to 1" << endl;
  cout << "  minimum image update: to set the minimum of the image update factor under which it stays constant (0 or a negative value" << endl;
  cout << "                        means no minimum thus allowing a 0 update)" << endl;
  cout << "  maximum image update: to set the maximum of the image update factor over which it stays constant (0 or a negative value means" << endl;
  cout << "                        no maximum)" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLEM::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLEM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the denominator threshold option
  key_word = "denominator threshold";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_dataSpaceDenominatorThreshold, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLEM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the minimum image update option
  key_word = "minimum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_minimumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLEM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the maximum image update option
  key_word = "maximum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_maximumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerMLEM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLEM::ReadOptionsList(const string& a_optionsList)
{
  // There are 4 floating point variables as options
  const int nb_options = 4;
  FLTNB options[nb_options];
  
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "MLEM configuration"))
  {
    Cerr("***** iOptimizerMLEM::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_initialValue = options[0];
  m_dataSpaceDenominatorThreshold = options[1];
  m_minimumImageUpdateFactor = options[2];
  m_maximumImageUpdateFactor = options[3];
  
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLEM::CheckSpecificParameters()
{
  // Check that initial image value is strictly positive
  if (m_initialValue<=0.)
  {
    Cerr("***** iOptimizerMLEM->CheckSpecificParameters() -> Provided initial image value (" << m_initialValue << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that denominator threshold value is strictly positive
  if (m_dataSpaceDenominatorThreshold<=0.)
  {
    Cerr("***** iOptimizerMLEM->CheckSpecificParameters() -> Provided data space denominator threshold (" << m_dataSpaceDenominatorThreshold << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that maximum image update factor is higher than the minimum
  if (m_minimumImageUpdateFactor>0. && m_maximumImageUpdateFactor>0. && m_maximumImageUpdateFactor<m_minimumImageUpdateFactor)
  {
    Cerr("***** iOptimizerMLEM->CheckSpecificParameters() -> Provided minimum/maximum (" << m_minimumImageUpdateFactor << "/" << m_maximumImageUpdateFactor << " are inconsistent !" << endl);
    return 1;
  }
  // Cannot deal with list-mode transmission data
  if (m_dataSpec==SPEC_TRANSMISSION && m_dataMode==MODE_LIST)
  {
    Cerr("***** iOptimizerMLEM->CheckSpecificParameters() -> Cannot reconstruct list-mode transmission data !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLEM::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerMLEM::InitializeSpecific() -> Use the MLEM optimizer" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      Cout("  --> Data space denominator threshold: " << m_dataSpaceDenominatorThreshold << endl);
      if (m_minimumImageUpdateFactor>0.) Cout("  --> Minimum image update factor: " << m_minimumImageUpdateFactor << endl);
      else Cerr("!!!!! The minimum update value is not set, if using subsets, voxels could be trapped in 0 value causing some negative bias !" << endl);
      if (m_maximumImageUpdateFactor>0.) Cout("  --> Maximum image update factor: " << m_maximumImageUpdateFactor << endl);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLEM::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
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

int iOptimizerMLEM::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                 FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                 FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Case for emission tomography
  if (m_dataSpec==SPEC_EMISSION)
  {
    // Truncate data to 0 if negative
    if (a_data<0.) a_data = 0.;
    // If the foward model is to close to zero, then ignore this data (set to 1 and not 0 because this line is taken into account in the sensitivity)
    if (a_forwardModel>m_dataSpaceDenominatorThreshold) *ap_backwardValues = a_data / a_forwardModel;
    else *ap_backwardValues = 1.;
  }
  // Case for transmission tomography (equivalent of log-converted MLEM from Nuyts et al 1998)
  else if (m_dataSpec==SPEC_TRANSMISSION)
  {
    // Subtract scatters
    a_data -= a_additiveCorrections;
    a_forwardModel -= a_additiveCorrections;
    // Truncate data and model to 1 if inferior
    if (a_data<1.) a_data = 1.;
    if (a_forwardModel<1.) a_forwardModel = 1.;
    // Log-convert the data and the model
    a_data = log(a_blankValue/a_data);
    a_forwardModel = log(a_blankValue/a_forwardModel);
    // If the foward model is to close to zero, then ignore this data (set to 1 and not 0 because this line is taken into account in the sensitivity)
    if (a_forwardModel>m_dataSpaceDenominatorThreshold) *ap_backwardValues = a_data / a_forwardModel;
    else *ap_backwardValues = 1.;
  }

  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerMLEM::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                  FLTNB a_sensitivity, FLTNB* ap_correctionValues, INTNB a_voxel )
{
  // Compute image update factor
  FLTNB image_update_factor = *ap_correctionValues / a_sensitivity;
  // Apply minimum image update factor
  if ( m_minimumImageUpdateFactor > 0. && image_update_factor < m_minimumImageUpdateFactor ) image_update_factor = m_minimumImageUpdateFactor;
  // Apply maximum image update factor
  if ( m_maximumImageUpdateFactor > 0. && image_update_factor > m_maximumImageUpdateFactor ) image_update_factor = m_maximumImageUpdateFactor;
  // Update image
  *ap_newImageValue = a_currentImageValue * image_update_factor;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
