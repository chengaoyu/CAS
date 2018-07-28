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
  \brief    Implementation of class iOptimizerOriginalAML
*/

#include "iOptimizerOriginalAML.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerOriginalAML::iOptimizerOriginalAML() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image for AML
  m_nbBackwardImages = 1;
  // AML does not accept penalties
  //m_penaltyEnergyFunctionDerivativesOrder = 0;
  // AML is only compatible with histogram data
  m_listmodeCompatibility = false;
  m_histogramCompatibility = true;
  // AML is only compatible with emission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = false;

  // --------------------------
  // Specific member parameters
  // --------------------------

  // Threshold applied to the denominator of the data space operation to avoid to high ratios
  m_dataSpaceDenominatorThreshold = -1.;
  // This is the AML bound parameter
  m_bound = 1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerOriginalAML::~iOptimizerOriginalAML()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerOriginalAML::ShowHelpSpecific()
{
  cout << "This optimizer is the AML algorithm derived from the AB-EMML of C. Byrne, Inverse Problems, 1998, vol. 14, pp. 1455-67." << endl;
  cout << "The bound B is taken to infinity, so only the bound A can be parameterized." << endl;
  cout << "This bound must be quantitative (same unit as the reconstructed image)." << endl;
  cout << "It is provided as a single value and thus assuming a uniform bound." << endl;
  cout << "Subsets can be used." << endl;
  cout << "With a negative or null bound, this algorithm implements equation 6 of A. Rahmim et al, Phys. Med. Biol., 2012, vol. 57, pp. 733-55." << endl;
  cout << "If a positive bound is provided, then we suppose that the bound A is taken to minus infinity. In that case, this algorithm implements" << endl;
  cout << "equation 22 of K. Van Slambrouck et al, IEEE TMI, Jan 2015, vol. 34, pp. 126-136." << endl;
  cout << "The following options can be used (in this particular order when provided as a list):" << endl;
  cout << "  initial image value: to set the uniform voxel value for the initial image" << endl;
  cout << "  denominator threshold: to set the threshold of the data space denominator under which the ratio is set to 1" << endl;
  cout << "  bound: to set the bound parameter that shift the Poisson law (quantitative, negative or null for standard AML and positive for infinite AML)." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOriginalAML::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOriginalAML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the denominator threshold option
  key_word = "denominator threshold";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_dataSpaceDenominatorThreshold, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOriginalAML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the bound value
  key_word = "bound";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_bound, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOriginalAML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOriginalAML::ReadOptionsList(const string& a_optionsList)
{
  // There are 3 floating point variables as options
  const int nb_options = 3;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "AML configuration"))
  {
    Cerr("***** iOptimizerOriginalAML::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_initialValue = options[0];
  m_dataSpaceDenominatorThreshold = options[1];
  m_bound = options[2];
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOriginalAML::CheckSpecificParameters()
{
  // Check that denominator threshold value is strictly positive
  if (m_dataSpaceDenominatorThreshold<=0.)
  {
    Cerr("***** iOptimizerOriginalAML->CheckSpecificParameters() -> Provided data space denominator threshold (" << m_dataSpaceDenominatorThreshold << ") must be strictly positive !" << endl);
    return 1;
  }
  // If a negative or null bound is provided (standard AML), then check that the initial image value is above the provided bound
  if (m_bound<=0. && m_initialValue<m_bound)
  {
    Cerr("***** iOptimizerOriginalAML::CheckSpecificParameters() -> The initial image value (" << m_initialValue << ") must be higher than or equal to the provided bound value ("
                                                                                  << m_bound << ") !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOriginalAML::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerOriginalAML::InitializeSpecific() -> Use the AML optimizer" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      Cout("  --> Data space denominator threshold: " << m_dataSpaceDenominatorThreshold << endl);
      if (m_bound<=0.) Cout("  --> Bound value: " << m_bound << endl);
      else Cout("  --> Bound value taken to minus infinity" << endl);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOriginalAML::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                          FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                          FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Line weight here is a simply 1
  *ap_weight = 1.;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOriginalAML::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                        FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                        FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Implementation of the standard AML (bound is negative or null)
  if (m_bound<=0.)
  {
    // Compute the bound projection
    FLTNB shift = ForwardProject(ap_Line) * m_bound;
    // Compute denominator (this is the shifted model)
    FLTNB denominator = a_forwardModel - shift;
    // If the denominator is to close to zero, then ignore this data
    if (denominator<=m_dataSpaceDenominatorThreshold) *ap_backwardValues = 1.;
    // Otherwise, do it normally
    else
    {
      // Compute numerator (this is the shifted data)
      FLTNB numerator = a_data - shift;
      // Truncate to 0 if negative
      if (numerator<0.) *ap_backwardValues = 0.;
      // Otherwise, update backward values
      else *ap_backwardValues = numerator / denominator;
    }
  }
  // Implementation of AML with the bound taken to minus infinity (when the provided bound is positive)
  else
  {
    *ap_backwardValues = (a_data - a_forwardModel) / ForwardProject(ap_Line);
  }
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOriginalAML::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                         FLTNB a_sensitivity, FLTNB* ap_correctionValues, INTNB a_voxel )
{
  // Compute image update factor
  FLTNB image_update_factor = *ap_correctionValues / a_sensitivity;
  // Update image for the standard AML (bound is negative or null)
  if (m_bound<=0.) *ap_newImageValue = m_bound + (a_currentImageValue - m_bound) * image_update_factor;
  // Update image for a bound taken to minus infinity (when the provided bound is positive)
  else *ap_newImageValue = a_currentImageValue + image_update_factor;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
