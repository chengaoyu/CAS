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
  \brief    Implementation of class iOptimizerNEGML
*/

#include "iOptimizerNEGML.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerNEGML::iOptimizerNEGML() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image for NEGML
  m_nbBackwardImages = 1;
  // NEGML does not accept penalties
  //m_penaltyEnergyFunctionDerivativesOrder = 0;
  // NEGML is only compatible with histogram data
  m_listmodeCompatibility = false;
  m_histogramCompatibility = true;
  // NEGML is only compatible with emission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = false;

  // --------------------------
  // Specific member parameters
  // --------------------------

  // This is the NEGML psi parameter
  m_psi = -1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerNEGML::~iOptimizerNEGML()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerNEGML::ShowHelpSpecific()
{
  cout << "This optimizer is the NEGML algorithm from K. Van Slambrouck et al, IEEE TMI, Jan 2015, vol. 34, pp. 126-136." << endl;
  cout << "Subsets can be used. This implementation only consider the psi parameter, but not the alpha image design parameter" << endl;
  cout << "which is supposed to be 1 for all voxels. It implements equation 17 of the reference paper." << endl;
  cout << "The following options can be used (in this particular order when provided as a list):" << endl;
  cout << "  initial image value: to set the uniform voxel value for the initial image" << endl;
  cout << "  psi: to set the psi parameter that sets the transition from Poisson to Gaussian statistics (must be positive)." << endl;
  cout << "       (if set to 0, then it is taken to infinity and implements equation 21 in the reference paper)." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerNEGML::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerNEGML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the psi value option
  key_word = "psi";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_psi, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerNEGML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerNEGML::ReadOptionsList(const string& a_optionsList)
{
  // There are 2 floating point variables as options
  const int nb_options = 2;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "NEGML configuration"))
  {
    Cerr("***** iOptimizerNEGML::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_initialValue = options[0];
  m_psi = options[1];
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerNEGML::CheckSpecificParameters()
{
  // Check that the psi value is positive (a zero value means we take psi to infinity)
  if (m_psi<0.)
  {
    Cerr("***** iOptimizerNEGML::CheckSpecificParameters() -> Provided psi value (" << m_psi << ") must be positive !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerNEGML::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerNEGML::InitializeSpecific() -> Use the NEGML optimizer" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      if (m_psi>0.) Cout("  --> Psi: " << m_psi << endl);
      else if (m_psi==0.) Cout("  --> Psi taken to infinity" << endl);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerNEGML::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                    FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                    FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // With psi taken to infinity, the weight is the forward projection of a uniform 1 image
  if (m_psi==0.) *ap_weight = ForwardProject(ap_Line);
  // Line weight here is a forward projection of alpha image divided by the max between psi and the forward model
  else *ap_weight = ForwardProject(ap_Line) / max( m_psi , a_forwardModel );
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerNEGML::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                  FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                  FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Compute numerator
  FLTNB numerator = a_data - a_forwardModel;
  // Compute denominator that will be strictly positive
  FLTNB denominator = 1.;
  if (m_psi>0.) denominator = max( m_psi , a_forwardModel );
  // Update backward values
  *ap_backwardValues = numerator / denominator;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerNEGML::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                   FLTNB a_sensitivity, FLTNB* ap_correctionValues, INTNB a_voxel )
{
  // Compute image update factor
  FLTNB image_update_factor = *ap_correctionValues / a_sensitivity;
  // Update image
  *ap_newImageValue = a_currentImageValue + image_update_factor;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
