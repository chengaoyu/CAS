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
  \brief    Implementation of class iOptimizerTemplate
*/

#include "iOptimizerTemplate.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerTemplate::iOptimizerTemplate() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // (inherited from vOptimizer)
  // ---------------------------

  // Initial value of the image voxels before iterating
  m_initialValue = 1.;
  // The number of backward images used by the optimizer (most of the time it is 1, but it can be more for special cases)
  m_nbBackwardImages = 1;
  // Specify here if the optimizer is compatible with list-mode and/or histogram data
  m_listmodeCompatibility = true;
  m_histogramCompatibility = true;
  // Specify here if the optimizer is compatible with emission and/or transmission data
  m_emissionCompatibility = false;
  m_transmissionCompatibility = false;

  // --------------------------
  // Specific member parameters
  // --------------------------

  // Set here the default values of the parameters specific to this optimizer

}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerTemplate::~iOptimizerTemplate()
{
  // Delete or free ONLY things that were built into this class
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerTemplate::ShowHelpSpecific()
{
  cout << "This optimizer is only a squeleton template to explain how to add an optimizer into CASToR." << endl;
  cout << "If you want to implement your own optimizer, start from here." << endl;
  cout << "The following options can be used (in this particular order when provided as a list):" << endl;
  cout << "  my parameter: blablabla" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerTemplate::ReadConfigurationFile(const string& a_configurationFile)
{
  // This function is designed to read parameters specific to the optimizer through a configuration file.
  // To do that, use the ReadDataASCIIFile() function, and take a look at other optimizers to see how it is done.
  // Do not check the parameters' values, the CheckSpecificParameters() function is designed to do that.
  // Return 1 if any problem.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerTemplate::ReadOptionsList(const string& a_optionsList)
{
  // This function is designed to read parameters specific to the optimizer through a list of options in a string.
  // To do that, use the ReadStringOption() function, and take a look at other optimizers to see how it is done.
  // Do not check the parameters' values, the CheckSpecificParameters() function is designed to do that.
  // Return 1 if any problem.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerTemplate::CheckSpecificParameters()
{
  // This function is designed to check that all parameters specific to the optimizer have been set properly.
  // Return 1 if any problem.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerTemplate::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=2) Cout("iOptimizerTemplate::InitializeSpecific() -> Use the template optimizer" << endl);

  // This function is designed to initialize everything that should be initialized before being able to launch the iteration.
  // Return 1 if any problem.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerTemplate::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                       FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                       FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // The purpose of this function is to assign the weight of the event. In most cases, it is simply 1 (note that all multiplicative
  // terms belonging to the system matrix are automatically taken into account during the back-projection, explaining why the weight
  // is simply 1 in most cases). The weight should be put into *ap_weight.
  // About the parameters:
  //   a_data: is the data of the event (in number of counts; 1 if list-mode data, and any value for histogram data)
  //   a_forwardModel: is the forward model (also in number of counts); this the forward projection of the image including all system matrix
  //                   related multiplicative terms as well as the additive terms; so it is directly comparable to the data
  //   ap_weight: is where the result should be put in
  //   a_multiplicativeCorrections: is the multiplicative corrections specific to the event; in other words, it does not contain the
  //                                multiplicative factors global to the image, which are given in a separate parameter as a_quantificationFactor
  //   a_additiveCorrections: is the additive terms (in number of counts); e.g. scatters and randoms for PET
  //   a_blankValue: is the blank measurement in transmission tomography
  //   a_quantificationFactor: is the quantification factor associated to the image where the event belongs to
  //   ap_Line: is the projection line associated to the event; in case a forward projection should be perform, it must be passed to the
  //            ForwardProject() function as an argument
  // Finally, note that this function is only called for histogram data. For list-mode data, the sensitivity is computed before launching the
  // the iterations, assuming weights of 1. So if the weights here for the specific optimizer are different than 1, then the optimizer cannot
  // handle list-mode data, so remember to set the m_listmodeCompatibility flag to false in the constructor, in such a case.

  *ap_weight = 1.;

  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerTemplate::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                     FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                     FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // The purpose of this function is to perform the data space operations specific to the optimizer (this function is called for each event)
  // in order to compute the correction terms that need to be back-projected. This value should be put into *ap_backwardValues. In case the
  // optimization requires multiple backward images, one must simply set the number of backward images into the constructor as m_nbBackwardImages.
  // Then, the multiple values to be back-projected can be put in ap_backwardValues[0], ap_backwardValues[1], etc.
  // About the parameters:
  //   a_data: is the data of the event (in number of counts; 1 if list-mode data, and any value for histogram data)
  //   a_forwardModel: is the forward model (also in number of counts); this the forward projection of the image including all system matrix
  //                   related multiplicative terms as well as the additive terms; so it is directly comparable to the data
  //   ap_backwardValues: is where the result should be put in (multiple values is accepted if m_nbBackwardImages>1)
  //   a_multiplicativeCorrections: is the multiplicative corrections specific to the event; in other words, it does not contain the
  //                                multiplicative factors global to the image, which are given in a separate parameter as a_quantificationFactor
  //   a_additiveCorrections: is the additive terms (in number of counts); e.g. scatters and randoms for PET
  //   a_blankValue: is the blank measurement in transmission tomography
  //   a_quantificationFactor: is the quantification factor associated to the image where the event belongs to
  //   ap_Line: is the projection line associated to the event; in case a forward projection should be perform, it must be passed to the
  //            ForwardProject() function as an argument
  // Look at the optimizers already implemented in order to get an idea of how to do it.

  *ap_backwardValues = 1.;

  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerTemplate::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                      FLTNB a_sensitivity, FLTNB* ap_correctionValues, INTNB a_voxel )
{
  // The purpose of this function is to perform the image space operations specific to the optimizer (this function is called for each voxel).
  // In other words, the new voxel value has to be computed based on the previous one, the sensitivity and the correction values which are the
  // back-projected correction terms computed in the DataSpaceSpecificOperations() function (there can be multiple correction values if
  // m_nbBackwardImages>1). Note that being here means that the sensitivity value is not 0, so no need to check it.

  *ap_newImageValue = a_currentImageValue;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
