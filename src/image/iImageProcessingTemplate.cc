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
  \brief    Implementation of class iImageProcessingTemplate
*/

#include "iImageProcessingTemplate.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iImageProcessingTemplate::iImageProcessingTemplate() : vImageProcessingModule()
{
  // Set the booleans that describe if this image processing module will affect each of the dynamic dimension.
  // The boolean members m_affectXXXDimensionFlag are used to avoid any misuse of an image processing module due to the indexation
  // on dynamic basis functions. Here is an example: if the module is acting on the time dynamics (1st pointer) and actually uses
  // the time information to do some processing along this dimension, then the m_affectTimeDimensionFlag must be set to 'true' in
  // the constructor. Doing that, the image processing manager will forbid any use of this processing module if time basis functions
  // are being used. So this insures that when entering this function, the first pointer of the a4p_image will be frames and not
  // generic time basis functions (i.e. there is an equivalence). The same reasoning applies to the respiratory and cardiac dimensions.
  m_affectTimeDimensionFlag = false;
  m_affectRespDimensionFlag = false;
  m_affectCardDimensionFlag = false;
  // Affect default values to the parameters specific to this module
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iImageProcessingTemplate::~iImageProcessingTemplate()
{
  // Delete or free all structures specific to this module that were allocated by this module
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageProcessingTemplate::ReadConfigurationFile(const string& a_configurationFile)
{
  // Implement here the reading of any options specific to this image processing module, through a configuration file

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageProcessingTemplate::ReadOptionsList(const string& a_optionsList)
{
  // Implement here the reading of any options specific to this image processing module, through a list of options separated by commas

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iImageProcessingTemplate::ShowHelp()
{
  // Here, display some help and guidance to how to use this image processing module and what it does.
  // Also describes the form of the configuration file or the list of options that parameterize your module.
  cout << "This image processing module is a template class dedicated to add your own custom image processing module." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageProcessingTemplate::CheckSpecificParameters()
{
  // Implement here all mandatory checks specific to this image processing module needed for a nice use

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageProcessingTemplate::InitializeSpecific()
{
  // Implement here the initialization of whatever member variables specifically used by this image processing module

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageProcessingTemplate::Process(FLTNB**** a4p_image)
{
  // Do the processing here.
  // The input parameter is the image to be processed.
  // For genericity purpose, the a4p_image is indexed using dynamic basis functions and not directly frames and gates indices.
  //   1st pointer: time basis functions
  //   2nd pointer: respiratory basis functions
  //   3rd pointer: cardiac basis functions
  //   4th pointer: voxels
  // The boolean members m_affectXXXDimensionFlag are used to avoid any misuse of an image processing module due to the indexation
  // on dynamic basis functions. Here is an example: if the module is acting on the time dynamics (1st pointer) and actually uses
  // the time information to do some processing along this dimension, then the m_affectTimeDimensionFlag must be set to 'true' in
  // the constructor. Doing that, the image processing manager will forbid any use of this processing module if time basis functions
  // are being used. So this insures that when entering this function, the first pointer of the a4p_image will be frames and not
  // generic time basis functions (i.e. there is an equivalence). The same reasoning applies to the respiratory and cardiac dimensions.
  // All informations about frames, durations, gates, etc, can be accessed via the object member mp_ImageDimensionsAndQuantification.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
