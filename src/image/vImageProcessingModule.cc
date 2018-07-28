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
  \brief    Implementation of class vImageProcessingModule
*/

#include "vImageProcessingModule.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vImageProcessingModule::vImageProcessingModule()
{
  // Affect default values
  mp_ImageDimensionsAndQuantification = NULL;
  m_verbose = 0;
  m_affectTimeDimensionFlag = false;
  m_affectRespDimensionFlag = false;
  m_affectCardDimensionFlag = false;
  m_checked = false;
  m_initialized = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vImageProcessingModule::~vImageProcessingModule()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vImageProcessingModule::CheckParameters()
{
  // Check image dimensions
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** vImageProcessingModule::CheckParameters() -> oImageDimensionsAndQuantification is null !" << endl);
    return 1;
  }
  // Check verbose level
  if (m_verbose<0)
  {
    Cerr("***** vImageProcessingModule::CheckParameters() -> Verbose level is negative !" << endl);
    return 1;
  }
  // Call the CheckSpecificParameters function of the child
  if (CheckSpecificParameters())
  {
    Cerr("***** vImageProcessingModule::CheckParameters() -> A problem occured while checking parameters specific to the image processing module !" << endl);
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

int vImageProcessingModule::Initialize()
{
  // Call the specific initialization function of the child
  if (InitializeSpecific())
  {
    Cerr("***** vImageProcessingModule::Initialize() -> A problem occured while initializing stuff specific to the image processing module !" << endl);
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
