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
  \ingroup  projector
  \brief    Implementation of class vProjector
*/

#include "vProjector.hh"
#include "vScanner.hh"
#include "vDataFile.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vProjector::vProjector()
{
  // Affect default values
  mp_Scanner = NULL;
  mp_ImageDimensionsAndQuantification = NULL;
  mp_sizeVox[0] = -1.;
  mp_sizeVox[1] = -1.;
  mp_sizeVox[2] = -1.;
  mp_nbVox[0] = -1;
  mp_nbVox[1] = -1;
  mp_nbVox[2] = -1;
  m_nbVoxXY = -1;
  mp_halfFOV[0] = -1.;
  mp_halfFOV[1] = -1.;
  mp_halfFOV[2] = -1.;
  m_sensitivityMode = false;
  m_applyTOF = -1;
  m_TOFnbSigmas = 3.;
  m_applyPOI = false;
  m_compatibleWithSPECTAttenuationCorrection = false;
  m_compatibleWithCompression = false;
  m_verbose = 0;
  m_checked = false;
  m_initialized = false;
  mp_mask = NULL;
  m_hasMask = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vProjector::~vProjector()
{
  ;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
{
  // Check that the parameter is not NULL
  if (ap_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** vProjector::SetImageDimensionsAndQuantification() -> Input image dimensions object is null !" << endl);
    return 1;
  }
  // Affect image dimensions
  mp_sizeVox[0] = ap_ImageDimensionsAndQuantification->GetVoxSizeX();
  mp_sizeVox[1] = ap_ImageDimensionsAndQuantification->GetVoxSizeY();
  mp_sizeVox[2] = ap_ImageDimensionsAndQuantification->GetVoxSizeZ();
  mp_nbVox[0] = ap_ImageDimensionsAndQuantification->GetNbVoxX();
  mp_nbVox[1] = ap_ImageDimensionsAndQuantification->GetNbVoxY();
  mp_nbVox[2] = ap_ImageDimensionsAndQuantification->GetNbVoxZ();
  m_nbVoxXY = mp_nbVox[0] * mp_nbVox[1];
  mp_halfFOV[0] = mp_sizeVox[0] * ((FLTNB)mp_nbVox[0]) / 2.;
  mp_halfFOV[1] = mp_sizeVox[1] * ((FLTNB)mp_nbVox[1]) / 2.;
  mp_halfFOV[2] = mp_sizeVox[2] * ((FLTNB)mp_nbVox[2]) / 2.;
  mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vProjector::ShowCommonHelp()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "------------------------------------------------------------------" << endl;
  cout << "-----  Common options for all projectors" << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout << "Currently a single common option is implemented, relevant only for TOF recon: " << endl;
  cout << "  the number of standard deviations for truncating the TOF distribution." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vProjector::ShowHelp()
{
  // Call the specific help function from the children
  ShowHelpSpecific();
  // Then, say if the child projector in use is compatible with SPECT attenuation correction or not
  if (m_compatibleWithSPECTAttenuationCorrection) cout << "This projector is compatible with SPECT attenuation correction." << endl;
  else cout << "This projector is NOT compatible with SPECT attenuation correction." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::ReadCommonOptionsList(const string& a_optionsList)
{
  // TODO : make this more explicit, here it is assumed that a single specific option can be provided 
  // (number of standard deviations for TOF truncation)
  if (a_optionsList!="")
  {
    FLTNB option[1];
    // Read the option for the number of standard deviations for the truncation of TOF weighting function
    if (ReadStringOption(a_optionsList, option, 1, ",", "Common options"))
    {
      Cerr("***** vProjector::ReadCommonOptionsList() -> Failed to correctly read the list of options !" << endl);
      return 1;
    }
    m_TOFnbSigmas = option[0];
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::CheckParameters()
{
  // Check scanner
  if (mp_Scanner==NULL)
  {
    Cerr("***** vProjector::CheckParameters() -> Please provide a valid scanner object !" << endl);
    return 1;
  }
  // Check image dimensions
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** vProjector::CheckParameters() -> Please provide a valid image dimensions and quantification object !" << endl);
    return 1;
  }
  if (mp_sizeVox[0]<=0. || mp_sizeVox[1]<=0. || mp_sizeVox[2]<=0.)
  {
    Cerr("***** vProjector::CheckParameters() -> One or more voxel sizes is negative or null !" << endl);
    return 1;
  }
  if (mp_nbVox[0]<=0 || mp_nbVox[1]<=0 || mp_nbVox[2]<=0)
  {
    Cerr("***** vProjector::CheckParameters() -> One or more number of voxels is negative or null !" << endl);
    return 1;
  }
  // Check TOF
  if ( m_applyTOF!=USE_TOFPOS && m_applyTOF!=USE_TOFBIN && m_applyTOF!=USE_NOTOF )
  {
    Cerr("***** vProjector::CheckParameters() -> TOF flag is incorrect or not set !" << endl);
    return 1;
  }
  // Check verbose level
  if (m_verbose<0)
  {
    Cerr("***** vProjector::CheckParameters() -> Verbose level is negative !" << endl);
    return 1;
  }
  // Check parameters of the child class
  if (CheckSpecificParameters())
  {
    Cerr("***** vProjector::CheckParameters() -> An error occurred while checking parameters of the child projector class !" << endl);
    return 1;
  }
  // Check the number of sigmas for TOF
  if (m_TOFnbSigmas<=0.)
  {
    Cerr("***** vProjector::CheckParameters() -> TOF number of standard deviations is not > 0 !" << endl);
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

int vProjector::Initialize()
{
  // First check that the parameters has been checked !
  if (!m_checked)
  {
    Cerr("***** vProjector::Initialize() -> Must call the CheckParameters() function before initialization !" << endl);
    return 1;
  }
  // Call the intialize function specific to the children
  if (InitializeSpecific())
  {
    Cerr("***** vProjector::Initialize() -> A problem occured while calling the specific initialization of the child projector !" << endl);
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

INTNB vProjector::EstimateMaxNumberOfVoxelsPerLine()
{
  return mp_ImageDimensionsAndQuantification->GetNbVoxXYZ();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::Project(int a_direction, oProjectionLine* ap_ProjectionLine, uint32_t* ap_index1, uint32_t* ap_index2, int a_nbIndices)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vProjector::Project() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // ---------------------------------------------------------------------------------------
  // First: Get cartesian coordinates from the scanner and average positions if compression.
  // Here, we differentiate the case with no compression (i.e. a_nbIndices==1) from the
  // compression case, because it will avoid to perform useless computation without
  // compression. However, it produces some duplication of code parts as a compromise.
  // ---------------------------------------------------------------------------------------

  // ______________________________________________________________________________________
  // Case1: no compression (i.e. a_nbIndices==1)
  if (a_nbIndices==1)
  {
    // Set indices of the line
    ap_ProjectionLine->SetIndex1(((int)(ap_index1[0])));
    ap_ProjectionLine->SetIndex2(((int)(ap_index2[0])));
    // Get positions and orientations from the scanner (mean depth of interaction intrinsicaly taken into account), taking POI into account if any
    if (mp_Scanner->GetPositionsAndOrientations( ((int)(ap_index1[0])), ((int)(ap_index2[0])),
                                                 ap_ProjectionLine->GetPosition1(), ap_ProjectionLine->GetPosition2(),
                                                 ap_ProjectionLine->GetOrientation1(), ap_ProjectionLine->GetOrientation2(),
                                                 ap_ProjectionLine->GetPOI1(), ap_ProjectionLine->GetPOI2() ))
    {
      Cerr("***** vProjector::Project() -> A problem occured while getting positions and orientations from scanner !" << endl);
      return 1;
    }
  }
  // ______________________________________________________________________________________
  // Case2: compression (i.e. a_nbIndices>1)
  else
  {
    // Set default indices of the line to -1 when compression
    ap_ProjectionLine->SetIndex1(-1);
    ap_ProjectionLine->SetIndex2(-1);
    // Buffer pointers for positions and orientations handled by the projection line
    FLTNB* position1 = ap_ProjectionLine->GetPosition1();
    FLTNB* position2 = ap_ProjectionLine->GetPosition2();
    FLTNB* orientation1 = ap_ProjectionLine->GetOrientation1();
    FLTNB* orientation2 = ap_ProjectionLine->GetOrientation2();
    FLTNB* buffer_position1 = ap_ProjectionLine->GetBufferPosition1();
    FLTNB* buffer_position2 = ap_ProjectionLine->GetBufferPosition2();
    FLTNB* buffer_orientation1 = ap_ProjectionLine->GetBufferOrientation1();
    FLTNB* buffer_orientation2 = ap_ProjectionLine->GetBufferOrientation2();
    // Zero the position and orientation
    for (int i=0; i<3; i++)
    {
      position1[i] = 0.;
      position2[i] = 0.;
      orientation1[i] = 0.;
      orientation2[i] = 0.;
    }
    // Loop on provided indices
    for (int l=0; l<a_nbIndices; l++)
    {
      // Get positions and orientations from scanner (mean depth of interaction intrinsicaly taken into account), taking POI into account if any
      if (mp_Scanner->GetPositionsAndOrientations( ((int)(ap_index1[l])), ((int)(ap_index2[l])),
                                                   buffer_position1, buffer_position2,
                                                   buffer_orientation1, buffer_orientation2,
                                                   ap_ProjectionLine->GetPOI1(), ap_ProjectionLine->GetPOI2() ))
      {
        Cerr("***** vProjector::Project() -> A problem occured while getting positions and orientations from scanner !" << endl);
        return 1;
      }
      // Add those contributions to the mean position and orientation
      for (int i=0; i<3; i++)
      {
        position1[i] += buffer_position1[i];
        position2[i] += buffer_position2[i];
        orientation1[i] += buffer_orientation1[i];
        orientation2[i] += buffer_orientation2[i];
      }
    }
    // Average positions and orientations
    for (int i=0; i<3; i++)
    {
      position1[i] /= ((FLTNB)a_nbIndices);
      position2[i] /= ((FLTNB)a_nbIndices);
      orientation1[i] /= ((FLTNB)a_nbIndices);
      orientation2[i] /= ((FLTNB)a_nbIndices);
    }
  }

  // --------------------------------------------------------------
  // Second: Modify the end points coordinates from common options,
  // random, offset, and LOR displacements.
  // -----------------------------------------------------------

  // Apply common options TODO

  // Apply LORs displacement TODO

  // Apply global image offset
  ap_ProjectionLine->ApplyOffset();

  // Apply bed position offset
  ap_ProjectionLine->ApplyBedOffset();

  // -----------------------------------------------------------
  // Third: project the line
  // -----------------------------------------------------------

  // Compute LOR length
  ap_ProjectionLine->ComputeLineLength();

  // Switch on different TOF options
  switch (m_applyTOF)
  {
    case USE_NOTOF:
      if (ProjectWithoutTOF( a_direction, ap_ProjectionLine ))
      {
        Cerr("***** vProjector::Project() -> A problem occured while projecting a line without time-of-flight !" << endl);
        return 1;
      }
      break;
    case USE_TOFPOS:
      if (ProjectWithTOFPos( a_direction, ap_ProjectionLine ))
      {
        Cerr("***** vProjector::Project() -> A problem occured while projecting a line with time-of-flight position !" << endl);
        return 1;
      }
      break;
    case USE_TOFBIN:
      if (ProjectWithTOFBin( a_direction, ap_ProjectionLine ))
      {
        Cerr("***** vProjector::Project() -> A problem occured while projecting a line with binned time-of-flight !" << endl);
        return 1;
      }
      break;
    // No default
  }

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    Cout("vProjector::Project() -> Exit function" << endl);
  }
  #endif

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
