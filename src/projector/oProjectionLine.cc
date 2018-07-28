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
  \brief    Implementation of class oProjectionLine
*/

#include "oProjectionLine.hh"
#include "sOutputManager.hh"
#include "vProjector.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oProjectionLine::oProjectionLine()
{
  // Default all values
  m_verbose = 0;
  m_checked = false;
  m_initialized = false;
  m_bedOffset = 0.;
  m_nbTOFBins = -1;
  m_currentTOFBin = 0;
  m_TOFResolution = -1.;
  m_TOFMeasurement = 0.;
  mp_POI1 = NULL;
  mp_POI2 = NULL;
  mp_POIResolution = NULL;
  m_useMatchedProjectors = false;
  mp_ForwardProjector = NULL;
  mp_BackwardProjector = NULL;
  m_length = 0.;
  m2p_allocatedNbVoxels = NULL;
  m2p_currentNbVoxels = NULL;
  m3p_voxelIndices = NULL;
  m3p_voxelWeights = NULL;
  mp_ImageDimensionsAndQuantification = NULL;
  m_computationStrategy = -1;
  mp_position1 = NULL;
  mp_position2 = NULL;
  mp_orientation1 = NULL;
  mp_orientation2 = NULL;
  mp_bufferPosition1 = NULL;
  mp_bufferPosition2 = NULL;
  mp_bufferOrientation1 = NULL;
  mp_bufferOrientation2 = NULL;
  m_threadNumber = -1;
  m_multiplicativeCorrection = 1.;
  m_index1 = -1;
  m_index2 = -1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oProjectionLine::~oProjectionLine() 
{
  // Go through the destructor only if the object was initialized
  if (m_initialized)
  {
    // Delete voxels lists
    for (int b=0; b<m_nbTOFBins; b++)
    {
      if (m3p_voxelIndices[FORWARD][b])
      {
        if (m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY) free(m3p_voxelIndices[FORWARD][b]);
        else delete m3p_voxelIndices[FORWARD][b];
      }
      if (m3p_voxelWeights[FORWARD][b])
      {
        if (m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY) free(m3p_voxelWeights[FORWARD][b]);
        else delete m3p_voxelWeights[FORWARD][b];
      }
      if (!m_useMatchedProjectors)
      {
        if (m3p_voxelIndices[BACKWARD][b])
        {
          if (m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY) free(m3p_voxelIndices[BACKWARD][b]);
          else delete m3p_voxelIndices[BACKWARD][b];
        }
        if (m3p_voxelWeights[BACKWARD][b])
        {
          if (m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY) free(m3p_voxelWeights[BACKWARD][b]);
          else delete m3p_voxelWeights[BACKWARD][b];
        }
      }
    }
    // Delete rest
    if (m3p_voxelIndices[FORWARD]) delete[] m3p_voxelIndices[FORWARD];
    if (m3p_voxelWeights[FORWARD]) delete[] m3p_voxelWeights[FORWARD];
    if (!m_useMatchedProjectors)
    {
      if (m3p_voxelIndices[BACKWARD]) delete[] m3p_voxelIndices[BACKWARD];
      if (m3p_voxelWeights[BACKWARD]) delete[] m3p_voxelWeights[BACKWARD];
    }
    if (m2p_allocatedNbVoxels[FORWARD]) delete m2p_allocatedNbVoxels[FORWARD];
    if (m2p_currentNbVoxels[FORWARD]) delete m2p_currentNbVoxels[FORWARD];
    if (!m_useMatchedProjectors)
    {
      if (m2p_allocatedNbVoxels[BACKWARD]) delete m2p_allocatedNbVoxels[BACKWARD];
      if (m2p_currentNbVoxels[BACKWARD]) delete m2p_currentNbVoxels[BACKWARD];
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oProjectionLine::CheckParameters()
{
  // Check mandatory parameters
  if (m_nbTOFBins<=0)
  {
    Cerr("***** oProjectionLine::CheckParameters() -> Forbidden number of TOF bins (" << m_nbTOFBins << ") !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (m_computationStrategy!=IMAGE_COMPUTATION_STRATEGY && m_computationStrategy!=FIXED_LIST_COMPUTATION_STRATEGY && m_computationStrategy!=ADAPTATIVE_LIST_COMPUTATION_STRATEGY)
  {
    Cerr("***** oProjectionLine::CheckParameters() -> Computation strategy incorrectly set !" << endl);
    return 1;
  }
  if (mp_POIResolution==NULL)
  {
    Cerr("***** oProjectionLine::CheckParameters() -> POI resolution not set !" << endl);
    return 1;
  }
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** oProjectionLine::CheckParameters() -> oImageDimensionsAndQuantification not set !" << endl);
    return 1;
  }
  if (m_threadNumber<0)
  {
    Cerr("***** oProjectionLine::CheckParameters() -> The thread number associated to this line is not set !" << endl);
    return 1;
  }
  // End
  m_checked = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oProjectionLine::Initialize()
{
  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** oProjectionLine::Initialize() -> Must call CheckParameters() before Initialize() !" << endl);
    return 1;
  }

  // Verbose
  if (m_verbose>=4) Cout("oProjectionLine::Initialize() -> Initialize this projection line for thread " << m_threadNumber << endl);

  // Default LOR length
  m_length = 0.;

  // Allocate pointers that depend on forward or backward projectors
  m2p_allocatedNbVoxels = new INTNB*[2];
  m2p_currentNbVoxels = new INTNB*[2];
  m3p_voxelIndices = new INTNB**[2];
  m3p_voxelWeights = new FLTNB**[2];

  // Allocate number of voxels and lists for forward projector
  m2p_allocatedNbVoxels[FORWARD] = new INTNB[m_nbTOFBins];
  m2p_currentNbVoxels[FORWARD] = new INTNB[m_nbTOFBins];
  m3p_voxelIndices[FORWARD] = new INTNB*[m_nbTOFBins];
  m3p_voxelWeights[FORWARD] = new FLTNB*[m_nbTOFBins];

  // Allocate number of voxels and lists for backward projector
  if (m_useMatchedProjectors)
  {
    // If one single projector, then simply link pointers
    m2p_allocatedNbVoxels[BACKWARD] = m2p_allocatedNbVoxels[FORWARD];
    m2p_currentNbVoxels[BACKWARD] = m2p_currentNbVoxels[FORWARD];
    m3p_voxelIndices[BACKWARD] = m3p_voxelIndices[FORWARD];
    m3p_voxelWeights[BACKWARD] = m3p_voxelWeights[FORWARD];
  }
  else
  {
    // If separated projectors, then allocate
    m2p_allocatedNbVoxels[BACKWARD] = new INTNB[m_nbTOFBins];
    m2p_currentNbVoxels[BACKWARD] = new INTNB[m_nbTOFBins];
    m3p_voxelIndices[BACKWARD] = new INTNB*[m_nbTOFBins];
    m3p_voxelWeights[BACKWARD] = new FLTNB*[m_nbTOFBins];
  }

  // --------------------------------------------------------------------------------------
  // Image computation strategy
  if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
  {
    // This image computation strategy is not compatible with system matrix projectors
    if (mp_ForwardProjector==NULL || mp_BackwardProjector==NULL)
    {
      Cerr("***** oProjectionLine::Initialize() -> Image computation strategy is not compatible with the use of system matrix projectors !" << endl);
      return 1;
    }
    // Verbose
    if (m_verbose>=5) Cout("  --> Choose the image computation strategy" << endl);
    // The number of voxels is held fixed to the image dimensions
    INTNB nb_voxels = mp_ImageDimensionsAndQuantification->GetNbVoxXYZ();
    // Set numbers and allocate
    for (int b=0; b<m_nbTOFBins; b++)
    {
      // Just to remember the size of the image
      m2p_allocatedNbVoxels[FORWARD][b] = nb_voxels;
      m2p_currentNbVoxels[FORWARD][b] = nb_voxels;
      // Voxel indices are useless with this computation strategy
      m3p_voxelIndices[FORWARD][b] = NULL;
      // Voxel weights
      m3p_voxelWeights[FORWARD][b] = new FLTNB[nb_voxels];
      // Do the same for backward projector if needed
      if (!m_useMatchedProjectors)
      {
        // Just to remember the size of the image
        m2p_allocatedNbVoxels[BACKWARD][b] = nb_voxels;
        m2p_currentNbVoxels[BACKWARD][b] = nb_voxels;
        // Voxel indices are useless with this computation strategy
        m3p_voxelIndices[BACKWARD][b] = NULL;
        // Voxel weights
        m3p_voxelWeights[BACKWARD][b] = new FLTNB[nb_voxels];
      }
    }
  }
  // --------------------------------------------------------------------------------------
  // Fixed list computation strategy
  else if (m_computationStrategy==FIXED_LIST_COMPUTATION_STRATEGY)
  {
    // Verbose
    if (m_verbose>=5) Cout("  --> Choose the fixed list computation strategy" << endl);
    // Set numbers and allocate
    for (int b=0; b<m_nbTOFBins; b++)
    {
      // For the system matrix case
      if (mp_ForwardProjector==NULL)
      {
        // Verbose
        if (m_verbose>=5) Cout("  --> System matrix for forward projection" << endl);
        // The number of allocated voxels is fixed with this strategy
        m2p_allocatedNbVoxels[FORWARD][b] = 0;
        m2p_currentNbVoxels[FORWARD][b] = 0;
        // Voxel indices set to NULL for the moment because the pointer will be copied during the execution
        m3p_voxelIndices[FORWARD][b] = NULL;
        // Voxel weights also set to NULL
        m3p_voxelWeights[FORWARD][b] = NULL;
      }
      // For the projector case
      else
      {
        // The number of allocated voxels is fixed with this strategy
        m2p_allocatedNbVoxels[FORWARD][b] = mp_ForwardProjector->EstimateMaxNumberOfVoxelsPerLine();
        m2p_currentNbVoxels[FORWARD][b] = 0;
        // Verbose
        if (m_verbose>=5) Cout("  --> Projector for forward projection using " << m2p_allocatedNbVoxels[FORWARD][b] << " allocated voxels" << endl);
        // Voxel indices
        m3p_voxelIndices[FORWARD][b] = new INTNB[m2p_allocatedNbVoxels[FORWARD][b]];
        // Voxel weights
        m3p_voxelWeights[FORWARD][b] = new FLTNB[m2p_allocatedNbVoxels[FORWARD][b]];
      }
      // Do the same for backward projector if needed
      if (!m_useMatchedProjectors)
      {
        // For the system matrix case
        if (mp_BackwardProjector==NULL)
        {
          // Verbose
          if (m_verbose>=5) Cout("  --> System matrix for backward projection" << endl);
          // The number of allocated voxels is fixed with this strategy
          m2p_allocatedNbVoxels[BACKWARD][b] = 0;
          m2p_currentNbVoxels[BACKWARD][b] = 0;
          // Voxel indices set to NULL for the moment because the pointer will be copied during the execution
          m3p_voxelIndices[BACKWARD][b] = NULL;
          // Voxel weights also set to NULL
          m3p_voxelWeights[BACKWARD][b] = NULL;
        }
        // For the projector case
        else
        {
          // The number of allocated voxels is fixed with this strategy
          m2p_allocatedNbVoxels[BACKWARD][b] = mp_BackwardProjector->EstimateMaxNumberOfVoxelsPerLine();
          m2p_currentNbVoxels[BACKWARD][b] = 0;
          // Verbose
          if (m_verbose>=5) Cout("  --> Projector for backward projection using " << m2p_allocatedNbVoxels[BACKWARD][b] << " allocated voxels" << endl);
          // Voxel indices
          m3p_voxelIndices[BACKWARD][b] = new INTNB[m2p_allocatedNbVoxels[BACKWARD][b]];
          // Voxel weights
          m3p_voxelWeights[BACKWARD][b] = new FLTNB[m2p_allocatedNbVoxels[BACKWARD][b]];
        }
      }
    }
  }
  // --------------------------------------------------------------------------------------
  // Adaptative list computation strategy
  else if (m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY)
  {
    // This image computation strategy is not compatible with system matrix projectors
    if (mp_ForwardProjector==NULL || mp_BackwardProjector==NULL)
    {
      Cerr("***** oProjectionLine::Initialize() -> Adaptative list computation strategy is not compatible with the use of system matrix projectors !" << endl);
      return 1;
    }
    // The number of voxels is set to the diagonal as a first guess
    INTNB nb_voxels = mp_ImageDimensionsAndQuantification->GetNbVoxDiagonal();
    // Verbose
    if (m_verbose>=5) Cout("  --> Choose the adaptative list computation strategy, starting with " << nb_voxels << " allocated voxels" << endl);
    // Set numbers and allocate
    for (int b=0; b<m_nbTOFBins; b++)
    {
      // The number of allocated voxels is fixed with this strategy
      m2p_allocatedNbVoxels[FORWARD][b] = nb_voxels;
      m2p_currentNbVoxels[FORWARD][b] = 0;
      // Voxel indices
      m3p_voxelIndices[FORWARD][b] = (INTNB*)malloc(nb_voxels*sizeof(INTNB));
      // Voxel weights
      m3p_voxelWeights[FORWARD][b] = (FLTNB*)malloc(nb_voxels*sizeof(FLTNB));
      // Do the same for backward projector if needed
      if (!m_useMatchedProjectors)
      {
        // The number of allocated voxels is fixed with this strategy
        m2p_allocatedNbVoxels[BACKWARD][b] = nb_voxels;
        m2p_currentNbVoxels[BACKWARD][b] = 0;
        // Voxel indices
        m3p_voxelIndices[BACKWARD][b] = (INTNB*)malloc(nb_voxels*sizeof(INTNB));
        // Voxel weights
        m3p_voxelWeights[BACKWARD][b] = (FLTNB*)malloc(nb_voxels*sizeof(FLTNB));
      }
    }
  }

  // Allocate position and orientation of end points (as well as buffers)
  mp_position1 = (FLTNB*)malloc(3*sizeof(FLTNB));
  mp_position2 = (FLTNB*)malloc(3*sizeof(FLTNB));
  mp_orientation1 = (FLTNB*)malloc(3*sizeof(FLTNB));
  mp_orientation2 = (FLTNB*)malloc(3*sizeof(FLTNB));
  mp_bufferPosition1 = (FLTNB*)malloc(3*sizeof(FLTNB));
  mp_bufferPosition2 = (FLTNB*)malloc(3*sizeof(FLTNB));
  mp_bufferOrientation1 = (FLTNB*)malloc(3*sizeof(FLTNB));
  mp_bufferOrientation2 = (FLTNB*)malloc(3*sizeof(FLTNB));

  // End
  m_initialized = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::ComputeLineLength()
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::ComputeLineLength() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10) Cout("oProjectionLine::ComputeLineLength() -> Compute length of the line" << endl);
  #endif

  m_length = sqrt( (mp_position1[0]-mp_position2[0])*(mp_position1[0]-mp_position2[0]) +
                   (mp_position1[1]-mp_position2[1])*(mp_position1[1]-mp_position2[1]) +
                   (mp_position1[2]-mp_position2[2])*(mp_position1[2]-mp_position2[2]) );
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

bool oProjectionLine::NotEmptyLine()
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::NotEmptyLine() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10) Cout("oProjectionLine::NotEmptyLine() -> Look if line is empty" << endl);
  #endif

  // Just look in all TOF bins whether there are some contributing voxels
  for (int b=0; b<m_nbTOFBins; b++) if (m2p_currentNbVoxels[FORWARD][b]!=0) return true;
  return false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::Reset()
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::Reset() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10) Cout("oProjectionLine::Reset() -> Reset buffers of the line" << endl);
  #endif

  // Set the length to zero
  m_length = 0.;
  // --------------------------------------------------------------------------------------
  // Image computation strategy
  if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
  {
    // Loop on the TOF bins
    for (int b=0; b<m_nbTOFBins; b++)
    {
      // Set all voxel weights to zero
      for (int v=0; v<m2p_allocatedNbVoxels[FORWARD][b]; v++) m3p_voxelWeights[FORWARD][b][v] = 0.;
      if (!m_useMatchedProjectors) for (int v=0; v<m2p_allocatedNbVoxels[BACKWARD][b]; v++) m3p_voxelWeights[BACKWARD][b][v] = 0.;
    }
  }
  // --------------------------------------------------------------------------------------
  // List computation strategy
  else if ( m_computationStrategy==FIXED_LIST_COMPUTATION_STRATEGY ||
            m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY )
  {
    // Loop on the TOF bins
    for (int b=0; b<m_nbTOFBins; b++)
    {
      // Only reset the number of current voxels, no need to reset matrices
      m2p_currentNbVoxels[FORWARD][b] = 0;
      if (!m_useMatchedProjectors) m2p_currentNbVoxels[BACKWARD][b] = 0;
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::ApplyOffset()
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::ApplyOffset() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10) Cout("oProjectionLine::ApplyOffset() -> Apply the global offset to the line end points" << endl);
  #endif

  // Offset position 1
  mp_position1[0] += mp_ImageDimensionsAndQuantification->GetOffsetX();
  mp_position1[1] += mp_ImageDimensionsAndQuantification->GetOffsetY();
  mp_position1[2] += mp_ImageDimensionsAndQuantification->GetOffsetZ();
  // Offset position 2
  mp_position2[0] += mp_ImageDimensionsAndQuantification->GetOffsetX();
  mp_position2[1] += mp_ImageDimensionsAndQuantification->GetOffsetY();
  mp_position2[2] += mp_ImageDimensionsAndQuantification->GetOffsetZ();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::ApplyBedOffset()
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::ApplyBedOffset() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10) Cout("oProjectionLine::ApplyBedOffset() -> Apply the bed position offset to the line end points" << endl);
  #endif

  // Offset position 1
  mp_position1[2] += m_bedOffset;
  // Offset position 2
  mp_position2[2] += m_bedOffset;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB oProjectionLine::GetVoxelIndex(int a_direction, int a_TOFBin, INTNB a_voxelInLine)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::GetVoxelIndex() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=12)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("oProjectionLine::GetVoxelIndex() -> Get voxel index of voxel number " << a_voxelInLine << " in TOF bin " << a_TOFBin << " of " << direction << " projector" << endl);
  }
  #endif

  // If image computation strategy, then simply return the voxel index itself
  if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
    return a_voxelInLine;
  // Otherwise, look up for it in the voxel indices table
  else
    return m3p_voxelIndices[a_direction][a_TOFBin][a_voxelInLine];
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::AddVoxelInTOFBin(int a_direction, int a_TOFBin, INTNB a_voxelIndex, FLTNB a_voxelWeight)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::AddVoxelInTOFBin() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=12)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("oProjectionLine::AddVoxelInTOFBin() -> Add voxel index " << a_voxelIndex << " of weight " << a_voxelWeight <<
         " into TOF bin " << a_TOFBin << " of " << direction << " projector" << endl);
  }
  #endif

  if (m_computationStrategy==FIXED_LIST_COMPUTATION_STRATEGY)
  {
    // No buffer overflow checks for this computation strategy
    m3p_voxelIndices[a_direction][a_TOFBin][m2p_currentNbVoxels[a_direction][a_TOFBin]] = a_voxelIndex;
    m3p_voxelWeights[a_direction][a_TOFBin][m2p_currentNbVoxels[a_direction][a_TOFBin]] = a_voxelWeight;
    m2p_currentNbVoxels[a_direction][a_TOFBin]++;
  }
  else if (m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY)
  {
    // Check if the number of contributing voxels is equal to the allocated size
    if (m2p_currentNbVoxels[a_direction][a_TOFBin]==m2p_allocatedNbVoxels[a_direction][a_TOFBin])
    {
      // We realloc for one more voxel
      m2p_allocatedNbVoxels[a_direction][a_TOFBin]++;
      m3p_voxelIndices[a_direction][a_TOFBin] = (INTNB*)realloc(m3p_voxelIndices[a_direction][a_TOFBin],m2p_allocatedNbVoxels[a_direction][a_TOFBin]*sizeof(INTNB));
      m3p_voxelWeights[a_direction][a_TOFBin] = (FLTNB*)realloc(m3p_voxelWeights[a_direction][a_TOFBin],m2p_allocatedNbVoxels[a_direction][a_TOFBin]*sizeof(FLTNB));
    }
    // Then register this voxel contribution
    m3p_voxelIndices[a_direction][a_TOFBin][m2p_currentNbVoxels[a_direction][a_TOFBin]] = a_voxelIndex;
    m3p_voxelWeights[a_direction][a_TOFBin][m2p_currentNbVoxels[a_direction][a_TOFBin]] = a_voxelWeight;
    m2p_currentNbVoxels[a_direction][a_TOFBin]++;
  }
  // Different way of doing for each computation strategy
  else if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
  {
    // Simply add the contribution to this voxel index
    m3p_voxelWeights[a_direction][a_TOFBin][a_voxelIndex] += a_voxelWeight;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::AddVoxelAllTOFBins(int a_direction, INTNB a_voxelIndex, FLTNB a_voxelWeight, HPFLTNB* a_tofWeights, INTNB a_tofBinFirst, INTNB a_tofBinLast)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::AddVoxelAllTOFBins() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=12)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("oProjectionLine::AddVoxelAllTOFBin() -> Add voxel index " << a_voxelIndex << " and weights for all TOF bins of " 
    << direction << " projector" << endl);
  }
  #endif

  if (m_computationStrategy==FIXED_LIST_COMPUTATION_STRATEGY)
  {
    for (INTNB t=a_tofBinFirst; t<=a_tofBinLast; t++)
    {
      // No buffer overflow checks for this computation strategy
      m3p_voxelIndices[a_direction][t][m2p_currentNbVoxels[a_direction][t]] = a_voxelIndex;
      m3p_voxelWeights[a_direction][t][m2p_currentNbVoxels[a_direction][t]] = a_voxelWeight * (FLTNB)(a_tofWeights[t]);
      m2p_currentNbVoxels[a_direction][t]++;
    }
  }
  else if (m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY)
  {
    for (INTNB t=a_tofBinFirst; t<=a_tofBinLast; t++)
    {
      // Check if the number of contributing voxels is equal to the allocated size
      if (m2p_currentNbVoxels[a_direction][t]==m2p_allocatedNbVoxels[a_direction][t])
      {
        // We realloc for one more voxel
        m2p_allocatedNbVoxels[a_direction][t]++;
        m3p_voxelIndices[a_direction][t] = (INTNB*)realloc(m3p_voxelIndices[a_direction][t],m2p_allocatedNbVoxels[a_direction][t]*sizeof(INTNB));
        m3p_voxelWeights[a_direction][t] = (FLTNB*)realloc(m3p_voxelWeights[a_direction][t],m2p_allocatedNbVoxels[a_direction][t]*sizeof(FLTNB));
      }
      // Then register this voxel contribution
      m3p_voxelIndices[a_direction][t][m2p_currentNbVoxels[a_direction][t]] = a_voxelIndex;
      m3p_voxelWeights[a_direction][t][m2p_currentNbVoxels[a_direction][t]] = a_voxelWeight * (FLTNB)(a_tofWeights[t]);
      m2p_currentNbVoxels[a_direction][t]++;
    }
  }
  // Different way of doing for each computation strategy
  else if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
  {
    // Simply add the contribution to this voxel index
    for (INTNB t=a_tofBinFirst; t<=a_tofBinLast; t++)
    {
      m3p_voxelWeights[a_direction][t][a_voxelIndex] += a_voxelWeight * (FLTNB)(a_tofWeights[t]);
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::AddVoxel(int a_direction, INTNB a_voxelIndex, FLTNB a_voxelWeight)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** oProjectionLine::AddVoxel() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=12)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("oProjectionLine::AddVoxel() -> Add voxel index " << a_voxelIndex << " of weight " << a_voxelWeight <<
         " in " << direction << " projector" << endl);
  }
  #endif

  // No TOF here
  int no_tof_bin = 0;

  if (m_computationStrategy==FIXED_LIST_COMPUTATION_STRATEGY)
  {
    // No buffer overflow checks for this computation strategy
    m3p_voxelIndices[a_direction][no_tof_bin][m2p_currentNbVoxels[a_direction][no_tof_bin]] = a_voxelIndex;
    m3p_voxelWeights[a_direction][no_tof_bin][m2p_currentNbVoxels[a_direction][no_tof_bin]] = a_voxelWeight;
    m2p_currentNbVoxels[a_direction][no_tof_bin]++;
  }
  else if (m_computationStrategy==ADAPTATIVE_LIST_COMPUTATION_STRATEGY)
  {
    // Check if the number of contributing voxels is equal to the allocated size
    if (m2p_currentNbVoxels[a_direction][no_tof_bin]==m2p_allocatedNbVoxels[a_direction][no_tof_bin])
    {
      // We realloc for one more voxel
      m2p_allocatedNbVoxels[a_direction][no_tof_bin]++;
      m3p_voxelIndices[a_direction][no_tof_bin] = (INTNB*)realloc(m3p_voxelIndices[a_direction][no_tof_bin],m2p_allocatedNbVoxels[a_direction][no_tof_bin]*sizeof(INTNB));
      m3p_voxelWeights[a_direction][no_tof_bin] = (FLTNB*)realloc(m3p_voxelWeights[a_direction][no_tof_bin],m2p_allocatedNbVoxels[a_direction][no_tof_bin]*sizeof(FLTNB));
    }
    // Then register this voxel contribution
    m3p_voxelIndices[a_direction][no_tof_bin][m2p_currentNbVoxels[a_direction][no_tof_bin]] = a_voxelIndex;
    m3p_voxelWeights[a_direction][no_tof_bin][m2p_currentNbVoxels[a_direction][no_tof_bin]] = a_voxelWeight;
    m2p_currentNbVoxels[a_direction][no_tof_bin]++;
  }
  // Different way of doing for each computation strategy
  else if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
  {
    // Simply add the contribution to this voxel index
    m3p_voxelWeights[a_direction][no_tof_bin][a_voxelIndex] += a_voxelWeight;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB oProjectionLine::ForwardProject(FLTNB* ap_image)
{
  // The forward projection value
  FLTNB value = 0.;
  // If image is NULL, then assume 1 everywhere
  if (ap_image==NULL)
  {
    // Loop over voxels inside the current TOF bin
    for (int v=0; v<m2p_currentNbVoxels[FORWARD][m_currentTOFBin]; v++)
      value += m3p_voxelWeights[FORWARD][m_currentTOFBin][v];
  }
  // Otherwise, project the image
  else
  {
    if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
    {
      // Loop over the whole image
      for (int v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++)
        value += m3p_voxelWeights[FORWARD][m_currentTOFBin][v] * ap_image[v];
    }
    else
    {
      // Loop over voxels inside the current TOF bin
      for (int v=0; v<m2p_currentNbVoxels[FORWARD][m_currentTOFBin]; v++)
        value += m3p_voxelWeights[FORWARD][m_currentTOFBin][v]
               * ap_image[ m3p_voxelIndices[FORWARD][m_currentTOFBin][v] ];
    }
  }
  // Apply multiplicative correction term
  value /= m_multiplicativeCorrection;
  // Return the value
  return value;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB oProjectionLine::ForwardProjectWithSPECTAttenuation(FLTNB* ap_attenuation, FLTNB* ap_image)
{
  #ifdef CASTOR_DEBUG
  // This function cannot be used by definition when using the image computation strategy
  if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
  {
    Cerr("***** oProjectionLine::ForwardProjectWithSPECTAttenuation() -> Cannot be used with an image computation strategy over the projection line !" << endl);
    Exit(EXIT_FAILURE); // Have to exit
  }
  // This function cannot be used with an incompatible projector
  if (mp_ForwardProjector!=NULL && !mp_ForwardProjector->GetCompatibilityWithSPECTAttenuationCorrection())
  {
    Cerr("***** oProjectionLine::ForwardProjectWithSPECTAttenuation() -> Cannot be used with an incompatible projector !" << endl);
    Exit(EXIT_FAILURE); // Have to exit
  }
  // Note that normally, the incompatibilities are checked before launching the reconstruction
  #endif

  // The forward projection value
  FLTNB value = 0.;

  // If image is NULL, then assume 1 everywhere in the emission object
  if (ap_image==NULL)
  {
    // Case without attenuation
    if (ap_attenuation == NULL)
    {
      // Loop over the element in the line
      for (int v=0; v<m2p_currentNbVoxels[FORWARD][m_currentTOFBin]; v++)
      {
        // Update forward projection
        value += m3p_voxelWeights[FORWARD][m_currentTOFBin][v];
      }
    }
    // Case with attenuation
    else
    {
      // Loop over the element in the line
      for (int v=0; v<m2p_currentNbVoxels[FORWARD][m_currentTOFBin]; v++)
      {
        FLTNB atn_sum = 0.0;
        // Compute the attenuation of an element belonging to the line
        for (int w=v; w<m2p_currentNbVoxels[FORWARD][m_currentTOFBin]; w++)
          atn_sum += ap_attenuation[ m3p_voxelIndices[FORWARD][m_currentTOFBin][w] ] * m3p_voxelWeights[FORWARD][m_currentTOFBin][w];
        // Take the inverse exponential to save a division after that (and convert from cm-1 to mm-1)
        atn_sum = std::exp( -atn_sum *0.1 );
        // Update forward projection
        value += m3p_voxelWeights[FORWARD][m_currentTOFBin][v] * atn_sum;
      }
    }
  }
  // Otherwise, project the image
  else
  {
    // Case without attenuation
    if (ap_attenuation == NULL)
    {
      // Loop over the element in the line
      for (int v=0; v<m2p_currentNbVoxels[FORWARD][m_currentTOFBin]; v++)
      {
        // Update projection
        value += m3p_voxelWeights[FORWARD][m_currentTOFBin][v] * ap_image[ m3p_voxelIndices[FORWARD][m_currentTOFBin][v] ];
      }
    }
    // Case with attenuation
    else
    {
      // Loop over the element in the line
      for (int v=0; v<m2p_currentNbVoxels[FORWARD][m_currentTOFBin]; v++)
      {
        FLTNB atn_sum = 0.0;
        // Compute the attenuation of an element belonging to the line
        for (int w=v; w<m2p_currentNbVoxels[FORWARD][m_currentTOFBin]; w++)
          atn_sum += ap_attenuation[ m3p_voxelIndices[FORWARD][m_currentTOFBin][w] ] * m3p_voxelWeights[FORWARD][m_currentTOFBin][w];
        // Take the inverse exponential to save a division after that (and convert from cm-1 to mm-1)
        atn_sum = exp( atn_sum * ((FLTNB)(-0.1)) );
        // Update projection
        value += m3p_voxelWeights[FORWARD][m_currentTOFBin][v] * ap_image[ m3p_voxelIndices[FORWARD][m_currentTOFBin][v] ] * atn_sum;
      }
    }
  }
  // Apply multiplicative correction term
  value /= m_multiplicativeCorrection;
  // Return the value
  return value;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::BackwardProject(FLTNB* ap_image, FLTNB a_value)
{
  // Apply multiplicative correction term
  a_value /= m_multiplicativeCorrection;
  // Image-based strategy
  if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
  {
    // Loop over the whole image
    for (int v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++)
      ap_image[v] += a_value * m3p_voxelWeights[BACKWARD][m_currentTOFBin][v];
  }
  // List-based strategies
  else
  {
    // Loop over voxels inside the current TOF bin
    for (int v=0; v<m2p_currentNbVoxels[BACKWARD][m_currentTOFBin]; v++)
      ap_image[ m3p_voxelIndices[BACKWARD][m_currentTOFBin][v] ] += a_value * m3p_voxelWeights[BACKWARD][m_currentTOFBin][v];
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oProjectionLine::BackwardProjectWithSPECTAttenuation(FLTNB* ap_attenuation, FLTNB* ap_image, FLTNB a_value)
{
  #ifdef CASTOR_DEBUG
  // This function cannot be used by definition when using the image computation strategy
  if (m_computationStrategy==IMAGE_COMPUTATION_STRATEGY)
  {
    Cerr("***** oProjectionLine::BackwardProjectWithSPECTAttenuation() -> Cannot be used with an image computation strategy over the projection line !" << endl);
    Exit(EXIT_FAILURE); // Have to exit
  }
  // This function cannot be used with an incompatible projector
  if (mp_BackwardProjector!=NULL && !mp_BackwardProjector->GetCompatibilityWithSPECTAttenuationCorrection())
  {
    Cerr("***** oProjectionLine::BackwardProjectWithSPECTAttenuation() -> Cannot be used with an incompatible projector !" << endl);
    Exit(EXIT_FAILURE); // Have to exit
  }
  // Note that normally, the incompatibilities are checked before launching the reconstruction
  #endif

  // Apply multiplicative correction term
  a_value /= m_multiplicativeCorrection;

  // Case where attenuation is not provided
  if (ap_attenuation == NULL)
  {
    // Loop over the element in the line
    for (int v=0; v<m2p_currentNbVoxels[BACKWARD][m_currentTOFBin]; v++)
    {
      // Update image
      ap_image[ m3p_voxelIndices[BACKWARD][m_currentTOFBin][v] ] += m3p_voxelWeights[BACKWARD][m_currentTOFBin][v] * a_value;
    }
  }
  // Case where attenuation is provided
  else
  {
    // Loop over the element in the line
    for (int v=0; v<m2p_currentNbVoxels[BACKWARD][m_currentTOFBin]; v++)
    {
      FLTNB atn_sum = 0.0;
      // Compute the attenuation of an element belonging to the line
      for (int w=v; w<m2p_currentNbVoxels[BACKWARD][m_currentTOFBin]; w++)
        atn_sum += ap_attenuation[ m3p_voxelIndices[BACKWARD][m_currentTOFBin][w] ] * m3p_voxelWeights[BACKWARD][m_currentTOFBin][w];
      // Take the inverse exponential to save a division after that (and convert from cm-1 to mm-1)
      atn_sum = exp( atn_sum * ((FLTNB)(-0.1)) );
      // Update image
      ap_image[ m3p_voxelIndices[BACKWARD][m_currentTOFBin][v] ] += m3p_voxelWeights[BACKWARD][m_currentTOFBin][v] * a_value * atn_sum;
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB oProjectionLine::ComputeLineIntegral(int a_direction)
{
  // Will be the cumulative integral
  FLTNB integral = 0.;
  // Simply loop over all voxels weights and sum them
  for (int v=0; v<m2p_currentNbVoxels[a_direction][m_currentTOFBin]; v++)
    integral += m3p_voxelWeights[a_direction][m_currentTOFBin][v];
  // Return the result
  return integral;
}

