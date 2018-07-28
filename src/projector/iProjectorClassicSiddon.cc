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
  \brief    Implementation of class iProjectorClassicSiddon
*/

#include <algorithm>
#include <cstring>
#include <cassert>

#include "iProjectorClassicSiddon.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorClassicSiddon::iProjectorClassicSiddon() : vProjector()
{
  // This projector is compatible with SPECT attenuation correction because
  // all voxels contributing to a given line are ordered from point 1 (focal)
  // to point 2 (scanner)
  m_compatibleWithSPECTAttenuationCorrection = true;
  // This projector is compatible with compression as it works only with the
  // cartesian coordinates of 2 end points of a line
  m_compatibleWithCompression = true;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorClassicSiddon::~iProjectorClassicSiddon()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorClassicSiddon::ReadConfigurationFile(const string& a_configurationFile)
{
  // No options for classic siddon
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorClassicSiddon::ReadOptionsList(const string& a_optionsList)
{
  // No options for classic siddon
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iProjectorClassicSiddon::ShowHelpSpecific()
{
  cout << "This projector is a simple line projector that computes the exact path length of a line through the voxels." << endl;
  cout << "It is implemented from the following published paper:" << endl;
  cout << "R. L. Siddon, \"Fast calculation of the exact radiological path for a three-dimensional CT array\", Med. Phys., vol. 12, pp. 252-5, 1985." << endl;
  cout << "All voxels are correctly ordered in the line, so this projector can be used with SPECT attenuation correction." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorClassicSiddon::CheckSpecificParameters()
{
  // Nothing to check for this projector
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorClassicSiddon::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=2) Cout("iProjectorClassicSiddon::InitializeSpecific() -> Use classic Siddon projector" << endl);
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB iProjectorClassicSiddon::EstimateMaxNumberOfVoxelsPerLine()
{
  // Find the maximum number of voxels along a given dimension
  INTNB max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxX();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxY()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxY();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxZ();
  // We should have at most 4 voxels in a given plane, so multiply by 4
  // (note: this is not true however it ensures no overflow and is already quite optimized for RAM usage !)
  max_nb_voxels_in_dimension *= 4;
  // Return the value
  return max_nb_voxels_in_dimension;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorClassicSiddon::ProjectWithoutTOF(int a_direction, oProjectionLine* ap_ProjectionLine )
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorClassicSiddon::ProjectWithoutTOF() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorClassicSiddon::Project without TOF -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Get event positions
  FLTNB* event1_position = ap_ProjectionLine->GetPosition1();
  FLTNB* event2_position = ap_ProjectionLine->GetPosition2();

  FLTNB event1[3] = { event1_position[0], event1_position[1], event1_position[2] };
  FLTNB event2[3] = { event2_position[0], event2_position[1], event2_position[2] };

  // **************************************
  // STEP 1: LOR length calculation
  // **************************************
  FLTNB length_LOR = ap_ProjectionLine->GetLength();

  // **************************************
  // STEP 2: Compute entrance voxel indexes
  // **************************************
  FLTNB alphaFirst[] = { 0.0, 0.0, 0.0 };
  FLTNB alphaLast[] = { 0.0, 0.0, 0.0 };

  FLTNB alphaMin = 0.0, alphaMax = 1.0;
  FLTNB delta_pos[] = {
    event2[ 0 ] - event1[ 0 ],
    event2[ 1 ] - event1[ 1 ],
    event2[ 2 ] - event1[ 2 ]
  };

  // Computation of alphaMin et alphaMax values (entrance and exit point of the ray)
  for( int i = 0; i < 3; ++i )
  {
    if( delta_pos[ i ] != 0.0 )
    {
      alphaFirst[i] = (-mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
      alphaLast[i]  = ( mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
      alphaMin = (std::max)(alphaMin,(std::min)(alphaFirst[i],alphaLast[i]));
      alphaMax = (std::min)(alphaMax,(std::max)(alphaFirst[i],alphaLast[i]));
    }
  }

  // if alphaMax is less than or equal to alphaMin no intersection
  // and return an empty buffer
  if( alphaMax <= alphaMin ) return 0;

  // Now we have to find the indices of the particular plane
  // (iMin,iMax), (jMin,jMax), (kMin,kMax)
  int iMin = 0, iMax = 0;
  int jMin = 0, jMax = 0;
  int kMin = 0, kMax = 0;

  // For the X-axis
  if( delta_pos[ 0 ] > 0.0 )
  {
    iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMin * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
    iMax = ::floor( 1 + ( event1[ 0 ] + alphaMax * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
  }
  else if( delta_pos[ 0 ] < 0.0 )
  {
    iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMax * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
    iMax = ::floor( 1 + ( event1[ 0 ] + alphaMin * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
  }
  if( delta_pos[ 0 ] == 0. )
  {
    iMin = 1, iMax = 0;
  }

  // For the Y-axis
  if( delta_pos[ 1 ] > 0. )
  {
    jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMin * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
    jMax = ::floor( 1 + ( event1[ 1 ] + alphaMax * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
  }
  else if( delta_pos[ 1 ] < 0. )
  {
    jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMax * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
    jMax = ::floor( 1 + ( event1[ 1 ] + alphaMin * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
  }
  else if( delta_pos[ 1 ] == 0. )
  {
    jMin = 1, jMax = 0;
  }

  // For the Z-axis
  if( delta_pos[ 2 ] > 0. )
  {
    kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMin * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
    kMax = ::floor( 1 + ( event1[ 2 ] + alphaMax * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
  }
  else if( delta_pos[ 2 ] < 0. )
  {
    kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMax * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
    kMax = ::floor( 1 + ( event1[ 2 ] + alphaMin * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
  }
  else if( delta_pos[ 2 ] == 0. )
  {
    kMin = 1, kMax = 0;
  }

  // Computing the last term n number of intersection
  int n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
    + ( kMax - kMin + 1 );

  // We create a buffer storing the merging data
  // We merge alphaMin, alphaMax, alphaX, alphaY and alphaZ
  FLTNB *alpha = new FLTNB[ mp_nbVox[ 0 ] + mp_nbVox[ 1 ] + mp_nbVox[ 2 ] + 3 ];
  FLTNB *tmpAlpha = new FLTNB[ mp_nbVox[ 0 ] + mp_nbVox[ 1 ] + mp_nbVox[ 2 ] + 3 ];
  FLTNB *alphaX = new FLTNB[ ( mp_nbVox[ 0 ] + 1 ) ];
  FLTNB *alphaY = new FLTNB[ ( mp_nbVox[ 1 ] + 1 ) ];
  FLTNB *alphaZ = new FLTNB[ ( mp_nbVox[ 2 ] + 1 ) ];

  INTNB iElement = iMax - iMin + 1;
  if( iElement > 0 )
  {
    FLTNB *idx = alphaX;
    if( delta_pos[ 0 ] > 0. )
    {
      for( int i = iMin; i <= iMax; ++i )
        *idx++ = ( ( (-mp_halfFOV[ 0 ]) + ( i - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    }
    else if( delta_pos[ 0 ] < 0. )
    {
      for( int i = iMax; i >= iMin; --i )
        *idx++ = ( ( (-mp_halfFOV[ 0 ]) + ( i - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    }
  }

  // For alphaY
  INTNB jElement = jMax - jMin + 1;
  if( jElement > 0 )
  {
    FLTNB *idx = alphaY;
    if( delta_pos[ 1 ] > 0. )
    {
      for( int j = jMin; j <= jMax; ++j )
        *idx++ = ( ( (-mp_halfFOV[ 1 ]) + ( j - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    }
    else if( delta_pos[ 1 ] < 0. )
    {
      for( int j = jMax; j >= jMin; --j )
        *idx++ = ( ( (-mp_halfFOV[ 1 ]) + ( j - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    }
  }

  // For alphaZ
  INTNB kElement = kMax - kMin + 1;
  if( kElement > 0 )
  {
    FLTNB *idx = alphaZ;
    if( delta_pos[ 2 ] > 0. )
    {
      for( int k = kMin; k <= kMax; ++k )
        *idx++ = ( ( (-mp_halfFOV[ 2 ]) + ( k - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    }
    else if( delta_pos[ 2 ] < 0. )
    {
      for( int k = kMax; k >= kMin; --k )
        *idx++ = ( ( (-mp_halfFOV[ 2 ]) + ( k - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    }
  }

  std::merge(
    alphaX, alphaX + iElement,
    tmpAlpha, tmpAlpha,
    alpha );

  ::memcpy( tmpAlpha, alpha, ( iElement ) * sizeof( FLTNB ) );

  std::merge(
    alphaY, alphaY + jElement,
    tmpAlpha, tmpAlpha + iElement,
    alpha );

  ::memcpy( tmpAlpha, alpha, ( iElement + jElement ) * sizeof( FLTNB ) );

  std::merge(
    alphaZ, alphaZ + kElement,
    tmpAlpha, tmpAlpha + iElement + jElement,
    alpha );

  // Computing some constants
  FLTNB const i1 = ( event1[ 0 ] - (-mp_halfFOV[ 0 ]) + mp_sizeVox[0] ) / mp_sizeVox[0];
  FLTNB const i2 = delta_pos[ 0 ] / mp_sizeVox[0];
  FLTNB const j1 = ( event1[ 1 ] - (-mp_halfFOV[ 1 ]) + mp_sizeVox[1] ) / mp_sizeVox[1];
  FLTNB const j2 = delta_pos[ 1 ] / mp_sizeVox[1];
  FLTNB const k1 = ( event1[ 2 ] - (-mp_halfFOV[ 2 ]) + mp_sizeVox[2] ) / mp_sizeVox[2];
  FLTNB const k2 = delta_pos[ 2 ] / mp_sizeVox[2];

  // Computing the index of the voxels
  FLTNB alphaMid = 0.0;
  INTNB i = 0, j = 0, k = 0; // indices of the voxel
  FLTNB *pAlpha = &alpha[ 1 ];
  FLTNB *pAlphaPrevious = &alpha[ 0 ];
  FLTNB coeff = 0.0;
  INTNB numVox = 0;
  // Loop over the number of crossed planes
  for( int nP = 0; nP < n - 1; ++nP, ++pAlpha, ++pAlphaPrevious )
  {
    alphaMid = ( *pAlpha + *pAlphaPrevious ) * 0.5;

    i = alphaMid * i2 + i1;
    if( i < 1 || i > mp_nbVox[ 0 ] ) continue;

    j = alphaMid * j2 + j1;
    if( j < 1 || j > mp_nbVox[ 1 ] ) continue;

    k = alphaMid * k2 + k1;
    if( k < 1 || k > mp_nbVox[ 2 ] ) continue; 

    numVox = ( i - 1 ) + ( j - 1 ) * mp_nbVox[0] + ( k - 1 ) * mp_nbVox[0] * mp_nbVox[1];

    // if this voxel is masked, skip
    if (m_hasMask && !mp_mask[numVox]) continue;

    coeff = length_LOR * ( *pAlpha - *pAlphaPrevious );

    ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff);
  }

  delete[] alpha;
  delete[] tmpAlpha;
  delete[] alphaX;
  delete[] alphaY;
  delete[] alphaZ;

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    Cout("iProjectorClassicSiddon::Project without TOF -> Exit function" << endl);
  }
  #endif

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorClassicSiddon::ProjectWithTOFPos(int a_direction, oProjectionLine* ap_ProjectionLine)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorClassicSiddon::ProjectWithTOFPos() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorClassicSiddon::Project with TOF measurement -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Simpler way now, hopefully it works
  FLTNB* event1_position = ap_ProjectionLine->GetPosition1();
  FLTNB* event2_position = ap_ProjectionLine->GetPosition2();

  FLTNB event1[3] = { event1_position[0], event1_position[1], event1_position[2] };
  FLTNB event2[3] = { event2_position[0], event2_position[1], event2_position[2] };

  // **************************************
  // STEP 1: LOR length calculation
  // **************************************
  FLTNB length_LOR = ap_ProjectionLine->GetLength();

  // **************************************
  // STEP 2: Compute entrance voxel indexes
  // **************************************
  FLTNB alphaFirst[] = { 0.0, 0.0, 0.0 };
  FLTNB alphaLast[] = { 0.0, 0.0, 0.0 };

  FLTNB alphaMin = 0.0, alphaMax = 1.0;
  FLTNB delta_pos[] = {
    event2[ 0 ] - event1[ 0 ],
    event2[ 1 ] - event1[ 1 ],
    event2[ 2 ] - event1[ 2 ]
  };

  // Get TOF info
  HPFLTNB tof_resolutionT = ap_ProjectionLine->GetTOFResolution();
  HPFLTNB tof_measurement = ap_ProjectionLine->GetTOFMeasurement();

  // TOF standard deviation and truncation
  HPFLTNB tof_sigma = tof_resolutionT * SPEED_OF_LIGHT * 0.5 / TWO_SQRT_TWO_LN_2;
  HPFLTNB tof_sigma_sqrt2 = sqrt(2.)*tof_sigma;
  HPFLTNB tof_half_span = tof_sigma * m_TOFnbSigmas;

  // convert delta time into delta length
  HPFLTNB tof_delta = tof_measurement * SPEED_OF_LIGHT * 0.5;

  // distance between the first event1 and the center of the Gaussian distribution along the LOR
  HPFLTNB lor_tof_center = length_LOR * 0.5 + tof_delta;

  // index along each axis of the first voxel falling inside the truncated Gaussian distribution
  HPFLTNB tof_edge_low[] = {0,0,0};
  // index along each axis of the last voxel falling inside the truncated Gaussian distribution
  HPFLTNB tof_edge_high[] = {0,0,0};
  HPFLTNB tof_center;
  INTNB tof_index;

  // low/high voxel edges (in absolute coordinates) for truncated TOF
  for (int ax=0;ax<3;ax++)
  {
    // absolute coordinate along each axis of the center of the TOF distribution
    tof_center = event1[ax] +  lor_tof_center * delta_pos[ax] / length_LOR;

    // absolute coordinate along each axis of the lowest voxel edge spanned by the TOF distribution, limited by the lowest FOV edge
    tof_edge_low[ax] = tof_center - tof_half_span * fabs(delta_pos[ax]) / length_LOR;
    tof_index = max( (INTNB)::floor( (tof_edge_low[ax] - (-mp_halfFOV[ax])) / mp_sizeVox[ax] ), (INTNB)0);
    // if low TOF edge above the highest FOV edge, return empty line
    if (tof_index>mp_nbVox[ax]-1) return 0;
    tof_edge_low[ax] = (HPFLTNB)tof_index *  mp_sizeVox[ax] - mp_halfFOV[ax];

    // absolute coordinate along each axis of the highest voxel edge spanned by the TOF distribution, limited by the highest FOV edge
    tof_edge_high[ax] = tof_center + tof_half_span * fabs(delta_pos[ax]) / length_LOR;
    tof_index = min( (INTNB)::floor( (tof_edge_high[ax] - (-mp_halfFOV[ax])) / mp_sizeVox[ax] ), mp_nbVox[ax]-1);
    // if high TOF edge below the lowest FOV edge, return empty line
    if (tof_index<0) return 0;
    tof_edge_high[ax] = (HPFLTNB)(tof_index+1) * mp_sizeVox[ax] - mp_halfFOV[ax];
  }

  // Computation of alphaMin et alphaMax values 
  // entrance and exit point of the ray at voxel grid (edges), with respect to the truncated TOF distribution, limited by the FOV,
  // absolute normalized distance from event1
  for( int i = 0; i < 3; ++i )
  {
    if( delta_pos[ i ] != 0.0 )
    {
      alphaFirst[i] = (tof_edge_low[i] - event1[i]) / (delta_pos[ i ]);
      alphaLast[i]  = (tof_edge_high[i] - event1[i]) / (delta_pos[ i ]);
      alphaMin = (std::max)(alphaMin,(std::min)(alphaFirst[i],alphaLast[i]));
      alphaMax = (std::min)(alphaMax,(std::max)(alphaFirst[i],alphaLast[i]));
    }
  }

  // if alphaMax is less than or equal to alphaMin no intersection
  // and return an empty buffer
  if( alphaMax <= alphaMin ) return 0;

  // Min/max indices of voxels intersected by the LOR along each axis, 0-N indices (same order as absolute coordinates) match the FOV
  // (iMin,iMax), (jMin,jMax), (kMin,kMax)
  int iMin = 0, iMax = 0;
  int jMin = 0, jMax = 0;
  int kMin = 0, kMax = 0;

  // For the X-axis
  if( delta_pos[ 0 ] > 0.0 )
  {
    iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMin * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
    iMax = ::floor( 1 + ( event1[ 0 ] + alphaMax * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
  }
  else if( delta_pos[ 0 ] < 0.0 )
  {
    iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMax * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
    iMax = ::floor( 1 + ( event1[ 0 ] + alphaMin * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
  }
  if( delta_pos[ 0 ] == 0. )
  {
    iMin = 1, iMax = 0;
  }

  // For the Y-axis
  if( delta_pos[ 1 ] > 0. )
  {
    jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMin * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
    jMax = ::floor( 1 + ( event1[ 1 ] + alphaMax * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
  }
  else if( delta_pos[ 1 ] < 0. )
  {
    jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMax * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
    jMax = ::floor( 1 + ( event1[ 1 ] + alphaMin * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
  }
  else if( delta_pos[ 1 ] == 0. )
  {
    jMin = 1, jMax = 0;
  }

  // For the Z-axis
  if( delta_pos[ 2 ] > 0. )
  {
    kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMin * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
    kMax = ::floor( 1 + ( event1[ 2 ] + alphaMax * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
  }
  else if( delta_pos[ 2 ] < 0. )
  {
    kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMax * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
    kMax = ::floor( 1 + ( event1[ 2 ] + alphaMin * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
  }
  else if( delta_pos[ 2 ] == 0. )
  {
    kMin = 1, kMax = 0;
  }

  // Computing the last term n number of intersection
  int n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
    + ( kMax - kMin + 1 );

  // We create a buffer storing the merging data
  // We merge alphaMin, alphaMax, alphaX, alphaY and alphaZ
  FLTNB *alpha = new FLTNB[ mp_nbVox[ 0 ] + mp_nbVox[ 1 ] + mp_nbVox[ 2 ] + 3 ];
  FLTNB *tmpAlpha = new FLTNB[ mp_nbVox[ 0 ] + mp_nbVox[ 1 ] + mp_nbVox[ 2 ] + 3 ];
  FLTNB *alphaX = new FLTNB[ ( mp_nbVox[ 0 ] + 1 ) ];
  FLTNB *alphaY = new FLTNB[ ( mp_nbVox[ 1 ] + 1 ) ];
  FLTNB *alphaZ = new FLTNB[ ( mp_nbVox[ 2 ] + 1 ) ];

  INTNB iElement = iMax - iMin + 1;
  if( iElement > 0 )
  {
    FLTNB *idx = alphaX;
    if( delta_pos[ 0 ] > 0. )
    {
      for( int i = iMin; i <= iMax; ++i )
        *idx++ = ( ( (-mp_halfFOV[ 0 ]) + ( i - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    }
    else if( delta_pos[ 0 ] < 0. )
    {
      for( int i = iMax; i >= iMin; --i )
        *idx++ = ( ( (-mp_halfFOV[ 0 ]) + ( i - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    }
  }

  // For alphaY
  INTNB jElement = jMax - jMin + 1;
  if( jElement > 0 )
  {
    FLTNB *idx = alphaY;
    if( delta_pos[ 1 ] > 0. )
    {
      for( int j = jMin; j <= jMax; ++j )
        *idx++ = ( ( (-mp_halfFOV[ 1 ]) + ( j - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    }
    else if( delta_pos[ 1 ] < 0. )
    {
      for( int j = jMax; j >= jMin; --j )
        *idx++ = ( ( (-mp_halfFOV[ 1 ]) + ( j - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    }
  }

  // For alphaZ
  INTNB kElement = kMax - kMin + 1;
  if( kElement > 0 )
  {
    FLTNB *idx = alphaZ;
    if( delta_pos[ 2 ] > 0. )
    {
      for( int k = kMin; k <= kMax; ++k )
        *idx++ = ( ( (-mp_halfFOV[ 2 ]) + ( k - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    }
    else if( delta_pos[ 2 ] < 0. )
    {
      for( int k = kMax; k >= kMin; --k )
        *idx++ = ( ( (-mp_halfFOV[ 2 ]) + ( k - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    }
  }

  std::merge(
    alphaX, alphaX + iElement,
    tmpAlpha, tmpAlpha,
    alpha );

  ::memcpy( tmpAlpha, alpha, ( iElement ) * sizeof( FLTNB ) );

  std::merge(
    alphaY, alphaY + jElement,
    tmpAlpha, tmpAlpha + iElement,
    alpha );

  ::memcpy( tmpAlpha, alpha, ( iElement + jElement ) * sizeof( FLTNB ) );

  std::merge(
    alphaZ, alphaZ + kElement,
    tmpAlpha, tmpAlpha + iElement + jElement,
    alpha );

  // Computing some constants
  FLTNB const i1 = ( event1[ 0 ] - (-mp_halfFOV[ 0 ]) + mp_sizeVox[0] ) / mp_sizeVox[0];
  FLTNB const i2 = delta_pos[ 0 ] / mp_sizeVox[0];
  FLTNB const j1 = ( event1[ 1 ] - (-mp_halfFOV[ 1 ]) + mp_sizeVox[1] ) / mp_sizeVox[1];
  FLTNB const j2 = delta_pos[ 1 ] / mp_sizeVox[1];
  FLTNB const k1 = ( event1[ 2 ] - (-mp_halfFOV[ 2 ]) + mp_sizeVox[2] ) / mp_sizeVox[2];
  FLTNB const k2 = delta_pos[ 2 ] / mp_sizeVox[2];

  // Computing the index of the voxels
  FLTNB alphaMid = 0.0;
  INTNB i = 0, j = 0, k = 0; // indices of the voxel
  FLTNB *pAlpha = &alpha[ 1 ];
  FLTNB *pAlphaPrevious = &alpha[ 0 ];
  FLTNB coeff = 0.0;
  INTNB numVox = 0;

  // tof : erf of previous voxel edge - Gaussian center
  FLTNB previous_edge_erf =  erf( (length_LOR * (*pAlphaPrevious) - lor_tof_center)/ tof_sigma_sqrt2 );
  FLTNB next_edge_erf;

  // Loop over the number of crossed planes
  for( int nP = 0; nP < n - 1; ++nP, ++pAlpha, ++pAlphaPrevious )
  {
    alphaMid = ( *pAlpha + *pAlphaPrevious ) * 0.5;

    i = alphaMid * i2 + i1;
    if( i < 1 || i > mp_nbVox[ 0 ] ) continue;

    j = alphaMid * j2 + j1;
    if( j < 1 || j > mp_nbVox[ 1 ] ) continue;

    k = alphaMid * k2 + k1;
    if( k < 1 || k > mp_nbVox[ 2 ] ) continue;

    numVox = ( i - 1 ) + ( j - 1 ) * mp_nbVox[0] + ( k - 1 ) * mp_nbVox[0] * mp_nbVox[1];

    // tof : erf of next voxel edge - Gaussian center
    next_edge_erf = erf( (length_LOR * (*pAlpha) - lor_tof_center) / tof_sigma_sqrt2 );
    // integration of the Gaussian is done on the LOR portion matching the whole current voxel along X axis
    coeff = 0.5 * fabs(previous_edge_erf - next_edge_erf) ;
    // keeping record of the previous edge, so as to save 1 erf computation
    previous_edge_erf = next_edge_erf;

    // if this voxel is masked, skip
    if (!m_hasMask || mp_mask[numVox]) ap_ProjectionLine->AddVoxelInTOFBin(a_direction, 0, numVox, coeff);

  }

  delete[] alpha;
  delete[] tmpAlpha;
  delete[] alphaX;
  delete[] alphaY;
  delete[] alphaZ;

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    Cout("iProjectorClassicSiddon::Project with TOF measurement-> Exit function" << endl);
  }
  #endif

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorClassicSiddon::ProjectWithTOFBin(int a_direction, oProjectionLine* ap_ProjectionLine)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorClassicSiddon::ProjectWithTOFBin() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorClassicSiddon::Project with TOF bins -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Simpler way now, hopefully it works
  FLTNB* event1_position = ap_ProjectionLine->GetPosition1();
  FLTNB* event2_position = ap_ProjectionLine->GetPosition2();

  FLTNB event1[3] = { event1_position[0], event1_position[1], event1_position[2] };
  FLTNB event2[3] = { event2_position[0], event2_position[1], event2_position[2] };

  // **************************************
  // STEP 1: LOR length calculation
  // **************************************
  FLTNB length_LOR = ap_ProjectionLine->GetLength();
  FLTNB length_LOR_half = length_LOR * 0.5;

  // **************************************
  // STEP 2: Compute entrance voxel indexes
  // **************************************
  FLTNB alphaFirst[] = { 0.0, 0.0, 0.0 };
  FLTNB alphaLast[] = { 0.0, 0.0, 0.0 };

  FLTNB alphaMin = 0.0, alphaMax = 1.0;
  FLTNB delta_pos[] = {
    event2[ 0 ] - event1[ 0 ],
    event2[ 1 ] - event1[ 1 ],
    event2[ 2 ] - event1[ 2 ]
  };

  // Get TOF info
  // TOF temporal resolution (FWHM)
  HPFLTNB tof_resolutionT = ap_ProjectionLine->GetTOFResolution();
  // TOF bin temporal size
  HPFLTNB tof_bin_sizeT = ap_ProjectionLine->GetTOFBinSize();
  INTNB tof_nb_bins = ap_ProjectionLine->GetNbTOFBins();
  INTNB tof_half_nb_bins = tof_nb_bins/2;

  // TOF spatial standard deviation (simple TOF Gaussian, not convolved with the TOF bin width, maybe TODO)
  HPFLTNB tof_sigma = tof_resolutionT * SPEED_OF_LIGHT * 0.5 / TWO_SQRT_TWO_LN_2;
  //HPFLTNB tof_sigma_sqrt2 = sqrt(2.)*tof_sigma;
  HPFLTNB tof_half_span = tof_sigma * m_TOFnbSigmas;

  // normalized Gaussian (integral=1) : precomputation of the coefficient that multiplies the exponential
  HPFLTNB gaussian_norm_coef = INV_SQRT_2_PI * 1./tof_sigma;
  
  // spatial size of the TOF bin in mm
  HPFLTNB tof_bin_size = tof_bin_sizeT * SPEED_OF_LIGHT * 0.5;
  
  // minimum and maximum TOF bins
  INTNB tof_bin_last = tof_half_nb_bins;
  INTNB tof_bin_first = -tof_half_nb_bins;

  // distance between the first event1 and the center of the Gaussian distribution of the first TOF bin along the LOR
  //HPFLTNB tof_bin_center_first = length_LOR_half + tof_bin_first * tof_bin_size;
  // distance between the first event1 and the center of the Gaussian distribution of the current TOF bin along the LOR
  HPFLTNB lor_tof_center = 0.;
  // normalization coefficient used for making the sum of TOF bin coefficients for a voxel equal the nonTOF voxel coefficient
  HPFLTNB tof_norm_coef = 0.;

  // tof help variables
  //HPFLTNB previous_edge_erf, next_edge_erf;
  HPFLTNB tof_weight;

  // Computation of alphaMin et alphaMax values
  // entrance and exit point of the ray at voxel grid (edges), limited by the FOV,
  // absolute normalized distance from event1
  for( int i = 0; i < 3; ++i )
  {
    if( delta_pos[ i ] != 0.0 )
    {
      alphaFirst[i] = (-mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
      alphaLast[i]  = (mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
      alphaMin = (std::max)(alphaMin,(std::min)(alphaFirst[i],alphaLast[i]));
      alphaMax = (std::min)(alphaMax,(std::max)(alphaFirst[i],alphaLast[i]));
    }
  }

  // if alphaMax is less than or equal to alphaMin no intersection
  // and return an empty buffer
  if( alphaMax <= alphaMin ) return 0;

  // Min/max indices of voxels intersected by the LOR along each axis, 0-N indices (same order as absolute coordinates) match the FOV
  // (iMin,iMax), (jMin,jMax), (kMin,kMax)
  int iMin = 0, iMax = 0;
  int jMin = 0, jMax = 0;
  int kMin = 0, kMax = 0;

  // For the X-axis
  if( delta_pos[ 0 ] > 0.0 )
  {
    iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMin * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
    iMax = ::floor( 1 + ( event1[ 0 ] + alphaMax * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
  }
  else if( delta_pos[ 0 ] < 0.0 )
  {
    iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMax * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
    iMax = ::floor( 1 + ( event1[ 0 ] + alphaMin * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
  }
  if( delta_pos[ 0 ] == 0. )
  {
    iMin = 1, iMax = 0;
  }

  // For the Y-axis
  if( delta_pos[ 1 ] > 0. )
  {
    jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMin * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
    jMax = ::floor( 1 + ( event1[ 1 ] + alphaMax * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
  }
  else if( delta_pos[ 1 ] < 0. )
  {
    jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMax * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
    jMax = ::floor( 1 + ( event1[ 1 ] + alphaMin * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
  }
  else if( delta_pos[ 1 ] == 0. )
  {
    jMin = 1, jMax = 0;
  }

  // For the Z-axis
  if( delta_pos[ 2 ] > 0. )
  {
    kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMin * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
    kMax = ::floor( 1 + ( event1[ 2 ] + alphaMax * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
  }
  else if( delta_pos[ 2 ] < 0. )
  {
    kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMax * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
    kMax = ::floor( 1 + ( event1[ 2 ] + alphaMin * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
  }
  else if( delta_pos[ 2 ] == 0. )
  {
    kMin = 1, kMax = 0;
  }

  // Computing the last term n number of intersection
  int n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
    + ( kMax - kMin + 1 );

  // We create a buffer storing the merging data
  // We merge alphaMin, alphaMax, alphaX, alphaY and alphaZ
  FLTNB *alpha = new FLTNB[ mp_nbVox[ 0 ] + mp_nbVox[ 1 ] + mp_nbVox[ 2 ] + 3 ];
  FLTNB *tmpAlpha = new FLTNB[ mp_nbVox[ 0 ] + mp_nbVox[ 1 ] + mp_nbVox[ 2 ] + 3 ];
  FLTNB *alphaX = new FLTNB[ ( mp_nbVox[ 0 ] + 1 ) ];
  FLTNB *alphaY = new FLTNB[ ( mp_nbVox[ 1 ] + 1 ) ];
  FLTNB *alphaZ = new FLTNB[ ( mp_nbVox[ 2 ] + 1 ) ];

  // temporary storage for TOF bins Gaussian integrals over the current voxel
  // allocation after potential returns
  HPFLTNB* tof_weights_temp = new HPFLTNB[tof_nb_bins];
  for (INTNB tof_bin=0; tof_bin<tof_nb_bins; tof_bin++) tof_weights_temp[tof_bin] = 0.;
  
  INTNB iElement = iMax - iMin + 1;
  if( iElement > 0 )
  {
    FLTNB *idx = alphaX;
    if( delta_pos[ 0 ] > 0. )
    {
      for( int i = iMin; i <= iMax; ++i )
        *idx++ = ( ( (-mp_halfFOV[ 0 ]) + ( i - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    }
    else if( delta_pos[ 0 ] < 0. )
    {
      for( int i = iMax; i >= iMin; --i )
        *idx++ = ( ( (-mp_halfFOV[ 0 ]) + ( i - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    }
  }

  // For alphaY
  INTNB jElement = jMax - jMin + 1;
  if( jElement > 0 )
  {
    FLTNB *idx = alphaY;
    if( delta_pos[ 1 ] > 0. )
    {
      for( int j = jMin; j <= jMax; ++j )
        *idx++ = ( ( (-mp_halfFOV[ 1 ]) + ( j - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    }
    else if( delta_pos[ 1 ] < 0. )
    {
      for( int j = jMax; j >= jMin; --j )
        *idx++ = ( ( (-mp_halfFOV[ 1 ]) + ( j - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    }
  }

  // For alphaZ
  INTNB kElement = kMax - kMin + 1;
  if( kElement > 0 )
  {
    FLTNB *idx = alphaZ;
    if( delta_pos[ 2 ] > 0. )
    {
      for( int k = kMin; k <= kMax; ++k )
        *idx++ = ( ( (-mp_halfFOV[ 2 ]) + ( k - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    }
    else if( delta_pos[ 2 ] < 0. )
    {
      for( int k = kMax; k >= kMin; --k )
        *idx++ = ( ( (-mp_halfFOV[ 2 ]) + ( k - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    }
  }

  std::merge(
    alphaX, alphaX + iElement,
    tmpAlpha, tmpAlpha,
    alpha );

  ::memcpy( tmpAlpha, alpha, ( iElement ) * sizeof( FLTNB ) );

  std::merge(
    alphaY, alphaY + jElement,
    tmpAlpha, tmpAlpha + iElement,
    alpha );

  ::memcpy( tmpAlpha, alpha, ( iElement + jElement ) * sizeof( FLTNB ) );

  std::merge(
    alphaZ, alphaZ + kElement,
    tmpAlpha, tmpAlpha + iElement + jElement,
    alpha );

  // Computing some constants
  FLTNB const i1 = ( event1[ 0 ] - (-mp_halfFOV[ 0 ]) + mp_sizeVox[0] ) / mp_sizeVox[0];
  FLTNB const i2 = delta_pos[ 0 ] / mp_sizeVox[0];
  FLTNB const j1 = ( event1[ 1 ] - (-mp_halfFOV[ 1 ]) + mp_sizeVox[1] ) / mp_sizeVox[1];
  FLTNB const j2 = delta_pos[ 1 ] / mp_sizeVox[1];
  FLTNB const k1 = ( event1[ 2 ] - (-mp_halfFOV[ 2 ]) + mp_sizeVox[2] ) / mp_sizeVox[2];
  FLTNB const k2 = delta_pos[ 2 ] / mp_sizeVox[2];

  // Computing the index of the voxels
  FLTNB alphaMid = 0.0;
  INTNB i = 0, j = 0, k = 0; // indices of the voxel
  FLTNB *pAlpha = &alpha[ 1 ];
  FLTNB *pAlphaPrevious = &alpha[ 0 ];
  FLTNB coeff = 0.0;
  INTNB numVox = 0;
  INTNB tof_bin_first_for_voxel = 0, tof_bin_last_for_voxel = 0;
  
  // Loop over the number of crossed planes
  for( int nP = 0; nP < n - 1; ++nP, ++pAlpha, ++pAlphaPrevious )
  {
    alphaMid = ( *pAlpha + *pAlphaPrevious ) * 0.5;

    i = alphaMid * i2 + i1;
    if( i < 1 || i > mp_nbVox[ 0 ] ) continue;

    j = alphaMid * j2 + j1;
    if( j < 1 || j > mp_nbVox[ 1 ] ) continue;

    k = alphaMid * k2 + k1;
    if( k < 1 || k > mp_nbVox[ 2 ] ) continue;

    numVox = ( i - 1 ) + ( j - 1 ) * mp_nbVox[0] + ( k - 1 ) * mp_nbVox[0] * mp_nbVox[1];
    
    // if this voxel is masked, skip
    if (m_hasMask && !mp_mask[numVox]) continue;

    coeff = length_LOR * ( *pAlpha - *pAlphaPrevious );

    // Compute the first and the last relevant TOF bin for this voxel
    tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(trunc( ((*pAlphaPrevious) * length_LOR - tof_half_span - length_LOR_half ) / tof_bin_size )));
    tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(trunc( ((*pAlpha) * length_LOR + tof_half_span - length_LOR_half) / tof_bin_size )));
    lor_tof_center = length_LOR_half + tof_bin_first_for_voxel * tof_bin_size;
      
    // if min/max TOF bins for this voxel do not fall into the total available TOF bins range, skip
    if(tof_bin_first_for_voxel>tof_bin_last || tof_bin_last_for_voxel<tof_bin_first) continue;

    // shift tof bin indices from -/+ to 0:nbBins range
    tof_bin_first_for_voxel += tof_half_nb_bins;
    tof_bin_last_for_voxel += tof_half_nb_bins;

    // first compute the normalization for TOF bin coefficients for the current voxel (simple sum)
    tof_norm_coef = 0.;
    for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
    {
      // erf ( voxel edge projected onto the LOR - TOF bin center along the LOR )
      //previous_edge_erf = erf( (length_LOR * (*pAlphaPrevious) - lor_tof_center) / tof_sigma_sqrt2 );
      //next_edge_erf = erf( (length_LOR * (*pAlpha) - lor_tof_center) / tof_sigma_sqrt2 );
      // integration of the Gaussian is done on the LOR portion matching the whole current voxel along X axis
      //tof_weight = 0.5 * fabs(previous_edge_erf - next_edge_erf);
      // approximation of the integration of the Gaussian over the voxel
      // TOF coefficient = value of the Gaussian at voxel center (projected onto the LOR) * voxel size (projected onto the LOR)
      HPFLTNB temp = (alphaMid * length_LOR - lor_tof_center) / tof_sigma;
      tof_weight =  exp(- 0.5 * temp * temp ) * gaussian_norm_coef * coeff;
      // add the weight to the sum for normalization
      tof_norm_coef += tof_weight;
      // save the weight temporarily
      tof_weights_temp[tof_bin] = tof_weight;
      // update TOF center along the LOR for the next TOF bin
      lor_tof_center += tof_bin_size;
    }

    // compute the final TOF bin coefficients for the current voxel
    if (tof_norm_coef>0.)
    {
      // first normalize TOF bin coefficients
      for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
      {
        tof_weights_temp[tof_bin] /= tof_norm_coef;
      }
      // then write all TOF bins for the current voxel
      ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, numVox, coeff, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
    }
  }

  delete[] alpha;
  delete[] tmpAlpha;
  delete[] alphaX;
  delete[] alphaY;
  delete[] alphaZ;
  delete[] tof_weights_temp;

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    Cout("iProjectorClassicSiddon::Project with TOF bins -> Exit function" << endl);
  }
  #endif

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
