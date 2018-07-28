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
  \brief    Implementation of class iProjectorDistanceDriven
*/

#include "iProjectorDistanceDriven.hh"
#include "sOutputManager.hh"

#include <cmath>
#include <algorithm>
#ifdef _WIN32
// Avoid compilation errors due to mix up between std::min()/max() and
// min max macros
#undef min
#undef max
#endif

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorDistanceDriven::iProjectorDistanceDriven() : vProjector()
{
  // This projector is not compatible with SPECT attenuation correction because
  // the voxels contributing to the line are not strictly ordered with respect to
  // their distance to point2 (due to interpolations at each plane crossed)
  m_compatibleWithSPECTAttenuationCorrection = false;
  // This projector is not compatible with compression as it works only with the
  // detection element indices
  m_compatibleWithCompression = false;
  // Default pointers and parameters
  m_toleranceX = 0.;
  m_toleranceY = 0.;
  m_toleranceZ = 0.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorDistanceDriven::~iProjectorDistanceDriven()
{
  ;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorDistanceDriven::ReadConfigurationFile(const string& a_configurationFile)
{
  // No options for distance driven
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorDistanceDriven::ReadOptionsList(const string& a_optionsList)
{
  // No options for distance driven
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iProjectorDistanceDriven::ShowHelpSpecific()
{
  cout << "This projector is a line projector based on computations of overlap between a pair of detection elements and voxels." << endl;
  cout << "It is implemented from the following published paper:" << endl;
  cout << "B. De Man and S. Basu, \"Distance-driven projection and backprojection in three dimensions\", Phys. Med. Biol., vol. 49, pp. 2463-75, 2004." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorDistanceDriven::CheckSpecificParameters()
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

int iProjectorDistanceDriven::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=2) Cout("iProjectorDistanceDriven::InitializeSpecific() -> Use Distance Driven projector" << endl);

  // Set the tolerance with respect to voxel sizes in each dimensions
  HPFLTNB tolerance_factor = 1.e-4;
  m_toleranceX = ((HPFLTNB)(mp_sizeVox[0])) * tolerance_factor;
  m_toleranceY = ((HPFLTNB)(mp_sizeVox[1])) * tolerance_factor;
  m_toleranceZ = ((HPFLTNB)(mp_sizeVox[2])) * tolerance_factor;

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB iProjectorDistanceDriven::EstimateMaxNumberOfVoxelsPerLine()
{
  // For this projector, we estimated that the number of voxels in a slice will be always higher than the number of contributing voxels for any line of response
  return mp_ImageDimensionsAndQuantification->GetNbVoxXY();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorDistanceDriven::ProjectWithoutTOF(int a_direction, oProjectionLine* ap_ProjectionLine )
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorDistanceDriven::ProjectWithoutTOF() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorDistanceDriven::Project without TOF -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Computing the coordinates of the points on the bound of the source/crystal
  // 1. The first point is on transaxial plan and the last point is on axial
  // plan.
  // Source spot size (CT) / Crystal 1 (PET):
  //
  //                  ------- P13 -------
  //                  |                  |
  //                  |                  |
  //                 P12      Center     P11
  //                  |         0        |
  //                  |                  |
  //                  ------- P14 -------
  //  Y
  //  <---------. X
  //            |
  //            |
  //            |
  //           \/ Z
  // --------------

  // Computing the position of the point2.
  // In Distance Driven the point p21 relies the point P11, p22 -> p12,
  // p23 -> p13 and p24 -> p14
  // Pixel (CT) / Crystal 2 (PET):
  //
  //                  ------- P23 -------
  //                  |                  |
  //                  |                  |
  //                 P22      Center     P21
  //                  |         0        |
  //                  |                  |
  //                  ------- P24 -------
  //  Y
  //  <---------. X
  //            |
  //            |
  //            |
  //           \/ Z
  // --------------

  // buffers storing point 1 and 2
  FLTNB p1_x[ 5 ];
  FLTNB p1_y[ 5 ];
  FLTNB p1_z[ 5 ];
  FLTNB p2_x[ 5 ];
  FLTNB p2_y[ 5 ];
  FLTNB p2_z[ 5 ];

  // Get the positions of the centers of the edges of crystal elements or source position
  if (mp_Scanner->GetEdgesCenterPositions(
        ap_ProjectionLine->GetIndex1(), // Index 1
        ap_ProjectionLine->GetIndex2(), // Index 2
        ap_ProjectionLine->GetPosition1(), // Line position for point 1
        ap_ProjectionLine->GetPosition2(), // Line position for point 2
        &p1_x[1], &p1_y[1], &p1_z[1], // Edges for point 1
        &p2_x[1], &p2_y[1], &p2_z[1] // Edges for point 2
    ))
  {
    Cerr("***** iProjectorDistanceDriven::ProjectWithoutTOF() -> A problem occured while getting the edges' center positions !" << endl);
    return 1;
  }

  // Point 1 focal (CT) or crystal (PET) position
  FLTNB* p1 = ap_ProjectionLine->GetPosition1();

  // Point 2 pixel (CT) or crystal (PET) position
  FLTNB* p2 = ap_ProjectionLine->GetPosition2();

  // Center position p1
  p1_x[ 0 ] = p1[ 0 ]; p1_y[ 0 ] = p1[ 1 ]; p1_z[ 0 ] = p1[ 2 ];

  // Center position p2
  p2_x[ 0 ] = p2[ 0 ]; p2_y[ 0 ] = p2[ 1 ]; p2_z[ 0 ] = p2[ 2 ];

  // storing the ray relying p1 and p2
  HPFLTNB r[ 5 ][ 3 ];

  // buffer storing the intersection
  HPFLTNB ri[ 4 ][ 3 ];

  // Take the position of ray relying p1 and p2 for all the 5 points
  for( int p = 0; p < 5; ++p )
  {
    r[ p ][ 0 ] = p2_x[ p ] - p1_x[ p ];
    r[ p ][ 1 ] = p2_y[ p ] - p1_y[ p ];
    r[ p ][ 2 ] = p2_z[ p ] - p1_z[ p ];
  }

  // Computing the square of distance for the center
  HPFLTNB const r2[ 3 ] = {
    r[ 0 ][ 0 ] * r[ 0 ][ 0 ],
    r[ 0 ][ 1 ] * r[ 0 ][ 1 ],
    r[ 0 ][ 2 ] * r[ 0 ][ 2 ]
  };

  // Find the first and last intersecting plane using the parametric
  // values alpha_min and alpha_max
  HPFLTNB alpha_min = 0.0, alpha_max = 1.0;

  // Buffer storing the alpha on each axis
  // 2 different buffers to get the min and the max
  HPFLTNB alpha_x_min[ 4 ], alpha_x_max[ 4 ];
  HPFLTNB alpha_y_min[ 4 ], alpha_y_max[ 4 ];
  HPFLTNB alpha_z_min[ 4 ], alpha_z_max[ 4 ];

  // Position of the projected points in plane
  // If the main direction is X:
  // pos[0] is the position of point1 in Y
  // pos[1] is the position of point2 in Y
  // pos[2] is the position of point3 in Z
  // pos[3] is the position of point4 in Z
  HPFLTNB pos[ 4 ];

  // Width of intersection
  // w1 normalizes the overlap 1
  // w2 normalizes the overlap 2
  HPFLTNB w1, w2;

  // Buffer storing the min. and max. indices bounding the overlap
  int index[ 4 ];

  // Buffer storing the distances
  // 50 HPFLTNB is largely enough!!!
  HPFLTNB distance_x[ 50 ];
  HPFLTNB distance_y[ 50 ];
  HPFLTNB distance_z[ 50 ];

  // Computing the angle of line vs. horizontal
  float const angle = std::acos( ( r[ 0 ][ 0 ] ) / ( std::sqrt( r2[ 0 ] + r2[ 1 ] ) ) ) * 180.0f / M_PI;

  // Condition on the largest component of r, taking the center
  if( ( angle >= 0.0 && angle <= 45.0 ) || ( angle >= 135.0 && angle <= 180.0 ) )
  {
    // Computing the parametric values alpha_min and alpha_max
    // Because the main direction is X, we use only the center
    alpha_x_min[ 0 ] = ( (-mp_halfFOV[ 0 ]) - p1_x[ 0 ] ) / r[ 0 ][ 0 ];
    alpha_x_min[ 1 ] = ( mp_halfFOV[ 0 ] - p1_x[ 0 ] ) / r[ 0 ][ 0 ];
    alpha_x_max[ 0 ] = alpha_x_min[ 0 ];
    alpha_x_max[ 1 ] = alpha_x_min[ 1 ];
    alpha_min = std::max( alpha_min, std::min( alpha_x_min[ 0 ], alpha_x_min[ 1 ] ) );
    alpha_max = std::min( alpha_max, std::max( alpha_x_max[ 0 ], alpha_x_max[ 1 ] ) );

    // Compute the parametric value for each line
    // For Y, we use only the point 1 and 2
    if( r[ 1 ][ 1 ] != 0.0 ) // point1
    {
      alpha_y_min[ 0 ] = ( (-mp_halfFOV[ 1 ]) - p1_y[ 1 ] ) / r[ 1 ][ 1 ];
      alpha_y_min[ 1 ] = ( mp_halfFOV[ 1 ] - p1_y[ 1 ] ) / r[ 1 ][ 1 ];
      alpha_y_max[ 0 ] = alpha_y_min[ 0 ];
      alpha_y_max[ 1 ] = alpha_y_min[ 1 ];
    }
    else
    {
      alpha_y_min[ 0 ] = 0.0;
      alpha_y_min[ 1 ] = 0.0;
      alpha_y_max[ 0 ] = 1.0;
      alpha_y_max[ 1 ] = 1.0;
    }

    if( r[ 2 ][ 1 ] != 0.0 ) // point2
    {
      alpha_y_min[ 2 ] = ( (-mp_halfFOV[ 1 ]) - p1_y[ 2 ] ) / r[ 2 ][ 1 ];
      alpha_y_min[ 3 ] = ( mp_halfFOV[ 1 ] - p1_y[ 2 ] ) / r[ 2 ][ 1 ];
      alpha_y_max[ 2 ] = alpha_y_min[ 2 ];
      alpha_y_max[ 3 ] = alpha_y_min[ 3 ];
    }
    else
    {
      alpha_y_min[ 2 ] = 0.0;
      alpha_y_min[ 3 ] = 0.0;
      alpha_y_max[ 2 ] = 1.0;
      alpha_y_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_y_min
    alpha_min = std::max( alpha_min,
      alpha_y_min[ std::minmax_element( alpha_y_min, alpha_y_min + 4 ).first - alpha_y_min ] );
    // Getting the maximum of alpha_y_max
    alpha_max = std::min( alpha_max,
      alpha_y_max[ std::minmax_element( alpha_y_max, alpha_y_max + 4 ).second - alpha_y_max ] );

    // For Z, we use only the point 3 and 4
    if( r[ 3 ][ 2 ] != 0.0 ) // point3
    {
      alpha_z_min[ 0 ] = ( (-mp_halfFOV[ 2 ]) - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_min[ 1 ] = ( mp_halfFOV[ 2 ] - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_max[ 0 ] = alpha_z_min[ 0 ];
      alpha_z_max[ 1 ] = alpha_z_min[ 1 ];
    }
    else
    {
      alpha_z_min[ 0 ] = 0.0;
      alpha_z_min[ 1 ] = 0.0;
      alpha_z_max[ 0 ] = 1.0;
      alpha_z_max[ 1 ] = 1.0;
    }

    if( r[ 4 ][ 2 ] != 0.0 ) // point4
    {
      alpha_z_min[ 2 ] = ( (-mp_halfFOV[ 2 ]) - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_min[ 3 ] = ( mp_halfFOV[ 2 ] - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_max[ 2 ] = alpha_z_min[ 2 ];
      alpha_z_max[ 3 ] = alpha_z_min[ 3 ];
    }
    else
    {
      alpha_z_min[ 2 ] = 0.0;
      alpha_z_min[ 3 ] = 0.0;
      alpha_z_max[ 2 ] = 1.0;
      alpha_z_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_z_min
    alpha_min = std::max( alpha_min,
      alpha_z_min[ std::minmax_element( alpha_z_min, alpha_z_min + 4 ).first
      - alpha_z_min ] );
    // Getting the maximum of alpha_z_max
    alpha_max = std::min( alpha_max,
      alpha_z_max[ std::minmax_element( alpha_z_max, alpha_z_max + 4 ).second
      - alpha_z_max ] );

    if( alpha_max <= alpha_min ) return 0;

    // Computing the first and the last plane
    int16_t i_min = 0, i_max = 0;
    if( r[ 0 ][ 0 ] > 0.0 )
    {
      i_min = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alpha_min * r[ 0 ][ 0 ] - p1_x[ 0 ] ) / mp_sizeVox[ 0 ] );
      i_max = ::floor( m_toleranceX + ( p1_x[ 0 ] + alpha_max * r[ 0 ][ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }
    else if( r[ 0 ][ 0 ] < 0.0 )
    {
      i_min = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alpha_max * r[ 0 ][ 0 ] - p1_x[ 0 ] ) / mp_sizeVox[ 0 ] );
      i_max = ::floor( m_toleranceX + ( p1_x[ 0 ] + alpha_min * r[ 0 ][ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }

    // Computing weight of normalization
    // Using the center
    HPFLTNB const factor( ::fabs( r[ 0 ][ 0 ] ) / ::sqrt( r2[ 0 ] + r2[ 1 ] + r2[ 2 ] ) );
    HPFLTNB weight( mp_sizeVox[ 0 ] / factor );

    // Computing the increment for each 4 lines to project
    // (the center is not projected)
    for( uint8_t p = 0; p < 4; ++p )
    {
      ri[ p ][ 0 ] = 1.0; //r[ p + 1 ][ 0 ] / r[ p + 1 ][ 0 ]
      ri[ p ][ 1 ] = r[ p + 1 ][ 1 ] / r[ p + 1 ][ 0 ];
      ri[ p ][ 2 ] = r[ p + 1 ][ 2 ] / r[ p + 1 ][ 0 ];
    }

    // Increment orientation
    int8_t incr_orient_trans = 0;
    int8_t idx_orient_trans = 0;

    //int8_t const incr_orient_trans = p2_y[ 1 ] < p2_y[ 2 ] ? 1 : -1;
    int8_t const incr_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : -1;

    // Index orientation
    //int8_t const idx_orient_trans = p2_y[ 1 ] < p2_y[ 2 ] ? 1 : 0;
    int8_t const idx_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : 0;

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = (-mp_halfFOV[ 0 ]) + ( mp_sizeVox[ 0 ] * 0.5 ) - p1_x[ 0 ];

    // Loop over the crossed X planes
    for( int16_t i = i_min; i < i_max; ++i )
    {
      // Computing the coordinates of crossed plane
      HPFLTNB const step = offset + i * mp_sizeVox[ 0 ];
      // in Y (point 1 and 2)
      pos[ 0 ] = p1_y[ 1 ] + step * ri[ 0 ][ 1 ];
      pos[ 1 ] = p1_y[ 2 ] + step * ri[ 1 ][ 1 ];
      // in Z (point 3 and 4)
      pos[ 2 ] = p1_z[ 3 ] + step * ri[ 2 ][ 2 ];
      pos[ 3 ] = p1_z[ 4 ] + step * ri[ 3 ][ 2 ];

      // Computing the factor w1 and w2 normalizing the overlaps
      w1 = ::fabs( pos[ 0 ] - pos[ 1 ] );
      w2 = ::fabs( pos[ 2 ] - pos[ 3 ] );
      HPFLTNB final_weight = 0.0;
      if( ( w1 * w2 ) != 0 ) final_weight = weight / ( w1 * w2 );

      // Computing the index min and max on each axis
      // In Y
      index[ 0 ] = ::floor( m_toleranceY + ( pos[ 0 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
      index[ 1 ] = ::floor( m_toleranceY + ( pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
      // In Z
      index[ 2 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );
      index[ 3 ] = ::floor( m_toleranceZ + ( pos[ 3 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );

      // apply mask if required
      int test_index1 = i + index[ 0 ] * mp_nbVox[ 0 ] + index[ 2 ] * m_nbVoxXY;
      int test_index2 = i + index[ 1 ] * mp_nbVox[ 0 ] + index[ 3 ] * m_nbVoxXY;
      bool test_index = (index[ 0 ] >= 0 && index[ 2 ] >= 0 && index[ 1 ] >= 0 && index[ 3 ] >= 0
                      && index[ 0 ] < mp_nbVox[ 1 ] && index[ 1 ] < mp_nbVox[ 1 ] && index[ 2 ] < mp_nbVox[ 2 ] && index[ 3 ] < mp_nbVox[ 2 ]);
      if (m_hasMask && test_index && !mp_mask[test_index1] && !mp_mask[test_index2]) continue;

      if( index[ 0 ] < index[ 1 ] )
      {
        incr_orient_trans = 1;
        idx_orient_trans = 1;
      }
      else if( index[ 0 ] > index[ 1 ] )
      {
        incr_orient_trans = -1;
        idx_orient_trans = 0;
      }

      // Getting the number of distance in Y
      int16_t const n_distance_y = ::abs( index[ 0 ] - index[ 1 ] ) + 2;
      // Computing the distances in Y
      distance_y[ 0 ] = pos[ 0 ];
      distance_y[ n_distance_y - 1 ] = pos[ 1 ];
      // Computing the rest if necessary
      if( n_distance_y > 2 )
      {
        for( int16_t d = 0; d < n_distance_y - 2; ++d )
          distance_y[ d + 1 ] = (-mp_halfFOV[ 1 ]) + ( index[ 0 ] + idx_orient_trans + ( incr_orient_trans * d ) ) * mp_sizeVox[ 1 ];
      }

      // Computing the final distance in Y
      for( int16_t d = 0; d < n_distance_y - 1; ++d )
        distance_y[ d ] = ::fabs( distance_y[ d + 1 ] - distance_y[ d ] );

      // Getting the number of distance in Z
      int16_t const n_distance_z = ::abs( index[ 2 ] - index[ 3 ] ) + 2;
      // Storing the positions in Y
      distance_z[ 0 ] = pos[ 2 ];
      distance_z[ n_distance_z - 1 ] = pos[ 3 ];
      // Storing the rest if necessary
      if( n_distance_z > 2 )
      {
        for( int16_t d = 0; d < n_distance_z - 2; ++d )
          distance_z[ d + 1 ] = (-mp_halfFOV[ 2 ]) + ( index[ 2 ] + idx_orient_axial + ( incr_orient_axial * d ) ) * mp_sizeVox[ 2 ];
      }

      // Computing the final distance in Z
      for( int16_t d = 0; d < n_distance_z - 1; ++d )
        distance_z[ d ] = ::fabs( distance_z[ d + 1 ] - distance_z[ d ] );

      // Loop over the overlap and store the elements
      int16_t index_y, index_z;
      for( int16_t jj = 0; jj < n_distance_z - 1; ++jj )
      {
        index_z = index[ 2 ] + jj * incr_orient_axial;
        if( index_z < 0 || index_z > mp_nbVox[ 2 ] - 1 ) continue;

        for( int16_t ii = 0; ii < n_distance_y - 1; ++ii )
        {
          index_y = index[ 0 ] + ii * incr_orient_trans;
          if( index_y < 0 || index_y > mp_nbVox[ 1 ] - 1 ) continue;

          ap_ProjectionLine->AddVoxel(a_direction, i + index_y * mp_nbVox[ 0 ] + index_z * m_nbVoxXY, final_weight * distance_z[ jj ] * distance_y[ ii ]);
        }
      }
    }
  }
  else
  {
    // Compute the parametric value for each line
    // For X, we use only the point 1 and 2
    if( r[ 1 ][ 0 ] != 0.0 ) // point1
    {
      alpha_x_min[ 0 ] = ( (-mp_halfFOV[ 0 ]) - p1_x[ 1 ] ) / r[ 1 ][ 0 ];
      alpha_x_min[ 1 ] = ( mp_halfFOV[ 0 ] - p1_x[ 1 ] ) / r[ 1 ][ 0 ];
      alpha_x_max[ 0 ] = alpha_x_min[ 0 ];
      alpha_x_max[ 1 ] = alpha_x_min[ 1 ];
    }
    else
    {
      alpha_x_min[ 0 ] = 0.0;
      alpha_x_min[ 1 ] = 0.0;
      alpha_x_max[ 0 ] = 1.0;
      alpha_x_max[ 1 ] = 1.0;
    }

    if( r[ 2 ][ 0 ] != 0.0 ) // point2
    {
      alpha_x_min[ 2 ] = ( (-mp_halfFOV[ 0 ]) - p1_x[ 2 ] ) / r[ 2 ][ 0 ];
      alpha_x_min[ 3 ] = ( mp_halfFOV[ 0 ] - p1_x[ 2 ] ) / r[ 2 ][ 0 ];
      alpha_x_max[ 2 ] = alpha_x_min[ 2 ];
      alpha_x_max[ 3 ] = alpha_x_min[ 3 ];
    }
    else
    {
      alpha_x_min[ 2 ] = 0.0;
      alpha_x_min[ 3 ] = 0.0;
      alpha_x_max[ 2 ] = 1.0;
      alpha_x_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_x_min
    alpha_min = std::max( alpha_min,
      alpha_x_min[ std::minmax_element( alpha_x_min, alpha_x_min + 4 ).first - alpha_x_min ] );
    // Getting the maximum of alpha_y_max
    alpha_max = std::min( alpha_max,
      alpha_x_max[ std::minmax_element( alpha_x_max, alpha_x_max + 4 ).second - alpha_x_max ] );

    // Computing the parametric values alpha_min and alpha_max
    // Because the main direction is Y, we use only the center
    alpha_y_min[ 0 ] = ( (-mp_halfFOV[ 1 ]) - p1_y[ 0 ] ) / r[ 0 ][ 1 ];
    alpha_y_min[ 1 ] = ( mp_halfFOV[ 1 ] - p1_y[ 0 ] ) / r[ 0 ][ 1 ];
    alpha_y_max[ 0 ] = alpha_y_min[ 0 ];
    alpha_y_max[ 1 ] = alpha_y_min[ 1 ];
    alpha_min = std::max( alpha_min,
      std::min( alpha_y_min[ 0 ], alpha_y_min[ 1 ] ) );
    alpha_max = std::min( alpha_max,
      std::max( alpha_y_max[ 0 ], alpha_y_max[ 1 ] ) );

    // For Z, we use only the point 3 and 4
    if( r[ 3 ][ 2 ] != 0.0 ) // point3
    {
      alpha_z_min[ 0 ] = ( (-mp_halfFOV[ 2 ]) - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_min[ 1 ] = ( mp_halfFOV[ 2 ] - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_max[ 0 ] = alpha_z_min[ 0 ];
      alpha_z_max[ 1 ] = alpha_z_min[ 1 ];
    }
    else
    {
      alpha_z_min[ 0 ] = 0.0;
      alpha_z_min[ 1 ] = 0.0;
      alpha_z_max[ 0 ] = 1.0;
      alpha_z_max[ 1 ] = 1.0;
    }

    if( r[ 4 ][ 2 ] != 0.0 ) // point4
    {
      alpha_z_min[ 2 ] = ( (-mp_halfFOV[ 2 ]) - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_min[ 3 ] = ( mp_halfFOV[ 2 ] - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_max[ 2 ] = alpha_z_min[ 2 ];
      alpha_z_max[ 3 ] = alpha_z_min[ 3 ];
    }
    else
    {
      alpha_z_min[ 2 ] = 0.0;
      alpha_z_min[ 3 ] = 0.0;
      alpha_z_max[ 2 ] = 1.0;
      alpha_z_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_z_min
    alpha_min = std::max( alpha_min, alpha_z_min[ std::minmax_element( alpha_z_min, alpha_z_min + 4 ).first - alpha_z_min ] );
    // Getting the maximum of alpha_z_max
    alpha_max = std::min( alpha_max, alpha_z_max[ std::minmax_element( alpha_z_max, alpha_z_max + 4 ).second - alpha_z_max ] );

    if( alpha_max <= alpha_min ) return 0;

    // Computing the first and the last plane
    int16_t j_min = 0, j_max = 0;
    if( r[ 0 ][ 1 ] > 0.0 )
    {
      j_min = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alpha_min * r[ 0 ][ 1 ] - p1_y[ 0 ] ) / mp_sizeVox[ 1 ] );
      j_max = ::floor( m_toleranceY + ( p1_y[ 0 ] + alpha_max * r[ 0 ][ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }
    else if( r[ 0 ][ 1 ] < 0.0 )
    {
      j_min = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alpha_max * r[ 0 ][ 1 ] - p1_y[ 0 ] ) / mp_sizeVox[ 1 ] );
      j_max = ::floor( m_toleranceY + ( p1_y[ 0 ] + alpha_min * r[ 0 ][ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }

    // Computing weight of normalization
    // Using the center
    HPFLTNB const factor( ::fabs( r[ 0 ][ 1 ] ) / ::sqrt( r2[ 0 ] + r2[ 1 ] + r2[ 2 ] ) );
    HPFLTNB weight( mp_sizeVox[ 1 ] / factor );

    // Computing the increment for each 4 lines to project
    // (the center is not projected)
    for( uint8_t p = 0; p < 4; ++p )
    {
      ri[ p ][ 0 ] = r[ p + 1 ][ 0 ] / r[ p + 1 ][ 1 ];
      ri[ p ][ 1 ] = 1.0;
      ri[ p ][ 2 ] = r[ p + 1 ][ 2 ] / r[ p + 1 ][ 1 ];
    }

    // Increment orientation and Index orientation
    int8_t incr_orient_trans = 0;
    int8_t idx_orient_trans = 0;

    //int8_t const incr_orient_trans = p2_x[ 1 ] < p2_x[ 2 ] ? 1 : -1;
    int8_t const incr_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : -1;
    //int8_t const idx_orient_trans = p2_x[ 1 ] < p2_x[ 2 ] ? 1 : 0;
    int8_t const idx_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : 0;

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = (-mp_halfFOV[ 1 ]) + ( mp_sizeVox[ 1 ] * 0.5 ) - p1_y[ 0 ];

    // Loop over the crossed Y planes
    for( int16_t j = j_min; j < j_max; ++j )
    {
      // Computing the coordinates of crossed plane
      HPFLTNB const step = offset + j * mp_sizeVox[ 1 ];
      // in X (point 1 and 2)
      pos[ 0 ] = p1_x[ 1 ] + step * ri[ 0 ][ 0 ];
      pos[ 1 ] = p1_x[ 2 ] + step * ri[ 1 ][ 0 ];
      // in Z (point 3 and 4)
      pos[ 2 ] = p1_z[ 3 ] + step * ri[ 2 ][ 2 ];
      pos[ 3 ] = p1_z[ 4 ] + step * ri[ 3 ][ 2 ];

      // Computing the factor w1 and w2 normalizing the overlaps
      w1 = ::fabs( pos[ 0 ] - pos[ 1 ] );
      w2 = ::fabs( pos[ 2 ] - pos[ 3 ] );
      HPFLTNB final_weight = 0.0;
      if( ( w1 * w2 ) != 0 ) final_weight = weight / ( w1 * w2 );

      // Computing the index min and max on each axis
      // In X
      index[ 0 ] = ::floor( m_toleranceX + ( pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
      index[ 1 ] = ::floor( m_toleranceX + ( pos[ 1 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
      // In Z
      index[ 2 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );
      index[ 3 ] = ::floor( m_toleranceZ + ( pos[ 3 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );

      // if voxel masked, skip
      int test_index1 = index[ 0 ] + j* mp_nbVox[ 0 ] + index[ 2 ] * m_nbVoxXY;
      int test_index2 = index[ 1 ] + j* mp_nbVox[ 0 ] + index[ 3 ] * m_nbVoxXY;
      bool test_index = (index[ 0 ] >= 0 && index[ 2 ] >= 0 && index[ 1 ] >= 0 && index[ 3 ] >= 0
                      && index[ 0 ] < mp_nbVox[ 0 ] && index[ 1 ] < mp_nbVox[ 0 ] && index[ 2 ] < mp_nbVox[ 2 ] && index[ 3 ] < mp_nbVox[ 2 ]);
      if (m_hasMask && test_index && !mp_mask[test_index1] && !mp_mask[test_index2]) continue;

      if( index[ 0 ] < index[ 1 ] )
      {
        incr_orient_trans = 1;
        idx_orient_trans = 1;
      }
      else if( index[ 0 ] > index[ 1 ] )
      {
        incr_orient_trans = -1;
        idx_orient_trans = 0;
      }

      // Getting the number of distance in X
      int16_t const n_distance_x = ::abs( index[ 0 ] - index[ 1 ] ) + 2;
      // Computing the distances in X
      distance_x[ 0 ] = pos[ 0 ];
      distance_x[ n_distance_x - 1 ] = pos[ 1 ];
      // Computing the rest if necessary
      if( n_distance_x > 2 )
      {
        for( int16_t d = 0; d < n_distance_x - 2; ++d )
          distance_x[ d + 1 ] = (-mp_halfFOV[ 0 ]) + ( index[ 0 ] + idx_orient_trans + ( incr_orient_trans * d ) ) * mp_sizeVox[ 0 ];
      }

      // Computing the final distance in X
      for( int16_t d = 0; d < n_distance_x - 1; ++d )
        distance_x[ d ] = ::fabs( distance_x[ d + 1 ] - distance_x[ d ] );

      // Getting the number of distance in Z
      int16_t const n_distance_z = ::abs( index[ 2 ] - index[ 3 ] ) + 2;
      // Storing the positions in Y
      distance_z[ 0 ] = pos[ 2 ];
      distance_z[ n_distance_z - 1 ] = pos[ 3 ];
      // Storing the rest if necessary
      if( n_distance_z > 2 )
      {
        for( int16_t d = 0; d < n_distance_z - 2; ++d )
          distance_z[ d + 1 ] = (-mp_halfFOV[ 2 ]) + ( index[ 2 ] + idx_orient_axial + ( incr_orient_axial * d ) ) * mp_sizeVox[ 2 ];
      }

      // Computing the final distance in Z
      for( int16_t d = 0; d < n_distance_z - 1; ++d )
        distance_z[ d ] = ::fabs( distance_z[ d + 1 ] - distance_z[ d ] );

      // Loop over the overlap and store the elements
      int16_t index_x, index_z;
      for( int16_t jj = 0; jj < n_distance_z - 1; ++jj )
      {
        index_z = index[ 2 ] + jj * incr_orient_axial;
        if( index_z < 0 || index_z > mp_nbVox[ 2 ] - 1 ) continue;

        for( int16_t ii = 0; ii < n_distance_x - 1; ++ii )
        {
          index_x = index[ 0 ] + ii * incr_orient_trans;
          if( index_x < 0 || index_x > mp_nbVox[ 0 ] - 1 ) continue;

          ap_ProjectionLine->AddVoxel(a_direction, index_x + j * mp_nbVox[ 0 ] + index_z * m_nbVoxXY, final_weight * distance_z[ jj ] * distance_x[ ii ]);
        }
      }
    }
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorDistanceDriven::ProjectWithTOFPos(int a_direction, oProjectionLine* ap_ProjectionLine)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorDistanceDriven::ProjectWithTOFPos() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorDistanceDriven::Project with TOF measurement -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Computing the coordinates of the points on the bound of the source/crystal
  // 1. The first point is on transaxial plan and the last point is on axial
  // plan.
  // Source spot size (CT) / Crystal 1 (PET):
  //
  //                  ------- P13 -------
  //                  |                  |
  //                  |                  |
  //                 P12      Center     P11
  //                  |         0        |
  //                  |                  |
  //                  ------- P14 -------
  //  Y
  //  <---------. X
  //            |
  //            |
  //            |
  //           \/ Z
  // --------------

  // Computing the position of the point2.
  // In Distance Driven the point p21 relies the point P11, p22 -> p12,
  // p23 -> p13 and p24 -> p14
  // Pixel (CT) / Crystal 2 (PET):
  //
  //                  ------- P23 -------
  //                  |                  |
  //                  |                  |
  //                 P22      Center     P21
  //                  |         0        |
  //                  |                  |
  //                  ------- P24 -------
  //  Y
  //  <---------. X
  //            |
  //            |
  //            |
  //           \/ Z
  // --------------

  // buffers storing point 1 and 2
  FLTNB p1_x[ 5 ];
  FLTNB p1_y[ 5 ];
  FLTNB p1_z[ 5 ];
  FLTNB p2_x[ 5 ];
  FLTNB p2_y[ 5 ];
  FLTNB p2_z[ 5 ];

  if (mp_Scanner->GetEdgesCenterPositions(
        ap_ProjectionLine->GetIndex1(), // Index 1
        ap_ProjectionLine->GetIndex2(), // Index 2
        ap_ProjectionLine->GetPosition1(), // Line position for point 1
        ap_ProjectionLine->GetPosition2(), // Line position for point 2
        &p1_x[1], &p1_y[1], &p1_z[1], // Edges for point 1
        &p2_x[1], &p2_y[1], &p2_z[1] // Edges for point 2
    ))
  {
    Cerr("***** iProjectorDistanceDriven::ProjectWithoutTOF() -> A problem occured while getting the edges' center positions !" << endl);
    return 1;
  }

  // Point 1 focal (CT) or crystal (PET) position
  FLTNB* p1 = ap_ProjectionLine->GetPosition1();

  // Point 2 pixel (CT) or crystal (PET) position
  FLTNB* p2 = ap_ProjectionLine->GetPosition2();

  // Center position p1
  p1_x[ 0 ] = p1[ 0 ]; p1_y[ 0 ] = p1[ 1 ]; p1_z[ 0 ] = p1[ 2 ];

  // Center position p2
  p2_x[ 0 ] = p2[ 0 ]; p2_y[ 0 ] = p2[ 1 ]; p2_z[ 0 ] = p2[ 2 ];

  // storing the ray relying p1 and p2
  HPFLTNB r[ 5 ][ 3 ];

  // buffer storing the intersection
  HPFLTNB ri[ 4 ][ 3 ];

  // Take the position of ray relying p1 and p2 for all the 5 points
  for( int p = 0; p < 5; ++p )
  {
    r[ p ][ 0 ] = p2_x[ p ] - p1_x[ p ];
    r[ p ][ 1 ] = p2_y[ p ] - p1_y[ p ];
    r[ p ][ 2 ] = p2_z[ p ] - p1_z[ p ];
  }

  // Computing the square of distance for the center
  HPFLTNB const r2[ 3 ] = {
    r[ 0 ][ 0 ] * r[ 0 ][ 0 ],
    r[ 0 ][ 1 ] * r[ 0 ][ 1 ],
    r[ 0 ][ 2 ] * r[ 0 ][ 2 ]
  };

  // Find the first and last intersecting plane using the parametric
  // values alpha_min and alpha_max
  HPFLTNB alpha_min = 0.0, alpha_max = 1.0;

  // Buffer storing the alpha on each axis
  // 2 different buffers to get the min and the max
  HPFLTNB alpha_x_min[ 4 ], alpha_x_max[ 4 ];
  HPFLTNB alpha_y_min[ 4 ], alpha_y_max[ 4 ];
  HPFLTNB alpha_z_min[ 4 ], alpha_z_max[ 4 ];

  // Position of the projected points in plane
  // If the main direction is X:
  // pos[0] is the position of point1 in Y
  // pos[1] is the position of point2 in Y
  // pos[2] is the position of point3 in Z
  // pos[3] is the position of point4 in Z
  HPFLTNB pos[ 4 ];

  // Width of intersection
  // w1 normalizes the overlap 1
  // w2 normalizes the overlap 2
  HPFLTNB w1, w2;

  // Buffer storing the min. and max. indices bounding the overlap
  int index[ 4 ];

  // Buffer storing the distances
  // 50 HPFLTNB is largely enough!!!
  HPFLTNB distance_x[ 50 ];
  HPFLTNB distance_y[ 50 ];
  HPFLTNB distance_z[ 50 ];

  FLTNB length_LOR = ap_ProjectionLine->GetLength();

  // Get TOF info
  HPFLTNB tof_resolutionT = ap_ProjectionLine->GetTOFResolution();
  HPFLTNB tof_measurement = ap_ProjectionLine->GetTOFMeasurement();

  // TOF standard deviation and truncation
  HPFLTNB tof_sigma = tof_resolutionT * SPEED_OF_LIGHT * 0.5 / TWO_SQRT_TWO_LN_2;
  HPFLTNB tof_half_span = tof_sigma * m_TOFnbSigmas;

  // convert delta time into delta length
  HPFLTNB tof_delta = tof_measurement * SPEED_OF_LIGHT * 0.5;

  // normalized Gaussian (integral=1) : precomputation of the coefficient that multiplies the exponential
  HPFLTNB gaussian_norm_coef = INV_SQRT_2_PI * 1./tof_sigma;
  
  // distance between the first event1 and the center of the Gaussian distribution along the LOR
  HPFLTNB lor_tof_center = length_LOR * 0.5 + tof_delta;

  // index along each axis of the first voxel falling inside the truncated Gaussian distribution
  HPFLTNB tof_edge_low[] = {0,0,0};
  // index along each axis of the last voxel falling inside the truncated Gaussian distribution
  HPFLTNB tof_edge_high[] = {0,0,0};
  HPFLTNB tof_center[] = {0,0,0};
  INTNB tof_index;

  // low/high voxel edges (in absolute coordinates) for truncated TOF
  for (int ax=0;ax<3;ax++)
  {
    // absolute coordinate along each axis of the center of the TOF distribution
    tof_center[ax] = p1[ax] +  lor_tof_center * r[0][ax] / length_LOR;

    // absolute coordinate along each axis of the lowest voxel edge spanned by the TOF distribution, limited by the lowest FOV edge
    tof_edge_low[ax] = tof_center[ax] - tof_half_span * fabs(r[0][ax]) / length_LOR;
    tof_index = max( (INTNB)::floor( (tof_edge_low[ax] - (-mp_halfFOV[ax])) / mp_sizeVox[ax] ), (INTNB)0);
    // if low TOF edge above the highest FOV edge, return empty line
    if (tof_index>mp_nbVox[ax]-1) return 0;
    tof_edge_low[ax] = (HPFLTNB)tof_index *  mp_sizeVox[ax] - mp_halfFOV[ax];

    // absolute coordinate along each axis of the highest voxel edge spanned by the TOF distribution, limited by the highest FOV edge
    tof_edge_high[ax] = tof_center[ax] + tof_half_span * fabs(r[0][ax]) / length_LOR;
    tof_index = min( (INTNB)::floor( (tof_edge_high[ax] - (-mp_halfFOV[ax])) / mp_sizeVox[ax] ), mp_nbVox[ax]-1);
    // if high TOF edge below the lowest FOV edge, return empty line
    if (tof_index<0) return 0;
    tof_edge_high[ax] = (HPFLTNB)(tof_index+1) * mp_sizeVox[ax] - mp_halfFOV[ax];
  }

  // Computing the angle of line vs. horizontal
  float const angle = std::acos( ( r[ 0 ][ 0 ] ) / ( std::sqrt( r2[ 0 ] + r2[ 1 ] ) ) ) * 180.0f / M_PI;

  // Condition on the largest component of r, taking the center
  if( ( angle >= 0.0 && angle <= 45.0 ) || ( angle >= 135.0 && angle <= 180.0 ) )
  {
    // Computing the parametric values alpha_min and alpha_max
    // taking into account TOF truncation
    // Because the main direction is X, we use only the center
    alpha_x_min[ 0 ] = ( tof_edge_low[ 0 ] - p1_x[ 0 ] ) / r[ 0 ][ 0 ];
    alpha_x_min[ 1 ] = ( tof_edge_high[ 0 ] - p1_x[ 0 ] ) / r[ 0 ][ 0 ];
    alpha_x_max[ 0 ] = alpha_x_min[ 0 ];
    alpha_x_max[ 1 ] = alpha_x_min[ 1 ];
    alpha_min = std::max( alpha_min, std::min( alpha_x_min[ 0 ], alpha_x_min[ 1 ] ) );
    alpha_max = std::min( alpha_max, std::max( alpha_x_max[ 0 ], alpha_x_max[ 1 ] ) );

    // Compute the parametric value for each line
    // For Y, we use only the point 1 and 2
    if( r[ 1 ][ 1 ] != 0 ) // point1
    {
      alpha_y_min[ 0 ] = ( tof_edge_low[ 1 ] - p1_y[ 1 ] ) / r[ 1 ][ 1 ];
      alpha_y_min[ 1 ] = ( tof_edge_high[ 1 ] - p1_y[ 1 ] ) / r[ 1 ][ 1 ];
      alpha_y_max[ 0 ] = alpha_y_min[ 0 ];
      alpha_y_max[ 1 ] = alpha_y_min[ 1 ];
    }
    else
    {
      alpha_y_min[ 0 ] = 0.0;
      alpha_y_min[ 1 ] = 0.0;
      alpha_y_max[ 0 ] = 1.0;
      alpha_y_max[ 1 ] = 1.0;
    }

    if( r[ 2 ][ 1 ] != 0 ) // point2
    {
      alpha_y_min[ 2 ] = ( tof_edge_low[ 1 ] - p1_y[ 2 ] ) / r[ 2 ][ 1 ];
      alpha_y_min[ 3 ] = ( tof_edge_high[ 1 ] - p1_y[ 2 ] ) / r[ 2 ][ 1 ];
      alpha_y_max[ 2 ] = alpha_y_min[ 2 ];
      alpha_y_max[ 3 ] = alpha_y_min[ 3 ];
    }
    else
    {
      alpha_y_min[ 2 ] = 0.0;
      alpha_y_min[ 3 ] = 0.0;
      alpha_y_max[ 2 ] = 1.0;
      alpha_y_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_y_min
    alpha_min = std::max( alpha_min,
      alpha_y_min[ std::minmax_element( alpha_y_min, alpha_y_min + 4 ).first - alpha_y_min ] );
    // Getting the maximum of alpha_y_max
    alpha_max = std::min( alpha_max,
      alpha_y_max[ std::minmax_element( alpha_y_max, alpha_y_max + 4 ).second - alpha_y_max ] );

    // For Z, we use only the point 3 and 4
    if( r[ 3 ][ 2 ] != 0 ) // point3
    {
      alpha_z_min[ 0 ] = ( tof_edge_low[ 2 ] - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_min[ 1 ] = ( tof_edge_high[ 2 ] - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_max[ 0 ] = alpha_z_min[ 0 ];
      alpha_z_max[ 1 ] = alpha_z_min[ 1 ];
    }
    else
    {
      alpha_z_min[ 0 ] = 0.0;
      alpha_z_min[ 1 ] = 0.0;
      alpha_z_max[ 0 ] = 1.0;
      alpha_z_max[ 1 ] = 1.0;
    }

    if( r[ 4 ][ 2 ] != 0 ) // point4
    {
      alpha_z_min[ 2 ] = ( tof_edge_low[ 2 ] - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_min[ 3 ] = ( tof_edge_high[ 2 ] - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_max[ 2 ] = alpha_z_min[ 2 ];
      alpha_z_max[ 3 ] = alpha_z_min[ 3 ];
    }
    else
    {
      alpha_z_min[ 2 ] = 0.0;
      alpha_z_min[ 3 ] = 0.0;
      alpha_z_max[ 2 ] = 1.0;
      alpha_z_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_z_min
    alpha_min = std::max( alpha_min,
      alpha_z_min[ std::minmax_element( alpha_z_min, alpha_z_min + 4 ).first
      - alpha_z_min ] );
    // Getting the maximum of alpha_z_max
    alpha_max = std::min( alpha_max,
      alpha_z_max[ std::minmax_element( alpha_z_max, alpha_z_max + 4 ).second
      - alpha_z_max ] );

    if( alpha_max <= alpha_min ) return 0;

    // Computing the first and the last plane
    int16_t i_min = 0, i_max = 0;
    if( r[ 0 ][ 0 ] > 0.0 )
    {
      i_min = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alpha_min * r[ 0 ][ 0 ] - p1_x[ 0 ] ) / mp_sizeVox[ 0 ] );
      i_max = ::floor( m_toleranceX + ( p1_x[ 0 ] + alpha_max * r[ 0 ][ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }
    else if( r[ 0 ][ 0 ] < 0.0 )
    {
      i_min = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alpha_max * r[ 0 ][ 0 ] - p1_x[ 0 ] ) / mp_sizeVox[ 0 ] );
      i_max = ::floor( m_toleranceX + ( p1_x[ 0 ] + alpha_min * r[ 0 ][ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }

    // Computing weight of normalization
    // Using the center
    HPFLTNB const factor( ::fabs( r[ 0 ][ 0 ] ) / ::sqrt( r2[ 0 ] + r2[ 1 ] + r2[ 2 ] ) );
    // removed because taken into account in TOF weight computation
    HPFLTNB weight( mp_sizeVox[ 0 ] / factor );
    HPFLTNB const factor_for_tof( length_LOR / r[ 0 ][ 0 ] );

    // Computing the increment for each 4 lines to project
    // (the center is not projected)
    for( uint8_t p = 0; p < 4; ++p )
    {
      ri[ p ][ 0 ] = 1.0; //r[ p + 1 ][ 0 ] / r[ p + 1 ][ 0 ]
      ri[ p ][ 1 ] = r[ p + 1 ][ 1 ] / r[ p + 1 ][ 0 ];
      ri[ p ][ 2 ] = r[ p + 1 ][ 2 ] / r[ p + 1 ][ 0 ];
    }

    // Increment orientation
    int8_t incr_orient_trans = 0;
    int8_t idx_orient_trans = 0;

    //int8_t const incr_orient_trans = p2_y[ 1 ] < p2_y[ 2 ] ? 1 : -1;
    int8_t const incr_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : -1;

    // Index orientation
    //int8_t const idx_orient_trans = p2_y[ 1 ] < p2_y[ 2 ] ? 1 : 0;
    int8_t const idx_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : 0;

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = (-mp_halfFOV[ 0 ]) + ( mp_sizeVox[ 0 ] * 0.5 ) - p1_x[ 0 ];

    // Loop over the crossed X planes
    for( int16_t i = i_min; i < i_max; ++i )
    {
      // Computing the coordinates of crossed plane
      HPFLTNB const step = offset + i * mp_sizeVox[ 0 ];
      // in Y (point 1 and 2)
      pos[ 0 ] = p1_y[ 1 ] + step * ri[ 0 ][ 1 ];
      pos[ 1 ] = p1_y[ 2 ] + step * ri[ 1 ][ 1 ];
      // in Z (point 3 and 4)
      pos[ 2 ] = p1_z[ 3 ] + step * ri[ 2 ][ 2 ];
      pos[ 3 ] = p1_z[ 4 ] + step * ri[ 3 ][ 2 ];

      // Computing the factor w1 and w2 normalizing the overlaps
      w1 = ::fabs( pos[ 0 ] - pos[ 1 ] );
      w2 = ::fabs( pos[ 2 ] - pos[ 3 ] );
      HPFLTNB final_weight = 0.0;
      if( ( w1 * w2 ) != 0 ) final_weight = weight / ( w1 * w2 );

      // Computing the index min and max on each axis
      // In Y
      index[ 0 ] = ::floor( m_toleranceY + ( pos[ 0 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
      index[ 1 ] = ::floor( m_toleranceY + ( pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
      // In Z
      index[ 2 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );
      index[ 3 ] = ::floor( m_toleranceZ + ( pos[ 3 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );

      // apply mask if required
      int test_index1 = i + index[ 0 ] * mp_nbVox[ 0 ] + index[ 2 ] * m_nbVoxXY;
      int test_index2 = i + index[ 1 ] * mp_nbVox[ 0 ] + index[ 3 ] * m_nbVoxXY;
      bool test_index = (index[ 0 ] >= 0 && index[ 2 ] >= 0 && index[ 1 ] >= 0 && index[ 3 ] >= 0
                      && index[ 0 ] < mp_nbVox[ 1 ] && index[ 1 ] < mp_nbVox[ 1 ] && index[ 2 ] < mp_nbVox[ 2 ] && index[ 3 ] < mp_nbVox[ 2 ]);
      if (m_hasMask && test_index && !mp_mask[test_index1] && !mp_mask[test_index2]) continue;

      if( index[ 0 ] < index[ 1 ] )
      {
        incr_orient_trans = 1;
        idx_orient_trans = 1;
      }
      else if( index[ 0 ] > index[ 1 ] )
      {
        incr_orient_trans = -1;
        idx_orient_trans = 0;
      }

      // Getting the number of distance in Y
      int16_t const n_distance_y = ::abs( index[ 0 ] - index[ 1 ] ) + 2;
      // Computing the distances in Y
      distance_y[ 0 ] = pos[ 0 ];
      distance_y[ n_distance_y - 1 ] = pos[ 1 ];
      // Computing the rest if necessary
      if( n_distance_y > 2 )
      {
        for( int16_t d = 0; d < n_distance_y - 2; ++d )
          distance_y[ d + 1 ] = (-mp_halfFOV[ 1 ]) + ( index[ 0 ] + idx_orient_trans + ( incr_orient_trans * d ) ) * mp_sizeVox[ 1 ];
      }

      // Computing the final distance in Y
      for( int16_t d = 0; d < n_distance_y - 1; ++d )
        distance_y[ d ] = ::fabs( distance_y[ d + 1 ] - distance_y[ d ] );

      // Getting the number of distance in Z
      int16_t const n_distance_z = ::abs( index[ 2 ] - index[ 3 ] ) + 2;
      // Storing the positions in Y
      distance_z[ 0 ] = pos[ 2 ];
      distance_z[ n_distance_z - 1 ] = pos[ 3 ];
      // Storing the rest if necessary
      if( n_distance_z > 2 )
      {
        for( int16_t d = 0; d < n_distance_z - 2; ++d )
          distance_z[ d + 1 ] = (-mp_halfFOV[ 2 ]) + ( index[ 2 ] + idx_orient_axial + ( incr_orient_axial * d ) ) * mp_sizeVox[ 2 ];
      }

      // Computing the final distance in Z
      for( int16_t d = 0; d < n_distance_z - 1; ++d )
        distance_z[ d ] = ::fabs( distance_z[ d + 1 ] - distance_z[ d ] );

      // TOF coefficient = value of the Gaussian at voxel center (projected onto the LOR) * voxel size (projected onto the LOR)
      HPFLTNB temp = (step * factor_for_tof - lor_tof_center) / tof_sigma;
      HPFLTNB tof_weight =  exp(- 0.5 * temp * temp ) * gaussian_norm_coef;

      // Loop over the overlap and store the elements
      int16_t index_y, index_z;
      for( int16_t jj = 0; jj < n_distance_z - 1; ++jj )
      {
        index_z = index[ 2 ] + jj * incr_orient_axial;
        if( index_z < 0 || index_z > mp_nbVox[ 2 ] - 1 ) continue;

        for( int16_t ii = 0; ii < n_distance_y - 1; ++ii )
        {
          index_y = index[ 0 ] + ii * incr_orient_trans;
          if( index_y < 0 || index_y > mp_nbVox[ 1 ] - 1 ) continue;

          ap_ProjectionLine->AddVoxelInTOFBin(a_direction, 0, i + index_y * mp_nbVox[ 0 ] + index_z * m_nbVoxXY, final_weight * distance_z[ jj ] * distance_y[ ii ] * tof_weight);
        }
      }
    }
  }
  else
  {
    // Compute the parametric value for each line
    // taking into account TOF truncation
    // For X, we use only the point 1 and 2
    if( r[ 1 ][ 0 ] != 0.0 ) // point1
    {
      alpha_x_min[ 0 ] = ( tof_edge_low[ 0 ] - p1_x[ 1 ] ) / r[ 1 ][ 0 ];
      alpha_x_min[ 1 ] = ( tof_edge_high[ 0 ] - p1_x[ 1 ] ) / r[ 1 ][ 0 ];
      alpha_x_max[ 0 ] = alpha_x_min[ 0 ];
      alpha_x_max[ 1 ] = alpha_x_min[ 1 ];
    }
    else
    {
      alpha_x_min[ 0 ] = 0.0;
      alpha_x_min[ 1 ] = 0.0;
      alpha_x_max[ 0 ] = 1.0;
      alpha_x_max[ 1 ] = 1.0;
    }

    if( r[ 2 ][ 0 ] != 0 ) // point2
    {
      alpha_x_min[ 2 ] = ( tof_edge_low[ 0 ] - p1_x[ 2 ] ) / r[ 2 ][ 0 ];
      alpha_x_min[ 3 ] = ( tof_edge_high[ 0 ] - p1_x[ 2 ] ) / r[ 2 ][ 0 ];
      alpha_x_max[ 2 ] = alpha_x_min[ 2 ];
      alpha_x_max[ 3 ] = alpha_x_min[ 3 ];
    }
    else
    {
      alpha_x_min[ 2 ] = 0.0;
      alpha_x_min[ 3 ] = 0.0;
      alpha_x_max[ 2 ] = 1.0;
      alpha_x_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_x_min
    alpha_min = std::max( alpha_min,
      alpha_x_min[ std::minmax_element( alpha_x_min, alpha_x_min + 4 ).first - alpha_x_min ] );
    // Getting the maximum of alpha_y_max
    alpha_max = std::min( alpha_max,
      alpha_x_max[ std::minmax_element( alpha_x_max, alpha_x_max + 4 ).second - alpha_x_max ] );

    // Computing the parametric values alpha_min and alpha_max
    // Because the main direction is Y, we use only the center
    alpha_y_min[ 0 ] = ( tof_edge_low[1] - p1_y[ 0 ] ) / r[ 0 ][ 1 ];
    alpha_y_min[ 1 ] = ( tof_edge_high[1] - p1_y[ 0 ] ) / r[ 0 ][ 1 ];
    alpha_y_max[ 0 ] = alpha_y_min[ 0 ];
    alpha_y_max[ 1 ] = alpha_y_min[ 1 ];
    alpha_min = std::max( alpha_min,
      std::min( alpha_y_min[ 0 ], alpha_y_min[ 1 ] ) );
    alpha_max = std::min( alpha_max,
      std::max( alpha_y_max[ 0 ], alpha_y_max[ 1 ] ) );

    // For Z, we use only the point 3 and 4
    if( r[ 3 ][ 2 ] != 0 ) // point3
    {
      alpha_z_min[ 0 ] = ( tof_edge_low[ 2 ] - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_min[ 1 ] = ( tof_edge_high[ 2 ]- p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_max[ 0 ] = alpha_z_min[ 0 ];
      alpha_z_max[ 1 ] = alpha_z_min[ 1 ];
    }
    else
    {
      alpha_z_min[ 0 ] = 0.0;
      alpha_z_min[ 1 ] = 0.0;
      alpha_z_max[ 0 ] = 1.0;
      alpha_z_max[ 1 ] = 1.0;
    }

    if( r[ 4 ][ 2 ] != 0 ) // point4
    {
      alpha_z_min[ 2 ] = ( tof_edge_low[ 2 ] - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_min[ 3 ] = ( tof_edge_high[ 2 ] - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_max[ 2 ] = alpha_z_min[ 2 ];
      alpha_z_max[ 3 ] = alpha_z_min[ 3 ];
    }
    else
    {
      alpha_z_min[ 2 ] = 0.0;
      alpha_z_min[ 3 ] = 0.0;
      alpha_z_max[ 2 ] = 1.0;
      alpha_z_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_z_min
    alpha_min = std::max( alpha_min, alpha_z_min[ std::minmax_element( alpha_z_min, alpha_z_min + 4 ).first - alpha_z_min ] );
    // Getting the maximum of alpha_z_max
    alpha_max = std::min( alpha_max, alpha_z_max[ std::minmax_element( alpha_z_max, alpha_z_max + 4 ).second - alpha_z_max ] );

    if( alpha_max <= alpha_min ) return 0;

    // Computing the first and the last plane
    int16_t j_min = 0, j_max = 0;
    if( r[ 0 ][ 1 ] > 0.0 )
    {
      j_min = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alpha_min * r[ 0 ][ 1 ] - p1_y[ 0 ] ) / mp_sizeVox[ 1 ] );
      j_max = ::floor( m_toleranceY + ( p1_y[ 0 ] + alpha_max * r[ 0 ][ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }
    else if( r[ 0 ][ 1 ] < 0.0 )
    {
      j_min = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alpha_max * r[ 0 ][ 1 ] - p1_y[ 0 ] ) / mp_sizeVox[ 1 ] );
      j_max = ::floor( m_toleranceY + ( p1_y[ 0 ] + alpha_min * r[ 0 ][ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }

    // Computing weight of normalization
    // Using the center
    HPFLTNB const factor( ::fabs( r[ 0 ][ 1 ] ) / ::sqrt( r2[ 0 ] + r2[ 1 ] + r2[ 2 ] ) );
    HPFLTNB weight( mp_sizeVox[ 1 ] / factor );
    HPFLTNB const factor_for_tof( length_LOR / r[ 0 ][ 1 ] );

    // Computing the increment for each 4 lines to project
    // (the center is not projected)
    for( uint8_t p = 0; p < 4; ++p )
    {
      ri[ p ][ 0 ] = r[ p + 1 ][ 0 ] / r[ p + 1 ][ 1 ];
      ri[ p ][ 1 ] = 1.0;
      ri[ p ][ 2 ] = r[ p + 1 ][ 2 ] / r[ p + 1 ][ 1 ];
    }

    // Increment orientation and Index orientation
    int8_t incr_orient_trans = 0;
    int8_t idx_orient_trans = 0;

    //int8_t const incr_orient_trans = p2_x[ 1 ] < p2_x[ 2 ] ? 1 : -1;
    int8_t const incr_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : -1;
    //int8_t const idx_orient_trans = p2_x[ 1 ] < p2_x[ 2 ] ? 1 : 0;
    int8_t const idx_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : 0;

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = (-mp_halfFOV[ 1 ]) + ( mp_sizeVox[ 1 ] * 0.5 ) - p1_y[ 0 ];

    // Loop over the crossed Y planes
    for( int16_t j = j_min; j < j_max; ++j )
    {
      // Computing the coordinates of crossed plane
      HPFLTNB const step = offset + j * mp_sizeVox[ 1 ];
      // in Y (point 1 and 2)
      pos[ 0 ] = p1_x[ 1 ] + step * ri[ 0 ][ 0 ];
      pos[ 1 ] = p1_x[ 2 ] + step * ri[ 1 ][ 0 ];
      // in Z (point 3 and 4)
      pos[ 2 ] = p1_z[ 3 ] + step * ri[ 2 ][ 2 ];
      pos[ 3 ] = p1_z[ 4 ] + step * ri[ 3 ][ 2 ];

      // Computing the factor w1 and w2 normalizing the overlaps
      w1 = ::fabs( pos[ 0 ] - pos[ 1 ] );
      w2 = ::fabs( pos[ 2 ] - pos[ 3 ] );
      HPFLTNB final_weight = 0.0;
      if( ( w1 * w2 ) != 0 ) final_weight = weight / ( w1 * w2 );

      // Computing the index min and max on each axis
      // In Y
      index[ 0 ] = ::floor( m_toleranceX + ( pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
      index[ 1 ] = ::floor( m_toleranceX + ( pos[ 1 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
      // In Z
      index[ 2 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );
      index[ 3 ] = ::floor( m_toleranceZ + ( pos[ 3 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );

      // if voxel masked, skip
      int test_index1 = index[ 0 ] + j* mp_nbVox[ 0 ] + index[ 2 ] * m_nbVoxXY;
      int test_index2 = index[ 1 ] + j* mp_nbVox[ 0 ] + index[ 3 ] * m_nbVoxXY;
      bool test_index = (index[ 0 ] >= 0 && index[ 2 ] >= 0 && index[ 1 ] >= 0 && index[ 3 ] >= 0
                      && index[ 0 ] < mp_nbVox[ 0 ] && index[ 1 ] < mp_nbVox[ 0 ] && index[ 2 ] < mp_nbVox[ 2 ] && index[ 3 ] < mp_nbVox[ 2 ]);
      if (m_hasMask && test_index && !mp_mask[test_index1] && !mp_mask[test_index2]) continue;

      if( index[ 0 ] < index[ 1 ] )
      {
        incr_orient_trans = 1;
        idx_orient_trans = 1;
      }
      else if( index[ 0 ] > index[ 1 ] )
      {
        incr_orient_trans = -1;
        idx_orient_trans = 0;
      }

      // Getting the number of distance in X
      int16_t const n_distance_x = ::abs( index[ 0 ] - index[ 1 ] ) + 2;
      // Computing the distances in X
      distance_x[ 0 ] = pos[ 0 ];
      distance_x[ n_distance_x - 1 ] = pos[ 1 ];
      // Computing the rest if necessary
      if( n_distance_x > 2 )
      {
        for( int16_t d = 0; d < n_distance_x - 2; ++d )
          distance_x[ d + 1 ] = (-mp_halfFOV[ 0 ]) + ( index[ 0 ] + idx_orient_trans + ( incr_orient_trans * d ) ) * mp_sizeVox[ 0 ];
      }

      // Computing the final distance in X
      for( int16_t d = 0; d < n_distance_x - 1; ++d )
        distance_x[ d ] = ::fabs( distance_x[ d + 1 ] - distance_x[ d ] );

      // Getting the number of distance in Z
      int16_t const n_distance_z = ::abs( index[ 2 ] - index[ 3 ] ) + 2;
      // Storing the positions in Y
      distance_z[ 0 ] = pos[ 2 ];
      distance_z[ n_distance_z - 1 ] = pos[ 3 ];
      // Storing the rest if necessary
      if( n_distance_z > 2 )
      {
        for( int16_t d = 0; d < n_distance_z - 2; ++d )
          distance_z[ d + 1 ] = (-mp_halfFOV[ 2 ]) + ( index[ 2 ] + idx_orient_axial + ( incr_orient_axial * d ) ) * mp_sizeVox[ 2 ];
      }

      // Computing the final distance in Z
      for( int16_t d = 0; d < n_distance_z - 1; ++d )
        distance_z[ d ] = ::fabs( distance_z[ d + 1 ] - distance_z[ d ] );

      // TOF coefficient = value of the Gaussian at voxel center (projected onto the LOR) * voxel size (projected onto the LOR)
      HPFLTNB temp = (step * factor_for_tof - lor_tof_center) / tof_sigma;
      HPFLTNB tof_weight =  exp(- 0.5 * temp * temp ) * gaussian_norm_coef;

      // Loop over the overlap and store the elements
      int16_t index_x, index_z;
      for( int16_t jj = 0; jj < n_distance_z - 1; ++jj )
      {
        index_z = index[ 2 ] + jj * incr_orient_axial;
        if( index_z < 0 || index_z > mp_nbVox[ 2 ] - 1 ) continue;

        for( int16_t ii = 0; ii < n_distance_x - 1; ++ii )
        {
          index_x = index[ 0 ] + ii * incr_orient_trans;
          if( index_x < 0 || index_x > mp_nbVox[ 0 ] - 1 ) continue;

          ap_ProjectionLine->AddVoxelInTOFBin(a_direction, 0, index_x + j * mp_nbVox[ 0 ] + index_z * m_nbVoxXY, final_weight * distance_z[ jj ] * distance_x[ ii ] * tof_weight);
        }
      }
    }
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorDistanceDriven::ProjectWithTOFBin(int a_direction, oProjectionLine* ap_ProjectionLine)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorDistanceDriven::ProjectWithTOFBin() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorDistanceDriven::Project with TOF bins -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Computing the coordinates of the points on the bound of the source/crystal
  // 1. The first point is on transaxial plan and the last point is on axial
  // plan.
  // Source spot size (CT) / Crystal 1 (PET):
  //
  //                  ------- P13 -------
  //                  |                  |
  //                  |                  |
  //                 P12      Center     P11
  //                  |         0        |
  //                  |                  |
  //                  ------- P14 -------
  //  Y
  //  <---------. X
  //            |
  //            |
  //            |
  //           \/ Z
  // --------------

  // Computing the position of the point2.
  // In Distance Driven the point p21 relies the point P11, p22 -> p12,
  // p23 -> p13 and p24 -> p14
  // Pixel (CT) / Crystal 2 (PET):
  //
  //                  ------- P23 -------
  //                  |                  |
  //                  |                  |
  //                 P22      Center     P21
  //                  |         0        |
  //                  |                  |
  //                  ------- P24 -------
  //  Y
  //  <---------. X
  //            |
  //            |
  //            |
  //           \/ Z
  // --------------


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
  HPFLTNB tof_weight;

  // length LOR, connecting central crystal points
  HPFLTNB length_LOR = ap_ProjectionLine->GetLength();
  HPFLTNB length_LOR_half = length_LOR * 0.5;

  INTNB tof_bin_first_for_voxel = 0, tof_bin_last_for_voxel = 0;

  // buffers storing point 1 and 2
  FLTNB p1_x[ 5 ];
  FLTNB p1_y[ 5 ];
  FLTNB p1_z[ 5 ];
  FLTNB p2_x[ 5 ];
  FLTNB p2_y[ 5 ];
  FLTNB p2_z[ 5 ];

  // Get the positions of the centers of the edges of crystal elements or source position
  if (mp_Scanner->GetEdgesCenterPositions(
        ap_ProjectionLine->GetIndex1(), // Index 1
        ap_ProjectionLine->GetIndex2(), // Index 2
        ap_ProjectionLine->GetPosition1(), // Line position for point 1
        ap_ProjectionLine->GetPosition2(), // Line position for point 2
        &p1_x[1], &p1_y[1], &p1_z[1], // Edges for point 1
        &p2_x[1], &p2_y[1], &p2_z[1] // Edges for point 2
    ))
  {
    Cerr("***** iProjectorDistanceDriven::ProjectWithoutTOF() -> A problem occured while getting the edges' center positions !" << endl);
    return 1;
  }

  // Point 1 focal (CT) or crystal (PET) position
  FLTNB* p1 = ap_ProjectionLine->GetPosition1();

  // Point 2 pixel (CT) or crystal (PET) position
  FLTNB* p2 = ap_ProjectionLine->GetPosition2();

  // Center position p1
  p1_x[ 0 ] = p1[ 0 ]; p1_y[ 0 ] = p1[ 1 ]; p1_z[ 0 ] = p1[ 2 ];

  // Center position p2
  p2_x[ 0 ] = p2[ 0 ]; p2_y[ 0 ] = p2[ 1 ]; p2_z[ 0 ] = p2[ 2 ];

  // storing the ray connecting p1 and p2
  HPFLTNB r[ 5 ][ 3 ];

  // buffer storing the intersection
  HPFLTNB ri[ 4 ][ 3 ];

  // Take the position of ray connecting p1 and p2 for all the 5 points
  for( int p = 0; p < 5; ++p )
  {
    r[ p ][ 0 ] = p2_x[ p ] - p1_x[ p ];
    r[ p ][ 1 ] = p2_y[ p ] - p1_y[ p ];
    r[ p ][ 2 ] = p2_z[ p ] - p1_z[ p ];
  }

  // Computing the square of distance for the center
  HPFLTNB const r2[ 3 ] = {
    r[ 0 ][ 0 ] * r[ 0 ][ 0 ],
    r[ 0 ][ 1 ] * r[ 0 ][ 1 ],
    r[ 0 ][ 2 ] * r[ 0 ][ 2 ]
  };

  // Find the first and last intersecting plane using the parametric
  // values alpha_min and alpha_max
  HPFLTNB alpha_min = 0.0, alpha_max = 1.0;

  // Buffer storing the alpha on each axis
  // 2 different buffers to get the min and the max
  HPFLTNB alpha_x_min[ 4 ], alpha_x_max[ 4 ];
  HPFLTNB alpha_y_min[ 4 ], alpha_y_max[ 4 ];
  HPFLTNB alpha_z_min[ 4 ], alpha_z_max[ 4 ];

  // Position of the projected points in plane
  // If the main direction is X:
  // pos[0] is the position of point1 in Y
  // pos[1] is the position of point2 in Y
  // pos[2] is the position of point3 in Z
  // pos[3] is the position of point4 in Z
  HPFLTNB pos[ 4 ];

  // Width of intersection
  // w1 normalizes the overlap 1
  // w2 normalizes the overlap 2
  HPFLTNB w1, w2;

  // Buffer storing the min. and max. indices bounding the overlap
  int index[ 4 ];

  // Buffer storing the distances
  // 50 HPFLTNB is largely enough!!!
  HPFLTNB distance_x[ 50 ];
  HPFLTNB distance_y[ 50 ];
  HPFLTNB distance_z[ 50 ];

  // Computing the angle of line vs. horizontal
  float const angle = std::acos( ( r[ 0 ][ 0 ] ) / ( std::sqrt( r2[ 0 ] + r2[ 1 ] ) ) ) * 180.0f / M_PI;

  // Condition on the largest component of r, taking the center
  if( ( angle >= 0.0 && angle <= 45.0 ) || ( angle >= 135.0 && angle <= 180.0 ) )
  {
    // Computing the parametric values alpha_min and alpha_max
    // Because the main direction is X, we use only the center
    alpha_x_min[ 0 ] = ( (-mp_halfFOV[ 0 ]) - p1_x[ 0 ] ) / r[ 0 ][ 0 ];
    alpha_x_min[ 1 ] = ( mp_halfFOV[ 0 ] - p1_x[ 0 ] ) / r[ 0 ][ 0 ];
    alpha_x_max[ 0 ] = alpha_x_min[ 0 ];
    alpha_x_max[ 1 ] = alpha_x_min[ 1 ];
    alpha_min = std::max( alpha_min, std::min( alpha_x_min[ 0 ], alpha_x_min[ 1 ] ) );
    alpha_max = std::min( alpha_max, std::max( alpha_x_max[ 0 ], alpha_x_max[ 1 ] ) );

    // Compute the parametric value for each line
    // For Y, we use only the point 1 and 2
    if( r[ 1 ][ 1 ] != 0.0 ) // point1
    {
      alpha_y_min[ 0 ] = ( (-mp_halfFOV[ 1 ]) - p1_y[ 1 ] ) / r[ 1 ][ 1 ];
      alpha_y_min[ 1 ] = ( mp_halfFOV[ 1 ] - p1_y[ 1 ] ) / r[ 1 ][ 1 ];
      alpha_y_max[ 0 ] = alpha_y_min[ 0 ];
      alpha_y_max[ 1 ] = alpha_y_min[ 1 ];
    }
    else
    {
      alpha_y_min[ 0 ] = 0.0;
      alpha_y_min[ 1 ] = 0.0;
      alpha_y_max[ 0 ] = 1.0;
      alpha_y_max[ 1 ] = 1.0;
    }

    if( r[ 2 ][ 1 ] != 0.0 ) // point2
    {
      alpha_y_min[ 2 ] = ( (-mp_halfFOV[ 1 ]) - p1_y[ 2 ] ) / r[ 2 ][ 1 ];
      alpha_y_min[ 3 ] = ( mp_halfFOV[ 1 ] - p1_y[ 2 ] ) / r[ 2 ][ 1 ];
      alpha_y_max[ 2 ] = alpha_y_min[ 2 ];
      alpha_y_max[ 3 ] = alpha_y_min[ 3 ];
    }
    else
    {
      alpha_y_min[ 2 ] = 0.0;
      alpha_y_min[ 3 ] = 0.0;
      alpha_y_max[ 2 ] = 1.0;
      alpha_y_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_y_min
    alpha_min = std::max( alpha_min,
      alpha_y_min[ std::minmax_element( alpha_y_min, alpha_y_min + 4 ).first - alpha_y_min ] );
    // Getting the maximum of alpha_y_max
    alpha_max = std::min( alpha_max,
      alpha_y_max[ std::minmax_element( alpha_y_max, alpha_y_max + 4 ).second - alpha_y_max ] );

    // For Z, we use only the point 3 and 4
    if( r[ 3 ][ 2 ] != 0.0 ) // point3
    {
      alpha_z_min[ 0 ] = ( (-mp_halfFOV[ 2 ]) - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_min[ 1 ] = ( mp_halfFOV[ 2 ] - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_max[ 0 ] = alpha_z_min[ 0 ];
      alpha_z_max[ 1 ] = alpha_z_min[ 1 ];
    }
    else
    {
      alpha_z_min[ 0 ] = 0.0;
      alpha_z_min[ 1 ] = 0.0;
      alpha_z_max[ 0 ] = 1.0;
      alpha_z_max[ 1 ] = 1.0;
    }

    if( r[ 4 ][ 2 ] != 0.0 ) // point4
    {
      alpha_z_min[ 2 ] = ( (-mp_halfFOV[ 2 ]) - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_min[ 3 ] = ( mp_halfFOV[ 2 ] - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_max[ 2 ] = alpha_z_min[ 2 ];
      alpha_z_max[ 3 ] = alpha_z_min[ 3 ];
    }
    else
    {
      alpha_z_min[ 2 ] = 0.0;
      alpha_z_min[ 3 ] = 0.0;
      alpha_z_max[ 2 ] = 1.0;
      alpha_z_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_z_min
    alpha_min = std::max( alpha_min,
      alpha_z_min[ std::minmax_element( alpha_z_min, alpha_z_min + 4 ).first
      - alpha_z_min ] );
    // Getting the maximum of alpha_z_max
    alpha_max = std::min( alpha_max,
      alpha_z_max[ std::minmax_element( alpha_z_max, alpha_z_max + 4 ).second
      - alpha_z_max ] );

    if( alpha_max <= alpha_min ) return 0;

    // temporary storage for TOF bins Gaussian integrals over the current voxel
    // allocation after potential returns
    HPFLTNB* tof_weights_temp = new HPFLTNB[tof_nb_bins];
    for (INTNB tof_bin=0; tof_bin<tof_nb_bins; tof_bin++) tof_weights_temp[tof_bin] = 0.;

    // Computing the first and the last plane
    int16_t i_min = 0, i_max = 0;
    if( r[ 0 ][ 0 ] > 0.0 )
    {
      i_min = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alpha_min * r[ 0 ][ 0 ] - p1_x[ 0 ] ) / mp_sizeVox[ 0 ] );
      i_max = ::floor( m_toleranceX + ( p1_x[ 0 ] + alpha_max * r[ 0 ][ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }
    else if( r[ 0 ][ 0 ] < 0.0 )
    {
      i_min = ::ceil( mp_nbVox[ 0 ] - m_toleranceX - ( mp_halfFOV[ 0 ] - alpha_max * r[ 0 ][ 0 ] - p1_x[ 0 ] ) / mp_sizeVox[ 0 ] );
      i_max = ::floor( m_toleranceX + ( p1_x[ 0 ] + alpha_min * r[ 0 ][ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
    }

    // Computing weight of normalization
    // Using the center
    HPFLTNB const factor( ::fabs( r[ 0 ][ 0 ] ) / ::sqrt( r2[ 0 ] + r2[ 1 ] + r2[ 2 ] ) );
    HPFLTNB const factor_for_tof( length_LOR / r[ 0 ][ 0 ] );
    HPFLTNB weight( mp_sizeVox[ 0 ] / factor );

    // Computing the increment for each 4 lines to project
    // (the center is not projected)
    for( uint8_t p = 0; p < 4; ++p )
    {
      ri[ p ][ 0 ] = 1.0; //r[ p + 1 ][ 0 ] / r[ p + 1 ][ 0 ]
      ri[ p ][ 1 ] = r[ p + 1 ][ 1 ] / r[ p + 1 ][ 0 ];
      ri[ p ][ 2 ] = r[ p + 1 ][ 2 ] / r[ p + 1 ][ 0 ];
    }

    // Increment orientation
    int8_t incr_orient_trans = 0;
    int8_t idx_orient_trans = 0;

    //int8_t const incr_orient_trans = p2_y[ 1 ] < p2_y[ 2 ] ? 1 : -1;
    int8_t const incr_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : -1;

    // Index orientation
    //int8_t const idx_orient_trans = p2_y[ 1 ] < p2_y[ 2 ] ? 1 : 0;
    int8_t const idx_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : 0;

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = (-mp_halfFOV[ 0 ]) + ( mp_sizeVox[ 0 ] * 0.5 ) - p1_x[ 0 ];

    // Loop over the crossed X planes
    for( int16_t i = i_min; i < i_max; ++i )
    {
      // Computing the coordinates of crossed plane
      HPFLTNB const step = offset + i * mp_sizeVox[ 0 ];
      // in Y (point 1 and 2)
      pos[ 0 ] = p1_y[ 1 ] + step * ri[ 0 ][ 1 ];
      pos[ 1 ] = p1_y[ 2 ] + step * ri[ 1 ][ 1 ];
      // in Z (point 3 and 4)
      pos[ 2 ] = p1_z[ 3 ] + step * ri[ 2 ][ 2 ];
      pos[ 3 ] = p1_z[ 4 ] + step * ri[ 3 ][ 2 ];

      // Computing the factor w1 and w2 normalizing the overlaps
      w1 = ::fabs( pos[ 0 ] - pos[ 1 ] );
      w2 = ::fabs( pos[ 2 ] - pos[ 3 ] );
      HPFLTNB final_weight = 0.0;
      if( ( w1 * w2 ) != 0 ) final_weight = weight / ( w1 * w2 );

      // Computing the index min and max on each axis
      // In Y
      index[ 0 ] = ::floor( m_toleranceY + ( pos[ 0 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
      index[ 1 ] = ::floor( m_toleranceY + ( pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
      // In Z
      index[ 2 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );
      index[ 3 ] = ::floor( m_toleranceZ + ( pos[ 3 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );
 
      // apply mask if required
      int test_index1 = i + index[ 0 ] * mp_nbVox[ 0 ] + index[ 2 ] * m_nbVoxXY;
      int test_index2 = i + index[ 1 ] * mp_nbVox[ 0 ] + index[ 3 ] * m_nbVoxXY;
      bool test_index = (index[ 0 ] >= 0 && index[ 2 ] >= 0 && index[ 1 ] >= 0 && index[ 3 ] >= 0
                      && index[ 0 ] < mp_nbVox[ 1 ] && index[ 1 ] < mp_nbVox[ 1 ] && index[ 2 ] < mp_nbVox[ 2 ] && index[ 3 ] < mp_nbVox[ 2 ]);
      if (m_hasMask && test_index && !mp_mask[test_index1] && !mp_mask[test_index2]) continue;

      if( index[ 0 ] < index[ 1 ] )
      {
        incr_orient_trans = 1;
        idx_orient_trans = 1;
      }
      else if( index[ 0 ] > index[ 1 ] )
      {
        incr_orient_trans = -1;
        idx_orient_trans = 0;
      }

      // Getting the number of distance in Y
      int16_t const n_distance_y = ::abs( index[ 0 ] - index[ 1 ] ) + 2;
      // Computing the distances in Y
      distance_y[ 0 ] = pos[ 0 ];
      distance_y[ n_distance_y - 1 ] = pos[ 1 ];
      // Computing the rest if necessary
      if( n_distance_y > 2 )
      {
        for( int16_t d = 0; d < n_distance_y - 2; ++d )
          distance_y[ d + 1 ] = (-mp_halfFOV[ 1 ]) + ( index[ 0 ] + idx_orient_trans + ( incr_orient_trans * d ) ) * mp_sizeVox[ 1 ];
      }

      // Computing the final distance in Y
      for( int16_t d = 0; d < n_distance_y - 1; ++d )
        distance_y[ d ] = ::fabs( distance_y[ d + 1 ] - distance_y[ d ] );

      // Getting the number of distance in Z
      int16_t const n_distance_z = ::abs( index[ 2 ] - index[ 3 ] ) + 2;
      // Storing the positions in Y
      distance_z[ 0 ] = pos[ 2 ];
      distance_z[ n_distance_z - 1 ] = pos[ 3 ];
      // Storing the rest if necessary
      if( n_distance_z > 2 )
      {
        for( int16_t d = 0; d < n_distance_z - 2; ++d )
          distance_z[ d + 1 ] = (-mp_halfFOV[ 2 ]) + ( index[ 2 ] + idx_orient_axial + ( incr_orient_axial * d ) ) * mp_sizeVox[ 2 ];
      }

      // Computing the final distance in Z
      for( int16_t d = 0; d < n_distance_z - 1; ++d )
        distance_z[ d ] = ::fabs( distance_z[ d + 1 ] - distance_z[ d ] );

      // Compute the first and the last relevant TOF bin for this voxel
      tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(trunc( (step * factor_for_tof - tof_half_span - length_LOR_half ) / tof_bin_size )));
      tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(trunc( (step * factor_for_tof + tof_half_span - length_LOR_half) / tof_bin_size )));
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
        // TOF coefficient = value of the Gaussian at voxel center (projected onto the LOR) * voxel size (projected onto the LOR)
        HPFLTNB temp = (step * factor_for_tof - lor_tof_center) / tof_sigma;
        tof_weight =  exp(- 0.5 * temp * temp ) * gaussian_norm_coef;
        // add the weight to the sum for normalization
        tof_norm_coef += tof_weight;
        // save the weight temporarily
        tof_weights_temp[tof_bin] = tof_weight;
        // update TOF center along the LOR for the next TOF bin
        lor_tof_center += tof_bin_size;
      }

      // compute and write the final TOF bin coefficients for current voxels
      if (tof_norm_coef>0.)
      {
        for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
        {
          // normalized tof coef
          tof_weights_temp[tof_bin] /= tof_norm_coef;
        }

        // Loop over the overlap and store the elements
        int16_t index_y, index_z;
        for( int16_t jj = 0; jj < n_distance_z - 1; ++jj )
        {
          index_z = index[ 2 ] + jj * incr_orient_axial;
          if( index_z < 0 || index_z > mp_nbVox[ 2 ] - 1 ) continue;

          for( int16_t ii = 0; ii < n_distance_y - 1; ++ii )
          {
            index_y = index[ 0 ] + ii * incr_orient_trans;
            if( index_y < 0 || index_y > mp_nbVox[ 1 ] - 1 ) continue;

            ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, i + index_y * mp_nbVox[ 0 ] + index_z * m_nbVoxXY, final_weight * distance_z[ jj ] * distance_y[ ii ], tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
          }
        }
      }
    }
    delete[] tof_weights_temp;
  }
  else
  {
    // Compute the parametric value for each line
    // For X, we use only the point 1 and 2
    if( r[ 1 ][ 0 ] != 0.0 ) // point1
    {
      alpha_x_min[ 0 ] = ( (-mp_halfFOV[ 0 ]) - p1_x[ 1 ] ) / r[ 1 ][ 0 ];
      alpha_x_min[ 1 ] = ( mp_halfFOV[ 0 ] - p1_x[ 1 ] ) / r[ 1 ][ 0 ];
      alpha_x_max[ 0 ] = alpha_x_min[ 0 ];
      alpha_x_max[ 1 ] = alpha_x_min[ 1 ];
    }
    else
    {
      alpha_x_min[ 0 ] = 0.0;
      alpha_x_min[ 1 ] = 0.0;
      alpha_x_max[ 0 ] = 1.0;
      alpha_x_max[ 1 ] = 1.0;
    }

    if( r[ 2 ][ 0 ] != 0.0 ) // point2
    {
      alpha_x_min[ 2 ] = ( (-mp_halfFOV[ 0 ]) - p1_x[ 2 ] ) / r[ 2 ][ 0 ];
      alpha_x_min[ 3 ] = ( mp_halfFOV[ 0 ] - p1_x[ 2 ] ) / r[ 2 ][ 0 ];
      alpha_x_max[ 2 ] = alpha_x_min[ 2 ];
      alpha_x_max[ 3 ] = alpha_x_min[ 3 ];
    }
    else
    {
      alpha_x_min[ 2 ] = 0.0;
      alpha_x_min[ 3 ] = 0.0;
      alpha_x_max[ 2 ] = 1.0;
      alpha_x_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_x_min
    alpha_min = std::max( alpha_min,
      alpha_x_min[ std::minmax_element( alpha_x_min, alpha_x_min + 4 ).first - alpha_x_min ] );
    // Getting the maximum of alpha_y_max
    alpha_max = std::min( alpha_max,
      alpha_x_max[ std::minmax_element( alpha_x_max, alpha_x_max + 4 ).second - alpha_x_max ] );

    // Computing the parametric values alpha_min and alpha_max
    // Because the main direction is Y, we use only the center
    alpha_y_min[ 0 ] = ( (-mp_halfFOV[ 1 ]) - p1_y[ 0 ] ) / r[ 0 ][ 1 ];
    alpha_y_min[ 1 ] = ( mp_halfFOV[ 1 ] - p1_y[ 0 ] ) / r[ 0 ][ 1 ];
    alpha_y_max[ 0 ] = alpha_y_min[ 0 ];
    alpha_y_max[ 1 ] = alpha_y_min[ 1 ];
    alpha_min = std::max( alpha_min,
      std::min( alpha_y_min[ 0 ], alpha_y_min[ 1 ] ) );
    alpha_max = std::min( alpha_max,
      std::max( alpha_y_max[ 0 ], alpha_y_max[ 1 ] ) );

    // For Z, we use only the point 3 and 4
    if( r[ 3 ][ 2 ] != 0.0 ) // point3
    {
      alpha_z_min[ 0 ] = ( (-mp_halfFOV[ 2 ]) - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_min[ 1 ] = ( mp_halfFOV[ 2 ] - p1_z[ 3 ] ) / r[ 3 ][ 2 ];
      alpha_z_max[ 0 ] = alpha_z_min[ 0 ];
      alpha_z_max[ 1 ] = alpha_z_min[ 1 ];
    }
    else
    {
      alpha_z_min[ 0 ] = 0.0;
      alpha_z_min[ 1 ] = 0.0;
      alpha_z_max[ 0 ] = 1.0;
      alpha_z_max[ 1 ] = 1.0;
    }

    if( r[ 4 ][ 2 ] != 0.0 ) // point4
    {
      alpha_z_min[ 2 ] = ( (-mp_halfFOV[ 2 ]) - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_min[ 3 ] = ( mp_halfFOV[ 2 ] - p1_z[ 4 ] ) / r[ 4 ][ 2 ];
      alpha_z_max[ 2 ] = alpha_z_min[ 2 ];
      alpha_z_max[ 3 ] = alpha_z_min[ 3 ];
    }
    else
    {
      alpha_z_min[ 2 ] = 0.0;
      alpha_z_min[ 3 ] = 0.0;
      alpha_z_max[ 2 ] = 1.0;
      alpha_z_max[ 3 ] = 1.0;
    }

    // Getting the minimum of alpha_z_min
    alpha_min = std::max( alpha_min, alpha_z_min[ std::minmax_element( alpha_z_min, alpha_z_min + 4 ).first - alpha_z_min ] );
    // Getting the maximum of alpha_z_max
    alpha_max = std::min( alpha_max, alpha_z_max[ std::minmax_element( alpha_z_max, alpha_z_max + 4 ).second - alpha_z_max ] );

    if( alpha_max <= alpha_min ) return 0;

    // temporary storage for TOF bins Gaussian integrals over the current voxel
    // allocation after potential returns
    HPFLTNB* tof_weights_temp = new HPFLTNB[tof_nb_bins];
    for (INTNB tof_bin=0; tof_bin<tof_nb_bins; tof_bin++) tof_weights_temp[tof_bin] = 0.;

    // Computing the first and the last plane
    int16_t j_min = 0, j_max = 0;
    if( r[ 0 ][ 1 ] > 0.0 )
    {
      j_min = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alpha_min * r[ 0 ][ 1 ] - p1_y[ 0 ] ) / mp_sizeVox[ 1 ] );
      j_max = ::floor( m_toleranceY + ( p1_y[ 0 ] + alpha_max * r[ 0 ][ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }
    else if( r[ 0 ][ 1 ] < 0.0 )
    {
      j_min = ::ceil( mp_nbVox[ 1 ] - m_toleranceY - ( mp_halfFOV[ 1 ] - alpha_max * r[ 0 ][ 1 ] - p1_y[ 0 ] ) / mp_sizeVox[ 1 ] );
      j_max = ::floor( m_toleranceY + ( p1_y[ 0 ] + alpha_min * r[ 0 ][ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[ 1 ] );
    }

    // Computing weight of normalization
    // Using the center
    HPFLTNB const factor( ::fabs( r[ 0 ][ 1 ] ) / ::sqrt( r2[ 0 ] + r2[ 1 ] + r2[ 2 ] ) );
    HPFLTNB const factor_for_tof( length_LOR / r[ 0 ][ 1 ] );
    HPFLTNB weight( mp_sizeVox[ 1 ] / factor );

    // Computing the increment for each 4 lines to project
    // (the center is not projected)
    for( uint8_t p = 0; p < 4; ++p )
    {
      ri[ p ][ 0 ] = r[ p + 1 ][ 0 ] / r[ p + 1 ][ 1 ];
      ri[ p ][ 1 ] = 1.0;
      ri[ p ][ 2 ] = r[ p + 1 ][ 2 ] / r[ p + 1 ][ 1 ];
    }

    // Increment orientation and Index orientation
    int8_t incr_orient_trans = 0;
    int8_t idx_orient_trans = 0;

    //int8_t const incr_orient_trans = p2_x[ 1 ] < p2_x[ 2 ] ? 1 : -1;
    int8_t const incr_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : -1;
    //int8_t const idx_orient_trans = p2_x[ 1 ] < p2_x[ 2 ] ? 1 : 0;
    int8_t const idx_orient_axial = p2_z[ 3 ] < p2_z[ 4 ] ? 1 : 0;

    // Computing an offset to get the correct position in plane
    HPFLTNB const offset = (-mp_halfFOV[ 1 ]) + ( mp_sizeVox[ 1 ] * 0.5 ) - p1_y[ 0 ];

    // Loop over the crossed Y planes
    for( int16_t j = j_min; j < j_max; ++j )
    {
      // Computing the coordinates of crossed plane
      HPFLTNB const step = offset + j * mp_sizeVox[ 1 ];
      // in Y (point 1 and 2)
      pos[ 0 ] = p1_x[ 1 ] + step * ri[ 0 ][ 0 ];
      pos[ 1 ] = p1_x[ 2 ] + step * ri[ 1 ][ 0 ];
      // in Z (point 3 and 4)
      pos[ 2 ] = p1_z[ 3 ] + step * ri[ 2 ][ 2 ];
      pos[ 3 ] = p1_z[ 4 ] + step * ri[ 3 ][ 2 ];

      // Computing the factor w1 and w2 normalizing the overlaps
      w1 = ::fabs( pos[ 0 ] - pos[ 1 ] );
      w2 = ::fabs( pos[ 2 ] - pos[ 3 ] );
      HPFLTNB final_weight = 0.0;
      if( ( w1 * w2 ) != 0 ) final_weight = weight / ( w1 * w2 );

      // Computing the index min and max on each axis
      // In X
      index[ 0 ] = ::floor( m_toleranceX + ( pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
      index[ 1 ] = ::floor( m_toleranceX + ( pos[ 1 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[ 0 ] );
      // In Z
      index[ 2 ] = ::floor( m_toleranceZ + ( pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );
      index[ 3 ] = ::floor( m_toleranceZ + ( pos[ 3 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[ 2 ] );

      int test_index1 = index[ 0 ] + j* mp_nbVox[ 0 ] + index[ 2 ] * m_nbVoxXY;
      int test_index2 = index[ 1 ] + j* mp_nbVox[ 0 ] + index[ 3 ] * m_nbVoxXY;
      bool test_index = (index[ 0 ] >= 0 && index[ 2 ] >= 0 && index[ 1 ] >= 0 && index[ 3 ] >= 0
                      && index[ 0 ] < mp_nbVox[ 0 ] && index[ 1 ] < mp_nbVox[ 0 ] && index[ 2 ] < mp_nbVox[ 2 ] && index[ 3 ] < mp_nbVox[ 2 ]);
      if (m_hasMask && test_index && !mp_mask[test_index1] && !mp_mask[test_index2]) continue;

      if( index[ 0 ] < index[ 1 ] )
      {
        incr_orient_trans = 1;
        idx_orient_trans = 1;
      }
      else if( index[ 0 ] > index[ 1 ] )
      {
        incr_orient_trans = -1;
        idx_orient_trans = 0;
      }

      // Getting the number of distance in X
      int16_t const n_distance_x = ::abs( index[ 0 ] - index[ 1 ] ) + 2;
      // Computing the distances in X
      distance_x[ 0 ] = pos[ 0 ];
      distance_x[ n_distance_x - 1 ] = pos[ 1 ];
      // Computing the rest if necessary
      if( n_distance_x > 2 )
      {
        for( int16_t d = 0; d < n_distance_x - 2; ++d )
          distance_x[ d + 1 ] = (-mp_halfFOV[ 0 ]) + ( index[ 0 ] + idx_orient_trans + ( incr_orient_trans * d ) ) * mp_sizeVox[ 0 ];
      }

      // Computing the final distance in X
      for( int16_t d = 0; d < n_distance_x - 1; ++d )
        distance_x[ d ] = ::fabs( distance_x[ d + 1 ] - distance_x[ d ] );

      // Getting the number of distance in Z
      int16_t const n_distance_z = ::abs( index[ 2 ] - index[ 3 ] ) + 2;
      // Storing the positions in Y
      distance_z[ 0 ] = pos[ 2 ];
      distance_z[ n_distance_z - 1 ] = pos[ 3 ];
      // Storing the rest if necessary
      if( n_distance_z > 2 )
      {
        for( int16_t d = 0; d < n_distance_z - 2; ++d )
          distance_z[ d + 1 ] = (-mp_halfFOV[ 2 ]) + ( index[ 2 ] + idx_orient_axial + ( incr_orient_axial * d ) ) * mp_sizeVox[ 2 ];
      }

      // Computing the final distance in Z
      for( int16_t d = 0; d < n_distance_z - 1; ++d )
        distance_z[ d ] = ::fabs( distance_z[ d + 1 ] - distance_z[ d ] );

      // Compute the first and the last relevant TOF bin for this voxel
      tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(trunc( (step * factor_for_tof - tof_half_span - length_LOR_half ) / tof_bin_size )));
      tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(trunc( (step * factor_for_tof + tof_half_span - length_LOR_half) / tof_bin_size )));
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
        // TOF coefficient = value of the Gaussian at voxel center (projected onto the LOR) * voxel size (projected onto the LOR)
        HPFLTNB temp = (step * factor_for_tof - lor_tof_center) / tof_sigma;
        tof_weight =  exp(- 0.5 * temp * temp ) * gaussian_norm_coef;
        // add the weight to the sum for normalization
        tof_norm_coef += tof_weight;
        // save the weight temporarily
        tof_weights_temp[tof_bin] = tof_weight;
        // update TOF center along the LOR for the next TOF bin
        lor_tof_center += tof_bin_size;
      }

      // compute and write the final TOF bin coefficients for current voxels
      if (tof_norm_coef>0.)
      {
        for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
        {
          // normalized tof coef
          tof_weights_temp[tof_bin] /= tof_norm_coef;
        }

        // Loop over the overlap and store the elements
        int16_t index_x, index_z;
        for( int16_t jj = 0; jj < n_distance_z - 1; ++jj )
        {
          index_z = index[ 2 ] + jj * incr_orient_axial;
          if( index_z < 0 || index_z > mp_nbVox[ 2 ] - 1 ) continue;

          for( int16_t ii = 0; ii < n_distance_x - 1; ++ii )
          {
            index_x = index[ 0 ] + ii * incr_orient_trans;
            if( index_x < 0 || index_x > mp_nbVox[ 0 ] - 1 ) continue;

            ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, index_x + j * mp_nbVox[ 0 ] + index_z * m_nbVoxXY, final_weight * distance_z[ jj ] * distance_x[ ii ], tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
          }
        }
      }
    }
    delete[] tof_weights_temp;
  }

  return 0;

}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
