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
  \brief    Implementation of class iProjectorIncrementalSiddonMulti
*/

#include "iProjectorIncrementalSiddonMulti.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorIncrementalSiddonMulti::iProjectorIncrementalSiddonMulti() : vProjector()
{
  m_nbLines = -1;
  // This projector is not compatible with SPECT attenuation correction because
  // the voxels contributing to the line are not strictly ordered with respect to
  // their distance to point 2 (due to the use of multiple lines that are
  // stack one after the other)
  m_compatibleWithSPECTAttenuationCorrection = false;
  // This projector is not compatible with compression as it works only with the
  // detection element indices
  m_compatibleWithCompression = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorIncrementalSiddonMulti::~iProjectorIncrementalSiddonMulti()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ReadConfigurationFile(const string& a_configurationFile)
{
  // Read the transaxial FWHM option
  string key_word = "number of lines";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_nbLines, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ReadOptionsList(const string& a_optionsList)
{
  // Read them
  if (ReadStringOption(a_optionsList, &m_nbLines, 1, ",", "Multi-Siddon configuration"))
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ReadConfigurationFile() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iProjectorIncrementalSiddonMulti::ShowHelpSpecific()
{
  cout << "This projector uses multiple ray-tracing for a single event in order to estimate the solid angle contribution." << endl;
  cout << "For each line of an event, the end points of the line are randomly chosen inside the detector element." << endl;
  cout << "The ray-tracing is performed with the incremental Siddon algorithm (see incrementalSiddon projector)." << endl;
  cout << "The only parameter of this projector is the number of lines to use per event:" << endl;
  cout << "  number of lines: the number of lines used per event" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::CheckSpecificParameters()
{
  // Send an error if less than 1 line
  if (m_nbLines<1)
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::CheckSpecificParameters() -> The provided number of lines is less than 1 !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=2) Cout("iProjectorIncrementalSiddonMulti::Initialize() -> Use incremental Siddon projector with " << m_nbLines << " lines per event" << endl);
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB iProjectorIncrementalSiddonMulti::EstimateMaxNumberOfVoxelsPerLine()
{
  // Find the maximum number of voxels along a given dimension
  INTNB max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxX();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxY()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxY();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxZ();
  // We should have at most 4 voxels in a given plane, so multiply by 4
  // (note: this is not true however it ensures no overflow and is already quite optimized for RAM usage !)
  max_nb_voxels_in_dimension *= 4;
  // Finally multiply by the number of lines
  max_nb_voxels_in_dimension *= ((INTNB)m_nbLines);
  // Return the value
  return max_nb_voxels_in_dimension;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ProjectWithoutTOF(int a_direction, oProjectionLine* ap_ProjectionLine )
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ProjectWithoutTOF() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorIncrementalSiddonMulti::Project without TOF -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Loop on the number of lines
  for (int line=0; line<m_nbLines; line++)
  {       
    // Get random positions from the scanner (mean depth of interaction intrinsicaly taken into account), taking POI into account if any
    if (mp_Scanner->GetRdmPositionsAndOrientations( ap_ProjectionLine->GetIndex1(), ap_ProjectionLine->GetIndex2(), 
                                                    ap_ProjectionLine->GetPosition1(), ap_ProjectionLine->GetPosition2(),
                                                    ap_ProjectionLine->GetOrientation1(), ap_ProjectionLine->GetOrientation2() ))
    {
      Cerr("***** vProjector::Project() -> A problem occured while getting positions and orientations from scanner !" << endl);
      return 1;
    }                                       
                                
    // Get end points position
    FLTNB* event1Float = ap_ProjectionLine->GetPosition1();
    FLTNB* event2Float = ap_ProjectionLine->GetPosition2();
    HPFLTNB event1[3] = { event1Float[0], event1Float[1], event1Float[2] };
    HPFLTNB event2[3] = { event2Float[0], event2Float[1], event2Float[2] };
  
  
  
  
  
  
  
  
    // **************************************
    // STEP 1: LOR length calculation
    // **************************************
    HPFLTNB length_LOR = ap_ProjectionLine->GetLength();
  
    // **************************************
    // STEP 2: Compute entrance voxel indexes
    // **************************************
    HPFLTNB alphaFirst[] = { 0.0, 0.0, 0.0 };
    HPFLTNB alphaLast[] = { 0.0, 0.0, 0.0 };
  
    HPFLTNB alphaMin = 0.0f, alphaMax = 1.0f;
    HPFLTNB delta_pos[] = {
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
    if( delta_pos[ 0 ] > 0.0f)
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMin * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMax * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    else if( delta_pos[ 0 ] < 0.0 )
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMax * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMin * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    if( delta_pos[ 0 ] == 0 )
    {
      iMin = 1, iMax = 0;
    }
  
    // For the Y-axis
    if( delta_pos[ 1 ] > 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMin * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMax * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] < 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMax * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMin * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] == 0 )
    {
      jMin = 1, jMax = 0;
    }
  
    // For the Z-axis
    if( delta_pos[ 2 ] > 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMin * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMax * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] < 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMax * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMin * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] == 0 )
    {
      kMin = 1, kMax = 0;
    }
  
    // Computing the last term n number of intersection
    INTNB n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
      + ( kMax - kMin + 1 );
  
    // Array storing the first alpha in X, Y and Z
    HPFLTNB alpha_XYZ[ 3 ] = { 1.0, 1.0, 1.0 };
  
    // Computing the first alpha in X
    if( delta_pos[ 0 ] > 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMin - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    else if( delta_pos[ 0 ] < 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMax - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
  
    // Computing the first alpha in Y
    if( delta_pos[ 1 ] > 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMin - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    else if( delta_pos[ 1 ] < 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMax - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
  
    // Computing the first alpha in Z
    if( delta_pos[ 2 ] > 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMin - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    else if( delta_pos[ 2 ] < 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMax - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
  
    // Computing the alpha updating
    HPFLTNB alpha_u[ 3 ] = {
      mp_sizeVox[0] / std::fabs( delta_pos[ 0 ] ),
      mp_sizeVox[1] / std::fabs( delta_pos[ 1 ] ),
      mp_sizeVox[2] / std::fabs( delta_pos[ 2 ] )
    };
  
    // Computing the index updating 
    INTNB index_u[ 3 ] = {
      (delta_pos[ 0 ] < 0) ? -1 : 1,
      (delta_pos[ 1 ] < 0) ? -1 : 1,
      (delta_pos[ 2 ] < 0) ? -1 : 1
    };
  
    // Check which alpha is the min/max and increment
    if( alpha_XYZ[ 0 ] == alphaMin ) alpha_XYZ[ 0 ] += alpha_u[ 0 ];
    if( alpha_XYZ[ 1 ] == alphaMin ) alpha_XYZ[ 1 ] += alpha_u[ 1 ];
    if( alpha_XYZ[ 2 ] == alphaMin ) alpha_XYZ[ 2 ] += alpha_u[ 2 ];
  
    // Computing the minimum value in the alpha_XYZ buffer
    HPFLTNB const min_alpha_XYZ = (std::min)( alpha_XYZ[ 0 ],
      (std::min)( alpha_XYZ[ 1 ], alpha_XYZ[ 2 ] ) );
  
    // Computing the first index of intersection
    HPFLTNB const alphaMid = ( min_alpha_XYZ + alphaMin ) / 2.0f;
    INTNB index_ijk[ 3 ] = {
      1 + (int)( ( event1[ 0 ] + alphaMid * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] ),
      1 + (int)( ( event1[ 1 ] + alphaMid * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] ),
      1 + (int)( ( event1[ 2 ] + alphaMid * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] )
    };
  
    INTNB const w = mp_nbVox[ 0 ];
    INTNB const wh = w * mp_nbVox[ 1 ];
  
    // Loop over the number of plane to cross
    HPFLTNB alpha_c = alphaMin;
    FLTNB coeff = 0.0;
    INTNB numVox = 0;
    for( int nP = 0; nP < n - 1; ++nP )
    {
      if( ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 1 ] )
       && ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 0 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 0 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff/((FLTNB)m_nbLines));
        }
  
        // Increment values
        alpha_c = alpha_XYZ[ 0 ];
        alpha_XYZ[ 0 ] += alpha_u[ 0 ];
        index_ijk[ 0 ] += index_u[ 0 ];
      }
      else if( ( alpha_XYZ[ 1 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 1 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 1 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 1 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff/((FLTNB)m_nbLines));
        }
  
        // Increment values
        alpha_c = alpha_XYZ[ 1 ];
        alpha_XYZ[ 1 ] += alpha_u[ 1 ];
        index_ijk[ 1 ] += index_u[ 1 ];
      }
      else if( ( alpha_XYZ[ 2 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 2 ] < alpha_XYZ[ 1 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 2 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 2 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff/((FLTNB)m_nbLines));
        }
  
        // Increment values
        alpha_c = alpha_XYZ[ 2 ];
        alpha_XYZ[ 2 ] += alpha_u[ 2 ];
        index_ijk[ 2 ] += index_u[ 2 ];
      }
    }                                      
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ProjectWithTOFPos(int a_Projector, oProjectionLine* ap_ProjectionLine)
{
  Cerr("***** iProjectorIncrementalSiddonMulti::ProjectWithTOFPos() -> Not yet implemented !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ProjectWithTOFBin(int a_Projector, oProjectionLine* ap_ProjectionLine)
{
  Cerr("***** iProjectorIncrementalSiddonMulti::ProjectWithTOFBin() -> Not yet implemented !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
