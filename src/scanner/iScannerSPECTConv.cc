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
  \ingroup  scanner
  \brief    Implementation of class iScannerSPECTConv
*/

#include "iScannerSPECTConv.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iScannerSPECTConv::iScannerSPECTConv() : vScanner() 
{
  // Set member variables to default values
  m_scannerType = SCANNER_SPECT_CONVERGENT;
  m_nbCrystals = -1;
  m_nbHeads = -1;
  m_acquisitionZoom = 1.;
  m_nbPixelsTrans = 0;
  m_pixelsSizeTrans = -1.;
  m_gapSizeTrans = -1.;
  m_nbPixelsAxial = 0;
  m_pixelsSizeAxial = -1.;
  m_gapSizeAxial = -1.;
  m_vPixelsSizeTrans = -1.;
  m_vPixelsSizeAxial = -1.;
  m_vNbPixelsTrans = 0;
  m_vNbPixelsAxial = 0;
  mp_nbOfBins[0] = 0;
  mp_nbOfBins[1] = 0;
  m_crystalDepth = -1.;
  mp_focalModelTrans = NULL;
  mp_nbCoefModelTrans = NULL;
  m2p_transFocalParameters = NULL;    
  mp_focalModelAxial = NULL;
  mp_nbCoefModelAxial = NULL;
  m2p_axialFocalParameters = NULL;  
  m_nbOfProjections = 0;
  mp_projectionAngles = NULL;
  mp_radius = NULL;
  mp_CORtoDetectorDistance = NULL; 
  m_rotDirection = GEO_ROT_CW;
  mp_crystalCentralPositionX = NULL; 
  mp_crystalCentralPositionY = NULL; 
  mp_crystalCentralPositionZ = NULL; 
  mp_crystalOrientationX = NULL; 
  mp_crystalOrientationY = NULL; 
  mp_crystalOrientationZ = NULL; 
  mp_crystalFocalPositionX = NULL; 
  mp_crystalFocalPositionY = NULL; 
  mp_crystalFocalPositionZ = NULL; 
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iScannerSPECTConv::~iScannerSPECTConv() 
{
  if (mp_focalModelTrans) delete[] mp_focalModelTrans;
  if (mp_nbCoefModelTrans) delete[] mp_nbCoefModelTrans;

  if (m2p_transFocalParameters)
  {
    for ( int i = 0; i < m_nbHeads; ++i )
      delete m2p_transFocalParameters[ i ];
    delete[] m2p_transFocalParameters;
  }

  if (mp_focalModelAxial) delete[] mp_focalModelAxial;
  if (mp_nbCoefModelAxial) delete[] mp_nbCoefModelAxial;

  if (m2p_axialFocalParameters)
  {
    for ( int i = 0; i < m_nbHeads; ++i )
      delete m2p_axialFocalParameters[ i ];
    delete[] m2p_axialFocalParameters;
  }
  
  if (mp_projectionAngles) delete[] mp_projectionAngles;
  if (mp_radius) delete[] mp_radius;
  if (mp_CORtoDetectorDistance) delete[] mp_CORtoDetectorDistance; 

  if (mp_crystalCentralPositionX) delete mp_crystalCentralPositionX;
  if (mp_crystalCentralPositionY) delete mp_crystalCentralPositionY;
  if (mp_crystalCentralPositionZ) delete mp_crystalCentralPositionZ;

  if (mp_crystalOrientationX) delete mp_crystalOrientationX;
  if (mp_crystalOrientationY) delete mp_crystalOrientationY;
  if (mp_crystalOrientationZ) delete mp_crystalOrientationZ;
 
  if (mp_crystalFocalPositionX) delete mp_crystalFocalPositionX;
  if (mp_crystalFocalPositionY) delete mp_crystalFocalPositionY;
  if (mp_crystalFocalPositionZ) delete mp_crystalFocalPositionZ;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      iScannerSPECTConv::DescribeSpecific()
  \brief   Implementation of the pure virtual eponym function that simply prints info about the scanner
*/
void iScannerSPECTConv::DescribeSpecific()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose==0) return;
  
  // Describe the scanner
  Cout("iScannerCT::DescribeSpecific() -> Here is some specific content of the SPECT convergent camera" << endl);
  
  
  Cout("  --> Total number of heads: " << m_nbHeads << endl);
  Cout("  --> Total number of crystals/pixels: " << m_nbCrystals << endl);
  Cout("  --> Total number of projections: " << m_nbOfProjections << endl);
  Cout("  --> Total number of transaxial bins: " << mp_nbOfBins[0] << endl);
  Cout("  --> Total number of axial bins: " << mp_nbOfBins[1] << endl);


  if( mp_projectionAngles ) 
  {
    int pr=0;
    int32_t nb_projections_by_head = m_nbOfProjections/m_nbHeads;

    Cout("  --> Projection angles (degree): "<< endl);
      
    for (int h=0 ; h<m_nbHeads ; h++)
    {
      Cout("  Head #"<<h<<": "<< endl);
      
      // Display ten by ten
      if( (nb_projections_by_head) > 10)
      {
        while ( pr< (h+1)*nb_projections_by_head-1
             && pr< ( m_nbOfProjections-1 ) )
        {
          for (int p=0 ; p<10 ; p++)
          {
            if(pr == (h+1)*nb_projections_by_head
            || pr == m_nbOfProjections)
              break;
            Cout(mp_projectionAngles[pr] << ", ");
            pr++;
          }
          Cout(endl);
        }
      }
      else
      {
        for (int p = h*nb_projections_by_head ; p<(h+1)*nb_projections_by_head ; p++)
          Cout(mp_projectionAngles[p] << ", ");
      }
    }
    Cout(endl);
  }

  if( mp_CORtoDetectorDistance ) 
  {
    int pr=0;
    int32_t nb_projections_by_head = m_nbOfProjections/m_nbHeads;
    
    Cout("  --> Distance between the center of rotation and the detector surface: "<< endl);
    
    for (int h=0 ; h<m_nbHeads ; h++)
    {
      Cout("  Head #"<<h<<": "<< endl);
      
      // Display ten by ten
      if( (nb_projections_by_head) > 10)
      {
        while ( pr< (h+1)*nb_projections_by_head-1
             && pr< ( m_nbOfProjections-1 ) )
        {
          for (int p=0 ; p<10 ; p++)
          {
            if(pr == (h+1)*nb_projections_by_head
            || pr == m_nbOfProjections)
              break;
            Cout(mp_CORtoDetectorDistance[pr] << ", ");
            pr++;
          }
          Cout(endl);
        }
      }
      else
      {
        for (int p = h*nb_projections_by_head ; p<(h+1)*nb_projections_by_head ; p++)
          Cout(mp_CORtoDetectorDistance[p] << ", ");
      }
    }
    Cout(endl);
  }


  if( mp_radius ) 
  {
    // Display ten by ten
    Cout("  --> Default radius for each head: " << endl);
    for (int h=0 ; h<m_nbHeads ; h++)
      Cout("      For head #"<<h<<": " << mp_radius[h] << endl);
  }
  
  
  Cout("  --> Total number of transaxial pixels as defined in the system file: " << m_nbPixelsTrans << endl);
  Cout("  --> Size of transaxial pixels as defined in the system file: " << m_pixelsSizeTrans << endl);
  Cout("  --> Gap size between each transaxial pixe as defined in the system file: " << m_gapSizeTrans << endl);
  Cout("  --> Total number of axial pixels as defined in the system file: " << m_nbPixelsAxial << endl);
  Cout("  --> Size of axial pixels as defined in the system file: " << m_pixelsSizeAxial << endl);
  Cout("  --> Gap size between each axial pixel as defined in the system file: " << m_gapSizeAxial << endl);
  Cout("  --> Number of transaxial virtual pixels (pixels actually used in reconstruction): " << m_vNbPixelsTrans << endl);
  Cout("  --> Number of axial virtual pixels (pixels actually used in reconstruction): " << m_vNbPixelsAxial << endl);
  Cout("  --> Transaxial size of virtual pixels (pixels actually used in reconstruction): " << m_vPixelsSizeTrans << endl);
  Cout("  --> Axial size of virtual pixels (pixels actually used in reconstruction): " << m_vPixelsSizeAxial << endl);
  
  Cout("  --> Focal models parameters: "<< endl);
  
  for (int h=0 ; h<m_nbHeads ; h++)
  {
    Cout("      For head #"<<h<<":" << endl);
    if( mp_focalModelTrans ) Cout("      Transaxial focal model: " << mp_focalModelTrans[h] << endl);
    if( mp_nbCoefModelTrans ) 
    { 
      Cout("      Number of coefficients: " << ( uint16_t ) mp_nbCoefModelTrans[h] << endl);
      Cout("      Coefficients: ");
      for( int p=0 ; p<mp_nbCoefModelTrans[h] ; p++)
        if( m2p_transFocalParameters[h] ) Cout( m2p_transFocalParameters[h][p] << ", " << endl);
    }

    if( mp_focalModelAxial ) Cout("      Axial focal model: " << mp_focalModelAxial[h] << endl);
    if( mp_nbCoefModelAxial ) 
    { 
      Cout("      Number of coefficients: " << ( uint16_t ) mp_nbCoefModelAxial[h] << endl);
      Cout("      Coefficients: ");
      for( int p=0 ; p<mp_nbCoefModelAxial[h] ; p++)
        if( m2p_axialFocalParameters[h] ) Cout( m2p_axialFocalParameters[h][p] << ", " << endl);
    }
    Cout(endl);
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerSPECTConv::Instantiate(bool a_scannerFileIsLUT)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerSPECTConv::Instantiate() -> Create scanner structure and read parameters from configuration file"<< endl); 

  // Get scanner manager
  sScannerManager* p_scannerManager; 
  p_scannerManager = sScannerManager::GetInstance();  

  // Get the number of heads in the scanner
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of detector heads", &m_nbHeads, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the number of SPECT heads !" << endl);
    return 1;
  }

  // Allocating memories that depend on number of heads
  mp_radius = new FLTNB[m_nbHeads];
  mp_focalModelTrans = new string[ m_nbHeads ];
  mp_focalModelAxial = new string[ m_nbHeads ];
  mp_nbCoefModelTrans = new uint8_t[ m_nbHeads ];
  mp_nbCoefModelAxial = new uint8_t[ m_nbHeads ];
  m2p_transFocalParameters = new FLTNB*[ m_nbHeads ];
  m2p_axialFocalParameters = new FLTNB*[ m_nbHeads ];

  // Get default radius
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "scanner radius", mp_radius, m_nbHeads, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the number of SPECT heads !" << endl);
    return 1;
  }
  for (int i = 0; i < m_nbHeads; i++) if (mp_radius[i]<=0.)
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> Scanner radius <= 0. ? really ? :) ... " << endl);
    return 1;
  } 
  // Detector dimensions
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "trans number of pixels", &m_nbPixelsTrans, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the transaxial number of pixels !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "trans pixel size", &m_pixelsSizeTrans, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the transaxial pixel size !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "axial number of pixels", &m_nbPixelsAxial, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the axial number of pixels !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "axial pixel size", &m_pixelsSizeAxial, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the axial pixel size !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "detector depth", &m_crystalDepth, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read detector depth !" << endl);
    return 1;
  } 

  // Head strings to get information from scanner configuration file
  string *head_name = new string[ m_nbHeads + 1 ];
  for ( int i = 0; i < m_nbHeads; i++ )
  {
    ostringstream oss( ostringstream::out );
    oss << "head" << i+1;
    head_name[ i ] = oss.str();
  }
  head_name[ m_nbHeads ] = "eof";

  // Read informations about SPECT heads
  for ( int i = 0; i < m_nbHeads; i++ )
  {
    if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "trans focal model", &mp_focalModelTrans[ i ], 1, KEYWORD_MANDATORY, head_name[ i ], head_name[ i + 1 ]))
    {
      Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read transaxial focal model of head " << i << " !" << endl);
      return 1;
    }
    if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "trans number of coef model", &mp_nbCoefModelTrans[ i ], 1, KEYWORD_MANDATORY, head_name[ i ], head_name[ i + 1 ]))
    {
      Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read transaxial number of coef model of head " << i << " !" << endl);
      return 1;
    }
    if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "axial focal model", &mp_focalModelAxial[ i ], 1, KEYWORD_MANDATORY, head_name[ i ], head_name[ i + 1 ]))
    {
      Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read axial focal model of head " << i << " !" << endl);
      return 1;
    }
    if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "axial number of coef model", &mp_nbCoefModelAxial[ i ], 1, KEYWORD_MANDATORY, head_name[ i ], head_name[ i + 1 ]))
    {
      Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read axial number of coef model of head " << i << " !" << endl);
      return 1;
    }
    m2p_transFocalParameters[i] = new FLTNB[mp_nbCoefModelTrans[ i ]];
    m2p_axialFocalParameters[i] = new FLTNB[mp_nbCoefModelAxial[ i ]];
    if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "trans parameters", m2p_transFocalParameters[i], mp_nbCoefModelTrans[ i ], KEYWORD_MANDATORY, head_name[ i ], head_name[ i + 1 ]))
    {
      Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read transaxial model parameters of head " << i << " !" << endl);
      return 1;
    }
    if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "axial parameters", m2p_axialFocalParameters[i], mp_nbCoefModelAxial[ i ], KEYWORD_MANDATORY, head_name[ i ], head_name[ i + 1 ]))
    {
      Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read axial model parameters of head " << i << " !" << endl);
      return 1;
    }
  }

  // Bed displacement
  m_defaultBedDisplacementInMm = 0.;
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "multiple bed displacement", &m_defaultBedDisplacementInMm, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the multiple bed displacement in the scanner header file !" << endl);
    return 1;
  }

  // Delete heads strings
  delete[] head_name;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn BuildLUT
  \param a_scannerFileIsLUT : boolean indicating if the file describing 
                              the SPECT camera is a generic file (0) or custom Look-up-table (1)
  \brief Call the functions to generate the LUT or read the user-made LUT depending on the user choice
  \return 0 if success, positive value otherwise
*/
int iScannerSPECTConv::BuildLUT(bool a_scannerFileIsLUT)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerSPECTConv::BuildLUT -> Build LUT for scanner elements coordinates and orientations"<< endl);
  
  // Initialize nb_crystals
  if (mp_nbOfBins[0] == 0 && mp_nbOfBins[1] == 0) // Number of bins not initialized in the datafile, then we read the number of crystals from the geom file
  {
    if (m_nbPixelsTrans==0 || m_nbPixelsAxial==0)
    {
     Cerr("***** iScannerSPECTConv::BuildLUT() -> Transaxial or axial number of pixels in the geometric file should be >0 !" << endl);
     return 1;
    }
    m_nbCrystals = m_nbPixelsTrans*m_nbPixelsAxial;
  }
  else // Initialization with number of bins
  {
    m_nbCrystals = mp_nbOfBins[0]*mp_nbOfBins[1];
  }
  
  mp_crystalCentralPositionX = new FLTNB[m_nbOfProjections*m_nbCrystals];
  mp_crystalCentralPositionY = new FLTNB[m_nbOfProjections*m_nbCrystals];
  mp_crystalCentralPositionZ = new FLTNB[m_nbOfProjections*m_nbCrystals];
  mp_crystalOrientationX = new FLTNB[m_nbOfProjections*m_nbCrystals];
  mp_crystalOrientationY = new FLTNB[m_nbOfProjections*m_nbCrystals];
  mp_crystalOrientationZ = new FLTNB[m_nbOfProjections*m_nbCrystals];
  mp_crystalFocalPositionX = new FLTNB[m_nbOfProjections*m_nbCrystals];
  mp_crystalFocalPositionY = new FLTNB[m_nbOfProjections*m_nbCrystals];
  mp_crystalFocalPositionZ = new FLTNB[m_nbOfProjections*m_nbCrystals]; 
    
  // Either generate the LUT from a generic file, or load the user precomputed LUT
  if (!a_scannerFileIsLUT)
  {
    if (ComputeLUT())
    {
     Cerr("***** iScannerSPECTConv::BuildLUT() -> A problem occured while generating scanner LUT !" << endl);
     return 1;
    }
  }
  else 
  {
    if (LoadLUT())
    {
      Cerr("***** iScannerSPECTConv::BuildLUT() -> A problem occured while loading scanner LUT !" << endl);
      return 1;
    }
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckParameters
  \brief Check that all parameters have been correctly initialized.
  \return 0 if success, positive value otherwise
  \todo Keep the check on crystal depth ?
*/
int iScannerSPECTConv::CheckParameters()
{  
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  
  // Check if all parameters have been correctly initialized. Return error otherwise
  if (m_nbCrystals<0)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Number of crystals has not been initialized !" <<endl);
    return 1;
  }
  if (m_nbHeads<0)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Number of heads in the SPECT system has not been initialized !" <<endl);
    return 1;
  }
  if (m_nbOfProjections<=0)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Number of projection angles has not been initialized !" <<endl);
    return 1;
  }
  if (m_nbPixelsTrans<=0 || m_nbPixelsAxial<=0 )
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Number of transaxial/axial pixels have not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (m_pixelsSizeTrans<=0 || m_pixelsSizeAxial<=0 )
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Transaxial/axial pixel sizes have not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (m_vNbPixelsTrans<=0 || m_vNbPixelsAxial<=0 )
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Transaxial/axial number of virtual pixels have not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (m_vPixelsSizeTrans<=0 || m_vPixelsSizeAxial<=0 )
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Transaxial/axial virtual pixel sizes have not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (m_gapSizeTrans<0 || m_gapSizeAxial<0 )
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Transaxial/axial pixel gap sizes have not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  // todo : maybe delete this as we don't used it
  if (m_crystalDepth<=0)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Crystal depth has not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (mp_focalModelTrans == NULL || mp_focalModelAxial == NULL)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Transaxial/axial focal models have not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_nbCoefModelTrans == NULL || mp_nbCoefModelAxial == NULL)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Number of coefficients for the transaxial/axial focal models have not correctly been initialized !" <<endl);
    return 1;
  }
  if (m2p_transFocalParameters == NULL || m2p_axialFocalParameters == NULL)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Parameters for the transaxial/axial focal models have not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_projectionAngles == NULL)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Projection angles have not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_radius == NULL)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Default scanner radius for each detector heads (scanner file) not correctly initialized !" <<endl);
    return 1;
  }
  if (mp_CORtoDetectorDistance == NULL)
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> Distances between center of rotation and each detector heads have not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_crystalCentralPositionX == NULL ||
      mp_crystalCentralPositionY == NULL ||
      mp_crystalCentralPositionZ == NULL  )
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> LUT elements (crystal central positions) have not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_crystalOrientationX == NULL ||
      mp_crystalOrientationY == NULL ||
      mp_crystalOrientationZ == NULL  )
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> LUT elements (crystal orientations) have not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_crystalFocalPositionX == NULL ||
      mp_crystalFocalPositionY == NULL ||
      mp_crystalFocalPositionZ == NULL  )
  {
    Cerr("***** iScannerSPECTConv::CheckParameters()-> LUT elements (crystal focal positions) have not correctly been initialized !" <<endl);
    return 1;
  }
  // No need to check mp_nbOfBins, Default values = 0

  // End
  m_allParametersChecked = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Initialize
  \brief Check general initialization and set several parameters to their default value
  \return 0 if success, positive value otherwise
*/
int iScannerSPECTConv::Initialize()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerSPECTConv::Initialize() -> Initialize remaining stuff for scanner to be ready"<< endl); 

  // Parameters checked ?
  if (!m_allParametersChecked)
  {
    Cerr("***** iScannerSPECTConv::Initialize() -> Parameters have not been checked !" << endl);
    return 1;
  }
  
  // Any initialization
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LoadLUT
  \brief Load a precomputed scanner LUT
  \details Read mandatory data from the header of the LUT. Then load the LUT elements for each crystal
  \todo Not yet implemented
  \return 0 if success, positive value otherwise
*/
int iScannerSPECTConv::LoadLUT()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // This has not been designed yet
  Cerr("iScannerSPECTConv::LoadLUT() -> Not yet implemented !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ComputeLUT
  \brief Computes the LUT of the scanner from a generic (.geom) file.
  \details Read mandatory data from the geom file. Then compute the LUT elements for each crystal from the geometry described in the file
           Compute the look-up-tables of the system containing the locations of the scanner elements center in space and their orientations
  \todo  center of rotation & head first angles : get this from acquisition header file
  \return 0 if success, positive value otherwise
*/
int iScannerSPECTConv::ComputeLUT()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if(m_verbose>=VERBOSE_DETAIL) Cout("iScannerSPECTConv::ComputeLUT() -> Start LUT generation" << endl);
  
  // Compute LUT according to SPECT camera geom files. Check the Castor presentation file for more details about the related  mandatory/optional fields to be read

  // ============================================================================================================
  // Parameters declarations
  // ============================================================================================================
  
  // todo center of rotation & head first angles : get this from acquisition header file
  FLTNB *CORx, *CORy, *CORz;
  FLTNB *head_angleX, *head_angleY, *head_angleZ;
  
  CORx = new FLTNB[m_nbHeads];
  CORy = new FLTNB[m_nbHeads]; 
  CORz = new FLTNB[m_nbHeads];
  head_angleX = new FLTNB[m_nbHeads];
  head_angleY = new FLTNB[m_nbHeads];
  head_angleZ = new FLTNB[m_nbHeads];
    
  //  for(int hId=0 ; hId<m_nbHeads ; hId++)
    // todo Initialisation COR and head_angle

  // ============================================================================================================
  // Generate the LUT
  // ============================================================================================================
  
  // oMatrix crystal_center_ref will be used to gather cartesian positions of the crystals center (on the surface of the crystals) 
  // in the reference head (directly above isocenter)
  oMatrix** crystal_center_ref = new oMatrix *[m_nbCrystals];
  // oMatrix crystal_center_out will be used to recover positions of each crystal after each rotation (each projection angles)
  oMatrix* crystal_center_out = new oMatrix(3,1);
  
  // same method for matrices recovering the positions of focal point
  oMatrix* focal_projection_position_mtx = new oMatrix(3,3);
  oMatrix* focal_projection_position_mtx_output = new oMatrix(3,3);
  
  for (int c = 0; c<m_nbCrystals; c++)
    crystal_center_ref[c]  = new oMatrix(3,1);


  for(int i=0 ; i<3 ; i++)
  {
    crystal_center_out->SetMatriceElt(i,0,0);
    
    for(int j=0 ; j<3 ; j++)
    {
      focal_projection_position_mtx->SetMatriceElt(i,j,0);
      focal_projection_position_mtx_output->SetMatriceElt(i,j,0);
    }
  }

  // Generation of the rotation matrix allowing to compute the position of all the rsectors. 
  oMatrix** rotation_mtx = new oMatrix*[m_nbOfProjections];

  for(int i=0; i<m_nbOfProjections; i++)
  {
    rotation_mtx[i] = new oMatrix(3,3);
    rotation_mtx[i]->SetMatriceElt(0,0, cos(mp_projectionAngles[i] * M_PI/180.) );
    rotation_mtx[i]->SetMatriceElt(1,0, sin(mp_projectionAngles[i] * M_PI/180.) );
    rotation_mtx[i]->SetMatriceElt(2,0,0);
    rotation_mtx[i]->SetMatriceElt(0,1, -sin(mp_projectionAngles[i] * M_PI/180.) );
    rotation_mtx[i]->SetMatriceElt(1,1, cos(mp_projectionAngles[i] * M_PI/180.) );
    rotation_mtx[i]->SetMatriceElt(2,1,0);
    rotation_mtx[i]->SetMatriceElt(0,2,0);
    rotation_mtx[i]->SetMatriceElt(1,2,0);
    rotation_mtx[i]->SetMatriceElt(2,2,1);
  }

  // Recover the trans/axial gap size from the geom file
  sScannerManager* p_scannerManager; 
  p_scannerManager = sScannerManager::GetInstance();
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "trans gap size", &m_gapSizeTrans, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::ComputeLUT() -> An error occurred while trying to read the transaxial gap size !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "axial gap size", &m_gapSizeAxial, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::ComputeLUT() -> An error occurred while trying to read the axial gap size !" << endl);
    return 1;
  }
  // Check for negative gap sizes
  if (m_gapSizeTrans<0. || m_gapSizeAxial<0.)
  {
    Cerr("***** iScannerSPECTConv::ComputeLUT() -> Crystal gap sizes cannot be negative !" << endl);
    return 1;
  }

  // Number of bins have been recovered from the datafile (>0)
  if (mp_nbOfBins[0]>0 && mp_nbOfBins[1]>0)
  {
    // In this case, we suppose that we have only one monolithic crystal, otherwise, we throw an error
    if (m_nbPixelsTrans>1 || m_nbPixelsAxial>1)
    {
      Cerr("***** iScannerSPECTConv::ComputeLUT() -> Number of bins provided in the datafile header while the scanner is not a monolithic crystal !" << endl);
      return 1;
    }
    // Check also that positive gap sizes is not possible
    if (m_gapSizeTrans>0. || m_gapSizeAxial>0.)
    {
      Cerr("***** iScannerSPECTConv::ComputeLUT() -> Positive crystal gap sizes has no sense with monolithic crystals !" << endl);
      return 1;
    }
    // Check that the zoom is not less than one
    if (m_acquisitionZoom<1.)
    {
      Cerr("***** iScannerSPECTConv::ComputeLUT() -> An acquisition zoom less than 1 has no sense !" << endl);
      return 1;
    }
    // Compute the size of the pseudo crystals from the number of bins and the size of the monolithic crystal
    m_vPixelsSizeTrans = m_pixelsSizeTrans/(((FLTNB)(mp_nbOfBins[0]))*m_acquisitionZoom);
    m_vPixelsSizeAxial = m_pixelsSizeAxial/(((FLTNB)(mp_nbOfBins[1]))*m_acquisitionZoom);
    // Set the number of crystals to be used in the reconstruction
    m_vNbPixelsTrans = mp_nbOfBins[0];
    m_vNbPixelsAxial = mp_nbOfBins[1];
  }
  // Otherwise, we recover these informations from the geom file, which means that we have an already discretized detector
  else
  {
    // Check that no zoom has been provided, because it makes sense only with monolithic detectors
    if (m_acquisitionZoom!=1.)
    {
      Cerr("***** iScannerSPECTConv::ComputeLUT() -> An acquisition zoom has no sense with pixelated detectors !" << endl);
      return 1;
    }
    m_vPixelsSizeTrans = m_pixelsSizeTrans;
    m_vPixelsSizeAxial = m_pixelsSizeAxial;
    m_vNbPixelsTrans = m_nbPixelsTrans;
    m_vNbPixelsAxial = m_nbPixelsAxial;
  }
  
  // Generate cartesian coordinates of the crystal centers for the SPECT at the reference position (directly above isocenter)
  // Starting position of the SPECT camera for all axis
  FLTNB x_start = -(((FLTNB)(m_vNbPixelsTrans))*m_vPixelsSizeTrans + ((FLTNB)(m_vNbPixelsTrans-1))*m_gapSizeTrans) / 2. + (m_vPixelsSizeTrans/2.);
  FLTNB z_start = -(((FLTNB)(m_vNbPixelsAxial))*m_vPixelsSizeAxial + ((FLTNB)(m_vNbPixelsAxial-1))*m_gapSizeAxial) / 2. + (m_vPixelsSizeAxial/2.);
  // Loop on crystals
  for ( uint32_t jj = 0; jj < m_vNbPixelsAxial; ++jj )
  {
    FLTNB Zcrist = z_start + jj * (m_vPixelsSizeAxial + m_gapSizeAxial);
    for ( uint32_t ii = 0; ii < m_vNbPixelsTrans; ++ ii )
    {
      FLTNB Xcrist = x_start + ii * (m_vPixelsSizeTrans + m_gapSizeTrans);
      crystal_center_ref[ii + m_vNbPixelsTrans * jj]->SetMatriceElt(0,0,Xcrist);
      crystal_center_ref[ii + m_vNbPixelsTrans * jj]->SetMatriceElt(1,0,     0); // This value will be set for each projection angles
      crystal_center_ref[ii + m_vNbPixelsTrans * jj]->SetMatriceElt(2,0,Zcrist);
    }
  }

  // ============================================================================================================
  // Loop over all the projection angles
  // Positions of the scanner elements are progressively stored in the LUT file.
  // ============================================================================================================
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Generate positions for each projection angle"<< endl); 
  // Loop on projections
  for (int a=0 ; a<m_nbOfProjections ; a++)
  {
    // Set y-position of the crystal reference matrix according to the distance between SPECT camera and center of rotation of the specific head.(mp_CORtoDetectorDistance)
    int hId = a*m_nbHeads/m_nbOfProjections;
    // Get surface crystal position (and not center crystal positions for SPECT), and set the reference crystal matrix
    FLTNB Ycrist = mp_CORtoDetectorDistance[hId];
    // Loop on crystals
    for (int c=0 ; c < m_nbCrystals ; c++)
    {
      crystal_center_ref[c]->SetMatriceElt(1,0,Ycrist);
      // Crystal indexation
      int cryID = a*m_nbCrystals + c;
      // Set cartesian coordinates of the crystal surface center positions
      rotation_mtx[a]->Multiplication(crystal_center_ref[c], crystal_center_out);
      // Get crystal surface positions
      mp_crystalCentralPositionX[cryID] = crystal_center_out->GetMatriceElt(0,0);
      mp_crystalCentralPositionY[cryID] = crystal_center_out->GetMatriceElt(1,0);
      mp_crystalCentralPositionZ[cryID] = crystal_center_out->GetMatriceElt(2,0);
      // Compute focal positions
      if (ComputeFocalPositions(crystal_center_ref[c]->GetMatriceElt(0,0), 
                                crystal_center_ref[c]->GetMatriceElt(1,0), 
                                crystal_center_ref[c]->GetMatriceElt(2,0), 
                                                                      hId, 
                                                                    cryID) )
      {
        Cerr("***** iScannerSPECTConv::ComputeLUT() -> An error occurred while computing the focal positions ! " << endl);
        return 1;
      }
      // Set elements of the focal matrix
      focal_projection_position_mtx->SetMatriceElt(0,0,mp_crystalFocalPositionX[cryID] );
      focal_projection_position_mtx->SetMatriceElt(1,0,mp_crystalFocalPositionY[cryID] );
      focal_projection_position_mtx->SetMatriceElt(2,0,mp_crystalFocalPositionZ[cryID] );
      focal_projection_position_mtx->SetMatriceElt(0,1,0);
      focal_projection_position_mtx->SetMatriceElt(1,1,0);
      focal_projection_position_mtx->SetMatriceElt(2,1,0);
      focal_projection_position_mtx->SetMatriceElt(0,2,0);
      focal_projection_position_mtx->SetMatriceElt(1,2,0);
      focal_projection_position_mtx->SetMatriceElt(2,2,1);
      // Rotate and compute projections of focal position for this crystal
      rotation_mtx[a]->Multiplication(focal_projection_position_mtx, focal_projection_position_mtx_output);
      mp_crystalFocalPositionX[cryID] = focal_projection_position_mtx_output->GetMatriceElt(0,0);
      mp_crystalFocalPositionY[cryID] = focal_projection_position_mtx_output->GetMatriceElt(1,0); 
      mp_crystalFocalPositionZ[cryID] = focal_projection_position_mtx_output->GetMatriceElt(2,0); 
      // Get crystal orientations
      mp_crystalOrientationX[cryID] = rotation_mtx[a]->GetMatriceElt(0,0);
      mp_crystalOrientationY[cryID] = rotation_mtx[a]->GetMatriceElt(0,1);
      mp_crystalOrientationZ[cryID] = 0;
    }
  }

  // Delete all temporary tables
  for (int i=0; i<m_nbOfProjections; i++) delete rotation_mtx[i];
  delete[] rotation_mtx;
  for (int c=0; c<m_nbCrystals; c++) delete crystal_center_ref[c];
  delete[] crystal_center_ref;
  delete crystal_center_out;
  delete focal_projection_position_mtx;
  delete focal_projection_position_mtx_output;
  delete[] head_angleX;
  delete[] head_angleY;
  delete[] head_angleZ;
  delete[] CORx;
  delete[] CORy;
  delete[] CORz;

  // End
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> LUT generation completed" << endl);
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ComputeFocalPositions
  \param a_posX : cartesian position of the crystal on the x-axis 
  \param a_posY : cartesian position of the crystal on the y-axis 
  \param a_posZ : cartesian position of the crystal on the z-axis 
  \param a_headID : head index of the crystal 
                   (required to select the related focal model parameters
                    and distance to center of rotation)
  \param a_cryID : crystal index in the LUT
  \brief Compute focal positions for a specific crystal ID
  \details Compute the focal positions using the implemented "constant", 
           "polynomial", "slanthole" focal models related to the gamma cameras
           The "custom" focal model is dedicated to user-made focal model
           and should be implemented by the user
  \return 0 if success, positive value otherwise
*/
int iScannerSPECTConv::ComputeFocalPositions(FLTNB a_posX, FLTNB a_posY, FLTNB a_posZ, int a_headID, int a_cryID)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_NORMAL)

  mp_crystalFocalPositionY[a_cryID] = -a_posY;

  // ===================================================================
  // AXIAL FOCAL POSITIONS
  // ===================================================================
  // constant
  if(mp_focalModelAxial[a_headID] == "constant")
  {
    mp_crystalFocalPositionZ[a_cryID] = a_posZ;
  }
  
  // polynomial
  else if(mp_focalModelAxial[a_headID] == "polynomial" ) 
  {
    if (mp_nbCoefModelAxial[a_headID] >= 1 && mp_nbCoefModelAxial[a_headID] <= 3)
    {
      FLTNB focal_posZ = 0;
      
      for(int coef=0 ; coef<mp_nbCoefModelAxial[a_headID] ; coef++)
        focal_posZ += m2p_axialFocalParameters[a_headID][coef]*abs(pow(a_posZ,coef));
        
      // Compute the projection of the focal position on the Y plane (thales) in order to be sure that the focal positions are not in the field of view (LOR would not go through the whole FOV) 
      mp_crystalFocalPositionZ[a_cryID] = -(2*mp_CORtoDetectorDistance[a_headID]-focal_posZ)*a_posZ/focal_posZ;
      
    }
    else if (mp_nbCoefModelAxial[a_headID] > 3) // Error (nb parameters >3)
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Error, the number of coeffs for the axial model should be <4 (max : polynom order 2 = 3 coeffs) ! " << endl);
      return 1;
    }
    else // Error (nb parameters == 0)
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Error, the number of coeffs for the axial model should be >0 ('axial number of coef model' in geom file) ! " << endl);
      return 1;
    }
  }

  // slanthole
  else if(mp_focalModelAxial[a_headID] == "slanthole" )
  {
    if (mp_nbCoefModelAxial[a_headID] == 1) 
    {
      // Compute the projection of the focal position on the Y plane
      mp_crystalFocalPositionZ[a_cryID] = -2.*m2p_axialFocalParameters[a_headID][0]*mp_CORtoDetectorDistance[a_headID] + a_posZ;
    } 
    else // Error (nb parameters != 1)
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Error, the number of coeffs for a slanthole model should be equal to 1 (slope) : 'axial number of coef model' in geom file ! " << endl);
      return 1;
    }
  }
  
  // custom
  else if(mp_focalModelAxial[a_headID] == "custom" )
  {
    if (mp_nbCoefModelAxial[a_headID] == 1) // Uni-Parameter
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Custom model should be implemented by the user (in iScannerSPECTConv::ComputeFocalPositions()) ! " << endl);
      return 1;
    } 
    else if (mp_nbCoefModelAxial[a_headID] > 1) // Multi-Parameters
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Custom model should be implemented by the user (in iScannerSPECTConv::ComputeFocalPositions()) ! " << endl);
      return 1;
    }
    else // Error (nb parameters == 0)
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Custom model should be implemented by the user (in iScannerSPECTConv::ComputeFocalPositions()) ! " << endl);
      return 1;
    }
  }
  else  // Error, unknown model
  {
    Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Error, current model " << mp_focalModelAxial[a_headID] << " is unknown !" << endl);
    Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Should be either 'constant' (parallel), 'polynomial', 'hyperbolic, or 'custom'" << endl);
    return 1;
  }

  // ===================================================================
  // TRANSAXIAL FOCAL POSITIONS
  // ===================================================================
  // constant
  if(mp_focalModelTrans[a_headID] == "constant")
  {
    mp_crystalFocalPositionX[a_cryID] = a_posX;
  }
  
  // polynomial
  else if(mp_focalModelTrans[a_headID] == "polynomial" ) 
  {
    if (mp_nbCoefModelTrans[a_headID] >= 1 && mp_nbCoefModelTrans[a_headID] <= 3)
    {
      FLTNB focal_posX = 0;
      
      for(int coef=0 ; coef<mp_nbCoefModelTrans[a_headID] ; coef++)
        focal_posX += m2p_transFocalParameters[a_headID][coef]*abs(pow(a_posX,coef));
        
      // Compute the projection of the focal position on the Y plane (thales) in order to be sure that the focal positions are not in the field of view (LOR would not go through the whole FOV) 
      mp_crystalFocalPositionX[a_cryID] = -(2*mp_CORtoDetectorDistance[a_headID]-focal_posX)*a_posX/focal_posX;
      
    } 
    else if (mp_nbCoefModelTrans[a_headID] > 3)  // Error (nb parameters >3)
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Error, the number of coeffs for the trans model should be <4 (max : polynom order 2 = 3 coeffs) ! " << endl);
      return 1;
    }
    else // Error (nb parameters == 0)
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Error, the number of coeffs for the axial model should be >0 ('trans number of coef model' in geom file) ! " << endl);
      return 1;
    }
  }
  
  // slanthole
  else if(mp_focalModelTrans[a_headID] == "slanthole" ) 
  {
    if (mp_nbCoefModelTrans[a_headID] == 1) 
    {
      // Compute the projection of the focal position on the Y plane
      mp_crystalFocalPositionX[a_cryID] = -2.*m2p_transFocalParameters[a_headID][0]*mp_CORtoDetectorDistance[a_headID] + a_posX;
    } 
    else // Error (nb parameters != 1
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Error, the number of coeffs for a slanthole model should be equal to 1 (slope) : 'trans number of coef model' in geom file ! " << endl);
      return 1;
    }
  }
  
  // custom
  else if(mp_focalModelTrans[a_headID] == "custom" )
  {
    if (mp_nbCoefModelTrans[a_headID] == 1) // Uni-Parameter
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Custom model should be implemented by the user (in iScannerSPECTConv::ComputeFocalPositions()) ! " << endl);
      return 1;
    } 
    else if (mp_nbCoefModelTrans[a_headID] > 1) // Multi-Parameters
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Custom model should be implemented by the user (in iScannerSPECTConv::ComputeFocalPositions()) ! " << endl);
      return 1;
    }
    else // Error (nb parameters == 0)
    {
      Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Custom model should be implemented by the user (in iScannerSPECTConv::ComputeFocalPositions()) ! " << endl);
      return 1;
    }
  }
  else  // Error, unknown model
  {
    Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Error, current model " << mp_focalModelTrans[a_headID] << " is unknown !" << endl);
    Cerr("***** iScannerSPECTConv::ComputeFocalPositions() -> Should be either 'constant' (parallel), 'polynomial', 'hyperbolic, or 'custom'" << endl);
    return 1;
  }
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerSPECTConv::GetPositionsAndOrientations( int a_index1, int a_index2,
                                                    FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                                    FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3],
                                                    FLTNB* ap_POI1, FLTNB* ap_POI2 )
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // SPECT indexes : 1st for the projection, 2nd for the crystal

  // First check projection angle existency
  if (a_index1<0 || a_index1>=m_nbOfProjections)
  {
    Cerr("***** iScannerSPECTConv::GetPositionsAndOrientations() -> Projection index (" << a_index1 << ") out of range [0:" << m_nbOfProjections-1 << "] !" << endl);
    return 1;
  }
  
  // Second check crystals existency
  if (a_index2<0 || a_index2>=m_nbCrystals)
  {
    Cerr("***** iScannerSPECTConv::GetPositionsAndOrientations() -> Crystal index (" << a_index2 << ") out of range [0:" << m_nbCrystals-1 << "] !" << endl);
    return 1;
  }

  // Get crystal index related to the projection
  int index = a_index1*m_nbCrystals + a_index2;
  
  // Get position of the focal point related to the crystal
  ap_Position1[0] = mp_crystalFocalPositionX[index];
  ap_Position1[1] = mp_crystalFocalPositionY[index];
  ap_Position1[2] = mp_crystalFocalPositionZ[index];

  // todo : Due to the current implementation of SPECT projection, POI 
  //        and DOI are not handled and ignored.
  //        An error is returned if POI are provided

  // Case when POI is not provided
  if (ap_POI2==NULL)
  {
    //FLTNB depth = mp_meanDepthOfInteraction[GetLayer(a_index2)] - mp_sizeCrystalDepth[GetLayer(a_index2)]/2;
    //ap_Position2[0] = mp_crystalCentralPositionX[a_index2] + depth*mp_crystalOrientationX[a_index2];
    //ap_Position2[1] = mp_crystalCentralPositionY[a_index2] + depth*mp_crystalOrientationY[a_index2];
    //ap_Position2[2] = mp_crystalCentralPositionZ[a_index2] + depth*mp_crystalOrientationZ[a_index2];
    ap_Position2[0] = mp_crystalCentralPositionX[index];
    ap_Position2[1] = mp_crystalCentralPositionY[index];
    ap_Position2[2] = mp_crystalCentralPositionZ[index];
  }
  // Case when POI[2] is negative (meaning we only have POI[0] or POI[1] specified and to be taken into account)
  else if (ap_POI2[2]<0.)
  {
    Cerr("***** iScannerSPECTConv::GetPositionsAndOrientations() -> POI management not implemented yet for SPECT !" << endl);
    return 1;
  }
  // Case when only the DOI is provided
  else if (ap_POI2[0]==0. && ap_POI2[1]==0.)
  {
    Cerr("***** iScannerSPECTConv::GetPositionsAndOrientations() -> POI management not implemented yet for SPECT !" << endl);
    return 1;
  }
  // Case when the full POI is taken into account
  else
  {
    Cerr("***** iScannerSPECTConv::GetPositionsAndOrientations() -> POI management not implemented yet for SPECT !" << endl);
    return 1;
  }

  // Get orientations
  ap_Orientation1[0] = -mp_crystalOrientationX[a_index2];
  ap_Orientation1[1] = -mp_crystalOrientationY[a_index2];
  ap_Orientation1[2] = -mp_crystalOrientationZ[a_index2];
  ap_Orientation2[0] = mp_crystalOrientationX[a_index2];
  ap_Orientation2[1] = mp_crystalOrientationY[a_index2];
  ap_Orientation2[2] = mp_crystalOrientationZ[a_index2];

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerSPECTConv::GetRdmPositionsAndOrientations(int a_index1, int a_index2,
                                                      FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                                      FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3] )
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // SPECT indexes : 1st for the projection, 2nd for the crystal

  // First check projection angle existency
  if (a_index1<0 || a_index1>=m_nbOfProjections)
  {
    Cerr("***** iScannerSPECTConv::GetPositionsAndOrientations() -> Projection index (" << a_index1 << ") out of range [0:" << m_nbOfProjections-1 << "] !" << endl);
    return 1;
  }
  
  // Second check crystals existency
  if (a_index2<0 || a_index2>=m_nbCrystals)
  {
    Cerr("***** iScannerSPECTConv::GetPositionsAndOrientations() -> Crystal index (" << a_index2 << ") out of range [0:" << m_nbCrystals-1 << "] !" << endl);
    return 1;
  }
  
  // Get crystal index related to the projection
  int index = a_index1*m_nbCrystals + a_index2;

  // Get position of the focal point related to the crystal
  ap_Position1[0] = mp_crystalFocalPositionX[index];
  ap_Position1[1] = mp_crystalFocalPositionY[index];
  ap_Position1[2] = mp_crystalFocalPositionZ[index];
  
  // Get instance of random number generator
  sRandomNumberGenerator* p_RNG = sRandomNumberGenerator::GetInstance(); 

  // Get random numbers for the first crystal
  FLTNB axial = (p_RNG->GenerateRdmNber()-0.5) * m_pixelsSizeAxial;
  FLTNB trans = (p_RNG->GenerateRdmNber()-0.5) * m_pixelsSizeTrans;
  // Do not consider random depth (position on the surface of the crystal)
  //FLTNB depth = (p_RNG->GenerateRdmNber()-0.5) * m_crystalDepth;
  
  ap_Position2[0] = mp_crystalCentralPositionX[index] + trans*mp_crystalOrientationY[index] + axial*mp_crystalOrientationX[index]*mp_crystalOrientationZ[index];
  ap_Position2[1] = mp_crystalCentralPositionY[index] + trans*mp_crystalOrientationX[index] + axial*mp_crystalOrientationY[index]*mp_crystalOrientationZ[index];
  ap_Position2[2] = mp_crystalCentralPositionZ[index] + axial*sqrt(1-mp_crystalOrientationZ[index]*mp_crystalOrientationZ[index]);

  // Get orientations
  ap_Orientation1[0] = -1.;
  ap_Orientation1[1] = -1.;
  ap_Orientation1[2] = -1.;
  ap_Orientation2[0] = mp_crystalOrientationX[a_index2];
  ap_Orientation2[1] = mp_crystalOrientationY[a_index2];
  ap_Orientation2[2] = mp_crystalOrientationZ[a_index2];
  
  return 0;
}
  
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerSPECTConv::GetPositionWithRandomDepth(int a_index1, int a_index2, FLTNB ap_Position1[3], FLTNB ap_Position2[3])
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  // This function was first implemented for PET testing purpose. Not implemented yet.
  Cerr("***** iScannerSPECTConv::GetPositionWithRandomDepth() -> This function was implemented for PET testing purpose. Not implemented for SPECT !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerSPECTConv::GetTwoCorners(int a_index1, int a_index2,
                                     FLTNB ap_CornerInf1[3], FLTNB ap_CornerSup1[3],
                                     FLTNB ap_CornerInf2[3], FLTNB ap_CornerSup2[3])
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  Cerr("***** iScannerSPECTConv::GetTwoCorners() -> Not implemented yet !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerSPECTConv::GetEdgesCenterPositions( int a_index1, int a_index2,
                                                FLTNB ap_pos_line_point1[3], FLTNB ap_pos_line_point2[3],
                                                FLTNB ap_pos_point1_x[4], FLTNB ap_pos_point1_y[4], FLTNB ap_pos_point1_z[4],
                                                FLTNB ap_pos_point2_x[4], FLTNB ap_pos_point2_y[4], FLTNB ap_pos_point2_z[4]
)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  Cerr("***** iScannerSPECTConv::GetEdgesCenterPositions() -> Not implemented yet !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetGeometricInfoFromDataFile
  \param a_pathToDF : string containing the path to datafile header
  \brief Recover geometric informations specific to the scanner class from the datafile header
  \details -Recover nb of bins and projections
           -Recover the projection angles from the header 
           (directly read from the datafile, or extrapolated from a first and last angle)
           -Recover the distance between the gamma camera detector surfaces and the center of rotation
           (directly read from the datafile, or extracted from the gamma camera configuratino file)
  \return 0 if success, positive value otherwise
*/
int iScannerSPECTConv::GetGeometricInfoFromDataFile(string a_pathToDF)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_DETAIL) Cout("iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Specific to SPECT" << endl);

  // This function is intended to be called after the scanner initialization, so that any field present in the datafile header, similar to
  // one in the scanner configuration file, may overload the scanner value.

  // Get the number of bins for monolithic detectors
  if (ReadDataASCIIFile(a_pathToDF , "Number of bins", mp_nbOfBins, 2, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error while reading number of bins in the header data file " << endl);
    return 1;
  }
  // Get the acquisition zoom for monolithic detectors
  if (ReadDataASCIIFile(a_pathToDF , "Zoom", &m_acquisitionZoom, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error while reading acquisition zoom in the header data file " << endl);
    return 1;
  }
  // Get the number of projections
  if (ReadDataASCIIFile(a_pathToDF , "Number of projections", &m_nbOfProjections, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error while reading number of projections in the header data file " << endl);
    return 1;
  }
  // Get the head rotation direction
  string rotation_direction = "";
  if (ReadDataASCIIFile(a_pathToDF , "Head rotation direction", &rotation_direction, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error while reading head rotation orientation in the header data file " << endl);
    return 1;
  }
  // Set and check the rotation direction
  if (SetRotDirection(rotation_direction) )
  {
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error occured while trying to initialize head rotation orientation " << endl);
    return 1;
  }

  // Allocate a one-dimensional vector to retrieve angles from the datafile, then copy it in the member variables
  FLTNB* angles = new FLTNB[m_nbOfProjections];
  FLTNB first_and_last_angles[2] = {-1.,-1.};

  // Two possible initializations for projection angles :
  // - All angles are provided with the keyword "Angles"
  // - The first and last angles are provided, and the intermediate angles are extrapolated from the number of projections
  
  // 'First/Last angles' tags : checking issue during data reading/conversion (==1)
  if (ReadDataASCIIFile(a_pathToDF, "First and last projection angles", first_and_last_angles, 2, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error while reading Angle mandatory field in the header data file '" << endl);
    return 1;
  }
  
  // Check for 'Projection angles' tag
  int rvalue = ReadDataASCIIFile(a_pathToDF, "Projection angles", angles, m_nbOfProjections, KEYWORD_OPTIONAL);
  
  // Error while reading "Angles" tag (==1)
  if(rvalue==1)
  {
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error while reading Angles field in the header data file !'" << endl);
    return 1;
  }

  // Check if information on projection angles has been provided
  if ( rvalue>=2 &&  // "Angles" tag not found
      (first_and_last_angles[0] <0 || first_and_last_angles[1] <0) ) // Tags first/last angles not found)
  {
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> No information on projection angles provided in the datafile !'" << endl);
    Cerr("                                                           This information should be provided using either the 'Angles' tag, or both 'First angles', 'Last angles' tags !'" << endl);
    return 1;
  }
  else if (rvalue>=2) // "Angles" tag not found, but first angle and last angle provided
  {
    // Put first and last angles in [0:360[
    while (first_and_last_angles[0]>=360.) first_and_last_angles[0] -= 360.;
    while (first_and_last_angles[0]<0.) first_and_last_angles[0] += 360.;
    while (first_and_last_angles[1]>=360.) first_and_last_angles[1] -= 360.;
    while (first_and_last_angles[1]<0.) first_and_last_angles[1] += 360.;
    // Rotation direction
    FLTNB dir = (m_rotDirection == GEO_ROT_CCW) ? -1. : 1.;
    // Compute angle increment
    FLTNB angle_increment = dir*(first_and_last_angles[1] - first_and_last_angles[0]);
    while (angle_increment>=360.) angle_increment -= 360.;
    while (angle_increment<0.) angle_increment += 360.;
    angle_increment /= ((FLTNB)(m_nbOfProjections-1));
    // Compute all projection angles
    for(int a=0 ; a<m_nbOfProjections ; a++)
    {
      angles[a] = first_and_last_angles[0] + dir * angle_increment * ((FLTNB)a);
      while (angles[a]>=360.) angles[a] -= 360.;
      while (angles[a]<0.) angles[a] += 360.;
    }
    // Set the last angle to be sure it is the one requested (rounding errors may change it a bit)
    angles[m_nbOfProjections-1] = first_and_last_angles[1];
  }
  // else : Angles have been recovered using ReadDataASCIIFile() above

  // Instanciate here the projection angles variable
  mp_projectionAngles = new FLTNB[m_nbOfProjections];
    
  for(int a=0 ; a<m_nbOfProjections ; a++)
    mp_projectionAngles[a] = angles[a];

  // Recover distance between the detectors and scanner center of rotation
  // Allocate with the number of projections by default
  mp_CORtoDetectorDistance = new FLTNB[m_nbOfProjections]; 
  
  // Initialize by default with the scanner radius
  for(int hId=0 ; hId<m_nbHeads ; hId++)
    for(int a=0 ; a<m_nbOfProjections/m_nbHeads ; a++)
      mp_CORtoDetectorDistance[a+hId*m_nbOfProjections/m_nbHeads] = mp_radius[hId];
  
  // Read datafile value if any. First case : we have a radius specific to each projections
  int read_flag = 0;
  read_flag = ReadDataASCIIFile(a_pathToDF, "Distance camera surface to COR", mp_CORtoDetectorDistance, m_nbOfProjections, KEYWORD_OPTIONAL);
  
  if (read_flag==1)
  {
    // Error during reading
    Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error while reading the distance between the camera detectors to the center of rotation in the header data file " << endl);
    return 1;
  } 
  // Check the second case : we have a global radius for each projection
  else if (read_flag==2)
  {
    read_flag = ReadDataASCIIFile(a_pathToDF, "Global distance camera surface to COR", mp_CORtoDetectorDistance, 1, KEYWORD_OPTIONAL);

    if(read_flag==0) 
    {
      // Field was found : Initialize the distance for each projection angle with the global value
      for(int a=1 ; a<m_nbOfProjections ; a++)
        mp_CORtoDetectorDistance[a] = mp_CORtoDetectorDistance[0];
    }
    else if(read_flag==1)
    {
      // Error during reading
      Cerr("***** iScannerSPECTConv::GetGeometricInfoFromDataFile() -> Error while reading the global distance between the camera detectors to the center of rotation in the header data file " << endl);
      return 1;
    }
    else if (read_flag==2)
    {
      // Initialization with default values from the scanner file
      for(int hId=0 ; hId<m_nbHeads ; hId++)
        for(int a=0 ; a<m_nbOfProjections/m_nbHeads ; a++)
          mp_CORtoDetectorDistance[a+hId*m_nbOfProjections/m_nbHeads] = mp_radius[hId];
    }
  }

  delete[] angles;
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerSPECTConv::GetSPECTSpecificParameters(uint16_t* ap_nbOfProjections, 
                                                  uint16_t* ap_nbHeads,
                                                  FLTNB* ap_acquisitionZoom,
                                                  uint16_t* ap_nbOfBins,
                                                  FLTNB*  ap_pixSizeXY, 
                                                  FLTNB*& ap_angles, 
                                                  FLTNB*& ap_CORtoDetectorDistance,
                                                  int*    ap_headRotDirection)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("iScannerSPECTConv::GetSPECTSpecificParameters ..." << endl);
  // Verify that all parameters have been correctly checked
  if (!m_allParametersChecked)
  {
    Cerr("***** iScannerSPECTConv::GetSPECTSpecificParameters() -> Parameters have not been checked !" << endl);
    return 1;
  }
  // Get them
  *ap_nbOfProjections = m_nbOfProjections;
  *ap_nbHeads = m_nbHeads;
  *ap_acquisitionZoom = m_acquisitionZoom;
  ap_nbOfBins[0] = mp_nbOfBins[0];
  ap_nbOfBins[1] = mp_nbOfBins[1];
  ap_pixSizeXY[0] = m_vPixelsSizeTrans;
  ap_pixSizeXY[1] = m_vPixelsSizeAxial;
  ap_angles = mp_projectionAngles;
  ap_CORtoDetectorDistance = mp_CORtoDetectorDistance;
  *ap_headRotDirection = m_rotDirection;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTAngles
  \param ap_projectionAngles : an array containing the projection angles
  \brief Set the projection angles with the array provided in parameter
  \return 0 if success, positive value otherwise
*/
int iScannerSPECTConv::PROJ_SetSPECTAngles(FLTNB* ap_projectionAngles)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if(m_verbose>=VERBOSE_NORMAL) Cout("iScannerSPECTConv::PROJ_SetSPECTAngles ..." << endl);

  // Check initialization of the number of projections
  if(m_nbOfProjections <= 0)
  {
    Cerr("***** iScannerSPECTConv::PROJ_SetSPECTAngles -> Error number of projection should be >0 ! '" << endl);
    return 1;
  }
  else
  {
    mp_projectionAngles = new FLTNB[m_nbOfProjections];
    
    for(int a=0 ; a<m_nbOfProjections ; a++)
      mp_projectionAngles[a] = ap_projectionAngles[a];
  }  
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTCORtoDetectorDistance
  \param a_distance
  \brief Set distance between the center of rotation and SPECT detectors if arg value>0, 
         Set with the geometric information in the scanner configuration file otherwise
  \return 0 if success, positive value otherwise
*/
int iScannerSPECTConv::PROJ_SetSPECTCORtoDetectorDistance(FLTNB a_distance)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerSPECTConv::PROJ_SetSPECTCORtoDetectorDistance ..." << endl);
  
  // Check initialization of the number of projections
  if(m_nbOfProjections <= 0)
  {
    Cerr("***** iScannerSPECTConv::PROJ_SetSPECTCORtoDetectorDistance -> Error number of projection should be >0 ! '" << endl);
    return 1;
  }
  else
  {
    mp_CORtoDetectorDistance = new FLTNB[m_nbOfProjections];
    
    if(a_distance>0)
    {
      // Set all radius to the provided distance
      for(int a=0 ; a<m_nbOfProjections ; a++)
        mp_CORtoDetectorDistance[a] = a_distance;
    }
    else
    {
      // Set all distance according to scanner geometric informations (default)
      for(int hId=0 ; hId<m_nbHeads ; hId++)
        for(int a=0 ; a<m_nbOfProjections/m_nbHeads ; a++)
          mp_CORtoDetectorDistance[a+hId*m_nbOfProjections/m_nbHeads] = mp_radius[hId];
    }
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetSPECTNbBins
  \brief Get the number of SPECT heads in the pointer provided in parameter
  \return 0 by default (no error)
*/
int iScannerSPECTConv::PROJ_GetSPECTNbBins(uint16_t* ap_nbOfBins) 
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  ap_nbOfBins = mp_nbOfBins;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTNbBins
  \param ap_nbOfBins
  \brief Set number of bins
  \return 0 by default (no error)
*/
int iScannerSPECTConv::PROJ_SetSPECTNbBins(uint16_t* ap_nbOfBins)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  for(int i=0 ; i<2 ; i++)
    mp_nbOfBins[i]=ap_nbOfBins[i]; 
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTNbProjections
  \param a_nbOfProjections
  \brief Set number of projections
  \return 0 by default (no error)
*/
int iScannerSPECTConv::PROJ_SetSPECTNbProjections(uint32_t a_nbOfProjections)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  m_nbOfProjections = a_nbOfProjections; 
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowHelp
  \brief Display help
  \todo Provide informations about SPECT system initialization ?
*/
void iScannerSPECTConv::ShowHelp()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  cout << "This scanner class is dedicated to the description of parallel, convergent and multi-convergent SPECT systems." << endl;
}
