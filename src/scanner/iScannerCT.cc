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
  \brief    Implementation of class iScannerCT
*/

#include "iScannerCT.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iScannerCT::iScannerCT() : vScanner() 
{
  // Set member variables to default values
  m_scannerType = SCANNER_CT;
  m_nbPixels = -1;
  m_nbPixelsTrans = 0;
  m_pixelsSizeTrans = -1.;
  m_gapSizeTrans = -1.;
  m_nbPixelsAxial = 0;
  m_pixelsSizeAxial = -1.;
  m_gapSizeAxial = -1.;
  m_nbOfProjections = 0;
  mp_projectionAngles = NULL;
  m_CORtoDetectorDistance = -1.;
  m_CORtoSourceDistance = -1.;
  m_rotDirection = GEO_ROT_CW;
  m_detectorDepth = -1.;
  m_spotSizeWidth = -1.;
  m_spotSizeDepth = -1.;
  mp_crystalCentralPositionX = NULL; 
  mp_crystalCentralPositionY = NULL; 
  mp_crystalCentralPositionZ = NULL; 
  mp_crystalOrientationX = NULL; 
  mp_crystalOrientationY = NULL; 
  mp_crystalOrientationZ = NULL; 
  mp_sourcePositionX = NULL; 
  mp_sourcePositionY = NULL; 
  mp_sourcePositionZ = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iScannerCT::~iScannerCT() 
{
  if (mp_projectionAngles) delete[] mp_projectionAngles;
  if (mp_crystalCentralPositionX) delete mp_crystalCentralPositionX;
  if (mp_crystalCentralPositionY) delete mp_crystalCentralPositionY;
  if (mp_crystalCentralPositionZ) delete mp_crystalCentralPositionZ;
  if (mp_crystalOrientationX) delete mp_crystalOrientationX;
  if (mp_crystalOrientationY) delete mp_crystalOrientationY;
  if (mp_crystalOrientationZ) delete mp_crystalOrientationZ;
  if (mp_sourcePositionX) delete mp_sourcePositionX;
  if (mp_sourcePositionY) delete mp_sourcePositionY;
  if (mp_sourcePositionZ) delete mp_sourcePositionZ;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      iScannerCT::DescribeSpecific()
  \brief   Implementation of the pure virtual eponym function that simply prints info about the scanner
*/
void iScannerCT::DescribeSpecific()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose==0) return;
  
  // Describe the scanner
  Cout("iScannerCT::DescribeSpecific() -> Here is some specific content of the CT scanner" << endl);
  Cout("  --> Total number of pixels: " << m_nbPixels << endl);
  Cout("  --> Total number of projections: " << m_nbOfProjections << endl);
  if( mp_projectionAngles ) 
  {
    // Display ten by ten
    Cout("  --> Projection angles: " << endl);
    if(m_nbOfProjections > 10)
    {
      int pr=0;
      while ( pr<m_nbOfProjections-1 )
      {
        for (int p=0 ; p<10 ; p++)
        {
          if(pr == m_nbOfProjections)
            break;
          Cout(mp_projectionAngles[pr] << ", ");
          pr++;
        }
        Cout(endl);
      }
    }
    else
    {
      for (int p=0 ; p<m_nbOfProjections-1 ; p++)
        Cout(mp_projectionAngles[p] << ", ");
      Cout(mp_projectionAngles[m_nbOfProjections-1] << endl);
    }
  }

  Cout("  --> Distance between the center of rotation and the detector surface: " << m_CORtoDetectorDistance << endl);
  Cout("  --> Distance between the center of rotation and the source surface: " << m_CORtoSourceDistance << endl);
  
  Cout("  --> Total number of transaxial pixels as defined in the system file: " << m_nbPixelsTrans << endl);
  Cout("  --> Size of transaxial pixels as defined in the system file: " << m_pixelsSizeTrans << endl);
  Cout("  --> Gap size between each transaxial pixe as defined in the system file: " << m_gapSizeTrans << endl);
  Cout("  --> Total number of axial pixels as defined in the system file: " << m_nbPixelsAxial << endl);
  Cout("  --> Size of axial pixels as defined in the system file: " << m_pixelsSizeAxial << endl);
  Cout("  --> Gap size between each axial pixel as defined in the system file: " << m_gapSizeAxial << endl);
  Cout("  --> Depth of detection pixels: " << m_detectorDepth << endl);
  Cout("  --> Width of the source, along the direction tangential to the scanner radius: " << m_spotSizeWidth << endl);
  Cout("  --> Depth of the source, along the direction of the scanner radius: " << m_spotSizeDepth << endl);
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::Instantiate(bool a_scannerFileIsLUT)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerCT::Instantiate() -> Create scanner structure and read parameters from configuration file" << endl); 

  // Get scanner manager
  sScannerManager* p_scannerManager = sScannerManager::GetInstance();

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
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "detector depth", &m_detectorDepth, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read detector depth !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "detector radius", &m_CORtoDetectorDistance, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read detector to COR distance as \"detector radius\" !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "source radius", &m_CORtoSourceDistance, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read source to COR distance as \"source radius\" !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "spot size width", &m_spotSizeWidth, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the spot size width !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "spot size height", &m_spotSizeDepth, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the spot size height !" << endl);
    return 1;
  }

  // Bed displacement
  m_defaultBedDisplacementInMm = 0.;
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "multiple bed displacement", &m_defaultBedDisplacementInMm, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** iScannerSPECTConv::Instantiate() -> An error occurred while trying to read the multiple bed displacement in the scanner header file !" << endl);
    return 1;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::BuildLUT(bool a_scannerFileIsLUT)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerCT::BuildLUT -> Build LUT for scanner elements coordinates and orientations"<< endl);
  
  // Initialize nb_pixels
  m_nbPixels = m_nbPixelsTrans * m_nbPixelsAxial;

  // Allocate LUT structures
  mp_crystalCentralPositionX = new FLTNB[m_nbOfProjections*m_nbPixels];
  mp_crystalCentralPositionY = new FLTNB[m_nbOfProjections*m_nbPixels];
  mp_crystalCentralPositionZ = new FLTNB[m_nbOfProjections*m_nbPixels];
  mp_crystalOrientationX = new FLTNB[m_nbOfProjections*m_nbPixels];
  mp_crystalOrientationY = new FLTNB[m_nbOfProjections*m_nbPixels];
  mp_crystalOrientationZ = new FLTNB[m_nbOfProjections*m_nbPixels];
  mp_sourcePositionX = new FLTNB[m_nbOfProjections];
  mp_sourcePositionY = new FLTNB[m_nbOfProjections];
  mp_sourcePositionZ = new FLTNB[m_nbOfProjections]; 

  // Either generate the LUT from a generic file, or load the user precomputed LUT
  if (!a_scannerFileIsLUT)
  {
    if (ComputeLUT())
    {
     Cerr("***** iScannerCT::BuildLUT() -> A problem occured while generating scanner LUT !" << endl);
     return 1;
    }
  }
  else 
  {
    if (LoadLUT())
    {
      Cerr("***** iScannerCT::BuildLUT() -> A problem occured while loading scanner LUT !" << endl);
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

int iScannerCT::CheckParameters()
{  
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  
  // Check if all parameters have been correctly initialized. Return error otherwise
  if (m_nbPixels<0)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Number of crystals has not been initialized !" <<endl);
    return 1;
  }
  if (m_nbOfProjections<=0)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Number of projection angles has not been initialized !" <<endl);
    return 1;
  }
  if (m_nbPixelsTrans<=0 || m_nbPixelsAxial<=0)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Number of transaxial/axial pixels have not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (m_pixelsSizeTrans<=0 || m_pixelsSizeAxial<=0)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Transaxial/axial pixel sizes have not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (m_gapSizeTrans<0. || m_gapSizeAxial<0.)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Transaxial/axial pixel gap sizes have not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (m_detectorDepth<=0)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Crystal depth has not correctly been initialized ! (should be >0)" <<endl);
    return 1;
  }
  if (mp_projectionAngles == NULL)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Projection angles have not correctly been initialized !" <<endl);
    return 1;
  }
  if (m_CORtoDetectorDistance<=0.)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Distance between center of rotation and detector surface has not correctly been initialized !" <<endl);
    return 1;
  }
  if (m_CORtoSourceDistance<=0.)
  {
    Cerr("***** iScannerCT::CheckParameters()-> Distance between center of rotation and source surface has not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_crystalCentralPositionX == NULL ||
      mp_crystalCentralPositionY == NULL ||
      mp_crystalCentralPositionZ == NULL  )
  {
    Cerr("***** iScannerCT::CheckParameters()-> LUT elements (crystal central positions) have not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_crystalOrientationX == NULL ||
      mp_crystalOrientationY == NULL ||
      mp_crystalOrientationZ == NULL  )
  {
    Cerr("***** iScannerCT::CheckParameters()-> LUT elements (crystal orientations) have not correctly been initialized !" <<endl);
    return 1;
  }
  if (mp_sourcePositionX == NULL ||
      mp_sourcePositionY == NULL ||
      mp_sourcePositionZ == NULL  )
  {
    Cerr("***** iScannerCT::CheckParameters()-> LUT elements (crystal focal positions) have not correctly been initialized !" <<endl);
    return 1;
  }

  // End
  m_allParametersChecked = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::Initialize()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerCT::Initialize() -> Initialize remaining stuff for scanner to be ready"<< endl); 

  // Parameters checked ?
  if (!m_allParametersChecked)
  {
    Cerr("***** iScannerCT::Initialize() -> Parameters have not been checked !" << endl);
    return 1;
  }
  
  // Any initialization
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::LoadLUT()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // This has not been designed yet
  Cerr("iScannerCT::LoadLUT() -> Not yet implemented !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::ComputeLUT()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if(m_verbose>=VERBOSE_DETAIL) Cout("iScannerCT::ComputeLUT() -> Start LUT generation" << endl);

  // Use multi-threading here as the rotation of the pixels and source can take some times if many projection angles are used
  int nb_threads = 1;
  if (mp_ID!=NULL) nb_threads = mp_ID->GetNbThreadsMax();
  #ifdef CASTOR_OMP
  omp_set_num_threads(nb_threads);
  #endif

  // ============================================================================================================
  // Generate the LUT
  // ============================================================================================================

  // oMatrix crystal_center_ref will be used to gather cartesian positions of the crystals center (on the surface of the crystals)
  // in the reference head (directly above isocenter = +y). Same for the source reference position but below the isocenter ;-)
  oMatrix** crystal_center_ref = new oMatrix *[m_nbPixels];
  for (int c = 0; c<m_nbPixels; c++)
    crystal_center_ref[c]  = new oMatrix(3,1);
  oMatrix* source_center_ref = new oMatrix(3,1);
  // oMatrix crystal_center_out will be used to recover positions of each crystal after each rotation (each projection angles).
  // Same for the source_center_out
  oMatrix** crystal_center_out = new oMatrix *[nb_threads];
  oMatrix** source_center_out = new oMatrix *[nb_threads];
  for (int th=0; th<nb_threads; th++)
  {
    crystal_center_out[th] = new oMatrix(3,1);
    source_center_out[th] = new oMatrix(3,1);
  }

  // Generation of the rotation matrix allowing to compute the position of all the projections. 
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
  sScannerManager* p_scannerManager = sScannerManager::GetInstance();
  if(ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "trans gap size", &m_gapSizeTrans, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerCT::ComputeLUT() -> An error occurred while trying to read the transaxial gap size !" << endl);
    return 1;
  } 

  if(ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "axial gap size", &m_gapSizeAxial, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerCT::ComputeLUT() -> An error occurred while trying to read the axial gap size !" << endl);
    return 1;
  }

  // Generate cartesian coordinates of the pixels centers for the CT at the reference position (directly above isocenter = +y)
  FLTNB x_start = -(((FLTNB)(m_nbPixelsTrans))*m_pixelsSizeTrans + ((FLTNB)(m_nbPixelsTrans-1))*m_gapSizeTrans) / 2. + (m_pixelsSizeTrans/2.);
  FLTNB z_start = -(((FLTNB)(m_nbPixelsAxial))*m_pixelsSizeAxial + ((FLTNB)(m_nbPixelsAxial-1))*m_gapSizeAxial) / 2. + (m_pixelsSizeAxial/2.);
  for(uint32_t jj = 0; jj < m_nbPixelsAxial; ++jj )
  {
    FLTNB Zcryst = z_start + jj * (m_pixelsSizeAxial + m_gapSizeAxial);
    for(uint32_t ii = 0; ii < m_nbPixelsTrans; ++ ii )
    {
      FLTNB Xcryst = x_start + ii * (m_pixelsSizeTrans + m_gapSizeTrans);
      crystal_center_ref[ii + m_nbPixelsTrans * jj]->SetMatriceElt(0,0,Xcryst);
      // Set y-position of the crystal reference matrix according to the distance between CT detector and center of rotation of the scanner (m_CORtoDetectorDistance)
      crystal_center_ref[ii + m_nbPixelsTrans * jj]->SetMatriceElt(1,0,m_CORtoDetectorDistance);
      crystal_center_ref[ii + m_nbPixelsTrans * jj]->SetMatriceElt(2,0,Zcryst);
    }
  }
  // Generate cartesian coordinates of the source for the CT at the reference position (directly above isocenter = +y)
  source_center_ref->SetMatriceElt(0,0,0.);
  // Set y-position of the source reference matrix according to the distance between CT source and center of rotation of the scanner (m_CORtoSourceDistance)
  source_center_ref->SetMatriceElt(1,0,-m_CORtoSourceDistance);
  source_center_ref->SetMatriceElt(2,0,0.); // The source is supposed to be centered along the Z axis

  // ============================================================================================================
  // Loop over all the projection angles
  // Positions of the scanner elements are progressively stored in the LUT file.
  // ============================================================================================================
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Generate positions for each projection angle"<< endl); 

  // Multi-threading loop over projection angles to compute the LUT
  int a;
  #pragma omp parallel for private(a) schedule(guided)
  for (a=0 ; a<m_nbOfProjections ; a++)
  {
    // Get the thread index
    int th = 0;
    #ifdef CASTOR_OMP
    th = omp_get_thread_num();
    #endif
    // Get surface crystal position (and not center crystal positions for SPECT), and set the reference crystal matrix
    for (int c=0 ; c < m_nbPixels ; c++)
    {
      // crystal indexation
      int cryID = a*m_nbPixels + c;

      // Compute rotation of the crystal center position from the reference and the rotation matrix associated to the projection
      rotation_mtx[a]->Multiplication(crystal_center_ref[c], crystal_center_out[th]);
      // Store crystal surface positions
      mp_crystalCentralPositionX[cryID] = crystal_center_out[th]->GetMatriceElt(0,0);
      mp_crystalCentralPositionY[cryID] = crystal_center_out[th]->GetMatriceElt(1,0);
      mp_crystalCentralPositionZ[cryID] = crystal_center_out[th]->GetMatriceElt(2,0);
      // Store crystal orientations
      mp_crystalOrientationX[cryID] = rotation_mtx[a]->GetMatriceElt(0,1);
      mp_crystalOrientationY[cryID] = rotation_mtx[a]->GetMatriceElt(0,0);
      mp_crystalOrientationZ[cryID] = 0;
    }
    // Compute rotation of the source center position from the reference and the rotation matrix associated to the projection
    rotation_mtx[a]->Multiplication(source_center_ref, source_center_out[th]);
    // Store crystal surface positions
    mp_sourcePositionX[a] = source_center_out[th]->GetMatriceElt(0,0);
    mp_sourcePositionY[a] = source_center_out[th]->GetMatriceElt(1,0);
    mp_sourcePositionZ[a] = source_center_out[th]->GetMatriceElt(2,0);
  }

  for(int i=0; i<m_nbOfProjections; i++)
    delete rotation_mtx[i];
  delete[] rotation_mtx;
  for (int c=0; c<m_nbPixels; c++)
    delete crystal_center_ref[c];
  delete[] crystal_center_ref;
  delete source_center_ref;
  for (int th=0; th<nb_threads; th++)
  {
    delete crystal_center_out[th];
    delete source_center_out[th];
  }
  delete[] crystal_center_out;
  delete[] source_center_out;

  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> LUT generation completed" << endl);

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::GetPositionsAndOrientations( int a_index1, int a_index2,
                                             FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                             FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3],
                                             FLTNB* ap_POI1, FLTNB* ap_POI2 )
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // CT indices : 1st for the projection, 2nd for the pixel

  // First check projection angle existency
  if (a_index1<0 || a_index1>=m_nbOfProjections)
  {
    Cerr("***** iScannerCT::GetPositionsAndOrientations() -> Projection index (" << a_index1 << ") out of range [0:" << m_nbOfProjections-1 << "] !" << endl);
    return 1;
  }
  
  // Second check pixels existency
  if (a_index2<0 || a_index2>=m_nbPixels)
  {
    Cerr("***** iScannerCT::GetPositionsAndOrientations() -> Pixel index (" << a_index2 << ") out of range [0:" << m_nbPixels-1 << "] !" << endl);
    return 1;
  }

  // Get crystal index related to the projection
  int index = a_index1*m_nbPixels + a_index2;
  
  // Get position of the source
  ap_Position1[0] = mp_sourcePositionX[a_index1];
  ap_Position1[1] = mp_sourcePositionY[a_index1];
  ap_Position1[2] = mp_sourcePositionZ[a_index1];

  // todo : Due to the current implementation of CT projection, POI 
  //        and DOI are not handled and ignored.
  //        An error is returned if POI are provided

  // Case when POI is not provided
  if (ap_POI2==NULL)
  {
    ap_Position2[0] = mp_crystalCentralPositionX[index];
    ap_Position2[1] = mp_crystalCentralPositionY[index];
    ap_Position2[2] = mp_crystalCentralPositionZ[index];
  }
  // Case when POI[2] is negative (meaning we only have POI[0] or POI[1] specified and to be taken into account)
  else if (ap_POI2[2]<0.)
  {
    Cerr("***** iScannerCT::GetPositionsAndOrientations() -> POI management not implemented yet for CT !" << endl);
    return 1;
  }
  // Case when only the DOI is provided
  else if (ap_POI2[0]==0. && ap_POI2[1]==0.)
  {
    Cerr("***** iScannerCT::GetPositionsAndOrientations() -> POI management not implemented yet for CT !" << endl);
    return 1;
  }
  // Case when the full POI is taken into account
  else
  {
    Cerr("***** iScannerCT::GetPositionsAndOrientations() -> POI management not implemented yet for CT !" << endl);
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

int iScannerCT::GetRdmPositionsAndOrientations(int a_index1, int a_index2,
                                                      FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                                      FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3] )
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  Cerr("***** iScannerCT::GetRdmPositionsAndOrientations() -> Not yet implemented for CT !" << endl);

  return 0;
}
  
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::GetPositionWithRandomDepth(int a_index1, int a_index2, FLTNB ap_Position1[3], FLTNB ap_Position2[3])
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  // This function was first implemented for PET testing purpose. Not implemented yet.
  Cerr("***** iScannerCT::GetPositionWithRandomDepth() -> This function was implemented for PET testing purpose. Not implemented for CT !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::GetTwoCorners(int a_index1, int a_index2,
                              FLTNB ap_CornerInf1[3], FLTNB ap_CornerSup1[3],
                              FLTNB ap_CornerInf2[3], FLTNB ap_CornerSup2[3])
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  Cerr("***** iScannerCT::GetTwoCorners() -> Not implemented yet !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::GetEdgesCenterPositions( int a_index1, int a_index2,
                                         FLTNB ap_pos_line_point1[3], FLTNB ap_pos_line_point2[3],
                                         FLTNB ap_pos_point1_x[4], FLTNB ap_pos_point1_y[4], FLTNB ap_pos_point1_z[4],
                                         FLTNB ap_pos_point2_x[4], FLTNB ap_pos_point2_y[4], FLTNB ap_pos_point2_z[4]
)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  if ( a_index1<0 || a_index1>=m_nbOfProjections
    || a_index2<0 || a_index2>=m_nbPixels )
  {
    Cerr("***** iScannerCT::GetPositionEdge() -> Crystal index or projection index out of range !" << endl);
    return 1;
  }

  // Computing global index
  int global_index = m_nbPixels * a_index1 + a_index2;

  // Getting the half size of axial/trans crystal 1 and 2 depending on the layer
  FLTNB half_spot_size_trans = m_spotSizeWidth / 2.0;
  FLTNB half_spot_size_axial = m_spotSizeDepth / 2.0;
  FLTNB half_pixel_trans = m_pixelsSizeTrans / 2.0;
  FLTNB half_pixel_axial = m_pixelsSizeAxial / 2.0;

  //////////////////////////////////////////////////////////////////////////
  // Computing coordinates depending on the orientation for point1
  // X-axis
  ap_pos_point1_x[ 0 ] = ap_pos_line_point1[ 0 ] - half_spot_size_trans * -mp_crystalOrientationY[ global_index ];
  ap_pos_point1_x[ 1 ] = ap_pos_line_point1[ 0 ] + half_spot_size_trans * -mp_crystalOrientationY[ global_index ];
  ap_pos_point1_x[ 2 ] = ap_pos_line_point1[ 0 ];
  ap_pos_point1_x[ 3 ] = ap_pos_line_point1[ 0 ];

  // Y-axis
  ap_pos_point1_y[ 0 ] = ap_pos_line_point1[ 1 ] - half_spot_size_trans * mp_crystalOrientationX[ global_index ];
  ap_pos_point1_y[ 1 ] = ap_pos_line_point1[ 1 ] + half_spot_size_trans * mp_crystalOrientationX[ global_index ];
  ap_pos_point1_y[ 2 ] = ap_pos_line_point1[ 1 ];
  ap_pos_point1_y[ 3 ] = ap_pos_line_point1[ 1 ];

  // Z-axis
  ap_pos_point1_z[ 0 ] = ap_pos_line_point1[ 2 ];
  ap_pos_point1_z[ 1 ] = ap_pos_line_point1[ 2 ];
  ap_pos_point1_z[ 2 ] = ap_pos_line_point1[ 2 ] - half_spot_size_axial;
  ap_pos_point1_z[ 3 ] = ap_pos_line_point1[ 2 ] + half_spot_size_axial;

  //////////////////////////////////////////////////////////////////////////
  // Computing coordinates depending on the orientation for point2
  // X-axis
  ap_pos_point2_x[ 0 ] = ap_pos_line_point2[ 0 ] + half_pixel_trans * mp_crystalOrientationY[ global_index ];
  ap_pos_point2_x[ 1 ] = ap_pos_line_point2[ 0 ] - half_pixel_trans * mp_crystalOrientationY[ global_index ];
  ap_pos_point2_x[ 2 ] = ap_pos_line_point2[ 0 ];
  ap_pos_point2_x[ 3 ] = ap_pos_line_point2[ 0 ];

  // Y-axis
  ap_pos_point2_y[ 0 ] = ap_pos_line_point2[ 1 ] + half_pixel_trans * -mp_crystalOrientationX[ global_index ];
  ap_pos_point2_y[ 1 ] = ap_pos_line_point2[ 1 ] - half_pixel_trans * -mp_crystalOrientationX[ global_index ];
  ap_pos_point2_y[ 2 ] = ap_pos_line_point2[ 1 ];
  ap_pos_point2_y[ 3 ] = ap_pos_line_point2[ 1 ];

  // Z-axis
  ap_pos_point2_z[ 0 ] = ap_pos_line_point2[ 2 ];
  ap_pos_point2_z[ 1 ] = ap_pos_line_point2[ 2 ];
  ap_pos_point2_z[ 2 ] = ap_pos_line_point2[ 2 ] - half_pixel_axial;
  ap_pos_point2_z[ 3 ] = ap_pos_line_point2[ 2 ] + half_pixel_axial;

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

// TODO SS: If possible, simplify all this ping-pong between scanner and datafile information. For example, this function should
//          be implemented inside the datafile, and then the information are passed to the scanner, instead of mixing everything.
int iScannerCT::GetGeometricInfoFromDataFile(string a_pathToDF)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_DETAIL) Cout("iScannerCT::GetGeometricInfoFromDataFile() -> Specific to CT" << endl);

  // This function is intended to be called after the scanner initialization, so that any field present in the datafile header, similar to
  // one in the scanner configuration file, may overload the scanner value.

  if (ReadDataASCIIFile(a_pathToDF , "Number of projections", &m_nbOfProjections, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerCT::GetGeometricInfoFromDataFile() -> Error while reading number of projections in the header data file " << endl);
    return 1;
  }
  
  // Allocate a one-dimensional vector to retrieve angles from the datafile, then copy it in the member variables
  FLTNB* angles = new FLTNB[m_nbOfProjections];
  FLTNB first_and_last_angles[2] = {-1.,-1.};

  // Two possible initializations for projection angles :
  // - All angles are provided with the keyword "Angles"
  // - The first and last angles are provided, and the intermediate angles are extrapolated from the number of projections and the rotation direction
  
  // 'First/Last angles' tags : checking issue during data reading/conversion (==1)
  if (ReadDataASCIIFile(a_pathToDF, "First and last projection angles", first_and_last_angles, 2, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** iScannerCT::GetGeometricInfoFromDataFile() -> Error while reading Angle mandatory field in the header data file '" << endl);
    return 1;
  }
  // Scanner rotation direction
  string rotation_direction = "";
  if (ReadDataASCIIFile(a_pathToDF , "Scanner rotation direction", &rotation_direction, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** iScannerCT::GetGeometricInfoFromDataFile() -> Error while reading head rotation orientation in the header data file " << endl);
    return 1;
  }
  if (SetRotDirection(rotation_direction) )
  {
    Cerr("***** iScannerCT::GetGeometricInfoFromDataFile() -> Error occured while trying to initialize scanner rotation orientation " << endl);
    return 1;
  }
  
  // Check for 'Projection angles' tag
  int rvalue = ReadDataASCIIFile(a_pathToDF, "Projection angles", angles, m_nbOfProjections, KEYWORD_OPTIONAL);
  
  // Error while reading "Angles" tag (==1)
  if(rvalue==1)
  {
    Cerr("***** iScannerCT::GetGeometricInfoFromDataFile() -> Error while reading Angles field in the header data file !'" << endl);
    return 1;
  }

  // Check if information on projection angles has been provided
  if ( rvalue>=2 &&  // "Angles" tag not found
      (first_and_last_angles[0] <0 || first_and_last_angles[1] <0) ) // Tags first/last angles not found)
  {
    Cerr("***** iScannerCT::GetGeometricInfoFromDataFile() -> No information on projection angles provided in the datafile !'" << endl);
    Cerr("                                                    This information should be provided using either the 'Angles' tag, or both 'First angles', 'Last angles' tags !'" << endl);
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
  for(int a=0 ; a<m_nbOfProjections ; a++) mp_projectionAngles[a] = angles[a];

  delete[] angles;
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerCT::GetCTSpecificParameters(uint16_t* ap_nbOfProjections, 
                                        FLTNB*& ap_angles, 
                                        int*    ap_detectorRotDirection)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("iScannerCT::GetCTSpecificParameters() -> Copy pointers inside scanner to the datafile to get information" << endl
                                   << "                                         depending on the datafile inside the scanner" << endl);
  // Verify that all parameters have been correctly checked
  if (!m_allParametersChecked)
  {
    Cerr("***** iScannerCT::GetCTSpecificParameters() -> Parameters have not been checked !" << endl);
    return 1;
  }
  // Get them
  *ap_nbOfProjections = m_nbOfProjections;
  ap_angles = mp_projectionAngles;
  *ap_detectorRotDirection = m_rotDirection;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iScannerCT::ShowHelp()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  cout << "This scanner class is dedicated to the description of CT systems." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
