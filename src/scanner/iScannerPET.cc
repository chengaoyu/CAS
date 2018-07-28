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
  \brief    Implementation of class iScannerPET
*/

#include "iScannerPET.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iScannerPET::iScannerPET() : vScanner() 
{
  // Set member variables to default values
  m_scannerType = SCANNER_PET;
  m_nbLayers = -1;
  m_nbCrystals = -1;
  m_minAngleDifference = -1.;
  mp_nbCrystalsInLayer = NULL;
  mp_crystalCentralPositionX = NULL;
  mp_crystalCentralPositionY = NULL;
  mp_crystalCentralPositionZ = NULL;
  mp_crystalOrientationX = NULL;
  mp_crystalOrientationY = NULL;
  mp_crystalOrientationZ = NULL;
  m_rotDirection = GEO_ROT_CW;
  mp_sizeCrystalTrans = NULL;
  mp_sizeCrystalAxial = NULL;
  mp_sizeCrystalDepth = NULL;
  mp_meanDepthOfInteraction = NULL;
  mp_positionMatrix_ref = NULL;
  mp_positionMatrix_out = NULL;
  mp_rotationMatrix = NULL;
  // Variables depending on the acquisition
  m_maxAxialDiffmm = -1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iScannerPET::~iScannerPET() 
{
  if (mp_sizeCrystalTrans != NULL) delete mp_sizeCrystalTrans; 
  if (mp_sizeCrystalAxial != NULL) delete mp_sizeCrystalAxial; 
  if (mp_sizeCrystalDepth != NULL) delete mp_sizeCrystalDepth; 
  if (mp_meanDepthOfInteraction != NULL)  delete mp_meanDepthOfInteraction;
  if (mp_crystalCentralPositionX != NULL) delete mp_crystalCentralPositionX; 
  if (mp_crystalCentralPositionY != NULL) delete mp_crystalCentralPositionY; 
  if (mp_crystalCentralPositionZ != NULL) delete mp_crystalCentralPositionZ; 
  if (mp_crystalOrientationX != NULL) delete mp_crystalOrientationX; 
  if (mp_crystalOrientationY != NULL) delete mp_crystalOrientationY; 
  if (mp_crystalOrientationZ != NULL) delete mp_crystalOrientationZ; 
  if (mp_nbCrystalsInLayer != NULL) delete mp_nbCrystalsInLayer;
  if (mp_positionMatrix_ref != NULL) delete mp_positionMatrix_ref;
  if (mp_positionMatrix_out != NULL) delete mp_positionMatrix_out;
  if (mp_rotationMatrix != NULL) delete mp_rotationMatrix;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      iScannerPET::DescribeSpecific()
  \brief   Implementation of the pure virtual eponym function that simply prints info about the scanner
*/
void iScannerPET::DescribeSpecific()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose==0) return;
  
  // Describe the scanner
  Cout("iScannerPET::DescribeSpecific() -> Here is some specific content of the PET scanner" << endl);
  Cout("  --> Number of layers: " << m_nbLayers << endl);
  Cout("  --> Total number of crystals: " << m_nbCrystals << endl);
  for (int l=0 ; l<m_nbLayers ; l++)
    if( mp_nbCrystalsInLayer ) Cout("  --> Nb crystals layer ["<<l+1<<"]: " << mp_nbCrystalsInLayer[l] << endl);

  for (int l=0 ; l<m_nbLayers ; l++)
  {
    if(l>0) Cout("  --> Layer "<<l<<":"<< endl);
    if( mp_sizeCrystalTrans ) Cout("  --> Crystal transaxial dimension (mm): " << mp_sizeCrystalTrans[l] << endl);
    if( mp_sizeCrystalAxial ) Cout("  --> Crystal axial dimension (mm): " << mp_sizeCrystalAxial[l] << endl);
    if( mp_sizeCrystalDepth ) Cout("  --> Crystal depth dimension (mm): " << mp_sizeCrystalDepth[l] << endl);
  }
  if(mp_meanDepthOfInteraction )
    for (int l=0 ; l<m_nbLayers ; l++)
      Cout("  --> Mean depth of interaction for layer ["<<l+1<<"]: " << mp_meanDepthOfInteraction[l] << endl);
  
  Cout("  --> (Transaxial) minimum angle difference (radian): " << m_minAngleDifference << endl);
  if(m_maxAxialDiffmm > 0.)Cout("  --> Maximum axial difference for LORs (mm): " << m_maxAxialDiffmm << endl);
}







// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int iScannerPET::Instantiate()
  \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                               the system is a generic file (0) or custom Look-up-table (1)
  \brief   Get mandatory informations from the scanner file 
           and allocate memory for the member variables
  \return  0 if success, positive value otherwise
*/
int iScannerPET::Instantiate(bool a_scannerFileIsLUT)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerPET::Instantiate() -> Create scanner structure and read parameters from configuration file"<< endl); 

  // Get scanner manager
  sScannerManager* p_scannerManager; 
  p_scannerManager = sScannerManager::GetInstance();  

  // Get the number of layers in the scanner
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of layers", &m_nbLayers, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerPET::Instantiate() -> An error occurred while trying to read the number of layers in the scanner header file !" << endl);
    return 1;
  }
  // Check the correct initialization of the number of layers
  if (m_nbLayers<=0)
  {
    Cerr("***** iScannerPET::Instantiate() -> Incorrect value for the number of layer (must be >0) !" << endl);
    return 1;
  }
  // Get the number of elements in the system
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of elements", &m_nbCrystals, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iScannerPET::Instantiate() -> An error occurred while trying to read the number of elements in the scanner header file !" << endl);
    return 1;
  }

  // Instanciation of layer dependent member variables
  mp_nbCrystalsInLayer = new int[m_nbLayers];
  mp_sizeCrystalTrans = new FLTNB[m_nbLayers]; 
  mp_sizeCrystalAxial = new FLTNB[m_nbLayers]; 
  mp_sizeCrystalDepth = new FLTNB[m_nbLayers];
  mp_meanDepthOfInteraction = new FLTNB[m_nbLayers];
  
  // Instanciation of number of crystals dependent member variables
  mp_crystalCentralPositionX = new FLTNB[m_nbCrystals];
  mp_crystalCentralPositionY = new FLTNB[m_nbCrystals];
  mp_crystalCentralPositionZ = new FLTNB[m_nbCrystals];
  mp_crystalOrientationX = new FLTNB[m_nbCrystals];
  mp_crystalOrientationY = new FLTNB[m_nbCrystals];
  mp_crystalOrientationZ = new FLTNB[m_nbCrystals];

  // Initialize layer size with default value (=0);
  for (int l=0 ; l<m_nbLayers ; l++)
  {
    mp_sizeCrystalTrans[l] = 0.;
    mp_sizeCrystalAxial[l] = 0.;
    mp_sizeCrystalDepth[l] = 0.;
  }
  
  // Size of crystals
  if (a_scannerFileIsLUT)
  {
    // For the moment, only the depth is mandatory as the others are not yet implemented
    if ( ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "crystals size trans", mp_sizeCrystalTrans, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
         ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "crystals size axial", mp_sizeCrystalAxial, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
         ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "crystals size depth", mp_sizeCrystalDepth, m_nbLayers, KEYWORD_MANDATORY) )
    {
      Cerr("***** iScannerPET::Instantiate() -> An error occurred while trying to read the crystals size in the scanner header file !" << endl);
      return 1;
    }
  }
  else
  {
    // Mandatory parameter
    if ( ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "crystals size trans", mp_sizeCrystalTrans, m_nbLayers, KEYWORD_MANDATORY) ||
         ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "crystals size axial", mp_sizeCrystalAxial, m_nbLayers, KEYWORD_MANDATORY) ||
         ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "crystals size depth", mp_sizeCrystalDepth, m_nbLayers, KEYWORD_MANDATORY) )
    {
      Cerr("***** iScannerPET::Instantiate() -> An error occurred while trying to read the crystals size in the scanner header file !" << endl);
      return 1;
    }
  }

  // ----------------------------------------------
  // Optional parameters
  // ----------------------------------------------
  
  // Initialize with default values
  for (int l=0 ; l<m_nbLayers ; l++) mp_meanDepthOfInteraction[l] = -1;
  m_minAngleDifference = 0.;
  m_defaultBedDisplacementInMm = 0.;
  
  // Check if value is provided in the scanner conf file
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "mean depth of interaction", mp_meanDepthOfInteraction, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
  {
    Cerr("***** iScannerPET::Instantiate() -> An error occurred while trying to read the mean depth of interaction in the scanner header file !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "min angle difference", &m_minAngleDifference, 1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
  {
    Cerr("***** iScannerPET::Instantiate() -> An error occurred while trying to read the min angle difference in the scanner header file !" << endl);
    return 1;
  }
  if (ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "multiple bed displacement", &m_defaultBedDisplacementInMm, 1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
  {
    Cerr("***** iScannerPET::Instantiate() -> An error occurred while trying to read the multiple bed displacement in the scanner header file !" << endl);
    return 1;
  }

  // Instanciation of objects for matrix calculation during reconstruction
  mp_positionMatrix_ref = new oMatrix(3,1);
  mp_positionMatrix_out = new oMatrix(3,1);
  mp_rotationMatrix = new oMatrix(3,3); 

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
  \todo default values to max ring diff
*/
int iScannerPET::BuildLUT(bool a_scannerFileIsLUT)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerPET::BuildLUT() -> Build LUT for scanner elements coordinates and orientations"<< endl); 
    
  // Either generate the LUT from a generic file, or load the user precomputed LUT
  if (!a_scannerFileIsLUT)
  {
    if (ComputeLUT())
    {
     Cerr("***** iScannerPET::BuildLUT() -> A problem occured while generating scanner LUT !" << endl);
     return 1;
    }
  }
  else 
  {
    if (LoadLUT())
    {
      Cerr("***** iScannerPET::BuildLUT() -> A problem occured while loading scanner LUT !" << endl);
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
  \brief Check if all parameters have been correctly initialized
  \return 0 if success, positive value otherwise
*/
int iScannerPET::CheckParameters()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Check if all parameters have been correctly initialized. Return error otherwise
  if (m_scannerType == -1)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Scanner type not initialized !" << endl);
    return 1;
  }
  if (m_verbose == -1)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Verbosity not initialized !" << endl);
    return 1;
  }
  if (m_nbLayers <= 0)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Incorrect value for the number of layer (must be >0) !" << endl);
    return 1;
  }
  if (m_nbCrystals <0)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Number of crystals not initialized !" << endl);
    return 1;
  }
  if (mp_nbCrystalsInLayer == NULL)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Number of crystals in layer(s) not initialized !" << endl);
    return 1;
  }
  if (mp_crystalCentralPositionX==NULL || mp_crystalCentralPositionY==NULL || mp_crystalCentralPositionZ==NULL)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Crystals central positions not initialized !" << endl);
    return 1;
  }
  if (mp_crystalOrientationX==NULL || mp_crystalOrientationY==NULL || mp_crystalOrientationZ==NULL)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Crystals orientations not initialized !" << endl);
    return 1;
  }
  if (mp_sizeCrystalTrans==NULL || mp_sizeCrystalAxial==NULL || mp_sizeCrystalDepth==NULL)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Crystals dimensions not initialized !" << endl);
    return 1;
  }
  if (mp_meanDepthOfInteraction == NULL)
  {
    Cout("***** iScannerPET::CheckParameters() -> Mean depth of interaction not initialized !" << endl);
    return 1;
  }
  if (m_minAngleDifference < 0)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Minimum angle difference not initialized !" << endl);
    return 1;
  }
  if (m_scannerType == -1)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Scanner type not initialized !" << endl);
    return 1;
  }
  if (m_scannerType == -1)
  {
    Cerr("***** iScannerPET::CheckParameters() -> Scanner type not initialized !" << endl);
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
/*
  \fn Initialize
  \brief Check general initialization and set several parameters to their default value
  \return 0 if success, positive value otherwise
*/
int iScannerPET::Initialize()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerPET::Initialize() -> Initialize remaining stuff for scanner to be ready"<< endl); 

  // Parameters checked ?
  if (!m_allParametersChecked)
  {
    Cerr("***** iScannerPET::Initialize() -> Parameters have not been checked !" << endl);
    return 1;
  }
  
  // Set mean depth of interaction to the center of crystal if not initialized by default.
  for(int l=0 ; l<m_nbLayers ; l++)
    if (mp_meanDepthOfInteraction[l]<0) mp_meanDepthOfInteraction[l] = mp_sizeCrystalDepth[l]/2;
    
  // Conversion of the minimal angle difference into radians
  m_minAngleDifference *= M_PI/180;
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LoadLUT
  \brief Load a precomputed scanner LUT
  \details Read mandatory data from the header of the LUT. 
           Then load the LUT elements for each crystal
  \return 0 if success, positive value otherwise
*/
int iScannerPET::LoadLUT()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("iScannerPET::LoadLUT() -> Start loading LUT from user-provided file" <<  endl);

  // Get scanner manager
  sScannerManager* p_scannerManager; 
  p_scannerManager = sScannerManager::GetInstance(); 

  if(ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of crystals in layer", mp_nbCrystalsInLayer, m_nbLayers, 1) )
  {
    Cerr("***** iScannerPET::LoadLUT() -> An error occurred while trying to read a mandatory parameter in the scanner header file  !" << endl);
    return 1;
  }

  // Open the scanner LUT file
  string scanner_lut_file = p_scannerManager->GetPathToScannerFile();
  //todo : maybe more practical to directly read the path from the header.
  //       Do not put restrictions on the file extension (.lut)
  scanner_lut_file = scanner_lut_file.substr(0, scanner_lut_file.find_last_of(".")).append(".lut");

  // Open file
  FILE* LUT_file = fopen(scanner_lut_file.c_str(), "rb");
  if (LUT_file==NULL)
  {
    Cerr("***** iScannerPET::LoadLUT() -> Input LUT file '" << scanner_lut_file << "' is missing or corrupted !" << endl);
    return 1;
  }

  // Read data for each index
  int nb_data_read = 0;
  for (int i=0; i<m_nbCrystals; i++)
  {
    FLTNBLUT buffer;
    // Read central crystal position X, then Y and Z
    nb_data_read += fread(&buffer,sizeof(FLTNBLUT),1,LUT_file);
    mp_crystalCentralPositionX[i] = ((FLTNB)buffer);
    nb_data_read += fread(&buffer,sizeof(FLTNBLUT),1,LUT_file);
    mp_crystalCentralPositionY[i] = ((FLTNB)buffer);
    nb_data_read += fread(&buffer,sizeof(FLTNBLUT),1,LUT_file);
    mp_crystalCentralPositionZ[i] = ((FLTNB)buffer);
    // Read crystal orientation X, then Y and Z
    nb_data_read += fread(&buffer,sizeof(FLTNBLUT),1,LUT_file);
    mp_crystalOrientationX[i] = ((FLTNB)buffer);
    nb_data_read += fread(&buffer,sizeof(FLTNBLUT),1,LUT_file);
    mp_crystalOrientationY[i] = ((FLTNB)buffer);
    nb_data_read += fread(&buffer,sizeof(FLTNBLUT),1,LUT_file);
    mp_crystalOrientationZ[i] = ((FLTNB)buffer);
  }

  // Close file
  fclose(LUT_file);

  // Check reading
  if (nb_data_read!=m_nbCrystals*6)
  {
    Cerr("***** iScannerPET::LoadLUT() -> Failed to read all data in input LUT file '" << scanner_lut_file << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> LUT successfully loaded" << endl);
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ComputeLUT
  \brief Compute the LUT of the scanner from a generic (.geom) file
  \details Read mandatory data from the geom file. Then compute the LUT elements for each crystal from the geometry described in the file
  \details Compute the look-up-tables of the system containing the locations of the scanner elements center in space and their orientations
  \todo Add option to inverse rsector if transaxial orientation is counter-clockwise.?
  \todo rotation for non-perfectly cylindric scanner (angles along the Z axis)
  \return 0 if success, positive value otherwise
*/
int iScannerPET::ComputeLUT()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_DETAIL) Cout("iScannerPET::ComputeLUT() -> Start LUT generation" << endl);
  
  // ============================================================================================================
  // Init local variables to recover data from scanner file
  // ============================================================================================================
  
  // Layer dependent variables
  int *nb_ang_rsector_lyr, // nb of angular position for the rsectors
      *nb_axl_rsector_lyr, // nb of axial subdivision of the rsector for a same angular position
      *nb_trs_mod_lyr,
      *nb_axl_mod_lyr,
      *nb_trs_submod_lyr,
      *nb_axl_submod_lyr,
      *nb_trs_crystal_lyr,
      *nb_axl_crystal_lyr;
  
  FLTNBLUT  *radius_lyr,
            *rsector_first_angle_lyr,
            *rsector_angular_span_lyr,
            *gap_axl_rsector_lyr,
            *gap_trs_mod_lyr,
            *gap_axl_mod_lyr,
            *gap_trs_submod_lyr,
            *gap_axl_submod_lyr,
            *gap_trs_crystal_lyr,
            *gap_axl_crystal_lyr;

  nb_ang_rsector_lyr = new int[m_nbLayers];
  nb_axl_rsector_lyr = new int[m_nbLayers];
  nb_trs_mod_lyr = new int[m_nbLayers];
  nb_axl_mod_lyr = new int[m_nbLayers];
  nb_trs_submod_lyr = new int[m_nbLayers];
  nb_axl_submod_lyr = new int[m_nbLayers];
  nb_trs_crystal_lyr = new int[m_nbLayers];
  nb_axl_crystal_lyr = new int[m_nbLayers];

  radius_lyr = new FLTNBLUT[m_nbLayers];
  gap_axl_rsector_lyr = new FLTNBLUT[m_nbLayers];
  gap_trs_mod_lyr = new FLTNBLUT[m_nbLayers];
  gap_axl_mod_lyr = new FLTNBLUT[m_nbLayers];
  gap_trs_submod_lyr = new FLTNBLUT[m_nbLayers];
  gap_axl_submod_lyr = new FLTNBLUT[m_nbLayers];
  gap_trs_crystal_lyr = new FLTNBLUT[m_nbLayers];
  gap_axl_crystal_lyr = new FLTNBLUT[m_nbLayers];
  rsector_first_angle_lyr = new FLTNBLUT[m_nbLayers];
  rsector_angular_span_lyr = new FLTNBLUT[m_nbLayers];

  string rotation_direction = "";

  // -------------------------------------------------------------------
  // Recover mandatory data from scanner file
  sScannerManager* p_scannerManager; 
  p_scannerManager = sScannerManager::GetInstance();

  if( ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of rsectors",              nb_ang_rsector_lyr, m_nbLayers, KEYWORD_MANDATORY) ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of crystals transaxial",  nb_trs_crystal_lyr, m_nbLayers, KEYWORD_MANDATORY) ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of crystals axial",       nb_axl_crystal_lyr, m_nbLayers, KEYWORD_MANDATORY) ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "scanner radius",                           radius_lyr, m_nbLayers, KEYWORD_MANDATORY) )
  {
    Cerr("***** iScannerPET::ComputeLUT() -> An error occurred while trying to read a mandatory parameter for scanner LUT generation in the scanner header file  !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------
  // Recover optionnal data from scanner file
  
  // Initialize defaulf values
  for(int l=0 ; l<m_nbLayers ; l++)
  {
    nb_axl_rsector_lyr[ l ] =1;
    nb_trs_mod_lyr[ l ] =1;
    nb_axl_mod_lyr[ l ] =1;
    nb_trs_submod_lyr[ l ] = 1;
    nb_axl_submod_lyr[ l ] = 1;
    gap_axl_rsector_lyr[ l ] = 0.;
    gap_trs_mod_lyr[ l ] = 0.;
    gap_axl_mod_lyr[ l ] = 0.;
    gap_trs_submod_lyr[ l ] = 0.;
    gap_axl_submod_lyr[ l ] = 0.;
    gap_trs_crystal_lyr[ l ] = 0.;
    gap_axl_crystal_lyr[ l ] = 0.;
    rsector_first_angle_lyr[ l ] = 0;
    rsector_angular_span_lyr[ l ] = 360.;
  }

  if( ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of rsectors axial",       nb_axl_rsector_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of modules transaxial",       nb_trs_mod_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of modules axial",            nb_axl_mod_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of submodules transaxial", nb_trs_submod_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "number of submodules axial",      nb_axl_submod_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "rsector gap axial",             gap_axl_rsector_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "module gap transaxial",             gap_trs_mod_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "module gap axial",                  gap_axl_mod_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "submodule gap transaxial",       gap_trs_submod_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "submodule gap axial",            gap_axl_submod_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "crystal gap transaxial",        gap_trs_crystal_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "crystal gap axial",             gap_axl_crystal_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "rsectors first angle",        rsector_first_angle_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "rsectors angular span",      rsector_angular_span_lyr, m_nbLayers, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR ||
      ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "rotation direction",         &rotation_direction,          1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
  {
    Cerr("***** iScannerPET::ComputeLUT() -> An error occurred while trying to read an optionnal parameter for scanner LUT generation in the scanner header file  !" << endl);
    return 1;
  }
  
  // Compute nb rings per layer
  uint16_t *p_nbRings = new uint16_t[m_nbLayers];
  for(int lyr=0 ; lyr<m_nbLayers ; lyr++)
    p_nbRings[lyr] = nb_axl_rsector_lyr[ lyr ]
                    * nb_axl_mod_lyr[ lyr ]
                    * nb_axl_submod_lyr[ lyr ]
                    * nb_axl_crystal_lyr[ lyr ];

  // Set rotation direction for geometry generation
  if (SetRotDirection(rotation_direction) )
  {
    Cerr("***** iScannerPET::ComputeLUT() -> Error occured while trying to initialize head rotation orientation " << endl);
    return 1;
  }

  // Z-SHIFT MANAGEMENT

  // Variables
  int rsector_nb_zshift;
  FLTNBLUT *zshift;
  
  int rvalue = ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "rsectors nbZShift", &rsector_nb_zshift, 1, KEYWORD_OPTIONAL);

  // reading error
  if( rvalue == KEYWORD_OPTIONAL_ERROR ) 
  {
    Cerr("***** iScannerPET::ComputeLUT() -> Error occurred when trying to read z-shift nb !" << endl);
    return 1;
  }
  // z-shift not found, or == 0
  else if( rvalue > 1 || rsector_nb_zshift==0)
  {
    rsector_nb_zshift = 1; // set to default value (=1)
    zshift = new FLTNBLUT[rsector_nb_zshift];
    zshift[0] = 0;
  }
  // (==0, zhift value provided)
  else
  {
    zshift = new FLTNBLUT[rsector_nb_zshift];

    if(ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "rsectors ZShift", zshift, rsector_nb_zshift, KEYWORD_MANDATORY) )
    {
      Cerr("***** iScannerPET::ComputeLUT() -> No data found about modules z-shift in the scanner header file, whereas z-shift is enabled !" << endl);
      return 1;
    }
  }
  
  // keep track of the number of crystals already indexed between layer 
  int nb_crystals_cur = 0;

  // Loop on layer, crystals indexed according to the layer they belong to
  for(int lyr=0 ; lyr<m_nbLayers ; lyr++)
  {
    if (lyr>0)
      nb_crystals_cur+=mp_nbCrystalsInLayer[lyr-1];

    int nb_ang_rsector = nb_ang_rsector_lyr[lyr],
        nb_axl_rsector = nb_axl_rsector_lyr[lyr],
        nb_trs_mod = nb_trs_mod_lyr[lyr],
        nb_axl_mod = nb_axl_mod_lyr[lyr],
        nb_trs_submod = nb_trs_submod_lyr[lyr],
        nb_axl_submod = nb_axl_submod_lyr[lyr],
        nb_trs_crystal = nb_trs_crystal_lyr[lyr],
        nb_axl_crystal = nb_axl_crystal_lyr[lyr];

        FLTNBLUT radius              = radius_lyr[lyr],
                 rsector_first_angle = rsector_first_angle_lyr[lyr],
                 angular_span        = rsector_angular_span_lyr[lyr],
                 gap_axl_rsector     = gap_axl_rsector_lyr[lyr],
                 gap_trs_mod       = gap_trs_mod_lyr[lyr],
                 gap_axl_mod       = gap_axl_mod_lyr[lyr],
                 gap_trs_submod    = gap_trs_submod_lyr[lyr],
                 gap_axl_submod    = gap_axl_submod_lyr[lyr],
                 gap_trs_crystal   = gap_trs_crystal_lyr[lyr],
                 gap_axl_crystal   = gap_axl_crystal_lyr[lyr]; 

    FLTNBLUT size_trs_submod = nb_trs_crystal*mp_sizeCrystalTrans[lyr] + (nb_trs_crystal-1)*gap_trs_crystal;
    FLTNBLUT size_axl_submod = nb_axl_crystal*mp_sizeCrystalAxial[lyr] + (nb_axl_crystal-1)*gap_axl_crystal;
    FLTNBLUT size_trs_mod    = nb_trs_submod*size_trs_submod + (nb_trs_submod-1)*gap_trs_submod;
    FLTNBLUT size_axl_mod    = nb_axl_submod*size_axl_submod + (nb_axl_submod-1)*gap_axl_submod;
    FLTNBLUT size_trs_rsector = nb_trs_mod*size_trs_mod + (nb_trs_mod-1)*gap_trs_mod;
    FLTNBLUT size_axl_rsector = nb_axl_mod*size_axl_mod + (nb_axl_mod-1)*gap_axl_mod;


    int nb_rsectors = nb_ang_rsector * nb_axl_rsector;
    int nb_mod = nb_axl_mod * nb_trs_mod;
    int nb_submod = nb_axl_submod * nb_trs_submod;
    int nb_crystal = nb_axl_crystal * nb_trs_crystal;

    // Check the number of crystals < nb elements provided by the scanner configuration file
    int nb_crystals = nb_rsectors*nb_mod*nb_submod*nb_crystal + nb_crystals_cur;
    if(m_nbCrystals<nb_crystals)
    {
      Cerr("***** iScannerPET::ComputeLUT() -> Computed number of crystals computed from the geom file ("<< nb_crystals
         <<") > not consistent with the total number of crystals (including potential layers) provided in the geom file ("<< m_nbCrystals <<") !" << endl);
      return 1;
    }
    
    mp_nbCrystalsInLayer[lyr] = nb_rsectors*nb_mod*nb_submod*nb_crystal;
    int number_crystals_in_ring = mp_nbCrystalsInLayer[lyr]/p_nbRings[lyr];
    
    // ============================================================================================================
    // Main part of the function: Generate the LUT
    // ============================================================================================================
  
    if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Generate positions for each crystal"<< endl); 
    
    // loop to nb_ang_rsector+1. crystal_center[0] will be used to gather position of the reference rsector angular position (directly above isocenter)
    oMatrix ******crystal_center = new oMatrix *****[nb_ang_rsector+1]; 
  
    for(int rsa = 0; rsa < nb_ang_rsector+1 ; rsa++) 
    {
      crystal_center[rsa] = new oMatrix ****[ nb_axl_rsector ];

      for (int i = 0; i<nb_axl_rsector ; i++)
      {
        crystal_center[rsa][i] = new oMatrix ***[ nb_mod ];
        
        for (int j = 0; j<nb_mod ; j++)
        {
          crystal_center[rsa][i][j] = new oMatrix **[ nb_submod ];
    
          for (int k = 0; k<nb_submod; k++)
          {
            crystal_center[rsa][i][j][k] = new oMatrix*[ nb_crystal ];
    
            for (int l = 0; l<nb_crystal; l++)
            {
              crystal_center[rsa][i][j][k][l]  = new oMatrix(3,1);
            }
          }
        }
      }
    }



    // Generation of the rotation matrix allowing to compute the position of all the rsectors. 
    oMatrix** rotation_mtx = new oMatrix*[nb_rsectors];
  
    for(int i=0; i<nb_ang_rsector; i++)
      rotation_mtx[i] = new oMatrix(3,3);
    
    // convert first angle and angular span in radians
    FLTNBLUT rsector_first_angle_rad = rsector_first_angle*M_PI/180.;  
    FLTNBLUT angular_span_rad = angular_span*M_PI/180.;


    int dir = (m_rotDirection == GEO_ROT_CCW) ? -1 : 1;
    
    
    for (int i = 0; i<nb_ang_rsector; i++)
    { 
      FLTNBLUT angle = remainderf(rsector_first_angle_rad + ((FLTNBLUT)i)*angular_span_rad/((FLTNBLUT)(nb_ang_rsector)), 2.*M_PI);
      
      rotation_mtx[i]->SetMatriceElt(0,0,cos(angle) );
      rotation_mtx[i]->SetMatriceElt(1,0,-dir*sin(angle) );
      rotation_mtx[i]->SetMatriceElt(2,0,0);
      rotation_mtx[i]->SetMatriceElt(0,1,dir*sin(angle) );
      rotation_mtx[i]->SetMatriceElt(1,1,cos(angle) );
      rotation_mtx[i]->SetMatriceElt(2,1,0);
      rotation_mtx[i]->SetMatriceElt(0,2,0);
      rotation_mtx[i]->SetMatriceElt(1,2,0);
      rotation_mtx[i]->SetMatriceElt(2,2,1);
    }

    // Compute scanner elements positions for the first rsector
    for (int r=0; r < nb_axl_rsector ; r++)
    {
      // Define the transaxial and axial edge start positions for the global rsector
      FLTNBLUT x_start_r = (-dir*size_trs_rsector) / 2.;
      FLTNBLUT z_start_r = -(nb_axl_rsector*size_axl_rsector + (nb_axl_rsector-1)*gap_axl_rsector) / 2.;

      z_start_r += r * (size_axl_rsector + gap_axl_rsector);

      for (int i=0; i < nb_mod ; i++)
      {
        // Define the transaxial and axial edge start positions for the rsector
        FLTNBLUT x_start_m = x_start_r;//-dir*(nb_trs_mod*size_trs_mod + (nb_trs_mod-1)*gap_trs_mod) / 2.;
        FLTNBLUT z_start_m = z_start_r;//-(nb_axl_mod*size_axl_mod + (nb_axl_mod-1)*gap_axl_mod) / 2.;
    
        // Define the transaxial and axial edge start positions for the i-Module in the rsector. 
        // Enumeration starting with the transaxial modules.
        x_start_m += dir*(i%nb_trs_mod) *  (size_trs_mod + gap_trs_mod);
        z_start_m += int(i/nb_trs_mod) * (size_axl_mod + gap_axl_mod);
    
        for (int j=0 ; j < nb_submod ; j++)
        {
          FLTNBLUT x_start_sm = x_start_m;
          FLTNBLUT z_start_sm = z_start_m;
          
          x_start_sm += dir*(j%nb_trs_submod) *  (size_trs_submod + gap_trs_submod);
          z_start_sm += int(j/nb_trs_submod) * (size_axl_submod + gap_axl_submod);
  
          for (int k=0 ; k < nb_crystal ; k++) 
          {
            // Define the transaxial and axial center positions for the j-SubModule (crystal) i-Module of the rsector.
            // Enumeration starting with the transaxial submodules.
            FLTNBLUT Xcrist = x_start_sm + dir* ( (k%nb_trs_crystal) * (mp_sizeCrystalTrans[lyr] + gap_trs_crystal) + mp_sizeCrystalTrans[lyr]/2 );
            FLTNBLUT Ycrist = radius + mp_sizeCrystalDepth[lyr]/2;
            FLTNBLUT Zcrist = z_start_sm + int(k/nb_trs_crystal) * (mp_sizeCrystalAxial[lyr] + gap_axl_crystal) + mp_sizeCrystalAxial[lyr]/2;
            
            crystal_center[0][r][i][j][k]->SetMatriceElt(0,0,Xcrist);
            crystal_center[0][r][i][j][k]->SetMatriceElt(1,0,Ycrist);
            crystal_center[0][r][i][j][k]->SetMatriceElt(2,0,Zcrist); 
          }
        }
      }
    }


    // ============================================================================================================
    // Loop over all the other rsectors.
    // Positions of the scanner elements are progressively stored in the LUT file.
    // ============================================================================================================
    for (int rsa=0 ; rsa<nb_ang_rsector ; rsa++)
      for (int rs=0 ; rs<nb_axl_rsector ; rs++)
        for (int m=0 ; m<nb_mod ; m++)
          for (int sm=0 ; sm<nb_submod ; sm++)
            for (int c=0 ; c<nb_crystal ; c++)
            {          
              // crystal indexation
              int cryID =                    rs * nb_axl_mod     * nb_axl_submod  * nb_axl_crystal * number_crystals_in_ring // nb indexed crystals in the rings covered by the previous (axial) rsectors
                        + int(m/nb_trs_mod)     * nb_axl_submod  * nb_axl_crystal * number_crystals_in_ring // nb indexed crystals in the rings covered by the previous (axial) modules
                        + int(sm/nb_trs_submod) * nb_axl_crystal * number_crystals_in_ring // nb indexed crystals in the rings covered by the previous (axial) submodules
                        + int(c/nb_trs_crystal) * number_crystals_in_ring // nb indexed crystals in the rings covered by the previous (axial) crystals
                        +                   rsa * nb_trs_mod    * nb_trs_submod * nb_trs_crystal // nb indexed crystals in the previous (transaxial) rsectors
                        +  m%nb_trs_mod         * nb_trs_submod * nb_trs_crystal // nb indexed crystals in the previous (transaxial) modules
                        + sm%nb_trs_submod      * nb_trs_crystal // nb indexed crystals in the previous (transaxial) submodules
                        +  c%nb_trs_crystal // previous crystals in the submodule
                        + nb_crystals_cur; // nb indexed crystals in the previous layer(s)
  
  
              rotation_mtx[rsa]->Multiplication(crystal_center[0][rs][m][sm][c], crystal_center[rsa+1][rs][m][sm][c]);
              
              mp_crystalCentralPositionX[cryID]  = (FLTNB)crystal_center[rsa+1][rs][m][sm][c]->GetMatriceElt(0,0);
              mp_crystalCentralPositionY[cryID]  = (FLTNB)crystal_center[rsa+1][rs][m][sm][c]->GetMatriceElt(1,0);
              mp_crystalCentralPositionZ[cryID]  = (FLTNB)crystal_center[rsa+1][rs][m][sm][c]->GetMatriceElt(2,0);
              mp_crystalCentralPositionZ[cryID] += (FLTNB)zshift[rsa%rsector_nb_zshift];
    
              // TODO Rotation for non-perfectly cylindric scanner (angles along the Z axis)
              mp_crystalOrientationX[cryID] = (FLTNB)rotation_mtx[rsa]->GetMatriceElt(0,1);
              mp_crystalOrientationY[cryID] = (FLTNB)rotation_mtx[rsa]->GetMatriceElt(0,0);
              mp_crystalOrientationZ[cryID] = 0.;
              
            }
          
    // ============================================================================================================
    // End
    // ============================================================================================================

    // Delete objects
    for (int rsa = 0; rsa<nb_ang_rsector+1 ; rsa++)
     for (int i = 0; i<nb_axl_rsector ; i++)
      for (int j = 0; j<nb_mod; j++)
       for (int k = 0; k<nb_submod; k++)
        for (int l = 0; l<nb_crystal; l++)
        {
         delete crystal_center[rsa][i][j][k][l];
        }


    for(int rsa = 0; rsa < nb_ang_rsector+1 ; rsa++)
     for(int i = 0; i < nb_axl_rsector ; i++)
      for (int j = 0; j<nb_mod; j++)
       for (int k = 0; k<nb_submod; k++)
       {
        delete[] crystal_center[rsa][i][j][k];
       }

    for(int rsa = 0; rsa < nb_ang_rsector+1 ; rsa++)
     for(int i = 0; i < nb_axl_rsector ; i++)
      for (int j = 0; j<nb_mod; j++)
      {
       delete[] crystal_center[rsa][i][j];
      }

    for(int rsa = 0; rsa < nb_ang_rsector+1 ; rsa++)
      for(int i = 0; i < nb_axl_rsector ; i++)
      {
       delete[] crystal_center[rsa][i];
      }

    for(int rsa = 0; rsa < nb_ang_rsector+1 ; rsa++)
      delete[] crystal_center[rsa];
            
            
    for(int i = 0; i < nb_ang_rsector ; i++)
      delete rotation_mtx[i];
    
    delete[] crystal_center;
    delete[] rotation_mtx;
    
  }

  // ============================================================================================================
  // Save LUT if required
  // ============================================================================================================

  if (p_scannerManager->SaveLUTFlag())
  {
    // Make file names
    string path_to_geom_file = p_scannerManager->GetPathToScannerFile();
    string path_to_LUT = path_to_geom_file.substr(0, path_to_geom_file.find_last_of("."));
    string path_to_header_LUT = path_to_LUT + ".ghscan";
    path_to_LUT.append(".glut");
    // Verbose
    if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Save LUT into file with extension '.glut'" << endl);
    // Open file
    ofstream LUT_file, header_LUT_file;
    LUT_file.open(path_to_LUT.c_str(), ios::binary | ios::out);
    // Write the corresponding crystal parameters in the binary LUT
    for (int i=0 ; i<m_nbCrystals ; i++)
    {
      LUT_file.write(reinterpret_cast<char*>(&mp_crystalCentralPositionX[i]), sizeof(FLTNB));
      LUT_file.write(reinterpret_cast<char*>(&mp_crystalCentralPositionY[i]), sizeof(FLTNB));
      LUT_file.write(reinterpret_cast<char*>(&mp_crystalCentralPositionZ[i]), sizeof(FLTNB));
      LUT_file.write(reinterpret_cast<char*>(&mp_crystalOrientationX[i]), sizeof(FLTNB));
      LUT_file.write(reinterpret_cast<char*>(&mp_crystalOrientationY[i]), sizeof(FLTNB));
      LUT_file.write(reinterpret_cast<char*>(&mp_crystalOrientationZ[i]), sizeof(FLTNB));
    }
    // Close file
    LUT_file.close(); 
    // ---------------------
    // --> Write header file
    // ---------------------
    if (m_verbose>=VERBOSE_DETAIL) Cout("  -->  Save header LUT file with extension '.ghscan'" << endl);
    // Open file
    header_LUT_file.open(path_to_header_LUT.c_str(), ios::out); 
    // Start writing
    string scanner_name = path_to_geom_file.substr(0, path_to_geom_file.find_last_of("."));
    header_LUT_file << "scanner name:" << "    " << GetFileFromPath(scanner_name) << endl; 
    header_LUT_file << "modality:" << "    " << "PET" << endl; 
    
    header_LUT_file << "scanner radius:" << "    " << radius_lyr[0];
    for (int lyr=1 ; lyr<m_nbLayers ; lyr++) 
      header_LUT_file << "," << radius_lyr[lyr] ; header_LUT_file << endl;
      
    header_LUT_file << "number of elements:" << "    " << m_nbCrystals << endl; 
    header_LUT_file << "number of layers:" << "    " << m_nbLayers << endl;
    header_LUT_file << "number of crystals in layer(s):" << "    " << mp_nbCrystalsInLayer[0]; // TODO nb_cry_in_layer
    for (int lyr=1 ; lyr<m_nbLayers ; lyr++) 
      header_LUT_file << ","<< mp_nbCrystalsInLayer[lyr] ; header_LUT_file << endl;
      
    header_LUT_file << "crystals size depth:" << "    " << mp_sizeCrystalDepth[0];
    for (int lyr=1 ; lyr<m_nbLayers ; lyr++) 
      header_LUT_file << ","<< mp_sizeCrystalDepth[lyr] ; header_LUT_file << endl;
      
    header_LUT_file << "crystals size transaxial:" << "    " << mp_sizeCrystalTrans[0];  
    for (int lyr=1 ; lyr<m_nbLayers ; lyr++) 
      header_LUT_file << ","<< mp_sizeCrystalTrans[lyr] ; header_LUT_file << endl;
      
    header_LUT_file << "crystals size axial:" << "    " << mp_sizeCrystalAxial[0]; 
    for (int lyr=1 ; lyr<m_nbLayers ; lyr++) 
      header_LUT_file << ","<< mp_sizeCrystalAxial[lyr] ; header_LUT_file << endl;

    // Default reconstruction parameters
    uint32_t def_dim_trs = 0, def_dim_axl = 0; // default number of voxels
    FLTNB def_FOV_trs = -1., def_FOV_axl = -1; // computed from modules width and gaps of the 1st layer
    
    if( ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "voxels number transaxial", &def_dim_trs, 1, KEYWORD_MANDATORY) ||
        ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "voxels number axial", &def_dim_axl, 1, KEYWORD_MANDATORY) ||
        ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "field of view transaxial", &def_FOV_trs, 1, KEYWORD_MANDATORY) ||
        ReadDataASCIIFile(p_scannerManager->GetPathToScannerFile(), "field of view axial", &def_FOV_axl, 1, KEYWORD_MANDATORY) )
        {
          Cerr("***** iScannerPET::ComputeLUT() -> Error occurred when trying to read transaxial/axial dimensions and voxel sizes from scanner geom file !" << endl);
          return 1;
        }
    
    header_LUT_file << "voxels number transaxial:" << "    " << def_dim_trs << endl;  
    header_LUT_file << "voxels number axial:" << "    " << def_dim_axl << endl; 
  
    header_LUT_file << "field of view transaxial:" << "    " << def_FOV_trs << endl;  
    header_LUT_file << "field of view axial:" << "    " << def_FOV_axl << endl; 
    
    header_LUT_file << "min angle difference:" << "    " << m_minAngleDifference << " #deg" << endl;
    
    header_LUT_file << "mean depth of interaction:" << "    " << mp_meanDepthOfInteraction[0];
    for (int lyr=1 ; lyr<m_nbLayers ; lyr++) 
      header_LUT_file << "," << mp_meanDepthOfInteraction[lyr] ;
    header_LUT_file << " #optional (default value : center of crystal ). Input value must correspond to the distance from the crystal surface, or negative value if default" << endl;
    
    header_LUT_file << "description:" << "    " << "LUT generated from geom file: " << GetFileFromPath(p_scannerManager->GetPathToScannerFile()) << endl;
    
    if(m_verbose>=2) Cout("iScannerPET::ComputeLUT() -> Header LUT file writing completed" << endl);
  }

  delete[] p_nbRings;
  delete[] nb_ang_rsector_lyr;
  delete[] nb_axl_rsector_lyr;
  delete[] nb_trs_mod_lyr;
  delete[] nb_axl_mod_lyr;
  delete[] nb_trs_submod_lyr;
  delete[] nb_axl_submod_lyr;
  delete[] nb_trs_crystal_lyr;
  delete[] nb_axl_crystal_lyr;

  delete[] radius_lyr;
  delete[] rsector_angular_span_lyr;
  delete[] rsector_first_angle_lyr;
  delete[] zshift;
  delete[] gap_axl_rsector_lyr;
  delete[] gap_trs_mod_lyr;
  delete[] gap_axl_mod_lyr;
  delete[] gap_trs_submod_lyr;
  delete[] gap_axl_submod_lyr;
  delete[] gap_trs_crystal_lyr;
  delete[] gap_axl_crystal_lyr;

  if (m_verbose>=VERBOSE_DETAIL) Cout("  --> LUT generation completed" << endl);

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetPositionsAndOrientations
  \param a_index1 : index of the 1st crystal 
  \param a_index2 : index of the 2nd crystal 
  \param ap_Position1[3] : x,y,z cartesian position of center of the 1st crystal
  \param ap_Position2[3] : x,y,z cartesian position of center of the 2nd crystal
  \param ap_Orientation1[3] : x,y,z components of the orientation vector related to the 1st crystal
  \param ap_Orientation2[3] : x,y,z components of the orientation vector related to the 2nd crystal
  \param ap_POI1 : x,y,z components of the Point Of Interation related to the 1st crystal
  \param ap_POI2 : x,y,z components of the Point Of Interation related to the 2nd crystal
  \brief Get the central positions and orientations of the scanner elements from their indices.
  \details This method is very general and is used by the vProjector. 
           From these positions and orientations, other methods can be used by specific projectors to get specific positions.
           Position calculations include POI and mean depth of interaction
  \todo some cases depending on POI are not implemented
  \return 0 if success, positive value otherwise
*/
int iScannerPET::GetPositionsAndOrientations( int a_index1, int a_index2,
                                               FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                               FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3],
                                               FLTNB* ap_POI1, FLTNB* ap_POI2 )
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  // First check crystals existency
  if (a_index1<0 || a_index1>=m_nbCrystals)
  {
    Cerr("***** iScannerPET::GetPositionsAndOrientations() -> Crystal index 1 (" << a_index1 << ") out of range [0:" << m_nbCrystals-1 << "] !" << endl);
    return 1;
  }
  if (a_index2<0 || a_index2>=m_nbCrystals)
  {
    Cerr("***** iScannerPET::GetPositionsAndOrientations() -> Crystal index 2 (" << a_index2 << ") out of range [0:" << m_nbCrystals-1 << "] !" << endl);
    return 1;
  }

  // The POI specification is as follows:
  //  - POI[0] and POI[1] express the local x and y coordinates: x is the tangential axis, and y is colinear to the global z axis.
  //    Their origins are in the middle of the scanner elements so they range from minus half to half the element dimensions.
  //  - POI[2] express the DOI along the local z axis which is colinear to the crystal orientation (which itself is oriented from
  //    the center of the scanner to the exterior).
  //    Its origin is at the entrance of the element and thus ranges from 0 to the local z size of the element (crystal depth).
  //  - If POI[2] has a negative value, this means that the mean depth of interaction should be considered instead, while still
  //    using POI[0] and POI[1].
  // Remember here that the crystal position is at its center.

  // Case when POI is not provided (we use the mean depth of interaction)
  if (ap_POI1==NULL)
  { 
    FLTNB depth = mp_meanDepthOfInteraction[GetLayer(a_index1)] - mp_sizeCrystalDepth[GetLayer(a_index1)]/2;
    ap_Position1[0] = mp_crystalCentralPositionX[a_index1] + depth*mp_crystalOrientationX[a_index1];
    ap_Position1[1] = mp_crystalCentralPositionY[a_index1] + depth*mp_crystalOrientationY[a_index1];
    ap_Position1[2] = mp_crystalCentralPositionZ[a_index1] + depth*mp_crystalOrientationZ[a_index1];
  }
  // Case when POI[2] is negative (meaning we only have POI[0] or POI[1] specified and to be taken into account)
  else if (ap_POI1[2]<0.)
  {
    // TODO
  }
  // Case when only the DOI is provided
  else if (ap_POI1[0]==0. && ap_POI1[1]==0.)
  {
    FLTNB depth = ap_POI1[2] - mp_sizeCrystalDepth[GetLayer(a_index1)]/2;
    ap_Position1[0] = mp_crystalCentralPositionX[a_index1] + depth*mp_crystalOrientationX[a_index1];
    ap_Position1[1] = mp_crystalCentralPositionY[a_index1] + depth*mp_crystalOrientationY[a_index1];
    ap_Position1[2] = mp_crystalCentralPositionZ[a_index1] + depth*mp_crystalOrientationZ[a_index1];
  }
  // Case when the full POI is taken into account
  else
  {
    // TODO
  }

  // Case when POI is not provided (we use the mean depth of interaction)
  if (ap_POI2==NULL)
  {
    FLTNB depth = mp_meanDepthOfInteraction[GetLayer(a_index2)] - mp_sizeCrystalDepth[GetLayer(a_index2)]/2;
    ap_Position2[0] = mp_crystalCentralPositionX[a_index2] + depth*mp_crystalOrientationX[a_index2];
    ap_Position2[1] = mp_crystalCentralPositionY[a_index2] + depth*mp_crystalOrientationY[a_index2];
    ap_Position2[2] = mp_crystalCentralPositionZ[a_index2] + depth*mp_crystalOrientationZ[a_index2];
  }
  // Case when POI[2] is negative (meaning we only have POI[0] or POI[1] specified and to be taken into account)
  else if (ap_POI2[2]<0.)
  {
    // TODO
  }
  // Case when only the DOI is provided
  else if (ap_POI2[0]==0. && ap_POI2[1]==0.)
  {
    FLTNB depth = ap_POI2[2] - mp_sizeCrystalDepth[GetLayer(a_index2)]/2;
    ap_Position2[0] = mp_crystalCentralPositionX[a_index2] + depth*mp_crystalOrientationX[a_index2];
    ap_Position2[1] = mp_crystalCentralPositionY[a_index2] + depth*mp_crystalOrientationY[a_index2];
    ap_Position2[2] = mp_crystalCentralPositionZ[a_index2] + depth*mp_crystalOrientationZ[a_index2];
  }
  // Case when the full POI is taken into account
  else
  {
    // TODO
  }

  // Get orientations
  ap_Orientation1[0] = mp_crystalOrientationX[a_index1];
  ap_Orientation1[1] = mp_crystalOrientationY[a_index1];
  ap_Orientation1[2] = mp_crystalOrientationZ[a_index1];
  ap_Orientation2[0] = mp_crystalOrientationX[a_index2];
  ap_Orientation2[1] = mp_crystalOrientationY[a_index2];
  ap_Orientation2[2] = mp_crystalOrientationZ[a_index2];
  
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetRdmPositionsAndOrientations
  \param a_index1 : index of the 1st crystal 
  \param a_index2 : index of the 2nd crystal 
  \param ap_Position1[3] : x,y,z cartesian position of center of the 1st crystal
  \param ap_Position2[3] : x,y,z cartesian position of center of the 2nd crystal
  \param ap_Orientation1[3] : x,y,z components of the orientation vector related to the 1st crystal
  \param ap_Orientation2[3] : x,y,z components of the orientation vector related to the 2nd crystal
  \brief Get random positions of the scanner elements and their orientations from their indices.
  \details - Find the scanner elements described by the two indexes passed in parameters. 
  \details - Compute random positions inside the elements described by the indexes passed in parameters
  \details - Find the scanner elements described by the two indexes passed in parameters.
  \details - Write the corresponding random cartesian coordinates in the positions parameters.
  \details Position calculations include POI and mean depth of interaction
  \todo fix the possibility to draw LOR outside the actual crystal volume (if mp_meanDepthOfInteraction != 0.5)
  \todo some cases depending on POI are not implemented
  \return 0 if success, positive value otherwise
*/
int iScannerPET::GetRdmPositionsAndOrientations(int a_index1, int a_index2,
                                                FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                                FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3] )
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  // First check crystals existency
  if (a_index1<0 || a_index1>=m_nbCrystals)
  {
    Cerr("***** iScannerPET::GetRdmPositionsAndOrientations() -> Crystal index 1 (" << a_index1 << ") out of range [0:" << m_nbCrystals-1 << "] !" << endl);
    return 1;
  }
  if (a_index2<0 || a_index2>=m_nbCrystals)
  {
    Cerr("***** iScannerPET::GetRdmPositionsAndOrientations() -> Crystal index 2 (" << a_index2 << ") out of range [0:" << m_nbCrystals-1 << "] !" << endl);
    return 1;
  }

  if (mp_sizeCrystalTrans[GetLayer(a_index1)]<=0. || mp_sizeCrystalAxial[GetLayer(a_index1)]<=0. ||
      mp_sizeCrystalTrans[GetLayer(a_index2)]<=0. || mp_sizeCrystalAxial[GetLayer(a_index2)]<=0. )
  {
    Cerr("***** iScannerPET::GetRdmPositionsAndOrientations() -> Crystal sizes are unknown or equal to 0. Crystal dimensions are mandatory for this function !" << endl);
    return 1;
  }

  sRandomNumberGenerator* p_RNG = sRandomNumberGenerator::GetInstance(); 

  // Get random numbers for the first crystal
  FLTNB rdm_axl1 = p_RNG->GenerateRdmNber();
  FLTNB rdm_trs1 = p_RNG->GenerateRdmNber();
  FLTNB rdm_depth1 = p_RNG->GenerateRdmNber();

  FLTNB size_crystalTrans1 = mp_sizeCrystalTrans[GetLayer(a_index1)];
  FLTNB size_crystalAxial1 = mp_sizeCrystalAxial[GetLayer(a_index1)];
  FLTNB size_crystalDepth1 = mp_sizeCrystalDepth[GetLayer(a_index1)];

  FLTNB axial1 = (rdm_axl1-0.5) * size_crystalAxial1;
  FLTNB trans1 = (rdm_trs1-0.5) * size_crystalTrans1;
  FLTNB depth1 = (rdm_depth1-0.5) * size_crystalDepth1;

  ap_Position1[0] = mp_crystalCentralPositionX[a_index1] 
                  + depth1 * mp_crystalOrientationX[a_index1] 
                  + trans1 * mp_crystalOrientationY[a_index1] 
                  + axial1 * mp_crystalOrientationX[a_index1] * mp_crystalOrientationZ[a_index1];
                  
  ap_Position1[1] = mp_crystalCentralPositionY[a_index1] 
                  + depth1 * mp_crystalOrientationY[a_index1] 
                  + trans1 * mp_crystalOrientationX[a_index1] 
                  + axial1 * mp_crystalOrientationY[a_index1] * mp_crystalOrientationZ[a_index1];
                  
  ap_Position1[2] = mp_crystalCentralPositionZ[a_index1] 
                  + depth1 * mp_crystalOrientationZ[a_index1] 
                  + axial1 * sqrt(1-mp_crystalOrientationZ[a_index1] * mp_crystalOrientationZ[a_index1]);

  // Get orientations
  ap_Orientation1[0] = mp_crystalOrientationX[a_index1];
  ap_Orientation1[1] = mp_crystalOrientationY[a_index1];
  ap_Orientation1[2] = mp_crystalOrientationZ[a_index1];

  // Get random numbers for the second crystal
  FLTNB rdm_axl2 = p_RNG->GenerateRdmNber();
  FLTNB rdm_trs2 = p_RNG->GenerateRdmNber();
  FLTNB rdm_depth2 = p_RNG->GenerateRdmNber();

  FLTNB size_crystalTrans2 = mp_sizeCrystalTrans[GetLayer(a_index2)];
  FLTNB size_crystalAxial2 = mp_sizeCrystalAxial[GetLayer(a_index2)];
  FLTNB size_crystalDepth2 = mp_sizeCrystalDepth[GetLayer(a_index2)];

  FLTNB axial2 = (rdm_axl2-0.5) * size_crystalAxial2;
  FLTNB trans2 = (rdm_trs2-0.5) * size_crystalTrans2;
  FLTNB depth2 = (rdm_depth2-0.5) * size_crystalDepth2;

  ap_Position2[0] = mp_crystalCentralPositionX[a_index2] 
                  + depth2 * mp_crystalOrientationX[a_index2] 
                  + trans2 * mp_crystalOrientationY[a_index2] 
                  + axial2 * mp_crystalOrientationX[a_index2] * mp_crystalOrientationZ[a_index2];
                  
  ap_Position2[1] = mp_crystalCentralPositionY[a_index2] 
                  + depth2 * mp_crystalOrientationY[a_index2] 
                  + trans2 * mp_crystalOrientationX[a_index2] 
                  + axial2 * mp_crystalOrientationY[a_index2] * mp_crystalOrientationZ[a_index2];
                  
  ap_Position2[2] = mp_crystalCentralPositionZ[a_index2] 
                  + depth2 * mp_crystalOrientationZ[a_index2] 
                  + axial2 * sqrt(1-mp_crystalOrientationZ[a_index2] * mp_crystalOrientationZ[a_index2]);

  // Get orientations
  ap_Orientation2[0] = mp_crystalOrientationX[a_index2];
  ap_Orientation2[1] = mp_crystalOrientationY[a_index2];
  ap_Orientation2[2] = mp_crystalOrientationZ[a_index2];

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetPositionWithRandomDepth
  \param a_index1 : index of the 1st crystal
  \param a_index2 : index of the 2nd crystal
  \param ap_Position1[3] : x,y,z cartesian position of the point related to the 1st index (see child function for more details)
  \param ap_Position2[3] : x,y,z cartesian position of the point related to the 2st index (see child function for more details)
  \brief Get the positions and orientations of scanner elements from their indices, with a random depth.
  \details Method for testing purposes. Does not include POI and mean depth of interaction
  \return 0 if success, positive value otherwise
*/
int iScannerPET::GetPositionWithRandomDepth(int a_index1, int a_index2, FLTNB ap_Position1[3], FLTNB ap_Position2[3])
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  // Get instance of random number generator
  sRandomNumberGenerator* p_RNG = sRandomNumberGenerator::GetInstance(); 
  
  // Shoot first random number
  FLTNB shoot1 = p_RNG->GenerateRdmNber();
  // Compute associated depth
  FLTNB depth1 =  (shoot1-0.5) * mp_sizeCrystalDepth[GetLayer(a_index1)];
  
  // Compute associated position
  ap_Position1[0] = mp_crystalCentralPositionX[a_index1] + depth1*mp_crystalOrientationX[a_index1];
  ap_Position1[1] = mp_crystalCentralPositionY[a_index1] + depth1*mp_crystalOrientationY[a_index1];
  ap_Position1[2] = mp_crystalCentralPositionZ[a_index1] + depth1*mp_crystalOrientationZ[a_index1];
  
  // Shoot second random number
  FLTNB shoot2 = p_RNG->GenerateRdmNber();
  // Compute associated depth
  FLTNB depth2 =  (shoot2-0.5) * mp_sizeCrystalDepth[GetLayer(a_index2)];
  // Compute associated position
  ap_Position2[0] = mp_crystalCentralPositionX[a_index2] + depth2*mp_crystalOrientationX[a_index2];
  ap_Position2[1] = mp_crystalCentralPositionY[a_index2] + depth2*mp_crystalOrientationY[a_index2];
  ap_Position2[2] = mp_crystalCentralPositionZ[a_index2] + depth2*mp_crystalOrientationZ[a_index2];
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetTwoCorners
  \param a_index1 : index of the 1st crystal 
  \param a_index2 : index of the 2nd crystal
  \param ap_CornerInf1[3]
  \param ap_CornerSup1[3]
  \param ap_CornerInf2[3]
  \param ap_CornerSup2[3]
  \brief Get the cartesian coordinaters of the two opposite corners of a scanner element.
  \todo Not implemented yet 
  \return 0 if success, positive value otherwise
*/
int iScannerPET::GetTwoCorners(int a_index1, int a_index2,
                               FLTNB ap_CornerInf1[3], FLTNB ap_CornerSup1[3],
                               FLTNB ap_CornerInf2[3], FLTNB ap_CornerSup2[3])
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  if (a_index1<0 || a_index1>=m_nbCrystals || a_index2<0 || a_index2>=m_nbCrystals)
  {
    Cerr("***** iScannerPET::GetTwoCorners() -> Crystal index out of range !" << endl);
    return 1;
  }

  if (mp_sizeCrystalTrans[GetLayer(a_index1)]<=0. || mp_sizeCrystalAxial[GetLayer(a_index1)]<=0. ||
      mp_sizeCrystalTrans[GetLayer(a_index2)]<=0. || mp_sizeCrystalAxial[GetLayer(a_index2)]<=0. )
  {
    Cerr("***** iScannerPET::GetRdmPositionsAndOrientations() -> Crystal sizes are unknown or equal to 0. Crystal dimensions are mandatory for this function !" << endl);
    return 1;
  }
  
  Cerr("***** iScannerPET::GetTwoCorners() -> Not implemented yet !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iScannerPET::GetEdgesCenterPositions( int a_index1, int a_index2,
                                          FLTNB ap_pos_line_point1[3], FLTNB ap_pos_line_point2[3],
                                          FLTNB ap_pos_point1_x[4], FLTNB ap_pos_point1_y[4], FLTNB ap_pos_point1_z[4],
                                          FLTNB ap_pos_point2_x[4], FLTNB ap_pos_point2_y[4], FLTNB ap_pos_point2_z[4]
)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  if( a_index1 < 0 || a_index1 >= m_nbCrystals
   || a_index2 < 0 || a_index2 >= m_nbCrystals )
  {
    Cerr("***** iScannerPET::GetEdgesCenterPositions() -> Crystal indices out of range !" << endl);
    return 1;
  }

  // Getting the half size of axial/trans crystal 1 and 2 depending on the layer
  FLTNB half_crystal_trans_1 = mp_sizeCrystalTrans[GetLayer(a_index1)] / 2.0;
  FLTNB half_crystal_axial_1 = mp_sizeCrystalAxial[GetLayer(a_index1)] / 2.0;
  FLTNB half_crystal_trans_2 = mp_sizeCrystalTrans[GetLayer(a_index2)] / 2.0;
  FLTNB half_crystal_axial_2 = mp_sizeCrystalAxial[GetLayer(a_index2)] / 2.0;

  //////////////////////////////////////////////////////////////////////////
  // Computing coordinates depending on the orientation for point1
  // X-axis
  ap_pos_point1_x[ 0 ] = ap_pos_line_point1[ 0 ] - half_crystal_trans_1 * mp_crystalOrientationY[ a_index1 ];
  ap_pos_point1_x[ 1 ] = ap_pos_line_point1[ 0 ] + half_crystal_trans_1 * mp_crystalOrientationY[ a_index1 ];
  ap_pos_point1_x[ 2 ] = ap_pos_line_point1[ 0 ];
  ap_pos_point1_x[ 3 ] = ap_pos_line_point1[ 0 ];

  // Y-axis
  ap_pos_point1_y[ 0 ] = ap_pos_line_point1[ 1 ] + half_crystal_trans_1 * mp_crystalOrientationX[ a_index1 ];
  ap_pos_point1_y[ 1 ] = ap_pos_line_point1[ 1 ] - half_crystal_trans_1 * mp_crystalOrientationX[ a_index1 ];
  ap_pos_point1_y[ 2 ] = ap_pos_line_point1[ 1 ];
  ap_pos_point1_y[ 3 ] = ap_pos_line_point1[ 1 ];

  // Z-axis
  ap_pos_point1_z[ 0 ] = ap_pos_line_point1[ 2 ];
  ap_pos_point1_z[ 1 ] = ap_pos_line_point1[ 2 ];
  ap_pos_point1_z[ 2 ] = ap_pos_line_point1[ 2 ] - half_crystal_axial_1;
  ap_pos_point1_z[ 3 ] = ap_pos_line_point1[ 2 ] + half_crystal_axial_1;

  //////////////////////////////////////////////////////////////////////////
  // Computing coordinates depending on the orientation for point2
  // X-axis
  ap_pos_point2_x[ 0 ] = ap_pos_line_point2[ 0 ] + half_crystal_trans_2 * mp_crystalOrientationY[ a_index2 ];
  ap_pos_point2_x[ 1 ] = ap_pos_line_point2[ 0 ] - half_crystal_trans_2 * mp_crystalOrientationY[ a_index2 ];
  ap_pos_point2_x[ 2 ] = ap_pos_line_point2[ 0 ];
  ap_pos_point2_x[ 3 ] = ap_pos_line_point2[ 0 ];

  // Y-axis
  ap_pos_point2_y[ 0 ] = ap_pos_line_point2[ 1 ] - half_crystal_trans_2 * mp_crystalOrientationX[ a_index2 ];
  ap_pos_point2_y[ 1 ] = ap_pos_line_point2[ 1 ] + half_crystal_trans_2 * mp_crystalOrientationX[ a_index2 ];
  ap_pos_point2_y[ 2 ] = ap_pos_line_point2[ 1 ];
  ap_pos_point2_y[ 3 ] = ap_pos_line_point2[ 1 ];

  // Z-axis
  ap_pos_point2_z[ 0 ] = ap_pos_line_point2[ 2 ];
  ap_pos_point2_z[ 1 ] = ap_pos_line_point2[ 2 ];
  ap_pos_point2_z[ 2 ] = ap_pos_line_point2[ 2 ] - half_crystal_axial_2;
  ap_pos_point2_z[ 3 ] = ap_pos_line_point2[ 2 ] + half_crystal_axial_2;

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetLayer
  \param a_idx : index of the crystal in the loaded LUT
  \brief Get the layer from which the 'a_idx' crystal belongs to
  \return layer index
*/
int iScannerPET::GetLayer(int a_idx)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_MAX)
  
  int layer=0;
  int sum_crystals = mp_nbCrystalsInLayer[layer];
  // Get the layer which the crystal belongs to
  while (a_idx >= sum_crystals)
  {
    layer++;
    sum_crystals += mp_nbCrystalsInLayer[layer];
  }
  return layer;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IsAvailableLOR
  \param a_elt1 : index of the 1st crystal
  \param a_elt2 : index of the 2nd crystal
  \brief Check if the LOR formed by the crystalf whose indices are passed 
         in parameters is available according to the scanner restrictions
  \details This function is dedicated to analytic projection and list-mode sensitivity image generation
           The PET implementation checks the LOR availability according to the minimal (transaxial) angle difference 
           and the maximal ring difference between two crystals
  \todo min angle difference (currently just system dependent) may be estimated from the FOV requested by the user ?
  \todo Restriction for ring_ID NOT safe for user using their own LUT !!!
        Perhaps throw warning in this case
  \return 1 if the LOR is available, 0 otherwise
*/
int iScannerPET::IsAvailableLOR(int a_elt1, int a_elt2)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  // Get absolute angle between the two crystals
  FLTNB abs_angle = mp_crystalOrientationX[a_elt1]*mp_crystalOrientationX[a_elt2] 
                  + mp_crystalOrientationY[a_elt1]*mp_crystalOrientationY[a_elt2];

  // Handle boundaries (arg of acos() outside range [-1. ; +1]) 
  abs_angle = (abs_angle>1.) ? 1 : abs_angle;
  abs_angle = (abs_angle<-1.) ? -1 : abs_angle;
   
  FLTNB angle_diff = acos(abs_angle);
 
  // Check restrictions (Axial constraint in mm and Transaxial constraint on min angle difference)
  // TODO : min angle difference (currently just system dependent) may be estimated from the FOV requested by the user ?  
  if ( ( m_maxAxialDiffmm <0. || abs(mp_crystalCentralPositionZ[a_elt1]-mp_crystalCentralPositionZ[a_elt2]) <= m_maxAxialDiffmm ) // a max axial diff restriction
    && ( angle_diff>=m_minAngleDifference || FLTNBIsEqual(angle_diff,m_minAngleDifference,.0001)) ) // Returns 1 if computed angle diff ~ min angle diff (E-4)
  {
    return 1;
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetGeometricInfoFromDataFile
  \brief Retrieve PET geometric informations from the datafile
  \details Recover the (axial) max ring difference
  \return 0 if success, positive value otherwise
*/
int iScannerPET::GetGeometricInfoFromDataFile(string a_pathToDF)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_DETAIL) Cout("iScannerPET::GetGeometricInfoFromDataFile() -> Specific to PET" << endl);
  
  // This function is intended to be called after the scanner initialization, so that any field present in the datafile header, similar to
  // one in the scanner configuration file, may overload the scanner value.

  // Get maximum axial difference in mm
  if (ReadDataASCIIFile(a_pathToDF, "Maximum axial difference mm", &m_maxAxialDiffmm, 1, KEYWORD_OPTIONAL)==1)
  {
    Cerr("***** iScannerPET::GetGeometricInfoFromDataFile() -> Error while reading max number of ring difference in the header data file '" << endl);
    return 1;
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetPETSpecificParameters
  \param ap_maxAxialDiffmm
  \brief Set pointers passed in argument with the related PET specific variables
         This function is used to recover these values in the datafile object
  \return 0 if success, positive value otherwise
*/
int iScannerPET::PROJ_GetPETSpecificParameters(FLTNB* ap_maxAxialDiffmm)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_NORMAL) Cout("iScannerPET::PROJ_GetPETSpecificParameters ..." << endl);
  
  // Verify that all parameters have been correctly checked
  if (!m_allParametersChecked)
  {
    Cerr("***** iScannerPET::PROJ_GetPETSpecificParameters() -> Parameters have not been checked !" << endl);
    return 1;
  }
  
  *ap_maxAxialDiffmm = m_maxAxialDiffmm;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowHelp
  \brief Display help
  \todo Provide informations about PET system initialization ?
*/
void iScannerPET::ShowHelp()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  cout << "This scanner class is dedicated to the description of PET systems." << endl;
}
