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
  \brief    Implementation of class vScanner
*/

#include "sScannerManager.hh"
#include "vScanner.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vScanner::vScanner() 
{
  // Initialize all members to default values
  m_scannerType = SCANNER_UNKNOWN;
  m_verbose = -1;
  mp_ID = NULL;
  m_allParametersChecked = false;
  mp_rotationMatrix = NULL; 
  mp_positionMatrix_ref = NULL;
  mp_positionMatrix_out = NULL;
  m_rotDirection = GEO_ROT_CW;
  m_defaultBedDisplacementInMm = 0.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vScanner::~vScanner() 
{
  if (mp_rotationMatrix) delete mp_rotationMatrix; 
  if (mp_positionMatrix_ref) delete mp_positionMatrix_ref;
  if (mp_positionMatrix_out) delete mp_positionMatrix_out;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int vScanner::Describe()
  \brief   A function used to describe the generic parts of the datafile
*/
void vScanner::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose==0) return;
  
  // Describe the datafile
  Cout("vScanner::Describe() -> Here is some generic content of the scanner" << endl);
  if (m_scannerType == 0)
  Cout("  --> Scanner type: " << GetScannerTypeString() << endl);

  // Call the specific function of the scanner
  DescribeSpecific();
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      inline int vScanner::GetScannerTypeString()
  \return  the string corresponding to the scanner type (modality) 
           as defined by macro SCANNER_TYPE in sScannerManager.hh
*/
string vScanner::GetScannerTypeString()
{
  if (m_scannerType == SCANNER_PET)
    return "SCANNER PET";
  else if (m_scannerType == SCANNER_SPECT_PINHOLE)
    return "SCANNER SPECT PINHOLE";
  else if (m_scannerType == SCANNER_SPECT_CONVERGENT)
    return "SCANNER SPECT CONVERGENT";
  else if (m_scannerType == SCANNER_CT)
    return "SCANNER CT";
  else if (m_scannerType == SCANNER_SINOGRAM)
    return "SCANNER SINOGRAM";
  
  // Default
  return "Unknown";
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ComputeLUT
  \brief Virtual function which should be implemented by the child classes.
         It computes the LUT of the scanner from a generic (.geom) file.
         The vScanner implementation throws error by default as it should be implemented by the child class
  \todo  Make the mother function virtual pure ? (should iScannerSinogram  also implement this function ?)
  \todo  iScannerCT implementation will probably consists in computing one projection, then computing
         the others on-the-fly during reconstruction (Often too much data to keep in memory for CT)
         Have to check if we offer the precomputation of the entire LUT in some situation, as for PET/SPECT
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::ComputeLUT() 
{
  Cerr("***** vScanner::ComputeLUT() -> Call to ComputeLUT() which is not implemented by the scanner child class !" << endl);
  return 1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LoadLUT
  \brief Virtual function which should be implemented by the child classes.
         Load a precomputed scanner LUT. 
         The vScanner implementation throws error by default as it should be implemented by the child class
  \todo  Make the mother function virtual pure ? (should iScannerSinogram  also implement this function ?)
  \todo  iScannerCT implementation will probably consists in computing one projection, then computing
         the others on-the-fly during reconstruction (Often too much data to keep in memory for CT)
         Have to check if we offer the precomputation of the entire LUT in some situation, as for PET/SPECT
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::LoadLUT() 
{
  Cerr("***** vScanner::ComputeLUT() -> Call to ComputeLUT() which is not implemented by the scanner child class !" << endl);
  return 1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
// ----- PET Specific Functions --- //
/*
  \fn IsAvailableLOR
  \param a_elt1 : index of the 1st scanner element
  \param a_elt2 : index of the 2nd scanner element
  \brief This function is implemented in child classes.
         Check if the LOR is available according to the scanner restrictions
  \details This function is related to analytic projection and list-mode sensitivity image generation
  \return 1 if the LOR is available, 0 otherwise (vScanner implementation returns 1 by default)
*/
int vScanner::IsAvailableLOR(int a_elt1, int a_elt2)
{
  Cerr("***** vScanner::IsAvailableLOR() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                    This function only works with PET scanner objects !!" << endl);
  return 1;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SetPETMaxRingDiff
  \param a_maxAxialDiffmm
  \brief Set the maximal axial difference in mm between 2 crystals forming a lor
  \details This function is surcharged by the PET scanner daughter class 
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::SetPETMaxAxialDiffmm(FLTNB a_maxAxialDiffmm)
{
  Cerr("***** vScanner::SetPETMaxAxialDiffmm() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                       This function only works with PET scanner objects !!" << endl);
  return 1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetPETSpecificParameters
  \param ap_maxRingDiff
  \brief Get geometric PET specific parameters to initialize the datafile
  \details This function is surcharged by the PET scanner daughter classes
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::PROJ_GetPETSpecificParameters(FLTNB* ap_maxRingDiff)
{
  Cerr("***** vScanner::PROJ_GetPETSpecificParameters() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                                   This function only works with PET scanner objects !!" << endl);
  return 1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vScanner::GetSPECTSpecificParameters(uint16_t* ap_nbOfProjections, 
                                         uint16_t* ap_nbHeads,
                                         FLTNB* ap_acquisitionZoom,
                                         uint16_t* ap_nbOfBins,
                                           FLTNB*  ap_pixSizeXY, 
                                           FLTNB*& ap_angles, 
                                           FLTNB*& ap_CORtoDetectorDistance,
                                              int* ap_headRotDirection)
{
  Cerr("***** vScanner::GetSPECTSpecificParameters() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                                This function only works with SPECT scanner objects !!" << endl);
  return 1;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
// ----- CT Specific Functions --- //
/*
  \fn GetCTSpecificParameters
  \param ap_nbOfProjections
  \param ap_angles
  \param ap_detectorRotDirection
  \brief Recover geometric CT specific parameters from the scanner to initialize the datafile
  \details This function is surcharged by the CT scanner daughter classes
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::GetCTSpecificParameters(uint16_t* ap_nbOfProjections, 
                                      FLTNB*& ap_angles, 
                                      int* ap_detectorRotDirection)
{
  Cerr("***** vScanner::GetCTSpecificParameters() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                             This function only works with CT scanner objects !!" << endl);
  return 1;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      SetRotDirection
  \param   a_rotDirection
  \brief   Set rotation direction of the system
  \details Set rotation direction of the scanner elements (head for SPECT, rsector/modules for PET)
           for the generation of the geometry
  \return  0 if success, positive value otherwise (unknown key)
*/
int vScanner::SetRotDirection( string a_rotDirection )
{
  if(     a_rotDirection == "CCW" ||
          a_rotDirection == "Ccw" ||
          a_rotDirection == "ccw" )
    m_rotDirection = GEO_ROT_CCW ;
    
  else if(a_rotDirection == ""   || // Default
          a_rotDirection == "CW" ||
          a_rotDirection == "Cw" ||
          a_rotDirection == "cw" )
    m_rotDirection = GEO_ROT_CW ;
    
  else
  {
    Cerr("***** vScanner::SetRotDirection -> Error while initializing rotation direction !" << endl);
    Cerr("     "<< a_rotDirection <<"' is unknown. Direction must be  'CW' (clockwise) or 'CCW' (counter-clockwise).");
    return 1;
  }
  
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTNbBins
  \param ap_nbOfBins
  \brief Set SPECT number of Bins 
  \details This function is surcharged by the SPECT scanner daughter classes
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::PROJ_SetSPECTNbBins(uint16_t* ap_nbOfBins)
{
  Cerr("***** vScanner::PROJ_SetSPECTNbBins() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                         This function only works with SPECT scanner objects !!" << endl);
  return 1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTNbProjections
  \param a_nbOfProjections
  \brief Set SPECT number of views
  \details This function is surcharged by the SPECT scanner daughter classes
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::PROJ_SetSPECTNbProjections(uint32_t a_nbOfProjections)
{
  Cerr("***** vScanner::PROJ_SetSPECTNbProjections() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                                This function only works with SPECT scanner objects !!" << endl);
  return 1;
}






// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTAngles
  \param ap_projectionAngles
  \brief Set SPECT projection angles
  \details This function is surcharged by the SPECT scanner daughter classes
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::PROJ_SetSPECTAngles(FLTNB* a2p_projectionAngles)
{
  Cerr("***** vScanner::SetSPECTAngles() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                    This function only works with SPECT scanner objects !!" << endl);
  return 1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTCORtoDetectorDistance
  \param a_CORtoDetectorDistance
  \brief Set distance between center of rotation and SPECT detectors
  \details This function is surcharged by the SPECT scanner daughter classes
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
int vScanner::PROJ_SetSPECTCORtoDetectorDistance(FLTNB a_CORtoDetectorDistance)
{
  Cerr("***** vScanner::SetSPECTCORtoDetectorDistance() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                                   This function only works with SPECT scanner objects !!" << endl);
  return 1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetSPECTNbProjections
  \brief return the total number of projections for a SPECT acquisition
  \details This function is surcharged by the SPECT scanner daughter classes
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
uint16_t vScanner::PROJ_GetSPECTNbProjections()
{
  Cerr("***** vScanner::GetSPECTNbProjections() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                           This function only works with SPECT scanner objects !!" << endl);
  return 1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetSPECTNbPixels
  \brief return the total number of pixels for a SPECT reconstruction
  \details This function is surcharged by the SPECT scanner daughter classes
           Returns an error by default.
  \return 1 (error) if not surcharged by a daughter class
*/
uint16_t vScanner::PROJ_GetSPECTNbPixels()
{
  Cerr("***** vScanner::GetSPECTNbPixels() -> This function is not implemented by the Instantiated scanner class !!" << endl);
  Cerr("                                      This function only works with SPECT scanner objects !!" << endl);
  return 1;
}

