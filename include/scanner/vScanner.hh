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
  \brief    Declaration of class vScanner
*/

#ifndef VSCANNER_HH
#define VSCANNER_HH 1

#include "gVariables.hh"
#include "sOutputManager.hh"
#include "oMatrix.hh"
#include "sRandomNumberGenerator.hh"

/**
 * @defgroup GEO_ROT system geometry rotation direction
 *
 *    \brief Rotation orientation for angles computations \n
 *           Defined in vScanner.hh
 * @{
 */
/** Constant corresponding to a clockwise rotation direction (=0) */
#define GEO_ROT_CW 0
/** Constant corresponding to a counter-clockwise rotation direction (=1) */
#define GEO_ROT_CCW 1
/** @} */

class sScannerManager;

/*!
  \class   vScanner
  \brief   Generic class for scanner objects.
  \details This class is designed to be a mother virtual class that should not be used on its own; only its children are used. \n
           It has a pure virtual function that is used to get the cartesian coordinates of any line described by two indices. \n
           Although new scanner classes are not supposed to be added, it uses the same auto-declaration class system than other 
           reconstruction classes for the initialization of its children.
*/
class vScanner
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   vScanner constructor. 
               Initialize the member variables to their default values.
    */
    vScanner();
    /*!
      \brief   vScanner destructor.
    */
    virtual ~vScanner();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      int vScanner::Describe()
      \brief   A function used to describe the generic parts of the datafile
    */
    void Describe();
    /*!
      \fn      vScanner::DescribeSpecific() = 0
      \brief   A pure virtual function used to describe the specific parts of the scanner
    */
    virtual void DescribeSpecific() = 0;
    /*!
      \fn      virtual int vScanner::Instantiate() = 0
      \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                                   the system is a generic file (0) or custom Look-up-table (1)
      \brief   This function is implemented in child classes.
               Read the mandatory fields from the scanner file and allocate memory for the member variables
      \return  0 if success, positive value otherwise
    */
    virtual int Instantiate(bool a_scannerFileIsLUT) = 0;
    /*!
      \fn      virtual int vScanner::BuildLUT() = 0
      \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                                   the system is a generic file (0) or custom Look-up-table (1)
      \brief   This function is implemented in child classes. \n
               Instantiate the scanner look-up-table (LUT)
      \return  0 if success, positive value otherwise
    */
    virtual int BuildLUT(bool a_scannerFileIsLUT) = 0;
    /*!
      \fn      virtual int vScanner::CheckParameters() = 0
      \brief   This function is implemented in child classes. \n
               Check that all parameters have been correctly initialized.
      \return  0 if success, positive value otherwise
    */
    virtual int CheckParameters() = 0;
    /*!
      \fn      virtual int vScanner::Initialize() = 0
      \brief   This function is implemented in child classes.  \n
               Check initialization and set several parameters to their default value
      \return  0 if success, positive value otherwise
    */
    virtual int Initialize() = 0;
    /*!
      \fn      virtual int vScanner::GetPositionsAndOrientations() = 0
      \param   a_index1 : 1st index of the event
      \param   a_index2 : 2nd index of the event
      \param   ap_Position1[3] : x,y,z cartesian position of the point related to the 1st index (see child function for more details)
      \param   ap_Position2[3] : x,y,z cartesian position of the point related to the 2st index (see child function for more details)
      \param   ap_Orientation1[3] : x,y,z components of the orientation vector related to the 1st index (see child function for more details)
      \param   ap_Orientation2[3] : x,y,z components of the orientation vector related to the 2nd index (see child function for more details)
      \param   ap_POI1 : x,y,z components of the Point Of Interation related to the 1st index (see child function for more details)
      \param   ap_POI2 : x,y,z components of the Point Of Interation related to the 2nd index (see child function for more details)
      \brief   This is a pure virtual method that must be implemented by children.  \n
               Get the central positions and orientations of the scanner elements from their indices.
      \return  0 if success, positive value otherwise
    */
    virtual int GetPositionsAndOrientations( int a_index1, int a_index2,
                                             FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                             FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3],
                                             FLTNB* ap_POI1 = NULL, FLTNB* ap_POI2 = NULL ) = 0;
    /*!
      \fn      virtual int vScanner::GetRdmPositionsAndOrientations() = 0
      \param   a_index1 : 1st index of the event
      \param   a_index2 : 2nd index of the event
      \param   ap_Position1[3] : x,y,z cartesian position of the point related to the 1st index (see child function for more details)
      \param   ap_Position2[3] : x,y,z cartesian position of the point related to the 2st index (see child function for more details)
      \param   ap_Orientation1[3] : x,y,z components of the orientation vector related to the 1st index (see child function for more details)
      \param   ap_Orientation2[3] : x,y,z components of the orientation vector related to the 2nd index (see child function for more details)
      \brief   This is a pure virtual method that must be implemented by children.  \n
               Get random positions of the scanner elements and their orientations from their indices.
      \return  0 if success, positive value otherwise
    */
    virtual int GetRdmPositionsAndOrientations( int a_index1, int a_index2,
                                                FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                                FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3] ) = 0;
    /*!
      \fn      virtual int vScanner::GetPositionWithRandomDepth() = 0
      \param   a_index1 : 1st index of the event
      \param   a_index2 : 2nd index of the event
      \param   ap_Position1[3] : x,y,z cartesian position of the point related to the 1st event (see child function for more details)
      \param   ap_Position2[3] : x,y,z cartesian position of the point related to the 2st event (see child function for more details)
      \brief   This is a pure virtual method that must be implemented by children.  \n
               Get the positions of scanner elements from their indices, with a random depth.
      \return  0 if success, positive value otherwise
    */
    virtual int GetPositionWithRandomDepth( int a_index1, int a_index2, FLTNB ap_Position1[3], FLTNB ap_Position2[3] ) = 0;
    /*!
      \fn      virtual int vScanner::GetTwoCorners() = 0
      \param   a_index1 : 1st index of the event
      \param   a_index2 : 2nd index of the event
      \param   ap_CornerInf1[3]
      \param   ap_CornerSup1[3]
      \param   ap_CornerInf2[3]
      \param   ap_CornerSup2[3]
      \brief   This is a pure virtual method that must be implemented by children.  \n
               Get the cartesian coordinaters of the two opposite corners of a scanner element.
      \todo    Not implemented yet 
      \return  0 if success, positive value otherwise
    */
    virtual int GetTwoCorners( int a_index1, int a_index2,
                               FLTNB ap_CornerInf1[3], FLTNB ap_CornerSup1[3],
                               FLTNB ap_CornerInf2[3], FLTNB ap_CornerSup2[3]) = 0;
    /*!
      \fn      virtual int vScanner::GetEdgesCenterPositions() = 0
      \param   int a_index1 : 1st index of the event
      \param   int a_index2 : 2nd index of the event
      \param   FLTNB ap_pos_line_point1[3] : current position of point 1 of the ray (supposed to be centered in the section of the detection element and at the desired depth of interaction)
      \param   FLTNB ap_pos_line_point2[3] : current position of point 1 of the ray (supposed to be centered in the section of the detection element and at the desired depth of interaction)
      \param   FLTNB ap_pos_point1_x[4] : the resulting X position of the 4 edges of point 1
      \param   FLTNB ap_pos_point1_y[4] : the resulting Y position of the 4 edges of point 1
      \param   FLTNB ap_pos_point1_z[4] : the resulting Z position of the 4 edges of point 1
      \param   FLTNB ap_pos_point2_x[4] : the resulting X position of the 4 edges of point 1
      \param   FLTNB ap_pos_point2_y[4] : the resulting Y position of the 4 edges of point 1
      \param   FLTNB ap_pos_point2_z[4] : the resulting Z position of the 4 edges of point 1
      \brief   This is a pure virtual method that must be implemented by children.  \n
               Get the cartesian coordinaters of the center of the 4 edges of the detection element. \n
               It is typically used for the Distance Driven projector
      \return  0 if success, positive value otherwise
    */
    virtual int GetEdgesCenterPositions( int a_index1, int a_index2,
      FLTNB ap_pos_line_point1[3], FLTNB ap_pos_line_point2[3],
      FLTNB ap_pos_point1_x[4], FLTNB ap_pos_point1_y[4], FLTNB ap_pos_point1_z[4],
      FLTNB ap_pos_point2_x[4], FLTNB ap_pos_point2_y[4], FLTNB ap_pos_point2_z[4] ) = 0;
    /*!
      \fn      virtual int vScanner::ComputeLUT()
      \brief   Virtual function which should be implemented by the child classes. \n
               It computes the LUT of the scanner from a generic (.geom) file. \n
               The vScanner implementation throws error by default as it should be implemented by the child class
      \todo    Make the mother function virtual pure ? (should iScannerSinogram  also implement this function ?) \n
      #         iScannerCT implementation will probably consists in computing one projection, then computing
      #         the others on-the-fly during reconstruction (Often too much data to keep in memory for CT) \n
      #         Have to check if we offer the precomputation of the entire LUT in some situation, as for PET/SPECT
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int ComputeLUT();
    /*!
      \fn      virtual int vScanner::LoadLUT()
      \brief   Virtual function which should be implemented by the child classes. \n
               Load a precomputed scanner LUT.  \n
               The vScanner implementation throws error by default as it should be implemented by the child class
      \todo    Make the mother function virtual pure ? (should iScannerSinogram  also implement this function ?) \n
      #         iScannerCT implementation will probably consists in computing one projection, then computing
      #         the others on-the-fly during reconstruction (Often too much data to keep in memory for CT) \n
      #         Have to check if we offer the precomputation of the entire LUT in some situation, as for PET/SPECT
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int LoadLUT();
    /*!
      \fn      virtual int vScanner::GetSystemNbElts() = 0
      \brief   This is a pure virtual method that must be implemented by children
      \return  the number of elements in the system
    */
    virtual int GetSystemNbElts() = 0;
    /*!
      \fn      virtual int vScanner::IsAvailableLOR()
      \param   a_elt1 : index of the 1st scanner element
      \param   a_elt2 : index of the 2nd scanner element
      \brief   This function is implemented in child classes. \n
               Check if the LOR is available according to the scanner restrictions
      \details This function is related to analytic projection and list-mode sensitivity image generation
      \return  1 if the LOR is available, 0 otherwise
    */
    virtual int IsAvailableLOR(int a_elt1, int a_elt2);
    /*!
      \fn      virtual void vScanner::ShowHelp() = 0
      \brief   This function is implemented in child classes \n
               Display help specific to the scanner class
    */
    virtual void ShowHelp() = 0;
    /*!
      \fn      virtual int vScanner::GetGeometricInfoFromDataFile() = 0
      \param   a_path : string containing the path to datafile header
      \brief   This function is implemented in child classes \n
               Recover geometric informations specific to the scanner class from the datafile header
      \return  0 if success. Positive value otherwise
    */
    virtual int GetGeometricInfoFromDataFile(string a_path) = 0;
    /*!
      \fn      inline int vScanner::GetScannerType()
      \return  the type (modality) of the system as defined by macro
               SCANNER_TYPE in sScannerManager.hh
    */
    inline int GetScannerType()
           {return m_scannerType;}
    /*!
      \fn      inline int vScanner::GetScannerTypeString()
      \return  the string corresponding to the scanner type (modality) 
               as defined by macro SCANNER_TYPE in sScannerManager.hh
    */
    string  GetScannerTypeString();
    /*!
      \fn      inline void vScanner::SetVerbose()
      \param   a_verboseLevel
      \brief   Set verbosity
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      inline void vScanner::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ID
      \brief   Set the pointer to the image dimensions and quantification object
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ID)
           {mp_ID = ap_ID;}
    /*!
      \fn      inline FLTNB vScanner::GetDefaultBedDisplacementInMm()
      \return  the default bed displacement between two successive bed positions
    */
    inline FLTNB GetDefaultBedDisplacementInMm()
           {return m_defaultBedDisplacementInMm;}
    /*!
      \fn      inline virtual FLTNB vScanner::GetDetectionElementSizeTrans() = 0
      \return  return the transaxial size of the detection element, function implemented by children
    */
    inline virtual FLTNB GetDetectionElementSizeTrans() = 0;
    /*!
      \fn      inline virtual FLTNB vScanner::GetDetectionElementSizeAxial() = 0
      \return  return the axial size of the detection element, function implemented by children
    */
    inline virtual FLTNB GetDetectionElementSizeAxial() = 0;

  // -------------------------------------------------------------------
  // PET dedicated functions (only overloaded by PET scanner class)
  public:
    /*!
      \fn      virtual int vScanner::SetPETMaxAxialDiffmm()
      \param   a_maxAxialDiffmm
      \brief   Set the maximal axial difference in mm between 2 crystals forming a lor
      \details This function is surcharged by the PET scanner daughter class  \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int SetPETMaxAxialDiffmm(FLTNB a_maxAxialDiffmm);
  
    /*!
      \fn      SetRotDirection
      \param   a_rotDirection
      \brief   Set rotation direction of the system
      \details Set rotation direction of the scanner elements (head for SPECT, rsector/modules for PET)
               for the generation of the geometry
      \return  0 if success, positive value otherwise (unknown key)
    */
    virtual int SetRotDirection( string a_rotDirection );

    /*!
      \fn      virtual int vScanner::PROJ_GetPETSpecificParameters()
      \param   ap_maxRingDiff
      \brief   Get geometric PET specific parameters to initialize the datafile
      \details This function is surcharged by the PET scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int PROJ_GetPETSpecificParameters(FLTNB* ap_maxRingDiff);


  // -------------------------------------------------------------------
  // SPECT dedicated functions (only surcharged by SPECT scanner classes)
  public:
    /*!
      \fn      virtual int vScanner::GetSPECTSpecificParameters()
      \param   ap_nbOfProjections
      \param   ap_nbHeads
      \param   ap_acquisitionZoom
      \param   ap_nbOfBins
      \param   ap_pixSizeXY
      \param   ap_angles
      \param   ap_CORtoDetectorDistance
      \param   ap_headRotDirection
      \brief   Recover geometric SPECT specific parameters from the scanner to initialize the datafile
      \details This function is surcharged by the SPECT scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int GetSPECTSpecificParameters( uint16_t* ap_nbOfProjections, 
                                            uint16_t* ap_nbHeads,
                                            FLTNB* ap_acquisitionZoom,
                                            uint16_t* ap_nbOfBins,
                                            FLTNB*  ap_pixSizeXY, 
                                            FLTNB*& ap_angles, 
                                            FLTNB*& ap_CORtoDetectorDistance,
                                            int* ap_headRotDirection );
    /*!
      \fn      virtual int vScanner::GetCTSpecificParameters()
      \param   ap_nbOfProjections
      \param   ap_angles
      \param   ap_detectorRotDirection
      \brief   Recover geometric CT specific parameters from the scanner to initialize the datafile
      \details This function is surcharged by the CT scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not overloaded by a daughter class
    */
    virtual int GetCTSpecificParameters( uint16_t* ap_nbOfProjections, 
                                         FLTNB*& ap_angles, 
                                         int* ap_detectorRotDirection );
     /*!
      \fn      virtual int vScanner::PROJ_SetSPECTNbBins()
      \param   ap_nbOfBins
      \brief   Set SPECT number of Bins 
      \details This function is surcharged by the SPECT scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int PROJ_SetSPECTNbBins( uint16_t* ap_nbOfBins );
    /*!
      \fn      virtual int vScanner::PROJ_SetSPECTNbProjections()
      \param   a_nbOfProjections
      \brief   Set SPECT number of views
      \details This function is surcharged by the SPECT scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int PROJ_SetSPECTNbProjections( uint32_t a_nbOfProjections );
    /*!
      \fn      virtual int vScanner::PROJ_SetSPECTAngles()
      \param   ap_projectionAngles
      \brief   Set SPECT projection angles
      \details This function is surcharged by the SPECT scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int PROJ_SetSPECTAngles( FLTNB* ap_projectionAngles );
    /*!
      \fn      virtual int vScanner::PROJ_SetSPECTCORtoDetectorDistance()
      \param   a_CORtoDetectorDistance
      \brief   Set distance between center of rotation and SPECT detectors
      \details This function is surcharged by the SPECT scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual int PROJ_SetSPECTCORtoDetectorDistance( FLTNB a_CORtoDetectorDistance );
    /*!
      \fn      virtual uint16_t vScanner::PROJ_GetSPECTNbProjections()
      \brief   return the total number of projections for a SPECT acquisition
      \details This function is surcharged by the SPECT scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual uint16_t PROJ_GetSPECTNbProjections(); // required for analytic projection (computations of indices for projection main loops in sScannerManager)
    /*!
      \fn      virtual uint16_t vScanner::PROJ_GetSPECTNbPixels()
      \brief   return the total number of pixels for a SPECT reconstruction
      \details This function is surcharged by the SPECT scanner daughter classes \n
               Returns an error by default.
      \return  1 (error) if not surcharged by a daughter class
    */
    virtual uint16_t PROJ_GetSPECTNbPixels(); // required for analytic projection (computations of index for projection main loops in sScannerManager)

    
  // -------------------------------------------------------------------
  // Private member functions
  private:


  // -------------------------------------------------------------------
  // Data members
  protected:
    int m_scannerType;                         /*!< System type. Default =-1 (unknown) */
    int m_verbose;                             /*!< Verbosity. Default =-1 */
    oImageDimensionsAndQuantification* mp_ID;  /*!< Image dimensions and quantification object */
    bool m_allParametersChecked;               /*!< Boolean indicating if all variables of the class have been checked. Default =false */
    oMatrix *mp_rotationMatrix;                /*!< Matrix strucure for rotation computations*/
    oMatrix *mp_positionMatrix_ref;            /*!< Matrix strucure for rotation computations*/
    oMatrix *mp_positionMatrix_out;            /*!< Matrix strucure for rotation computations*/
    int m_rotDirection;                        /*!< Rotation Direction for the generation of the geometry (clockwise by default) */
    FLTNB m_defaultBedDisplacementInMm;        /*!< The default displacement of the scanner between two successive bed positions. Default = 0. */
};


// ----------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ----------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_SCANNER(CLASS) \
  static vScanner *make_scanner() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_SCANNER(NAME,CLASS)                                                             \
  class NAME##ScannerCreator                                                                  \
  {                                                                                           \
    public:                                                                                   \
      NAME##ScannerCreator()                                                                  \
        { sAddonManager::GetInstance()->mp_listOfScannerTypes[#NAME] = CLASS::make_scanner; } \
  };                                                                                          \
  static NAME##ScannerCreator ScannerCreator##NAME;

#endif
