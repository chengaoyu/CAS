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
  \brief    Declaration of class iScannerSPECTConv
*/

#ifndef ISCANNERSPECTConv_HH
#define ISCANNERSPECTConv_HH 1

#include "gVariables.hh"
#include "vScanner.hh"
#include "sAddonManager.hh"

/*!
  \class   iScannerSPECTConv
  \brief   This class is used to represent any SPECT camera 
           with parallel/convergent collimator
  \details Inherits from vScanner
*/
class iScannerSPECTConv : public vScanner
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iScannerSPECTConv constructor. 
               Initialize the member variables to their default values.
    */
    iScannerSPECTConv();
    /*!
      \brief   iScannerSPECTConv destructor. 
    */
    ~iScannerSPECTConv();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_SCANNER(iScannerSPECTConv)
    /*!
      \fn      iScannerSPECTConv::DescribeSpecific()
      \brief   Implementation of the pure virtual eponym function that simply prints info about the scanner
    */
    void DescribeSpecific();
    /*!
      \fn      int iScannerSPECTConv::Instantiate()
      \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                                   the system is a generic file (0) or custom Look-up-table (1)
      \brief   Get mandatory informations from the scanner file 
               and allocate memory for the member variables
      \return  0 if success, positive value otherwise
    */
    int Instantiate(bool a_scannerFileIsLUT);
    /*!
      \fn      int iScannerSPECTConv::CheckParameters()
      \brief   Check that all parameters have been correctly initialized.
      \return  0 if success, positive value otherwise
      \todo    Keep the check on crystal depth ?
    */
    int CheckParameters();
    /*!
      \fn      int iScannerSPECTConv::Initialize()
      \brief   Check general initialization and set several parameters to their default value
      \return  0 if success, positive value otherwise
    */
    int Initialize();
    /*!
      \fn      int iScannerSPECTConv::BuildLUT()
      \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                                    the SPECT camera is a generic file (0) or custom Look-up-table (1)
      \brief   Call the functions to generate the LUT or read the user-made LUT depending on the user choice
      \return  0 if success, positive value otherwise
    */
    int BuildLUT( bool a_scannerFileIsLUT );
    /*!
      \fn      intiScannerSPECTConv::GetPositionsAndOrientations()
      \param   a_index1 : 1st index of the event (projection angle)
      \param   a_index2 : 2nd index of the event (crystal index in the gamma camera)
      \param   ap_Position1[3] : x,y,z cartesian position of the focal point
      \param   ap_Position2[3] : x,y,z cartesian position of the crystal
      \param   ap_Orientation1[3] : return -1 by default (no orientation components required for the focal point)
      \param   ap_Orientation2[3] : x,y,z components of the orientation vector related to the crystal
      \param   ap_POI1 : ignored in SPECT (no POI component required for the focal point)
      \param   ap_POI2 : x,y,z components of the Point Of Interation related to the crystal.
                         Currently ignored (POI management not implemented in SPECT)
      \brief   Get the central positions and orientations of the scanner elements from their indices.
      \todo    Point Of Interactions management is NOT implemented and ignored
      \return  0 if success, positive value otherwise
    */
    int GetPositionsAndOrientations( int a_index1, int a_index2,
                                     FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                     FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3],
                                     FLTNB* ap_POI1 = NULL, FLTNB* ap_POI2 = NULL );
    /*!
      \fn      int iScannerSPECTConv::GetRdmPositionsAndOrientations()
      \param   a_index1 : 1st index of the event (projection angle)
      \param   a_index2 : 2nd index of the event (crystal index in the gamma camera)
      \param   ap_Position1[3] : x,y,z cartesian position of the focal point
      \param   ap_Position2[3] : x,y,z cartesian position of the crystal
      \param   ap_Orientation1[3] : return -1 by default (no orientation components required for the focal point)
      \param   ap_Orientation2[3] : x,y,z components of the orientation vector related to the crystal
      \brief   Get the focal point and random positions on the crystal surface and its orientations from the event indices.
      \details - Computed the LUT index described by the projection angle and crystal index passed in parameters. \n 
               - Compute random positions on the surface of the crystal \n
               - Write the corresponding random cartesian coordinates in the positions parameters.
      \todo    This implementation has to be checked and adapted,
      #         as the current implementation of SPECT projection assumes crystal position located on the center of the crystal surface
      \return  0 if success, positive value otherwise
    */
    int GetRdmPositionsAndOrientations( int a_index1, int a_index2,
                                        FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                        FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3] );
    /*!
      \fn      int iScannerSPECTConv::GetPositionWithRandomDepth()
      \param   a_index1 :
      \param   a_index2 : index of the crystal
      \param   ap_Position1[3] :
      \param   ap_Position2[3] : x,y,z cartesian position of the point related to the crystal
      \brief   Get the positions and orientations of scanner elements from their indices, with a random depth.
      \todo    Not yet implemented
      \return  0 if success, positive value otherwise
    */
    int GetPositionWithRandomDepth( int a_index1, int a_index2, FLTNB ap_Position1[3], FLTNB ap_Position2[3] );
    /*!
      \fn      int iScannerSPECTConv::GetTwoCorners()
      \param   a_index1 : index of the projection angle
      \param   a_index2 : index of the  crystal 
      \param   ap_CornerInf1[3]
      \param   ap_CornerSup1[3]
      \param   ap_CornerInf2[3]
      \param   ap_CornerSup2[3]
      \brief   Get the cartesian coordinaters of the two opposite corners of a scanner element.
      \todo    Not implemented yet 
      \return  0 if success, positive value otherwise
    */
    int GetTwoCorners( int a_index1, int a_index2,
                       FLTNB ap_CornerInf1[3], FLTNB ap_CornerSup1[3],
                       FLTNB ap_CornerInf2[3], FLTNB ap_CornerSup2[3] );
    /*!
      \fn      int iScannerSPECTConv::GetEdgesCenterPositions()
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
      \brief   Implementation of the pure virtual function from vScanner.  \n
               Not yet implemented for the SPECT convergent systems. Simply return 1. \n
               It is typically used for the Distance Driven projector
      \return  0 if success, positive value otherwise
    */
    int GetEdgesCenterPositions( int a_index1, int a_index2,
      FLTNB ap_pos_line_point1[3], FLTNB ap_pos_line_point2[3],
      FLTNB ap_pos_point1_x[4], FLTNB ap_pos_point1_y[4], FLTNB ap_pos_point1_z[4],
      FLTNB ap_pos_point2_x[4], FLTNB ap_pos_point2_y[4], FLTNB ap_pos_point2_z[4] );
    /*!
      \fn      int iScannerSPECTConv::GetGeometricInfoFromDataFile()
      \param   a_pathToDF : string containing the path to datafile header
      \brief   Recover geometric informations specific to the scanner class from the datafile header
      \details -Recover nb of bins and projections \n
               -Recover the projection angles from the header
               (directly read from the datafile, or extrapolated from a first and last angle) \n
               -Recover the distance between the gamma camera detector surfaces and the center of rotation \n
               (directly read from the datafile, or extracted from the gamma camera configuratino file)
      \return  0 if success, positive value otherwise
    */
    int GetGeometricInfoFromDataFile( string a_pathToDF );
    /*!
      \fn      inline int iScannerSPECTConv::GetSystemNbElts()
      \brief   Get the number of elements in the system. For a SPECT system, 
               returns the number of crystal in one gamma camera head
      \return  a number of crystals
    */
    inline int GetSystemNbElts()
           {return m_nbCrystals;}
    /*!
      \fn      void iScannerSPECTConv::ShowHelp()
      \brief   Display help
      \todo    Provide informations about SPECT system initialization ?
    */
    void ShowHelp();
    /*!
      \fn      inline FLTNB iScannerSPECTConv::GetDetectionElementSizeTrans()
      \return  return the transaxial size of the detection element
      \todo    take the multi layer into account
    */
    inline FLTNB GetDetectionElementSizeTrans()
           {return m_vPixelsSizeTrans;}
    /*!
      \fn      inline FLTNB iScannerSPECTConv::GetDetectionElementSizeAxial()
      \return  return the axial size of the detection element
      \todo    take the multi layer into account
    */
    inline FLTNB GetDetectionElementSizeAxial()
           {return m_vPixelsSizeTrans;}

  // -------------------------------------------------------------------
  // Functions dedicated to Analytic Projection
  public:
    /*!
      \fn      inline int iScannerSPECTConv::IsAvailableLOR()
      \param   a_elt1 : index of the 1st crystal
      \param   a_elt2 : index of the 2nd crystal
      \brief   Check if the LOR formed by the crystalf whose indices are passed 
               in parameters is available according to the scanner restrictions \n
               It is dedicated to PET and return 1 by default for SPECT
      \return  1 if the LOR is available, 0 otherwise
    */
    inline int IsAvailableLOR(int a_elt1, int a_elt2)
           {return 1;};
    /*!
      \fn      int iScannerSPECTConv::GetSPECTSpecificParameters()
      \param   ap_nbOfProjections : total number of views
      \param   ap_nbHeads : number of heads in the system
      \param   ap_acquisitionZoom : zoom used during acquisition for monolithic detectors
      \param   ap_nbOfBins : 2 elements array containing transaxial/axial number of pixels
      \param   ap_pixSizeXY : 2 elements array containing transaxial/axial pixel sizes
      \param   ap_angles : Array containing angles of each projection view
      \param   ap_CORtoDetectorDistance : Radius (distance between FOV center and detector)
      \param   ap_headRotDirection : head rotation direction
      \brief   Set pointers passed in argument with the related SPECT specific variables \n
               This function is used to recover these values in the datafile object
      \return  0 if success, positive value otherwise
    */
    int GetSPECTSpecificParameters( uint16_t* ap_nbOfProjections, 
                                    uint16_t* ap_nbHeads,
                                    FLTNB* ap_acquisitionZoom, 
                                    uint16_t* ap_nbOfBins,
                                    FLTNB*  ap_pixSizeXY, 
                                    FLTNB*& ap_angles, 
                                    FLTNB*& ap_CORtoDetectorDistance,
                                    int* ap_headRotDirection);
    /*!
      \fn      int iScannerSPECTConv::PROJ_SetSPECTAngles()
      \param   ap_projectionAngles : an array containing the projection angles
      \brief   Set the projection angles with the array provided in parameter
      \return  0 if success, positive value otherwise
    */
    int PROJ_SetSPECTAngles(FLTNB* ap_projectionAngles);
    /*!
      \fn      int iScannerSPECTConv::PROJ_SetSPECTCORtoDetectorDistance()
      \param   a_distance
      \brief   Set distance between the center of rotation and SPECT detectors if arg value>0, \n
               Set with the geometric information in the scanner configuration file otherwise
      \return  0 if success, positive value otherwise
    */
    int PROJ_SetSPECTCORtoDetectorDistance(FLTNB a_distance);
    /*!
      \fn      inline uint16_t iScannerSPECTConv::PROJ_GetSPECTNbProjections()
      \return  the number of projection angles
    */
    inline uint16_t PROJ_GetSPECTNbProjections()
           {return m_nbOfProjections;}
    /*!
      \fn      inline uint16_t iScannerSPECTConv::PROJ_GetSPECTNbPixels()
      \return  the total number of pixels in a projection
    */
    inline uint16_t PROJ_GetSPECTNbPixels()
           {return m_vNbPixelsTrans*m_vNbPixelsAxial;}
    /*!
      \fn      int iScannerSPECTConv::PROJ_GetSPECTNbBins()
      \brief   Get the number of SPECT heads in the pointer provided in parameter
      \return  0 by default (no error)
    */
    int PROJ_GetSPECTNbBins(uint16_t* ap_nbOfBins);
    /*!
      \fn      int iScannerSPECTConv::PROJ_SetSPECTNbBins()
      \param   ap_nbOfBins
      \brief   Set number of bins
      \return  0 by default (no error)
    */
    int PROJ_SetSPECTNbBins(uint16_t* ap_nbOfBins);
    /*!
      \fn      int iScannerSPECTConv::PROJ_SetSPECTNbProjections()
      \param   a_nbOfProjections
      \brief   Set number of projections
      \return  0 by default (no error)
    */
    int PROJ_SetSPECTNbProjections(uint32_t a_nbOfProjections);

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      int iScannerSPECTConv::LoadLUT()
      \brief   Load a precomputed scanner LUT
      \details Read mandatory data from the header of the LUT. Then load the LUT elements for each crystal
      \todo    Not yet implemented
      \return  0 if success, positive value otherwise
    */
    int LoadLUT();
    /*!
      \fn      int iScannerSPECTConv::ComputeLUT()
      \brief   Computes the LUT of the scanner from a generic (.geom) file.
      \details Read mandatory data from the geom file. Then compute the LUT elements for each crystal from the geometry described in the file \n
               Compute the look-up-tables of the system containing the locations of the scanner elements center in space and their orientations
      \todo    center of rotation & head first angles : get this from acquisition header file
      \return  0 if success, positive value otherwise
    */
    int ComputeLUT();
    /*!
      \fn      int iScannerSPECTConv::ComputeFocalPositions()
      \param   a_posX : cartesian position of the crystal on the x-axis 
      \param   a_posY : cartesian position of the crystal on the y-axis 
      \param   a_posZ : cartesian position of the crystal on the z-axis 
      \param   a_headID : head index of the crystal 
                         (required to select the related focal model parameters
                          and distance to center of rotation)
      \param   a_cryID : crystal index in the LUT
      \brief   Compute focal positions for a specific crystal ID
      \details Compute the focal positions using the implemented "constant", 
               "polynomial", "slanthole" focal models related to the gamma cameras \n
               The "custom" focal model is dedicated to user-made focal model
               and should be implemented by the user
      \return  0 if success, positive value otherwise
    */
    int ComputeFocalPositions( FLTNB a_posX, FLTNB a_posY, FLTNB a_posZ, int a_headID, int a_cryID );
  
  
  // -------------------------------------------------------------------
  // Data members
  private:
    int m_nbCrystals;                  /*!< Total number of crystal in the scanner */
    int m_nbHeads;                     /*!< Total number of SPECT heads*/
    
    uint16_t mp_nbOfBins[2];           /*!< 2 dimensionnal array containing the number of transaxial bins (X,Y). Default values : 2,2*/
    uint16_t m_nbOfProjections;        /*!< Total number of projection angles*/
    
    FLTNB* mp_projectionAngles;        /*!< Array containing all the projection angles ('m_nbOfProjections' elements)*/
    FLTNB* mp_CORtoDetectorDistance;   /*!< Distance between the center of rotation and the detector surface of each head. One value for each projection angle.*/
    FLTNB* mp_radius;                  /*!< Default radius (distance between the center of rotation and detector surface) for each head for the system*/
    
    uint32_t m_nbPixelsTrans;          /*!< Total number of transaxial pixels as defined in the system file*/
    FLTNB m_pixelsSizeTrans;           /*!< Size of transaxial pixels as defined in the system file*/
    FLTNB m_gapSizeTrans;              /*!< Gap size between each transaxial pixe as defined in the system filel*/

    uint32_t m_nbPixelsAxial;          /*!< Total number of axial pixels as defined in the system file*/
    FLTNB m_pixelsSizeAxial;           /*!< Size of axial pixels as defined in the system file*/
    FLTNB m_gapSizeAxial;              /*!< Gap size between each axial pixel as defined in the system file*/

    FLTNB m_acquisitionZoom;                      /*!< The zoom used during the acquisition to limit the area of detection with monolithic crystals */

    uint32_t m_vNbPixelsTrans;         /*!< Number of trans virtual pixels (pixels actually used in reconstruction)*/
    uint32_t m_vNbPixelsAxial;         /*!< Number of axial virtual pixels (pixels actually used in reconstruction)*/
    
    FLTNB m_vPixelsSizeTrans;          /*!< Trans size of virtual pixels (pixels actually used in reconstruction)*/
    FLTNB m_vPixelsSizeAxial;          /*!< Axial size of virtual pixels (pixels actually used in reconstruction)*/

    FLTNB m_crystalDepth;              /*!< Depth of crystals*/
 
    string* mp_focalModelTrans;        /*!< Transaxial focal model (should be 'constant', 'polynomial', 'slanthole', or 'custom'). Specific to each head*/
    uint8_t* mp_nbCoefModelTrans;      /*!< Number of coefficients of the transaxial focal model. Specific to each head*/
    FLTNB** m2p_transFocalParameters;  /*!< Parameters of the transaxial focal model. Specific to each head*/

    string* mp_focalModelAxial;        /*!< Axial focal model (should be 'constant', 'polynomial', 'slanthole', or 'custom'). Specific to each head*/
    uint8_t* mp_nbCoefModelAxial;      /*!< Number of coefficients of the axial focal model. Specific to each head*/
    FLTNB** m2p_axialFocalParameters;  /*!< Parameters of the axial focal model. Specific to each head*/

    FLTNB* mp_crystalCentralPositionX; /*!< Cartesian coordinate on X-axis of the center of each crystal, at each projection*/
    FLTNB* mp_crystalCentralPositionY; /*!< Cartesian coordinate on Y-axis of the center of each crystal, at each projection*/
    FLTNB* mp_crystalCentralPositionZ; /*!< Cartesian coordinate on Z-axis of the center of each crystal, at each projection*/

    FLTNB* mp_crystalOrientationX;     /*!< X-axis orientation of each crystal, at each projection*/
    FLTNB* mp_crystalOrientationY;     /*!< Y-axis orientation of each crystal, at each projection*/
    FLTNB* mp_crystalOrientationZ;     /*!< Z-axis orientation of each crystal, at each projection*/
 
    FLTNB* mp_crystalFocalPositionX;   /*!< X-axis position of the focal point for each crystal, at each projection*/
    FLTNB* mp_crystalFocalPositionY;   /*!< Y-axis position of the focal point for each crystal, at each projection*/
    FLTNB* mp_crystalFocalPositionZ;   /*!< Z-axis position of the focal point for each crystal, at each projection*/
};


// Class for automatic insertion (set here the visible scanner type name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_SCANNER(SPECT_CONVERGENT,iScannerSPECTConv)

#endif
