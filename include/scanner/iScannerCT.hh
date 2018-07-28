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
  \brief    Declaration of class iScannerCT
*/

#ifndef ISCANNERCT_HH
#define ISCANNERCT_HH 1

#include "gVariables.hh"
#include "vScanner.hh"
#include "sAddonManager.hh"

/*!
  \class   iScannerCT
  \brief   This class is used to represent any CT camera 
           with either a CBCT ascii description or a LUT file for any
           style of detector
  \details Inherits from vScanner
*/
class iScannerCT : public vScanner
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iScannerCT constructor. 
               Initialize the member variables to their default values.
    */
    iScannerCT();
    /*!
      \brief   iScannerCT destructor.
    */
    ~iScannerCT();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_SCANNER(iScannerCT)
    /*!
      \fn      iScannerCT::DescribeSpecific()
      \brief   Implementation of the pure virtual eponym function that simply prints info about the scanner
    */
    void DescribeSpecific();
    /*!
      \fn      int iScannerCT::Instantiate()
      \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                                   the system is a generic file (0) or custom Look-up-table (1)
      \brief   Get mandatory informations from the scanner file 
               and allocate memory for the member variables
      \return  0 if success, positive value otherwise
    */
    int Instantiate(bool a_scannerFileIsLUT);
    /*!
      \fn      int iScannerCT::CheckParameters()
      \brief   Check that all parameters have been correctly initialized.
      \return  0 if success, positive value otherwise
      \todo    Keep the check on crystal depth ?
    */
    int CheckParameters();
    /*!
      \fn      int iScannerCT::Initialize()
      \brief   Check general initialization and set several parameters to their default value
      \return  0 if success, positive value otherwise
    */
    int Initialize();
    /*!
      \fn      int iScannerCT::BuildLUT()
      \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                                    the CT camera is a generic file (0) or custom Look-up-table (1)
      \brief   Call the functions to generate the LUT or read the user-made LUT depending on the user choice
      \return  0 if success, positive value otherwise
    */
    int BuildLUT( bool a_scannerFileIsLUT );
    /*!
      \fn      initScannerCT::GetPositionsAndOrientations()
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
      \fn      int iScannerCT::GetRdmPositionsAndOrientations()
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
      \fn      int iScannerCT::GetPositionWithRandomDepth()
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
      \fn      int iScannerCT::GetTwoCorners()
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
      \fn      int iScannerCT::GetEdgesCenterPositions()
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
               Get the cartesian coordinaters of the center of the 4 edges of the detection element and the source. \n
               It is typically used for the Distance Driven projector
      \return  0 if success, positive value otherwise
    */
    int GetEdgesCenterPositions( int a_index1, int a_index2,
      FLTNB ap_pos_line_point1[3], FLTNB ap_pos_line_point2[3],
      FLTNB ap_pos_point1_x[4], FLTNB ap_pos_point1_y[4], FLTNB ap_pos_point1_z[4],
      FLTNB ap_pos_point2_x[4], FLTNB ap_pos_point2_y[4], FLTNB ap_pos_point2_z[4] );
    /*!
      \fn      int iScannerCT::GetGeometricInfoFromDataFile()
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
      \fn      inline int iScannerCT::GetSystemNbElts()
      \brief   Get the number of elements in the system. For a CT system, 
               returns the number of pixels in the detector
      \return  a number of crystals
    */
    inline int GetSystemNbElts()
           {return m_nbPixels;}
    /*!
      \fn      void iScannerCT::ShowHelp()
      \brief   Display help
      \todo    Provide informations about CT system initialization ?
    */
    void ShowHelp();
    /*!
      \fn      inline FLTNB iScannerCT::GetDetectionElementSizeTrans()
      \return  return the transaxial size of the detection element
    */
    inline FLTNB GetDetectionElementSizeTrans()
           {return m_pixelsSizeTrans;}
    /*!
      \fn      inline FLTNB iScannerCT::GetDetectionElementSizeAxial()
      \return  return the axial size of the detection element
    */
    inline FLTNB GetDetectionElementSizeAxial()
           {return m_pixelsSizeAxial;}
    /*!
      \fn      int iScannerCT::GetCTSpecificParameters()
      \param   ap_nbOfProjections : total number of views
      \param   ap_angles : Array containing angles of each projection view
      \param   ap_detectorRotDirection : detector rotation direction
      \brief   Set pointers passed in argument with the related CT specific variables \n
               This function is used to recover these values in the datafile object
      \return  0 if success, positive value otherwise
    */
    int GetCTSpecificParameters( uint16_t* ap_nbOfProjections, 
                                 FLTNB*& ap_angles, 
                                 int* ap_detectorRotDirection );

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      int iScannerCT::LoadLUT()
      \brief   Load a precomputed scanner LUT
      \details Read mandatory data from the header of the LUT. Then load the LUT elements for each pixel
      \todo    Not yet implemented
      \return  0 if success, positive value otherwise
    */
    int LoadLUT();
    /*!
      \fn      int iScannerCT::ComputeLUT()
      \brief   Computes the LUT of the scanner from a generic (.geom) file.
      \details Read mandatory data from the geom file. Then compute the LUT elements for each crystal from the geometry described in the file \n
               Compute the look-up-tables of the system containing the locations of the scanner elements center in space and their orientations
      \todo    center of rotation & head first angles : get this from acquisition header file
      \return  0 if success, positive value otherwise
    */
    int ComputeLUT();  
  
  // -------------------------------------------------------------------
  // Data members
  private:
    int m_nbPixels;             /*!< Total number of pixels in the detector */
    
    uint16_t m_nbOfProjections;      /*!< Total number of projection angles*/
    FLTNB* mp_projectionAngles;      /*!< Array containing all the projection angles ('m_nbOfProjections' elements)*/
    FLTNB m_CORtoDetectorDistance;   /*!< Distance between the center of rotation and the detector surface. One value for each projection angle.*/
    FLTNB m_CORtoSourceDistance;     /*!< Distance between the center of rotation and the source surface. One value for each projection angle.*/
    
    uint32_t m_nbPixelsTrans;        /*!< Total number of transaxial pixels as defined in the system file*/
    FLTNB m_pixelsSizeTrans;         /*!< Size of transaxial pixels as defined in the system file*/
    FLTNB m_gapSizeTrans;            /*!< Gap size between each transaxial pixe as defined in the system filel*/

    uint32_t m_nbPixelsAxial;        /*!< Total number of axial pixels as defined in the system file*/
    FLTNB m_pixelsSizeAxial;         /*!< Size of axial pixels as defined in the system file*/
    FLTNB m_gapSizeAxial;            /*!< Gap size between each axial pixel as defined in the system file*/

    FLTNB m_detectorDepth;           /*!< Depth of detection pixels*/
    FLTNB m_spotSizeWidth;           /*!< Width of the source, along the direction tangential to the scanner radius */
    FLTNB m_spotSizeDepth;           /*!< Depth of the source, along the direction of the scanner radius */

    // For the LUT
    FLTNB* mp_crystalCentralPositionX; /*!< Cartesian coordinate on X-axis of the center of each crystal, at each projection*/
    FLTNB* mp_crystalCentralPositionY; /*!< Cartesian coordinate on Y-axis of the center of each crystal, at each projection*/
    FLTNB* mp_crystalCentralPositionZ; /*!< Cartesian coordinate on Z-axis of the center of each crystal, at each projection*/

    FLTNB* mp_crystalOrientationX;     /*!< X-axis orientation of each crystal, at each projection*/
    FLTNB* mp_crystalOrientationY;     /*!< Y-axis orientation of each crystal, at each projection*/
    FLTNB* mp_crystalOrientationZ;     /*!< Z-axis orientation of each crystal, at each projection*/
 
    FLTNB* mp_sourcePositionX;   /*!< X-axis position of the focal point for each crystal, at each projection*/
    FLTNB* mp_sourcePositionY;   /*!< Y-axis position of the focal point for each crystal, at each projection*/
    FLTNB* mp_sourcePositionZ;   /*!< Z-axis position of the focal point for each crystal, at each projection*/
};

// Class for automatic insertion (set here the visible scanner type name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_SCANNER(CT,iScannerCT)

#endif
