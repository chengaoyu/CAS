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
  \brief    Declaration of class iScannerPET
*/

#ifndef ISCANNERPET_HH
#define ISCANNERPET_HH 1

#include "gVariables.hh"
#include "vScanner.hh"
#include "sAddonManager.hh"

/*!
  \class   iScannerPET
  \brief   This class is used to represent any cylindrical PET scanner
  \details Inherits from vScanner
*/
class iScannerPET : public vScanner
{
  // Constructor & Destructor
  public:
    /*!
      \brief   iScannerPET constructor. 
               Initialize the member variables to their default values.
    */
    iScannerPET();

    /*!
      \brief   iScannerPET destructor. 
    */
    ~iScannerPET();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_SCANNER(iScannerPET)
    /*!
      \fn      iScannerPET::DescribeSpecific()
      \brief   Implementation of the pure virtual eponym function that simply prints info about the scanner
    */
    void DescribeSpecific();
    /*!
      \fn      int iScannerPET::Instantiate()
      \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                                   the system is a generic file (0) or custom Look-up-table (1)
      \brief   Get mandatory informations from the scanner file 
               and allocate memory for the member variables
      \return  0 if success, positive value otherwise
    */
    int Instantiate(bool a_scannerFileIsLUT);
    /*!
      \fn      int iScannerPET::CheckParameters()
      \brief   Check if all parameters have been correctly initialized
      \return  0 if success, positive value otherwise
    */
    int CheckParameters();
    /*!
      \fn      int iScannerPET::Initialize()
      \brief   Check general initialization and set several parameters to their default value
      \return  0 if success, positive value otherwise
    */
    int Initialize();
    /*!
      \fn      int iScannerPET::BuildLUT()
      \param   a_scannerFileIsLUT : boolean indicating if the file describing 
                                    the SPECT camera is a generic file (0) or custom Look-up-table (1)
      \brief   Call the functions to generate the LUT or read the user-made LUT depending on the user choice
      \return  0 if success, positive value otherwise
      \todo    default values to max ring diff
    */
    int BuildLUT(bool a_scannerFileIsLUT);
    /*!
      \fn      int iScannerPET::GetPositionsAndOrientations()
      \param   a_index1 : index of the 1st crystal 
      \param   a_index2 : index of the 2nd crystal 
      \param   ap_Position1[3] : x,y,z cartesian position of center of the 1st crystal
      \param   ap_Position2[3] : x,y,z cartesian position of center of the 2nd crystal
      \param   ap_Orientation1[3] : x,y,z components of the orientation vector related to the 1st crystal
      \param   ap_Orientation2[3] : x,y,z components of the orientation vector related to the 2nd crystal
      \param   ap_POI1 : x,y,z components of the Point Of Interation related to the 1st crystal
      \param   ap_POI2 : x,y,z components of the Point Of Interation related to the 2nd crystal
      \brief   Get the central positions and orientations of the scanner elements from their indices.
      \details This method is very general and is used by the vProjector. 
               From these positions and orientations, other methods can be used by specific projectors to get specific positions.
               Position calculations include POI and mean depth of interaction
      \todo    some cases depending on POI are not implemented
      \return  0 if success, positive value otherwise
    */
    int GetPositionsAndOrientations( int a_index1, int a_index2,
                                     FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                     FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3],
                                     FLTNB* ap_POI1 = NULL, FLTNB* ap_POI2 = NULL );
    /*!
      \fn      int iScannerPET::GetRdmPositionsAndOrientations()
      \param   a_index1 : index of the 1st crystal 
      \param   a_index2 : index of the 2nd crystal 
      \param   ap_Position1[3] : x,y,z cartesian position of center of the 1st crystal
      \param   ap_Position2[3] : x,y,z cartesian position of center of the 2nd crystal
      \param   ap_Orientation1[3] : x,y,z components of the orientation vector related to the 1st crystal
      \param   ap_Orientation2[3] : x,y,z components of the orientation vector related to the 2nd crystal
      \brief   Get random positions of the scanner elements and their orientations from their indices.
      \details - Find the scanner elements described by the two indexes passed in parameters. 
      \details - Compute random positions inside the elements described by the indexes passed in parameters
      \details - Find the scanner elements described by the two indexes passed in parameters.
      \details - Write the corresponding random cartesian coordinates in the positions parameters.
      \details Position calculations include POI and mean depth of interaction
      \todo    fix the possibility to draw LOR outside the actual crystal volume (if mp_meanDepthOfInteraction != 0.5)
      \todo    some cases depending on POI are not implemented
      \return  0 if success, positive value otherwise
    */
    int GetRdmPositionsAndOrientations( int a_index1, int a_index2,
                                        FLTNB ap_Position1[3], FLTNB ap_Position2[3],
                                        FLTNB ap_Orientation1[3], FLTNB ap_Orientation2[3] );
    /*!
      \fn      int iScannerPET::GetPositionWithRandomDepth()
      \param   a_index1 : index of the 1st crystal
      \param   a_index2 : index of the 2nd crystal
      \param   ap_Position1[3] : x,y,z cartesian position of the point related to the 1st index (see child function for more details)
      \param   ap_Position2[3] : x,y,z cartesian position of the point related to the 2st index (see child function for more details)
      \brief   Get the positions and orientations of scanner elements from their indices, with a random depth.
      \details Method for testing purposes. Does not include POI and mean depth of interaction
      \return  0 if success, positive value otherwise
    */
    int GetPositionWithRandomDepth( int a_index1, int a_index2, FLTNB ap_Position1[3], FLTNB ap_Position2[3] );
    /*!
      \fn      int iScannerPET::GetTwoCorners()
      \param   a_index1 : index of the 1st crystal 
      \param   a_index2 : index of the 2nd crystal
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
      \fn      int iScannerPET::GetEdgesCenterPositions()
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
               Get the cartesian coordinaters of the center of the 4 edges of the detection elements. \n
               It is typically used for the Distance Driven projector
      \return  0 if success, positive value otherwise
    */
    int GetEdgesCenterPositions( int a_index1, int a_index2,
      FLTNB ap_pos_line_point1[3], FLTNB ap_pos_line_point2[3],
      FLTNB ap_pos_point1_x[4], FLTNB ap_pos_point1_y[4], FLTNB ap_pos_point1_z[4],
      FLTNB ap_pos_point2_x[4], FLTNB ap_pos_point2_y[4], FLTNB ap_pos_point2_z[4] );
    /*!
      \fn      inline int iScannerPET::GetSystemNbElts()
      \return  the number of crystals in the PET system
    */
    inline int GetSystemNbElts()
           {return m_nbCrystals;};

    /*!
      \fn      int iScannerPET::IsAvailableLOR()
      \param   a_elt1 : index of the 1st crystal
      \param   a_elt2 : index of the 2nd crystal
      \brief   Check if the LOR formed by the crystalf whose indices are passed 
               in parameters is available according to the scanner restrictions
      \details This function is dedicated to analytic projection and list-mode sensitivity image generation \n
               The PET implementation checks the LOR availability according to the minimal (transaxial) angle difference 
               and the maximal ring difference between two crystals
      \todo    min angle difference (currently just system dependent) may be estimated from the FOV requested by the user ?
      \todo    Restriction for ring_ID NOT safe for user using their own LUT !!! \n
      #         Perhaps throw warning in this case
      \return  1 if the LOR is available, 0 otherwise
    */
    int IsAvailableLOR(int a_elt1, int a_elt2);
    /*!
      \fn      int iScannerPET::GetGeometricInfoFromDataFile()
      \brief   Retrieve PET geometric informations from the datafile
      \details Recover the (axial) max ring difference
      \return  0 if success, positive value otherwise
    */
    int GetGeometricInfoFromDataFile(string a_pathToDataFilename);
    /*!
      \fn      inline int iScannerPET::SetPETMaxAxialDiff()
      \param   a_maxAxialDiffmm
      \brief   Set the maximal axial difference in mm between 2 crystals forming a lor (if any)
    */
    inline int SetPETMaxAxialDiffmm(FLTNB a_maxAxialDiffmm)
           {m_maxAxialDiffmm = a_maxAxialDiffmm; return 0;}
    /*!
      \fn      void iScannerPET::ShowHelp()
      \brief   Display help
      \todo    Provide informations about PET system initialization ?
    */
    void ShowHelp();
    /*!
      \fn      inline FLTNB iScannerPET::GetDetectionElementSizeTrans()
      \return  return the transaxial size of the detection element
      \todo    take the multi layer into account
    */
    inline FLTNB GetDetectionElementSizeTrans()
           {return mp_sizeCrystalTrans[0];}
    /*!
      \fn      inline FLTNB iScannerPET::GetDetectionElementSizeAxial()
      \return  return the axial size of the detection element
      \todo    take the multi layer into account
    */
    inline FLTNB GetDetectionElementSizeAxial()
           {return mp_sizeCrystalAxial[0];}

  // -------------------------------------------------------------------
  // Functions dedicated to Analytic Projection
    /*!
      \fn      int iScannerPET::PROJ_GetPETSpecificParameters()
      \param   ap_maxAxialDiffmm
      \brief   Set pointers passed in argument with the related PET specific variables
               This function is used to recover these values in the datafile object
      \return  0 if success, positive value otherwise
    */
    int PROJ_GetPETSpecificParameters(FLTNB* ap_maxAxialDiffmm);


  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      int iScannerPET::LoadLUT()
      \brief   Load a precomputed scanner LUT
      \details Read mandatory data from the header of the LUT. \n
               Then load the LUT elements for each crystal
      \return  0 if success, positive value otherwise
    */
    int LoadLUT();
    /*!
      \fn      int iScannerPET::ComputeLUT()
      \brief   Compute the LUT of the scanner from a generic (.geom) file
      \details Read mandatory data from the geom file. Then compute the LUT elements for each crystal from the geometry described in the file
      \details Compute the look-up-tables of the system containing the locations of the scanner elements center in space and their orientations
      \todo    Add option to inverse rsector if transaxial orientation is counter-clockwise. ?
      \todo    rotation for non-perfectly cylindric scanner (angles along the Z axis)
      \return  0 if success, positive value otherwise
    */
    int ComputeLUT();
    /*!
      \fn      int iScannerPET::GetLayer()
      \param   a_idx : index of the crystal in the loaded LUT
      \brief   Get the layer from which the 'a_index' crystal belongs to
      \return  layer index
    */
    int GetLayer(int a_idx);


  // -------------------------------------------------------------------
  // Data members
  private:
    int m_nbLayers;                    /*!< Number of crystal layers in the scanner */
    int m_nbCrystals;                  /*!< Total number of crystal in the scanner */
    int* mp_nbCrystalsInLayer;         /*!< Number of crystals (for each layer) */
      
    FLTNB* mp_crystalCentralPositionX; /*!< Cartesian coordinate on X-axis of the center of each crystal*/
    FLTNB* mp_crystalCentralPositionY; /*!< Cartesian coordinate on Y-axis of the center of each crystal*/
    FLTNB* mp_crystalCentralPositionZ; /*!< Cartesian coordinate on Z-axis of the center of each crystal*/

    FLTNB* mp_crystalOrientationX;     /*!< X-axis orientation of each crystal*/
    FLTNB* mp_crystalOrientationY;     /*!< Y-axis orientation of each crystal*/
    FLTNB* mp_crystalOrientationZ;     /*!< Z-axis orientation of each crystal*/

    FLTNB* mp_sizeCrystalTrans;        /*!< Transaxial size of a crystal (for each layer)*/
    FLTNB* mp_sizeCrystalAxial;        /*!< Axial size of a crystal (for each layer)*/
    FLTNB* mp_sizeCrystalDepth;        /*!< Depth of a crystal (for each layer)*/
    FLTNB* mp_meanDepthOfInteraction;  /*!< Mean depth of interaction in a crystal (for each layer)*/
    
    FLTNB m_minAngleDifference;        /*!< Minimal transaxial angle difference for the system*/

    FLTNB m_maxAxialDiffmm;            /*!< Maximal axial difference in mm for a specific acquisition*/
};


// Class for automatic insertion (set here the visible scanner type name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_SCANNER(PET,iScannerPET)

#endif
