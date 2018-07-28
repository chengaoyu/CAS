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
  \ingroup  image
  \brief    Declaration of class vDeformation
*/

#ifndef VDEFORMATION_HH
#define VDEFORMATION_HH 1

#include "gVariables.hh"

// These definitions are used to simply discriminate between forward and backward deformation
#define FORWARD_DEFORMATION  0
#define BACKWARD_DEFORMATION 1

class oImageSpace;
class oImageDimensionsAndQuantification;

/*!
  \class   vDeformation
  \brief   This is the mother class of image-based transformation class
  \details This class is a virtual one, in the sense that it cannot be used on its own
           because several pure virtual functions belong to it. Its children are
           implementations of actual deformation models. \n
           Everywhere in the code, this parent class should be used instead of any of its children. \n
           It can be used during the projection/reconstruction process by the oDeformationManager through the
           use of the Deformation functions that cannot be overloaded: \n
           - ApplyDeformationsToForwardImage() : Initial deformation of the forward image matrice \n
           - PerformDeformation() : Deformation of the forward image during the reconstruction process \n
           - ApplyDeformationsToBackwardImage() : Final deformation of the backward image \n
           - PerformSensitivityDeformation() : Forward or backward deformation during the sensitivity generation process \n
           
           All children must implement the following pure virtual functions: \n
            - ReadAndCheckConfigurationFile(): read specific options from a configuration file \n
            - ReadAndCheckOptionsList(): read specific options from a string \n
            - ShowHelp(): print helps about the deformation model specifications \n
            - Initialize(): initialize specific stuff of the model (if required) \n
            - CheckParameters(): Check the initialization of the parameters (if required) \n
            - ApplyDeformations(): implement the actual deformation
*/
class vDeformation
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      vDeformation::vDeformation
      \brief   Constructor of vDeformation. Simply set all data members to default values.
    */
    vDeformation();
    /*!
      \fn      vDeformation::~vDeformation
      \brief   Destructor of vDeformation.
    */
    virtual ~vDeformation();


  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      inline void vDeformation::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the image dimensions in use
    */
    inline void SetImageDimensionsAndQuantification( oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification )
           {mp_ID = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      inline void vDeformation::SetVerbose()
      \param   a_verboseLevel
      \brief   Set the verbose level
    */
    inline void SetVerbose( int a_verbose )
           {m_verbose = a_verbose;}
    /*!
      \fn      inline void vDeformation::SetNbTransformations()
      \param   a_nbTransformations
      \brief   Set the number of transformation in the data to be performed on the dataset
               (equal to the number of gates in one frame)
    */
    inline void SetNbTransformations( int a_nbTransformations )
           {m_nbTransformations = a_nbTransformations;}
    /*!
      \fn      virtual int vDeformation::CheckParameters()
      \brief   This function is used to check parameters after the latter
               have been all set using Set functions.
      \return  0 if success, positive value otherwise.
    */
    virtual int CheckParameters();
    /*!
      \fn      virtual int vDeformation::CheckSpecificParameters() = 0
      \brief   This function is used to check the parameters of the child functions before initialization if required.
      \details It could be overloaded by the child if needed. Default implementation is empty and return 0.
      \return  0 if success, other value otherwise.
    */
    virtual int CheckSpecificParameters() = 0;
    /*!
      \fn      virtual int vDeformation::ReadAndCheckConfigurationFile() = 0
      \param   const string& a_configurationFile : ASCII file containing informations about a deformation model
      \brief   This function is used to read options from a configuration file. \n
               It is pure virtual so must be implemented by children.
      \return  0 if success, other value otherwise.
    */
    virtual int ReadAndCheckConfigurationFile( const string& a_fileOptions ) = 0;
    /*!
      \fn      virtual int vDeformation::ReadAndCheckOptionsList() = 0
      \param   const string& a_optionsList : a list of parameters separated by commas
      \brief   This function is used to read parameters from a string. \n
               It is pure virtual so must be implemented by children.
      \return  0 if success, other value otherwise.
    */
    virtual int ReadAndCheckOptionsList( const string& a_listOptions ) = 0;
    /*!
      \fn      virtual int vDeformation::Initialize() = 0
      \brief   This function is used to initialize specific data related to the child deformation model. \n
               It is pure virtual so must be implemented by children.
      \return  0 if success, other value otherwise.
    */
    virtual int Initialize() = 0;
    /*!
      \fn      virtual void vDeformation::ShowHelp() = 0
      \brief   This function is used to print out specific help about the deformation and its options. \n
               It is pure virtual so must be implemented by children.
    */
    virtual void ShowHelp() = 0;


  // -----------------------------------------------------------------------------------------
  // Reconstruction deformation function
  public:
    /*!
      \fn      virtual int vDeformation::PerformDeformation()
      \param   ap_Image : required to access oImageSpace image matrices
      \param   a_defIdx : index of the deformation
      \param   a_fr : frame index
      \param   a_rimg : respiratory image index
      \param   a_cimg : cardiac image index
      \brief   Apply deformations during reconstruction
      \details 1. Recover all the data of the multithreaded backward image matrices in the first one (thread index 0) \n
               2. Perform backward deformation of the backward image to the reference position with defIdx-1 \n
               3. Add coefficients of the backward image matrice to the temporary backup image matrice & reset backward image \n
               4. Apply forward deformation of the forward image (backward deformation to the reference position with defIdx-1, and forward deformation with defIdx)
      \return  0 if success, positive value otherwise
    */
    virtual int PerformDeformation( oImageSpace* ap_Image, int a_defIdx, int a_fr, int a_rimg, int a_cimg );
    /*!
      \fn      virtual int vDeformation::PerformHistoSensitivityDeformation()
      \param   ap_Image : required to access oImageSpace image matrices
      \param   a_defIdx : index of the deformation
      \param   fr : frame index
      \param   rimg : respiratory image index
      \param   cimg : cardiac image index
      \brief   Apply deformations on the sensitivity image during reconstruction in histogram mode
      \details 1. Recover all the data of the multithreaded sensitivity image matrice in the first one (thread index 0) \n
               2. Perform backward deformation of the sensitivity image to the reference position with defIdx-1 \n
               3. Add coefficients of the sensitivity image matrice to the temporary backup image matrice & reset sensitivity image
      \return  0 if success, positive value otherwise
    */
    
    virtual int PerformDeformationBis(oImageSpace* ap_Image, int a_defIdx, int a_fr, int a_rimg, int a_cimg);
    
    virtual int PerformHistoSensitivityDeformation( oImageSpace* ap_Image, int a_defIdx, int fr, int rimg, int cimg );
    /*!
      \fn      virtual int vDeformation::ApplyDeformationsToBackwardImage()
      \param   ap_Image : required to access the backward image and its deformation backup matrice
      \param   a_fr : frame index
      \param   a_defIdx : index of the deformation
      \brief   Apply backward transformation of the backward image to the reference position
      \details Loop on frames \n
               Recover any potential data stored in the backup matrice m2p_defTmpBackwardImage
      \todo    Bed management for patient motion correction
      \return  0 if success, positive value otherwise
    */
    virtual int ApplyDeformationsToBackwardImage( oImageSpace* ap_Image, int a_fr, int a_defIdx );
    /*!
      \fn      virtual int vDeformation::ApplyDeformationsToHistoSensitivityImage()
      \param   ap_Image : required to access the backward image and its deformation backup matrice
      \param   a_fr : frame index
      \param   a_defIdx : index of the deformation
      \brief   Apply backward transformations of the sensitivity image to the reference position (histogram mode)
      \details Loop on frames \n
               Recover any potential data stored in the backup matrice m4p_defTmpSensitivityImage
      \todo    Bed management for patient motion correction
      \return  0 if success, positive value otherwise
    */
    virtual int ApplyDeformationsToHistoSensitivityImage( oImageSpace* ap_Image, int a_fr, int a_defIdx );
    /*!
      \fn      virtual int vDeformation::PerformSensitivityDeformation()
      \param   ap_Image : required to access oImageSpace image matrices
      \param   a_defDirection : a direction for the deformation to perform (forward or backward)
      \param   a_defIdx : index of the deformation
      \param   fr : frame index
      \param   rg : respiratory gate index
      \param   cg : cardiac gate index
      \brief   Apply image deformations during sensitivity image generation for list-mode
      \details Depending on the deformation direction (forward or backward): \n
               Forward : Perform forward deformation of the forward image to the deformation index position \n
               Backward: Perform backward deformation of the backward image to the reference position
      \return  0 if success, positive value otherwise
    */
    virtual int PerformSensitivityDeformation( oImageSpace* ap_Image, int a_defDirection, int a_defIdx, int fr, int rg, int cg );
    /*!
      \fn      virtual int vDeformation::ApplyDeformations() = 0
      \param   ap_inputImage : input image to deform
      \param   ap_outputImage : image in which the output of the deformation should be recovered
      \param   a_direction : a direction for the deformation to perform (forward or backward)
      \param   a_defIdx : index of the deformation
      \brief   This function prepares the deformation to perform  \n
               It is a virtual pure deformation function to be implemented in the child class
      \return  0 if success, other value otherwise.
    */
    virtual int ApplyDeformations( FLTNB* ap_inputImage, FLTNB* ap_outputImage, int a_direction, int a_defIdx ) = 0;
    
    
    /*
      \fn Tlerp
      \param ap_inputImage : input image matrix 
      \param ap_outputImage : output image matrix
      \param iov : index of the voxel to interpolate in the output image
      \param iiv : index of the input image central voxel for interpolation
      \param dx : x-axis difference between output voxel cartesian position after transformation and center of estimated vox position
      \param dy : y-axis difference between output voxel cartesian position after transformation and center of estimated vox position
      \param dz : z-axis difference between output voxel cartesian position after transformation and center of estimated vox position
      \brief This function performs a trilinear interpolation for a specific voxel
      \todo : perhaps use padded image in order to avoid if statements
      \return 0 if success, other value otherwise.
    */
    int Tlerp(HPFLTNB *ap_inputImage, HPFLTNB *ap_outputImage, uint32_t iov, uint32_t iiv, FLTNB dX, FLTNB dY, FLTNB dZ);
    
    
    
    
  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    oImageDimensionsAndQuantification* mp_ID; /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    int m_verbose;                            /*!< The verbose level */
    int m_nbTransformations;                  /*!< Number of transformations to be performed on the dataset (corresponding to the number of gates subsets of the data)*/
    int m_checked;                            /*!< Boolean indicating whether the parameters were checked or not */
    int m_initialized;                        /*!< Boolean indicating whether the manager was initialized or not */
};


// ----------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ----------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_DEFORMATION(CLASS) \
  static vDeformation *make_deformation() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_DEFORMATION(NAME,CLASS)                                                             \
  class NAME##DeformationCreator                                                                  \
  {                                                                                               \
    public:                                                                                       \
      NAME##DeformationCreator()                                                                  \
        { sAddonManager::GetInstance()->mp_listOfDeformations[#NAME] = CLASS::make_deformation; } \
  };                                                                                              \
  static NAME##DeformationCreator DeformationCreator##NAME;

#endif
