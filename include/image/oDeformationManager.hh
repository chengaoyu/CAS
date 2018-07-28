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
  \brief    Declaration of class oDeformationManager
*/

#ifndef ODEFORMATIONMANAGER_HH
#define ODEFORMATIONMANAGER_HH 1

#include "gVariables.hh"
#include "vDeformation.hh"
#include "vDataFile.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/**
 * @defgroup DEF_TYPE Deformation type
 *
 *    \brief Nature of motion (respiratory, cardiac, patient motion, ...) \n
 *           Defined in oDeformationManager.hh
 * @{
 */
/** Constant corresponding to a deformation for respiratory motion correction (=0) */
#define DEF_RESP_MOT 0
/** Constant corresponding to a deformation for cardiac motion correction (=1) */
#define DEF_CARD_MOT 1
/** Constant corresponding to a deformation for involuntary patien motion correction (=2) */
#define DEF_IPAT_MOT 2
/** Constant corresponding to a deformation for dual respiratory and cardiac motion correction (=3) */
#define DEF_DUAL_MOT 3
/** @} */

class vDataFile;



/*!
  \class   oDeformationManager
  \brief   This class is designed to manage the image-based deformation part of the reconstruction
  \details As each manager class, it is created in the main program, all parameters are
           then set, checked, and the manager is initialized. \n
           The manager is then used by the algorithm itself, where the function PerformDeformation()
           is called each time the data belongs to a new "gate". \n
           The deformation functions are not multithreaded, therefore the threads are synchronized  
           (inside the main loop of the main algorithm) before performing any image transformation.
*/
class oDeformationManager
{
  // Constructor & Destructor
  public:
    /*!
      \fn      oDeformationManager::oDeformationManager
      \brief   Constructor of oDeformationManager. Simply set all data members to default values.
    */
    oDeformationManager();
    
    /*!
      \fn      oDeformationManager::~oDeformationManager
      \brief   Destructor of oDeformationManager. Free memory from all allocated tabs.
    */
    ~oDeformationManager();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      oDeformationManager::CheckParameters
      \brief   This function is used to check parameters after the latter
               have been all set using Set functions.
      \return  0 if success, positive value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      oDeformationManager::Initialize
      \brief   Set the flags for the different motion types and instanciate/initialize deformation objects
               through the ParseOptionsAndInitializeDeformations() private function.
      \return  0 if success, positive value otherwise.
    */
    int Initialize();
    /*!
      \fn      oDeformationManager::InstantiateImageForDeformation
      \param   oImageSpace* ap_Image : required to call oImageSpace instanciation functions
      \brief   If deformation is enabled, ask the Image Space to Instantiate the temporary backward image for deformation \n
               If reconstruction is in histogram mode, the temporal sensitivity image is instanciated as well
    */
    void InstantiateImageForDeformation(oImageSpace* ap_Image);
    /*!
      \fn      oDeformationManager::DeallocateImageForDeformation
      \param   oImageSpace* ap_Image : required to call oImageSpace deallocation functions
      \brief   If deformation is enabled, ask the Image Space to free memory of the temporary backward image for deformation \n
               If reconstruction is in histogram mode, the temporal sensitivity image is deallocated as well
    */ 
    void DeallocateImageForDeformation(oImageSpace* ap_Image);
    /*!
      \fn      oDeformationManager::InitImageForDeformation
      \param   oImageSpace* ap_Image : required to call oImageSpace initialization functions
      \brief   If deformation is enabled, ask the Image Space to initialize the temporary backward image for deformation \n
               If reconstruction is in histogram mode, the temporal sensitivity image is initialized as well
    */
    void InitImageForDeformation(oImageSpace* ap_Image);


  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      oDeformationManager::SetVerbose
      \param   a_verboseLevel
      \brief   Set the verbose level
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      oDeformationManager::SetDataMode
      \param   a_dataMode : (histogram/list-mode)
      \brief   Set the mode of reconstruction 
      \details required to enable/disable sensitivity image deformation in histogram mode
    */
    inline void SetDataMode(int a_dataMode)
           {m_dataMode = a_dataMode;}
    /*!
      \fn      oDeformationManager::SetImageDimensionsAndQuantification
      \param   ap_ImageDimensionsAndQuantification
      \brief   Set the image dimensions in use
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification) 
                                                   {mp_ID = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      oDeformationManager::SetOptions
      \param   a_options
      \param   a_nbGates
      \brief   Set the motion options contained in the provided string, and the related number of gates
    */
    inline void SetOptions(const string& a_options) 
           {m_options = a_options;}
           
    /*!
      \fn      oDeformationManager::SetNbTransformations
      \param   a_nbTransformations
      \brief   Set the total number of transformations/deformations
    */
    inline void SetNbTransformations(int a_nbTransformations) 
           {m_nbTransformations = a_nbTransformations;}

    /*!
      \fn      oDeformationManager::SetMotionType
      \param   a_motionType
      \brief   Set the nature of motion correction (Deformation type macro)
    */
    void SetMotionType(int a_motionType);
    
    /*!
      \fn      oDeformationManager::UseDeformationResp
      \brief   Indicate if the respiratory motion deformation is enabled
      \return  true if enabled, false otherwise
    */
    inline bool UseDeformationResp()
           {return m_UseDeformationResp;}
    /*!
      \fn      oDeformationManager::UseDeformationCard
      \brief   Indicate if the cardiac motion deformation is enabled
      \return  true if enabled, false otherwise
    */
    inline bool UseDeformationCard()
           {return m_UseDeformationCard;}
    /*!
      \fn      oDeformationManager::UseDeformationInv
      \brief   Indicate if the involuntary patient motion deformation is enabled
      \return  true if enabled, false otherwise
    */
    inline bool UseDeformationInv()
           {return m_UseDeformationIPat;}
    /*!
      \fn      oDeformationManager::GetNbSensImagesRespDeformation
      \param   a_value : default number of cardiac images
      \brief   return the required number of respiratory images in the sensitivity image depending on the respiratory deformation model
      \return  a number of images
    */
    inline int GetNbSensImagesRespDeformation(int a_value) 
           {if (UseDeformationResp()
            ||  UseDeformationInv() ) return 0;
            else                      return a_value;}
    /*!
      \fn      oDeformationManager::GetNbSensImagesCardDeformation
      \param   a_value : default number of cardiac images
      \brief   return the required number of cardiac images in the sensitivity image depending on the cardiac deformation model
      \return  a number of images
    */
    inline int GetNbSensImagesCardDeformation(int a_value) 
           {if (UseDeformationCard()) return 0; else return a_value;}


  // -------------------------------------------------------------------
    // Deformation functions
    /*!
      \fn      oDeformationManager::ApplyDeformationForSensitivityGeneration
      \param   oImageSpace* ap_Image : required to access oImageSpace image matrices
      \param   int a_defDirection : direction of the deformation (forward/backward)
      //\param   int a_defType : Nature of the motion (Respiratory/Cardiac/Involuntary Patient)
      \param   int idx
      \param   int fr
      \param   int rg
      \param   int cg
      \brief   Apply deformations during the list-mode sensitivity image generation 
      \details Perform deformation on the forward_image or the backward_image matrices corresponding to the current fr, rg, cg (if any), and depending on the defDirection.
      \todo    Some changes required if we merge respiratory/cardiac motion objects
      \todo    Check and implement patient motion
      \return  0 if success, positive value otherwise
    */
    //int ApplyDeformationForSensitivityGeneration(oImageSpace* ap_Image, int a_defDirection, int fr, int rg, int cg);
    int ApplyDeformationForSensitivityGeneration(oImageSpace* ap_Image, int a_defDirection, int idx, int fr, int rg, int cg);
    /*!
      \fn      oDeformationManager::ApplyDeformationsToBackwardImage
      \param   oImageSpace* ap_Image : required to access oImageSpace image matrices
      \brief   Apply final backward deformations on the backward image
      \details Call the eponym function for the deformation object, as well as >ApplyDeformationsToHistoSensitivityImage() if data mode is histogram. \n
               Then reinitialize the temporary backup deformation images (the backward image, and the sensitivity image if data mode is histogram) 
      \return  0 if success, positive value otherwise
    */
    int ApplyDeformationsToBackwardImage(oImageSpace* ap_Image);
    /*!
      \fn      oDeformationManager::PerformDeformation
      \param   oImageSpace* ap_Image : required to access oImageSpace image matrices
      \brief   Apply deformations during reconstruction
      \details Call the eponym function for the deformation object, 
               as well as PerformHistoSensitivityDeformation() if data mode is histogram.
      \todo    why the check on frames ?
      \return  0 if success, positive value otherwise
    */
    int PerformDeformation(oImageSpace* ap_Image);
    /*!
      \fn      oDeformationManager::TestDeformationOnImage
      \param   ap_inputImage : input image to deform
      \param   ap_outputImage : image in which the output of the deformation should be recovered
      \param   a_direction : a direction for the deformation to perform (forward or backward)
      \param   a_defIdx : index of the deformation
      \brief   Apply deformation specified by arguments on provided input image, for testing purposes
      \return  0 if success, positive value otherwise
    */
    int TestDeformationOnImage(FLTNB* ap_inputImage, FLTNB* ap_outputImage, int a_direction, int a_defIdx);

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      oDeformationManager::ParseOptionsAndInitializeDeformations
      \brief   Parse respiratory/cardiac/involuntary patient motion options contained in the previously provided
               strings. This function is called inside the Initialize() function.
      \details Manage the options reading and initialize specific vDeformation \n
               Options are a string containing first the name of the deformation,
               then either a ':' and a configuration file specific to the deformation
               - or - as many ',' as needed parameters for this deformation. \n
               Specific pure virtual functions of the vDeformation are used to read parameters and initialize them.
      \todo    Some cleaning if we merge respiratory and cardiac motion objects
      \return  0 if success, positive value otherwise
    */
    int ParseOptionsAndInitializeDeformations();


  // -------------------------------------------------------------------
  // Data members
  private:
    // Image dimensions
    oImageDimensionsAndQuantification* mp_ID; /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    // Options for each deformation type
    string m_options;                     /*!< The string containing options for the motion correction */

    // Deformation objects and associated bool
    vDeformation* mp_Deformation;         /*!< Deformation object for image-based motion */
    bool m_UseDeformationResp;                /*!< Flag indicating that transformation for respiratory motion is enabled */
    bool m_UseDeformationCard;                /*!< Flag indicating that transformation for cardiac motion is enabled */
    bool m_UseDeformationIPat;                /*!< Flag indicating that transformation for involuntary patient motion is enabled */
    
    // Variable indicating the current gate/index of the motion
    int m_curMotIdx;                          /*!< Current gate for the motion */
    int16_t* mp_curPatMotIdx;                 /*!< Current index for the ip motion by frame */

    // Number of gates for cyclic motion
    int m_nbTransformations;                  /*!< Number of image-based transformations/deformations */

    // Verbose level
    int m_verbose;                            /*!< The verbose level */
    // Data mode
    int m_dataMode;                           /*!< Data mode (list-mode (=0), histogram (=1)). Recovered from the datafile */
    // Has been checked ?
    bool m_checked;                           /*!< Boolean indicating whether the parameters were checked or not */
    // Has been initialized ?
    bool m_initialized;                       /*!< Boolean indicating whether the manager was initialized or not */
};

#endif
