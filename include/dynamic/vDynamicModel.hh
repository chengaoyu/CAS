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
  \ingroup  dynamic
  \brief    Declaration of class vDynamicModel
*/

#ifndef VDYNAMICMODEL_HH
#define VDYNAMICMODEL_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oOptimizerManager.hh"

class oImageSpace;

/*!
  \class   vDynamicModel
  \brief   This is the mother class of dynamic model classes
  \details This class is a virtual one, in the sense that it cannot be used on its own \n
           because several pure virtual functions belong to it. 
           Its children are implementations of actual dynamic models. \n
           Everywhere in the code, this parent class should be used instead of any of its children. \n
           It can be used during the reconstruction process by the oDynamicModelManager through the
           use of the EstimateModelParameters() and EstimateImageWithModel() functions
           
           All children must implement the following pure virtual functions: \n
            - ReadAndCheckConfigurationFile(): read specific options from a configuration file \n
            - ReadAndCheckOptionsList(): read specific options from a string \n
            - ShowHelp(): print helps about the projector specifications \n
            - Initialize(): initialize specific data of the projector (if required) \n
            - CheckParameters(): Check the initialization of the parameters (if required) \n
            
            - EstimateModelParameters() : (virtual only)  \n
                                          Estimate any temporal functions or coefficients  \n
                                          related to the dynamic model (if required) \n
            - EstimateImageWithModel(): Fit the dynamic model to the series of dynamic images
*/
class vDynamicModel
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      vDynamicModel::vDynamicModel
      \brief   Constructor of vDynamicModel. Simply set all data members to default values.
    */
    vDynamicModel();
    /*!
      \fn      vDynamicModel::~vDynamicModel
      \brief   Destructor of vDynamicModel.
    */
    virtual ~vDynamicModel();


  // -----------------------------------------------------------------------------------------
  // Public member functions related to the initialization of the model
  public:
    /*!
      \fn      vDynamicModel::SetImageDimensionsAndQuantification
      \param   ap_ImageDimensionsAndQuantification
      \brief   Set the image dimensions in use
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
                           {mp_ID = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      vDynamicModel::SetVerbose
      \param   a_verboseLevel
      \brief   Set the verbose level
    */
    inline void SetVerbose(int a_verbose) 
               {m_verbose = a_verbose;}
    /*!
      \fn      vDynamicModel::CheckParameters
      \brief   This function is used to check parameters after the latter
               have been all set using Set functions.
      \return  0 if success, positive value otherwise.
    */
    virtual int CheckParameters();
    /*!
      \fn      vDynamicModel::CheckSpecificParameters
      \brief   This function is used to check the parameters of the child functions before initialization if required.
      \details It could be overloaded by the child if needed. Default implementation is empty and return 0.
      \return  0 if success, other value otherwise.
    */
    virtual int CheckSpecificParameters() = 0;
    /*!
      \fn      vDynamicModel::ReadAndCheckConfigurationFile
      \param   const string& a_configurationFile : ASCII file containing informations about a dynamic model
      \brief   This function is used to read options from a configuration file. \n
               It is pure virtual so must be implemented by children.
      \return  0 if success, other value otherwise.
    */
    virtual int ReadAndCheckConfigurationFile(string a_fileOptions) = 0;
    /*!
      \fn      vDynamicModel::ReadAndCheckOptionsList
      \param   const string& a_optionsList : a list of parameters separated by commas
      \brief   This function is used to read parameters from a string. \n
               It is pure virtual so must be implemented by children.
      \return  0 if success, other value otherwise.
    */
    virtual int ReadAndCheckOptionsList(string a_listOptions) = 0;
    /*!
      \fn      vDynamicModel::Initialize
      \brief   This function is used to initialize specific data related to the child deformation model. \n
               It is pure virtual so must be implemented by children.
      \return  0 if success, other value otherwise.
    */
    virtual int Initialize() = 0;
    /*!
      \fn      vDynamicModel::ShowHelp
      \brief   This function is used to print out specific help about the deformation and its options. It is
               pure virtual so must be implemented by children.
    */
    virtual void ShowHelp() = 0;


  // -----------------------------------------------------------------------------------------
  // Public member functions called by the main iterative algorithm class
  public:
    /*!
      \fn      vDynamicModel::EstimateModelParameters
      \param   ap_ImageS : pointer to the ImageSpace
      \param   a_ite : index of the actual iteration
      \param   a_sset : index of the actual subset
      \brief   This function is pure virtual so must be implemented by children. \n
               It can be used to estimate any temporal functions or coefficients 
               related to the dynamic model, if required
      \return  0 if success, other value otherwise.
    */
    virtual int EstimateModelParameters(oImageSpace* ap_Image, int a_ite, int a_sset) = 0;
    /*!
      \fn      vDynamicModel::EstimateImageWithModel
      \param   ap_ImageS : pointer to the ImageSpace
      \param   a_ite : index of the actual iteration
      \param   a_sset : index of the actual subset
      \brief   This function is pure virtual so must be implemented by children. \n
               It is used to fit the dynamic model to the series of dynamic images
      \return  0 if success, other value otherwise.
    */
    virtual int EstimateImageWithModel(oImageSpace* ap_Image, int a_ite, int a_sset) = 0;
    /*!
      \fn      vDynamicModel::SaveParametricImages
      \param   a_iteration : current iteration index
      \param   a_subset : current number of subsets (or -1 by default)
      \brief   This function is pure virtual so must be implemented by children \n
               Call SaveParametricImages() function of the dynamic model object,
               in order to write on disk any parametric image of the model
      \return  0 if success, positive value otherwise
    */
    virtual int SaveParametricImages(int a_iteration, int a_subset = -1);
    /*!
      \fn      vDynamicModel::ApplyOutputFOVMaskingOnParametricImages
      \brief   Mask the outside of the transaxial FOV based on the m_fovOutPercent
      \details Similar to the eponym function in ImageSpace, but on parametric images
    */
    virtual int ApplyOutputFOVMaskingOnParametricImages();
    /*!
      \fn      vDynamicModel::ComputeOutputParImage
      \brief   Compute output image using the m2p_parametricImages matrix Store the result in the m2p_outputParImages matrix
    */
    virtual void ComputeOutputParImage();

  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    oImageDimensionsAndQuantification* mp_ID; /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    int m_verbose;                            /*!< The verbose level */
    int m_nbTimeBF;                           /*!< Number of time basis functions in the model*/
    bool m_checked;                           /*!< Boolean indicating whether the parameters were checked or not */
    bool m_initialized;                       /*!< Boolean indicating whether the manager was initialized or not */
    bool m_saveParImageFlag;                  /*!<Flag indicating if parametric images should be written on disk (default=true) */
    uint16_t m_nbModelParam;                         /*!<Nb parameters in the model */
    FLTNB** m2p_parametricImages; /*!< Image matrix containing the parametric images \n
                                       2 pointers:  \n
                                       1: Parametric image related to the dynamic model basis functions. \n
                                       2: 3D voxels */
                                       
    FLTNB** m2p_modelTACs;       /*!< Vector containing the Model temporal basis functions \n
                                       2 pointers:  \n
                                       1: index of the temporal function \n
                                       2: coefficient of the functions for each time points of a dynamic acquisition */

    FLTNB** m2p_outputParImages; /*!< Image matrix to gather the parametric image before writing on disk \n
                                      By default it will point directly to the parametric m2p_parametricImages. \n
                                      They are allocated if post-processing are enabled before writing the image (i.e FOV masking) \n
                                       2 pointers:  \n
                                       1: Parametric image related to the dynamic model basis functions. \n
                                       2: 3D voxels */
                                  
};

// ----------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ----------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_DYNAMICMODEL(CLASS) \
  static vDynamicModel *make_dynamic_model() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_DYNAMICMODEL(NAME,CLASS)                                                               \
  class NAME##DynamicModelCreator                                                                    \
  {                                                                                                  \
    public:                                                                                          \
      NAME##DynamicModelCreator()                                                                    \
        { sAddonManager::GetInstance()->mp_listOfDynamicModels[#NAME] = CLASS::make_dynamic_model; } \
  };                                                                                                 \
  static NAME##DynamicModelCreator DynamicModelCreator##NAME;

#endif
