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
  \brief    Declaration of class iDynamicModelTemplate
*/

#ifndef IPATLAKMODEL_HH
#define IPATLAKMODEL_HH 1

// Mother class of dynamic model
#include "vDynamicModel.hh"
// Required to automatically add the class in the CASToR code
#include "sAddonManager.hh"

/*!
  \class   iDynamicModelTemplate
  \brief   This class is a child of the vDynamicModel class implementing a template squeleton
  \details Use this class to implement your own custom deformation model.
*/
class iDynamicModelTemplate : public vDynamicModel
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      iDynamicModelTemplate::iDynamicModelTemplate
      \brief   Constructor of iDynamicModelTemplate. Simply set all data members to default values.
    */
    iDynamicModelTemplate();
    /*!
      \fn      iDynamicModelTemplate::~iDynamicModelTemplate
      \brief   Destructor of iDynamicModelTemplate
    */
    ~iDynamicModelTemplate();


  // -----------------------------------------------------------------------------------------
  // Public member functions related to the initialization of the model
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_DYNAMICMODEL(iDynamicModelTemplate)
    /*!
      \fn      iDynamicModelTemplate::CheckSpecificParameters
      \brief   This function is used to check whether all member variables
               have been correctly initialized or not.
      \return  0 if success, positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iDynamicModelTemplate::ReadAndCheckConfigurationFile
      \param   a_configurationFile : ASCII file containing informations about a dynamic model
      \brief   This function is used to read options from a configuration file.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckConfigurationFile(string a_fileOptions);
    /*!
      \fn      iDynamicModelTemplate::ReadAndCheckOptionsList
      \param   a_optionsList : a list of parameters separated by commas
      \brief   This function is used to read parameters from a string.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckOptionsList(string a_listOptions);
    /*!
      \fn      iDynamicModelTemplate::Initialize
      \brief   This function is used to initialize the model parametric images and basis functions
      \todo    Read Interfile for parametric images initialization
      \return  0 if success, other value otherwise.
    */
    int Initialize();
    /*!
      \fn      iDynamicModelTemplate::ShowHelp
      \brief   This function is used to print out specific help about the deformation model and its options.
    */
    void ShowHelp();


  // -----------------------------------------------------------------------------------------
  // Public member functions called by the main iterative algorithm class
    /*!
      \fn      iDynamicModelTemplate::EstimateModelParameters
      \param   ap_ImageS : pointer to the ImageSpace
      \param   a_ite : index of the actual iteration (not used)
      \param   a_sset : index of the actual subset (not used)
      \brief   Estimate the model parametric images
      \return  0 if success, other value otherwise.
    */
    int EstimateModelParameters(oImageSpace* ap_Image, int a_ite, int a_sset);
    /*!
      \fn      iDynamicModelTemplate::EstimateImageWithModel
      \param   ap_ImageS : pointer to the ImageSpace
      \param   a_ite : index of the actual iteration (not used)
      \param   a_sset : index of the actual subset (not used)
      \brief   Estimate image using model parametric images and basis functions
      \return  0 if success, other value otherwise.
    */
    int EstimateImageWithModel(oImageSpace* ap_Image, int a_ite, int a_sset);
    /*!
      \fn      iDynamicModelTemplate::SaveParametricImages
      \param   a_ite : index of the actual iteration
      \brief   Write parametric images on disk
      \return  0 if success, other value otherwise.
    */
    int SaveParametricImages(int a_iteration, int a_subset = -1);


  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    FLTNB** m2p_parametricImages; /*!< 2 dimensions (or more if required) array containing parametric image. 
                                       It must be a 2 dimensions array [nb of basis functions][nb of voxels] */
};

// Class for automatic insertion (set here the visible dynamic model's name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_DYNAMICMODEL(template,iDynamicModelTemplate)

#endif
