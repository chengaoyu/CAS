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
  \brief    Declaration of class iDeformationTemplate
*/

#ifndef IDEFORMATIONTEMPLATE_HH
#define IDEFORMATIONTEMPLATE_HH 1

// (Optional) may be required if one want access to the image matrices
#include "oImageSpace.hh"
// Required to automatically add the class in the CASToR code
#include "sAddonManager.hh"
// Mother class of deformation model
#include "vDeformation.hh"

/*!
  \class   iDeformationTemplate
  \brief   This class is a child of the vDeformation class implementing a template squeleton
  \details Use this class to implement your own custom deformation model.
*/
class iDeformationTemplate : public vDeformation
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      iDeformationTemplate::iDeformationTemplate
      \brief   Constructor of iDeformationTemplate. Simply set all data members to default values.
    */
    iDeformationTemplate();
    /*!
      \fn      iDeformationTemplate::~iDeformationTemplate
      \brief   Destructor of iDeformationTemplate. Free memory from all allocated tabs.
    */
    ~iDeformationTemplate();


  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_DEFORMATION(iDeformationTemplate)
    /*!
      \fn      iDeformationTemplate::ReadAndCheckConfigurationFile
      \param   a_configurationFile : ASCII file containing informations about a dynamic model
      \brief   This function is an implementation of the pure virtual mother function. It is used to read options from a configuration file.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckConfigurationFile(const string& a_fileOptions);
    /*!
      \fn      iDeformationTemplate::ReadAndCheckOptionsList
      \param   a_optionsList : a list of parameters separated by commas
      \brief   This function is an implementation of the pure virtual mother function. It is used to read options from a list of options.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckOptionsList(const string& a_listOptions);
    /*!
      \fn      iDeformationTemplate::CheckSpecificParameters
      \brief   This function is an implementation of the pure virtual mother function. It is used to
               check parameters of the child deformation model before initialization.
      \return  0 if success, other value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iDeformationTemplate::Initialize
      \brief   This function is an implementation of the pure virtual mother function. It is used to
               initialize specific stuff to the child deformation model.
      \return  0 if success, other value otherwise.
    */
    int Initialize();
    /*!
      \fn      iDeformationTemplate::ShowHelp
      \brief   This function is used to print out specific help about the deformation model and its options.
    */
    void ShowHelp();
    /*!
      \fn      iDeformationTemplate::ApplyDeformations
      \param   ap_inputImage : input image to deform
      \param   ap_outputImage : image in which the output of the deformation should be recovered
      \param   a_direction : a direction for the deformation to perform (forward or backward)
      \param   a_defIdx : index of the deformation
      \brief   This function is an implementation of the pure virtual mother function. The actual deformation should be implemented here 
      \return  0 if success, other value otherwise.
    */
    int ApplyDeformations(FLTNB* ap_inputImage, FLTNB* ap_outputImage, int a_direction, int a_defIdx);


  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:


  // -----------------------------------------------------------------------------------------
  // Data members
  private:
  
};

// Class for automatic insertion (set here the visible image deformation's name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_DEFORMATION(template,iDeformationTemplate)

#endif
