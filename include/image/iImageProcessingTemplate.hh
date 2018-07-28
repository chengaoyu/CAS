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
  \brief    Declaration of class iImageProcessingTemplate
*/

#ifndef IIMAGEPROCESSINGTEMPLATE_HH
#define IIMAGEPROCESSINGTEMPLATE_HH 1

#include "gVariables.hh"
#include "gOptions.hh"
#include "oImageSpace.hh"
#include "vImageProcessingModule.hh"
#include "sAddonManager.hh"

/*!
  \class   iImageProcessingTemplate
  \brief   This class is a template of an image processing module to serve as an example
  \details This class is a child of vImageProcessingModule. It is a template of an image
           processing module that implements all mandatory parts to be compilable but does not
           do nothing. It can be used as a starting point to implement ones own image processing
           module. Many explanations are provided within the corpus of all functions to guide
           through the implementation of ones own module. In a few words, one simply has to
           implement few pure virtual functions inherited from the vImageProcessingModule
           abstract class, while observing a few rules.
*/
class iImageProcessingTemplate : public vImageProcessingModule
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public iImageProcessingTemplate::iImageProcessingTemplate()
      \brief   The constructor of iImageProcessingTemplate
      \details This is the default and unique constructor. It does not take any parameter and
               its role is:  \n
               (i) to affect default values to parameters specific to this module; \n
               (ii) to affect the values of the three booleans m_affectTimeDimensionFlag,
               m_affectRespDimensionFlag and m_affectCardDimensionFlag \n
               These booleans are inherited from the abstract class vImageProcessingModule, that specifies
               if this specific module will affect or operate on the different dynamic dimensions.
               This is used by the oImageProcessingManager to check that dynamic basis functions are not used
               at the same time as a module that work on these dimensions too.
    */
    iImageProcessingTemplate();
    /*!
      \fn      public iImageProcessingTemplate::~iImageProcessingTemplate()
      \brief   The destructor of iImageProcessingTemplate
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built specifically
               by this module.
    */
    ~iImageProcessingTemplate();


  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do NOT add semi-colon at the end of the line)
    FUNCTION_IMAGE_PROCESSING_MODULE(iImageProcessingTemplate)
    /*!
      \fn      public int iImageProcessingTemplate::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to a child module, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vImageProcessingModule. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_configurationFile);
    /*!
      \fn      public int iImageProcessingTemplate::ReadOptionsList()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to a child module, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vImageProcessingModule. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_optionsList);
    /*!
      \fn      public int iImageProcessingTemplate::ShowHelp()
      \brief   A function used to show help about the child module
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the module, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vImageProcessingModule.
    */
    void ShowHelp();
    /*!
      \fn      public int iImageProcessingTemplate::Process()
      \param   FLTNB**** a4p_image
      \brief   A function used to actually perform the processing
      \details This function implements processing specific to the child module. The provided parameter
               is the image to be processed. First pointer for the time basis functions/frames dimension,
               second pointer for the respiratory basis functions/gates dimension, third pointer for the
               cardiac basis functions/gates dimension, and fourth pointer for the voxels. It is the
               implementation of the pure virtual function inherited from the abstract class
               vImageProcessingModule.
      \return  An integer reflecting the processing success; 0 if success, another value otherwise.
    */
    int Process(FLTNB**** a4p_image);


  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int iImageProcessingTemplate::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child module
      \details This function is used to check that all parameters specific to the module are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vImageProcessingModule.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iImageProcessingTemplate::InitializeSpecific()
      \brief   A private function used to initialize everything specific to the child module
      \details This function is used to initialize everything specific to the module that should be
               initialized. It is called by the Initialize() function of the mother class. It is
               the implementation of the pure virtual function inherited from the abstract mother
               class vImageProcessingModule.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int InitializeSpecific();


  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    // Declare here any members useful for this specific module
};


// Class for automatic insertion (set here the visible image processing module's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_IMAGE_PROCESSING_MODULE(template,iImageProcessingTemplate)

#endif
