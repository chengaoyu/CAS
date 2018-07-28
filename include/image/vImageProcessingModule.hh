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
  \brief    Declaration of class vImageProcessingModule
*/

#ifndef VIMAGEPROCESSINGMODULE_HH
#define VIMAGEPROCESSINGMODULE_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oImageSpace.hh"

/*!
  \class   vImageProcessingModule
  \brief   This abstract class is the generic image processing module class used by the oImageProcessingManager
  \details This abstract class is the base of all implemented image processing modules inheriting from it. \n
           It is used by the oImageProcessingManager that instantiate a collection of children objects based on the
           provided options. It implements two public functions:  \n
           (i) CheckParameters() which checks the mandatory common parameters and calls the pure virtual 
           CheckSpecificParameters() function implemented by each child; \n
           (ii) Initialize() which initializes some common stuff and calls the pure virtual InitializeSpecific()
           function implemented by each child. \n
           It also specifies other pure virtual functions dedicated to the reading of options and help associated 
           to each child, and the main Process() function which actually implements the specific processing of 
           each child module. As an example of a child module, see the iImageProcessingTemplate child class 
           that illustrates how a specific image processing module should be implemented.
*/
class vImageProcessingModule
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public vImageProcessingModule::vImageProcessingModule()
      \brief   The constructor of vImageProcessingModule
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    vImageProcessingModule();
    /*!
      \fn      virtual public vImageProcessingModule::~vImageProcessingModule()
      \brief   The destructor of vImageProcessingModule
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were build by this class. \n
               It is virtual, so that it is automatically called when a child object is deleted.
    */
    virtual ~vImageProcessingModule();

  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int vImageProcessingModule::CheckParameters()
      \brief   A public function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized. At the end, it calls the pure virtual
               CheckSpecificParameters() function implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int vImageProcessingModule::Initialize()
      \brief   A public function used to initialize the module
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. At the end, it calls the pure virtual InitializeSpecific()
               function implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int Initialize();


  // -----------------------------------------------------------------------------------------
  // Pure virtual public member functions that need to be implemented by children
  public:
    /*!
      \fn      public virtual int vImageProcessingModule::ReadConfigurationFile() = 0
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to a child module, from
               a configuration file. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadConfigurationFile(const string& a_configurationFile) = 0;
    /*!
      \fn      public virtual int vImageProcessingModule::ReadOptionsList() = 0
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to a child module, from
               a list of options. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadOptionsList(const string& a_optionsList) = 0;
    /*!
      \fn      public virtual int vImageProcessingModule::ShowHelp() = 0
      \brief   A function used to show help about the child module
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the module, and how to set them through the use
               of a configuration file or a list of options. It is pure virtual so is implemented by
               children.
    */
    virtual void ShowHelp() = 0;
    /*!
      \fn      public virtual int vImageProcessingModule::Process() = 0
      \param   FLTNB**** a4p_image
      \brief   A function used to actually perform the processing
      \details This function implements processing specific to the child module. The provided parameter
               is the image to be processed. First pointer for the time basis functions/frames dimension,
               second pointer for the respiratory basis functions/gates dimension, third pointer for the
               cardiac basis functions/gates dimension, and fourth pointer for the voxels. It is pure
               virtual so is implented by children.
      \return  An integer reflecting the processing success; 0 if success, another value otherwise.
    */
    virtual int Process(FLTNB**** a4p_image) = 0;


  // -----------------------------------------------------------------------------------------
  // Pure virtual private member functions that need to be implemented by children
  private:
    /*!
      \fn      private virtual int vImageProcessingModule::CheckSpecificParameters() = 0
      \brief   A private function used to check the parameters settings specific to the child module
      \details This function is used to check that all parameters specific to the module are correctly set
               within allowed values. It is called by the CheckParameters() function. It is pure virtual
               so is implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    virtual int CheckSpecificParameters() = 0;
    /*!
      \fn      private virtual int vImageProcessingModule::InitializeSpecific() = 0
      \brief   A private function used to initialize everything specific to the child module
      \details This function is used to initialize everything specific to the module that should be
               initialized. It is called by the Initialize() function. It is pure virtual so is
               implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    virtual int InitializeSpecific() = 0;


  // -----------------------------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public void vImageProcessingModule::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the member m_verboseLevel to the provided value
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      public void vImageProcessingModule::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the member mp_ImageDimensionsAndQuantification to the provided value
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      public bool vImageProcessingModule::GetAffectTimeDimensionFlag()
      \brief   Return the boolean value of m_affectTimeDimensionFlag member
      \return  m_affectTimeDimensionFlag
    */
    inline bool GetAffectTimeDimensionFlag()
           {return m_affectTimeDimensionFlag;}
    /*!
      \fn      public bool vImageProcessingModule::GetAffectRespDimensionFlag()
      \brief   Return the boolean value of m_affectRespDimensionFlag member
      \return  m_affectRespDimensionFlag
    */
    inline bool GetAffectRespDimensionFlag()
           {return m_affectRespDimensionFlag;}
    /*!
      \fn      public bool vImageProcessingModule::GetAffectCardDimensionFlag()
      \brief   Return the boolean value of m_affectCardDimensionFlag member
      \return  m_affectCardDimensionFlag
    */
    inline bool GetAffectCardDimensionFlag()
           {return m_affectCardDimensionFlag;}


  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification;  /*!< The image dimensions */
    bool m_affectTimeDimensionFlag;         /*!< A boolean that specify if the module is affecting the time frame dynamic dimension */
    bool m_affectRespDimensionFlag;         /*!< A boolean that specify if the module is affecting the respiratory dynamic dimension */
    bool m_affectCardDimensionFlag;         /*!< A boolean that specify if the module is affecting the cardiac dynamic dimension */
    bool m_checked;                         /*!< A boolean that says if the function CheckParameters() has been called */
    bool m_initialized;                     /*!< A boolean that says if the function Initialize() has been called */
    int m_verbose;                          /*!< The verbose level associated to this class */
};


// ----------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ----------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_IMAGE_PROCESSING_MODULE(CLASS) \
  static vImageProcessingModule *make_image_processing_module() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_IMAGE_PROCESSING_MODULE(NAME,CLASS)                                                                       \
  class NAME##ImageProcessingModuleCreator                                                                              \
  {                                                                                                                     \
    public:                                                                                                             \
      NAME##ImageProcessingModuleCreator()                                                                              \
        { sAddonManager::GetInstance()->mp_listOfImageProcessingModules[#NAME] = CLASS::make_image_processing_module; } \
  };                                                                                                                    \
  static NAME##ImageProcessingModuleCreator ImageProcessingModuleCreator##NAME;

#endif
