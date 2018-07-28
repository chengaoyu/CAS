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
  \ingroup management
  \brief Declaration of class sAddonManager
*/

#ifndef SADDONMANAGER_HH
#define SADDONMANAGER_HH 1

#include "gVariables.hh"
#include "vProjector.hh"
#include "vOptimizer.hh"
#include "vImageConvolver.hh"
#include "vImageProcessingModule.hh"
#include "vPenalty.hh"
#include "vDeformation.hh"
#include "vDynamicModel.hh"

/*!
  \class   sAddonManager
  \brief   This class is designed to manage the automatic declaration of 'addon' classes
  \details An 'addon' class is any class deriving from an abstract class used for the implementation of specific modules. \n 
           As an example, the class vProjector is an abstract class that cannot be used on its own, but all its children \n 
           are implementing actual projectors that can be used. \n 
           This addon manager allows to automatically declare the existency of an 'addon' class into the CASToR project. \n 
           Based on specific macros to be added in the 'addon' class, and what is defined in this addon manager, \n
           there is nothing more to do to be able to use the 'addon' inside CASToR. \n 
           This class is a singleton. For all things to work properly, additional macros are defined inside \n
           the header file of abstract classes. \n 
           The philosophy of this 'addon' mechanism is taken from the GATE software.
*/
class sAddonManager
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public :
    /*!
      \fn      public static sAddonManager* sAddonManager::GetInstace()
      \brief   Instantiate the singleton if not already done, then return the pointer to its instance
      \return  mp_Instance
    */
    static sAddonManager *GetInstance()
    {
      if (mp_Instance == NULL) mp_Instance = new sAddonManager;
      return mp_Instance;
    };
    /*!
      \fn      public sAddonManager::~sAddonManager()
      \brief   The destructor of sAddonManager
      \details This is the default and unique destructor. It simply replace the instance pointer by NULL.
    */
    ~sAddonManager() {mp_Instance = NULL;}


  // -----------------------------------------------------------------------------------------
  // Public lists and methods for the different objects based on abstract classes
  public:

    // The function's pointer to build the map of all functions dedicated to the creation of vProjector objects
    typedef vProjector *(*maker_projector) ();
    // The map indexed by the name of the vProjector children, and which points on the corresponding function that
    // call the constructor of each child. The corresponding function is created by a macro defined in the vProjector
    // header file and used in the children classes. The pointer to this function is automatically added into the map
    // by another macro defined in the vProjector header file (describring a static creator class) and used in the
    // children classes.
    std::map<string,maker_projector> mp_listOfProjectors; /*!< The map of all vProjector children's constructors */
    /*!
      \fn      public void sAddonManager::ShowHelpProjector()
      \brief   Show help about all implemented projectors
      \details Loop over all entries of mp_listOfProjectors and call the associated ShowHelp() function of vProjector.
    */
    void ShowHelpProjector();

    // The function's pointer to build the map of all functions dedicated to the creation of vOptimizer objects
    typedef vOptimizer *(*maker_optimizer) ();
    // The map indexed by the name of the vOptimizer children, and which points on the corresponding function that
    // call the constructor of each child. The corresponding function is created by a macro defined in the vOptimizer
    // header file and used in the children classes. The pointer to this function is automatically added into the map
    // by another macro defined in the vOptimizer header file (describring a static creator class) and used in the
    // children classes.
    std::map<string,maker_optimizer> mp_listOfOptimizers;
    /*!
      \fn      public void sAddonManager::ShowHelpOptimizer()
      \brief   Show help about all implemented optimizers
      \details Loop over all entries of mp_listOfOptimizers and call the associated ShowHelp() function of vOptimizer.
    */
    void ShowHelpOptimizer();

    // The function's pointer to build the map of all functions dedicated to the creation of vPenalty objects
    typedef vPenalty *(*maker_penalty) ();
    // The map indexed by the name of the vPenalty children, and which points on the corresponding function that
    // call the constructor of each child. The corresponding function is created by a macro defined in the vPenalty
    // header file and used in the children classes. The pointer to this function is automatically added into the map
    // by another macro defined in the vPenalty header file (describring a static creator class) and used in the
    // children classes.
    std::map<string,maker_penalty> mp_listOfPenalties;
    /*!
      \fn      public void sAddonManager::ShowHelpPenalty()
      \brief   Show help about all implemented penalties
      \details Loop over all entries of mp_listOfPenalties and call the associated ShowHelp() function of vPenalty.
    */
    void ShowHelpPenalty();

    // The function's pointer to build the map of all functions dedicated to the creation of vImageConvolver objects
    typedef vImageConvolver *(*maker_image_convolver) ();
    // The map indexed by the name of the vImageConvolver children, and which points on the corresponding function that
    // call the constructor of each child. The corresponding function is created by a macro defined in the vImageConvolver
    // header file and used in the children classes. The pointer to this function is automatically added into the map
    // by another macro defined in the vImageConvolver header file (describring a static creator class) and used in the
    // children classes.
    std::map<string,maker_image_convolver> mp_listOfImageConvolvers;
    /*!
      \fn      public void sAddonManager::ShowHelpImageConvolver()
      \brief   Show help about all implemented image convolvers
      \details Loop over all entries of mp_listOfImageConvolvers and call the associated ShowHelp() function of vImageConvolver.
    */
    void ShowHelpImageConvolver();

    // The function's pointer to build the map of all functions dedicated to the creation of vImageProcessingModule objects
    typedef vImageProcessingModule *(*maker_image_processing_module) ();
    // The map indexed by the name of the vImageProcessingModule children, and which points on the corresponding function that
    // call the constructor of each child. The corresponding function is created by a macro defined in the vImageProcessingModule
    // header file and used in the children classes. The pointer to this function is automatically added into the map
    // by another macro defined in the vImageProcessingModule header file (describring a static creator class) and used in the
    // children classes.
    std::map<string,maker_image_processing_module> mp_listOfImageProcessingModules;
    /*!
      \fn      public void sAddonManager::ShowHelpImageProcessingModule()
      \brief   Show help about all implemented image processing modules
      \details Loop over all entries of mp_listOfImageProcessingModules and call the associated ShowHelp() function of vImageProcessingModule.
    */
    void ShowHelpImageProcessingModule();
    
    // The function's pointer to build the map of all functions dedicated to the creation of vScanner objects
    typedef vScanner *(*maker_scanner) ();
    // The map indexed by the name of the vScanner children, and which points on the corresponding function that
    // call the constructor of each child. The corresponding function is created by a macro defined in the vScanner
    // header file and used in the children classes. The pointer to this function is automatically added into the map
    // by another macro defined in the vScanner header file (describring a static creator class) and used in the
    // children classes.
    std::map<string,maker_scanner> mp_listOfScannerTypes;
    /*!
      \fn      public void sAddonManager::ShowHelpScanner()
      \brief   Show help about all implemented scanners
      \details Loop over all entries of mp_listOfScannerTypes and call the associated ShowHelp() function of vScanner.
    */
    void ShowHelpScanner();

    // The function's pointer to build the map of all functions dedicated to the creation of vDynamicModel objects
    typedef vDynamicModel *(*maker_dynamic_model) ();
    // The map indexed by the name of the vDynamicModel children, and which points on the corresponding function that
    // call the constructor of each child. The corresponding function is created by a macro defined in the vDynamicModel
    // header file and used in the children classes. The pointer to this function is automatically added into the map
    // by another macro defined in the vDynamicModel header file (describring a static creator class) and used in the
    // children classes.
    std::map<string,maker_dynamic_model> mp_listOfDynamicModels;
    /*!
      \fn      public void sAddonManager::ShowHelpDynamicModel()
      \brief   Show help about all implemented dynamic models
      \details Loop over all entries of mp_listOfDynamicModels and call the associated ShowHelp() function of vDynamicModel.
    */
    void ShowHelpDynamicModel();

    // The function's pointer to build the map of all functions dedicated to the creation of vDeformation objects
    typedef vDeformation *(*maker_deformation) ();
    // The map indexed by the name of the vDeformation children, and which points on the corresponding function that
    // call the constructor of each child. The corresponding function is created by a macro defined in the vDeformation
    // header file and used in the children classes. The pointer to this function is automatically added into the map
    // by another macro defined in the vDeformation header file (describring a static creator class) and used in the
    // children classes.
    std::map<string,maker_deformation> mp_listOfDeformations;
    /*!
      \fn      public void sAddonManager::ShowHelpDeformation()
      \brief   Show help about all implemented deformations
      \details Loop over all entries of mp_listOfDeformations and call the associated ShowHelp() function of vDeformation.
    */
    void ShowHelpDeformation();


  // -------------------------------------------------------------------
  // Private constructor
  private:
    /*!
      \fn      private sAddonManager::sAddonManager()
      \brief   The constructor of sAddonManager
      \details This is the default and unique constructor. It does nothing.
    */
    sAddonManager();


  // -------------------------------------------------------------------
  // Private data members
  private:
    static sAddonManager *mp_Instance; /*!< Pointer to this singleton object */   
};

#endif
