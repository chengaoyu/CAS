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
  \brief    Declaration of class oImageProcessingManager
*/

#ifndef OIMAGEPROCESSINGMANAGER_HH
#define OIMAGEPROCESSINGMANAGER_HH 1

#include "gVariables.hh"
#include "oImageSpace.hh"
#include "vImageProcessingModule.hh"

/*!
  \class   oImageProcessingManager
  \brief   This class is designed to manage the different image processing modules and to apply them
  \details This manager class is supposed to be created and initialized in the main program.
           To do so, the following steps must be used: (i) The empty constructor is called which affect
           all members with default values. (ii) All parameters are set through the use of SetXXX()
           functions. (iii) The CheckParameters() function is called to check that everything mandatory
           has been set. (iv) The Initialize() function is called to initialize everything. (v) Now
           the action functions of the manager can be called to apply the different image processing modules.
           In a few words, based on supplied options, the manager will create children of the
           abstract vImageProcessingModule class which are specific image processing modules. As an example,
           see the iImageProcessingTemplate child class that illustrates how a specific image processing
           module should be implemented.
*/
class oImageProcessingManager
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public oImageProcessingManager::oImageProcessingManager()
      \brief   The constructor of oImageProcessingManager
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    oImageProcessingManager();
    /*!
      \fn      public oImageProcessingManager::~oImageProcessingManager()
      \brief   The destructor of oImageProcessingManager
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were build by this class.
    */
    ~oImageProcessingManager();

  // -----------------------------------------------------------------------------------------
  // Public member functions for initialization
  public:
    /*!
      \fn      public int oImageProcessingManager::CheckParameters()
      \brief   A function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int oImageProcessingManager::Initialize()
      \brief   A function used to initialize the manager and all image processing modules it manages
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. In a few words, it parses the options, then creates and
               initializes all image processing modules based on the provided options.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int Initialize();
    /*!
      \fn      public static void oImageProcessingManager::ShowCommonHelp()
      \brief   This function does not take any parameter and is used to display some help about
               the syntax of the options describing the image processing module that should be used.
               It is static so that it can be called without any object initialization.
    */
    static void ShowCommonHelp();


  // -----------------------------------------------------------------------------------------
  // Public member functions for actions
  public:
    /*!
      \fn      public int oImageProcessingManager::ConvolveForward()
      \param   oImageSpace* ap_ImageSpace
      \brief   A function used to apply image processing modules onto the forward image of the
               oImageSpace
      \details Based on the different mp_applyForward of all managed image processing modules,
               it will apply them or not onto the forward image of the provided oImageSpace.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ApplyProcessingForward(oImageSpace* ap_ImageSpace);
    /*!
      \fn      public int oImageProcessingManager::ConvolveIntra()
      \param   oImageSpace* ap_ImageSpace
      \brief   A function used to apply image processing modules onto the current image of
               the oImageSpace
      \details Based on the different mp_applyIntra of all managed image processing modules, it
               will apply them or not onto the current estimated image of the provided oImageSpace.
               The processed image is put back as the current estimate for the next update.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ApplyProcessingIntra(oImageSpace* ap_ImageSpace);
    /*!
      \fn      public int oImageProcessingManager::ConvolvePost()
      \param   oImageSpace* ap_ImageSpace
      \brief   A function used to apply convolvers onto the output image of the oImageSpace
      \details Based on the different mp_applyPost of all managed image processing modules,
               it will apply them or not onto the output image of the provided oImageSpace.
               This function can be used when one wants to apply processing onto the image
               as a post-processing step, right before being saved. The processed image is
               only used to be saved and is not put back into the iterative process.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ApplyProcessingPost(oImageSpace* ap_ImageSpace);


  // -----------------------------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline void oImageProcessingManager::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the member m_verboseLevel to the provided value
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      public inline void oImageProcessingManager::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the member mp_ImageDimensionsAndQuantification to the provided value
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification) 
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      public inline void oImageProcessingManager::SetOptions()
      \param   vector<string> a_options
      \brief   Set the member m_options to the provided value
    */
    inline void SetOptions(vector<string> a_options)
           {m_options = a_options;}


  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int oImageProcessingManager::ParseOptionsAndInitializeImageProcessingModules()
      \brief   A function used to parse options and initialize image processing modules
      \details This function first parses the options contained in the member m_options. Each
               string of the vector describes an image processing module to be used. Based on a
               specific syntax, the options are parsed to get the name of the processing module,
               its associated parameters and the steps of application. Based on this, the image
               processing modules are initialized. This function is private because it is called
               by the Initialize() function.
      \return  An integer reflecting the parsing and initialization status; 0 if no problem,
               another value otherwise.
    */
    int ParseOptionsAndInitializeImageProcessingModules();


  // -----------------------------------------------------------------------------------------
  // Data members
  private:
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification; /*!< The image dimensions */
    vector<string> m_options;              /*!< A vector containing strings of options, each string is associated to an image processing module */
    int m_nbImageProcessingModules;        /*!< The number of image processing modules managed by this manager */
    vImageProcessingModule** 
      m2p_ImageProcessingModules;          /*!< The actual image processing modules (as many as m_nbImageProcessingModules) */
    bool* mp_applyForward;                 /*!< As many booleans as m_nbImageProcessingModules specifying if each module should be apply within the ConvolveForward function */
    bool* mp_applyIntra;                   /*!< As many booleans as m_nbImageProcessingModules specifying if each module should be apply within the ConvolveIntra function */
    bool* mp_applyPost;                    /*!< As many booleans as m_nbImageProcessingModules specifying if each module should be apply within the ConvolvePost function */
    bool m_checked;                        /*!< A boolean that says if the function CheckParameters() has been called */
    bool m_initialized;                    /*!< A boolean that says if the function Initialize() has been called */
    int m_verbose;                         /*!< The verbose level associated to this class */
};

#endif
