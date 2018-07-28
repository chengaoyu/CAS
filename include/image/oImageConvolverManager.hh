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
  \brief    Declaration of class oImageConvolverManager
*/

#ifndef OIMAGECONVOLVERMANAGER_HH
#define OIMAGECONVOLVERMANAGER_HH 1

#include "gVariables.hh"
#include "oImageSpace.hh"
#include "vImageConvolver.hh"

/*!
  \class   oImageConvolverManager
  \brief   This class is designed to manage the different image convolvers and to apply them
  \details This manager class is supposed to be created and initialized in the main program.
           To do so, the following steps must be used: \n
           (i) The empty constructor is called which affect all members with default values. \n
           (ii) All parameters are set through the use of SetXXX() functions. \n
           (iii) The CheckParameters() function is called to check that everything mandatory
           has been set. \n 
           (iv) The Initialize() function is called to initialize everything. \n
           (v) Now the action functions of the manager can be called to apply the different image convolvers. \n
           In a few words, based on supplied options, the manager will create children of the
           abstract vImageConvolver class which are specific image convolver modules. As an example,
           see the iImageConvolverTemplate child class that illustrates how a specific image convolver
           module should be implemented.
*/
class oImageConvolverManager
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public oImageConvolverManager::oImageConvolverManager()
      \brief   The constructor of oImageConvolverManager
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    oImageConvolverManager();
    /*!
      \fn      public oImageConvolverManager::~oImageConvolverManager()
      \brief   The destructor of oImageConvolverManager
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were build by this class.
    */
    ~oImageConvolverManager();


  // -----------------------------------------------------------------------------------------
  // Public member functions for initialization
  public:
    /*!
      \fn      public int oImageConvolverManager::CheckParameters()
      \brief   A function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int oImageConvolverManager::Initialize()
      \brief   A function used to initialize the manager and all image convolvers it manages
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. In a few words, it parses the options, then creates and
               initializes all image convolvers based on the provided options.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int Initialize();
    /*!
      \fn      public static void oImageConvolverManager::ShowCommonHelp()
      \brief   This function does not take any parameter and is used to display some help about
               the syntax of the options describing the image convolvers that should be used.
               It is static so that it can be called without any object initialization.
    */
    static void ShowCommonHelp();


  // -----------------------------------------------------------------------------------------
  // Public member functions for actions
  public:
    /*!
      \fn      public int oImageConvolverManager::ConvolveForward()
      \param   oImageSpace* ap_ImageSpace
      \brief   A function used to apply convolvers onto the forward image of the oImageSpace
      \details Based on the different mp_applyForward of all managed convolvers, it will apply
               them or not onto the forward image of the provided oImageSpace. The convolvers
               are applied on all dynamic dimensions. This function is typically used when
               image-based PSF modeling is part of the iterative reconstruction process, so
               that the convolution is applied onto the image to be forward projected.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ConvolveForward(oImageSpace* ap_ImageSpace);
    /*!
      \fn      public int oImageConvolverManager::ConvolveBackward()
      \param   oImageSpace* ap_ImageSpace
      \brief   A function used to apply convolvers onto the backward images of the oImageSpace
      \details Based on the different mp_applyBackward of all managed convolvers, it will apply
               them or not onto the backward images of the provided oImageSpace. The convolvers
               are applied on all dynamic dimensions. This function is typically used when
               image-based PSF modeling is part of the iterative reconstruction process, so
               that the convolution is applied onto the correction images after the backward
               projection. Note that in the case of histogram-based reconstruction, the
               convolution will also be applied to the sensitivity image that is computed
               in synchronization with the backward images containing the correction terms.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ConvolveBackward(oImageSpace* ap_ImageSpace);
    /*!
      \fn      public int oImageConvolverManager::ConvolveSensitivity()
      \param   oImageSpace* ap_ImageSpace
      \brief   A function used to apply convolvers onto the sensitivity image of the oImageSpace
      \details Based on the different mp_applyBackward of all managed convolvers, it will apply
               them or not onto the sensitivity image of the provided oImageSpace. The convolvers
               are applied on all dynamic dimensions. This function is typically used when
               image-based PSF modeling is part of the iterative listmode-based reconstruction
               process so that the convolution is applied onto the first-computed sensitivity image.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ConvolveSensitivity(oImageSpace* ap_ImageSpace);
    /*!
      \fn      public int oImageConvolverManager::ConvolveIntra()
      \param   oImageSpace* ap_ImageSpace
      \brief   A function used to apply convolvers onto the current image of the oImageSpace
      \details Based on the different mp_applyIntra of all managed convolvers, it will apply
               them or not onto the current estimated image of the provided oImageSpace. The
               convolvers are applied on all dynamic dimensions. This function can be used when
               one wants to apply convolution onto the current estimated image right after it has
               been updated within the iterative process. The convolved image is thus put back as
               the current estimate for the next update.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ConvolveIntra(oImageSpace* ap_ImageSpace);
    /*!
      \fn      public int oImageConvolverManager::ConvolvePost()
      \param   oImageSpace* ap_ImageSpace
      \brief   A function used to apply convolvers onto the output image of the oImageSpace
      \details Based on the different mp_applyPost of all managed convolvers, it will apply
               them or not onto the output image of the provided oImageSpace. The convolvers are
               applied on all dynamic dimensions. This function can be used when one wants to
               apply convolution onto the image as a post-processing step, right before being
               saved. The convolved image is only used to be saved and is not put back into the
               iterative process.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ConvolvePost(oImageSpace* ap_ImageSpace);


  // -----------------------------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline void oImageConvolverManager::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the member m_verboseLevel to the provided value
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      public inline void oImageConvolverManager::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the member mp_ImageDimensionsAndQuantification to the provided value
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification) 
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      public inline void oImageConvolverManager::SetOptions()
      \param   vector<string> a_options
      \brief   Set the member m_options to the provided value
    */
    inline void SetOptions(vector<string> a_options)
           {m_options = a_options;}


  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int oImageConvolverManager::ParseOptionsAndInitializeImageConvolvers()
      \brief   A function used to parse options and initialize image convolvers
      \details This function first parses the options contained in the member m_options. Each
               string of the vector describes an image convolver to be used. Based on a specific
               syntax, the options are parsed to get the name of the convolver module, its
               associated parameters and the steps of application. Based on this, the image
               convolvers are initialized. This function is private because it is called by the
               Initialize() function.
      \return  An integer reflecting the parsing and initialization status; 0 if no problem,
               another value otherwise.
    */
    int ParseOptionsAndInitializeImageConvolvers();


  // -----------------------------------------------------------------------------------------
  // Data members
  private:
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification; /*!< The image dimensions */
    vector<string> m_options;              /*!< A vector containing strings of options, each string is associated to a convolver */
    int m_nbImageConvolvers;               /*!< The number of convolvers managed by this manager */
    vImageConvolver** m2p_ImageConvolvers; /*!< The actual convolvers (as many as m_nbImageConvolvers) */
    bool* mp_applyForward;                 /*!< As many booleans as m_nbImageConvolvers specifying if each convolver should be apply within the ConvolveForward function */
    bool* mp_applyBackward;                /*!< As many booleans as m_nbImageConvolvers specifying if each convolver should be apply within the ConvolveBackward function */
    bool* mp_applyIntra;                   /*!< As many booleans as m_nbImageConvolvers specifying if each convolver should be apply within the ConvolveIntra function */
    bool* mp_applyPost;                    /*!< As many booleans as m_nbImageConvolvers specifying if each convolver should be apply within the ConvolvePost function */
    bool m_checked;                        /*!< A boolean that says if the function CheckParameters() has been called */
    bool m_initialized;                    /*!< A boolean that says if the function Initialize() has been called */
    int m_verbose;                         /*!< The verbose level associated to this class */
};

#endif
