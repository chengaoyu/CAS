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
  \brief    Declaration of class vImageConvolver
*/

#ifndef VIMAGECONVOLVER_HH
#define VIMAGECONVOLVER_HH 1

#include "gVariables.hh"
#include "gOptions.hh"

class oImageSpace;
class oImageDimensionsAndQuantification;

/*!
  \class   vImageConvolver
  \brief   This abstract class is the generic image convolver class used by the oImageConvolverManager
  \details This abstract class is the base of all implemented image convolvers inheriting from it. It
           is used by the oImageConvolverManager that instantiate a collection of children objects based on the
           provided options. It implements four main public functions: \n
           (i) CheckParameters() which checks the mandatory common parameters and calls the pure virtual 
           CheckSpecificParameters() function implemented by each child; \n
           (ii) Initialize() which initializes some common stuff and calls the pure virtual BuildConvolutionKernel()
           function implemented by each child; \n
           (iii) ApplyConvolution() which actually applies the convolution onto the provided image; \n
           (iv) ApplyConvolutionTranspose() which applies the transpose of the convolution. \n
           It also specifies other pure virtual functions dedicated to the reading of options and help associated to
           each child, and the Convolve() and ConvolveTranspose() functions which actually implement the specific
           convolving of each child module. As an example of a child module, see the iImageProcessingTemplate child
           class that illustrates how a specific image processing module should
           be implemented.
*/
class vImageConvolver
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public vImageConvolver::vImageConvolver()
      \brief   The constructor of vImageConvolver
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    vImageConvolver();
    /*!
      \fn      virtual public vImageConvolver::~vImageConvolver()
      \brief   The destructor of vImageConvolver
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were build by this class.
               It is virtual, so that it is automatically called when a child object is deleted.
    */
    virtual ~vImageConvolver();


  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int vImageConvolver::CheckParameters()
      \brief   A public function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized. At the end, it calls the pure virtual
               CheckSpecificParameters() function implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int vImageConvolver::Initialize()
      \brief   A public function used to initialize the module
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. At the beginning, it calls the pure virtual BuildConvolutionKernel()
               function implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int Initialize();
    /*!
      \fn      public int vImageConvolver::ApplyConvolution()
      \param   FLTNB* ap_image
      \brief   A public function used to apply the convolution module on the provided image
      \details This function is the first main action function used to apply the convolution on the provided
               image. It copy the provided image into the padded image buffer using the CopyToPaddedImage()
               private function, and then apply the convolution using the private Convolve() function.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ApplyConvolution(FLTNB* ap_image);
    /*!
      \fn      public int vImageConvolver::ApplyConvolutionTranspose()
      \param   FLTNB* ap_image
      \brief   A public function used to apply the transpose convolution module on the provided image
      \details This function is the second main action function used to apply the transpose of the convolution
               on the provided image. It copy the provided image into the padded image buffer using the
               CopyToPaddedImage() private function, and then apply the transpose of the convolution using
               the private ConvolveTranspose() function.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int ApplyConvolutionTranspose(FLTNB* ap_image);


  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private void vImageConvolver::CopyToPaddedImage()
      \param   FLTNB* ap_inputImage
      \brief   A private function used to copy the provided image into the padded buffer
      \details This function is copies the provided image into the padded image buffer. It is used right
               before apply the convolution or its transpose.
    */
    void CopyToPaddedImage(FLTNB* ap_inputImage);
    /*!
      \fn      private virtual int vImageConvolver::Convolve()
      \param   FLTNB* ap_outputImage
      \brief   A private function used to apply the convolution on the padded image to the provided output image
      \details This function is used to apply the convolution on the padded image to the provided output image.
               It is virtual so it can be overloaded. The current implementation is valid only for stationary
               kernels. So for spatially variant kernels, it must be overloaded.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    virtual int Convolve(FLTNB* ap_outputImage);
    /*!
      \fn      private virtual int vImageConvolver::ConvolveTranspose()
      \param   FLTNB* ap_outputImage
      \brief   A private function used to apply the transpose convolution on the padded image to the provided output image
      \details This function is used to apply the transpose of the convolution on the padded image to the provided
               output image. It is virtual so it can be overloaded. The current implementation is valid only for
               stationary kernels as it simply call the Convolve() function. So for spatially variant kernels, it
               must be overloaded.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    virtual int ConvolveTranspose(FLTNB* ap_outputImage);


  // -----------------------------------------------------------------------------------------
  // Pure virtual public member functions that need to be implemented by children
  public:
    /*!
      \fn      public virtual int vImageConvolver::ReadConfigurationFile() = 0
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to a child convolver, from
               a configuration file. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadConfigurationFile(const string& a_fileOptions) = 0;
    /*!
      \fn      public virtual int vImageConvolver::ReadOptionsList() = 0
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to a child convolver, from
               a list of options. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadOptionsList(const string& a_listOptions) = 0;
    /*!
      \fn      public virtual int vImageConvolver::ShowHelp() = 0
      \brief   A function used to show help about the child module
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the module, and how to set them through the use
               of a configuration file or a list of options. It is pure virtual so is implemented by
               children.
    */
    virtual void ShowHelp() = 0;


  // -----------------------------------------------------------------------------------------
  // Pure virtual private member functions that need to be implemented by children
  private:
    /*!
      \fn      private virtual int vImageConvolver::CheckSpecificParameters() = 0
      \brief   A private function used to check the parameters settings specific to the child convolver
      \details This function is used to check that all parameters specific to the convolver are correctly set
               within allowed values. It is called by the CheckParameters() function. It is pure virtual so
               is implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    virtual int CheckSpecificParameters() = 0;
    /*!
      \fn      private virtual int vImageConvolver::BuildConvolutionKernel() = 0
      \brief   A private function used to build the convolution kernel specific to the child convolver
      \details This function is used to build the convolution kernels associated to the child convolver.
               It is called by the Initialize() function. It is pure virtual so is implemented by
               children. To be the most generic possible, one can build has many convolution kernels
               as desired in order to implement spatially variant convolutions. The number of kernels
               should be specified, the kernels' dimensions allocated and specified, and same for the
               actual kernels' values.
      \return  An integer reflecting the building status; 0 if no problem, another value otherwise.
    */
    virtual int BuildConvolutionKernel() = 0;


  // -----------------------------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public void vImageConvolver::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the member m_verboseLevel to the provided value
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      public void vImageConvolver::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the member mp_ImageDimensionsAndQuantification to the provided value
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}

  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    // Standards
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification; /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    int m_verbose;                         /*!< The verbose level */
    // Booleans
    bool m_checked;        /*!< Boolean that says if the parameters were checked or not */
    bool m_initialized;    /*!< Boolean that says if the convolver was initialized or not */
    bool m_stationary;     /*!< Boolean that says if the kernel is stationary or not */
    // Padded image
    FLTNB* mp_paddedImage; /*!< The actual padded buffer image */
    INTNB m_offsetX;       /*!< The offset of the padded image along X */
    INTNB m_offsetY;       /*!< The offset of the padded image along Y */
    INTNB m_offsetZ;       /*!< The offset of the padded image along Z */
    INTNB m_dimPadX;       /*!< The number of voxels of the padded image along X */
    INTNB m_dimPadY;       /*!< The number of voxels of the padded image along Y */
    INTNB m_dimPadZ;       /*!< The number of voxels of the padded image along Z */
    INTNB m_dimPadXY;      /*!< The number of voxels of the padded image in a slice */
    INTNB m_dimPadXYZ;     /*!< The total number of voxels of the padded image */
    // Convolution kernel
    INTNB m_nbKernels;     /*!< The number of kernels (1 if stationary, more otherwise */
    INTNB* mp_dimKernelX;  /*!< The dimension of each kernel along X */
    INTNB* mp_dimKernelY;  /*!< The dimension of each kernel along Y */
    INTNB* mp_dimKernelZ;  /*!< The dimension of each kernel along Z */
    FLTNB** m2p_kernel;    /*!< The actual kernels, first pointer for the number of kernels, second pointer for the kernel values */
};


// ----------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ----------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_IMAGE_CONVOLVER(CLASS) \
  static vImageConvolver *make_image_convolver() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_IMAGE_CONVOLVER(NAME,CLASS)                                                                \
  class NAME##ImageConvolverCreator                                                                      \
  {                                                                                                      \
    public:                                                                                              \
      NAME##ImageConvolverCreator()                                                                      \
        { sAddonManager::GetInstance()->mp_listOfImageConvolvers[#NAME] = CLASS::make_image_convolver; } \
  };                                                                                                     \
  static NAME##ImageConvolverCreator ImageConvolverCreator##NAME;

#endif
