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
  \brief    Declaration of class iImageConvolverStationaryGaussian
*/

#ifndef IIMAGECONVOLVERSTATIONARYGAUSSIAN_HH
#define IIMAGECONVOLVERSTATIONARYGAUSSIAN_HH 1

#include "gVariables.hh"
#include "gOptions.hh"
#include "oImageSpace.hh"
#include "vImageConvolver.hh"
#include "sAddonManager.hh"

/*!
  \class   iImageConvolverStationaryGaussian
  \brief   This class is an image convolver module implementing stationary gaussian filtering
  \details This class is a child of vImageConvolver. It implements a stationary gaussian
           convolution filter, parameterized by a transaxial and an axial FWHM, as well as the number
           of sigmas to be included in the convolution kernel.
*/
class iImageConvolverStationaryGaussian : public vImageConvolver
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public iImageConvolverStationaryGaussian::iImageConvolverStationaryGaussian()
      \brief   The constructor of iImageConvolverStationaryGaussian
      \details This is the default and unique constructor. It does not take any parameter and
               its role is to affect default values to parameters specific to this module.
    */
    iImageConvolverStationaryGaussian();
    /*!
      \fn      public iImageConvolverStationaryGaussian::~iImageConvolverStationaryGaussian()
      \brief   The destructor of iImageConvolverStationaryGaussian
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built specifically
               by this module. Here it does nothing.
    */
    ~iImageConvolverStationaryGaussian();


  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-colon at the end of the line)
    FUNCTION_IMAGE_CONVOLVER(iImageConvolverStationaryGaussian)
    /*!
      \fn      public int iImageConvolverStationaryGaussian::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to a child module, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vImageConvolver. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_fileOptions);
    /*!
      \fn      public int iImageConvolverStationaryGaussian::ReadOptionsList()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to a child module, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vImageConvolver. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_listOptions);
    /*!
      \fn      public int iImageConvolverStationaryGaussian::ShowHelp()
      \brief   A function used to show help about the child module
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the module, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vImageConvolverModule.
    */
    void ShowHelp();


  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int iImageConvolverStationaryGaussian::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child module
      \details This function is used to check that all parameters specific to the module are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vImageConvolver.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iImageConvolverStationaryGaussian::BuildConvolutionKernel()
      \brief   A private function used to build the convolution kernel specific to the child convolver
      \details This function is used to build the convolution kernels associated to the child convolver.
               It is called by the Initialize() function. It is the implementation of the pure virtual
               function inherited from the abstract mother class vImageConvolver. In this case, the
               kernel is stationary so m_nbKernels = 1. We keep the implementation of the convolution
               from the mother class which already implements stationary convolutions.
      \return  An integer reflecting the building status; 0 if no problem, another value otherwise.
    */
    int BuildConvolutionKernel();


  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    FLTNB m_transFWHM;    /*!< The transaxial FWHM in mm */
    FLTNB m_axialFWHM;    /*!< The axial FWHM in mm */
    FLTNB m_nbSigmas;     /*!< The number of sigmas of the Gaussian distribution included in the kernel */
    INTNB m_dimKernelXY;  /*!< The number of voxels in a slice of the kernel */
    INTNB m_dimKernelXYZ; /*!< The total number of voxels of the kernel */
};


// Class for automatic insertion (set here the visible image convolver's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_IMAGE_CONVOLVER(gaussian,iImageConvolverStationaryGaussian)

#endif
