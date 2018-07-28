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
  \brief    Implementation of class iImageConvolverStationaryGaussian
*/

#include "vImageConvolver.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oImageSpace.hh"
#include "iImageConvolverStationaryGaussian.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iImageConvolverStationaryGaussian::iImageConvolverStationaryGaussian() : vImageConvolver()
{
  m_nbSigmas = -1.;
  m_transFWHM = -1.;
  m_axialFWHM = -1.;
  m_dimKernelXY = -1;
  m_dimKernelXYZ = -1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iImageConvolverStationaryGaussian::~iImageConvolverStationaryGaussian()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iImageConvolverStationaryGaussian::ShowHelp()
{
  cout << "This convolver is based on a classic stationary Gaussian kernel." << endl;
  cout << "It can be anisotropic so that the transaxial and axial FWHM can be different." << endl;
  cout << "One also have to choose the number of sigmas of the Gaussian distribution that will be included in the kernel." << endl;
  cout << "The following options can be used (in this particular order when provided as a list, separated by commas):" << endl;
  cout << "  transaxial FWHM: to set the transaxial FWHM (in mm)" << endl;
  cout << "  axial FWHM: to set the axial FWHM (in mm)" << endl;
  cout << "  number of sigmas: to set the number of sigmas included in the kernel (recommendations: at least 3. and maximum 5.)" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryGaussian::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the transaxial FWHM option
  key_word = "transaxial FWHM";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_transFWHM, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iImageConvolverStationaryGaussian::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the axial FWHM option
  key_word = "axial FWHM";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_axialFWHM, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iImageConvolverStationaryGaussian::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the number of sigma option
  key_word = "number of sigmas";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_nbSigmas, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iImageConvolverStationaryGaussian::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryGaussian::ReadOptionsList(const string& a_optionsList)
{
  // There are 3 floating point variables as options
  FLTNB options[3];
  // Read them
  if (ReadStringOption(a_optionsList, options, 3, ",", "Stationary Gaussian convolver configuration"))
  {
    Cerr("***** iImageConvolverStationaryGaussian::ReadConfigurationFile() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_transFWHM = options[0];
  m_axialFWHM = options[1];
  m_nbSigmas = options[2];
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryGaussian::CheckSpecificParameters()
{
  // Check number of sigmas
  if (m_nbSigmas<=0.)
  {
    Cerr("***** iImageConvolverStationaryGaussian::CheckSpecificParameters() -> Number of sigmas included in the kernel must be strictly positive !" << endl);
    return 1;
  }
  // Check the transaxial FWHM
  if (m_transFWHM<0.)
  {
    Cerr("***** iImageConvolverStationaryGaussian::CheckSpecificParameters() -> Transaxial FWHM is negative !" << endl);
    return 1;
  }
  // Check the axial FWHM
  if (m_axialFWHM<0.)
  {
    Cerr("***** iImageConvolverStationaryGaussian::CheckSpecificParameters() -> Axial FWHM is negative !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryGaussian::BuildConvolutionKernel()
{
  // Verbose
  if (m_verbose>=2)
  {
    Cout("iImageConvolverStationaryGaussian::BuildConvolutionKernel() -> Compute convolution kernel" << endl);
    Cout("  --> Transaxial FWHM: " << m_transFWHM << " mm" << endl);
    Cout("  --> Axial FWHM: " << m_axialFWHM << " mm" << endl);
    Cout("  --> Included sigmas: " << m_nbSigmas << endl);
  }

  // --------------------------------
  // Stationary kernel
  // --------------------------------

  // Set the number of kernels to 1
  m_nbKernels = 1;
  // Specify that it is stationary (we do not overload the Convolve() function here)
  m_stationary = true;
  // Compute kernel half size based on the number of sigmas included in the kernel
  FLTNB two_sqrt_two_ln_two = 2.*sqrt(2.*log(2.));
  FLTNB half_kern_floatX = m_transFWHM * m_nbSigmas / (two_sqrt_two_ln_two*mp_ImageDimensionsAndQuantification->GetVoxSizeX());
  FLTNB half_kern_floatY = m_transFWHM * m_nbSigmas / (two_sqrt_two_ln_two*mp_ImageDimensionsAndQuantification->GetVoxSizeY());
  FLTNB half_kern_floatZ = m_axialFWHM * m_nbSigmas / (two_sqrt_two_ln_two*mp_ImageDimensionsAndQuantification->GetVoxSizeZ());
  // For each half size, if the decimal part is below 0.5, then it will be taken into account into the central kernel's voxel,
  // otherwise, we had a voxel to the half size
  INTNB half_kern_dimX = ((INTNB)(half_kern_floatX));
  if (half_kern_floatX-((FLTNB)half_kern_dimX)>0.5) half_kern_dimX++;
  INTNB half_kern_dimY = ((INTNB)(half_kern_floatY));
  if (half_kern_floatY-((FLTNB)half_kern_dimY)>0.5) half_kern_dimY++;
  INTNB half_kern_dimZ = ((INTNB)(half_kern_floatZ));
  if (half_kern_floatZ-((FLTNB)half_kern_dimZ)>0.5) half_kern_dimZ++;
  // Allocate kernel dimensions
  mp_dimKernelX = (INTNB*)malloc(1*sizeof(INTNB));
  mp_dimKernelY = (INTNB*)malloc(1*sizeof(INTNB));
  mp_dimKernelZ = (INTNB*)malloc(1*sizeof(INTNB));
  // Compute kernel dimensions
  mp_dimKernelX[0] = half_kern_dimX * 2 + 1;
  mp_dimKernelY[0] = half_kern_dimY * 2 + 1;
  mp_dimKernelZ[0] = half_kern_dimZ * 2 + 1;
  m_dimKernelXY = mp_dimKernelX[0] * mp_dimKernelY[0];
  m_dimKernelXYZ = m_dimKernelXY * mp_dimKernelZ[0];
  // Verbose
  if (m_verbose>=2) Cout ("  --> Kernel dimensions: [" << mp_dimKernelX[0] << ";" << mp_dimKernelY[0] << ";" << mp_dimKernelZ[0] << "]" << endl);
  // Allocate kernel
  m2p_kernel = (FLTNB**)malloc(1*sizeof(FLTNB*));
  m2p_kernel[0] = (FLTNB*)malloc(m_dimKernelXYZ*sizeof(FLTNB));

  // --------------------------------
  // Compute the kernel
  // --------------------------------

  // Variables
  FLTNB sigma_x = m_transFWHM / two_sqrt_two_ln_two;
  FLTNB sigma_y = m_transFWHM / two_sqrt_two_ln_two;
  FLTNB sigma_z = m_axialFWHM / two_sqrt_two_ln_two;
  FLTNB mu_x = (FLTNB)(half_kern_dimX);
  FLTNB mu_y = (FLTNB)(half_kern_dimY);
  FLTNB mu_z = (FLTNB)(half_kern_dimZ);
  // Compute kernel by integration of the Gaussian distribution (this method is not exactly exact... but is closer to it compared to single Gaussian points)
  FLTNB sum_kernel = 0.;
  for (int z=0, index=0; z<mp_dimKernelZ[0]; z++)
  {
    FLTNB kern_Z = 0.5 * ( erf(  ( (((FLTNB)z)-mu_z+0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeZ()) / (sqrt(2.)*sigma_z)  )
                         - erf(  ( (((FLTNB)z)-mu_z-0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeZ()) / (sqrt(2.)*sigma_z)  ) );
    for (int y=0; y<mp_dimKernelY[0]; y++)
    {
      FLTNB kern_ZY = 0.5 * ( erf(  ( (((FLTNB)y)-mu_y+0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeY()) / (sqrt(2.)*sigma_y)  )
                            - erf(  ( (((FLTNB)y)-mu_y-0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeY()) / (sqrt(2.)*sigma_y)  ) )
                          * kern_Z;
      for (int x=0; x<mp_dimKernelX[0]; x++, index++)
      {
        m2p_kernel[0][index] = 0.5 * ( erf(  ( (((FLTNB)x)-mu_x+0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeX()) / (sqrt(2.)*sigma_x)  )
                                     - erf(  ( (((FLTNB)x)-mu_x-0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeX()) / (sqrt(2.)*sigma_x)  ) )
                                   * kern_ZY;
        sum_kernel += m2p_kernel[0][index];
      }
    }
  }
  // Normalize kernel
  for (int v=0; v<m_dimKernelXYZ; v++) m2p_kernel[0][v] /= sum_kernel;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
