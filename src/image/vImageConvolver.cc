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
  \brief    Implementation of class vImageConvolver
*/

#include "vImageConvolver.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oImageSpace.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vImageConvolver::vImageConvolver() 
{
  // Standards
  mp_ImageDimensionsAndQuantification = NULL;
  m_verbose = -1;
  // Booleans
  m_checked = false;
  m_initialized = false;
  m_stationary = false;
  // Padded image
  mp_paddedImage = NULL;
  m_offsetX = -1;
  m_offsetY = -1;
  m_offsetZ = -1;
  m_dimPadX = -1;
  m_dimPadY = -1;
  m_dimPadZ = -1;
  m_dimPadXY = -1;
  m_dimPadXYZ = -1;
  // Convolution kernel
  m_nbKernels = -1;
  mp_dimKernelX = NULL;
  mp_dimKernelY = NULL;
  mp_dimKernelZ = NULL;
  m2p_kernel = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vImageConvolver::~vImageConvolver() 
{
  // Simply free what should be
  if (mp_paddedImage) free(mp_paddedImage);
  if (mp_dimKernelX) free(mp_dimKernelX);
  if (mp_dimKernelY) free(mp_dimKernelY);
  if (mp_dimKernelZ) free(mp_dimKernelZ);
  if (m2p_kernel)
  {
    for (int k=0; k<m_nbKernels; k++) if (m2p_kernel[k]) free(m2p_kernel[k]);
    free(m2p_kernel);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

// Check global parameters then check the specific parameters
int vImageConvolver::CheckParameters()
{
  // Check verbose
  if (m_verbose<0)
  {
    Cerr("***** vImageConvolver::CheckParameters() -> The verbose level must be set positive !" << endl);
    return 1;
  }
  // Check oImageDimensionsAndQuantification
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** vImageConvolver::CheckParameters() -> The image dimensions must be set through a oImageDimensionsAndQuantification !" << endl);
    return 1;
  }
  // Finally check the specific parameters from the child convolver
  if (CheckSpecificParameters())
  {
    Cerr("***** vImageConvolver::CheckParameters() -> A problem occured while checking specific parameters of the child convolver !" << endl);
    return 1;
  }
  // All set
  m_checked = true;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

// Check that the kernel has been successfully built, then allocate the padded image
int vImageConvolver::Initialize()
{
  // Check if parameters have been checked
  if (!m_checked)
  {
    Cerr("***** vImageConvolver::Initialize() -> Parameters have not been checked ! Please call CheckParameters() before." << endl);
    return 1;
  }

  // Build the convolution kernel
  if (BuildConvolutionKernel())
  {
    Cerr("***** vImageConvolver::Initialize() -> A problem occured while building the convolution kernel !" << endl);
    return 1;
  }

  // ---------------------------------------------------------------------------------------
  // Check that all kernel dimensions are odd and find the maximum size along each dimension
  // ---------------------------------------------------------------------------------------

  // Check that the number of kernels is not null or negative
  if (m_nbKernels<1)
  {
    Cerr("***** vImageConvolver::Initialize() -> The number of kernels is negative or null !" << endl);
    return 1;
  }
  // Check that the kernel dimensions have been allocated
  if (mp_dimKernelX==NULL || mp_dimKernelY==NULL || mp_dimKernelZ==NULL)
  {
    Cerr("***** vImageConvolver::Initialize() -> The kernel dimensions are not allocated !" << endl);
    return 1;
  }
  // Check that the kernel has been allocated
  if (m2p_kernel==NULL)
  {
    Cerr("***** vImageConvolver::Initialize() -> The kernel is not allocated !" << endl);
    return 1;
  }
  // In case the kernel is stationary and the user forgot to explicitely specify it, we force it here
  if (m_nbKernels==1) m_stationary = true;

  // Maximum kernel sizes along each dimension
  INTNB max_kern_size_Z = 0;
  INTNB max_kern_size_Y = 0;
  INTNB max_kern_size_X = 0;

  // Loop over all kernels
  for (int k=0; k<m_nbKernels; k++)
  {
    // Check oddity along Z
    if (mp_dimKernelZ[k]%2==0)
    {
      Cerr("***** vImageConvolver::Initialize() -> Dimension of " << k << "th convolution kernel along Z (" << mp_dimKernelZ[k] << " is not odd !" << endl);
      return 1;
    }
    // Check oddity along Y
    if (mp_dimKernelY[k]%2==0)
    {
      Cerr("***** vImageConvolver::Initialize() -> Dimension of " << k << "th convolution kernel along Y (" << mp_dimKernelY[k] << " is not odd !" << endl);
      return 1;
    }
    // Check oddity along X
    if (mp_dimKernelX[k]%2==0)
    {
      Cerr("***** vImageConvolver::Initialize() -> Dimension of " << k << "th convolution kernel along X (" << mp_dimKernelX[k] << " is not odd !" << endl);
      return 1;
    }
    // Check maximum along Z
    if (mp_dimKernelZ[k]>max_kern_size_Z) max_kern_size_Z = mp_dimKernelZ[k];
    // Check maximum along Y
    if (mp_dimKernelY[k]>max_kern_size_Y) max_kern_size_Y = mp_dimKernelY[k];
    // Check maximum along X
    if (mp_dimKernelX[k]>max_kern_size_X) max_kern_size_X = mp_dimKernelX[k];
  }

  // Set the offsets
  m_offsetZ = max_kern_size_Z / 2;
  m_offsetY = max_kern_size_Y / 2;
  m_offsetX = max_kern_size_X / 2;

  // In the convolution we assume an intensity conservation along the Z axis, so we check if the offset is higher or equal to the image dimension.
  // In this case, the convolution will still work (it is implemented so), but computation time will be lost for nothing, so we throw a warning.
  if (m_offsetZ>=mp_ImageDimensionsAndQuantification->GetNbVoxZ())
  {
    FLTNB over_computation_factor = ((FLTNB)(max_kern_size_Z)) / ((FLTNB)(mp_ImageDimensionsAndQuantification->GetNbVoxZ()*2-1));
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!! vImageConvolver::Initialize() -> The maximum half size of the convolution kernels is higher or equal to the number of slices of the image." << endl);
    Cerr("!!!!!                                  This won't affect the results, but a convolution running time factor of " << over_computation_factor << " will be lost for nothing." << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
  }

  // ---------------------------------------------------------------------------------------
  // Build the padded image
  // ---------------------------------------------------------------------------------------

  // Compute the padded image dimensions
  m_dimPadZ = mp_ImageDimensionsAndQuantification->GetNbVoxZ() + max_kern_size_Z - 1;
  m_dimPadY = mp_ImageDimensionsAndQuantification->GetNbVoxY() + max_kern_size_Y - 1;
  m_dimPadX = mp_ImageDimensionsAndQuantification->GetNbVoxX() + max_kern_size_X - 1;
  m_dimPadXY = m_dimPadX * m_dimPadY;
  m_dimPadXYZ = m_dimPadXY * m_dimPadZ;

  // Allocate the padded image and fill it with 0.
  mp_paddedImage = (FLTNB*)calloc(m_dimPadXYZ,sizeof(FLTNB));

  // All set
  m_initialized = true;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vImageConvolver::ApplyConvolution(FLTNB* ap_image)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vImageConvolver::ApplyConvolution() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  // Copy image into padded buffer image
  CopyToPaddedImage(ap_image);
  // Then convolve
  if (Convolve(ap_image))
  {
    Cerr("***** vImageConvolver::ApplyConvolution() -> A problem occured while actually convolving !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vImageConvolver::ApplyConvolutionTranspose(FLTNB* ap_image)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vImageConvolver::ApplyConvolutionTranspose() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  // Copy image into padded buffer image
  CopyToPaddedImage(ap_image);
  // Then convolve
  if (ConvolveTranspose(ap_image))
  {
    Cerr("***** vImageConvolver::ApplyConvolutionTranspose() -> A problem occured while actually convolving !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vImageConvolver::CopyToPaddedImage(FLTNB* ap_inputImage)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vImageConvolver::CopyToPaddedImage() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  // Copy the image into the padded buffer, centered with respect to the pad
  INTNB z;
  #pragma omp parallel for private(z) schedule(guided)
  for (z=0; z<mp_ImageDimensionsAndQuantification->GetNbVoxZ(); z++)
  {
    // Pad coordinate Z
    INTNB padz = z + m_offsetZ;
    for (INTNB y=0; y<mp_ImageDimensionsAndQuantification->GetNbVoxY(); y++)
    {
      // Pad coordinate Y
      INTNB pady = y + m_offsetY;
      for (INTNB x=0; x<mp_ImageDimensionsAndQuantification->GetNbVoxX(); x++)
      {
        // Pad coordinate X
        INTNB padx = x + m_offsetX;
        // Global original coordinate
        INTNB coord_orig = z*mp_ImageDimensionsAndQuantification->GetNbVoxXY() + y*mp_ImageDimensionsAndQuantification->GetNbVoxX() + x;
        // Global pad coordinate
        INTNB coord_pad = padz*m_dimPadXY + pady*m_dimPadX + padx;
        // Affect buffer
        mp_paddedImage[coord_pad] = ap_inputImage[coord_orig];
        // Zero image
        ap_inputImage[coord_orig] = 0.;
      }
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vImageConvolver::Convolve(FLTNB* ap_outputImage)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vImageConvolver::Convolve() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  // This function is designed to work universally and exclusively on a stationary kernel.
  // Note that we suppose that a stationary kernel is also symmetric (but not necessarily isotropic).
  // For asymetric and non-stationary kernels, have to overload this function in the child convolver.
  // Note also that we choose the out-to-in implementation (outer loop on the output image) so that the
  // parallelism is operated on the output and is thus thread-safe.
  // Finally note that transaxially, we do not implement an image intensity preservation policy as we suppose that
  // the distance between the object and the image border is taller than the kernel half-width (it should be so).
  if (m_stationary)
  {
    INTNB stationary_kernel = 0;
    // Convolve (OpenMP on Z dimension)
    INTNB z;
    #pragma omp parallel for private(z) schedule(guided)
    for (z=0; z<mp_ImageDimensionsAndQuantification->GetNbVoxZ(); z++)
    {
      // As a first step, we virtually go through the kernel along Z to compute the actual kernel integral that will
      // be really used in the convolution. For each virtual slice contributing to the convolved slice, we test if it
      // is inside or outside the image; in the former case we count it into the integral and in the latter we don't.
      FLTNB kernel_integral = 0.;
      for (INTNB zk=0; zk<mp_dimKernelZ[stationary_kernel]; zk++)
      {
        // Index of the virtual slice which will (or not) contribute to the convolved slice
        INTNB virtualZ = z + zk - m_offsetZ;
        // Test if it is not in the image, then we move onto the next virtual slice
        if (virtualZ<0 || virtualZ>=mp_ImageDimensionsAndQuantification->GetNbVoxZ()) continue;
        // Otherwise, we update the kernel integral
        INTNB kernZ_base = zk * mp_dimKernelX[stationary_kernel] * mp_dimKernelY[stationary_kernel];
        for (INTNB xyk=0; xyk<mp_dimKernelX[stationary_kernel]*mp_dimKernelY[stationary_kernel]; xyk++)
          kernel_integral += m2p_kernel[stationary_kernel][kernZ_base+xyk];
      }
      // Some precomputation for maximum efficiency
      INTNB indexZ_base = z*mp_ImageDimensionsAndQuantification->GetNbVoxXY();
      for (INTNB y=0; y<mp_ImageDimensionsAndQuantification->GetNbVoxY(); y++)
      {
        // Some precomputation for maximum efficiency
        INTNB indexZY_base = indexZ_base + y*mp_ImageDimensionsAndQuantification->GetNbVoxX();
        for (INTNB x=0; x<mp_ImageDimensionsAndQuantification->GetNbVoxX(); x++)
        {
          // Output image index
          INTNB index_img = indexZY_base + x;
          // Inner loops on convolution kernel dimensions
          for (INTNB zk=0, index_kern=0; zk<mp_dimKernelZ[stationary_kernel]; zk++)
          {
            // Some precomputation for maximum efficiency
            INTNB indexZ_pad = (z + zk) * m_dimPadXY;
            for (INTNB yk=0; yk<mp_dimKernelY[stationary_kernel]; yk++)
            {
              // Some precomputation for maximum efficiency
              INTNB indexZY_pad = indexZ_pad + (y + yk) * m_dimPadX;
              for (INTNB xk=0; xk<mp_dimKernelX[stationary_kernel]; xk++, index_kern++)
              {
                // Padded image index
                INTNB index_pad = indexZY_pad + x + xk;
                // Apply contribution
                ap_outputImage[index_img] += mp_paddedImage[index_pad] * m2p_kernel[stationary_kernel][index_kern] / kernel_integral;
              }
            }
          }
        }
      }
    }
    // Normal end
    return 0;
  }
  // Otherwise the convolve function has to be implemented by the child non-stationary convolver
  else
  {
    Cerr("***** vImageConvolver::Convolve() -> To use a non-stationary kernel, you should overload this function to implement your own in your child convolver !" << endl);
    return 1;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vImageConvolver::ConvolveTranspose(FLTNB* ap_outputImage)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vImageConvolver::ConvolveTranspose() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  // This function is designed to work universally on a stationary kernel.
  // Note that we suppose a stationary kernel is also symmetric (otherwise it is absurd)
  if (m_stationary)
  {
    return Convolve(ap_outputImage);
  }
  // Otherwise the convolve function has to be implemented by the child non-stationary convolver
  else
  {
    Cerr("***** vImageConvolver::ConvolveTranspose() -> To use a non-stationary kernel, you should overload this function to implement your own in your child convolver !" << endl);
    return 1;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
