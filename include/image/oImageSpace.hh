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
  \brief    Declaration of class oImageSpace
*/

#ifndef OIMAGESPACE_HH
#define OIMAGESPACE_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oDeformationManager.hh"
#include "oInterfileIO.hh"

class vOptimizer;

/*!
  \class   oImageSpace
  \brief   This class holds all the matrices in the image domain that can be used in the algorithm: image, forward-image, correction, additional image, sensitivity image.
  \details Image matrices are public and can be directly accessed from each classes. \n
           It also includes many functions for initializating, reseting, or deforming the images. \n
           Mandatory Reconstruction image matrices : \n
           Image:                  4 pointers: 1: dynamic frames, 2: respiratory gates, 3: cardiac gates, 4: 3D voxels \n
           ForwardImage:           4 pointers: 1: dynamic frames, 2: respiratory gates, 3: cardiac gates, 4: 3D voxels \n
           BackwardImage:          6 pointers: 1: number of images (for optimizer), 2: threads, 3: dynamic frames, 4: respiratory gates, 5: cardiac gates, 6: 3D voxels \n
           SensitivityImage:       5 pointers: 1: threads, 2: dynamic frames, 3: respiratory gates, 4: cardiac gates, 5: 3D voxels \n
           OutputImage:            4 pointers: 1: dynamic frames, 2: respiratory gates, 3: cardiac gates, 4: 3D voxels \n
           VisitedVoxelsImage:     1 pointer: 1: 3D voxels \n
           Optional Reconstruction image matrices : \n
           MultiModalImage:        5 pointers: 1: number of multimodal images, 2: dynamic frames, 3: respiratory gates, 4: cardiac gates, 5: 3D voxels \n
           DefTmpBackwardImage:    6 pointers: 1: number of images (for optimizer), 2: threads, 3: dynamic frames, 4: respiratory gates, 5: cardiac gates, 6: 3D voxels \n
           DefTmpSensitivityImage: 5 pointers: 1: threads, 2: dynamic frames, 3: respiratory gates, 4: cardiac gates, 5: 3D voxels \n
           Projection image matrices : \n
           Attenuation:            4 pointers: 1: dynamic frames, 2: respiratory gates, 3: cardiac gates, 4: 3D voxels: \n
           ProjectionImage:        2 pointers: 1: number of projections, 2: 2D pixels 
*/
class oImageSpace
{
  // Constructor & Destructor
  public:
    /*!
      \brief   oImageSpace constructor. 
               Initialize the member variables to their default values.
    */
    oImageSpace();
    /*!
      \brief   oImageSpace destructor. 
    */
    ~oImageSpace();


  // -------------------------------------------------------------------
  // Public member functions
  public:

  // Let these pointers public to simplify getting them in other classes
  FLTNB****   m4p_image; /*!< Dynamic array for the reconstructed image, containing a 4-vectors voxellized image matrix.
                              4 pointers: 
                              1: dynamic frames, 
                              2: respiratory gates, 
                              3: cardiac gates, 
                              4: 3D voxels */
                              
  FLTNB****   m4p_forwardImage; /*!< Dynamic array for the forward image (image data before the projection step), containing a 4-vectors voxellized image matrix.
                                      4 pointers: 
                                      1: dynamic frames, 
                                      2: respiratory gates, 
                                      3: cardiac gates, 
                                      4: 3D voxels */
                                      
  FLTNB****** m6p_backwardImage;  /*!< Dynamic array for the backward image (image-based correction factors computed after the projection/backprojection steps),
                                      containing a 6-vectors voxellized image matrix.
                                       6 pointers: 
                                       1: number of images (dedicated to the optimization algorithm)
                                       2: threads, 
                                       3: dynamic frames, 
                                       4: respiratory gates, 
                                       5: cardiac gates, 
                                       6: 3D voxels */
                                       
  FLTNB*****  m5p_sensitivity; /*!< Dynamic array for the sensitivity image (image-based normalization factors), containing a 5-vectors voxellized image matrix.
                                    5 pointers: 
                                    1: threads, 
                                    2: dynamic frames, 
                                    3: respiratory gates, 
                                    4: cardiac gates, 
                                    5: 3D voxels */
                                     
  FLTNB*****   m5p_multiModalImage; /*!< Dynamic array for multimodal images to be used in reconstruction, containing a 5-vectors voxellized image matrix.
                                        5 pointers:
                                        1: number of multimodal images
                                        2: dynamic frames, 
                                        3: respiratory gates, 
                                        4: cardiac gates, 
                                        5: 3D voxels */

  FLTNB*      mp_maskImage; /*!< Array for a mask image to be used in reconstruction, for instance a background mask.
                                        1 pointer:
                                        1: 3D voxels */
                                        
  /*FLTNB*      mp_visitedVoxelsImage; !< 1-D Dynamic array, containing binary information regarding which 3D voxels have been visited during the projection steps.
                                          1 pointer:
                                          3D voxels */
  FLTNB****   m4p_visitedVoxelsImage; /*!< Dynamic array, containing binary information regarding which 3D voxels have been visited during the projection steps.
                                          4 pointers: 
                                          1: dynamic frames, 
                                          2: respiratory gates, 
                                          3: cardiac gates, 
                                          4: 3D voxels */
  FLTNB**     m2p_miscellaneousImage; /*!< A bunch of images that can be allocated on demand using the AllocateMiscellaneousImage() function which return
                                           a pointer to the newly created image. They can thus be used anywhere in the code for any purpose. */
                                     
  FLTNB****   m4p_attenuation; /*!< Dynamic array for the attenuation images, used in analytical projection or sensitivity list-mode generation, containing a 4-vectors voxellized image matrix.
                                    4 pointers: 
                                    1: dynamic frames, 
                                    2: respiratory gates, 
                                    3: cardiac gates, 
                                    4: 3D voxels */
                                    
  FLTNB****   m4p_outputImage; /*!< Dynamic array for the output image to be written on disk, containing a 4-vectors voxellized image matrix
                                    Required in the case the reconstruction uses intrinsic temporal basis functions.
                                    In any other case, m4p_outputImage just refers to the address of the image matrix containing the regular image.
                                    4 pointers: 
                                    1: dynamic frames, 
                                    2: respiratory gates, 
                                    3: cardiac gates, 
                                    4: 3D voxels */

  FLTNB****   m4p_refDynForwardImage;  /*!< Buffer dynamic array for the forward image (image-based correction factors computed after the projection/backprojection steps), 
                                            containing a 5-vectors voxellized image matrix.
                                             Required for image-based deformation, in order to store the forward image of the reference position from which deformation are performed
                                             5 pointers: 
                                             1: number of images (dedicated to the optimization algorithm)
                                             2: dynamic frames, 
                                             3: respiratory gates, 
                                             4: cardiac gates, 
                                             5: 3D voxels */
                                             
  FLTNB*****  m5p_refDynBackwardImage;  /*!< Buffer dynamic array for the backward image (image-based correction factors computed after the projection/backprojection steps), 
                                            containing a 6-vectors voxellized image matrix.
                                             Required when image-based deformation, in order to store the correction factors in the reference position
                                             6 pointers: 
                                             1: number of images (dedicated to the optimization algorithm)
                                             2: threads, 
                                             3: dynamic frames, 
                                             4: respiratory gates, 
                                             5: cardiac gates, 
                                             6: 3D voxels */
                                             
  FLTNB****   m4p_refDynSensitivityImage;  /*!< Buffer dynamic array for the sensitivity image (image-based sensitivity factors), containing a 5-vectors voxellized image matrix.
                                             Required when image-based deformation is enabled during histogram-based reconstruction, in order to store the sensitivity factors in the reference position
                                             5 pointers:
                                             1: threads, 
                                             2: dynamic frames, 
                                             3: respiratory gates, 
                                             4: cardiac gates, 
                                             5: 3D voxels */
                                             
  FLTNB**     m2p_projectionImage;  /*!< Dynamic array for a projection image (currently dedicated to SPECT analytical projection), containing a 2-vectors voxellized image matrix
                                    Required in the case the reconstruction uses intrinsic temporal basis functions. In any other case, m4p_outputImage just refers to the address of the image matrix containing the regular image.
                                    2 pointers:
                                    1: number of projections
                                    2: 2D pixels */


  // -------------------------------------------------------------------
  /*!
    \fn      oImageSpace::InstantiateImage
    \brief   Allocate memory for the main image matrices
  */
  void InstantiateImage();
  /*!
    \fn      oImageSpace::DeallocateImage
    \brief   Free memory for the main image matrices
  */
  void DeallocateImage();
  /*!
    \fn      oImageSpace::InstantiateForwardImage
    \brief   Allocate memory for the forward image matrices
    \details The dimensions are taken from the dynamic basis (number of time basis, respiratory and cardiac basis functions, as opposed to the number of frames/gates)
  */
  void InstantiateForwardImage();
  /*!
    \fn      oImageSpace::DeallocateForwardImage
    \brief   Free memory for the forward image matrices
    \details The dimensions are taken from the dynamic basis (number of time basis, respiratory and cardiac basis functions, as opposed to the number of frames/gates)
  */
  void DeallocateForwardImage();
  /*!
    \fn      oImageSpace::InstantiateBackwardImageFromDynamicBasis
    \param   a_nbBackwardImages : number of backward images required for the optimization algorithm
    \brief   Allocate memory for the backward image matrices and set the number of backward images for the whole class
    \details The dimensions are taken from the dynamic basis (number of time basis, respiratory and cardiac basis functions, as opposed to the number of frames/gates)
  */
  void InstantiateBackwardImageFromDynamicBasis(int a_nbBackwardImages);
  /*!
    \fn      oImageSpace::DeallocateBackwardImageFromDynamicBasis
    \brief   Free memory for the backward image matrices
    \details The dimensions are taken from the dynamic basis (number of time basis, respiratory and cardiac basis functions, as opposed to the number of frames/gates)
  */
  void DeallocateBackwardImageFromDynamicBasis();
  /*!
    \fn      oImageSpace::InstantiateBackwardImageFromDynamicBins
    \brief   Allocate memory for the backward image matrices and initialize them
    \details The dimensions are taken from the dynamic bins (number of times frames, respiratory and cardiac gates, as opposed to the number of basis functions)
  */
  void InstantiateBackwardImageFromDynamicBins();
  /*!
    \fn      oImageSpace::DeallocateBackwardImageFromDynamicBins
    \brief   Free memory of the backward image matrices
    \details The dimensions are taken from the dynamic bins (number of times frames, respiratory and cardiac gates, as opposed to the number of basis functions)
  */
  void DeallocateBackwardImageFromDynamicBins();
  /*!
    \fn      oImageSpace::InstantiateSensitivityImage
    \param   a_pathToSensitivityImage : path to the sensitivity image
    \brief   Allocate the sensitivity image matrices
    \details Sensitivity image initialization depends on the reconstruction mode :
             - list-mode: Sensitivity image has been computed before reconstruction, it is loaded from the path provided in parameter
             - histogram: Sensitivity image is calculated on the fly during reconstruction.  
             First dimension (thread) is only used in histogram mode, as the on-the-fly sensitivity image computation must be thread safe
  */
  void InstantiateSensitivityImage(const string& a_pathToSensitivityImage);
  /*!
    \fn      oImageSpace::DeallocateSensitivityImage
    \brief   Free memory for the sensitivity image matrices
    \details Sensitivity image deallocation depends on the reconstruction mode (multithreaded in histogram but not in list-mode)
  */
  void DeallocateSensitivityImage();
  /*!
    \fn      oImageSpace::AllocateMiscellaneousImage
    \brief   Allocate a new miscellaneous image on m2p_miscellaneousImages and return the pointer to this image
    \return  The pointer to the newly allocated image
  */
  FLTNB* AllocateMiscellaneousImage();
  /*!
    \fn      oImageSpace::DeallocateMiscellaneousImage
    \brief   Deallocate all allocated miscellaneous images
  */
  void DeallocateMiscellaneousImage();
  /*!
    \fn      oImageSpace::InitMultiModalImage
    \param   a_pathToMultiModalImage : path to multimodal image
    \brief   Memory allocation and initialization for the multimodal image matrices
    \details Nothing is performed if the path provided in parameter is empty; the image is directly read
    \return  0 if success, positive value otherwise
  */
  int InitMultiModalImage(const vector<string>& a_pathToMultiModalImage);
  /*!
    \fn      oImageSpace::DeallocateMultiModalImage
    \brief   Free memory for the multimodal image
  */
  void DeallocateMultiModalImage();
  /*!
    \fn      oImageSpace::InitMaskImage
    \param   a_pathToImage : path to the mask image
    \brief   Memory allocation and initialization for the mask image
    \details Nothing is performed if the path provided in parameter is empty; the image is directly read
    \return  0 if success, positive value otherwise
  */
  int InitMaskImage(const string& a_pathToImage);
  /*!
    \fn      oImageSpace::DeallocateMaskImage
    \brief   Free memory for the mask image
  */
  void DeallocateMaskImage();


  /*!
    \fn      oImageSpace::InstantiateOutputImage
    \brief   Instanciate Image matrices dedicated to output writing on disk
    \details Additionnal output image matrix is needed if the reconstruction uses intrinsic temporal basis functions
             In this case, the image matrices are defined in the temporal image basis functions space, therefore requiring an additional step to recover the images in the regular image-space.
             If no intrinsic temporal basis functions are used, m4p_outputImage just refers to the address of the image matrix containing the regular image.
  */
  void InstantiateOutputImage();
  /*!
    \fn      oImageSpace::InstantiateRefImageForDeformation
    \brief   Memory allocation for the reference forward/backward image required for image-based deformation. This function is called from the Deformation Manager
  */
  /*!
    \fn      oImageSpace::InstantiateRefImagesForDeformation
    \brief   Allocate memory for the buffer sensitivity image required for image-based deformation. This function is called from the Deformation Manager
  */
  void InstantiateRefImagesForDeformation();
  /*!
    \fn      oImageSpace::InstantiateSensImageForDeformation
    \brief   Memory allocation for the buffer sensitivity image required for image-based deformation (only for PET HISTOGRAM mode)
  */
  //void InstantiateSensImageForDeformation();
  /*!
    \fn      oImageSpace::InstantiateVisitedVoxelsImage
    \brief   Memory allocation and initialization for the image matrix containing binary information regarding which 3D voxels have been visited during the projection steps.
  */
  void InstantiateVisitedVoxelsImage();
  /*!
    \fn      oImageSpace::DeallocateBwdImageForDeformation
    \brief   Free memory for the buffer backward image required for image-based deformation
  */
  //void DeallocateBwdImageForDeformation();
  /*!
    \fn      oImageSpace::DeallocateRefImagesForDeformation
    \brief   Free memory for the buffer sensitivity image required for image-based deformation. This function is called from the Deformation Manager
  */
  void DeallocateRefImagesForDeformation();
  /*!
    \fn      oImageSpace::DeallocateOutputImage
    \brief   Free memory for the Image matrices dedicated to output writing on disk
  */
  void DeallocateOutputImage();
  /*!
    \fn      oImageSpace::DeallocateVisitedVoxelsImage
    \brief   Free memory for the image matrix containing binary information regarding which 3D voxels have been visited during the projection steps.
  */
  void DeallocateVisitedVoxelsImage();
  /*!
    \fn      oImageSpace::InitImage
    \param   a_pathToInitialImage : path to an existing image
    \param   a_value : value to initialize each voxel with, if an input image is not provided
    \brief   Initialize the main image, either using:
             - an existing image at path 'a_pathToInitialImage'
             - initialize each voxel with 'a_value'.
    \return  0 if success, positive value otherwise
  */
  int InitImage(const string& a_pathToInitialImage, FLTNB a_value);
  /*!
    \fn      oImageSpace::InitBackwardImage
    \brief   Initialize each voxel of the backward images to 0, also for sensitivity if not loaded (estimated on the fly).
  */
  void InitBackwardImage();
  /*!
    \fn      oImageSpace::InitSensitivityImage()
    \param   a_pathToSensitivityImage : path to the sensitivity image (should be provided only in list-mode reconstruction)
    \brief   Initialization for the sensitivity image matrices
    \details Sensitivity image initialization depends on the reconstruction mode :
             - list-mode: Sensitivity image has been computed before reconstruction, it is loaded from the path provided in parameter
             - histogram: Sensitivity image is calculated on the fly during reconstruction.  
             First dimension (thread) is only used in histogram mode, as the on-the-fly sensitivity image computation must be thread safe
    \return  0 if success, positive value otherwise
  */
  int InitSensitivityImage(const string& a_pathToSensitivityImage);
  /*!
    \fn      oImageSpace::InitBwdImageForDeformation
    \brief   Initialize the buffer backward image dedicated to image-based deformation
    \details this function is used by the vDeformation class, in order to reset the backward image 
             after its content has been recovered in the buffer backward image (typically, after a deformation step)
  */
  //void InitBwdImageForDeformation();
  /*!
    \fn      oImageSpace::InitSensImageForDeformation
    \brief   Initialize the buffer sensitivity image dedicated to image-based deformation, if required (histogram mode, as sensitivity is not loaded)
  */
  //void InitSensImageForDeformation();
  /*!
    \fn      oImageSpace::InitRefImageForDeformation
    \brief   Initialize the references image dedicated to image-based deformation. This function is called from the Deformation Manager
  */
  void InitRefImagesForDeformation();
  /*!
    \fn      oImageSpace::LoadInitialImage
    \param   a_pathToImage : path to an existing image
    \brief   Load the initial image provided by the user in the corresponding matrix
    \return  0 if success, positive value otherwise
  */
  int LoadInitialImage(const string& a_pathToImage);
  /*!
    \fn      oImageSpace::InstantiateOutputImage
    \brief   Compute output image using the m4p_image matrix and the time/respiratory/cardiac basis functions. Store the result in the m4p_outputImage matrix
    \details If time/respiratory/cardiac basis functions have been initialized, this function has no effect.
  */
  void ComputeOutputImage();
  /*!
    \fn      oImageSpace::ApplyOutputFOVMasking
    \brief   Mask the outside of the transaxial FOV based on the m_fovOutPercent
    \details Do this on m4p_outputImage
  */
  int ApplyOutputFOVMasking();
  /*!
    \fn      oImageSpace::ApplyOutputMaskImage
    \brief   Mask the outside of the provided input mask image
    \details Do this on m4p_outputImage
  */
  int ApplyOutputMaskImage();
  /*!
    \fn      oImageSpace::ApplyMaskToSensitivity
    \brief   Apply the mask to the sensitivity image (only for the first thread, the image must be reduced beforehand)
    \details Do this on m5p_sensitivityImage
  */
  int ApplyMaskToSensitivity();
  /*!
    \fn      oImageSpace::ApplyOutputFlip()
    \brief   Just flip the output image
    \details Do this on m4p_outputImage
  */
  int ApplyOutputFlip();
  /*!
    \fn      oImageSpace::SaveOutputImage
    \param   a_iteration : current iteration index
    \param   a_subset : current number of subsets (or -1 by default)
    \brief   Call the interfile function to write output image on disk
    \return  0 if success, positive value otherwise
  */
  int SaveOutputImage(int a_iteration, int a_subset = -1);
  /*!
    \fn      oImageSpace::SaveDebugImage
    \param   a_name : output name of the image
    \brief   Just a debug function dedicated to write any kind of image on disk in raw format, for debugging purposes
  */
  void SaveDebugImage(const string& a_name);
  /*!
    \fn SaveSensitivityImage
    \param a_pathToSensitivityImage : path to the sensitivity image (should be provided only in list-mode reconstruction)
    \brief Call the interfile function to write the sensitivity image on disk
    \return  0 if success, positive value otherwise
  */
  int SaveSensitivityImage(const string& a_pathToSensitivityImage);
  /*!
    \fn      oImageSpace::PrepareForwardImage
    \brief   Copy current image matrix in the forward-image buffer matrix
  */
  void PrepareForwardImage();
  /*!
    \fn      oImageSpace::Reduce
    \brief   Merge parallel results into the matrix of the backward image matrix of the first thread. Also for MPI.
  */
  void Reduce();
  /*!
    \fn      oImageSpace::CleanNeverVisitedVoxels
    \brief   Based on the visitedVoxelsImage, clean the never visited voxels in the image.
             This function must be called at the end of each iteration
  */
  void CleanNeverVisitedVoxels();


  // -------------------------------------------------------------------
  // Functions for List-Mode Sensitivity Generation process
  
  /*!
    \fn      oImageSpace::LMS_InstantiateImage
    \brief   Allocate memory for the main image matrices (for list-mode sensitivity generation)
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_InstantiateImage();
  /*!
    \fn      oImageSpace::LMS_InstantiateForwardImage
    \brief   Allocate memory for the forward image matrices (for list-mode sensitivity generation)
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_InstantiateForwardImage();
  /*!
    \fn      oImageSpace::LMS_InstantiateSensitivityImage
    \brief   Allocate memory for the sensitivity image matrices (for list-mode sensitivity generation)
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_InstantiateSensitivityImage();
  /*!
    \fn      oImageSpace::LMS_DeallocateImage
    \brief   Free memory for the main image matrices (for list-mode sensitivity generation)
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_DeallocateImage();
  /*!
    \fn      oImageSpace::LMS_DeallocateForwardImage
    \brief   Free memory for the forward image matrices (for list-mode sensitivity generation)
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_DeallocateForwardImage();
  /*!
    \fn      oImageSpace::LMS_DeallocateSensitivityImage
    \brief   Free memory for the sensitivity image matrices (for list-mode sensitivity generation)
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_DeallocateSensitivityImage();
  /*!
    \fn      oImageSpace::LMS_DeallocateAttenuationImage
    \brief   Free memory for the Attenuation image matrices (for analytical projection or list-mode sensitivity generation)
  */
  void LMS_DeallocateAttenuationImage();
  /*!
    \fn      oImageSpace::InitAttenuationImage
    \param   a_pathToAtnImage : path to an existing image
    \brief   Memory allocation and initialisation for the attenuation image using either :
             - an image provided by the user (a_pathToAtnImage or a_pathToCTImage)
             - the default value (=a_value)
    \return  0 if success, positive value otherwise
    \todo    If fr/rg/cg > 1 and only one attenuation image is provided, initialize each redundant matrice with pointers to the first one (deal with that here)
  */
  int InitAttenuationImage(const string& a_pathToAtnImage);
  /*!
    \fn      oImageSpace::LMS_CopyAtnToImage
    \brief   Copy the attenuation image contained in the 'm2p_attenuation' matrix inside the m2p_image matrix.
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_CopyAtnToImage();
  /*!
    \fn      oImageSpace::LMS_CopyAtnToForwardImage
    \param   a_use1stMotion
    \param   a_use2ndMotion
    \brief   Copy the attenuation image contained in the 'm2p_attenuation' matrix inside the m4p_forwardImage matrix.
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_CopyAtnToForwardImage(bool a_use1stMotion, bool a_use2ndMotion);
  /*!
    \fn      oImageSpace::LMS_CopyBackwardToSensitivityImage()
    \brief   Copy the backward image containing the result of the sensitivity back-projection, into the sensitivity matrix.
    \details This function is dedicated to list-mode sensitivity (LMS) generation, it is useful to apply convolution and to save the image.
  */
  void LMS_CopyBackwardToSensitivity();
  /*!
    \fn      oImageSpace::LMS_PrepareForwardImage
    \brief   Copy current image in forward-image buffer (for list-mode sensitivity generation)
    \details This function is dedicated to list-mode sensitivity (LMS) generation.
  */
  void LMS_PrepareForwardImage();
  /*!
    \fn      oImageSpace::ReduceBackwardImage
    \param   a_imageIndex
    \param   a_timeIndex
    \param   a_respIndex
    \param   a_cardIndex
    \brief   Merge parallel results into the backward image matrix of the first thread for the specific image / time / respiratory / cardiac provided indices
  */
  void ReduceBackwardImage(int a_imageIndex, int a_timeIndex, int a_respIndex, int a_cardIndex);
  /*!
    \fn      oImageSpace::ApplyMaskToBackwardImage
    \brief   Apply the mask to the backward image matrix of the first thread for the specific image / time / respiratory / cardiac provided indices
    \details Do this on m6p_backwardImage
  */
  int ApplyMaskToBackwardImage(int a_imageIndex, int a_timeIndex, int a_respIndex, int a_cardIndex);
  /*!
    \fn      oImageSpace::LMS_SaveSensitivityImage
    \param   a_pathToSensitivityImage : path to the sensitivity image (should be provided only in list-mode reconstruction)
    \param   ap_DeformationManager : Pointer to the deformation manager objet (required to retrieve the number of gates in the sensitivity image) 
    \brief   Call the interfile function to write the sensitivity image on disk
    \details If image deformation is enabled for respiratory/cardiac gated data, the gated images are summed up into one image and normalize
    \return  0 if success, positive value otherwise
  */
  int LMS_SaveSensitivityImage(const string& a_pathToSensitivityImage, oDeformationManager* ap_DeformationManager);


  // -------------------------------------------------------------------
  // Functions for ANALYTICAL PROJECTION process
  
  /*!
    \fn      oImageSpace::PROJ_InstantiateProjectionImage
    \param   a_nbProjs : a number of projection slices in the projection
    \param   a_nbPixels : a total number of pixels in the projection slices
    \brief   Instanciate and initialize projection image for analytical projection
    \details This function is currently only dedicated to SPECT projection, and used by the analytical projection script
    \todo    support for PET (requires projection format)
  */
  void PROJ_InstantiateProjectionImage(int a_nbProjs, int a_nbPixels);
  /*!
    \fn      oImageSpace::PROJ_DeallocateProjectionImage
    \param   a_nbProjs : a number of projection slices in the projection
    \brief   Free memory for the projection image for analytical projection
    \details This function is currently only dedicated to SPECT projection, and used by the analytical projection script
    \todo    support for PET (requires projection format)
  */
  void PROJ_DeallocateProjectionImage(int a_nbProjs);
  /*!
    \fn      oImageSpace::PROJ_InitImage
    \param   a_pathToInitialImage : path to the image to project
    \brief   Load the initial image for the analytical projection
    \return  0 if success, positive value otherwise
  */
  int PROJ_InitImage(const string& a_pathToInitialImage);
  /*!
    \fn      oImageSpace::PROJ_LoadInitialImage
    \param   a_pathToImage : path to the image to project
    \brief   Load the initial image for the analytical projection
    \return  0 if success, positive value otherwise
  */
  int PROJ_LoadInitialImage(const string& a_pathToImage);
  /*!
    \fn      oImageSpace::PROJ_SaveProjectionImage
    \brief   Save an image of the projected data (for analytic projection)
    \return  0 if success, positive value otherwise
    \todo    Support for PET projeted data
    #         Implement interfile file recovery in a scanner class function
  */
  int PROJ_SaveProjectionImage();


  // -------------------------------------------------------------------
  // Public Get and Set functions
  public:
    /*!
      \fn      oImageSpace::SetVerbose
      \param   a_verboseLevel
      \brief   set verbosity
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      oImageSpace::SetImageDimensionsAndQuantification
      \param   ap_ImageDimensionsAndQuantification
      \brief   set the pointer to the oImageDimensionsAndQuantification object
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ID = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      oImageSpace::GetNbBackwardImages
      \return  the number of backward images required for the selected optimizer
    */
    inline int GetNbBackwardImages()
           {return m_nbBackwardImages;}
    /*!
      \fn      oImageSpace::IsLoadedSensitivity
      \return  boolean indicating if the sensitivity is preloaded for the reconstruction (true) or computed on the fly (false)
    */
    inline bool IsLoadedSensitivity()
           {return m_loadedSensitivity;}
    /*!
      \fn      oImageSpace::IsLoadedMultiModal
      \return  boolean indicating if one or several multimodal images have been loaded (true) or not (false)    */
    inline bool IsLoadedMultiModal()
           {return m_loadedMultiModal;}
    /*!
      \fn      oImageSpace::IsLoadedMask
      \return  boolean indicating if a mask image has been loaded (true) or not (false)
    */
    inline bool IsLoadedMask()
           {return m_loadedMask;}
    /*!
      \fn      oImageSpace::Checked
      \brief   Simply check that the image dimensions and verbosity has been set
      \return  true if image dimensions and verbosity are set
    */
    inline bool Checked()
           {return mp_ID!=NULL && m_verbose!=-1;}
    /*!
      \fn      oImageSpace::GetNbMiscellaneousImages
      \return  The number of allocated miscellaneous images
    */
    inline int GetNbMiscellaneousImages()
           {return m_nbMiscellaneousImages;}


  // -------------------------------------------------------------------
  // Data members
  private:
    oImageDimensionsAndQuantification* mp_ID; /*!< Pointer to the oImageDimensionsAndQuantification object */
    int m_verbose;                            /*!< Verbosity */
    bool m_loadedSensitivity;                 /*!< Flag indicating if reconstruction is listmode, otherwise histogram because sensitivity is computed on-the-fly */
    bool m_loadedMultiModal;                  /*!< Flag indicating if multimodal images have been loaded */
    bool m_loadedMask;                        /*!< Flag indicating if a mask image has been loaded */
    int m_nbBackwardImages;                   /*!< Number of backward images in the related image matrix (for the optimizer) */
    int m_nbMiscellaneousImages;              /*!< Number of allocated miscellaneous images */
};

#endif
