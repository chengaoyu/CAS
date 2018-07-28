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
  \ingroup main_programs
  \brief  CASToR: Customizable and Advanced Software for Tomographic Reconstruction \n
          This is the main reconstruction program of the CASToR software.
          It reads/parses/checks the command-line options, initialize each manager classes with the correct set of options, and launch the reconstruction algorithm
          (currently managed by oIterativeAlgorithm).
          For list-mode datafiles, if a sensitivity image is not provided using "-sens", it is generated before reconstruction (managed by oSensitivityGenerator) .
*/


#include "gVariables.hh"
#include "gOptions.hh"
#include "oIterativeAlgorithm.hh"
#include "oSensitivityGenerator.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "iDataFilePET.hh"
#include "iDataFileSPECT.hh"
#include "iDataFileCT.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include "sRandomNumberGenerator.hh"
#include "iScannerPET.hh"
#include "sAddonManager.hh"
#include "sChronoManager.hh"

// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        H E L P     F U N C T I O N S
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================


/*!
  \fn      ShowHelp()
  \brief   Display main command line options for castor-recon
*/
void ShowHelp()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "Usage:  castor-recon.exe  -df file.cdh  -(f/d)out output  -it iter  [settings]" << endl;
  cout << endl;
  cout << "[Main options]:" << endl;
  cout << "  -df file.cdh         : Give an input CASTOR datafile header (no default)." << endl;
//MULTIFILE  cout << "  -df file.cdh         : Give an input CASTOR datafile header or several ones separated by commas (no default)." << endl;
//MULTIBED  cout << "                         Can use this option multiple times to specify multiple bed positions." << endl;
  cout << "  -fout name           : Give the root name for all output files (no default, alternative to -dout)" << endl;
  cout << "  -dout name           : Give the name of the output directory where all output files will be written (no default, alternative to -fout)" << endl;
  cout << "  -it  list            : Give the sequence of iterations:subsets separated by commas (no default)." << endl;
  cout << "  -dim x,y,z           : Give the number of voxels in each dimension (default: those of the scanner)." << endl;
  cout << "  -fov x,y,z           : Give the size of the field-of-view in each dimension, in mm (default: those of the scanner, or calculated from" << endl;
  cout << "                         the voxel sizes provided using '-vox')." << endl;
  cout << "  -vox x,y,z           : Give the size of the voxels in each dimension, in mm (default: those of the scanner, or calculated from the fov" << endl;
  cout << "                         if a fov value is provided using '-fov')." << endl;
  cout << "  -vb                  : Give the general verbosity level, from 0 (no verbose) to 5 and above (at the event level) (default: 1)." << endl;
  cout << endl;
  cout << "[Specific options]:" << endl;
  cout << "  -help-dim            : Print out specific help about dimensions settings." << endl; // managed by main
  cout << "  -help-in             : Print out specific help about input settings." << endl; // managed by main
  cout << "  -help-out            : Print out specific help about output settings." << endl; // managed by main
  cout << "  -help-algo           : Print out specific help about reconstruction optimizer algorithm settings." << endl; // managed by main
  cout << "  -help-proj           : Print out specific help about projection operators." << endl; // managed by main
  cout << "  -help-dynamic        : Print out specific help about dynamic methodologies settings." << endl; // managed by main
  cout << "  -help-imgp           : Print out specific help about image processing modules." << endl; // managed by main
  cout << "  -help-comp           : Print out specific help about computing settings." << endl; // managed by main
  cout << "  -help-corr           : Print out specific help about all corrections that can be disabled." << endl;
  cout << "  -help-misc           : Print out specific help about miscellaneous and verbose settings." << endl; // managed by main
  cout << endl;
  cout << "[Implemented Modules]:" << endl;
  cout << "  -help-scan           : Show the list of all scanners from the configuration directory." << endl; // managed by sScannerManager
  cout << "  -help-projm          : Show the list and description of all implemented projectors." << endl; // managed by oProjectorManager
  cout << "  -help-opti           : Show the list and description of all implemented optimizer algorithms." << endl; // managed by oOptimizerManager
  //cout << "  -help-pnlt           : Show the list and description of all implemented penalties for optimizers." << endl; // managed by oOptimizerManager
  //cout << "  -help-motion         : Show the list and description of all implemented image-based deformation models." << endl; // managed by oImageDeformationManager
  //cout << "  -help-dynm           : Show the list and description of all implemented dynamic models." << endl; // managed by oDynamicModelManager
  cout << "  -help-conv           : Show the list and description of all implemented image convolvers." << endl; // managed by oImageConvolverManager
  cout << "  -help-proc           : Show the list and description of all implemented image processing modules." << endl; // managed by oImageProcessingManager
  cout << endl;
  cout << endl;
  cout << "  --help,-h,-help      : Print out this help page." << endl; // managed by main
  cout << endl;
  #ifdef CASTOR_MPI
  cout << "  Compiled with MPI" << endl;
  #endif
  #ifdef CASTOR_OMP
  cout << "  Compiled with OpenMP" << endl;
  #endif
  #ifdef CASTOR_GPU
  cout << "  Compiled for GPU" << endl;
  #endif
  #if defined(CASTOR_OMP) || defined(CASTOR_MPI) || defined(CASTOR_GPU)
  cout << endl;
  #endif
  #ifdef BUILD_DATE
  cout << "  Build date: " << BUILD_DATE << endl;
  cout << endl;
  #endif
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpInput()
  \brief   Display command line options related to input settings for castor-recon
*/
void ShowHelpInput()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Input settings]:" << endl;
  cout << endl;
//MULTIBED  cout << "  -df-inv              : Invert the order of provided datafiles corresponding to the different bed positions." << endl;
//MULTIBED  cout << endl;
//MULTIBED  cout << "  -df file.Cdf         : Give an input CASTOR datafile (no default). Can use this option multiple times to specify multiple bed positions." << endl;
//MULTIFILE  cout << "  -df file.Cdf         : Give an input CASTOR datafile header, or several ones separated by commas (no default)." << endl;
//MULTIFILE  cout << "                         Several datafiles can be provided only in histogram mode, when these multiple datafiles represent different time" << endl;
//MULTIFILE  cout << "                         frames (identified by the time start and duration inside the header) or different respiratory or cardiac gates" << endl;
//MULTIFILE  cout << "                         (identified by the respiratory or cardiac gate index provided in the header)." << endl;
  cout << "  -df file.Cdf         : Give an input CASTOR datafile header (no default)." << endl;
  cout << endl;
  cout << "  -img file.hdr        : Give an input image as the initialization of the algorithm (default: uniform value)." << endl;
  cout << endl;
  cout << "  -sens file.hdr       : Provide the sensitivity image (default: sensitivity image is computed before reconstruction)." << endl;
  cout << "                         The image file should integrate all sensitivity images if more than one are required. If dual-gating is enabled and if it" << endl;
  cout << "                         requires sensitivity images for each gate, the image should integrate nb_resp_gates * nb_card_gates sensitivity images" << endl;
  cout << "                         (all cardiac-gated based sensitivity images for each one of the respiratory gates)." << endl;
  cout << endl;
  cout << "  -multimodal file.hdr : Provide additional images, from other modalities (anatomical, functional), or processed, for use in constrained recon " << endl;
  cout << "                         To provide several additional images for a single reconstructed image, the option should be used several times, new option for new image file " << endl;
  cout << endl;
  cout << "  -mask file.hdr       : Provide a mask image. The mask image is used in the projection, which requires the application of this same mask to the sensitivity, to avoid border issues. " << endl;
  cout << endl;
  cout << "  -norm file.cdh       : For list-mode data, provide a normalization data file for sensitivity computation (default: use the scanner LUT and" << endl;
  cout << "                         assume all LORs with a weight of 1.). This restricts the computation of the sensitivity to the LORs provided in the" << endl;
  cout << "                         normalization file and associated normalization factors and/or attenuation factors." << endl;
  cout << "                         For dynamic reconstructions with multiple frames, as many normalization files as frames can be supplied, their names" << endl;
  cout << "                         separated by commas. This is useful when dead-times correction is included in the normalization factors." << endl;
//MULTIBED  cout << "                         Can also use this option multiple times when multiple bed positions are reconstructed at once." << endl;
  cout << endl;
  cout << "  -atn file.hdr        : Give an input attenuation image (unit has to be cm-1) for sensitivity image generation or SPECT attenuation correction." << endl;
  cout << endl;
  cout << "  -help-in             : Print out this help." << endl;
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpOutput()
  \brief   Display command line options related to output settings for castor-recon
*/
void ShowHelpOutput()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Output settings]:" << endl;
  cout << endl;
  cout << "  -fout name           : Give the root name for all output files. All output files will be written as 'name_suffix.ext'." << endl;
  cout << "                         So the provided name should not end with '.' or '/' character. (no default, alternative to -dout)" << endl;
  cout << "  -dout name           : Give the name of the output directory where all output files will be written. All files will also" << endl;
  cout << "                         be prefixed by the name of the directory. The provided name should not end with '.' or '/' character." << endl;
  cout << "                         (no default, alternative to -fout)" << endl;
  cout << endl;
  cout << "  -oit list            : Give the sequence of output iterations as a list of 'a:b' pairs separated by commas. This will output one" << endl;
  cout << "                         iteration over 'a' until 'b' is reached, then it goes to the next pair of setting." << endl;
  cout << "                         Set '-1' to save only the last iteration. (default: save all iterations)" << endl;
  cout << endl;
  cout << "  -onbp prec            : By default, numbers are displayed using scientific format. This option allows to customize the format and precision" << endl;
  cout << "                        : The format is format,precision. f is the format (s=scientific, f=fixed), p is the precision" << endl;
  cout << "                          eg: -onbp f,5 --> fixed with 5 digits precision, -onbp -->  scientifici with max precision." << endl;
  cout << endl;
  cout << "  -omd                 : (M)erge (D)ynamic images. Indicate if a dynamic serie of 3D images should be written on disk in one file" << endl;
  cout << "                         instead of a serie of 3D images associated with an interfile metaheader" << endl;
  cout << endl;
/* TO BE IMPLEMENTED
  cout << "  -otb                 : Flag to say that we want to save the time basis images too (default: only the frames are saved)." << endl;
  cout << endl;
*/
  cout << "  -osens               : Flag to say that we want to save the sensitivity image of each subset/iteration, when in histogram mode." << endl;
  cout << endl;
  cout << "  -osub                : Flag to say that we want to save the image after each subset." << endl;
  cout << endl;
  cout << "  -olut                : If a scanner LUT (geometry information) is computed from a .geom file, it will be save on disk in the scanner repository" << endl;
  cout << endl;
  cout << "  -sens-histo          : If input file is a histogram, compute the global sensitivity from it (and still proceed to normal reconstruction after)" << endl;
  cout << endl;
  cout << "  -sens-only           : For list-mode data, exit directly after computing and saving the sensitivity." << endl;
  cout << endl;
  cout << "  -help-out            : Print out this help." << endl;
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpDimensions()
  \brief   Display command line options related to image dimensions for castor-recon
*/
void ShowHelpDimensions()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Dimensions options]:" << endl;
  cout << endl;
  cout << "  -dim x,y,z           : Give the number of voxels in each dimension (default: those of the scanner)." << endl;
  cout << endl;
  cout << "  -fov x,y,z           : Give the size of the field-of-view in each dimension, in mm (default: those of the scanner)." << endl;
  cout << endl;
  cout << "  -vox x,y,z           : Give the size of the voxels in each dimension, in mm (default: those of the scanner, or calculated from the fov if a fov value is provided using '-fov')." << endl;
  cout << endl;
  cout << "  -off x,y,z           : Give the offset of the field-of-view in each dimension, in mm (default: 0.,0.,0.)." << endl;
  cout << "                         (note this has no effect when using a pre-computed system matrix)" << endl;
  cout << endl;
  cout << "  -fov-out percent     : Give the percentage of the eliptical transaxial FOV to be kept while saving the image (default: no making)" << endl;
  cout << endl;
  cout << "  -slice-out value     : Give the number of axial slices to be masked at each border of the axial field-of-view (default: 0)" << endl;
  cout << endl;
  cout << "  -flip-out value      : Flip the image before saving it (not done in the computation); specify the axis (e.g. 'X', 'XY', 'YZ') (default: no flip)" << endl;
  cout << endl;
  cout << "  -help-dim            : Print out this help." << endl;
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpAlgo()
  \brief   Display command line options related to the Optimization algorithm module for castor-recon
*/
void ShowHelpAlgo()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "[Algorithm settings]:" << endl;
  cout << endl;
  cout << "  -opti param          : Give the algorithm to be used, along with a configuration file (algo:file.conf) or the list of parameters" << endl;
  cout << "                         associated to the algorithm (algo,param1,param2,...). If the algorithm only is specified, the default" << endl;
  cout << "                         configuration file is used if any. By default, the MLEM algorithm is used. For specific help, use option -help-opti." << endl;
  cout << endl;
  cout << "  -opti-fom            : Flag to say that we want to compute and print figures-of-merit (likelihood and RMSE) in the data-space." << endl;
  cout << endl;
  cout << "  -opti-stat           : Flag to say that we want to compute and print basic statistics about the image update." << endl;
  cout << endl;
/* NOT YET IMPLEMENTED
  cout << "  -pnlt param          : Give the penalty to be used with the algorithm (if the latter allows to do so), along with a configuration file (penalty:file.conf)" << endl;
  cout << "                         or the list of parameters associated to the penalty (penalty,param1,param2,...). If the penalty only is specified, the default" << endl;
  cout << "                         configuration file is used if any. By default, no penalty is used. For specific help, use option -help-penalty." << endl;
  cout << endl;
*/
  cout << "  -help-opti           : Print out specific help about algorithm settings." << endl; // managed by oOptimizerManager
  cout << endl;
/* NOT YET IMPLEMENTED
  cout << "  -help-pnlt           : Print out specific help about penalty settings." << endl; // managed by oOptimizerManager
  cout << endl;
*/
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpProj()
  \brief   Display command line options related to the Projector module for castor-recon
*/
void ShowHelpProj()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "[Projector settings]:" << endl;
  cout << endl;
  cout << "  -proj  param         : Give the projector to be used for both forward and backward projections, along with a configuration file (proj:file.conf)" << endl;
  cout << "                         or the list of parameters associated to the projector (proj,param1,param2,...). If the projector only is specified, the" << endl;
  cout << "                         default configuration file is used. By default, the Siddon projector is used. For specific help, use option -help-proj." << endl;
  cout << endl;
  cout << "  -projF param         : Give the projector to be used for forward projections. See option -proj for details." << endl;
  cout << endl;
  cout << "  -projB param         : Give the projector to be used for backward projections. See option -proj for details." << endl;
  cout << endl;
  cout << "  -proj-common         : Give parameterization of different common projector-related options." << endl;
  cout << endl;
/* TO BE IMPLEMENTED
  cout << "  -ignore-POI          : Flag to say that we want to ignore any potential POI information (default: use it if present in the datafile)." << endl;
  cout << endl;
  * */
  cout << "  -ignore-TOF          : Flag to say that we want to ignore any potential TOF information (default: use it if present in the datafile)." << endl;
  cout << endl;
  cout << "  -help-projm           : Print out specific help about projector settings." << endl; // managed by oProjectorManager
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpImgp()
  \brief   Display command line options related to the Image Processing module for castor-recon
*/
void ShowHelpImgp()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "[Image processing settings]:" << endl;
  cout << endl;
  cout << "  -conv  param;when    : Give an image convolver model to be used within the algorithm, along with a configuration file (conv:file.conf::when) or the" << endl;
  cout << "                         list of parameters associated to the convolver (conv,param1,param2,...::when). If the convolver only is specified, its default" << endl;
  cout << "                         configuration file is used. By default, no convolver is applied. Multiple convolvers can be combined simply by repeating this" << endl;
  cout << "                         this option. The mandatory 'when' parameter specifies when the convolver is applied (psf, sieve, forward, backward, post, intra)." << endl;
  cout << "                         For more specific help, use option -help-conv." << endl;
  cout << endl;
  cout << "  -help-conv           : Print out specific help about the image convolver settings." << endl; // managed by oImageConvolverManager
  cout << endl;
  cout << "  -proc  param;when    : Give an image processing module to be used within the algorithm, along with a configuration file (proc:file.conf;when) or the" << endl;
  cout << "                         list of parameters associated to the module (proc,param1,param2,...;when). If the module only is specified, its default" << endl;
  cout << "                         configuration file is used. By default, no image processing module is applied. Multiple modules can be combined simply by" << endl;
  cout << "                         repeating this this option. The mandatory 'when' parameter specifies when the module is applied (forward, backward, post, intra)." << endl;
  cout << "                         For more specific help, use option -help-proc." << endl;
  cout << endl;
  cout << "  -help-proc           : Print out specific help about the image processing module settings." << endl; // managed by oImageProcessingManager
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpDynamic()
  \brief   Display command line options related to the dynamic features for castor-recon
*/
void ShowHelpDynamic()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Dynamic settings]:" << endl;
  cout << endl;
  cout << "  -frm list            : Give the framing for the reconstruction where 'list' is a list of all frames durations in seconds separated" << endl;
  cout << "                         by commas. To specify a dead frame, add 'df' concatenated to the frame duration. Milliseconds maximum precision." << endl;
  cout << "                         (default: 1 frame of the whole input file duration)" << endl;
  cout << endl;
  cout << "  -g      path_to_file : Provide text file for dynamic data splitting associated to respiratory/cardiac gating or involuntary patient motion correction" << endl;
  cout << endl;
/* TO BE VALIDATED
  cout << "  -time-basis nb:file  : Give the number of time basis functions related to the frames given by option -frm." << endl;
  cout << "                         The file should contain one line for each time-basis function, and as many numbers as frames." << endl;
  cout << "  -resp-basis nb:file  : Give the number of intrinsic respiratory basis functions related to the respiratory gates" << endl;
  cout << "                         The file should contain one line for each respiratory basis function, and as many numbers as respiratory gates." << endl;
  cout << "  -card-basis nb:file  : Give the number of intrinsic cardiac basis functions related to the cardiac gates" << endl;
  cout << "                         The file should contain one line for each cardiac basis function, and as many numbers as cardiac gates." << endl;
  cout << endl;  
  cout << " -dynamic-model param  : Dynamic model applied to either the frames of a dynamic acquisition, respiratory-gated frames, cardiac-gated frames, or simultaneously between these datasets." << endl;
  cout << "                         Select the dynamic model to be used, along with a configuration file (model:file) or the list of parameters associated to the model (model_name,param1,param2,...)." << endl;  
  cout << endl; 
  cout << "  -rm param            : Provide an image-based deformation model to be used for respiratory motion correction, along with a configuration file (deformation:file)" << endl;
  cout << "                         or the list of parameters associated to the projector (deformation,param1,param2,...)." << endl;
  cout << endl; 
  cout << "  -cm param            : Provide an image-based deformation model to be used for cardiac motion correction, along with a configuration file (deformation:file)" << endl;
  cout << "                         or the list of parameters associated to the projector (deformation,param1,param2,...)." << endl;
  cout << endl; 
  cout << "  -rcm param           : Provide an image-based deformation model to be used for both respiratory and cardiac motion corrections, along with a configuration file (deformation:file)" << endl;
  cout << "                         or the list of parameters associated to the projector (deformation,param1,param2,...)." << endl;
  cout << endl; 
  cout << "  -im param            : Provide an image-based deformation model to be used for involuntary motion correction, along with a configuration file (deformation:file)" << endl;
  cout << "                         or the list of parameters associated to the projector (deformation,param1,param2,...)." << endl;
  cout << endl;
  cout << "  -qdyn file           : Provide a text file containing quantitative factors specific to dynamic frames, respiratory or cardiac gates." << endl;
  cout << "                         The file should provide factors with the keywords 'FRAME_QUANTIFICATION_FACTORS' and 'GATE_QUANTIFICATION_FACTORS' " << endl;
  cout << "                         The number of factors must be consistent with the number of frames/gates " << endl;
  cout << "                         If the data contains several frames and gates, the gate quantification factors should be entered on a separate line for each frame " << endl;
  cout << endl;
  cout << "  -help-dynm           : Print out specific help about dynamic model." << endl; // managed by oDynamicModelManager
  cout << endl;
  cout << "  -help-motion         : Print out specific help about deformation in image space for motion correction." << endl; // managed by oImageDeformationManager
  cout << endl;
  */
  cout << "  -help-dynamic        : Print out this help." << endl;
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpComputation()
  \brief   Display command line options related to the computation settings for castor-recon
*/
void ShowHelpComputation()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Computation settings]:" << endl;
  cout << endl;
  #ifdef CASTOR_OMP
  cout << "  -th param            : Set the number of threads for parallel computation (default: 1). If 0 is given, then the maximum number of available threads is automatically determined." << endl;
  cout << "                         Can also give two parameters separated by a comma (e.g. 16,4), to distinguish between the number of threads for projection and image operations respectively." << endl;
  cout << endl;
  #endif
  #ifdef CASTOR_GPU
  cout << "  -gpu                 : Flag to say that we want to use the GPU device (default: use the CPU only)." << endl;
  cout << endl;
  #endif
  cout << "  -load size           : Give a number between 0 and 100 to specify the number of events which will be loaded in the RAM." << endl;
  cout << "                         It will be interpreted as a percent ratio  (default: 100 means all data buffered. 0 means no data buffered.)" << endl;
  cout << endl;
  cout << "  -proj-comp           : Give the strategy for projection line computation. Here are the three different strategies that can be used:" << endl;
  cout << "                     1 : Image-based system matrix elements storage: The voxels weights are added in a matrix representing the whole image, so" << endl;
  cout << "                         the addition of a new line to the previous ones is straightforward only by adding the weights to the corresponding voxels." << endl;
  cout << "                         As it is insanely long, it can possibly be used for example with extremely complex projectors that makes use of huge number" << endl;
  cout << "                         of ray tracings for a single event, where the list of contributions can become longer than the number of voxels in the image." << endl;
  cout << "                         This strategy is not compatible with SPECT reconstruction including attenuation correction." << endl;
  cout << "                     2 : Fixed-size list storage of system matrix elements: The voxels are added one by one in two separated lists, one containing voxel" << endl;
  cout << "                         indices and the other voxel weights. When a voxel is added to the oProjectionLine, it is simply pilled-up to the list. The list" << endl;
  cout << "                         has a fixed size which is provided by the EstimateMaxNumberOfVoxelsPerLine() function from the vProjector class. There are no" << endl;
  cout << "                         ckecks at all for possible buffer overflows. This is the fastest strategy and default one." << endl;
  cout << "                     3 : Adaptative-size list storage of system matrix elements: This is the same as the fixed-size strategy except that the size can be" << endl;
  cout << "                         upgraded if the current number of contributing voxels exceed the list's size. The first allocated size corresponds to the diagonal" << endl;
  cout << "                         of the image. During the first iteration, this size will be upgraded until it will reach a size suitable for all events. Thus it" << endl;
  cout << "                         is a bit slower than the fixed-list strategy during the first iteration, but is optimal with respect to RAM usage." << endl;
  cout << endl;
  cout << "  -rng-seed            : Give a seed for the random number generator (should be >=0)" << endl;
  cout << endl;
  cout << "  -rng-extra           : Give the number of additional random number generators, if needed (in addition to one rng per thread) " << endl;
  cout << endl;
  cout << "  -help-comp           : Print out this help." << endl;
  cout << endl;
  #ifdef CASTOR_MPI
  cout << "  Compiled with MPI" << endl;
  #endif
  #ifdef CASTOR_OMP
  cout << "  Compiled with OpenMP" << endl;
  #endif
  #ifdef CASTOR_GPU
  cout << "  Compiled for GPU" << endl;
  #endif
  #if defined(CASTOR_OMP) || defined(CASTOR_MPI) || defined(CASTOR_GPU)
  cout << endl;
  #endif
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpMiscellaneous()
  \brief   Display command line options related to miscellaneous settings for castor-recon
*/
void ShowHelpMiscellaneous()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Miscellaneous settings]:" << endl;
  cout << endl;
  cout << "  -vb                  : Give the general verbosity level, from 0 (no verbose) to 5 and above (at the event level) (default: 1)." << endl;
  cout << "  -vb-algo             : Give the verbose level specific to the algorithm (default: same as general verbose level)." << endl;
  cout << "  -vb-opti             : Give the verbose level specific to the optimizer (default: same as general verbose level)." << endl;
  cout << "  -vb-proj             : Give the verbose level specific to the projector (default: same as general verbose level)." << endl;
  cout << "  -vb-conv             : Give the verbose level specific to the image convolver (default: same as general verbose level)." << endl;
  cout << "  -vb-proc             : Give the verbose level specific to the image processing (default: same as general verbose level)." << endl;
  cout << "  -vb-scan             : Give the verbose level specific to the scanner (default: same as general verbose level)." << endl;
  cout << "  -vb-data             : Give the verbose level specific to the data and image management (default: same as general verbose level)." << endl;
  cout << "  -vb-defo             : Give the verbose level specific to the deformation (default: same as general verbose level)." << endl;
  cout << "  -vb-dyna             : Give the verbose level specific to the dynamic model (default: same as general verbose level)." << endl;
  cout << "  -vb-sens             : Give the verbose level specific to the sensitivity computation (default: same as general verbose level)." << endl;
  cout << endl;
  cout << "  -conf                : Give the path to the CASToR config directory (default: located through the CASTOR_CONFIG environment variable)." << endl;
  cout << endl;
  cout << "  -help-misc           : Print out this help." << endl;
  cout << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpCorrection()
  \brief   Display command line options related to correction settings for castor-recon
*/
void ShowHelpCorrection()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Correction settings]:" << endl;
  cout << endl;
  cout << "  -ignore-corr list    : Give the list of corrections that should be ignored, separated by commas (default: all corrections applied if present)." << endl;
  cout << "                         Here is a list of all corrections that can be ignored:" << endl;
  cout << "                         - attn: to ignore the attenuation correction (emission only)" << endl;
  cout << "                         - norm: to ignore the normalization correction (emission only)" << endl;
  cout << "                         - rand: to ignore the random correction (PET only)" << endl;
  cout << "                         - scat: to ignore the scatter correction" << endl;
  cout << "                         - deca: to ignore the decay correction (emission only)" << endl;
  cout << "                         - brat: to ignore the branching ratio correction (emission only)" << endl;
  cout << "                         - fdur: to ignore the frame duration correction (emission only)" << endl;
  cout << "                         - cali: to ignore the calibration correction (emission only)" << endl;
  cout << endl;
}

// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        M A I N     P R O G R A M
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================

int main(int argc, char** argv)
{
  // ============================================================================================================
  // MPI stuff
  // ============================================================================================================
  int mpi_rank = 0;
  int mpi_size = 1;
  #ifdef CASTOR_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  #endif

  // No argument, then show help
  if (argc==1)
  {
    ShowHelp();
    Exit(EXIT_SUCCESS);
  }

  // ============================================================================================================
  // Parameterized variables with their default values
  // ============================================================================================================

  // --------------------------------------------------------------------------------
  // Dimensions settings
  // --------------------------------------------------------------------------------

  // Initialization of the voxel and field-of-view values in each spatial dimensions. default values depend of the scanner.
  INTNB nb_voxX=-1, nb_voxY=-1, nb_voxZ=-1;
  FLTNB fov_sizeX=-1., fov_sizeY=-1., fov_sizeZ=-1.;
  FLTNB vox_sizeX=-1., vox_sizeY=-1., vox_sizeZ=-1.;
  // Initialization offset values of the field-of-view in each spatial dimensions
  FLTNB offsetX = 0., offsetY = 0., offsetZ = 0.;
  FLTNB fov_out = 0.;
  INTNB slice_out = 0;
  // Output flip
  string flip_out = "";

  // --------------------------------------------------------------------------------
  // Dynamic settings
  // --------------------------------------------------------------------------------

  // Frame list descriptor
  string frame_list = "";
  // Number of respiratory gates in the data. default : 1
  int nb_resp_gates = 1;
  // Number of cardiac gates in the data. default : 1
  int nb_card_gates = 1;
  // Path to file containing the respiratory/cardiac gate splitting of the data
  string path_to_4D_data_splitting_file = "";
  // String gathering the name of the used dynamic model applied to the frames/resp gates/card gates of the dynamic acquisition (default : none), corresponding file, and corresponding string gathering parameters 
  string dynamic_model_options = "";
  // String gathering the name of the used deformation model for respiratory motion (default : none), corresponding file, and corresponding string gathering parameters 
  string resp_motion_options = "";
  // String gathering the name of the used deformation model for cardiac motion (default : none), corresponding file, and corresponding string gathering parameters  
  string card_motion_options = "";
  // String gathering the name of the used deformation model for both respiratory and cardiac motions (default : none), corresponding file, and corresponding string gathering parameters  
  string double_motion_options = "";
  // String gathering the name of the used deformation model for involuntary motion (default : none), corresponding file, and corresponding string gathering parameters  
  string ipat_motion_options = "";
  // Number of time basis functions and associated file
  int nb_time_basis = 1;
  string path_to_time_basis_coef = "";
  // Number of respiratory basis functions and associated file
  int nb_resp_basis = 1;
  string path_to_resp_basis_coef = "";
  // Number of cardiac basis functions and associated file
  int nb_card_basis = 1;
  string path_to_card_basis_coef = "";
  // Number of cardiac basis functions and associated file
  string path_to_dynamic_quantification_file= "";

  // --------------------------------------------------------------------------------
  // Input settings
  // --------------------------------------------------------------------------------

  // Number of bed positions
  int nb_beds = 0;
  // Vector (for bed positions) containing string which gather the path to the data filenames provided by the user. no default 
  vector<string> path_to_data_filename;
  // Invert bed datafile names order (bed positions)
  bool invert_datafile_bed_order_flag = false;
  // Path to an image used as initialization. no default (uniform image is used in this case)
  string path_to_initial_img = "";
  // Path to an image used as initialization for the sensitivity image. no default (sensitivity image is computed in this case)
  string path_to_sensitivity_img;
  // Path to a normalization data file for sensitivity computation with list-mode data
  vector<string> path_to_normalization_filename;
  // Path to an image used for the attenuation. default : uniform image
  string path_to_attenuation_img;
  // Paths to additional images
  vector<string> path_to_multimodal_img;
  // Path to an image used as mask
  string path_to_mask_img = "";
  // Number of resp and card atn images for sensitivity image generation
  int nb_atn_resp_imgs = 1;
  int nb_atn_card_imgs = 1;
  // Ignored corrections (by default, if present, there are applied)
  string ignored_corrections = "";

  // --------------------------------------------------------------------------------
  // Output settings
  // --------------------------------------------------------------------------------

  // Output directory name.
  string path_dout = "";
  // Or root name
  string path_fout = "";
  // Iterations to be saved
  string output_iterations = "";
  // Merge output dynamic images on disk into one file or not
  bool merge_dynamic_imgs_flag = false;
  // Write scanner LUT generated by a geom file on disk
  bool save_LUT_flag = false;
  // Save the sensitivity image in histogram mode for each subset/iteration
  bool save_sens_histo = false;
  // Save the image after each subset
  bool save_subset_image = false;
  // Save the time basis functions (flag). default : no saving
  //bool save_time_basis_functions_flag = false; //TODO
  // Exit directly after computing the sensitivity
  bool exit_after_sensitivity = false;
  // Compute sensitivity with a histogram file as input
  bool sensitivity_from_histogram = false;
  // Precision for output number display
  string onb_prec = "s,0";
  
  // --------------------------------------------------------------------------------
  // Projector settings
  // --------------------------------------------------------------------------------

  // String gathering the name of the projector for forward/backward projection with specific options (default: Native Siddon for both) 
  string options_projector  = "incrementalSiddon";
  string options_projectorF = "incrementalSiddon";
  string options_projectorB = "incrementalSiddon";
  string options_projector_common = "";
  // Default projector computation strategy
  int projector_computation_strategy = FIXED_LIST_COMPUTATION_STRATEGY;
  // POI & TOF flags for the use of such information in the reconstruction 
  bool ignore_POI = false;
  bool ignore_TOF = false;

  // --------------------------------------------------------------------------------
  // Algorithm settings
  // --------------------------------------------------------------------------------

  // Numbers of iterations and associated numbers of subsets (default : 1 for both) 
  string nb_iterations_subsets = "";
  // String gathering the name of the optimizer with specific options (default: MLEM)
  string options_optimizer = "MLEM";
  // Boolean to say if we want to compute FOM in the data-space
  bool optimizer_fom = false;
  // Boolean to say if we want to compute image update basic statistics
  bool optimizer_stat = false;
  // String gathering the name of the penalty with specific options (default: no penalty)
  string options_penalty = "";

  // --------------------------------------------------------------------------------
  // Image convolvers and processing modules
  // --------------------------------------------------------------------------------

  // String vector gathering the name of the image convolvers with specific options (default: no image convolver)
  vector<string> options_image_convolver;
  // String vector gathering the name of the image processing modules with specific options (default: no image processing module)
  vector<string> options_image_processing;

  // --------------------------------------------------------------------------------
  // Computation settings
  // --------------------------------------------------------------------------------

  // Using GPU (flag) ->NOTE : default : only CPU
  bool gpu_flag = false;
  // Number of threads
  string nb_threads = "1";

  // --------------------------------------------------------------------------------
  // Miscellaneous settings
  // --------------------------------------------------------------------------------

  // General verbose level
  int verbose_general = 1;
  // Specific verbose levels
  int verbose_algo = -1;
  int verbose_opti = -1;
  int verbose_proj = -1;
  int verbose_conv = -1;
  int verbose_proc = -1;
  int verbose_scan = -1;
  int verbose_data = -1;
  int verbose_defo = -1;
  int verbose_dyna = -1;
  int verbose_sens = -1;
  // Path to config directory
  string path_to_config_dir = "";
  // Initial seed for random number generator
  int64_t random_generator_seed = -1;
  int nb_extra_random_generators = 0;

  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================

  // TODO : replace remaining atoi, atof, etc, by ConvertFromString()

  // Must manually increment the option index when an argument is needed after an option
  for (int i=1; i<argc; i++)
  {
    // Get the option as a string
    string option = (string)argv[i];

    // --------------------------------------------------------------------------------
    // Miscellaneous settings
    // --------------------------------------------------------------------------------

    // Show help
    if (option=="-h" || option=="--help" || option=="-help")
    {
      ShowHelp();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for integrated scanners (managed by scanner manager)
    else if (option=="-help-scan")
    {
      if(sScannerManager::GetInstance()->ShowScannersDescription())
        Cerr("***** castor-recon() -> An error occured when trying to output the available scanners from the scanner repository !'" << endl;);
      Exit(EXIT_SUCCESS);
    }
    // Specific help for optimizer settings (managed by optimizer children)
    else if (option=="-help-opti")
    {
      sAddonManager::GetInstance()->ShowHelpOptimizer();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for penalty settings (managed by penalty children)
    else if (option=="-help-pnlt")
    {
      sAddonManager::GetInstance()->ShowHelpPenalty();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for projector settings (managed by projector children)
    else if (option=="-help-projm")
    {
      // Call specific showHelp function from vProjector children
      sAddonManager::GetInstance()->ShowHelpProjector();
      // Call the static showHelp function from vProjector
      vProjector::ShowCommonHelp();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for image convolver settings (managed by image convolver children)
    else if (option=="-help-conv")
    {
      // Call specific showHelp function from vImageConvolver children
      sAddonManager::GetInstance()->ShowHelpImageConvolver();
      // Call the static showHelp function from oImageConvolverManager
      oImageConvolverManager::ShowCommonHelp();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for image processing settings (managed by image processing children)
    else if (option=="-help-proc")
    {
      // Call specific showHelp function from vImageProcessingModule children
      sAddonManager::GetInstance()->ShowHelpImageProcessingModule();
      // Call the static showHelp function from oImageProcessingManager
      oImageProcessingManager::ShowCommonHelp();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for dynamic model (managed by dynamic model children)
    else if (option=="-help-dynm")
    {
      sAddonManager::GetInstance()->ShowHelpDynamicModel();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for image deformation settings (managed by projector children)
    else if (option=="-help-motion")
    {
      sAddonManager::GetInstance()->ShowHelpDeformation();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for input settings (managed by main)
    else if (option=="-help-in")
    {
      ShowHelpInput();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for output settings (managed by main)
    else if (option=="-help-out")
    {
      ShowHelpOutput();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for dimensions settings (managed by main)
    else if (option=="-help-dim")
    {
      ShowHelpDimensions();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for optimizer settings (managed by main)
    else if (option=="-help-algo")
    {
      ShowHelpAlgo();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for projector settings (managed by main)
    else if (option=="-help-proj")
    {
      ShowHelpProj();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for image processing settings (managed by main)
    else if (option=="-help-imgp")
    {
      ShowHelpImgp();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for computation settings (managed by main)
    else if (option=="-help-comp")
    {
      ShowHelpComputation();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for miscellaneous settings (managed by main)
    else if (option=="-help-misc")
    {
      ShowHelpMiscellaneous();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for dynamic settings (managed by main)
    else if (option=="-help-dynamic")
    {
      ShowHelpDynamic();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for correction settings (managed by main)
    else if (option=="-help-corr")
    {
      ShowHelpCorrection();
      Exit(EXIT_SUCCESS);
    }
    // General verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_general))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_general << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Algorithm verbosity level
    else if (option=="-vb-algo")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_algo))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_algo << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Optimizer verbosity level
    else if (option=="-vb-opti")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_opti))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_opti << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Projector verbosity level
    else if (option=="-vb-proj")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_proj))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_proj << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Image convolver verbosity level
    else if (option=="-vb-conv")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_conv))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_conv << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Image processing verbosity level
    else if (option=="-vb-proc")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_proc))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_proc << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Scanner verbosity level
    else if (option=="-vb-scan")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_scan))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_scan << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Data and image verbosity level
    else if (option=="-vb-data")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_data))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_data << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Deformation verbosity level
    else if (option=="-vb-defo")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_defo))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_defo << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Dynamic verbosity level
    else if (option=="-vb-dyna")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_dyna))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_dyna << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Sensitivity verbosity level
    else if (option=="-vb-sens")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_sens))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided verbosity level '" << verbose_sens << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // RNG seed
    else if (option=="-rng-seed")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if(ConvertFromString(argv[i+1], &random_generator_seed))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided number '" << random_generator_seed << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // number of extra RNG
    else if (option=="-rng-extra")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if(ConvertFromString(argv[i+1], &nb_extra_random_generators))
      {
        Cerr("***** castor-recon() -> Exception when trying to read provided number '" << nb_extra_random_generators << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Path to config directory
    else if (option=="-conf")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_config_dir = (string)argv[i+1];
      i++;
    }

    // --------------------------------------------------------------------------------
    // Dimensions settings
    // --------------------------------------------------------------------------------

    // Dimensions: number of voxels
    else if (option=="-dim")
    { 
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      INTNB input[3];
      if (ReadStringOption(argv[i+1], input, 3, ",", option))
      {
        Cerr("***** castor-recon() -> Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
        Exit(EXIT_FAILURE);
      }
      nb_voxX = input[0];
      nb_voxY = input[1];
      nb_voxZ = input[2];
      i++;
    }
    // Dimensions: size of the field-of-view in mm
    else if (option=="-fov")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      FLTNB input[3];
      if (ReadStringOption(argv[i+1], input, 3, ",", option))
      {
        Cerr("***** castor-recon() -> Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
        Exit(EXIT_FAILURE);
      }
      fov_sizeX = input[0];
      fov_sizeY = input[1];
      fov_sizeZ = input[2];
      i++;
    }
    // Dimensions: size of the voxels in mm
    else if (option=="-vox")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      FLTNB input[3];
      if (ReadStringOption(argv[i+1], input, 3, ",", option))
      {
        Cerr("***** castor-recon() -> Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
        Exit(EXIT_FAILURE);
      }
      vox_sizeX = input[0];
      vox_sizeY = input[1];
      vox_sizeZ = input[2];
      i++;
    }
    // Dimensions: offset of the image in mm
    else if (option=="-off")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      FLTNB input[3];
      if (ReadStringOption(argv[i+1], input, 3, ",", option))
      {
        Cerr("***** castor-recon() -> Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
        Exit(EXIT_FAILURE);
      }
      offsetX = input[0];
      offsetY = input[1];
      offsetZ = input[2];
      i++;
    }
    // Size of the output transaxial FOV in percentage
    else if (option=="-fov-out")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      fov_out = atof(argv[i+1]);
      i++;
    }
    // Number of slices to be masked at output
    else if (option=="-slice-out")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      slice_out = ((INTNB)(atoi(argv[i+1])));
      i++;
    }
    // Output flip
    else if (option=="-flip-out")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      flip_out = (string)(argv[i+1]);
      i++;
    }

    // --------------------------------------------------------------------------------
    // Dynamic settings
    // --------------------------------------------------------------------------------

    // Dimensions: number of frames
    else if (option=="-frm") // TODO: more checks to be done in the option reading
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      frame_list = (string)argv[i+1];
      i++;
    }    
    // Dimensions: number of respiratory gates
    else if (option=="-g")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      // Path to the dynamic file
      path_to_4D_data_splitting_file = ((string)argv[i+1]);
      // Recover number of respiratory gates from file
      if (ReadDataASCIIFile(path_to_4D_data_splitting_file, "nb_respiratory_gates", &nb_resp_gates, 1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
      {
        Cerr("***** castor-recon() -> Error when trying to read the number of respiratory gates in the file " << path_to_4D_data_splitting_file << " for option " << option << endl);
        Exit(EXIT_FAILURE);
      }
      // Recover number of cardiac gates from file
      if (ReadDataASCIIFile(path_to_4D_data_splitting_file, "nb_cardiac_gates", &nb_card_gates, 1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
      {
        Cerr("***** castor-recon() -> Error when trying to read the number of cardiac gates in the file " << path_to_4D_data_splitting_file << " for option " << option << endl);
        Exit(EXIT_FAILURE);
      }
      // Check for wrong initialization
      if (nb_resp_gates<1 || nb_card_gates <1)
      {
        Cerr("***** castor-recon() -> Incorrect initialization of the number of gates for the option: " << option << ". This number should be >= 1" << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Time basis functions
    else if (option=="-time-basis")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string input = argv[i+1];
      size_t column_pos = input.find_first_of(":");
      if (column_pos==string::npos)
      {
        Cerr("***** castor-recon() -> Incorrect argument after option " << option << ", ':' sign is missing !" << endl);
        Exit(EXIT_FAILURE);
      }
      string str = input.substr(0, column_pos);
      nb_time_basis = atoi(str.c_str());
      path_to_time_basis_coef = input.substr(column_pos+1);
      i++;
    }
    // Respiratory basis functions
    else if (option=="-resp-basis")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string input = argv[i+1];
      size_t column_pos = input.find_first_of(":");
      if (column_pos==string::npos)
      {
        Cerr("***** castor-recon() -> Incorrect argument after option " << option << ", ':' sign is missing !" << endl);
        Exit(EXIT_FAILURE);
      }
      string str = input.substr(0, column_pos);
      nb_resp_basis = atoi(str.c_str());
      path_to_resp_basis_coef = input.substr(column_pos+1);
      i++;
    }
    // Cardiac basis functions
    else if (option=="-card-basis")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string input = argv[i+1];
      size_t column_pos = input.find_first_of(":");
      if (column_pos==string::npos)
      {
        Cerr("***** castor-recon() -> Incorrect argument after option " << option << ", ':' sign is missing !" << endl);
        Exit(EXIT_FAILURE);
      }
      string str = input.substr(0, column_pos);
      nb_card_basis = atoi(str.c_str());
      path_to_card_basis_coef = input.substr(column_pos+1);
      i++;
    }
    // Dynamic model applied to the frames/respiratory gates/cardiac gates of the dynamic acquisition
    else if (option=="-dynamic-model") 
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      dynamic_model_options = (string)argv[i+1];
      i++;
    }
    // Data for respiratory motion correction based on deformation
    else if (option=="-rm")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      resp_motion_options = (string)argv[i+1];
      i++;
    }
    // Data for cardiac motion correction based on deformation
    else if (option=="-cm")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      card_motion_options =  (string)argv[i+1];
      i++;
    }
    // Data for respiratory and cardiac motion corrections based on deformation
    else if (option=="-rcm")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      double_motion_options =  (string)argv[i+1];
      i++;
    }
    // Data for involuntary motion correction based on deformation
    else if (option=="-im")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      ipat_motion_options = (string)argv[i+1];
      i++;
    }
    // Data for involuntary motion correction based on deformation
    else if (option=="-qdyn")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_dynamic_quantification_file = (string)argv[i+1];
      i++;
    }

    // --------------------------------------------------------------------------------
    // Input settings
    // --------------------------------------------------------------------------------

    // DataFiles
    else if (option=="-df") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string file_name = (string)argv[i+1];
      path_to_data_filename.push_back(file_name);
      nb_beds++;
      i++;
    }
    // Invert datafile order
    else if (option=="-df-inv")
    {
      invert_datafile_bed_order_flag = true;
    }
    // Image for the initialisation of the algorithm : What should we do in case of multiple beds ? or multiple frames ? TODO
    else if (option=="-img")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_initial_img = argv[i+1];
      i++;
    }
    // Sensitivity image
    else if (option=="-sens")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
              
      path_to_sensitivity_img = argv[i+1];
      i++;
    }
    // Normalization data files
    else if (option=="-norm")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string file_name = (string)argv[i+1];
      path_to_normalization_filename.push_back(file_name);
      i++;
    }

    // Image for the attenuation
    else if (option=="-atn")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      
      path_to_attenuation_img = (string)argv[i+1];
      
      if (path_to_attenuation_img != "") // parse
      {
        // Retrieve nb of gated images from header if required
        Intf_fields IF;
        IntfKeyInitFields(&IF);
        if(IntfReadHeader(path_to_attenuation_img, &IF, 0) )
        {
          Cerr("***** castor-recon() -> An error occurred while trying to read the interfile header of attenuation file " << path_to_attenuation_img << " !" << endl);  
          Exit(EXIT_FAILURE);
        }
        nb_atn_resp_imgs = IF.nb_resp_gates;
        nb_atn_card_imgs = IF.nb_card_gates;
      }
      i++;
    }
    // Multimodal images
    else if (option=="-multimodal")
    {
      if (i>=argc-1)
      {
        Cerr("***** Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string file_name = (string)argv[i+1];
      path_to_multimodal_img.push_back(file_name);
      i++;
    }
    // Mask image
    else if (option=="-mask")
    {
      if (i>=argc-1)
      {
        Cerr("***** Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }

      path_to_mask_img = argv[i+1];
      i++;
    }

    // --------------------------------------------------------------------------------
    // Output settings
    // --------------------------------------------------------------------------------

    // Name of the output directory
    else if (option=="-dout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_dout = argv[i+1];
      i++;
    }
    // Base name of the output files
    else if (option=="-fout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_fout = argv[i+1];
      i++;
    }
    // List of iterations to be outputed
    else if (option=="-oit")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      output_iterations = (string)argv[i+1];
      i++;
    }
    // Output number precision
    else if (option=="-onbp")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      onb_prec = argv[i+1];
      i++;
    }
    // Flag to say that we want to save time basis functions too
    else if (option=="-omd")
    {
      merge_dynamic_imgs_flag = true;
    }
    // Flag to say that we want to save time basis functions too
    else if (option=="-olut")
    {
      save_LUT_flag = true;
    }
    // Flag to say that we want to save the sensitivity image for each subset/iteration in histogram mode
    else if (option=="-osens")
    {
      save_sens_histo = true;
    }
    // Flag to say that we want to save the image after each subset
    else if (option=="-osub")
    {
      save_subset_image = true;
    }
    // Flag to say that we want to save time basis functions too
    else if (option=="-otb")
    {
      Cerr("!!!!! castor-recon() -> Warning: option -otb not yet implemented !" << endl);
    }
    // Flag to say that we exit directly after computing and saving the sensitivity
    else if (option=="-sens-only")
    {
      exit_after_sensitivity = true;
    }
    // Flag to say that we still want to compute global sensitivity from the input histogram file
    else if (option=="-sens-histo")
    {
      sensitivity_from_histogram = true;
    }

    // --------------------------------------------------------------------------------
    // Algorithm settings
    // --------------------------------------------------------------------------------

    // List of iterations/subsets
    else if (option=="-it") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      nb_iterations_subsets = (string)argv[i+1];
      i++;
    }
    // Optimizer settings
    else if (option=="-opti")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      options_optimizer = (string)argv[i+1];
      i++;
    }
    // Optimizer FOM
    else if (option=="-opti-fom")
    {
      optimizer_fom = true;
    }
    // Optimizer image update stat
    else if (option=="-opti-stat")
    {
      optimizer_stat = true;
    }
    // Penalty settings
    else if (option=="-penalty")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      options_penalty = (string)argv[i+1];
      i++;
    }
    // Image convolver settings
    else if (option=="-conv")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string convolver = (string)argv[i+1];
      options_image_convolver.push_back(convolver);
      i++;
    }
    // Image processing settings
    else if (option=="-proc")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string module = (string)argv[i+1];
      options_image_processing.push_back(module);
      i++;
    }

    // --------------------------------------------------------------------------------
    // Projection settings
    // --------------------------------------------------------------------------------

    // Projector settings
    else if (option=="-proj")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      options_projectorF = (string)argv[i+1];
      options_projectorB = (string)argv[i+1];
      i++;
    }
    // Forward projector settings
    else if (option=="-projF")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      options_projectorF = (string)argv[i+1];
      i++;
    }
    // Backward projector settings
    else if (option=="-projB")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      options_projectorB = (string)argv[i+1];
      i++;
    }
    // Common projector settings
    else if (option=="-proj-common")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      options_projector_common = (string)argv[i+1];
      i++;
    }
    // Ignore TOF flag
    else if (option=="-ignore-TOF")
    {
      ignore_TOF = true;
    }
    // Ignore POI flag
    else if (option=="-ignore-POI")
    {
      ignore_POI = true;
    }
    // Projection line computation strategy
    else if (option=="-proj-comp")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      projector_computation_strategy = atoi(argv[i+1]);
      i++;
    }

    // --------------------------------------------------------------------------------
    // Quantification settings
    // --------------------------------------------------------------------------------

    // Corrections settings
    else if (option=="-ignore-corr")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      ignored_corrections = (string)argv[i+1];
      i++;
    }

    // --------------------------------------------------------------------------------
    // Computation settings
    // --------------------------------------------------------------------------------

    // Flag to say that we want to use the GPU
    #ifdef CASTOR_GPU
    else if (option=="-gpu")
    {
      gpu_flag = 1;
    }
    #endif
    // Number of threads
    #ifdef CASTOR_OMP
    else if (option=="-th")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      nb_threads = (string)argv[i+1];
      i++;
    }
    #else
    else if (option=="-th")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!! castor-recon() -> Option -th is available only if the code is compiled using the CASTOR_OMP environment variable set to 1 !" << endl);
      Cerr("!!!!!                   The execution will continue BUT WITH ONLY ONE THREAD !" << endl);
      Cerr("!!!!!                   We strongly advice to compile CASToR with OpenMP !" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      i++;
    }
    #endif

    // --------------------------------------------------------------------------------
    // Unknown option!
    // --------------------------------------------------------------------------------

    else
    {
      Cerr("***** castor-recon() -> Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Affect specific verbose leves if not set
  if (verbose_algo==-1) verbose_algo = verbose_general;
  if (verbose_opti==-1) verbose_opti = verbose_general;
  if (verbose_proj==-1) verbose_proj = verbose_general;
  if (verbose_conv==-1) verbose_conv = verbose_general;
  if (verbose_proc==-1) verbose_proc = verbose_general;
  if (verbose_scan==-1) verbose_scan = verbose_general;
  if (verbose_data==-1) verbose_data = verbose_general;
  if (verbose_defo==-1) verbose_defo = verbose_general;
  if (verbose_dyna==-1) verbose_dyna = verbose_general;
  if (verbose_sens==-1) verbose_sens = verbose_general;

  // ============================================================================================================
  // Some checks
  // ============================================================================================================

  // Data files
  if (nb_beds < 1)
  {
    Cerr("***** castor-recon() -> Please provide at least one data filename !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Output files
  if (path_fout.empty() && path_dout.empty())
  {
    Cerr("***** castor-recon() -> Please provide an output option for output files (-fout or -dout) !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Check that only one option has been provided
  if (!path_fout.empty() && !path_dout.empty())
  {
    Cerr("***** castor-recon() -> Please provide either output option -fout or -dout but not both !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Check if gated reconstruction is enabled but no file describing the gating of the data has been provided
  // TODO: maybe do it in the dynamic data manager if possible to clean this main as much as possible
  if ( (nb_resp_gates>1 || nb_card_gates>1) && path_to_4D_data_splitting_file.empty() )
  {
    Cerr("***** castor-recon() -> gating is enabled, but no file describing the splitting of the data has been provided (-gating option) !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Singletons initialization: create here all needed singletons
  // ============================================================================================================
  
  if (verbose_general>=5) Cout("----- Singletons initializations ... -----" << endl); 

  // Get user endianness (interfile I/O)
  GetUserEndianness();
  
  // ----------------------------------------------------------------------------------------
  // Create sOutputManager
  // ----------------------------------------------------------------------------------------
  sOutputManager* p_outputManager = sOutputManager::GetInstance();  
  // Set verbose level
  p_outputManager->SetVerbose(verbose_general);
  // Set MPI rank
  p_outputManager->SetMPIRank(mpi_rank);
  // Set output dynamic image policy
  p_outputManager->SetMergeDynImagesFlag(merge_dynamic_imgs_flag);
  // Set datafile name(s) in order to be able to recover them from interfile classes
  p_outputManager->SetDataFileName(path_to_data_filename);
  // Set output number precision
  p_outputManager->SetOutNbPrec(onb_prec);
  
  
  // Set path to the config directory
  if (p_outputManager->CheckConfigDir(path_to_config_dir))
  {
    Cerr("***** castor-recon() -> A problem occured while checking for the config directory path !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize output directory and base name
  if (p_outputManager->InitOutputDirectory(path_fout, path_dout))
  {
    Cerr("***** castor-recon() -> A problem occured while initializing output directory !" << endl);
    Exit(EXIT_FAILURE);
  }
  // General call verbose
  #ifdef CASTOR_VERSION
  Cout("castor-recon() -> Launch reconstruction from CASToR version " << CASTOR_VERSION << "." << endl);
  #endif
  // Log command line
  if (p_outputManager->LogCommandLine(argc,argv))
  {
    Cerr("***** castor-recon() -> A problem occured while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ----------------------------------------------------------------------------------------
  // Create oImageDimensionsAndQuantification
  // ----------------------------------------------------------------------------------------

  // Have to create it and set the number of threads before the scanner manager is used to compute the LUT,
  // because some scanners use multi-threading for better efficiency while computing the LUT.
  // However, the dimensions must not be set here because if not provided in the command line, they will be
  // read from the scanner configuration file right below.
  oImageDimensionsAndQuantification* p_ImageDimensionsAndQuantification = new oImageDimensionsAndQuantification();
  if (p_ImageDimensionsAndQuantification->SetNbThreads(nb_threads))
  {
    Cerr("***** castor-recon() -> A problem occured while setting the number of threads !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ----------------------------------------------------------------------------------------
  // Create sScannerManager
  // ----------------------------------------------------------------------------------------
  sScannerManager* p_ScannerManager = sScannerManager::GetInstance();  
  p_ScannerManager->SetVerbose(verbose_scan);
  p_ScannerManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_ScannerManager->SetSaveLUTFlag(save_LUT_flag);

  if (verbose_general>=5) Cout("----- Geometry Initialization ... -----" << endl);

  // TODO: put all that stuff into only one function taking the path_to_data_filename[0] as the only parameter

  // Get system name from the dataFile
  string scanner_name = "";
  if (ReadDataASCIIFile(path_to_data_filename[0], "Scanner name", &scanner_name, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-recon() -> A problem occured while trying to find the system name in the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->FindScannerSystem(scanner_name) )
  {
    Cerr("***** castor-recon() -> A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->BuildScannerObject() )
  {
    Cerr("***** castor-recon() -> A problem occurred during scanner object construction ! !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->InstantiateScanner() )
  {
    Cerr("***** castor-recon() -> A problem occurred while creating Scanner object !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->GetGeometricInfoFromDataFile(path_to_data_filename[0]))
  {
    Cerr("***** castor-recon() -> A problem occurred while retrieving scanner fields from the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->BuildLUT() )
  {
    Cerr("***** castor-recon() -> A problem occurred while generating/reading the LUT !" << endl);
    Exit(EXIT_FAILURE);
  } 
  // Check the scanner manager parameters and initialize the scanner
  if (p_ScannerManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking scanner manager parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing scanner !" << endl);
    Exit(EXIT_FAILURE);
  }

  // If no number of voxels provided, then get the default ones from the scanner
  if (nb_voxX<=0 || nb_voxY<=0 || nb_voxZ<=0)
  {
    if (ReadDataASCIIFile(p_ScannerManager->GetPathToScannerFile(), "voxels number transaxial", &nb_voxX, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** castor-recon() -> A problem occured while reading for default number of transaxial voxels !" << endl);
      Exit(EXIT_FAILURE);
    }
    if (ReadDataASCIIFile(p_ScannerManager->GetPathToScannerFile(), "voxels number transaxial", &nb_voxY, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** castor-recon() -> A problem occured while reading for default number of transaxial voxels !" << endl);
      Exit(EXIT_FAILURE);
    }
    if (ReadDataASCIIFile(p_ScannerManager->GetPathToScannerFile(), "voxels number axial", &nb_voxZ, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** castor-recon() -> A problem occured while reading for default number of axial voxels !" << endl);
      Exit(EXIT_FAILURE);
    }
  }
  // If no FOV nor VOX size provided, then get the default one
  if ( (fov_sizeX<=0 || fov_sizeY<=0 || fov_sizeZ<=0) && (vox_sizeX<=0 || vox_sizeY<=0 || vox_sizeZ<=0) )
  {
    if (ReadDataASCIIFile(p_ScannerManager->GetPathToScannerFile(), "field of view transaxial", &fov_sizeX, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** castor-recon() -> A problem occured while reading for default transaxial FOV size !" << endl);
      Exit(EXIT_FAILURE);
    }
    if (ReadDataASCIIFile(p_ScannerManager->GetPathToScannerFile(), "field of view transaxial", &fov_sizeY, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** castor-recon() -> A problem occured while reading for default transaxial FOV size !" << endl);
      Exit(EXIT_FAILURE);
    }
    if (ReadDataASCIIFile(p_ScannerManager->GetPathToScannerFile(), "field of view axial", &fov_sizeZ, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** castor-recon() -> A problem occured while reading for default axial FOV size !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  if (verbose_general>=5) Cout("----- Geometry Initialization  OK -----" << endl);

  // ----------------------------------------------------------------------------------------
  // Create oImageDimensionsAndQuantification
  // ----------------------------------------------------------------------------------------
  
  if (verbose_general>=5) Cout("----- Image dimensions initialization ... -----" << endl);

  p_ImageDimensionsAndQuantification->SetNbBeds(nb_beds);
  p_ImageDimensionsAndQuantification->SetNbVoxX(nb_voxX);
  p_ImageDimensionsAndQuantification->SetNbVoxY(nb_voxY);
  p_ImageDimensionsAndQuantification->SetNbVoxZ(nb_voxZ);
  p_ImageDimensionsAndQuantification->SetVoxSizeX(vox_sizeX);
  p_ImageDimensionsAndQuantification->SetVoxSizeY(vox_sizeY);
  p_ImageDimensionsAndQuantification->SetVoxSizeZ(vox_sizeZ);
  p_ImageDimensionsAndQuantification->SetFOVSizeX(fov_sizeX);
  p_ImageDimensionsAndQuantification->SetFOVSizeY(fov_sizeY);
  p_ImageDimensionsAndQuantification->SetFOVSizeZ(fov_sizeZ);
  p_ImageDimensionsAndQuantification->SetFOVOutMasking(fov_out,slice_out);
  p_ImageDimensionsAndQuantification->SetOffsetX(offsetX);
  p_ImageDimensionsAndQuantification->SetOffsetY(offsetY);
  p_ImageDimensionsAndQuantification->SetOffsetZ(offsetZ);
  p_ImageDimensionsAndQuantification->SetMPIRankAndSize(mpi_rank, mpi_size);
  p_ImageDimensionsAndQuantification->SetVerbose(verbose_data);
  p_ImageDimensionsAndQuantification->SetIgnoredCorrections(ignored_corrections);
  p_ImageDimensionsAndQuantification->SetFrames(frame_list);
  p_ImageDimensionsAndQuantification->SetNbTimeBasisFunctions(nb_time_basis);
  p_ImageDimensionsAndQuantification->SetTimeBasisFunctionsFile(path_to_time_basis_coef);
  p_ImageDimensionsAndQuantification->SetNbMultiModalImages(path_to_multimodal_img.size());
  if (p_ImageDimensionsAndQuantification->SetFlipOut(flip_out))
  {
    Cerr("***** castor-recon() -> A problem occured while setting the output flip option !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (resp_motion_options=="" && double_motion_options=="")
  {
    p_ImageDimensionsAndQuantification->SetNbRespGates(nb_resp_gates);
    p_ImageDimensionsAndQuantification->SetNbRespBasisFunctions(nb_resp_basis);
    p_ImageDimensionsAndQuantification->SetRespBasisFunctionsFile(path_to_resp_basis_coef);
  }
  else
  {
    if (path_to_resp_basis_coef!="")
    {
      Cerr("***** castor-recon() -> Cannot use both respiratory motion correction and respiratory basis functions, it has no sense !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Set only one gate here because we correct for motion, so we reconstruct only one image
    p_ImageDimensionsAndQuantification->SetNbRespGates(1);
  }
  if (card_motion_options=="" && double_motion_options=="")
  {
    p_ImageDimensionsAndQuantification->SetNbCardGates(nb_card_gates);
    p_ImageDimensionsAndQuantification->SetNbCardBasisFunctions(nb_card_basis);
    p_ImageDimensionsAndQuantification->SetCardBasisFunctionsFile(path_to_card_basis_coef);
  }
  else
  {
    if (path_to_card_basis_coef!="")
    {
      Cerr("***** castor-recon() -> Cannot use both cardiac motion correction and cardiac basis functions, it has no sense !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Set only one gate here because we correct for motion, so we reconstruct only one image
    p_ImageDimensionsAndQuantification->SetNbCardGates(1);
  }
  if (p_ImageDimensionsAndQuantification->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking image dimensions parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ImageDimensionsAndQuantification->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing image dimensions !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialization of DynamicDataManager class, related 4D data splitting management 
  if (p_ImageDimensionsAndQuantification->InitDynamicData(path_to_4D_data_splitting_file,
                                                            !resp_motion_options.empty(), 
                                                            !card_motion_options.empty(), 
                                                          !double_motion_options.empty(), 
                                                            !ipat_motion_options.empty(), 
                                                                           nb_resp_gates, 
                                                                           nb_card_gates) )
  {
    Cerr("***** castor-recon() -> A problem occured while initializing Dynamic data manager's class !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Get the number of events in the data from the header (in order to check consistency between
  // the number of events in the datafile and in the dynamic files, if any gating is enabled)
  int64_t nb_events = 0;
  ReadDataASCIIFile(path_to_data_filename[0], "Number of events", &nb_events, 1, KEYWORD_MANDATORY);
  // Check dynamic parameters
  if (p_ImageDimensionsAndQuantification->CheckDynamicParameters(nb_events) )
  {
    Cerr("***** castor-recon() -> A problem occured while checking Dynamic data manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize dynamic specific quantitative factors
  if (p_ImageDimensionsAndQuantification->SetDynamicSpecificQuantificationFactors(path_to_dynamic_quantification_file) )
  {
    Cerr("***** castor-recon() -> A problem occured while initializing specific dynamic quantification factors!" << endl);
    Exit(EXIT_FAILURE);
  } 

  if (verbose_general>=5) Cout("----- Image dimensions initialization OK -----" << endl);

  // ----------------------------------------------------------------------------------------
  // Create the image space
  // ----------------------------------------------------------------------------------------

  oImageSpace* p_ImageSpace = new oImageSpace();
  p_ImageSpace->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_ImageSpace->SetVerbose(verbose_data);

  // ----------------------------------------------------------------------------------------
  // Create sChronoManager
  // ----------------------------------------------------------------------------------------

  if (verbose_general>=5) Cout("----- Profiling manager initialization ... -----" << endl);

  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();
  p_ChronoManager->SetNbThreads( p_ImageDimensionsAndQuantification->GetNbThreadsForProjection(),
                                 p_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation() );
  p_ChronoManager->SetVerbose(0);
  if (p_ChronoManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking parameters for the chrono manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ChronoManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing the chrono manager !" << endl);
    Exit(EXIT_FAILURE);
  }

  if (verbose_general>=5) Cout("----- Profiling manager initialization OK -----" << endl);

  // ----------------------------------------------------------------------------------------
  // Random Number Generator initialization: (we first require to know the number of threads to use from p_ImageDimensionsAndQuantification)
  // ----------------------------------------------------------------------------------------
  
  if (verbose_general>=5) Cout("----- Random number generator initialization ... -----" << endl);
  
  sRandomNumberGenerator* p_RandomNumberGenerator = sRandomNumberGenerator::GetInstance(); 
  p_RandomNumberGenerator->SetVerbose(verbose_general);
  // Use a user-provided seed to initialize the RNG if one has been provided. Use random number otherwise.
  if (random_generator_seed>=0) p_RandomNumberGenerator->Initialize(random_generator_seed, p_ImageDimensionsAndQuantification->GetNbThreadsMax(), nb_extra_random_generators);
  else p_RandomNumberGenerator->Initialize(p_ImageDimensionsAndQuantification->GetNbThreadsMax(), nb_extra_random_generators);
  
  if (verbose_general >=5) Cout("----- Random number generator initialization OK -----" << endl);

  // ----------------------------------------------------------------------------------------
  // Create vDataFile
  // ----------------------------------------------------------------------------------------

  if (verbose_general>=5) Cout("----- DataFile initialization ... -----" << endl);
  
  vDataFile** p_DataFile = new vDataFile*[nb_beds];

  if (p_ScannerManager->GetScannerType() == SCANNER_PET)
  {
    // Create specific data file
    for (int i=0 ; i<nb_beds ; i++)
    {
      p_DataFile[i] = new iDataFilePET();
      (dynamic_cast<iDataFilePET*>(p_DataFile[i]))->SetIgnoreTOFFlag(ignore_TOF);
    }
  }
  else if (p_ScannerManager->GetScannerType() == SCANNER_SPECT_CONVERGENT)
  {
    // Create specific data file
    for (int i=0 ; i<nb_beds ; i++)
    {
      p_DataFile[i] = new iDataFileSPECT(); 
    }
  }
  else if (p_ScannerManager->GetScannerType() == SCANNER_CT)
  {
    // Create specific data file
    for (int i=0 ; i<nb_beds ; i++)
    {
      p_DataFile[i] = new iDataFileCT(); 
    }
  }
  // Unknown scanner
  else
  {
    Cerr("***** castor-recon() -> Unknown scanner type (" << p_ScannerManager->GetScannerType() << ") for datafile construction ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }

  // Load raw data in memory and do other stuff if needed.
  for (int bed=0 ; bed<nb_beds ; bed++)
  {
    // If it is asked to revert the order of the bed positions, then we proceed here
    if (invert_datafile_bed_order_flag) p_DataFile[bed]->SetHeaderDataFileName(path_to_data_filename.at(nb_beds-1-bed));
    else p_DataFile[bed]->SetHeaderDataFileName(path_to_data_filename.at(bed));
    p_DataFile[bed]->SetBedIndex(bed);
    p_DataFile[bed]->SetVerbose(verbose_data);
    p_DataFile[bed]->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
    p_DataFile[bed]->SetIgnorePOIFlag(ignore_POI);
    if (p_DataFile[bed]->ReadInfoInHeader())
    {
      Cerr("***** castor-recon() -> A problem occurred during datafile header reading ! Abort." << endl);
      Exit(EXIT_FAILURE);
    }
    if (p_DataFile[bed]->CheckParameters())
    {
      Cerr("***** castor-recon() -> A problem occurred while checking datafile parameters ! Abort." << endl);
      Exit(EXIT_FAILURE);
    }
    if (p_DataFile[bed]->ComputeSizeEvent())
    {
      Cerr("***** castor-recon() -> A problem occurred in datafile initialization ! Abort." << endl);
      Exit(EXIT_FAILURE);
    }
    if (p_DataFile[bed]->InitializeMappedFile())
    {
      Cerr("***** castor-recon() -> A problem occurred in datafile initialization ! Abort." << endl);
      Exit(EXIT_FAILURE);
    }
    if (p_DataFile[bed]->PrepareDataFile())
    {
      Cerr("***** castor-recon() -> A problem occured in datafile preparation ! Abort." << endl);
      Exit(EXIT_FAILURE);
    }
  }
  // Check consistency between all datafiles; to do that we compare the first with all the others
  for (int bed=1; bed<nb_beds; bed++)
  {
    if (p_DataFile[0]->CheckConsistencyWithAnotherBedDataFile(p_DataFile[bed]))
    {
      int bed_index_problem = bed + 1;
      if (invert_datafile_bed_order_flag) bed_index_problem = nb_beds - bed;
      Cerr("***** castor-recon() -> A problem occured while checking consistency between first bed and bed " << bed_index_problem << " !" << endl);
      Exit(EXIT_FAILURE);
    }
  }
  // And finally, deal with the multiple bed positions
  if (p_ImageDimensionsAndQuantification->DealWithBedPositions(p_DataFile))
  {
    Cerr("***** castor-recon() -> A problem occured while dealing with the bed positions !" << endl);
    Exit(EXIT_FAILURE);
  }

  if (verbose_general>=5) Cout("----- DataFile initialization OK -----" << endl);

  // ----------------------------------------------------------------------------------------
  // Create Projector Manager
  // ----------------------------------------------------------------------------------------

  // Verbose
  if (verbose_general>=5) Cout("----- Projector initialization ... -----" << endl);
  // Create object
  oProjectorManager* p_ProjectorManager = new oProjectorManager();
  // Set all parameters
  p_ProjectorManager->SetScanner(p_ScannerManager->GetScannerObject());
  p_ProjectorManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_ProjectorManager->SetDataFile(p_DataFile[0]);
  p_ProjectorManager->SetComputationStrategy(projector_computation_strategy);
  p_ProjectorManager->SetOptionsForward(options_projectorF);
  p_ProjectorManager->SetOptionsBackward(options_projectorB);
  p_ProjectorManager->SetOptionsCommon(options_projector_common);
  
  if (!path_to_mask_img.empty())
  {
    if (p_ImageSpace->InitMaskImage(path_to_mask_img))
    {
      Cerr("The option for using a mask image for projection is set, but the mask image is not provided or could not be loaded !" << endl);
      Exit(EXIT_FAILURE);
    }
    p_ProjectorManager->ProcessAndSetMask(p_ImageSpace->mp_maskImage);
    p_ImageSpace->DeallocateMaskImage();
  }

  p_ProjectorManager->SetVerbose(verbose_proj);
  // Check parameters
  if (p_ProjectorManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking projector manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize projector manager
  if (p_ProjectorManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing projector manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Check specific requirements for SPECT with attenuation correction
  if (p_ProjectorManager->CheckSPECTAttenuationCompatibility(path_to_attenuation_img))
  {
    Cerr("***** castor-recon() -> A problem occured while checking projector's compatibility with SPECT and attenuation correction !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (verbose_general>=5) Cout("----- Projector initialization OK -----" << endl);
  
  // ----------------------------------------------------------------------------------------
  // Create Optimizer Manager
  // ----------------------------------------------------------------------------------------

  // Verbose
  if (verbose_general>=5) Cout("----- Optimizer initialization ... -----" << endl);
  // Create object
  oOptimizerManager* p_OptimizerManager = new oOptimizerManager();
  // Set all parameters
  p_OptimizerManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_OptimizerManager->SetImageSpace(p_ImageSpace);
  p_OptimizerManager->SetDataMode(p_DataFile[0]->GetDataMode());
  p_OptimizerManager->SetDataType(p_DataFile[0]->GetDataType());
  p_OptimizerManager->SetDataSpec(p_DataFile[0]->GetDataSpec());
  p_OptimizerManager->SetOptionsOptimizer(options_optimizer);
  p_OptimizerManager->SetOptimizerFOMFlag(optimizer_fom);
  p_OptimizerManager->SetOptimizerImageStatFlag(optimizer_stat);
  p_OptimizerManager->SetOptionsPenalty(options_penalty);
  p_OptimizerManager->SetNbTOFBins(p_ProjectorManager->GetNbTOFBins());
  p_OptimizerManager->SetVerbose(verbose_opti);
  // Check parameters
  if (p_OptimizerManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking optimizer manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize optimizer manager
  if (p_OptimizerManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing optimizer manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (verbose_general>=5) Cout("----- Optimizer initialization OK -----" << endl);
  
  // ----------------------------------------------------------------------------------------
  // Create Image Convolver Manager
  // ----------------------------------------------------------------------------------------

  // Verbose
  if (verbose_general>=5) Cout("----- Image Convolver initialization (if any) ... -----" << endl);
  // Create object
  oImageConvolverManager* p_ImageConvolverManager = new oImageConvolverManager();
  // Set all parameters
  p_ImageConvolverManager->SetVerbose(verbose_conv);
  p_ImageConvolverManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_ImageConvolverManager->SetOptions(options_image_convolver);
  // Check parameters
  if (p_ImageConvolverManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking image convolver manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize image convolver manager
  if (p_ImageConvolverManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing image convolver manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (verbose_general>=5) Cout("----- Image Convolver initialization OK -----" << endl);

  // ----------------------------------------------------------------------------------------
  // Create Image Processing Manager
  // ----------------------------------------------------------------------------------------

  // Verbose
  if (verbose_general>=5) Cout("----- Image Processing initialization (if any) ... -----" << endl);
  // Create object
  oImageProcessingManager* p_ImageProcessingManager = new oImageProcessingManager();
  // Set all parameters
  p_ImageProcessingManager->SetVerbose(verbose_proc);
  p_ImageProcessingManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_ImageProcessingManager->SetOptions(options_image_processing);
  // Check parameters
  if (p_ImageProcessingManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking image processing manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize image processing manager
  if (p_ImageProcessingManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing image processing manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (verbose_general>=5) Cout("----- Image Processing initialization OK -----" << endl);

  // ----------------------------------------------------------------------------------------
  // Create Dynamic Model Manager 
  // ----------------------------------------------------------------------------------------

  // Verbose
  if (verbose_general>=5) Cout("----- Dynamic model initialization (if any) ... -----" << endl);
  // Create object
  oDynamicModelManager* p_DynamicModelManager = new oDynamicModelManager();
  // Set all parameters
  p_DynamicModelManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_DynamicModelManager->SetOptions(dynamic_model_options);
  p_DynamicModelManager->SetVerbose(verbose_dyna);
  // Check parameters
  if (p_DynamicModelManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking dynamic model manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize optimizer manager
  if (p_DynamicModelManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing dynamic model manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (verbose_general>=5) Cout("----- Dynamic model initialization OK -----" << endl);
  
  // ----------------------------------------------------------------------------------------
  // Create Deformation Manager 
  // ----------------------------------------------------------------------------------------

  // Verbose
  if (verbose_general>=5) Cout("----- Image deformation initialization (if any) ... -----" << endl);
  // Create object
  oDeformationManager* p_DeformationManager = new oDeformationManager();
  // Set all parameters
  p_DeformationManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_DeformationManager->SetDataMode(p_DataFile[0]->GetDataMode()); // required to know if sensitivity image deformation should be enabled
  
  if(resp_motion_options != "")
  {
    p_DeformationManager->SetOptions(resp_motion_options);
    p_DeformationManager->SetNbTransformations(nb_resp_gates-1);
    p_DeformationManager->SetMotionType(DEF_RESP_MOT);
  }
  else if (card_motion_options != "")
  {
    p_DeformationManager->SetOptions(card_motion_options);
    p_DeformationManager->SetNbTransformations(nb_card_gates-1);
    p_DeformationManager->SetMotionType(DEF_CARD_MOT);
  }
  else if (double_motion_options != "")
  {
    p_DeformationManager->SetOptions(double_motion_options);
    p_DeformationManager->SetNbTransformations(nb_resp_gates*nb_card_gates-1);
    p_DeformationManager->SetMotionType(DEF_DUAL_MOT);
  }
  else if (ipat_motion_options != "")
  {
    p_DeformationManager->SetOptions(ipat_motion_options);
    p_DeformationManager->SetNbTransformations(p_ImageDimensionsAndQuantification->GetNbIPatMotionSubsets() - 1);
    p_DeformationManager->SetMotionType(DEF_IPAT_MOT);
  }
  else
    p_DeformationManager->SetNbTransformations(0);

  p_DeformationManager->SetVerbose(verbose_defo);
  // Check parameters
  if (p_DeformationManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occured while checking image deformation manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize optimizer manager
  if (p_DeformationManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occured while initializing image deformation manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (verbose_general >=5) Cout("----- Image deformation initialization OK -----" << endl);

  // ============================================================================================================
  // Check for supported/implemented combinations of options, after creating managers/input data and before launching the algorithm
  // ============================================================================================================

  // ============================================================================================================
  // Algorithm initialization: create here sensitivity computation and launch algorithm
  // ============================================================================================================
  
  // --------------------------------------------------------------------------------------------
  // Create sensitivity if one of the following conditions is met
  //  1. Input file is a list-mode file and a sensitivity image was not provided
  //  2. Input file is a histogram and the -sens-histo option is provided (i.e. we explicitely
  //     ask to compute the global sensitivity from the histogram datafile)
  // --------------------------------------------------------------------------------------------
  if ( (path_to_sensitivity_img.empty() && p_DataFile[0]->GetDataMode()==MODE_LIST) ||
       (p_DataFile[0]->GetDataMode()==MODE_HISTOGRAM && sensitivity_from_histogram)
     )
  {
    // Verbose
    if (verbose_general>=5) Cout("----- Image Sensitivity generation initialization ... -----" << endl);
    // Create object
    oSensitivityGenerator* p_Sensitivity = new oSensitivityGenerator();
    // Set parameters
    p_Sensitivity->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
    p_Sensitivity->SetImageSpace(p_ImageSpace);
    p_Sensitivity->SetScanner(p_ScannerManager->GetScannerObject());
    p_Sensitivity->SetProjectorManager(p_ProjectorManager);
    p_Sensitivity->SetImageConvolverManager(p_ImageConvolverManager);
    p_Sensitivity->SetDeformationManager(p_DeformationManager);
    p_Sensitivity->SetGPUflag(gpu_flag);
    p_Sensitivity->SetPathToAttenuationImage(path_to_attenuation_img);
    p_Sensitivity->SetPathToMaskImage(path_to_mask_img);
    p_Sensitivity->SetNumberOfAtnGateImages(nb_atn_resp_imgs, nb_atn_card_imgs);
    p_Sensitivity->SetPathToNormalizationFileName(path_to_normalization_filename,invert_datafile_bed_order_flag);
    p_Sensitivity->SetDataFile(p_DataFile);
    p_Sensitivity->SetComputeFromHistogramFlag(sensitivity_from_histogram);
    p_Sensitivity->SetVerbose(verbose_sens);
    // Check parameters
    if (p_Sensitivity->CheckParameters())
    {
      Cerr("***** castor-recon() -> A problem occured while checking parameters of the sensitivity generator !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Initialize the sensitivity generator
    if (p_Sensitivity->Initialize())
    {
      Cerr("***** castor-recon() -> A problem occured while initializing the sensitivity generator !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Set the sensitivity mode ON for the projector manager
    p_ProjectorManager->SetSensitivityModeOn();
    // Launch the computation
    if (p_Sensitivity->Launch())
    {
      Cerr("***** castor-recon() -> A problem occured while computing the sensitivity !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Set the sensitivity mode OFF for the projector manager
    p_ProjectorManager->SetSensitivityModeOff();
    // Get the path to the sensitivity image (will be given to the algorithm if input file is a list-mode)
    if (p_DataFile[0]->GetDataMode()==MODE_LIST) path_to_sensitivity_img = p_Sensitivity->GetPathToSensitivityImage();
    else path_to_sensitivity_img = "";
    // Delete the generator
    delete p_Sensitivity;
    // Exit now if asked for
    if (exit_after_sensitivity)
    {
      if (verbose_general>=VERBOSE_LIGHT) Cout("castor-recon() -> Asked to exit after sensitivity computation." << endl);
      // Delete objects in the inverse order in which they were created
      delete p_ImageSpace;
      delete p_DeformationManager;
      delete p_DynamicModelManager;
      delete p_ImageProcessingManager;
      delete p_ImageConvolverManager;
      delete p_OptimizerManager;
      delete p_ProjectorManager;
      for (int i=0 ; i<nb_beds ; i++) delete p_DataFile[i]; 
      delete[] p_DataFile;
      delete p_ImageDimensionsAndQuantification;
      // And exit
      return 0;
    }
  }
  
  // ----------------------------------------------------------------------------------------
  // Create algorithm
  // ----------------------------------------------------------------------------------------

  // If the number of events in a datafile is below the number of threads for projections, then we must reduce this number of threads, otherwise
  // the datafile reading cannot work (and it has no sense anyway). This is done by the following function
  p_ImageDimensionsAndQuantification->CheckNumberOfProjectionThreadsConsistencyWithDataFileSize(p_DataFile);
  // Verbose
  if (verbose_general>=5) Cout("----- Iterative reconstruction algorithm initialization ... -----" << endl);
  // Create object
  oIterativeAlgorithm* p_Algorithm = new oIterativeAlgorithm();
  // Set parameters
  p_Algorithm->SetImageConvolverManager(p_ImageConvolverManager);
  p_Algorithm->SetImageProcessingManager(p_ImageProcessingManager);
  p_Algorithm->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_Algorithm->SetImageSpace(p_ImageSpace);
  p_Algorithm->SetProjectorManager(p_ProjectorManager);
  p_Algorithm->SetDynamicModelManager(p_DynamicModelManager);
  p_Algorithm->SetDeformationManager(p_DeformationManager);
  p_Algorithm->SetOptimizerManager(p_OptimizerManager);
  p_Algorithm->SetDataFile(p_DataFile);
  p_Algorithm->SetGPUflag(gpu_flag);
  p_Algorithm->SetVerbose(verbose_algo);
  p_Algorithm->SetNbBeds(nb_beds);
  p_Algorithm->SetPathInitImage(path_to_initial_img);
  p_Algorithm->SetPathToAttenuationImage(path_to_attenuation_img);
  p_Algorithm->SetPathToSensitivityImage(path_to_sensitivity_img);
  p_Algorithm->SetPathToMultiModalImage(path_to_multimodal_img);
  p_Algorithm->SetPathToMaskImage(path_to_mask_img);
  p_Algorithm->SetSaveSensitivityHistoFlag(save_sens_histo);
  p_Algorithm->SetSaveSubsetImageFlag(save_subset_image);
  if (p_Algorithm->SetNbIterationsAndSubsets(nb_iterations_subsets))
  {
    Cerr("***** castor-recon() -> Error while setting the numbers of iterations and subsets !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_Algorithm->SetOutputIterations(output_iterations))
  {
    Cerr("***** castor-recon() -> Error while setting the selected output iterations !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_Algorithm->Iterate())
  {
    Cerr("***** castor-recon() -> Error while performing the reconstruction" << endl);
    Exit(EXIT_FAILURE);
  }

  // Display profiling results
  p_ChronoManager->Display();

  // ============================================================================================================
  // End
  // ============================================================================================================

  // Delete objects in the inverse order in which they were created
  if (verbose_general>=5) Cout("----- Deleting CASToR objects ... -----" << endl);
  delete p_Algorithm;
  delete p_ImageSpace;
  delete p_DeformationManager;
  delete p_DynamicModelManager;
  delete p_ImageProcessingManager;
  delete p_ImageConvolverManager;
  delete p_OptimizerManager;
  delete p_ProjectorManager;
  for (int i=0 ; i<nb_beds ; i++) delete p_DataFile[i];
  delete[] p_DataFile;
  delete p_ImageDimensionsAndQuantification;
  if (verbose_general>=5) Cout("----- CASToR objects successfully deleted -----" << endl);
    
  // Ending
  if (verbose_general>=1) Cout(endl);
  #ifdef CASTOR_MPI
  MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}
