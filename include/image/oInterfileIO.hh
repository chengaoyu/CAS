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
  \brief    This group of functions manages Interfile image file format
  \details  These functions can read and write interfile header and image files. It contains an interfile dictionnary.
*/

#ifndef OINTERFILEIO_HH
#define OINTERFILEIO_HH 1

#include "oImageDimensionsAndQuantification.hh"
#include "sScannerManager.hh"

// ---------------------------------------------------------------------
/**
 * @defgroup INTF Interfile I/O related variables
 *
 * \brief This group contains variables and structures related to Interfile I/O. \n
 *        Those are either variables associated to certain Interfile key,
 *        type of data, and structures dedicated to recover those keys. \n
 *        Defined in oInterfileIO.hh
 */


/**
 * \ingroup INTF
 * @defgroup INTF_ENDIANNESS Data endianness
 *   \brief Variables related to Big and Little endian
 */
 /*@{*/
/** Variable related to Big Endian data (=0) */
#define INTF_BIG_ENDIAN 0
/** Variable related to Little Endian data (=1) */
#define INTF_LITTLE_ENDIAN 1
/** @} */


/**
 * \ingroup INTF
 * @defgroup INTF_IMG_TYPE Image data type according to v3.3 (+PET)
 *   \brief Variables related to the Interfile type of data
 */
/*@{*/
/** Variable corresponding to a static image (=0, default) */
#define INTF_IMG_STATIC 0
/** Variable corresponding to a dynamic image (=1) */
#define INTF_IMG_DYNAMIC 1
/** Variable corresponding to a (static) PET image (=2) */
#define INTF_IMG_PET 2
/** Variable corresponding to a (static) SPECT image (=3) */
#define INTF_IMG_SPECT 3
/** Variable corresponding to a gated image (=4) */
#define INTF_IMG_GATED 4
/** Variable corresponding to a gated SPECT image (=5) */
#define INTF_IMG_GSPECT 5
/** Variable corresponding to an unknown data type (=6) */
#define INTF_IMG_UNKNOWN 6
/** @} */


/**
 * \ingroup INTF
 * @defgroup INTF_LERP Flag for input image interpolation
 *   \brief These flags are passed as last parameter to the Interfile
 *          file reading function.\n
 *          Depending on the flag, the input image will be interpolated
 *          if its dimensions are different to the reconstruction dimensions,
 *          or an error will be thrown.\n
 *          Current implemented interpolation is linear (LERP)
 */
/*@{*/
/** Interpolation forbidden (for input image reading) */
#define INTF_LERP_DISABLED false
/** Interpolation allowed (for input image reading) */
#define INTF_LERP_ENABLED true
/** @} */


/**
 * \ingroup INTF
 * @defgroup INTF_NUMBER_FORMAT String related to the possible pixel/voxel data 
 *
 *  \brief "ASCII" and "BIT" are not supported in the current implementation
 *         "unsigned integer" read by default
 */
/*@{*/
/** String corresponding to an unsigned bit type (="bit") */
#define BIT_str "bit"
/** String corresponding to an unsigned interger type (="") */
#define UINT32_str "unsigned integer"
/** String corresponding to a signed integer type (="") */
#define INT32_str "signed integer"
/** String corresponding to a float type (="short float / float") */
#define FLT32_str "short float"
#define FLT32_str2 "float"
/** String corresponding to a double type (="long float) */
#define FLT64_str "long float"
/** String corresponding to a long double type (="long long float) */
#define LONGDOUBLE_str "long long float"
/** String corresponding to an ASCII type (="ASCII") */
#define ASCII_str "ASCII"
/** @} */


/**
 * \ingroup INTF
 * @defgroup PATIENT_ORIENTATION Global patient orientation
 *
 *   \brief Combination of Patient rotation, orientation anf slice orientation. Not used in the current implementation
 */
/*@{*/
/** String corresponding to Transverse / Supine / Head-in orientation (=0) (default) */
#define INTF_SUPINE_HEADIN_TRANSAXIAL 0
/** String corresponding to Sagittal / Supine / Head-in orientation (=1) (default) */
#define INTF_SUPINE_HEADIN_SAGITTAL 1
/** String corresponding to Coronal / Supine / Head-in orientation (=2) (default) */
#define INTF_SUPINE_HEADIN_CORONAL 2 
/** String corresponding to Transverse / Supine / Feet-in orientation (=3) (default) */
#define INTF_SUPINE_FEETIN_TRANSAXIAL 3
/** String corresponding to Sagittal / Supine / Feet-in orientation (=4) (default) */
#define INTF_SUPINE_FEETIN_SAGITTAL 4
/** String corresponding to Coronal / Supine / Feet-in orientation (=5) (default) */
#define INTF_SUPINE_FEETIN_CORONAL 5
/** String corresponding to Transverse / Prone / Head-in orientation (=6) (default) */
#define INTF_PRONE_HEADIN_TRANSAXIAL 6
/** String corresponding to Sagittal / Prone / Head-in orientation (=7) (default) */
#define INTF_PRONE_HEADIN_SAGITTAL 7
/** String corresponding to Coronal / Prone / Head-in orientation (=8) (default) */
#define INTF_PRONE_HEADIN_CORONAL 8
/** String corresponding to Transverse / Prone / Feet-in orientation (=9) (default) */
#define INTF_PRONE_FEETIN_TRANSAXIAL 9
/** String corresponding to Sagittal / Prone / Feet-in orientation (=10) (default) */
#define INTF_PRONE_FEETIN_SAGITTAL 10
/** String corresponding to Coronal / Prone / Feet-in orientation (=11) (default) */
#define INTF_PRONE_FEETIN_CORONAL 11
/** @} */


/**
 * \ingroup INTF
 * @defgroup INTF_SLICE_ORIENTATION Image slice orientation
 *
 *   \brief  Variable referring to the slice orientation (transverse/coronal/sagittal)
 *           in the image \n
 *           Not used in the current implementation
 */
/*@{*/
/** String corresponding to transverse slice orientation (=0) (default) */
#define INTF_TRANSVERSE 0
/** String corresponding to coronal slice orientation (=1)*/
#define INTF_CORONAL 1
/** String corresponding to sagittal slice orientation (=2) */
#define INTF_SAGITTAL 2
/** String corresponding to other slice orientation (=3) */
#define INTF_OTHER 3
/** @} */


/**
 * \ingroup INTF
 * @defgroup INTF_PATIENT_ROTATION Image patient position
 *
 *   \brief  Variable referring to the patient position (supine/prone)
 *           in the image \n
 *           Not used in the current implementation
 */
/*@{*/
/** String corresponding to supine patient position (=0) (default) */
#define INTF_SUPINE 0
/** String corresponding to prone patient position (=1)*/
#define INTF_PRONE 1
/** @} */


/**
 * \ingroup INTF
 * @defgroup INTF_PATIENT_ORIENTATION Image patient orientation
 *
 *   \brief  Variable referring to the slice orientation (head-in/feet-in)
 *           in the image \n
 *           Not used in the current implementation
 */
/*@{*/
/** String corresponding to head-in patient orientation (=0) (default) */
#define INTF_HEADIN 0
/** String corresponding to feet-in patient orientation (=1)*/
#define INTF_FEETIN 1
/** @} */



/*!
  \class   Intf_key
  \brief   Interfile key elements. \n
           This structure is used to recover and process the elements
           of an Interfile key ( key := value #anycomment) \n
           Declared in oInterfileIO.hh
*/
struct Intf_key
{
  string korig;  /*!< original line w/ comment */
  string kcase;  /*!< the whole recovered key */
  string klcase; /*!< all lower cases, no space */
  string kvalue; /*!< value of the key, no space */
};
/** @} */


// ---------------------------------------------------------------------
/*!
  \class   Intf_fields
  \brief   Interfile fields. \n
           This structure contains all the Interfile keys currently managed by CASToR \n
           Declared in oInterfileIO.hh
*/
struct Intf_fields
{
  string path_to_image; /*!< Absolute or relative path to the image file. Will be converted to absolute path in the struct.\n
  Associated with 'name of data file' intf key. */

  string originating_system; /*!< Name of the system as used and defined in CASTOR.\n
  Associated to 'originating system' intf key. */

  uint8_t endianness; /*!< Indicate if bytes follow a Little Endian or Big Endian ordering.\n
  Associated with 'imagedata byte order' intf key. */
  
  uint32_t data_offset; /*!< Data offset in the image file, as defined by the fields  "data offset in bytes" or "data starting block".\n
  Associated with 'imagedata byte order' intf key.*/
  
  string nb_format; /*!< Type of each pixel/voxel data. Assumed to be similar for each slice. \n
  Associated with 'number format' intf key.*/
  
  uint8_t nb_dims; /*!< Total number of dimensions in the image data.\n
  =3 by default. Computed from the other interfile keys*/
  
  uint32_t mtx_size[7]; /*!< Dimensions of the images (x,y,z) and non-spatial dimensions, such as number of frames.\n
  Support for sinogram reading (and sinogram dimensions) is not yet implemented.\n
  Associated with the 'matrix size[x]' and 'number of time frames' intf key.*/
                        
  FLTNB vox_offset[3]; /*!< Image position offset (x,y,z).\n
  Associated with the 'first pixel offset (mm) [x]' intf key.*/
  
//  float vox_size[3]; /*!< Voxel dimensions in mm.\n
  FLTNB vox_size[3]; /*!< Voxel dimensions in mm.\n
  Associated with 'scaling factor (mm/pixel) [x]' intf key.*/
  
//  float slice_thickness_mm ;  /*!< Read from the key 'slice thickness (pixels)'.\n
  FLTNB slice_thickness_mm ;  /*!< Read from the key 'slice thickness (pixels)'.\n
                                   Conversely to the key name, the value is assumed to be given in mm.*/

  int nb_bed_positions; /*! Number of bed positions */
  FLTNB bed_relative_position; /*!< Bed relative position in mm */
  bool bed_position_provided; /*!< True if the relative bed position has been provided from the datafile */

  uint16_t ctr_to_ctr_separation;  /*!< Recovered from the key 'centre-centre slice separation (pixels)'. Not currently implemented.\n
                                        Used to compute gap between two slices (gap = ctr-ctr slice separation - slice_thickness_mm).*/
                                        
  uint16_t nb_time_frames ; /*!< Number of time frames in a dynamic series of images.\n
  Associated with 'number of time frames' intf key.*/
  
  uint16_t nb_resp_gates ; /*!< Number of respiratory gates in a dynamic series of images.\n
  Associated with 'number of respiratory gates' intf key (CASToR-specific key).*/
  
  uint16_t nb_card_gates ; /*!< Number of cardiac gates in a dynamic series of images.\n
  Associated with 'number of cardiac gates' intf key (CASToR-specific key).*/
  
  uint32_t nb_total_imgs; /*!< Total number of images in the associated data files. \n
  Associated with 'total number of images' intf key.*/
  
  uint8_t nb_bytes_pixel; /*!< number of bytes for each pixel/voxel of the image data.\n
  Associated with 'number of bytes per pixel' intf key.*/
  
  int8_t slice_orientation; /*!< slice orientation : transverse (=0, default), sagittal (=1), coronal (=2).\n
  Associated with 'slice orientation' intf key.\n
  NOT used in the current implementation. */
  
  int8_t pat_rotation;      /*!< patient rotation : supine (=0, default), prone (=1).\n
  Associated with 'patient rotation' intf key.\n
  NOT used in the current implementation. */
  
  int8_t pat_orientation;   /*!< slice orientation : head-in (=0, default), feet-in (=1).\n
  Associated with 'patient orientation' intf key.\n
  NOT used in the current implementation. */
  
  FLTNB rescale_slope; /*!< multiplicative calibration values. Initialized from 'rescale slope'.\n
  Associated with 'imagedata byte order' intf key.*/
  
  FLTNB quant_units; /*!< multiplicative calibration values. Initialized from 'quantification units'.\n
  Associated with 'imagedata byte order' intf key.*/
  
  FLTNB rescale_intercept; /*!< additive (intercept) calibration values.\n
  Associated with 'uantification units' and 'rescale intercept' intf keys.*/
  
  uint32_t cmtx_size[3]; /*!< (C)astor reconstruction dimensions of the images (x,y,z).\n
  Initialized only if interpolation is required on a recovered image.*/
  
  FLTNB cvox_size[3]; /*!< (C)astor voxel size of the reconstructed images (x,y,z).\n
  Initialized only if interpolation is required on a recovered image.*/

  FLTNB cvox_offset[3]; 
  /*!< (C)astor image position offset (x,y,z).\n
  Initialized only if interpolation is required on a recovered image.*/
  
  bool is_mtx_size_different; /*!< Bool indicating the original matric/voxel sizes are different to the reconstruction size, \n
                                   hence requiring interpolation if the image is read.*/

  int data_type; /*!< Type of projeted/reconstructed dataset (Static|Dynamic|Gated|Tomographic|Curve|ROI|GSPECT|Other).\n
  (Curve|ROI|Other) are considered as 'Static'.\n
  Associated with 'type of data' intf key.*/
  
//  float study_duration; /*!< Acquisition duration.\n
  FLTNB study_duration; /*!< Acquisition duration.\n
  Associated with 'study duration (sec)' intf key.*/
  
  vector<FLTNB> image_duration; /*!< "image duration (sec)". Duration of the frame. \n
  Should be specific to each time frame.\n
  Associated with 'image duration (sec)' intf key.*/

  vector<FLTNB> image_start_time; /*!< "image start time (sec)". Start time of the frame. \n
  Should be specific to each time frame.\n
  Associated with 'image start time (sec)' intf key.*/
  
  uint32_t nb_time_windows; /*!< Number of time windows. \n
  This key is not used in the current implementation.\n
  To define a number of gates, 'number of respiratory gates' and 'number of cardiac gates' must be used.\n
  Associated with 'number of time windows' intf key.*/ 
  
  string process_status; /*!< "acquired" or "reconstructed" (default).\n
  Associated with 'process status' intf key.*/
  
  // --- Fields related to an image in a group  --- // // todo
  
  uint32_t nb_img_in_frame_groups; /*!< number of images in a frame.\n
  Associated with 'number of images this frame group' intf key.*/
  
//  float image_pause; /*!< Pause between time windows ?\n
  //FLTNB image_pause; 
  /*!< Pause between time windows ?\n
  Should be specific to each time frame.\n
  Not used in the current implementation.\n
  Associated with 'pause between images (sec)' intf key.*/
  
  vector<FLTNB> frame_group_pause; /*!< ""pause between frame groups (sec)".\n
  Associated with 'pause between frame groups (sec) ' intf key.*/ 

  // --- SPECT and projection related data --- //
  uint16_t nb_detector_heads; /*!< Number of detector heads in the system.\n
  Associated with 'number of detector heads' intf key.*/
  
  uint32_t nb_energy_windows; /*!< Total number of energy windows.  \n
  Associated with 'number of energy windows' intf key.*/
  
  uint16_t nb_projections; /*!< Number of projection in acquired data. \n
  Associated with 'number of projections' intf key.*/
  
//  float extent_rotation; /*!< Angular span ex: 180, 360.\n
  FLTNB extent_rotation; /*!< Angular span ex: 180, 360.\n
  Associated with 'extent of rotation' intf key.*/
  
  string direction_rotation; /*!< Direction of rotation : CCW (counter-clockwise), CW (clockwise).\n
  Associated with 'direction of rotation' intf key.*/
  
//  float first_angle; /*!< Angle of the first view.\n
  FLTNB first_angle; /*!< Angle of the first view.\n
  Associated with 'start angle' intf key.*/
  
  string projection_angles; /*!< All projection angles (for each view).\n
  Associated with 'projection_angles' intf key.*/
  
  string radius; /*!< Distance between center of rotation and detector, for each view.\n
  Associated with 'Center of rotation to detector distance' and 'Radius' intf keys.*/
};
/** @} */



#include "gVariables.hh"
#include "sOutputManager.hh"
#include "gOptions.hh"

  // -------------------------------------------------------------------
  // ----- "PUBLIC" FUNCTIONS WHICH COULD AND SHOULD BE ANYWHERE IN THE CODE WHERE INTERFILE IS REQUIRED -----

    /*!
      \fn IntfKeyGetValueFromFile(const string& a_pathToHeader, const string& a_key, T* ap_return, int a_nbElts, int a_mandatoryFlag)
      \param a_pathToHeader : path to the interfile header
      \param a_key : the key to recover
      \param T* ap_return : template array in which the data will be recovered
      \param int a_nbElts : number of elements to recover
      \param int a_mandatoryFlag : flag indicating if the data to recover if mandatory (true) or optionnal (false)
      \brief Look for "a_nbElts" elts in the "a_pathToHeader" interfile header matching the "a_keyword" key passed as parameter and return the corresponding value(s) in the "ap_return" templated array.
      \details If more than one elements are to be recovered, the function first check the key has a correct Interfile kay layout (brackets and commas :  {,,})\n
               Depending on the mandatoryFlag, the function will return an error (flag > 0) or a warning (flag = 0) if the key is not found
      \return 0 if success, and positive value otherwise (1 if error, 2 if key not found).
    */
    template <typename T> int IntfKeyGetValueFromFile(const string& a_pathToHeader, const string& a_key, T* ap_return, int a_nbElts, int a_mandatoryFlag);

    /*!
      \fn IntfKeyGetValueFromFile(const string& a_pathToHeader, const string& a_key, T* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences)
      \param a_pathToHeader : path to the interfile header
      \param a_key : the key to recover
      \param T* ap_return : template array in which the data will be recovered
      \param int a_nbElts : number of elements to recover
      \param int a_mandatoryFlag : flag indicating if the data to recover if mandatory (true) or optionnal (false)
      \param int a_nbOccurences : number of occurences of the field before recovering the value
      \brief Look for "a_nbElts" elts in the "a_pathToHeader" interfile header matching the "a_keyword" key passed as parameter and return the corresponding value(s) in the "ap_return" templated array.\n
             Parameter "a_nbOccurences" can be set to recover a specific occurrence of a recurring value
      \details If more than one elements are to be recovered, the function first check the key has a correct Interfile kay layout (brackets and commas :  {,,})\n
               Depending on the mandatoryFlag, the function will return an error (flag > 0) or a warning (flag = 0) if the key is not found
      \return 0 if success, and positive value otherwise (1 if error, 2 if key not found).
    */
    template <typename T> int IntfKeyGetRecurringValueFromFile(const string& a_pathToHeader, const string& a_key, T* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);

    /*!
      \fn IntfReadProjectionImage(const string& a_pathToHeaderFile, FLTNB* ap_ImgMatrix, Intf_fields* ap_IF, int vb, bool a_lerpFlag)
      \param a_pathToHeaderFile : path to the header file
      \param ap_ImgMatrix : 1 dimensional image matrix which will recover the image.
      \param ap_IF : Intf_fields structure containing image metadata
      \param vb : Verbosity level
      \param a_lerpFlag : if true, enable linear interpolation of the image if img dimensions differ from the reconstruction dimensions
      \brief Main function which reads a projection Interfile 3D projection image and store its content in the provided ap_ImgMatrix \n
      \details Call the main functions dedicated to Interfile reading : IntfReadHeader(), and IntfReadImage()
      \return 0 if success, positive value otherwise.
    */
    int IntfReadProjectionImage( const string& a_pathToHeaderFile, 
                                FLTNB* ap_ImgMatrix,
                          Intf_fields* ap_IF, 
                                   int vb, 
                                  bool a_lerpFlag);

    /*!
      \fn      IntfCheckDimensionsConsistency
      \param   Intf_fields ap_ImgFields1: Intf_fields structure containing image metadata for first image
      \param   Intf_fields ap_ImgFields1: Intf_fields structure containing image metadata for second image
      \brief   This function checks that the dimensions of the two interfile image fields are the same
      \details It also checks the originating systems
      \return  0 if success, positive value otherwise
    */
    int IntfCheckDimensionsConsistency(Intf_fields ImgFields1, Intf_fields ImgFields2);

    /*!
      \fn      IntfLoadImageFromScratch
      \param   const string& a_pathToHeaderFile: path to the header image file
      \param   Intf_fields* ap_ImgFields: Intf_fields structure containing image metadata
      \param   int vb: Verbosity level
      \brief   Main function dedicated to load an Interfile 3D image from scratch and return the pointer to this image
      \details Call the main functions dedicated to Interfile reading : IntfReadHeader(), then IntfReadImage()
      \return  The pointer to the allocated memory containing the image data
    */
    FLTNB* IntfLoadImageFromScratch( const string& a_pathToHeaderFile, 
                                     Intf_fields* ap_ImgFields,
                                     int vb );
    /*!
      \fn      IntfWriteImageFromIntfFields
      \param   const string& a_pathToImg: output base name (including the path)
      \param   FLTNB* ap_ImgMatrix: the image to be written
      \param   Intf_fields Img_fields: Intf_fields structure containing image metadata
      \param   int vb: Verbosity level
      \brief   Main function dedicated to write a 3D interfile image from a FLTNB matrix and associated interfile fields
      \details Call the main functions dedicated to Interfile reading : IntfWriteHeaderMainData(), then IntfWriteImage()
      \return  0 upon success, positive value otherwise
    */
    int IntfWriteImageFromIntfFields( const string& a_pathToImg,
                                      FLTNB* ap_ImgMatrix,
                                      Intf_fields Img_fields,
                                      int vb );

    /*!
      \fn IntfReadImage(const string& a_pathToHeaderFile, FLTNB* ap_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb, bool a_lerpFlag)
      \param a_pathToHeaderFile : path to the header file
      \param ap_ImgMatrix : 1 dimensional image matrix which will recover the image.
      \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
      \param vb : Verbosity level
      \param a_lerpFlag : if true, enable linear interpolation of the image if img dimensions differ from the reconstruction dimensions
      \brief Main function dedicated to Interfile 3D image loading
      \details Call the main functions dedicated to Interfile reading : IntfReadHeader(), IntfCheckConsistency(), then IntfReadImage()
      \return 0 if success, positive value otherwise.
    */
    int IntfReadImage(   const string& a_pathToHeaderFile, 
                             FLTNB* ap_ImgMatrix,
    oImageDimensionsAndQuantification* ap_ID, 
                                   int vb, 
                                  bool a_lerpFlag);
                                  
    /*!
      \fn IntfReadImage(const string& a_pathToHeaderFile, FLTNB**** a4p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb, bool a_lerpFlag)
      \param a_pathToHeaderFile : path to the main header file
      \param a4p_ImgMatrix : 4 dimensional image matrix which will recover the image.
      \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
      \param vb : Verbosity level
      \param a_lerpFlag : if true, enable linear interpolation of the image if img dimensions differ from the reconstruction dimensions
      \brief Main function dedicated to Interfile 5D (1D+1D+1D time + 3D) image loading
      \details Check is the main header file is a metaheader associated to several image files, or a unique interfile header\n
               Depending on the type of file input (metaheader or unique file), read the group of image files or the unique provided image file
      \todo : Check if image 3D dimensions are different from an image to another ?
              (very unlikely, but it would cause segfault if interpolation is enabled)
      \return 0 if success, positive value otherwise.
    */
    int IntfReadImage(   const string& a_pathToHeaderFile, 
                          FLTNB**** a4p_ImgMatrix, 
    oImageDimensionsAndQuantification* ap_ID, 
                                   int vb, 
                                  bool a_lerpFlag);
    
    /*!
      \fn IntfReadImgDynCoeffFile(const string& a_pathToHeaderFile, FLTNB** a2p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int a_nbFbases, int vb, bool a_lerpFlag)
      \param a_pathToHeaderFile : path to the header file
      \param a2p_ImgMatrix : 2 dimensional image matrix which will recover the image..
      \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
      \param a_nbFbases : Number of basis functions
      \param vb : verbosity
      \param a_lerpFlag : if true, enable linear interpolation of the image if img dimensions differ from the reconstruction dimensions
      \brief Function dedicated to Interfile image reading for dynamic coefficients images
      \details The total number of basis functions should be provided in parameters\n
               Check is the main header file is a metaheader associated to several image files, or a unique interfile header\n
               Depending on the type of file input (metaheader or unique file), read the group of image files or the unique provided image file
      \return 0 if success, positive value otherwise.
    */
    int IntfReadImgDynCoeffFile(const string& a_pathToHeaderFile, 
                                   FLTNB** a2p_ImgMatrix, 
           oImageDimensionsAndQuantification* ap_ID,
                                          int a_nbFbases, 
                                          int vb, 
                                         bool a_lerpFlag);
  
  
  // -------------------------------------------------------------------
  // ----- FUNCTIONS DEDICATED TO IMAGE DATA READING/WRITING -----

    /*!
      \fn int IntfCheckConsistency(Intf_fields* ap_IF, oImageDimensionsAndQuantification* ap_ID, int vb, int a_lerpFlag)
      \param Intf_fields* ap_IF : Structure containing the interfile fields read in a interfile header
      \param oImageDimensionsAndQuantification* ap_ImageDimensions : Provide the Image dimensions object containing reconstruction dimensions
      \param int vb : Verbosity level
      \param a_lerpFlag : if true, enable linear interpolation of the image if img dimensions differ from the reconstruction dimensions
      \brief Check if the mandatory fields have been initialize in the ap_IF structure, and check consistencies with the reconstruction parameters
      \details This function also checks if the matrix size of the original image is different to the reconstruction matrix size.\n
               In this case a boolean is set up to enable interpolation during image reading
      \todo Add check for all mandatory parameters, and temporal image dimensions
      \todo Float comparison ?
      \return 0 if success, positive value otherwise.
    */
    int IntfCheckConsistency(Intf_fields* ap_IF, oImageDimensionsAndQuantification* ap_ID, int vb, int a_lerpFlag);

    /*!
      \fn IntfGetPixelTypeAndReadData(Intf_fields a_IF, ifstream* ap_iFile, FLTNB* ap_outImgMatrix, FLTNB* ap_inImgMatrix, uint32_t* a_offset, int a_nbVox, int vb)
      \param a_IF : Interfile fields recovered from the header
      \param ap_iFile : Ifstream pointing to an image file
      \param ap_outImgMtx : 3D image matrix with reconstruction dimensions/voxel size
      \param ap_inImgMtx : 3D image matrix with original dimensions/voxel size
      \param ap_offset : Offset indicating the beginning of the data to read in the image file
      \param a_nbVox : A number of voxels in the 3D image matrix with reconstruction dimensions/voxel size
      \param vb : Verbosity level
      \brief The purpose of this function is to call the templated ReadData() function with the data type corresponding to the interfile image 
      \details It uses "number format" and "number of bytes per pixel" fields to identify the correct type\n
               ASCII and bit images NOT supported
      \return 0 if success, positive value otherwise.
    */
    int IntfGetPixelTypeAndReadData(Intf_fields a_IF, 
                                      ifstream* ap_iFile, 
                                      FLTNB* ap_outImgMatrix, 
                                      FLTNB* ap_inImgMatrix, 
                                      uint32_t* a_offset, 
                                            int a_nbVox, 
                                            int vb);

    /*!
      \fn IntfReadData(Intf_fields a_IF, ifstream* ap_iFile, FLTNB* ap_outImgMatrix, FLTNB* ap_inImgMatrix, uint32_t* a_offset, int a_nbVox, int vb, T* bytes)
      \param a_IF : Interfile fields recovered from the header
      \param ap_iFile : Ifstream pointing to an image file
      \param ap_outImgMtx : 3D image matrix with reconstruction dimensions/voxel size
      \param ap_inImgMtx : 3D image matrix with original dimensions/voxel size
      \param ap_offset : Offset indicating the beginning of the data to read in the image file
      \param a_nbVox : A number of voxels in the 3D image matrix with reconstruction dimensions/voxel size
      \param vb : Verbosity level
      \param T* bytes : Buffer of templated size, to recover any voxel value
      \brief Templated function which read an image voxel by voxel and store it in the ap_outImgMtx image matrix
      \details Call an interpolation function if the original image dimensions/voxel sizes are different to the reconstruction dimensions/voxel sizes\n
               Manage endianness and the optionnal calibration with the rescale slope/intercept interfile keys
      \todo Re-orient image data from the original orientation to the default orientation (transaxial/supine/headin)
      \todo Perhaps check the conversion from original image type to FLTNB type has been correctly performed (maybe time-consuming)
      \return 0 if success, positive value otherwise.
    */
    template <class T>
    int IntfReadData(Intf_fields a_IF, 
                       ifstream* ap_iFile, 
                       FLTNB* ap_outImgMatrix, 
                       FLTNB* ap_inImgMatrix, 
                       uint32_t* a_offset, 
                             int a_nbVox, 
                             int vb, 
                              T* bytes);


    /*!
      \fn ImageInterpolation
      \param ap_iImg : Image matrix to interpolate, with dimensions of the original interfile
      \param ap_oImg : Image matrix to recover the interpolated image to the reconstruction dimensions 
      \param ap_iDimVox[3] : X,Y,Z dimensions of the image to interpolate
      \param ap_oDimVox[3] : X,Y,Z dimensions for the reconstruction
      \param ap_iSizeVox[3] : X,Y,Z voxel size of the image to interpolate
      \param ap_oSizeVox[3] : X,Y,Z voxel size for the reconstruction
      \param ap_iOffVox[3] : X,Y,Z offset position of the image to interpolate (mm)
      \param ap_oOffVox[3] : X,Y,Z offset position for the reconstruction (mm)
      \brief Trilinear interpolation
      \return 0 if success, positive value otherwise.
    */
    int ImageInterpolation(FLTNB *ap_iImg,               FLTNB *ap_oImg, 
                     const uint32_t ap_iDimVox[3], const uint32_t ap_oDimVox[3],
                     const FLTNB ap_iSizeVox[3],   const FLTNB ap_oSizeVox[3],
                     const FLTNB ap_iOffVox[3],    const FLTNB ap_oOffVox[3] );
          

    /*!
      \fn IntfWriteImgFile(const string& a_pathToImg, FLTNB* ap_ImgMatrix, const Intf_fields& ap_IF, int vb)
      \param a_pathToImg : path to image basename
      \param ap_ImgMatrix : 1 dimensional image matrix which contains the image to write
      \param ap_IntfF : Intf_fields structure containing image metadata
      \param vb : verbosity
      \brief Main function dedicated to Interfile 3D image writing. \n
             Recover image information from a provided Intf_fields structure.
      \details Call the main functions dedicated to Interfile header and data writing :
               IntfWriteHeaderMainData() and then IntfWriteImage()
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteImgFile(const string& a_pathToImg, FLTNB* ap_ImgMatrix, const Intf_fields& ap_IntfF, int vb);

    /*!
      \fn IntfWriteImgFile(const string& a_pathToImg, FLTNB* ap_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb)
      \param a_pathToImg : path to image basename
      \param ap_ImgMatrix : 1 dimensional image matrix which contains the image to write
      \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
      \param vb : verbosity
      \brief Main function dedicated to Interfile 3D image writing
      \details Call the main functions dedicated to Interfile header and data writing : \n
               IntfWriteHeaderMainData() and then IntfWriteImage()
      \todo Get metadata from a Intf_fields object ? \n
           (useful to transfer keys from read images to written images)
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteImgFile(const string& a_pathToImg, FLTNB* ap_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb);

    /*!
      \fn IntfWriteProjFile(const string& a_pathToImg, FLTNB** a2p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, Intf_fields a_Imgfields, int vb)
      \param a_pathToImg : string containing the path to the image basename
      \param a2p_ImgMatrix : 2 dimensional image matrix which contains the image to write.
      \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
      \param a_Imgfields: Structure containing information about the projected data
      \param vb : verbosity
      \brief Function dedicated to Interfile image writing for projected data
      \details Call the main functions dedicated to Interfile header and data writing \n
               Currently work for SPECT projected data \n
               The total number of projections should be provided in parameters \n
               Depending on the output writing mode (stored in sOutputManager), \n
      \todo Get metadata from a Intf_fields object ? \n
           (useful to transfer keys from read images to written images)
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteProjFile(const string& a_pathToImg, FLTNB** a2p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, Intf_fields a_Imgfields, int vb);

    /*!
      \fn IntfWriteImgDynCoeffFile(const string& a_pathToImg, FLTNB** a2p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int a_nbFbases, int vb)
      \param a_pathToImg : string containing the path to the image basename
      \param a2p_ImgMatrix : 2 dimensional image matrix which contains the image to write.
      \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
      \param a_nbParImgs : Number of parametric images
      \param vb : verbosity
      \brief Function dedicated to Interfile image writing for dynamic coefficients images
      \details Call the main functions dedicated to Interfile header and data writing \n
               The total number of basis functions should be provided in parameters \n
               Depending on the output writing mode (stored in sOutputManager), \n
               One image with one file and one will be created for the whole dynamic image \n
               or a metaheader and associated multiple 3D image raw file/headers will be generated
      \todo Get metadata from a Intf_fields object ? \n
           (useful to transfer keys from read images to written images)
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteImgDynCoeffFile(const string& a_pathToImg, FLTNB** a2p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int a_nbParImgs, int vb);

    /*!
      \fn IntfWriteImgFile(const string& a_pathToImg, FLTNB**** a4p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb)
      \param a_pathToImg : string containing the path to the image basename
      \param a4p_ImgMatrix : 4 dimensional image matrix which contains the image to write.
      \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
      \param vb : verbosity
      \brief Main function dedicated to Interfile 6D (dynamic dims + 3D ) image writing
      \details Call the main functions dedicated to Interfile header and data writing \n
               Depending on the output writing mode (stored in sOutputManager), \n
               One image with one file and one will be created for the whole dynamic image \n
               or a metaheader and associated multiple 3D image raw file/headers will be generated \n
      \todo Get metadata from a Intf_fields object ? \n
           (useful to transfer keys from read images to written images)
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteImgFile(const string& a_pathToImg, FLTNB**** a4p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb);

    /*!
      \fn IntfWriteImage(const string& a_pathToImg, FLTNB* ap_outImgMtx, uint32_t a_dim, int vb)
      \param a_pathToImg : th to the directory where the image will be written
      \param ap_outImgMtx : Matrix containing the image to write
      \param a_dim : Number of voxels in the 3D image
      \param vb : Verbosity level
      \brief Write Interfile raw data whose path is provided in parameter, using image matrix provided in parameter.
      \brief For 1 dimensional image matrices
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteImage(const string& a_pathToImg, FLTNB* ap_outImgMtx, uint32_t a_dim, int vb);

    /*!
      \fn IntfWriteImage(vector<string> ap_pathToImgs, FLTNB** a2p_outImgMtx, uint32_t ap_dim[2], int vb)
      \param ap_pathToImgs: List of string containing the paths to the image to write
      \param a2p_outImgMtx : 2 dimensional image matrix (1D temporal + 3D)
      \param ap_imgDim[2] : Dimensions of ap_outImgMtx
      \param vb : Verbosity level
      \brief Write Interfile raw data whose path is provided in parameter, using the image matrix provided in parameter.
      \brief For 2 dimensional image matrices
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteImage(vector<string> ap_pathToImgs, FLTNB** a2p_outImgMtx, uint32_t ap_dim[2], int vb);

    /*!
      \fn IntfWriteImage(vector<string> ap_pathToImgs, FLTNB**** a4p_outImgMtx, uint32_t ap_dim[4], int vb)
      \param vector<string> ap_pathToImgs: List of string containing the paths to the image to write
      \param FLTNB**** a4p_outImgMtx : 4 dimensional image matrix (3D temporal + 3D)
      \param int ap_imgDim[4] : Dimensions of ap_outImgMtx
      \param int vb : Verbosity level
      \brief Write Interfile raw data whose path is provided in parameter, using the image matrix provided in parameter.
      \brief For 4 dimensional image matrices
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteImage(vector<string> ap_pathToImgs, FLTNB**** a4p_outImgMtx, uint32_t ap_dim[4], int vb);

    /*!
      \fn IntfWriteData(ofstream* ap_oFile, FLTNB* ap_outImgMatrix, int a_nbVox, int vb)
      \param ap_oFile : Ofstream pointing to an image file
      \param ap_outImgMtx : 3D image matrix with reconstruction dimensions/voxel size containing the image data
      \param a_nbVox : A number of voxels in the 3D image matrix with reconstruction dimensions/voxel size
      \param vb : Verbosity level
      \brief Write the content of the image matrix in the file pointed by ofstream
      \todo  keep original orientations ? Would require a loop on voxels and a reindexing of voxel indexes
      \todo  check writing error
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteData(ofstream* ap_oFile, FLTNB* ap_outImgMatrix, int a_nbVox, int vb);


  // -------------------------------------------------------------------
  // ----- FUNCTIONS DEDICATED TO INTERFILE HEADER KEYs DECODING/PARSING/READING -----

    /*!
      \fn IntfIsMHD(string a_pathToFile, vector<string> &ap_lPathToImgs)
      \param a_pathToFile : String containing path to an Interfile header
      \param ap_lPathToImgs : pointer to a list of strings containing the path to the different images
      \brief Check if the string in argument contains the path to a Interfile metaheader
      \details Check for '!total number of data sets' key, and the associated image file names \n
               If the file is metaheader, the names of the header files of the images will be returned in ap_lPathToImgs list \n
               It not, the file is a single header, that we add to the list as solo file
      \return 1 if we have a metaheader, \n
              0 if not, \n
              negative value if error.
    */
    int IntfIsMHD(string a_pathToFile, vector<string> &ap_lPathToImgs);

    /*!
      \fn IntfWriteMHD(const string& a_pathToMhd, const vector<string> &ap_lPathToImgs, Intf_fields a_IntfF, oImageDimensionsAndQuantification* ap_ID, int vb)
      \param a_path : path to the Meta interfile header to write
      \param ap_lPathToImgs : pointer to a list of strings containing the path to the different images
      \param a_IntfF : Structure containing interfile keys of the image to write
      \param ap_ID : oImageDimensionsAndQuantification object containing additional infos about reconstruction
      \param vb : verbosity
      \brief Write an Interfile meta header at the path provided in parameter, using the field stack provided in parameter.
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteMHD(const string& a_pathToMhd, const vector<string> &ap_lPathToImgs, Intf_fields a_IntfF, oImageDimensionsAndQuantification* ap_ID, int vb);

    /*!
      \fn IntfReadHeader(const string& a_pathToHeaderFile, Intf_fields* ap_IntfFields, int vb)
      \param const string& a_pathToHeaderFile
      \param Intf_fields* ap_IF
      \param int vb : Verbosity level
      \brief Read an Interfile header
      \details Initialize all mandatory fields from the Intf_fields structure passed in parameter with the related interfile key from the image interfile header
      \todo Check everything work for Interfile 3.3 version, extended interfile, etc.
      \todo Several keys from Intf_fields still needs some management
      \todo Deal with sinogram in interfile format (matrix size keys)
      \todo orientation related keys
      \todo time windows, energy windows keys
      \todo check start angle key
      \todo Initializations & checks at the end (ok for slice_thickness, check slice_spacing, etc..)
      \todo Specification of this function for projection data (currently we consider the interfile is a 3D image)
      \return 0 if success, positive value otherwise.
    */
    int IntfReadHeader(const string& a_pathToHeaderFile, Intf_fields* ap_IntfFields, int vb);

    /*!
      \fn IntfWriteHeaderMainData(const string& a_path, Intf_fields a_IntfF, oImageDimensionsAndQuantification* ap_ID, int awmode, int vb)
      \param a_path : path to the interfile header to write
      \param ap_IntfF : Reference to a Intf_fields structure containing interfile keys of the image to write
      \param int awmode : writing mode for dynamic images (one merged file or separate files for each frame/gate)
      \param vb : verbosity
      \brief Write main infos of an Interfile header at the path provided in parameter, \n
             using the interfile keys provided by the field structure parameter.
      \todo Still a lot of fields to handle, and awmode
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteHeaderMainData(const string& a_path, const Intf_fields& ap_IntfF, int vb);

    /*!
      \fn IntfWriteHeaderImgData(ofstream &ap_ofile, Intf_fields a_IntfF, oImageDimensionsAndQuantification* ap_ID, int vb)
      \param ap_ofile : pointer to the ofstream linked to the image header file 
      \param ap_IntfF : Reference to a Intf_fields structure containing interfile keys of the image to write
      \param verbosity : verbosity
      \brief Write basic image data info in the file provided in parameter
      \todo Data about positionning/offsets ?
      \todo include dectector heads (a set of fields for each head)
      \return 0 if success, positive value otherwise.
    */
    int IntfWriteHeaderImgData(ofstream &ap_ofile, const Intf_fields& ap_IntfF, int vb);

    /*!
      \fn void IntfKeyInitFields(Intf_fields* ap_IF)
      \param ap_IF : Structure containing Interfile keys
      \brief Init the file of an Interfile fields structure passed in parameter to their default values
    */
    void IntfKeyInitFields(Intf_fields* ap_IF);
    
    /*!
      \fn IntfKeySetFieldsOutput(Intf_fields* ap_IF, oImageDimensionsAndQuantification* ap_ID)
      \param ap_IF : Structure containing Interfile keys
      \param ap_ID : oImageDimensionsAndQuantification object containing additional infos about reconstruction
      \brief Init the keys of the Interfile header of an image to be written on disk
      \details Init the keys of the Interfile structure passed in parameter for output writing
               using the ImageDimensions object containing information about the reconstruction
    */
    void IntfKeySetFieldsOutput(Intf_fields* ap_IF, oImageDimensionsAndQuantification* ap_ID);

    /*!
      \fn IntfKeyPrintFields(Intf_fields a_IF)
      \param ap_IF
      \brief Print all the keys of the Intf_fields structure passed in parameter, as well as their values for debugging purposes
    */
    void IntfKeyPrintFields(Intf_fields a_IF);

    /*!
      \fn IntfRecoverKey(Intf_key* ap_Key, const string& a_line)
      \param ap_Key : Structure to recover the parsed key components (key, value,..)
      \param a_line : String to process
      \brief Process the line passed in parameter and write the key information in the ap_Key Intf_key member structure
      \details .korig : Get original line without comments \n
               .kcase : Get key without spaces and without comments \n
               .klcase: Same as kcase, in lower case \n
               .kvalue: Value of the key, without spaces
      \todo Check that IntfToLowerCase() doesn't cause issue with some characters or some ASCII file format (unicode, etc..)
      \return 0 if success, positive value otherwise.
    */
    int IntfRecoverKey(Intf_key* ap_Key, const string& a_line);

    /*!
      \fn IntfCheckKeyMatch(Intf_key ap_Key, const string& a_field)
      \param ap_Key : Structure containing the parsed key components (key, value,..)
      \param a_line : String containing an interfile key
      \brief Check if the key matches the string passed in parameter
      \todo Be sure it is appropriate to use Intf_key.klcase for each key comparison
      \return 1 if success, 0 otherwise (not found).
    */
    int IntfCheckKeyMatch(Intf_key a_Key, const string& a_field);

    /*!
      \fn int IntfKeyIsArray(Intf_key ap_Key)
      \param ap_Key
      \brief Check if the key passed in parameter is an array (contains brackets '{' and '}' )
      \return 1 if success, 0 otherwise (not array).
    */
    int IntfKeyIsArray(Intf_key ap_Key);

    /*!
      \fn IntfKeyGetArrayNbElts(Intf_key ap_Key)
      \param ap_Key
      \brief Return the number of elts in an Interfile array Key
      \return the number of elements in the array key, or negative value if error
    */
    int IntfKeyGetArrayNbElts(Intf_key ap_Key);    

    /*!
      \fn IntfKeyGetMaxArrayKey(Intf_key ap_Key)
      \param ap_Key
      \brief Return the maximum value from an array key (key value contains brackets '{,,}' )
      \return the max value in the array key.
    */
    int IntfKeyGetMaxArrayKey(Intf_key ap_Key);

    /*!
      \fn IntfKeyGetArrayElts(Intf_key a_Key, T* ap_return)
      \param ap_Key
      \param T* ap_return : Templated parameter in which the elts will be returned
      \brief Get all the elements in an array key in a templated array passed in parameter. \n
             It assumes the return variable has been instanciated with a correct number of elements.
      \return 0 if success, positive value otherwise
    */
    template <typename T> int IntfKeyGetArrayElts(Intf_key a_Key, T* ap_return);

    /*!
      \fn IntfKeyGetEndianStr(int a_val)
      \param a_val : 
      \brief return the endian string corresponding to the value passed in parameter (see module INTF_ENDIANNESS).
      \return "BIG_ENDIAN" if 0, \n
              "LITTLE_ENDIAN" if 1, \n
              "UNKNOWN" otherwise
    */
    string IntfKeyGetEndianStr(int a_val);

    /*!
      \fn IntfKeyGetModalityStr(int a_modalityIdx)
      \param a_modalityIdx
      \brief Convert the integer provided in parameter to the string related \n
             to the corresponding modality as defined by the scanner objects
      \todo Add other modalities as we implement them
      \return string corresponding to the modality
    */
    string IntfKeyGetModalityStr(int a_modalityIdx);

    /*!
      \fn IntfKeyGetInputImgDataType(const string& a_str)
      \param a_str : input string
      \brief Get the image data type corresponding to the image metadata passed in parameter
      \return int value corresponding to the image data type (see module INTF_IMG_TYPE).
    */
    int IntfKeyGetInputImgDataType(const string& a_str);

    /*!
      \fn IntfKeyGetOutImgDataType(oImageDimensionsAndQuantification* ap_ID)
      \param ap_ID
      \brief Get the image data type corresponding to the image metadata passed in parameter
      \return int value corresponding to the image data type (see module INTF_IMG_TYPE).
    */
    int IntfKeyGetOutputImgDataType(oImageDimensionsAndQuantification* ap_ID);

    /*!
      \fn IntfKeyGetPixTypeStr()
      \brief Return the string corresponding to the nb of bytes in the type FLTNB
      \return string corresponding to the pixel data type
    */
    string IntfKeyGetPixTypeStr();


  // -------------------------------------------------------------------
  // ----- NOT IMPLEMENTED KEY FUNCTION -----

    /*!
      \fn IntfKeyGetPatOrientation(Intf_fields ap_IF)
      \param ap_IF
      \brief Get the complete patient orientation from an Intf_fields structure according to 
             the values of keys 'slice orientation', 'patient rotation', and 'patient orientation'
      \return int value corresponding to the patient orientation (see module PATIENT_ORIENTATION).
      \todo NOT CURRENTLY IMPLEMENTED
    */
    int IntfKeyGetPatOrientation(Intf_fields ap_IF);
    
    /*!
      \fn IntfGetVoxIdxSHTOrientation(Intf_fields a_IF, int a_voxId)
      \param ap_IF
      \param a_voxId : index of the voxel in a 1-D image vector
      \brief Compute a voxel index corresponding to the default orientation (Sup/Hin/Trans) \n
             from the orientation informations contained in the Intf_fields passed in parameter
      \todo NOT CURRENTLY IMPLEMENTED
      \return new voxel index
    */  
    int IntfGetVoxIdxSHTOrientation(Intf_fields a_IF, int a_voxId);


  // -------------------------------------------------------------------
  // ----- SOME UTILITY FUNCTIONS -----

    /*!
      \fn IntfAllocInterpImg(FLTNB** a2p_img, Intf_fields a_IF)
      \param ap_img : pointer to 1 dimensional image matrix to recover the image to interpolate
      \param ap_IF : Structure containing the interfile fields read in a interfile header
      \brief Allocate memory for an image matrix to recover an image to interpolate
    */
    void IntfAllocInterpImg(FLTNB** a2p_img, Intf_fields a_IF);
    
    /*!
      \fn GetUserEndianness()
      \brief Check user/host computer endianness and write it to the global variable User_Endianness
      \details This function should be called once during the initialization step of the algorithm (currently, the singleton initialization)
      \todo Maybe better compute it in a preprocessor macro
    */
    void GetUserEndianness();

    /*!
      \fn IntfEraseSpaces(string* input_str)
      \param string* input_str
      \brief Erase space, blank characters (\(t,r,n)), and '!' before and after the characters in the string passed in parameter
    */
    void IntfEraseSpaces(string* input_str);
    
    /*!
      \fn void IntfToLowerCase(string* ap_str)
      \param string* ap_str : original string
      \brief Set all characters of the string passed in parameter to lower case
      \todo May be issues with non ASCII characters, file decoding, etc..
    */
    void IntfToLowerCase(string* ap_str);

    /*!
      \fn void IntfKeyGetPatientNameTag()
      \brief Recover datafile name(s) stored in sOutputManager in one string
    */
    string IntfKeyGetPatientNameTag();

    /*!
      \fn SwapBytes(Type *ap_type)
      \param T *ap_type : Variable of type T
      \brief Use std::reverse to swap the bits of a variable of any type
      \details Used for Little/Big Endian conversion
    */
    template <class Type> void SwapBytes(Type *ap_type);

#endif
