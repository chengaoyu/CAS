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
  \brief    Implementation of class Interfile management functions
*/

#include "oInterfileIO.hh"
#include <iomanip>

#ifdef _WIN32
// Avoid compilation errors due to mix up between std::min()/max() and
// min max macros
#undef min
#undef max
#endif

int User_Endianness = -1;  /*!< Global variable recovering endianness of user system */

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyGetValueFromFile
  \param a_pathToHeader : path to the interfile header
  \param a_key : the key to recover
  \param T* ap_return : template array in which the data will be recovered
  \param int a_nbElts : number of elements to recover
  \param int a_mandatoryFlag : flag indicating if the data to recover if mandatory (true) or optionnal (false)
  \brief Look for "a_nbElts" elts in the "a_pathToHeader" interfile header matching the "a_keyword" key 
         passed as parameter and return the corresponding value(s) in the "ap_return" templated array.
  \details If more than one elements are to be recovered, the function first check that
           the key has a correct Interfile key layout (brackets and commas :  {,,})
           Depending on the mandatoryFlag, the function will return an error (flag > 0)
           or a warning (flag = 0) if the key is not found
  \return 0 if success, and positive value otherwise (1 if error, 2 if key not found).
*/
template<class T>
int IntfKeyGetValueFromFile(const string& a_pathToHeader, 
                            const string& a_key, 
                                       T* ap_return, 
                                      int a_nbElts, 
                                      int a_mandatoryFlag)
{
  ifstream input_file(a_pathToHeader.c_str(), ios::in);
  string line;
  
  // Check file
  if (input_file)
  {
    while(!input_file.eof())
    {
      getline(input_file, line);
  
      if(line.empty())
      continue;
      
      Intf_key Key;
      
      // Process the Key at this line
      if (IntfRecoverKey(&Key, line) ) 
      {
        Cerr("***** IntfKeyGetValueFromFile()-> Error : couldn't correctly read interfile key '"<< line << "' !" << endl);
        return 1;
      }
      //if (line.find(a_keyword) != string::npos)
      if(IntfCheckKeyMatch(Key, a_key)) 
      { 
        
        if(a_nbElts == 1) // 1 elt is required, just return the result
        {
          if (ConvertFromString(Key.kvalue, &ap_return[0]))
          {
            Cerr("***** IntfKeyGetValueFromFile()-> Exception when trying to read tag '" << a_key << "' in file '" << a_pathToHeader << "'." << endl);
            return 1;
          }
          
          return 0;
        }
        else // Several elements required
        {
          // First check we have an array
          if (!IntfKeyIsArray(Key))
          {
            Cerr("***** IntfKeyGetValueFromFile() -> " << a_nbElts << " requested for interfile key " << a_key << " , but this key is not an array !" << endl);
            return 1;
          }
          else
          {
            // Check the number of elements in the key.
            if(IntfKeyGetArrayNbElts(Key) != a_nbElts)
            {
              Cerr("***** IntfKeyGetValueFromFile() -> Nb of elements to recover (=" << a_nbElts << ") does not correspond to the number of elements found in the key '" 
                                                                                                 << a_key << "' (" << IntfKeyGetArrayNbElts(Key) << ") !" << endl);
              return 1;
            }
            
            // Read array key
            if (IntfKeyGetArrayElts(Key, ap_return) )
            {
              Cerr("***** IntfKeyGetValueFromFile() -> " << a_nbElts << " requested for interfile key " << a_key << " , but this key is not an array !" << endl);
              return 1;
            }
            
            return 0;
          }
        }
      }
    }
    
    // Tag not found, throw an error message if the tag is mandatory
    if (a_mandatoryFlag > 0) 
    {
      Cerr("***** IntfKeyGetValueFromFile()-> Error when reading Interfile '" << a_pathToHeader << "'. Key '" << a_key << "' was not found." << endl);
      return KEYWORD_MANDATORY_NOT_FOUND;
    }
    else
    {
      return KEYWORD_OPTIONAL_NOT_FOUND;
    }
    
  }
  else
  {
    Cerr("***** IntfKeyGetValueFromFile()-> Couldn't find or read data-file '"<< a_pathToHeader << "' !" << endl);
    return 1;
  }
}

// Templated functions definitions
template int IntfKeyGetValueFromFile<string>(const string& a_pathToHeader, const string& a_key, string* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<int>(const string& a_pathToHeader, const string& a_key, int* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<int64_t>(const string& a_pathToHeader, const string& a_key, int64_t* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<float>(const string& a_pathToHeader, const string& a_key, float* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<double>(const string& a_pathToHeader, const string& a_key, double* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<long double>(const string& a_pathToHeader, const string& a_key, long double* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<uint8_t>(const string& a_pathToHeader, const string& a_key, uint8_t* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<uint16_t>(const string& a_pathToHeader, const string& a_key, uint16_t* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<uint32_t>(const string& a_pathToHeader, const string& a_key, uint32_t* ap_return, int a_nbElts, int a_mandatoryFlag);
template int IntfKeyGetValueFromFile<bool>(const string& a_pathToHeader, const string& a_key, bool* ap_return, int a_nbElts, int a_mandatoryFlag);


/*
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
template <typename T> int IntfKeyGetRecurringValueFromFile(const string& a_pathToHeader, 
                                                           const string& a_key, 
                                                              T* ap_return, 
                                                             int a_nbElts, 
                                                             int a_mandatoryFlag, 
                                                        uint16_t a_nbOccurrences)
{
  ifstream input_file(a_pathToHeader.c_str(), ios::in);
  string line;
  uint16_t nb_occurences_cur =0;
  
  // Check file
  if (input_file)
  {
    while(!input_file.eof())
    {
      getline(input_file, line);
  
      if(line.empty())
      continue;
      
      Intf_key Key;
      
      // Process the Key at this line
      if (IntfRecoverKey(&Key, line) ) 
      {
        Cerr("***** IntfKeyGetValueFromFile()-> Error : couldn't correctly read interfile key '"<< line << "' !" << endl);
        return 1;
      }


      
      //if (line.find(a_keyword) != string::npos)
      if(IntfCheckKeyMatch(Key, a_key)) 
      { 
        // Check if we reached the correct number of occurence of the key
        // Skip otherwise
        if(nb_occurences_cur < a_nbOccurrences)
        {
          nb_occurences_cur++;
          continue;
        }
        
        if(a_nbElts == 1) // 1 elt is required, just return the result
        {
          if (ConvertFromString(Key.kvalue, &ap_return[0]))
          {
            Cerr("***** IntfKeyGetValueFromFile()-> Exception when trying to read tag '" << a_key << "' in file '" << a_pathToHeader << "'." << endl);
            return 1;
          }
          
          return 0;
        }
        else // Several elements required
        {
          // First check we have an array
          if (!IntfKeyIsArray(Key))
          {
            Cerr("***** IntfKeyGetValueFromFile() -> " << a_nbElts << " requested for interfile key " << a_key << " , but this key is not an array !" << endl);
            return 1;
          }
          else
          {
            // Check the number of elements in the key.
            if(IntfKeyGetArrayNbElts(Key) != a_nbElts)
            {
              Cerr("***** IntfKeyGetValueFromFile() -> Nb of elements to recover (=" << a_nbElts << ") does not correspond to the number of elements found in the key '" 
                                                                                                 << a_key << "' (" << IntfKeyGetArrayNbElts(Key) << ") !" << endl);
              return 1;
            }
            
            // Read array key
            if (IntfKeyGetArrayElts(Key, ap_return) )
            {
              Cerr("***** IntfKeyGetValueFromFile() -> " << a_nbElts << " requested for interfile key " << a_key << " , but this key is not an array !" << endl);
              return 1;
            }
            
            return 0;
          }
        }
      }
    }
    
    // Tag not found, throw an error message if the tag is mandatory
    if (a_mandatoryFlag > 0) 
    {
      Cerr("***** IntfKeyGetValueFromFile()-> Error when reading Interfile '" << a_pathToHeader << "'. Key '" << a_key << "' was not found." << endl);
      return KEYWORD_MANDATORY_NOT_FOUND;
    }
    else
    {
      return KEYWORD_OPTIONAL_NOT_FOUND;
    }
    
  }
  else
  {
    Cerr("***** IntfKeyGetValueFromFile()-> Couldn't find or read data-file '"<< a_pathToHeader << "' !" << endl);
    return 1;
  }
}

// Templated functions definitions
template int IntfKeyGetRecurringValueFromFile<string>(const string& a_pathToHeader, const string& a_key, string* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<int>(const string& a_pathToHeader, const string& a_key, int* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<int64_t>(const string& a_pathToHeader, const string& a_key, int64_t* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<float>(const string& a_pathToHeader, const string& a_key, float* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<double>(const string& a_pathToHeader, const string& a_key, double* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<long double>(const string& a_pathToHeader, const string& a_key, long double* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<uint8_t>(const string& a_pathToHeader, const string& a_key, uint8_t* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<uint16_t>(const string& a_pathToHeader, const string& a_key, uint16_t* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<uint32_t>(const string& a_pathToHeader, const string& a_key, uint32_t* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);
template int IntfKeyGetRecurringValueFromFile<bool>(const string& a_pathToHeader, const string& a_key, bool* ap_return, int a_nbElts, int a_mandatoryFlag, uint16_t a_nbOccurrences);




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int IntfReadProjectionImage( const string& a_pathToHeaderFile, 
                            FLTNB* ap_ImgMatrix,
                      Intf_fields* ap_IF, 
                               int vb, 
                              bool a_lerpFlag)
{    
  if(vb >= 3) Cout("IntfReadProjectionImage()-> Read Interfile header : "<< a_pathToHeaderFile << endl);

  // Init Interfile Key structure
  IntfKeyInitFields(ap_IF);

  // Recover image infos from the header file
  if(IntfReadHeader(a_pathToHeaderFile, ap_IF, vb) )
  {
    Cerr("***** IntfReadProjectionImage()-> Error : while trying to read the interfile header '"<< a_pathToHeaderFile << "' !" << endl);
    return 1;
  }

  // Specific changes for projection data. //TODO : deal with that in IntfReadHeader, specification regarding the nature of the data
  ap_IF->mtx_size[2] = ap_IF->nb_total_imgs;

  int nb_tot_pixels = ap_IF->mtx_size[0]
                     *ap_IF->mtx_size[1]
                     *ap_IF->mtx_size[2];
    
  // Read image data
  ifstream img_file(ap_IF->path_to_image.c_str(), ios::binary | ios::in);
  if(img_file)
  {
    // Get the position of the beginning of the image data
    uint32_t offset = ap_IF->data_offset;

    // Call the right IntfReadData() function according to the pixel type, and read data
    IntfGetPixelTypeAndReadData(*ap_IF, &img_file, ap_ImgMatrix, NULL, &offset, nb_tot_pixels, vb);
  }
  else
  {
    Cerr("***** IntfReadProjectionImage()-> Error occurred while trying to read the image file at the path:  '"<< ap_IF->path_to_image << "' !" << endl);
    return 1;
  }
                    
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int IntfCheckDimensionsConsistency(Intf_fields ImgFields1, Intf_fields ImgFields2)
{
/*
  // Check that the originating systems are the same
  if (ImgFields1.originating_system != ImgFields2.originating_system)
  {
    Cerr("***** IntfCheckDimensionsConsistency() -> Originating systems are not the same !" << endl);
    return 1;
  }
*/
  // Check that the numbers of dimensions are the same
  if (ImgFields1.nb_dims != ImgFields2.nb_dims)
  {
    Cerr("***** IntfCheckDimensionsConsistency() -> Numbers of dimensions are not the same !" << endl);
    return 1;
  }
  // Loop over all dimensions
  for (int dim=0; dim<((int)ImgFields1.nb_dims); dim++)
  {
    // Check the size of this dimension
    if (ImgFields1.mtx_size[dim] != ImgFields2.mtx_size[dim])
    {
      Cerr("***** IntfCheckDimensionsConsistency() -> The sizes of the dimension " << dim+1 << " are not the same !" << endl);
      return 1;
    }
  }
  // For the first 3 dimensions, check the voxel size
  for (int dim=0; dim<std::min(3,((int)ImgFields1.nb_dims)); dim++)
  {
    // Check voxel sizes
    if (ImgFields1.vox_size[dim] != ImgFields2.vox_size[dim])
    {
      Cerr("***** IntfCheckDimensionsConsistency() -> Voxel sizes of dimension " << dim+1 << " are not the same !" << endl);
      return 1;
    }
  }
  // Check the slice thickness
  if (ImgFields1.slice_thickness_mm != ImgFields2.slice_thickness_mm)
  {
    Cerr("***** IntfCheckDimensionsConsistency() -> Slice thicknesses are not the same !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB* IntfLoadImageFromScratch( const string& a_pathToHeaderFile, 
                                 Intf_fields* ap_ImgFields,
                                 int vb )
{
  if (vb>=2) Cout("IntfLoadImageFromScratch() -> Read Interfile image '" << a_pathToHeaderFile << "'" << endl);

  // Init Interfile Key structure
  IntfKeyInitFields(ap_ImgFields);

  // Recover image infos from the header file
  if (IntfReadHeader(a_pathToHeaderFile, ap_ImgFields, vb))
  {
    Cerr("***** IntfLoadImageFromScratch() -> A problem occured while trying to read the interfile header '" << a_pathToHeaderFile << "' !" << endl);
    return NULL;
  }

  // Error if number of dimensions is more than 3
  if (ap_ImgFields->nb_dims>3)
  {
    Cerr("***** IntfLoadImageFromScratch() -> Cannot handle a number of dimensions higher than 3 !" << endl);
    return NULL;
  }

  // Check that the numbers of voxels have been read successfully (they are initialized to 0 by the IntfKeyInitFields functions
  for (int d=0; d<ap_ImgFields->nb_dims; d++) if (ap_ImgFields->mtx_size[d]==0)
  {
    Cerr("***** IntfLoadImageFromScratch() -> Number of voxels for dimension #" << d << " is 0, so has not been read correctly in header file '" << a_pathToHeaderFile << "' !" << endl);
    return NULL;
  }
  // Same for voxel sizes, initialized to -1 by default
  for (int d=0; d<ap_ImgFields->nb_dims; d++) if (ap_ImgFields->vox_size[d]<=0.)
  {
    Cerr("***** IntfLoadImageFromScratch() -> Voxel size for dimension #" << d << " is negative, so has not been read correctly in header file '" << a_pathToHeaderFile << "' !" << endl);
    return NULL;
  }

  // Compute total number of voxels
  INTNB dim_tot = ap_ImgFields->mtx_size[0] * ap_ImgFields->mtx_size[1] * ap_ImgFields->mtx_size[2];

  // Allocate image
  FLTNB* p_image = (FLTNB*)malloc(dim_tot*sizeof(FLTNB));

  // Open image data file
  ifstream img_file(ap_ImgFields->path_to_image.c_str(), ios::binary | ios::in);
  if (!img_file)
  {
    Cerr("***** IntfLoadImageFromScratch() -> Input image file '" << ap_ImgFields->path_to_image << "' is missing or corrupted !" << endl);
    return NULL;
  }

  // Get the position of the beginning of the image data
  uint32_t offset = ap_ImgFields->data_offset;

  // Call the right IntfReadData() function according to the pixel type, and read data
  IntfGetPixelTypeAndReadData(*ap_ImgFields, &img_file, p_image, NULL, &offset, dim_tot, vb);

  // Return the pointer to the image in memory
  return p_image;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int IntfWriteImageFromIntfFields(const string& a_pathToImg, FLTNB* ap_ImgMatrix, Intf_fields Img_fields, int vb)
{
  if (vb>=2) Cout("IntfWriteImageFromIntfFields() -> Write 3D image with output base name '" << a_pathToImg << "'" << endl);

  // Write Interfile header
  if (IntfWriteHeaderMainData(a_pathToImg, Img_fields, vb) )
  {
    Cerr("***** IntfWriteImageFromIntfFields() -> Error : while trying to write the interfile header '"<< a_pathToImg << "' !" << endl);
    return 1;
  }

  // Binary data file
  string path_to_image = a_pathToImg;
  path_to_image.append(".img");

  // Total number of voxels
  uint32_t dim_tot = Img_fields.mtx_size[0] * Img_fields.mtx_size[1] * Img_fields.mtx_size[2];

  // Read Interfile image
  if (IntfWriteImage(path_to_image, ap_ImgMatrix, dim_tot, vb) )
  {
    Cerr("***** IntfWriteImageFromIntfFields() -> Error : while trying to write the interfile image '"<< path_to_image << "' !" << endl);
    return 1;
  }
  
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int IntfReadImage(   const string& a_pathToHeaderFile, 
                         FLTNB* ap_ImgMatrix,
oImageDimensionsAndQuantification* ap_ID, 
                               int a_verbose, 
                              bool a_lerpFlag)
{
  // Verbose
  if (a_verbose>=VERBOSE_DETAIL) Cout("oInterfileIO::IntfReadImage() -> Read image from interfile header '" << a_pathToHeaderFile << "'" << endl);

  // Init Interfile Key structure
  Intf_fields Img_fields;
  IntfKeyInitFields(&Img_fields);

  // Recover image infos from the header file
  if (IntfReadHeader(a_pathToHeaderFile, &Img_fields, a_verbose) )
  {
    Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while trying to read the interfile header '" << a_pathToHeaderFile << "' !" << endl);
    return 1;
  }

  // Check consistency between image interfile data and the algorithm data
  if (IntfCheckConsistency(&Img_fields, ap_ID, a_verbose, a_lerpFlag) )
  {
    Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while checking consistencies between reconstruction parameters and interfile keys from header '"
         << a_pathToHeaderFile << "' !" << endl);
    return 1;
  }

  // If interpolation required, allocate matrix with original image dimensions.
  FLTNB* pimg_erp = NULL;
  IntfAllocInterpImg(&pimg_erp, Img_fields);
    
  // Read image data
  ifstream img_file(Img_fields.path_to_image.c_str(), ios::binary | ios::in);
  if (img_file)
  {
    // Get the position of the beginning of the image data
    uint32_t offset = Img_fields.data_offset;
    // Call the right IntfReadData() function according to the pixel type, and read data
    if (IntfGetPixelTypeAndReadData(Img_fields, &img_file, ap_ImgMatrix, pimg_erp, &offset, ap_ID->GetNbVoxXYZ(), a_verbose))
    {
      Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while reading data and eventually interpolating from the IntfGetPixelTypeAndReadData() function !" << endl);
      return 1;
    }
  }
  else
  {
    Cerr("***** oInterfileIO::IntfReadImage() -> An error occurred while trying to open the image file '" << Img_fields.path_to_image << "' !" << endl);
    // If interpolation was required, deallocate image matrix
    if (pimg_erp) delete[] pimg_erp; 
    return 1;
  }
  
  // If interpolation was required, deallocate image matrix
  if (pimg_erp) delete[] pimg_erp; 
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int IntfReadImage(   const string& a_pathToHeaderFile, 
                      FLTNB**** a4p_ImgMatrix, 
oImageDimensionsAndQuantification* ap_ID, 
                               int a_verbose, 
                              bool a_lerpFlag)
{
  if (a_verbose >= 5) Cout("oInterfileIO::IntfReadImage() -> Read image from interfile header '" << a_pathToHeaderFile << "'" << endl);
  
  // Init Interfile Key structure
  Intf_fields Img_fields;
  IntfKeyInitFields(&Img_fields);
  
  // List containing single/multiple interfile headers of the images to read
  vector<string> lpath_to_headers;

  // First check if the provided file is a metaheader with multiple files
  // and recover the list of header file
  if (IntfIsMHD(a_pathToHeaderFile, lpath_to_headers) < 0 ) // error
  {
    Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while trying to read the interfile header '" << a_pathToHeaderFile << "' in IntfIsMHD() function !" << endl);
    return 1;
  }

  // --- Case 1: We have one single data file containing all the data --- //
  if (lpath_to_headers.size() == 1)
  {
    // Recover image infos from either metaheader or single header file
    if (IntfReadHeader(a_pathToHeaderFile, &Img_fields, a_verbose) )
    {
      Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while trying to read the interfile header '" << a_pathToHeaderFile << "' !" << endl);
      return 1;
    }
    // Check consistency between image interfile data and the algorithm data
    if (IntfCheckConsistency(&Img_fields, ap_ID, a_verbose, a_lerpFlag) )
    {
      Cerr("***** oInterfileIO::IntfReadImage() -> A error occured while checking consistencies between reconstruction parameters and interfile keys in the header '" << a_pathToHeaderFile << "' !" << endl);
      return 1;
    }
    // If interpolation required, instanciate matrix with original image dimensions.
    FLTNB* pimg_erp = NULL;
    IntfAllocInterpImg(&pimg_erp, Img_fields);
    // Open file
    ifstream img_file(Img_fields.path_to_image.c_str(), ios::binary | ios::in);
    if (img_file)
    {
      // Get the position of the beginning of the image data
      uint32_t offset = Img_fields.data_offset;
      // Call the right IntfReadData() function according to the pixel type, and read data.
      // offset is updated with the called functions within the loops
      // TODO : RAJOUTER DES CHECKS SUR LES DIMENSIONS DYNAMIQUES DANS CHECKCONSISTENCY (gérer différences sensisitivityGenerator & reconstruction) !!!
      for (int d1=0 ; d1<Img_fields.nb_time_frames ; d1++)
        for (int d2=0 ; d2<Img_fields.nb_resp_gates ; d2++)
          for (int d3=0 ; d3<Img_fields.nb_card_gates ; d3++)
          {
            if (IntfGetPixelTypeAndReadData(Img_fields, &img_file, a4p_ImgMatrix[d1][d2][d3], pimg_erp, &offset, ap_ID->GetNbVoxXYZ(), a_verbose))
            {
              Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while reading data and eventually interpolating from the IntfGetPixelTypeAndReadData() function !" << endl);
              // If interpolation was required, deallocate image matrix with original image dimensions
              if (pimg_erp) delete[] pimg_erp;
              return 1;
            }
          }
    }
    else
    {
      Cerr("***** oInterfileIO::IntfReadImage() -> An error occurred while trying to read the image file '"<< Img_fields.path_to_image << "' !" << endl);
      // If interpolation was required, deallocate image matrix with original image dimensions
      if (pimg_erp) delete[] pimg_erp; 
      return 1;
    }
    // If interpolation was required, deallocate image matrix with original image dimensions
    if (pimg_erp) delete[] pimg_erp; 
  }

  // --- Case 2 : We have a number of data file equal to the number of image to load --- //
  else 
  {
    // Recover image infos from the metaheader
    if (IntfReadHeader(a_pathToHeaderFile, &Img_fields, a_verbose) )
    {
      Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while trying to read the interfile metaheader '" << a_pathToHeaderFile << "' !" << endl);
      return 1;
    }
    // Get dimensions from ImageDimensions object
    int dims[3];
    dims[0] = ap_ID->GetNbTimeFrames();
    dims[1] = ap_ID->GetNbRespGates();
    dims[2] = ap_ID->GetNbCardGates();
    // Check we have a number of file corresponding to the number of images to load
    if (lpath_to_headers.size() != (uint32_t)(dims[0]*dims[1]*dims[2]))
    {
      Cerr("***** oInterfileIO::IntfReadImage() -> Number of interfile images (" << lpath_to_headers.size() << 
           ") not consistent with the number of images to load (" << dims[0]*dims[1]*dims[2]<< ") !" << endl);
      return 1;
    }
    // If interpolation required, instanciate matrix with original image dimensions.
    FLTNB* pimg_erp = NULL;
    IntfAllocInterpImg(&pimg_erp, Img_fields);
    // Loop on dynamic image files
    for(int d1=0 ; d1<dims[0] ; d1++)
      for(int d2=0 ; d2<dims[1] ; d2++)
        for(int d3=0 ; d3<dims[2] ; d3++)
        {
          int idx_img = d1*dims[1]*dims[2] + d2*dims[2] + d3;
          // Recover image infos from the specific interfile header
          // Todo : Check if image 3D dimensions are different from an image to another ?
          //        (very unlikely, but it would cause segfault if interpolation is enabled)
          if (IntfReadHeader(lpath_to_headers[idx_img], &Img_fields, a_verbose) )
          {
            Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while trying to read the interfile header '" << lpath_to_headers[idx_img] << "' !" << endl);
            return 1;
          }
          // Check consistency between image interfile data and the algorithm data
          if (IntfCheckConsistency(&Img_fields, ap_ID, a_verbose, a_lerpFlag) )
          {
            Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while checking consistencies between reconstruction parameters and interfile keys in the header '"
                                               << lpath_to_headers[idx_img] << "' !" << endl);
            return 1;
          }
          // Open file
          ifstream img_file(Img_fields.path_to_image.c_str(), ios::binary | ios::in);
          if (img_file)
          {
            // Get the position of the beginning of the image data (assuming that offset is the same for each file (if any) )
            uint32_t offset = Img_fields.data_offset; 
            // Call the right IntfReadData() function according to the pixel type, and read data (offset not updated)
            if (IntfGetPixelTypeAndReadData(Img_fields, &img_file, a4p_ImgMatrix[d1][d2][d3], pimg_erp, &offset, ap_ID->GetNbVoxXYZ(), a_verbose))
            {
              Cerr("***** oInterfileIO::IntfReadImage() -> An error occured while reading data and eventually interpolating from the IntfGetPixelTypeAndReadData() function !" << endl);
              // If interpolation was required, deallocate image matrix with original image dimensions
              if (pimg_erp) delete[] pimg_erp;
              return 1;
            }
          }
          else
          {
            Cerr("***** oInterfileIO::IntfReadImage() -> An error occurred while trying to read the image file '"<< Img_fields.path_to_image << "' !" << endl);
            // If interpolation was required, deallocate image matrix with original image dimensions
            if (pimg_erp) delete[] pimg_erp;
            return 1;
          }
        }
    // If interpolation was required, deallocate image matrix with original image dimensions
    if (pimg_erp) delete[] pimg_erp;
  }

  /* If gating is enabled with separate 3D gated images, the image_duration key may be in each header
     The following check will be broken in this case.
  // Check some dynamic recovered infos
  if(Img_fields.nb_time_frames > 1 &&
     Img_fields.nb_time_frames != Img_fields.image_duration.size())
  {
    Cerr("***** IntfReadImage()-> Error : the number of provided image duration:  '"<< Img_fields.image_duration.size() 
      << "' does not match the number of frames '"<< Img_fields.nb_time_frames <<"' !" << endl);
    return 1;
  }
  */
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int IntfReadImgDynCoeffFile(const string& a_pathToHeaderFile, 
                               FLTNB** a2p_ImgMatrix, 
       oImageDimensionsAndQuantification* ap_ID,
                                      int a_nbFbasis, 
                                      int a_verbose, 
                                     bool a_lerpFlag)
{
  if (a_verbose >= 5) Cout("IntfReadImgDynCoeffFile" << endl);

  // Init Interfile Key structure
  Intf_fields Img_fields;
  IntfKeyInitFields(&Img_fields);
  
  // List containing single/multiple interfile headers of the images to read
  vector<string> lpath_to_headers;

  // First check if the provided file is a metaheader with multiple files
  // and recover the list of image header file
  if(IntfIsMHD(a_pathToHeaderFile, lpath_to_headers) <0) // error
  {
    Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> an error occured while trying to read the interfile header '" << a_pathToHeaderFile << "' !" << endl);
    return 1;
  }

  // --- Case 1: We have one single data file containing all the data --- //
  if(lpath_to_headers.size() == 1)
  {
    // Recover image infos from either metaheader or single header file
    if (IntfReadHeader(a_pathToHeaderFile, &Img_fields, a_verbose) )
    {
      Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> An error occured while trying to read the interfile header '" << a_pathToHeaderFile << "' !" << endl);
      return 1;
    }
    // Check consistency between image interfile data and the algorithm data
    if (IntfCheckConsistency(&Img_fields, ap_ID, a_verbose, a_lerpFlag) )
    {
      Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> An error occured while checking consistencies between reconstruction parameters and interfile keys in the header '" << a_pathToHeaderFile << "' !" << endl);
      return 1;
    }
    // If interpolation required, instanciate matrix with original image dimensions.
    FLTNB* pimg_erp = NULL;
    IntfAllocInterpImg(&pimg_erp, Img_fields);  
    // Open file
    ifstream img_file(Img_fields.path_to_image.c_str(), ios::binary | ios::in);
    if (img_file)
    {
      // Get the position of the beginning of the image data
      uint32_t offset = Img_fields.data_offset;
      // Call the right IntfReadData() function according to the pixel type, and read data
      // offset is updated with the called functions within the loops
      for (int bf=0 ; bf<a_nbFbasis ; bf++)
      {
        if (IntfGetPixelTypeAndReadData(Img_fields, &img_file, a2p_ImgMatrix[bf], pimg_erp, &offset, ap_ID->GetNbVoxXYZ(), a_verbose))
        {
          Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> An error occured while reading from file stream and eventually interpolating !" << endl);
          // If interpolation was required, deallocate image matrix with original image dimensions
          if (pimg_erp) delete[] pimg_erp;
          return 1;
        }
      }
    }
    else
    {
      Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> An error occurred while trying to read the image file '" << Img_fields.path_to_image << "' !" << endl);
      // If interpolation was required, deallocate image matrix with original image dimensions
      if (pimg_erp) delete[] pimg_erp; 
      return 1;
    }
    // If interpolation was required, deallocate image matrix with original image dimensions
    if (pimg_erp) delete[] pimg_erp; 
  }

  // --- Case 2: We have a number of data file equal to the number of image to load --- //
  else 
  {
    // Recover image infos from the metaheader
    if (IntfReadHeader(a_pathToHeaderFile, &Img_fields, a_verbose) )
    {
      Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> An error occured while trying to read the interfile metaheader '" << a_pathToHeaderFile << "' !" << endl);
      return 1;
    }
    // Check we have a number of file corresponding to the number of images to load
    if (lpath_to_headers.size() != (uint32_t)a_nbFbasis)
    {
      Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> Number of interfile images (" << lpath_to_headers.size() << 
           ") not consistent with the number of parametric images (" << a_nbFbasis<< ") !" << endl);
      return 1;
    }
    // If interpolation required, instanciate matrix with original image dimensions.
    FLTNB* pimg_erp = NULL;
    IntfAllocInterpImg(&pimg_erp, Img_fields); 
    // Loop
    for (int bf=0 ; bf<a_nbFbasis ; bf++)
    {
      // Recover image infos from the specific interfile header
      if (IntfReadHeader(lpath_to_headers[bf], &Img_fields, a_verbose) )
      {
        Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> an error occured while trying to read interfile header '" << a_pathToHeaderFile << "' !" << endl);
        // If interpolation was required, deallocate image matrix with original image dimensions
        if (pimg_erp) delete[] pimg_erp; 
        return 1;
      }
      // Check consistency between image interfile data and the algorithm data
      if (IntfCheckConsistency(&Img_fields, ap_ID, a_verbose, a_lerpFlag) )
      {
        Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> An error occured while checking consistencies between reconstruction parameters and interfile keys in the header '"<< a_pathToHeaderFile << "' !" << endl);
        // If interpolation was required, deallocate image matrix with original image dimensions
        if (pimg_erp) delete[] pimg_erp; 
        return 1;
      }
      // Open file stream
      ifstream img_file(Img_fields.path_to_image.c_str(), ios::binary | ios::in);
      if (img_file)
      {
        // Get the position of the beginning of the image data (assuming that offset is the same for each file (if any) )
        uint32_t offset = Img_fields.data_offset; 
        // Call the right IntfReadData() function according to the pixel type, and read data (offset not updated)
        if (IntfGetPixelTypeAndReadData(Img_fields, &img_file, a2p_ImgMatrix[bf], pimg_erp, &offset, ap_ID->GetNbVoxXYZ(), a_verbose))
        {
          Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> An error occured while reading image stream and eventually interpolating !" << endl);
          // If interpolation was required, deallocate image matrix with original image dimensions
          if (pimg_erp) delete[] pimg_erp;
          return 1;
        }
      }
      else
      {
        Cerr("***** oInterfileIO::IntfReadImgDynCoeffFile() -> An error occurred while trying to read the image file '" << Img_fields.path_to_image << "' !" << endl);
        // If interpolation was required, deallocate image matrix with original image dimensions
        if (pimg_erp) delete[] pimg_erp; 
        return 1;
      }
    }
    // If interpolation was required, deallocate image matrix with original image dimensions
    if (pimg_erp) delete[] pimg_erp; 
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfAllocInterpImg
  \param a2p_img : Pointer to 1 dimensional image matrix to recover the image to interpolate
  \param ap_IF : Structure containing the interfile fields read in a interfile header
  \brief Allocate memory for an image matrix to recover an image to interpolate
*/
void IntfAllocInterpImg(FLTNB **a2p_img, Intf_fields a_IF)
{
  if(a_IF.is_mtx_size_different == true)
  {
    uint32_t nb_vox = a_IF.mtx_size[0] *
                      a_IF.mtx_size[1] *
                      a_IF.mtx_size[2];

    // Allocate image and
    *a2p_img = new FLTNB[ nb_vox ];
    ::memset(*a2p_img, 0, sizeof(FLTNB) * nb_vox);
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfCheckConsistency
  \param ap_IF : Structure containing the interfile fields read in a interfile header
  \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
  \param vb : Verbosity level
  \param a_lerpFlag : if true, enable linear interpolation of the image if img dimensions differ from the reconstruction dimensions
  \brief Check if the mandatory fields have been initialize in the ap_IF structure, and check consistencies with the reconstruction parameters
  \details This function also checks if the matrix size of the original image is different to the reconstruction matrix size.
           In this case a boolean is set up to enable interpolation during image reading
  \todo Add check for all mandatory parameters, and temporal image dimensions
  \todo Float comparison ?
  \return 0 if success, positive value otherwise.
*/
int IntfCheckConsistency(Intf_fields* ap_IF, oImageDimensionsAndQuantification* ap_ID, int vb, int a_lerpFlag)
{
  if(vb >= 5) Cout("IntfCheckConsistency()" << endl);
  
  // Check if main keys have been initialized
  if(ap_IF->path_to_image.size()<1 ||
     ap_IF->nb_format == "" ||
     ap_IF->mtx_size[0]==0 || ap_IF->mtx_size[1]==0 || ap_IF->mtx_size[2]==0 )
    {
      Cerr("***** IntfCheckConsistency()-> Error : some mandatory keys not initialized. Cannot read the interfile image !" << endl);
      // Print the related key 
      if(ap_IF->path_to_image.size()<1)
        Cerr("                               Error when trying to read path to image data" << endl);
      if(ap_IF->nb_format == "")
        Cerr("                               Error when trying to read data voxel type " << endl);
      if(ap_IF->mtx_size[0]==0 || ap_IF->mtx_size[1]==0 || ap_IF->mtx_size[2]==0 )
        Cerr("                               Error when trying to read matrix size (image dimensions) : x= " << ap_IF->mtx_size[0] << ", y= " << ap_IF->mtx_size[1] << ", z= " << ap_IF->mtx_size[2]<< endl);
      return 1;
    }
  
  // If voxel size not found, throw warning
  if( ap_IF->vox_size[0]<0 || ap_IF->vox_size[1]<0 || ap_IF->vox_size[2]<0)
    {
      if(vb == 5)
        Cerr("***** IntfCheckConsistency()-> WARNING : No information found about voxel size ('scaling factor (mm/pixel)' tags.). Missing voxel sizes will be set to 1mm !" << endl);

      if(ap_IF->vox_size[0]<0) ap_IF->vox_size[0] = 1.;
      if(ap_IF->vox_size[1]<0) ap_IF->vox_size[1] = 1.;
      if(ap_IF->vox_size[2]<0) ap_IF->vox_size[2] = 1.;
      if(vb == 5)
        Cerr("                                         Voxel sizes : x= " << ap_IF->vox_size[0] << ", y= " << ap_IF->vox_size[1] << ", z= " << ap_IF->vox_size[2]<< endl);
    }
    
    
  // If original dimensions don't match reconstructions dimensions/voxel sizes,
  // recover this data in Intf_fields object (if a_lerpFlag==true) or throw error (if a_lerpFlag==false)
  if( ((INTNB)(ap_IF->mtx_size[0])) != ap_ID->GetNbVoxX()  ||
      ((INTNB)(ap_IF->mtx_size[1])) != ap_ID->GetNbVoxY()  ||
      ((INTNB)(ap_IF->mtx_size[2])) != ap_ID->GetNbVoxZ()  ||
      !FLTNBIsEqual(ap_IF->vox_size[0], ap_ID->GetVoxSizeX(), pow(0.1, std::numeric_limits<FLTNB>::digits10+1) ) ||
      !FLTNBIsEqual(ap_IF->vox_size[1], ap_ID->GetVoxSizeY(), pow(0.1, std::numeric_limits<FLTNB>::digits10+1) ) ||
      !FLTNBIsEqual(ap_IF->vox_size[2], ap_ID->GetVoxSizeZ(), pow(0.1, std::numeric_limits<FLTNB>::digits10+1) ) )
  {
    if(a_lerpFlag)
    {
      ap_IF->cmtx_size[0] = (uint32_t)(ap_ID->GetNbVoxX());
      ap_IF->cmtx_size[1] = (uint32_t)(ap_ID->GetNbVoxY());
      ap_IF->cmtx_size[2] = (uint32_t)(ap_ID->GetNbVoxZ());
      ap_IF->cvox_size[0] = ap_ID->GetVoxSizeX();
      ap_IF->cvox_size[1] = ap_ID->GetVoxSizeY();
      ap_IF->cvox_size[2] = ap_ID->GetVoxSizeZ();
      ap_IF->cvox_offset[0] = ap_ID->GetOffsetX();
      ap_IF->cvox_offset[1] = ap_ID->GetOffsetY();
      ap_IF->cvox_offset[2] = ap_ID->GetOffsetZ();
      ap_IF->is_mtx_size_different = true; // Set this boolean to true to indicate that an interpolation step will be required
    }
    else
    {
      Cerr("***** IntfCheckConsistency()-> Error : Image dimensions don't match reconstructions dimensions/voxel sizes" << endl);
      Cerr("                                       and linear interpolation is disabled (a_lerpFlag is false) !" << endl);
      Cerr("*****                                  Recovered image dimensions (x;y;z): "<< ap_IF->mtx_size[0] <<" ; "<< ap_IF->mtx_size[1] << " ; " << ap_IF->mtx_size[2] << endl);
      Cerr("*****                                  Reconstruction dimensions (x;y;z) : "<< ap_ID->GetNbVoxX() <<" ; "<< ap_ID->GetNbVoxY() << " ; " << ap_ID->GetNbVoxZ() << endl);
      Cerr("*****                                  Image voxel sizes (x;y;z)         : "<< ap_IF->vox_size[0] <<" ; "<< ap_IF->vox_size[1] << " ; " << ap_IF->vox_size[2] << endl);
      Cerr("*****                                  Reconstruction voxel sizes (x;y;z): "<< ap_ID->GetVoxSizeX() <<" ; "<< ap_ID->GetVoxSizeY() << " ; " << ap_ID->GetVoxSizeZ() << endl);
      return 1;  
    }
  }
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int IntfGetPixelTypeAndReadData(Intf_fields a_IF, 
                                  ifstream* ap_iFile, 
                                  FLTNB* ap_outImgMatrix, 
                                  FLTNB* ap_inImgMatrix, 
                                  uint32_t* ap_offset, 
                                        int a_nbVox, 
                                        int a_verbose)
{
  if (a_verbose >= 5) Cout("oInterfileIO::IntfGetPixelTypeAndReadData() " << endl);
  
  // To check if an error occurred in one of the called functions
  int error = 0;
  
  // Check the image data voxel size according to the Interfile fields 'number_format' and 'number of bytes per pixel'
  if (a_IF.nb_format == FLT32_str || a_IF.nb_format == FLT32_str2) // Float data
  {
    float flt;
    error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &flt);
  }
  else if (a_IF.nb_format == FLT64_str) // Double data
  {
    double db;
    error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &db);
  }
  else if (a_IF.nb_format == LONGDOUBLE_str) //long double data
  {
    if (a_IF.nb_bytes_pixel != sizeof(long double))
    {
      Cerr("***** oInterfileIO::IntfGetPixelTypeAndReadData() -> The long double format is not the same for this image file and for this platform !" << endl);
      return 1;
    }
    long double db;
    error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &db);
  }
  else if (a_IF.nb_format == INT32_str) // Signed integer data
  {
    if(a_IF.nb_bytes_pixel == 1)
    {
      int8_t s_int;
      error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &s_int);
    }
    else if(a_IF.nb_bytes_pixel == 2)
    {
      int16_t s_int;
      error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &s_int);
    }
    else if(a_IF.nb_bytes_pixel == 8)
    {
      int64_t s_int;
      error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &s_int);
    }
    else // default : 4 bytes
    {
      int32_t s_int;
      error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &s_int);
    }
  }
  else // Unsigned integer data (default)
  {
    if(a_IF.nb_bytes_pixel == 1)
    {
      uint8_t u_int;
      error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &u_int);
    }
    else if(a_IF.nb_bytes_pixel == 2)
    {
      uint16_t u_int;
      error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &u_int);
    }
    else if(a_IF.nb_bytes_pixel == 8)
    {
      uint64_t u_int;
      error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &u_int);
    }
    else // default : 4 bytes
    {
      uint32_t u_int;
      error = IntfReadData(a_IF, ap_iFile, ap_outImgMatrix, ap_inImgMatrix, ap_offset, a_nbVox, a_verbose, &u_int);
    }
  }
  
  if (error!=0)
  {
    Cerr("***** oInterfileIO::IntfGetPixelTypeAndReadData() -> An error occurred when trying to read data through the IntfReadData() function !" << endl);
    return 1;
  }
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

template <class T>
int IntfReadData(Intf_fields a_IF, 
                   ifstream* ap_iFile, 
                   FLTNB* ap_outImgMatrix, 
                   FLTNB* ap_inImgMatrix, 
                   uint32_t* a_offset, 
                         int a_nbVox, 
                         int a_verbose, 
                          T* bytes)
{
  if (a_verbose >= 5) Cout("oInterfileIO::IntfReadData() -> Read data from file stream and eventually interpolate" << endl);
  
  int nb_vox = 0; // number of voxels in the matrix which will recover the image data
  FLTNB* pimg = NULL; // pointer to the matrix which will recover the image data

  // Check if interpolation is required
  if(a_IF.is_mtx_size_different)
  {
    // Set the number of voxels in the image to the original image dimensions as they differ from reconstruction dimensions
    nb_vox = a_IF.mtx_size[0]*a_IF.mtx_size[1]*a_IF.mtx_size[2];
    pimg = ap_inImgMatrix; // recover the image into the image matrix initialized with interfile dimensions
  }
  else
  {
    nb_vox = a_nbVox;
    pimg = ap_outImgMatrix; // directly recover the image into the CASToR image matrix 
  }

  // Allocate temporary reading buffer
  bytes = (T*)malloc(nb_vox*sizeof(T));
  // Seek to the provided position
  ap_iFile->seekg(*a_offset);
  // Read the data
  ap_iFile->read((char*)bytes, nb_vox*sizeof(T));
  // Check reading status
  if (!ap_iFile->good())
  {
    if (ap_iFile->eof()) Cerr("***** oInterfileIO::IntfReadData() -> Not enough data in the file stream !" << endl);
    else Cerr("***** oInterfileIO::IntfReadData() -> An error occured while reading from file stream !" << endl);
    free(bytes);
    return 1;
  }

  // Loop over all voxels
  for (int v=0; v<nb_vox; v++)
  {
    // Todo Re-orient image data from the original orientation to the default orientation (transaxial/supine/headin)
    //int voxel_idx = IntfGetVoxIdxSHTOrientation(a_IF, v);
    int voxel_idx = v; // For now no orientations
    
    // Swap data if endianness is not the same as user processor
    if (a_IF.endianness != User_Endianness) SwapBytes(&(bytes[v]));
    // Cast to CASToR image format
    FLTNB buffer = (FLTNB)(bytes[v]);
    // Calibrate data using rescale slope and intercept if needed.
    // Then write in the image matrix
    pimg[voxel_idx] = buffer * a_IF.rescale_slope * a_IF.quant_units
                    + a_IF.rescale_intercept; 

    // set the offset to the next position
    *a_offset += sizeof(T);
  }

  free(bytes);

  // END OF NEW CODE

  // Call interpolation function if required
  if (a_IF.is_mtx_size_different 
   || a_IF.vox_offset[0] != a_IF.vox_offset[0] 
   || a_IF.vox_offset[1] != a_IF.vox_offset[1] 
   || a_IF.vox_offset[2] != a_IF.vox_offset[2] )
  {
    if (ImageInterpolation(ap_inImgMatrix, ap_outImgMatrix, 
                           a_IF.mtx_size,   a_IF.cmtx_size,
                           a_IF.vox_size,   a_IF.cvox_size,
                           a_IF.vox_offset, a_IF.cvox_offset) )
    {
      Cerr("***** oInterfileIO::IntfReadData() -> An error occurred while interpolating the input image to the reconstruction dimensions !" << endl);
      return 1;
    }
  }

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
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
                 const FLTNB ap_iOffVox[3],    const FLTNB ap_oOffVox[3] )
{
  // Assuming the two images are registered
  // todo : deal with any potential offset if the input data contains following keys
  // "origin (mm) [x], offset [x], first pixel offset (mm) [x] "

/* ORIGINAL FUNCTIONS IMPLEMENTED ONLY WITH FLOATS
  float const posOldImage[] = {-(float)ap_iDimVox[0]*ap_iSizeVox[0]/2 ,
                               -(float)ap_iDimVox[1]*ap_iSizeVox[1]/2 ,
                               -(float)ap_iDimVox[2]*ap_iSizeVox[2]/2 };
                               
  float const posNewImage[] = {-(float)ap_oDimVox[0]*ap_oSizeVox[0]/2 ,
                               -(float)ap_oDimVox[1]*ap_oSizeVox[1]/2 ,
                               -(float)ap_oDimVox[2]*ap_oSizeVox[2]/2 };
  
  // Padding the old image in each direction
  uint32_t const iDimVoxPad[] 
  {
    ap_iDimVox[ 0 ] + 2,
    ap_iDimVox[ 1 ] + 2,
    ap_iDimVox[ 2 ] + 2
  };

  uint32_t const nElts = iDimVoxPad[ 0 ] *
                         iDimVoxPad[ 1 ] * 
                         iDimVoxPad[ 2 ];

  // Allocating a new buffer storing the padded image
  float *pPadIData = new float[ nElts ];
  ::memset( pPadIData, 0, sizeof( float ) * nElts );
  
  for( uint32_t k = 0; k < ap_iDimVox[ 2 ]; ++k )
    for( uint32_t j = 0; j < ap_iDimVox[ 1 ]; ++j )
      for( uint32_t i = 0; i < ap_iDimVox[ 0 ]; ++i )
      {
        pPadIData[ ( i + 1 ) + ( j + 1 ) * iDimVoxPad[ 0 ]
          + ( k + 1 ) * iDimVoxPad[ 0 ] * iDimVoxPad[ 1 ] ] =
          ap_iImg[ i + j * ap_iDimVox[ 0 ]
            + k * ap_iDimVox[ 0 ] * ap_iDimVox[ 1 ] ];
      }
    

  // Computing the bounds in each direction depending on the pad
  float const boundMin[] 
  {
    posOldImage[ 0 ] - ap_iSizeVox[ 0 ] / 2.0f,
    posOldImage[ 1 ] - ap_iSizeVox[ 1 ] / 2.0f,
    posOldImage[ 2 ] - ap_iSizeVox[ 2 ] / 2.0f
  };
  
  float const boundMax[] 
  {
    posOldImage[ 0 ] + (float)ap_iDimVox[ 0 ] * ap_iSizeVox[ 0 ]
      + ap_iSizeVox[ 0 ] / 2.0f,
    posOldImage[ 1 ] + (float)ap_iDimVox[ 1 ] * ap_iSizeVox[ 1 ]
      + ap_iSizeVox[ 1 ] / 2.0f,
    posOldImage[ 2 ] + (float)ap_iDimVox[ 2 ] * ap_iSizeVox[ 2 ]
      + ap_iSizeVox[ 2 ] / 2.0f
  };
  
  // Computing and storing the position of the center of the pixel in each
  // direction for the initial image
  float **pOldCoordCenter = new float*[ 3 ];
  
  for( uint32_t i = 0; i < 3; ++i )
  {
    pOldCoordCenter[ i ] = new float[ iDimVoxPad[ i ] ];
    // Set the values
    for( uint32_t j = 0; j < iDimVoxPad[ i ]; ++j )
    {
      pOldCoordCenter[ i ][ j ] = posOldImage[ i ] - ap_iSizeVox[ i ] / 2.0f
        + j * ap_iSizeVox[ i ];
    }
  }

  // Computing and storing the position of the center of the pixel in each
  // direction of the resampled image
  float **pNewCoordCenter = new float*[ 3 ];
  for( uint32_t i = 0; i < 3; ++i )
  {
    pNewCoordCenter[ i ] = new float[ ap_oDimVox[ i ] ];
    // Set the values
    for( uint32_t j = 0; j < ap_oDimVox[ i ]; ++j )
    {
      pNewCoordCenter[ i ][ j ] = posNewImage[ i ] + ap_oSizeVox[ i ] / 2.0f
        + j * ap_oSizeVox[ i ];
    }
  }

  // Defining parameters
  float const invSizeX = 1.0f / ap_iSizeVox[ 0 ];
  float const invSizeY = 1.0f / ap_iSizeVox[ 1 ];
  float const invSizeZ = 1.0f / ap_iSizeVox[ 2 ];

  // Allocating memory for the 8 pixels kernel
  float *pKernelData = new float[ 8 ];

  // Loop over the elements of the new images
  for( uint32_t k = 0; k < ap_oDimVox[ 2 ]; ++k )
  {
    // Get the coordinate in Z
    float const z = pNewCoordCenter[ 2 ][ k ];
    if( z < boundMin[ 2 ] || z > boundMax[ 2 ] ) continue;

    // Find the bin in the old image
    int32_t const zBin = ( z - boundMin[ 2 ] ) * invSizeZ;

    // Computing the 2 z composants
    float const zComposantI0 = invSizeZ * ( pOldCoordCenter[ 2 ][ zBin + 1 ] - z );
    float const zComposantI1 = invSizeZ * ( z - pOldCoordCenter[ 2 ][ zBin ] );

    for( uint32_t j = 0; j < ap_oDimVox[ 1 ]; ++j )
    {
      // Get the coordinate in Y
      float const y = pNewCoordCenter[ 1 ][ j ];
      if( y < boundMin[ 1 ] || y > boundMax[ 1 ] ) continue;

      // Find the bin in the old image
      int32_t const yBin = ( y - boundMin[ 1 ] ) * invSizeY;

      // Computing the 2 y composants
      float const yComposantI0 = invSizeY * ( pOldCoordCenter[ 1 ][ yBin + 1 ]
        - y );
      float const yComposantI1 = invSizeY * ( y
        - pOldCoordCenter[ 1 ][ yBin ] );

      for( uint32_t i = 0; i < ap_oDimVox[ 0 ]; ++i )
      {
        // Get the coordinate in X
        float const x = pNewCoordCenter[ 0 ][ i ];
        if( x < boundMin[ 0 ] || x > boundMax[ 0 ] ) continue;

        // Find the bin in the old image
        int32_t const xBin = ( x - boundMin[ 0 ] ) * invSizeX;

        // Computing the 2 x composants
        float const xComposantI0 = invSizeX * (
          pOldCoordCenter[ 0 ][ xBin + 1 ] - x );
        float const xComposantI1 = invSizeX * ( x
          - pOldCoordCenter[ 0 ][ xBin ] );

        // Fill the buffer storing the data
        for( uint32_t kk = 0; kk < 2; ++kk )
        {
          for( uint32_t jj = 0; jj < 2; ++jj )
          {
            for( uint32_t ii = 0; ii < 2; ++ii )
            {
              pKernelData[ ii + jj * 2 + kk * 2 * 2 ] =
                pPadIData[
                  ( xBin + ii ) +
                  ( yBin + jj ) * iDimVoxPad[ 0 ] +
                  ( zBin + kk ) * iDimVoxPad[ 0 ] * iDimVoxPad[ 1 ]
                ];
            }
          }
        }

        // Computing the interpolation
        // In X
        float const xInterpVal0 = pKernelData[ 0 ] * xComposantI0 +
          pKernelData[ 1 ] * xComposantI1;

        float const xInterpVal1 = pKernelData[ 2 ] * xComposantI0 +
          pKernelData[ 3 ] * xComposantI1;

        float const xInterpVal2 = pKernelData[ 4 ] * xComposantI0 +
          pKernelData[ 5 ] * xComposantI1;

        float const xInterpVal3 = pKernelData[ 6 ] * xComposantI0 +
          pKernelData[ 7 ] * xComposantI1;

        // In Y
        float const yInterpVal0 = xInterpVal0 * yComposantI0 +
          xInterpVal1 * yComposantI1;

        float const yInterpVal1 = xInterpVal2 * yComposantI0 +
          xInterpVal3 * yComposantI1;

        // Final interpolation
        float const interpValTot = yInterpVal0 * zComposantI0 +
          yInterpVal1 * zComposantI1;

        ap_oImg[ i + j * ap_oDimVox[ 0 ]
          + k * ap_oDimVox[ 0 ] * ap_oDimVox[ 1 ] ] = interpValTot;
      }
    }
  }

  // Freeing the memory
  for( uint32_t i = 0; i < 3; ++i )
  {
    delete[] pOldCoordCenter[ i ];
    delete[] pNewCoordCenter[ i ];
  }
  delete[] pOldCoordCenter;
  delete[] pNewCoordCenter;
  delete[] pPadIData;
  delete[] pKernelData;
  
  return 0;
*/



  FLTNB const posOldImage[] = {-((FLTNB)(ap_iDimVox[0]))*ap_iSizeVox[0]*((FLTNB)0.5) ,
                               -((FLTNB)(ap_iDimVox[1]))*ap_iSizeVox[1]*((FLTNB)0.5) ,
                               -((FLTNB)(ap_iDimVox[2]))*ap_iSizeVox[2]*((FLTNB)0.5) };
                               
  FLTNB const posNewImage[] = {-((FLTNB)(ap_oDimVox[0]))*ap_oSizeVox[0]*((FLTNB)0.5) ,
                               -((FLTNB)(ap_oDimVox[1]))*ap_oSizeVox[1]*((FLTNB)0.5) ,
                               -((FLTNB)(ap_oDimVox[2]))*ap_oSizeVox[2]*((FLTNB)0.5) };
                               
                               /*
  FLTNB const posOldImage[] = {-((FLTNB)(ap_iDimVox[0]))*ap_iSizeVox[0]*((FLTNB)0.5) + ap_iOffVox[0] ,
                               -((FLTNB)(ap_iDimVox[1]))*ap_iSizeVox[1]*((FLTNB)0.5) + ap_iOffVox[1] ,
                               -((FLTNB)(ap_iDimVox[2]))*ap_iSizeVox[2]*((FLTNB)0.5) + ap_iOffVox[2] };
                               
  FLTNB const posNewImage[] = {-((FLTNB)(ap_oDimVox[0]))*ap_oSizeVox[0]*((FLTNB)0.5) + ap_oOffVox[0],
                               -((FLTNB)(ap_oDimVox[1]))*ap_oSizeVox[1]*((FLTNB)0.5) + ap_oOffVox[1],
                               -((FLTNB)(ap_oDimVox[2]))*ap_oSizeVox[2]*((FLTNB)0.5) + ap_oOffVox[2]};*/
  
  // Padding the old image in each direction
  uint32_t const iDimVoxPad[] 
  {
    ap_iDimVox[ 0 ] + 2,
    ap_iDimVox[ 1 ] + 2,
    ap_iDimVox[ 2 ] + 2
  };

  uint32_t const nElts = iDimVoxPad[ 0 ] *
                         iDimVoxPad[ 1 ] * 
                         iDimVoxPad[ 2 ];

  // Allocating a new buffer storing the padded image
  FLTNB *pPadIData = new FLTNB[ nElts ];
  ::memset( pPadIData, 0, sizeof( FLTNB ) * nElts );
  
  for( uint32_t k = 0; k < ap_iDimVox[ 2 ]; ++k )
    for( uint32_t j = 0; j < ap_iDimVox[ 1 ]; ++j )
      for( uint32_t i = 0; i < ap_iDimVox[ 0 ]; ++i )
      {
        pPadIData[ ( i + 1 ) + ( j + 1 ) * iDimVoxPad[ 0 ]
          + ( k + 1 ) * iDimVoxPad[ 0 ] * iDimVoxPad[ 1 ] ] =
          ap_iImg[ i + j * ap_iDimVox[ 0 ]
            + k * ap_iDimVox[ 0 ] * ap_iDimVox[ 1 ] ];
      }
    

  // Computing the bounds in each direction depending on the pad
  FLTNB const boundMin[] 
  {
    posOldImage[ 0 ] - ap_iSizeVox[ 0 ] * ((FLTNB)0.5),
    posOldImage[ 1 ] - ap_iSizeVox[ 1 ] * ((FLTNB)0.5),
    posOldImage[ 2 ] - ap_iSizeVox[ 2 ] * ((FLTNB)0.5)
  };
  
  FLTNB const boundMax[] 
  {
    posOldImage[ 0 ] + ((FLTNB)ap_iDimVox[ 0 ]) * ap_iSizeVox[ 0 ]
      + ap_iSizeVox[ 0 ] * ((FLTNB)0.5),
    posOldImage[ 1 ] + ((FLTNB)ap_iDimVox[ 1 ]) * ap_iSizeVox[ 1 ]
      + ap_iSizeVox[ 1 ] * ((FLTNB)0.5),
    posOldImage[ 2 ] + ((FLTNB)ap_iDimVox[ 2 ]) * ap_iSizeVox[ 2 ]
      + ap_iSizeVox[ 2 ] * ((FLTNB)0.5)
  };
  
  // Computing and storing the position of the center of the pixel in each
  // direction for the initial image
  FLTNB **pOldCoordCenter = new FLTNB*[ 3 ];
  
  for( uint32_t i = 0; i < 3; ++i )
  {
    pOldCoordCenter[ i ] = new FLTNB[ iDimVoxPad[ i ] ];
    // Set the values
    for( uint32_t j = 0; j < iDimVoxPad[ i ]; ++j )
    {
      pOldCoordCenter[ i ][ j ] = posOldImage[ i ] - ap_iSizeVox[ i ] / 2.0
        + j * ap_iSizeVox[ i ];
    }
  }

  // Computing and storing the position of the center of the pixel in each
  // direction of the resampled image
  FLTNB **pNewCoordCenter = new FLTNB*[ 3 ];
  for( uint32_t i = 0; i < 3; ++i )
  {
    pNewCoordCenter[ i ] = new FLTNB[ ap_oDimVox[ i ] ];
    // Set the values
    for( uint32_t j = 0; j < ap_oDimVox[ i ]; ++j )
    {
      pNewCoordCenter[ i ][ j ] = posNewImage[ i ] + ap_oSizeVox[ i ] / 2.0
        + j * ap_oSizeVox[ i ];
    }
  }

  // Defining parameters
  FLTNB const invSizeX = 1.0 / ap_iSizeVox[ 0 ];
  FLTNB const invSizeY = 1.0 / ap_iSizeVox[ 1 ];
  FLTNB const invSizeZ = 1.0 / ap_iSizeVox[ 2 ];

  // Allocating memory for the 8 pixels kernel
  FLTNB *pKernelData = new FLTNB[ 8 ];

  // Loop over the elements of the new images
  for( uint32_t k = 0; k < ap_oDimVox[ 2 ]; ++k )
  {
    // Get the coordinate in Z
    FLTNB const z = pNewCoordCenter[ 2 ][ k ];
    if( z < boundMin[ 2 ] || z > boundMax[ 2 ] ) continue;

    // Find the bin in the old image
    int32_t const zBin = ( z - boundMin[ 2 ] ) * invSizeZ;

    // Computing the 2 z composants
    FLTNB const zComposantI0 = invSizeZ * ( pOldCoordCenter[ 2 ][ zBin + 1 ] - z );
    FLTNB const zComposantI1 = invSizeZ * ( z - pOldCoordCenter[ 2 ][ zBin ] );

    for( uint32_t j = 0; j < ap_oDimVox[ 1 ]; ++j )
    {
      // Get the coordinate in Y
      FLTNB const y = pNewCoordCenter[ 1 ][ j ];
      if( y < boundMin[ 1 ] || y > boundMax[ 1 ] ) continue;

      // Find the bin in the old image
      int32_t const yBin = ( y - boundMin[ 1 ] ) * invSizeY;

      // Computing the 2 y composants
      FLTNB const yComposantI0 = invSizeY * ( pOldCoordCenter[ 1 ][ yBin + 1 ]
        - y );
      FLTNB const yComposantI1 = invSizeY * ( y
        - pOldCoordCenter[ 1 ][ yBin ] );

      for( uint32_t i = 0; i < ap_oDimVox[ 0 ]; ++i )
      {
        // Get the coordinate in X
        FLTNB const x = pNewCoordCenter[ 0 ][ i ];
        if( x < boundMin[ 0 ] || x > boundMax[ 0 ] ) continue;

        // Find the bin in the old image
        int32_t const xBin = ( x - boundMin[ 0 ] ) * invSizeX;

        // Computing the 2 x composants
        FLTNB const xComposantI0 = invSizeX * (
          pOldCoordCenter[ 0 ][ xBin + 1 ] - x );
        FLTNB const xComposantI1 = invSizeX * ( x
          - pOldCoordCenter[ 0 ][ xBin ] );

        // TODO : check for potential segfaults 
        // Fill the buffer storing the data
        for( uint32_t kk = 0; kk < 2; ++kk )
        {
          for( uint32_t jj = 0; jj < 2; ++jj )
          {
            for( uint32_t ii = 0; ii < 2; ++ii )
            {
              pKernelData[ ii + jj * 2 + kk * 2 * 2 ] =
                pPadIData[
                  ( xBin + ii ) +
                  ( yBin + jj ) * iDimVoxPad[ 0 ] +
                  ( zBin + kk ) * iDimVoxPad[ 0 ] * iDimVoxPad[ 1 ]
                ];
            }
          }
        }

        // Computing the interpolation
        // In X
        FLTNB const xInterpVal0 = pKernelData[ 0 ] * xComposantI0 +
          pKernelData[ 1 ] * xComposantI1;

        FLTNB const xInterpVal1 = pKernelData[ 2 ] * xComposantI0 +
          pKernelData[ 3 ] * xComposantI1;

        FLTNB const xInterpVal2 = pKernelData[ 4 ] * xComposantI0 +
          pKernelData[ 5 ] * xComposantI1;

        FLTNB const xInterpVal3 = pKernelData[ 6 ] * xComposantI0 +
          pKernelData[ 7 ] * xComposantI1;

        // In Y
        FLTNB const yInterpVal0 = xInterpVal0 * yComposantI0 +
          xInterpVal1 * yComposantI1;

        FLTNB const yInterpVal1 = xInterpVal2 * yComposantI0 +
          xInterpVal3 * yComposantI1;

        // Final interpolation
        FLTNB const interpValTot = yInterpVal0 * zComposantI0 +
          yInterpVal1 * zComposantI1;

        ap_oImg[ i + j * ap_oDimVox[ 0 ]
          + k * ap_oDimVox[ 0 ] * ap_oDimVox[ 1 ] ] = interpValTot;
      }
    }
  }

  // Freeing the memory
  for( uint32_t i = 0; i < 3; ++i )
  {
    delete[] pOldCoordCenter[ i ];
    delete[] pNewCoordCenter[ i ];
  }
  delete[] pOldCoordCenter;
  delete[] pNewCoordCenter;
  delete[] pPadIData;
  delete[] pKernelData;
  
  return 0;


}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteImgFile
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
int IntfWriteImgFile(const string& a_pathToImg, FLTNB* ap_ImgMatrix, const Intf_fields& ap_IntfF, int vb)
{
  if (vb>=3) Cout("IntfWriteImgFile (with Intf_fields)" << endl);

  // Write Interfile header
  if(IntfWriteHeaderMainData(a_pathToImg, ap_IntfF, vb) )
  {
    Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile header for'"<< a_pathToImg << "' !" << endl);
    return 1;
  }
  
  string path_to_image = a_pathToImg;
  path_to_image.append(".img");
  
  // Read Interfile image
  if(IntfWriteImage(path_to_image, ap_ImgMatrix, ap_IntfF.mtx_size[0] * ap_IntfF.mtx_size[1] * ap_IntfF.mtx_size[2], vb) )
  {
    Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile image '"<< path_to_image << "' !" << endl);
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteImgFile
  \param a_pathToImg : path to image basename
  \param ap_ImgMatrix : 1 dimensional image matrix which contains the image to write
  \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
  \param vb : verbosity
  \brief Main function dedicated to Interfile 3D image writing
  \details Call the main functions dedicated to Interfile header and data writing :
           IntfWriteHeaderMainData() and then IntfWriteImage()
  \todo Get metadata from a Intf_fields object ?
       (useful to transfer keys from read images to written images)
  \return 0 if success, positive value otherwise.
*/
int IntfWriteImgFile(const string& a_pathToImg, FLTNB* ap_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb)
{
  if (vb>=3) Cout("IntfWriteImgFile (3D image)" << endl);
    
  Intf_fields Img_fields;
  IntfKeySetFieldsOutput(&Img_fields, ap_ID);

  // Write Interfile header
  if(IntfWriteHeaderMainData(a_pathToImg, Img_fields, vb) )
  {
    Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile header for'"<< a_pathToImg << "' !" << endl);
    return 1;
  }
  
  string path_to_image = a_pathToImg;
  path_to_image.append(".img");
  
  // Read Interfile image
  if(IntfWriteImage(path_to_image, ap_ImgMatrix, ap_ID->GetNbVoxXYZ(), vb) )
  {
    Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile image '"<< path_to_image << "' !" << endl);
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteProjFile
  \param a_pathToImg : string containing the path to the image basename
  \param a2p_ImgMatrix : 2 dimensional image matrix which contains the image to write.
  \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
  \param a_Imgfields: Structure containing information about the projected data
  \param vb : verbosity
  \brief Function dedicated to Interfile image writing for projected data
  \details Call the main functions dedicated to Interfile header and data writing
           Currently work for SPECT projected data
           The total number of projections should be provided in parameters
           Depending on the output writing mode (stored in sOutputManager),
  \todo Get metadata from a Intf_fields object ?
       (useful to transfer keys from read images to written images)
  \return 0 if success, positive value otherwise.
*/
int IntfWriteProjFile(const string& a_pathToImg, FLTNB** a2p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, Intf_fields a_Imgfields, int vb)
{
  if (vb>=3) Cout("IntfWriteProjFile ..." << endl);

  //Set image names
  vector<string> lpath_to_images;
  lpath_to_images.push_back(a_pathToImg);
  
  if(IntfWriteHeaderMainData(a_pathToImg, a_Imgfields, vb) )
  {
    Cerr("***** IntfWriteProjFile()-> Error : while trying to write the interfile header for'"<< lpath_to_images[0] << "' !" << endl);
    return 1;
  }
      
  // Write interfile image

  // Set dimensions
  uint32_t dims[2];
  dims[0] = a_Imgfields.nb_projections;
  dims[1] = a_Imgfields.mtx_size[0]*a_Imgfields.mtx_size[1];
  
  if(IntfWriteImage(lpath_to_images, a2p_ImgMatrix, dims, vb) )
  {
    Cerr("***** IntfWriteFile()-> Error : while trying to write the interfile image '"<< a_pathToImg << "' !" << endl);
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteImgDynCoeffFile
  \param a_pathToImg : string containing the path to the image basename
  \param a2p_ImgMatrix : 2 dimensional image matrix which contains the image to write.
  \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
  \param a_nbParImgs : Number of parametric images
  \param vb : verbosity
  \brief Function dedicated to Interfile image writing for dynamic coefficients images
  \details Call the main functions dedicated to Interfile header and data writing
           The total number of basis functions should be provided in parameters
           Depending on the output writing mode (stored in sOutputManager),
           One image with one file and one will be created for the whole dynamic image
           or a metaheader and associated multiple 3D image raw file/headers will be generated
  \todo Get metadata from a Intf_fields object ?
       (useful to transfer keys from read images to written images)
  \return 0 if success, positive value otherwise.
*/
int IntfWriteImgDynCoeffFile(const string& a_pathToImg, FLTNB** a2p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int a_nbParImgs, int vb)
{
  if (vb>=3) Cout("IntfWriteImgDynCoeffFile  ..." << endl);
    
  Intf_fields Img_fields; //TODO keep refs to Intf_fields in case we need it later (transfer keys from read images to written images)
  IntfKeySetFieldsOutput(&Img_fields, ap_ID);
  Img_fields.nb_time_frames = a_nbParImgs;
  Img_fields.nb_resp_gates = Img_fields.nb_card_gates = 1;
  
  // Get input from user about how dynamic files should be written 
  // (one global image file, or separate image file for each frame/gates).
  sOutputManager* p_outputMgr = sOutputManager::GetInstance();
  
  //Set image names
  vector<string> lpath_to_images;
  
  if(p_outputMgr->MergeDynImages() == true)   // Case one file
  {
    lpath_to_images.push_back(a_pathToImg);
    
    if(IntfWriteHeaderMainData(lpath_to_images[0], Img_fields, vb) )
    {
      Cerr("***** IntfWriteImgDynCoeffFile()-> Error : while trying to write the interfile header for'"<< lpath_to_images[0] << "' !" << endl);
      return 1;
    }
  }
  else   // Case several files
  {
    int idx_file = 0;
    
    for(int bf=0 ; bf<a_nbParImgs ; bf++)
    {
      //lpath_to_images[idx_file] = a_pathToImg;
      lpath_to_images.push_back(a_pathToImg);
      
      // Add a suffix for basis functions
      stringstream ss; ss << bf + 1;
      lpath_to_images[idx_file].append("_par").append(ss.str());
      
      // Write header specific to this image
      string path_to_hdr = lpath_to_images[idx_file] + ".hdr";
      string path_to_img = lpath_to_images[idx_file] + ".img";

      // As the content will be appended, make sure there is no existing file with such name
      // remove it otherwise
      ifstream fcheck(path_to_hdr.c_str());
      if(fcheck.good())
      {
        fcheck.close();
        #ifdef _WIN32
        string dos_instruction = "del " + path_to_hdr;
        system(dos_instruction.c_str());
        #else
        remove(path_to_hdr.c_str());
        #endif
      }
            
      ofstream ofile(path_to_hdr.c_str(), ios::out | ios::app);
      ofile << setprecision(std::numeric_limits<FLTNB>::digits10 +1);
      ofile << "!INTERFILE := " << endl; 
      ofile << "!name of data file := " << GetFileFromPath(path_to_img) << endl;
      ofile << "imagedata byte order := " << IntfKeyGetEndianStr(User_Endianness) << endl;
      ofile << "!data offset in bytes := " << 0 << endl;
      ofile << "patient name:= " << GetFileFromPath(path_to_img) << endl;
      ofile << "!type of data := Dynamic" << endl; // /!\ Dynamic by default in order to be readable from ImageJ OpenMed plugin
    
      if(IntfWriteHeaderImgData(ofile, Img_fields, vb) )
      {
        Cerr("***** IntfWriteImgDynCoeffFile()-> Error : while trying to write the interfile header for'"<< lpath_to_images[idx_file] << "' !" << endl);
        return 1;
      }
      
      ofile << "!END OF INTERFILE := " << endl; 
      
      ofile.close();
      idx_file++;
    }
    
    // Write meta header
    if(IntfWriteMHD(a_pathToImg, lpath_to_images, Img_fields, ap_ID, vb) )
    {
      Cerr("***** IntfWriteImgDynCoeffFile()-> Error : while trying to write the interfile header '"<< a_pathToImg << "' !" << endl);
      return 1;
    }
  }
  
  // Write interfile image

  // Set dimensions
  uint32_t dims[2];
  dims[0] = a_nbParImgs;
  dims[1] = ap_ID->GetNbVoxXYZ();
  
  if(IntfWriteImage(lpath_to_images, a2p_ImgMatrix, dims, vb) )
  {
    Cerr("***** IntfWriteImgDynCoeffFile()-> Error : while trying to write the interfile image '"<< a_pathToImg << "' !" << endl);
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteImgFile
  \param a_pathToImg : string containing the path to the image basename
  \param a4p_ImgMatrix : 4 dimensional image matrix which contains the image to write.
  \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
  \param vb : verbosity
  \brief Main function dedicated to Interfile 6D (dynamic dims + 3D ) image writing
  \details Call the main functions dedicated to Interfile header and data writing 
           Depending on the output writing mode (stored in sOutputManager),
           One image with one file and one will be created for the whole dynamic image
           or a metaheader and associated multiple 3D image raw file/headers will be generated
  \todo Get metadata from a Intf_fields object ?
       (useful to transfer keys from read images to written images)
  \return 0 if success, positive value otherwise.
*/
int IntfWriteImgFile(const string& a_pathToImg, FLTNB**** a4p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb)
{
  if (vb>=3) Cout("IntfWriteImgFile (static/dynamic image) ..." << endl);
    
  Intf_fields Img_fields; //TODO keep refs to Intf_fields in case we need it later (transfer keys from read images to written images)
  IntfKeySetFieldsOutput(&Img_fields, ap_ID);


  // Get input from user about how dynamic files should be written 
  // (one global image file, or separate image file for each frame/gates).
  sOutputManager* p_outputMgr = sOutputManager::GetInstance();
  
  //Set image names
  vector<string> lpath_to_images;  

  if(Img_fields.nb_dims <=3 || p_outputMgr->MergeDynImages() == true)   // One file to write
  {
    lpath_to_images.push_back(a_pathToImg);
    
    if(IntfWriteHeaderMainData(a_pathToImg, Img_fields, vb) )
    {
      Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile header for'"<< lpath_to_images[0] << "' !" << endl);
      return 1;
    }
  }
  else   // Case several files
  {
     
    int idx_file = 0;
    
    for(int fr=0 ; fr<ap_ID->GetNbTimeFrames() ; fr++)
      for(int rg=0 ; rg<ap_ID->GetNbRespGates() ; rg++)
        for(int cg=0 ; cg<ap_ID->GetNbCardGates() ; cg++)
        {
          
          //lpath_to_images[idx_file] = a_pathToImg;
          
          lpath_to_images.push_back(a_pathToImg);
          
          // Add a suffix for time frames if more than 1
          if (ap_ID->GetNbTimeFrames() > 1)
          {
            stringstream ss; ss << fr + 1;
            lpath_to_images[idx_file].append("_fr").append(ss.str());
          }
          
          // Add a suffix for respiratory gates if more than 1
          if (ap_ID->GetNbRespGates() > 1)
          {
            stringstream ss; ss << rg + 1;
            lpath_to_images[idx_file].append("_rg").append(ss.str());
          }
          // Add a suffix for cardiac gates if more than 1
          if (ap_ID->GetNbCardGates() > 1)
          {
            stringstream ss; ss << cg + 1;
            lpath_to_images[idx_file].append("_cg").append(ss.str());
          }

          // Write header specific to this image
          string path_to_hdr = lpath_to_images[idx_file];
          string path_to_img = lpath_to_images[idx_file];

          // As the content will be appended, make sure there is no existing file with such name
          // remove it otherwise
          ifstream fcheck(path_to_hdr.append(".hdr").c_str());
          if(fcheck.good())
          {
            fcheck.close();
            #ifdef _WIN32
            string dos_instruction = "del " + path_to_hdr;
            system(dos_instruction.c_str());
            #else
            remove(path_to_hdr.c_str());
            #endif
          }
    
          ofstream ofile(path_to_hdr.c_str(), ios::out | ios::app);
          ofile << setprecision(std::numeric_limits<FLTNB>::digits10 +1);
          ofile << "!INTERFILE := " << endl; 
          
          
          ofile << "!name of data file := " << GetFileFromPath(path_to_img).append(".img") << endl;
          ofile << "imagedata byte order := " << IntfKeyGetEndianStr(User_Endianness) << endl;
          ofile << "!data offset in bytes := " << 0 << endl;
          ofile << "patient name:= " << GetFileFromPath(path_to_img) << endl;
          ofile << "!type of data := Dynamic" << endl; // /!\ Dynamic by default in order to be readable from ImageJ OpenMed plugin
      
          if(IntfWriteHeaderImgData(ofile, Img_fields, vb) )
          {
            Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile header for'"<< lpath_to_images[idx_file] << "' !" << endl);
            return 1;
          }

          // Write image duration related data (this information will be replicated on each gated image) 
          ofile << "image duration (sec) := " << (ap_ID->GetFrameTimeStopInSec(0, fr) - ap_ID->GetFrameTimeStartInSec(0, fr) )<< endl;
          
          // same start time for all gates
          ofile << "image start time (sec) := " << ap_ID->GetFrameTimeStartInSec(0, fr) << endl; 

          ofile << "!END OF INTERFILE := " << endl; 

          ofile.close();

          idx_file++;
          
        }
    
    // Write meta header
    if(IntfWriteMHD(a_pathToImg, lpath_to_images, Img_fields, ap_ID, vb) )
    {
      Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile header '"<< a_pathToImg << "' !" << endl);
      return 1;
    }
     
  }

  // Write interfile image

  // Set dimensions
  uint32_t dims[4];
  dims[0] = ap_ID->GetNbTimeFrames();
  dims[1] = ap_ID->GetNbRespGates();
  dims[2] = ap_ID->GetNbCardGates();
  dims[3] = ap_ID->GetNbVoxXYZ();

  if(IntfWriteImage(lpath_to_images, a4p_ImgMatrix, dims, vb) )
  {
    Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile image '"<< a_pathToImg << "' !" << endl);
    return 1;
  }
  
  return 0;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteImgFile
  \param a_pathToImg : string containing the path to the image basename
  \param a4p_ImgMatrix : 4 dimensional image matrix which contains the image to write.
  \param ap_ID : Provide the Image dimensions object containing reconstruction dimensions
  \param vb : verbosity
  \brief Main function dedicated to Interfile 6D (dynamic dims + 3D ) image writing
  \details Call the main functions dedicated to Interfile header and data writing 
           Depending on the output writing mode (stored in sOutputManager),
           One image with one file and one will be created for the whole dynamic image
           or a metaheader and associated multiple 3D image raw file/headers will be generated
  \todo Get metadata from a Intf_fields object ?
       (useful to transfer keys from read images to written images)
  \return 0 if success, positive value otherwise.

int IntfWriteSensitivityImg(const string& a_pathToImg, FLTNB**** a4p_ImgMatrix, oImageDimensionsAndQuantification* ap_ID, int vb)
{
  if(vb >= 2) Cout("IntfWriteSensitivityImg ..." << endl);
    
  Intf_fields Img_fields; //TODO keep refs to Intf_fields in case we need it later (transfer keys from read images to written images)
  IntfKeySetFieldsOutput(&Img_fields, ap_ID);

  // Get input from user about how dynamic files should be written 
  // (one global image file, or separate image file for each frame/gates).
  sOutputManager* p_outputMgr = sOutputManager::GetInstance();
  
  //Set image names
  vector<string> lpath_to_images;  
  
  // By default, sensitivity image are written in one file One file to write
  {
    lpath_to_images.push_back(a_pathToImg);
    
    if(IntfWriteHeaderMainData(a_pathToImg, Img_fields, ap_ID, vb) )
    {
      Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile header for'"<< lpath_to_images[0] << "' !" << endl);
      return 1;
    }
  }

  // Write interfile image

  // Set dimensions
  uint32_t dims[4];
  dims[0] = ap_ID->GetNbTimeFrames();
  dims[1] = ap_ID->GetNbRespGates();
  dims[2] = ap_ID->GetNbCardGates();
  dims[3] = ap_ID->GetNbVoxXYZ();

  if(IntfWriteImage(lpath_to_images, a4p_ImgMatrix, dims, vb) )
  {
    Cerr("***** IntfWriteImgFile()-> Error : while trying to write the interfile image '"<< a_pathToImg << "' !" << endl);
    return 1;
  }
  
  return 0;
}
*/


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfIsMHD
  \param a_pathToFile : String containing path to an Interfile header
  \param ap_lPathToImgs : pointer to a list of strings containing the path to the different images
  \brief Check if the string in argument contains the path to a Interfile metaheader
  \details Check for '!total number of data sets' key, and the associated image file names
           If the file is metaheader, the names of image file header will be returned in ap_lPathToImgs list
           It not, the file is a single header, that we add to the list as solo file
  \return 1 if we have a metaheader, 0 if not, negative value if error.
*/
int IntfIsMHD(string a_pathToFile, vector<string> &ap_lPathToImgs)
{
  // Create file object and check existence
  ifstream hfile(a_pathToFile.c_str(), ios::in);
  string line;
  
  if(hfile)
  {
    // Read until the end of file and check if we have a '!total number of data sets' key
    while(!hfile.eof())
    {
      getline(hfile, line);

      // Check if empty line.
      if(line.empty()) continue;
      
      // Initiate a Key which will recover the potential interfile data for each line
      Intf_key Key;
      
      // Process the Key at this line. Recover the Key, its counterpart in lower case, and the value in the structure Intf_key
      if (IntfRecoverKey(&Key, line) ) 
      {
        Cerr("***** IntfIsMHD()-> Error : couldn't correctly read interfile key '"<< line << "' !" << endl);
        return -1;
      }

      // --- Check if a total number of data set if provided
      if(IntfCheckKeyMatch(Key, "total number of data set")  ||
         IntfCheckKeyMatch(Key, "total number of data sets") ||
         IntfCheckKeyMatch(Key, "!total number of data set") ||
         IntfCheckKeyMatch(Key, "!total number of data sets")) 
      {
        //  File is MHD
        
        uint16_t nb_datasets = 0;
        if(ConvertFromString(Key.kvalue , &nb_datasets) )
        {
          Cerr("***** IntfIsMHD()-> Error when trying to read data from 'total number of data set' key. Recovered value = " << Key.kvalue << endl);
          return -1;
        }
    
        // Get the image file path from the following lines
        int file_idx = 0;
        stringstream ss;
        ss << (file_idx+1);
        string file_key = "";
        string file_keyp = "%";
        string file_key_space = "";
        string file_keyp_space = "%";
        file_key_space.append("data set [").append(ss.str()).append("]");
        file_keyp_space.append("data set [").append(ss.str()).append("]");
        file_key.append("data set[").append(ss.str()).append("]");
        file_keyp.append("data set[").append(ss.str()).append("]");
        
        string path_to_mhd_dir = GetPathOfFile(a_pathToFile);

        while(!hfile.eof())
        {
          // Read the following lines...
          getline(hfile, line);
          
          // Check if empty line.
          if(line.empty()) continue;
          
          // Process the Key at this line. Recover the Key, its counterpart in lower case, and the value in the structure Intf_key
          if (IntfRecoverKey(&Key, line) ) 
          {
            Cerr("***** IntfIsMHD()-> Error : couldn't correctly read interfile key '"<< line << "' !" << endl);
            return -1;
          }

          if(IntfCheckKeyMatch(Key, file_key)  || 
             IntfCheckKeyMatch(Key, file_keyp) ||
             IntfCheckKeyMatch(Key, file_key_space) ||
             IntfCheckKeyMatch(Key, file_keyp_space))
          {
            ap_lPathToImgs.push_back(GetPathOfFile(a_pathToFile));
            if(IntfKeyIsArray(Key))
            {
              string elts_str[3];
              
              // First, check we have the correct number of elts
              int nb_elts = IntfKeyGetArrayNbElts(Key);
              if(nb_elts<0)
              {
                Cerr("***** IntfIsMHD()-> Error : when trying to read Interfile array key: '"<< line << "' !" << endl);
                return -1;
              }
              else if(nb_elts != 3)
              {
                Cerr("***** IntfIsMHD()-> Error : when trying to read Interfile array key: '"<< line << "' !" << endl);
                Cerr("                    3 elements are expected following the format {xxx, path_to_the_img_file, xxx} (xxx = ignored data)" << endl);
                return -1;
              }
              
              if(IntfKeyGetArrayElts(Key, elts_str) )
              {
                Cerr("***** IntfIsMHD()-> Error : when trying to read Interfile array key: '"<< line << "' !" << endl);
                return -1;
              }
              ap_lPathToImgs.at(file_idx).append(elts_str[1]);
            }
            else
              ap_lPathToImgs.at(file_idx).append(Key.kvalue);
            
            // Compute next key to recover
            file_idx++;
            stringstream ss;
            ss << (file_idx+1);
            file_key = "";
            file_keyp = "%";
            file_key.append("data set [").append(ss.str()).append("]");
            file_keyp.append("data set [").append(ss.str()).append("]");
            
          }
        }
    
        // Check we have a correct nb of images
        if(nb_datasets != ap_lPathToImgs.size())
        {
          // error if not found
          Cerr("***** IntfIsMHD()-> The number of recovered file in the metaheader ("<<ap_lPathToImgs.size()<<")");
          Cerr("does not correspond to the expected number of datasets ("<<nb_datasets<<")!"<< endl);
          return -1;
        }
        else
          // Normal end, file is metaheader
          return 1;
      }
    }
  }
  else
  {
    Cerr("***** IntfIsMHD()-> Error : couldn't read header file '"<< a_pathToFile << "' !" << endl);
    return -1;
  }

  // If we reach this, the file is NOT considered as metaheader the 'total number of data set' key was not found
  // Add the to the list as unique header file
  ap_lPathToImgs.push_back(a_pathToFile);
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteMHD
  \param a_path : path to the Meta interfile header to write
  \param ap_lPathToImgs : pointer to a list of strings containing the path to the different images
  \param a_IntfF : Structure containing interfile keys of the image to write
  \param ap_ID : oImageDimensionsAndQuantification object containing additional infos about reconstruction
  \param vb : verbosity
  \brief Write an Interfile meta header at the path provided in parameter, using the field stack provided in parameter.
  \todo keys for multi bed
  \return 0 if success, positive value otherwise.
*/
int IntfWriteMHD(    const string& a_pathToMhd, 
             const vector<string>& ap_lPathToImgs, 
                       Intf_fields a_IntfF, 
oImageDimensionsAndQuantification* ap_ID, 
                               int vb)
{  
  // Get file and check existence
  string path_to_mhd_file = a_pathToMhd;
  path_to_mhd_file.append(".mhd");
  ofstream ofile(path_to_mhd_file.c_str(), ios::out);
  ofile << setprecision(std::numeric_limits<FLTNB>::digits10 +1);

  if(ofile)
  {
    sScannerManager* p_scanMgr = sScannerManager::GetInstance(); 

    // todo fields for multi-beds ?
    
    ofile << "!INTERFILE := " << endl; 
    //ofile << "% MetaHeader written by CASToRv1.0 := " << end;
    ofile << endl << "!GENERAL DATA := " << endl; 
    ofile << "data description:=image" << endl;
    ofile << "!originating system := " << p_scanMgr->GetScannerName() << endl;

    // Number of bed positions
    if (ap_ID->GetNbBeds()>1) ofile << "number of bed positions := " << ap_ID->GetNbBeds() << endl;
    // Relative position of each bed
    else ofile << "horizontal bed relative position (mm) := " << ap_ID->GetBedPosition(0) << endl;

    ofile << endl << "%DATA MATRIX DESCRIPTION:=" << endl; 
    ofile << "number of time frames := " << a_IntfF.nb_time_frames << endl;
    ofile << "number of time windows := " << a_IntfF.nb_resp_gates *
                                                   a_IntfF.nb_card_gates << endl;
    ofile << "number of respiratory gates := " << a_IntfF.nb_resp_gates << endl;
    ofile << "number of cardiac gates := " << a_IntfF.nb_card_gates << endl;
  
    ofile << endl << "%DATA SET DESCRIPTION:="<< endl;
    
    uint16_t nb_imgs = a_IntfF.nb_time_frames * a_IntfF.nb_resp_gates * a_IntfF.nb_card_gates ; 
    ofile << "!total number of data sets:=" << nb_imgs << endl;
    
    // Check we have the right number of strings
    if(ap_lPathToImgs.size() != nb_imgs)
    {
      Cerr("***** IntfWriteMHD()-> Error : nb of provided string inconsistent with the expected number of dynamic images: '"<< nb_imgs << "' !" << endl);
      ofile.close();
      return 1;
    }
  
    for(int ds=0 ; ds<nb_imgs ; ds++)
    {
      string path_to_header = ap_lPathToImgs.at(ds);
      ofile << "%data set ["<<ds+1<<"]:={0," << GetFileFromPath(path_to_header).append(".hdr") << ",UNKNOWN}"<< endl;
    }
  }
  else
  {
    Cerr("***** IntfWriteMHD()-> Error : couldn't find output Metaheader interfile '"<< path_to_mhd_file << "' !" << endl);
    return 1;
  }
    
  ofile.close();
  

    
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfReadHeader
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
int IntfReadHeader(const string& a_pathToHeaderFile, Intf_fields* ap_IF, int vb)
{
  if (vb >= 4)
  {
    Cout("--------------------------------------------------------------- " << endl);
    Cout("IntfReadHeader()-> Start reading header interfile " << a_pathToHeaderFile << endl);
    Cout("--------------------------------------------------------------- " << endl);
  }


  // Get file and check existence
  ifstream input_file(a_pathToHeaderFile.c_str(), ios::in);
  string line;
  
  if(input_file)
  {
    // Read until the end of file
    while(!input_file.eof())
    {
      getline(input_file, line);

      // Check if empty line.
      if(line.empty())
        continue;
      
      // Initiate a Key which will recover the potential interfile data for each line
      Intf_key Key;
      
      // Process the Key at this line. Recover the Key, its counterpart in lower case, and the value in the structure Intf_key
      if (IntfRecoverKey(&Key, line) ) 
      {
        Cerr("***** ReadIntfHeader()-> Error : couldn't correctly read interfile key '"<< line << "' !" << endl);
        return 1;
      }
      
      if (vb >=10) Cout("ReadIntfHeader()-> Key  " << Key.kcase << endl);


      // --- From there, find key Id and process result --- 
      
      // Check version. Warning message if not 3.3 / Castor keys ?
      if(IntfCheckKeyMatch(Key, "version of keys")) 
      {
        //if(Key.kvalue != 3.3 ||
        //   Key.kvalue != CASToRv1.0)
        //  Cout("***** ReadIntfHeader()-> Warning : Unexpected version of Interfile keys found: " << Key.kvalue << " (supported version = 3.3)" << endl);
      }
      
      
      // --- path to image file
      // "name of data file" / Intf_fields.path_to_image
      else if(IntfCheckKeyMatch(Key, "name of data file")) 
      {
        // Check key isn't empty, throw error otherwise
        if ( Key.kvalue.size()>0 ) 
        {
          // Look for path separator
          if ( Key.kvalue.find(OS_SEP) != string::npos ) 
          { 
            // Assuming we have an absolute path 
            ap_IF->path_to_image = Key.kvalue;
          }
          else
          {
            // Assuming we just have the image name 
            // -> use the relative path of the header file, image should be in the same dir
            string header_file_directory = GetPathOfFile(a_pathToHeaderFile); 
            ap_IF->path_to_image = header_file_directory.append(Key.kvalue); 
          } 
        }
        else
        {
          Cerr("***** ReadIntfHeader()-> Error : path to interfile image file is empty !" << endl);
          return 1;
        }
      }
      
      
      // --- Little or big endian
      // "imagedata byte order" / Intf_fields.endianness
      else if(IntfCheckKeyMatch(Key, "imagedata byte order")) 
      {
        string endianness;
        IntfToLowerCase(&Key.kvalue); // whole string in low caps
        
        endianness = Key.kvalue;
        
        if (endianness == "littleendian")
          ap_IF->endianness = 1;
        else // Big endian by default
          ap_IF->endianness = INTF_BIG_ENDIAN;
      }
      
      
      // --- Data offset using "data starting block"
      // "data starting block" & "data offset in bytes" / Intf_fields.data_offset
      else if(IntfCheckKeyMatch(Key, "data starting block")) 
      {
        // Throw warning indicating that the data offset has already been set
        if(ap_IF->data_offset>0 )
        {
          Cerr("***** ReadIntfHeader()-> Warning : data_offset value was already set to " << ap_IF->data_offset << endl);
          Cerr("*****                              The header may contain both 'data offset in bytes' and 'data starting  block' fields " << endl);
        }
        
        if(ConvertFromString(Key.kvalue , &ap_IF->data_offset) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'data starting block' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        
        ap_IF->data_offset*= 2048L;
      }
      

      // --- Data offset using "data offset in bytes"
      // "data starting block" & "data offset in bytes" / Intf_fields.data_offset
      else if(IntfCheckKeyMatch(Key, "data offset in bytes")) 
      {
        // Throw warning indicating that the data offset has already been set
        if(ap_IF->data_offset>0 )
        {
          Cerr("***** ReadIntfHeader()-> Warning : data_offset value was already set to " << ap_IF->data_offset << endl);
          Cerr("*****                              The header may contain both 'data offset in bytes' and 'data starting  block' fields " << endl);
        }
    
        if(ConvertFromString(Key.kvalue , &ap_IF->data_offset) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'data offset in bytes' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      
      
      // --- Type of each pixel/voxel data
      // "number format" / Intf_fields.nb_format
      else if(IntfCheckKeyMatch(Key, "number format"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->nb_format) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number format' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }

      // --- Originating system (the scanner)
      // "originating system" / Intf_fields.originating_system
      else if(IntfCheckKeyMatch(Key, "originating system"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->originating_system) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'originating system' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }

      // --- Bed relative position in mm
      // "bed relative position (mm)" / Intf_fields.bed_relative_position
      else if(IntfCheckKeyMatch(Key, "horizontal bed relative position (mm)"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->bed_relative_position) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'bed relative position (mm)' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      
      // --- Width
      // "matrix size [1]" / Intf_fields.mtx_size[0]
      else if(IntfCheckKeyMatch(Key, "matrix size [1]"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->mtx_size[0]) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'matrix size [1]' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }

      // --- Height
      // "matrix size [2]" / Intf_fields.mtx_size[1]
      else if(IntfCheckKeyMatch(Key, "matrix size [2]"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->mtx_size[1]) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'matrix size [2]' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      
      // --- Slices : Tomographic image
      // "matrix size [3]" / Intf_fields.mtx_size[2]
      else if(IntfCheckKeyMatch(Key, "matrix size [3]"))
      {
        // Check if we have an array key
        if (IntfKeyIsArray(Key)) 
        {
          //TODO : perhaps more convenient to use vector to deal with data whose matrix size [3] has more than 1 element (eg: sinograms)
          //       It will depend how we decide to handle sinogram reading
          //       For now, throw error by default as there is no implementation yet

          // Throw error by default
          Cerr("***** ReadIntfHeader()-> Array reading for 'matrix size [3]' key not implemented yet. Single value expected." << endl);
          return 1;
        }
        else
        {
          if(ConvertFromString(Key.kvalue , &ap_IF->mtx_size[2]) )
          {
            Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'matrix size [3]' key. Recovered value = " << Key.kvalue  << endl);
            return 1;
          }
        }
      }
      
      // --- Read nb frames or other.
      // "matrix size [4]" / Intf_fields.mtx_size[3]
      // Consistency with "number of time frames" : Priority to number of time frames if this has already been initialized
      // TODO : this may cause issue for implementation of sinogram management 
      else if(IntfCheckKeyMatch(Key, "matrix size [4]"))
      {
        // Could be sinogram or dynamic related (nb time frames).
        // Dynamic reconstructions : consistency with "number of frames"

        if (IntfKeyIsArray(Key)) 
        {
          //TODO : perhaps more convenient to use vector to deal with data whose matrix size [4] has more than 1 element (eg: sinograms)
          //       It will depend how we decide to handle sinogram reading
          //       For now, throw error by default as there is no implementation yet
          
          // Throw error by default
          Cerr("***** ReadIntfHeader()-> Array reading for 'matrix size [4]' key not implemented yet. Single value expected." << endl);
          return 1;
        }
        else
        {
          // Check if number time frames has already been initialized
          if(ap_IF->nb_time_frames>1)
          {
            ap_IF->mtx_size[3] = ap_IF->nb_time_frames;
            Cerr("***** ReadIntfHeader()-> WARNING : both 'number of time frames' and 'matrix size [4]' keys have been provided");
            Cerr("                         'number of time frames' selected by default" << endl);
          }
          else if( ConvertFromString(Key.kvalue , &ap_IF->mtx_size[3]) )
          {
            Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'matrix size [4]' key. Recovered value = " << Key.kvalue << endl);
            return 1;
          }
        }
      }


      // --- Related to sinogram (not implemented yet)
      // "matrix size [5]" / Intf_fields.mtx_size[4]
      else if(IntfCheckKeyMatch(Key, "matrix size [5]"))
      {
        // Could be sinogram or dynamic related
        if (IntfKeyIsArray(Key)) 
        {
          //TODO : perhaps more convenient to use vector to deal with data whose matrix size [5] has more than 1 element (eg: sinograms)
          //       It will depend how we decide to handle sinogram reading
          //       For now, throw error by default as there is no implementation yet
          
          // Throw error by default
          Cerr("***** ReadIntfHeader()-> Array reading for 'matrix size [5]' key not implemented yet. Single value expected." << endl);
          return 1;
        }
        else
        {
          if(ConvertFromString(Key.kvalue , &ap_IF->mtx_size[4]) )
          {
            Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'matrix size [5]' key. Recovered value = " << Key.kvalue << endl);
            return 1;
          }
        }
      }
      
      // --- Related to sinogram (not implemented yet)
      // "matrix size [6]" / Intf_fields.mtx_size[5]
      else if(IntfCheckKeyMatch(Key, "matrix size [6]"))
      {
        // Could be sinogram or dynamic related
        if (IntfKeyIsArray(Key)) 
        {
          //TODO : perhaps more convenient to use vector to deal with data whose matrix size [6] has more than 1 element (eg: sinograms)
          //       It will depend how we decide to handle sinogram reading
          //       For now, throw error by default as there is no implementation yet
          
          // Throw error by default
          Cerr("***** ReadIntfHeader()-> Array reading for 'matrix size [6]' key not implemented yet. Single value expected." << endl);
          return 1;
        }
        else
        {
          if(ConvertFromString(Key.kvalue , &ap_IF->mtx_size[5]) )
          {
            Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'matrix size [6]' key. Recovered value = " << Key.kvalue << endl);
            return 1;
          }
        }
      }
      
      // --- Related to sinogram (not implemented yet)
      // "matrix size [7]" / Intf_fields.mtx_size[6]
      else if(IntfCheckKeyMatch(Key, "matrix size [7]"))
      {
        // Could be sinogram related
        if (IntfKeyIsArray(Key)) 
        {
          //TODO : perhaps more convenient to use vector to deal with data whose matrix size [7] has more than 1 element (eg: sinograms)
          //       It will depend how we decide to handle sinogram reading
          //       For now, throw error by default as there is no implementation yet
          
          // Throw error by default
          Cerr("***** ReadIntfHeader()-> Array reading for 'matrix size [7]' key not implemented yet. Single value expected." << endl);
          return 1;
        }
        else
        {
          if( ConvertFromString(Key.kvalue , &ap_IF->mtx_size[6]) )
          {
            Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'matrix size [7]' key. Recovered value = " << Key.kvalue << endl);
            return 1;
          }
        }
        
      }
      
      // --- Size voxel width
      // "scaling factor (mm/pixel) [1]" / Intf_fields.vox_size[0]
      else if(IntfCheckKeyMatch(Key, "scaling factor (mm/pixel) [1]") 
           || IntfCheckKeyMatch(Key, "scale factor (mm/pixel) [1]"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->vox_size[0]) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'scaling factor (mm/pixel) [1]' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }

      // --- Size voxel height
      // "scaling factor (mm/pixel) [2]" / Intf_fields.vox_size[1]
      else if(IntfCheckKeyMatch(Key, "scaling factor (mm/pixel) [2]") 
           || IntfCheckKeyMatch(Key, "scale factor (mm/pixel) [2]"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->vox_size[1]) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'scaling factor (mm/pixel) [2]' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }

      // --- Size voxel axial
      // "scaling factor (mm/pixel) [3]" / Intf_fields.vox_size[2]
      // Prefered over slice thickness by default
      else if(IntfCheckKeyMatch(Key, "scaling factor (mm/pixel) [3]") 
           || IntfCheckKeyMatch(Key, "scale factor (mm/pixel) [3]")) 
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->vox_size[2]) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'scaling factor (mm/pixel) [3]' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }


      // --- Voxel offset width
      // "first pixel offset (mm) [1]" / Intf_fields.vox_offset[0]
      else if(IntfCheckKeyMatch(Key, "first pixel offset (mm) [1]") 
           || IntfCheckKeyMatch(Key, "origin (mm) [1]") 
           || IntfCheckKeyMatch(Key, "offset [1]")) 
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->vox_offset[0]) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from one of the following key :" << endl);
          Cerr("*****                    'first pixel offset (mm) [1]'" << Key.kvalue << endl);
          Cerr("*****                    'origin (mm) [1]'"<< endl);
          Cerr("*****                    'offset [1]'" << Key.kvalue << endl);
          Cerr("*****                    Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }


      // --- Voxel offset height
      // "first pixel offset (mm) [2]" / Intf_fields.vox_offset[1]
      else if(IntfCheckKeyMatch(Key, "first pixel offset (mm) [2]") 
           || IntfCheckKeyMatch(Key, "origin (mm) [2]")
           || IntfCheckKeyMatch(Key, "offset [2]")) 
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->vox_offset[1]) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from one of the following key :" << endl);
          Cerr("*****                    'first pixel offset (mm) [2]'" << Key.kvalue << endl);
          Cerr("*****                    'origin (mm) [2]'"<< endl);
          Cerr("*****                    'offset [2]'" << Key.kvalue << endl);
          Cerr("*****                    Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      
      
      // --- Voxel offset axial
      // "first pixel offset (mm) [3]" / Intf_fields.vox_offset[2]
      else if(IntfCheckKeyMatch(Key, "first pixel offset (mm) [3]") 
           || IntfCheckKeyMatch(Key, "origin (mm) [3]")
           || IntfCheckKeyMatch(Key, "offset [3]")) 
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->vox_offset[2]) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from one of the following key :" << endl);
          Cerr("*****                    'first pixel offset (mm) [3]'" << Key.kvalue << endl);
          Cerr("*****                    'origin (mm) [3]'"<< endl);
          Cerr("*****                    'offset [3]'" << Key.kvalue << endl);
          Cerr("*****                    Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }


      // --- Size thickness mm
      // "slice thickness (pixels)" / Intf_fields.slice_thickness_mm
      // We assuming the slice thickness is given in mm regardless of the name 
      else if(IntfCheckKeyMatch(Key, "slice thickness (pixels)"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->slice_thickness_mm) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'slice thickness (pixels)' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }

      // --- centre-centre slice separation (pixels). (used to compute slice spacing)
      // "process status" / Intf_fields.ctr_to_ctr_separation
      else if(IntfCheckKeyMatch(Key, "centre-centre slice separation (pixels)") ||
              IntfCheckKeyMatch(Key, "center-center slice separation (pixels)") )
      {
        if(ConvertFromString(Key.kvalue ,&ap_IF->ctr_to_ctr_separation) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'centre-centre slice separation (pixels)' key.");
          Cerr(" Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        Cerr("***** ReadIntfHeader()-> WARNING : 'centre-centre slice separation (pixels)' has no use in the current implementation !"<< endl);
      }
      
      // --- Number of time frames
      // "number of time frames" & "number of frame groups" / Intf_fields.nb_time_frames
      // Consistency with matrix size 4 : priority to nb_time_frames
      else if(IntfCheckKeyMatch(Key, "number of time frames") ||
              IntfCheckKeyMatch(Key, "number of frame groups")) 
      {
        if( ConvertFromString(Key.kvalue , &ap_IF->nb_time_frames) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of time frames' or 'number of frame groups' keys.");
          Cerr("Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      
      
      // --- Number of respiratory gates
      // "number of respiratory gates" / Intf_fields.nb_resp_gates
      else if(IntfCheckKeyMatch(Key, "number of respiratory gates")) 
      {
        if( ConvertFromString(Key.kvalue , &ap_IF->nb_resp_gates) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of respiratory gates' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      
      
      // --- Number of cardiac gates
      // "number of cardiac gates" / Intf_fields.nb_card_gates
      else if(IntfCheckKeyMatch(Key, "number of cardiac gates")) 
      {
        if( ConvertFromString(Key.kvalue , &ap_IF->nb_card_gates) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of cardiac gates' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        } 
      }
      
      // --- Total number of images
      // "total number of images" / Intf_fields.nb_total_imgs
      else if(IntfCheckKeyMatch(Key, "total number of images") ) 
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->nb_total_imgs) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'total number of images' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      
      // --- Number of bytes for each pixel/voxel of the image data
      // "number of bytes per pixel" / Intf_fields.nb_bytes_pixel
      else if(IntfCheckKeyMatch(Key, "number of bytes per pixel"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->nb_bytes_pixel) ) 
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of bytes per pixel' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }

      // --- Slice orientations
      // "slice orientation" / Intf_fields.slice_orientation
      // TODO : This information is recovered, but not used for the current implementation
      else if(IntfCheckKeyMatch(Key, "slice orientation"))
      {
        string slice_orientation;
        if( ConvertFromString(Key.kvalue , &slice_orientation) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'slice orientation' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        
        if (slice_orientation == "transverse")  
          ap_IF->slice_orientation = INTF_TRANSVERSE;
        else if (slice_orientation == "sagittal")
          ap_IF->slice_orientation = INTF_SAGITTAL;
        else if (slice_orientation == "coronal")
          ap_IF->slice_orientation = INTF_CORONAL;
      }

      // --- Patient rotation
      // "patient rotation" / Intf_fields.pat_rotation
      // TODO : This information is recovered, but not used for the current implementation
      else if(IntfCheckKeyMatch(Key, "patient rotation"))
      {
        string pat_rotation;
        if( ConvertFromString(Key.kvalue , &pat_rotation) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'patient rotation' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        
        if (pat_rotation == "supine")  
          ap_IF->pat_rotation = INTF_SUPINE;
        else if (pat_rotation == "prone")    
          ap_IF->pat_rotation = INTF_PRONE;
      }
      
      
      // --- Patient orientation
      // "patient orientation" / Intf_fields.pat_orientation
      // TODO : This information is recovered, but not used for the current implementation
      else if(IntfCheckKeyMatch(Key, "patient orientation"))
      {
        string pat_orientation;
        if( ConvertFromString(Key.kvalue , &pat_orientation) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'patient orientation' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        
        if (pat_orientation == "head_in")  
          ap_IF->pat_orientation = INTF_HEADIN;
        else if (pat_orientation == "feet_in")    
          ap_IF->pat_orientation = INTF_FEETIN;
      }


      // --- Multiplicative calibration value
      // "rescale slope" / Intf_fields.rescale_slope
      else if(IntfCheckKeyMatch(Key, "rescale slope")||
              IntfCheckKeyMatch(Key, "data rescale slope"))
      { 
        double rs= 1.;
        if( ConvertFromString(Key.kvalue , &rs) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'rescale slope' or 'data rescale slope' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        ap_IF->rescale_slope = rs; // factorization as this variable can also be modified by quantification units
        
        if (ap_IF->rescale_slope == 0.) // Only checking initialization, no need for FLTNBisEqual()
        {
          Cerr("***** ReadIntfHeader()-> Error : field 'resclale slope' units should be >0!" << endl);
          return 1;
        }
      }

      // Several global scale factors.
      // --- Multiplicative calibration value
      // "quantification units" / Intf_fields.rescale_slope
      // Note : throw warning if bad conversion, as this key is sometimes 
      // used with a string gathering the data unit (eg : kbq/cc)
      else if(IntfCheckKeyMatch(Key, "quantification units"))
      { 
        double qu= 1.;
        if(ConvertFromString(Key.kvalue , &qu) )
        {
          Cerr("***** ReadIntfHeader()-> WARNING : Error when trying to read numeric value from 'quantification units' key. Actual value = " << Key.kvalue << endl);
          Cerr("*****                              This key will be ignored" << endl);
          qu= 1.;
        }
        ap_IF->quant_units = qu; // factorization as this variable can also be modified by rescale slope
        
        if (ap_IF->quant_units == 0.) // Only checking initialization, no need for FLTNBisEqual()
        {
          Cerr("***** ReadIntfHeader()-> Error : field 'quantification units' should be >0!" << endl);
          return 1;
        } 
      }
      
      // --- Additive calibration value
      // "rescale intercept" / Intf_fields.rescale_intercept
      else if(IntfCheckKeyMatch(Key, "rescale intercept") ||
              IntfCheckKeyMatch(Key, "data rescale intercept"))
      { 
        if( ConvertFromString(Key.kvalue , &ap_IF->rescale_intercept) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'rescale intercept' or 'data rescale intercept' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        
        if (ap_IF->rescale_intercept == 0.) // Only checking initialization, no need for FLTNBisEqual()
        {
          Cerr("***** ReadIntfHeader()-> Error : field 'resclale intercept' units should be >0!" << endl);
          return 1;
        }
      }
      
      
      // --- Type of projeted/reconstructed dataset (Static|Dynamic|Gated|Tomographic|Curve|ROI|GSPECT|Other)
      //     Curve|ROI|Other are considered as Static
      // "type of data" / Intf_fields.data_type
      else if(IntfCheckKeyMatch(Key, "type of data"))
      {
        string data_str;
        if(ConvertFromString(Key.kvalue ,&data_str) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'type of data' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        
        ap_IF->data_type = IntfKeyGetInputImgDataType(data_str); 
      }
      
      
      // --- Acquisition duration
      // "study duration (sec)" / Intf_fields.study_duration
      else if(IntfCheckKeyMatch(Key, "study duration (sec)"))
      {
        if(ConvertFromString(Key.kvalue ,&ap_IF->study_duration) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'study duration (sec)' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      

      // --- Image duration
      // "image_duration (sec)" / Intf_fields.image_duration
      // Note : check after reading all headers that the number of image durations matches the actual number of frames
      else if(IntfCheckKeyMatch(Key, "image duration (sec)"))
      {
        FLTNB duration;
        if(ConvertFromString(Key.kvalue ,&duration) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'study duration (sec)' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        ap_IF->image_duration.push_back(duration);
      }



      // --- Image start time
      // "image start time (sec)" / Intf_fields.image_start_time
      // Note : check after reading all headers that the number of image start time matches the actual number of frames
      else if(IntfCheckKeyMatch(Key, "image start time (sec)"))
      {
        FLTNB stime;
        if(ConvertFromString(Key.kvalue ,&stime) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'study duration (sec)' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        ap_IF->image_start_time.push_back(stime);
      }
      
      
      
      // --- Pause between frame groups
      // "pause between frame groups (sec)" / Intf_fields.pause_duration
      // Note : check after reading all headers that the number of pauses matches the actual number of frames
      else if(IntfCheckKeyMatch(Key, "pause between frame groups (sec)"))
      {
        FLTNB pause_duration;
        if(ConvertFromString(Key.kvalue ,&pause_duration) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'study duration (sec)' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        ap_IF->frame_group_pause.push_back(pause_duration);
      }
      
      
      // TODO : this key is not used in the current implementation
      // --- number of time windows
      // "number of time windows " / Intf_fields.nb_time_windows
      else if(IntfCheckKeyMatch(Key, "number of time windows "))
      {
        if(ConvertFromString(Key.kvalue ,&ap_IF->nb_time_windows) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of time windows ' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        Cerr("***** ReadIntfHeader()-> WARNING : 'number of time windows' has no use in the current implementation !"<< endl);
      }
      
      
      // --- Process status : acquired or reconstructed
      // "process status" / Intf_fields.process_status
      else if(IntfCheckKeyMatch(Key, "process status"))
      {
        if(ConvertFromString(Key.kvalue ,&ap_IF->process_status) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'process status' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      }
      

      // --- Number of energy windows
      // "number of energy windows" / Intf_fields.nb_energy_windows
      // TODO : not implemented yet.
      else if(IntfCheckKeyMatch(Key, "number of energy windows"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->nb_energy_windows) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of energy windows' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
        //Cerr("***** ReadIntfHeader()-> WARNING : 'number of energy windows' has no use in the current implementation !"<< endl);
      }
      
      
      // --- Number of detector heads (SPECT systems)
      // "number of detector heads" / Intf_fields.nb_detector_heads
      else if(IntfCheckKeyMatch(Key, "number of detector heads"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->nb_detector_heads) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of detector heads' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      } 
      

      // --- Number of projections (for one head in SPECT)
      // "number of projections" / Intf_fields.nb_energy_windows
      else if(IntfCheckKeyMatch(Key, "number of projections"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->nb_projections) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of projections' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      } 
      
      
      // --- Angular span of the projections ex: 180, 360
      // "extent of rotation " / Intf_fields.extent_rotation
      else if(IntfCheckKeyMatch(Key, "extent of rotation"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->extent_rotation) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'extent of rotation' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      } 

      // --- Direction of rotation : CCW (counter-clockwise), CW (clockwise)
      // "direction of rotation" / Intf_fields.direction_rotation
      else if(IntfCheckKeyMatch(Key, "direction of rotation"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->direction_rotation) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'direction_rotation' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      } 
      
      // --- Angle of the first view
      // "start angle" / Intf_fields.first_angle
      // TODO : This may be wrong. First angle may be sum of values of :
      // "start angle" + "first projection angle in data set"
      else if(IntfCheckKeyMatch(Key, "start angle"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->first_angle) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'start angle' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      } 


      // --- All projection angles as an array key
      // "projection_angles" / Intf_fields.projection_angles
      else if(IntfCheckKeyMatch(Key, "projection_angles"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->projection_angles) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'projection_angles' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      } 
      

      // --- All distance between center of rotation and detectors, as unique value similar for each angles
      // "Center of rotation to detector distance" / Intf_fields.radius
      else if(IntfCheckKeyMatch(Key, "Center of rotation to detector distance"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->radius) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'Center of rotation to detector distance' key.");
          Cerr(" Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      } 
      
      
      // --- All distance between center of rotation and detectors, as an array key
      // "Radius" / Intf_fields.radius
      else if(IntfCheckKeyMatch(Key, "Radius"))
      {
        if(ConvertFromString(Key.kvalue , &ap_IF->radius) )
        {
          Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'Radius' key. Recovered value = " << Key.kvalue << endl);
          return 1;
        }
      } 
      
      // Return error if data compression
      else if(IntfCheckKeyMatch(Key, "data compression")) 
      {
        if(! (Key.kvalue == "none" || Key.kvalue.empty() ) )
        {
          Cerr("***** ReadIntfHeader()-> Error : compressed interfile images not handled by the current implementation !" << endl);
          return 1;
        }
      }
      
      // Return error if data encoded
      else if(IntfCheckKeyMatch(Key, "data encode")) 
      {
        if(! (Key.kvalue == "none" || Key.kvalue.empty() ) )
        {
          Cerr("***** ReadIntfHeader()-> Error : encoded interfile images not handled by the current implementation !" << endl);
          return 1;
        }
      }

  
      if(IntfCheckKeyMatch(Key, "number of slices") ) 
      {
        if(!(ap_IF->mtx_size[2] >0)) // Check if not been initialized elsewhere using other field (matrix size [3])
        {
          if(ConvertFromString(Key.kvalue , &ap_IF->mtx_size[2]) )
          {
            Cerr("***** ReadIntfHeader()-> Error when trying to read data from 'number of slices' key. Recovered value = " << Key.kvalue << endl);
            return 1;
          }
          ap_IF->nb_dims++;
        }
        continue;
      }
      
    }
  }
  else
  {
    Cerr("***** ReadIntfHeader()-> Error : couldn't find or read header interfile '"<< a_pathToHeaderFile << "' !" << endl);
    return 1;
  }
  
  // Recover slice number in ap_IF->mtx_size[2] if provided by "total number of images" or "number of images" tags, and not "matrix size [3]" tag. 
  if(ap_IF->mtx_size[2]==0 && ap_IF->nb_total_imgs>0)
    ap_IF->mtx_size[2] = ap_IF->nb_total_imgs;
  
  
  // Compute slice thickness if not initialized
  ap_IF->vox_size[2] = (ap_IF->vox_size[2] < 0) ?
                      ((ap_IF->vox_size[0]+ap_IF->vox_size[1])/2.) * ap_IF->slice_thickness_mm:
                        ap_IF->vox_size[2];
  
  // Compute nb dimensions for the image
  ap_IF->nb_dims = 3;
  
  if(ap_IF->mtx_size[3] >1 ||
     ap_IF->nb_time_frames>1)
    ap_IF->nb_dims++;
  if( ap_IF->nb_resp_gates >1)
    ap_IF->nb_dims++;
  if( ap_IF->nb_card_gates >1)
    ap_IF->nb_dims++;

  // If there are any frames, check nb frame durations matches nb (frames) group pauses
  if(ap_IF->image_duration.size()    > 0 &&
     ap_IF->frame_group_pause.size() > 0 &&
     ap_IF->image_duration.size() != ap_IF->frame_group_pause.size())
  {
    Cerr("***** ReadIntfHeader()-> Error : nb of recovered frame duration ('"<< ap_IF->image_duration.size()
      << ") does not match the nb of recovered pauses between frames ('"<< ap_IF->frame_group_pause.size() << ") !" << endl);
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteHeaderMainData
  \param a_path : path to the interfile header to write
  \param ap_IntfF : Structure containing interfile keys of the image to write
  \param vb : verbosity
  \brief Write main infos of an Interfile header at the path provided in parameter,
         using the interfile keys provided by the field structure parameter.
  \todo Recover in 'patient name' the name of the datafile used to reconstruct the images [ bed ]?
  \return 0 if success, positive value otherwise.
*/
int IntfWriteHeaderMainData(const string& a_path, const Intf_fields& ap_IntfF, int vb)
{
  if (vb >= 4)
  {
    Cout("--------------------------------------------------------------- " << endl);
    Cout("IntfWriteHeaderMainData()-> Start writing header interfile " << a_path << endl);
    Cout("--------------------------------------------------------------- " << endl);
  }
  
  string path_to_header, path_to_image;
  path_to_header = a_path;

  path_to_header.append(".hdr");
  path_to_image = GetFileFromPath(a_path).append(".img");
  
  // Get file and check existence
  ofstream ofile(path_to_header.c_str(), ios::out);
  ofile << setprecision(std::numeric_limits<FLTNB>::digits10 +1);
            
  if(ofile)
  {
    sScannerManager* p_scanMgr = sScannerManager::GetInstance(); 
    
    ofile << "!INTERFILE := " << endl; 

    // PET/SPECT/CT according to the scanner type defined in scanner manager
    string imaging_modality;
    if (p_scanMgr->GetScannerType() == SCANNER_PET)
      imaging_modality = "PET";
    else if (p_scanMgr->GetScannerType() == SCANNER_SPECT_CONVERGENT || p_scanMgr->GetScannerType() == SCANNER_SPECT_PINHOLE)
      imaging_modality = "SPECT";
    else if(p_scanMgr->GetScannerType() == SCANNER_CT)
      imaging_modality = "CT";
    else
      imaging_modality = "UNKOWN";

    ofile << "!imaging modality := " << imaging_modality << endl;
    
    // 3D SPECT use 3.3
    if(ap_IntfF.data_type == INTF_IMG_SPECT)
      ofile << "!version of keys := " << "3.3" << endl;
    else
      ofile << "!version of keys := " << "CASToRv" << CASTOR_VERSION << endl;
    // CASToR version
    ofile << "CASToR version := " << CASTOR_VERSION << endl;

    ofile << endl << "!GENERAL DATA := " << endl; 
    ofile << "!originating system := " << p_scanMgr->GetScannerName() << endl;

    // If more than one bed, then say it
    if (ap_IntfF.nb_bed_positions>1) ofile << "number of bed positions := " << ap_IntfF.nb_bed_positions << endl;
    // Otherwise, write the relative position
    else if (ap_IntfF.bed_position_provided) ofile << "horizontal bed relative position (mm) := " << ap_IntfF.bed_relative_position << endl;

    ofile << "!data offset in bytes := " << 0 << endl;
    ofile << "!name of data file := " << path_to_image << endl;

    //string patient_name = IntfKeyGetPatientNameTag();
    ofile << "patient name := " << GetFileFromPath(a_path) << endl;

    
    ofile << endl << "!GENERAL IMAGE DATA " << endl; 

    string data_type_str = (ap_IntfF.data_type == INTF_IMG_SPECT) ? 
                           "Tomographic" :
                           "Static";
    //ofile << "!type of data := " << data_type_str << endl; // such as Static|Dynamic|Gated|Tomographic|Curve|ROI|GSPECT|Other
    ofile << "!type of data := Dynamic" << endl; // /!\ Dynamic by default in order to be readable from ImageJ OpenMed plugin
     
    uint32_t nb_imgs = (ap_IntfF.process_status == "acquired") ?
                        ap_IntfF.nb_projections : // projeted data
                        ap_IntfF.nb_total_imgs  ; // reconstructed data
    ofile << "!total number of images := " << nb_imgs << endl;
    
    ofile << "imagedata byte order := " << IntfKeyGetEndianStr(ap_IntfF.endianness) << endl;

    // Depending on the mode of output writing, we should write a header for each dynamic image, or append everything in one single header
    // output dynamic files as one single file instead of 3D images with a meta header
    if( ap_IntfF.nb_time_frames>1 )
    {
      ofile << "number of time frames := " << ap_IntfF.nb_time_frames << endl;
      ofile << "!number of frame groups := " << ap_IntfF.nb_time_frames << endl;
    }
    else
      ofile << "!number of frame groups :=1 " << endl;
    
    
    if( ap_IntfF.nb_energy_windows>1 )
    ofile << "number of time windows := " << ap_IntfF.nb_resp_gates *
                                                   ap_IntfF.nb_card_gates << endl;
    if(ap_IntfF.nb_resp_gates > 1)
      ofile << "number of respiratory gates := " << ap_IntfF.nb_resp_gates << endl;
    if(ap_IntfF.nb_card_gates > 1)
      ofile << "number of cardiac gates := " << ap_IntfF.nb_card_gates << endl;

    if(ap_IntfF.process_status == "acquired")
    {
      ofile << "!number of projections := " << ap_IntfF.nb_projections << endl;
      ofile << "!extent of rotation := " << ap_IntfF.extent_rotation << endl;
    }
    ofile << "process status := " << ap_IntfF.process_status << endl;

      
    if (ap_IntfF.data_type == INTF_IMG_DYNAMIC)        // Dynamic study
      ofile <<  endl << "!DYNAMIC STUDY (General) :=" << endl;
    else if(p_scanMgr->GetScannerType() == 2)         // SPECT study
    {
      if( ap_IntfF.data_type == INTF_IMG_GATED)
        ofile << endl << "!GSPECT STUDY (General) :=" << endl;
      else
      {
        ofile << endl << "!SPECT STUDY (General) :=" << endl;
        ofile << endl << "!SPECT STUDY ("<< ap_IntfF.process_status <<" data) :="<< endl;
      }
    }
    else if (ap_IntfF.data_type == INTF_IMG_GATED )
      ofile << endl << "!GATED STUDY (General) :=" << endl;
    else                                              // Standard study
      ofile << endl << "!STATIC STUDY (General) :=" << endl;


    // Patient position
    // (not implemented as these informations are most of the time missing, and could be misleading)
    // output_header << "patient orientation := " << ap_IntfF.pat_orientation << endl;
    // output_header << "patient rotation := " << ap_IntfF.pat_rotation << endl;
    // output_header << "slice orientation := " << ap_IntfF.slice_orientation << endl;

 
    // Loop on dynamic images (if any) and write data
    
    if (ap_IntfF.nb_resp_gates>1)
      ofile << "!number of frame groups :=" << ap_IntfF.nb_time_frames << endl;
    
    // loop on frames
    for(int fr=0 ; fr<ap_IntfF.nb_time_frames ; fr++)
    {
      if( ap_IntfF.nb_time_frames>1 )
      {
        ofile << "!Dynamic Study (each frame group) :=" << endl;
        ofile << "!frame group number := " << fr+1 << endl;
      }
      
      // loop on respiratory gates
      for(int rg=0 ; rg<ap_IntfF.nb_resp_gates ; rg++)
      {
        if( ap_IntfF.nb_resp_gates>1 )
        {
          ofile << "!Respiratory Gated Study (each time window) :=" << endl;
          ofile << "!time window number := " << rg+1 << endl;
        }
        
          // loop on cardiac gates
          for(int cg=0 ; cg<ap_IntfF.nb_card_gates ; cg++)
          {
            if( ap_IntfF.nb_card_gates>1 )
            {
              ofile << "!Cardiac Gated Study (each time window) :=" << endl;
              ofile << "!time window number := " << cg+1 << endl;
            }

              // Write image specific data
              if(IntfWriteHeaderImgData(ofile, ap_IntfF, vb) )
              {
                Cerr("***** IntfWriteHeaderMainData()-> Error : while trying to write the interfile header '"<< path_to_header << "' !" << endl);
                return 1;
              }

            if( ap_IntfF.nb_card_gates>1 )
              ofile << "!number of images in this time window :=" 
                          << ap_IntfF.mtx_size[2] << endl;
             
          }
        
        if( ap_IntfF.nb_resp_gates>1 )
          ofile << "!number of images in this time window :="
                      << ap_IntfF.nb_card_gates *
                         ap_IntfF.mtx_size[2] << endl;
         
      }
      
      if (ap_IntfF.nb_time_frames>1)
      {
        ofile << "!number of images this frame group :=" 
                    << ap_IntfF.nb_resp_gates *
                       ap_IntfF.nb_card_gates *
                       ap_IntfF.mtx_size[2] << endl;
                                                                  
        ofile << "!image duration (sec) := " << ap_IntfF.image_duration[fr]  << endl;
        if(fr+1 == ap_IntfF.nb_time_frames) // Check if we reached the last frame
          ofile << "pause between frame groups (sec) " << 0.0 << endl;
        else
          ofile << "pause between frame groups (sec) " << ap_IntfF.frame_group_pause[fr] << endl;
      }
    }

    ofile << "!END OF INTERFILE := " << endl; 
  }
  else
  {
    Cerr("***** IntfWriteHeaderMainData()-> Error : couldn't find output header interfile '"<< a_path << "' !" << endl);
    return 1;
  }

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteHeaderImgData
  \param ap_ofile : pointer to the ofstream linked to the image header file 
  \param ap_IntfF : reference to a structure containing interfile keys of the image to write
  \param verbosity : verbosity
  \brief Write basic image data info in the file provided in parameter
  \todo Data about positionning/offsets ?
  \todo : include dectector heads (a set of fields for each head)
  \return 0 if success, positive value otherwise.
*/
int IntfWriteHeaderImgData(ofstream &ap_ofile, const Intf_fields& ap_IntfF, int vb)
{
  if(ap_ofile)
  {
    if(ap_IntfF.process_status == "acquired") // projection data (SPECT)
    {
      // Todo : include dectector heads (a set of fields for each head)
      // That means a lot of modifications as angles, radius, ... should be provided for each head
      ap_ofile << "!direction of rotation := " << ap_IntfF.direction_rotation << endl;
      ap_ofile << "start angle := " << ap_IntfF.first_angle << endl;
      ap_ofile << "projection angles := " << ap_IntfF.projection_angles << endl;

      // If one common radius for each projection (no "{}"), write it in the 'Radius' key
      if( ap_IntfF.radius.find("{}") != string::npos )
        ap_ofile << "Center of rotation to detector distance := " << ap_IntfF.radius << endl;
      else
        ap_ofile << "Radius := " << ap_IntfF.radius << endl;
        
      ap_ofile << "!matrix size [1] := " << ap_IntfF.mtx_size[0] << endl;
      ap_ofile << "!matrix size [2] := " << ap_IntfF.mtx_size[1] << endl;
      ap_ofile << "!number format := " << ap_IntfF.nb_format << endl;
      ap_ofile << "!number of bytes per pixel := " << sizeof(FLTNB) << endl;
      ap_ofile << "scaling factor (mm/pixel) [1] := " << ap_IntfF.vox_size[0] << endl;
      ap_ofile << "scaling factor (mm/pixel) [2] := " << ap_IntfF.vox_size[1] << endl;
      ap_ofile << "!data offset in bytes := " << 0 << endl;
    }
    else
    {
      ap_ofile << "number of dimensions := " << 3 << endl;
      ap_ofile << "!matrix size [1] := " << ap_IntfF.mtx_size[0] << endl;
      ap_ofile << "!matrix size [2] := " << ap_IntfF.mtx_size[1] << endl;
      ap_ofile << "!matrix size [3] := " << ap_IntfF.mtx_size[2] << endl;
      ap_ofile << "!number format := " << ap_IntfF.nb_format << endl;
      ap_ofile << "!number of bytes per pixel := " << sizeof(FLTNB) << endl;
      ap_ofile << "scaling factor (mm/pixel) [1] := " << ap_IntfF.vox_size[0] << endl;
      ap_ofile << "scaling factor (mm/pixel) [2] := " << ap_IntfF.vox_size[1] << endl;
      ap_ofile << "scaling factor (mm/pixel) [3] := " << ap_IntfF.vox_size[2] << endl;
      ap_ofile << "first pixel offset (mm) [1] := " << ap_IntfF.vox_offset[0] << endl;
      ap_ofile << "first pixel offset (mm) [2] := " << ap_IntfF.vox_offset[1] << endl;
      ap_ofile << "first pixel offset (mm) [3] := " << ap_IntfF.vox_offset[2] << endl;
      
      if(ap_IntfF.rescale_intercept != 1. ||
         ap_IntfF.rescale_slope != 0.     ||
         ap_IntfF.quant_units != 1.       )
      {
        ap_ofile << "data rescale offset := "  << ap_IntfF.rescale_intercept << endl;
        ap_ofile << "data rescale slope := " << ap_IntfF.rescale_slope << endl; 
        ap_ofile << "quantification units := " << ap_IntfF.quant_units << endl; 
      }
      // Todo, we may need keys for image positionning in space 
      //ap_ofile << "origin (mm) [1] := "  << ap_ID->Get??? << endl;
      //ap_ofile << "origin (mm) [2] := "  << ap_ID->Get??? << endl;
      //ap_ofile << "origin (mm) [3] := "  << ap_ID->Get??? << endl;
      
      
      // The following keys were used in interfiles format 3.3
      // Seem to have been replaced by matrix size [3]/scaling factor (mm/pixel) [3]
      // but we may need them for some compatibilities issues
      //  ap_ofile << "!number of slices := " << ap_IntfF.mtx_size[2] << endl;
      //  ap_ofile << "slice thickness (pixels) := " <<
      //    ap_IntfF.vox_size[2] /( (ap_IntfF.vox_size[0]+ap_IntfF.vox_size[1])/2.) << endl;
      //  ap_ofile << "centre-centre slice separation (pixels) := "  << endl;
    }
    
  }
  else
  {
    Cerr("***** IntfWriteHeaderImgData()-> Error : couldn't open provided interfile header file !" << endl);
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteImage()
  \param a_pathToImg : th to the directory where the image will be written
  \param ap_outImgMtx : Matrix containing the image to write
  \param a_dim : Number of voxels in the 3D image
  \param vb : Verbosity level
  \brief Write Interfile raw data whose path is provided in parameter, using image matrix provided in parameter.
  \brief For 1 dimensional image matrices
  \return 0 if success, positive value otherwise.
*/
int IntfWriteImage(const string& a_pathToImg, FLTNB* ap_outImgMtx, uint32_t a_dim, int vb)
{
  if(vb >= 5) Cout("IntfWriteImage()*" << endl);
  
  ofstream img_file(a_pathToImg.c_str(), ios::binary | ios::out);

  if(img_file)
  {
    // Call to writing function.
    if(IntfWriteData(&img_file, ap_outImgMtx, a_dim, vb) )
    {
      Cerr("***** IntfWriteImage()-> Error occurred when writing the image file '"<< a_pathToImg << "' !" << endl);
      return 1;
    }
  }
  else
  {
    Cerr("***** IntfWriteImage()-> Error occurred while trying to write the image file '"<< a_pathToImg << "' !" << endl);
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteImage
  \param ap_pathToImgs: List of string containing the paths to the image to write
  \param a2p_outImgMtx : 2 dimensional image matrix (1D temporal + 3D)
  \param ap_imgDim[2] : Dimensions of ap_outImgMtx
  \param vb : Verbosity level
  \brief Write Interfile raw data whose path is provided in parameter, using the image matrix provided in parameter.
  \brief For 2 dimensional image matrices
  \return 0 if success, positive value otherwise.
*/
int IntfWriteImage(vector<string> ap_pathToImgs, FLTNB** a2p_outImgMtx, uint32_t ap_dim[2], int vb)
{
  if(vb >= 5) Cout("IntfWriteImage()**" << endl);
  
  if(ap_pathToImgs.size() == 1) // We have one single data file containing all the data
  {
    // As the content will be appended, make sure there is no existing file with such name
    // remove it otherwise
    ifstream fcheck(ap_pathToImgs.at(0).append(".img").c_str());
    if(fcheck.good())
    {
      fcheck.close();
      #ifdef _WIN32
      string dos_instruction = "del " + ap_pathToImgs.at(0);
      system(dos_instruction.c_str());
      #else
      remove(ap_pathToImgs.at(0).c_str());
      #endif
    }
    
    // Add app to the options as we will write the image in several step.
    ofstream img_file(ap_pathToImgs.at(0).c_str(), ios::binary | ios::out | ios::app);
    
    if(img_file)
    {
      // Call to writing function for each image matrix
      for(uint32_t d1=0 ; d1<ap_dim[0] ; d1++)
        if( IntfWriteData(&img_file, a2p_outImgMtx[d1], ap_dim[1], vb) )
        {
          Cerr("***** IntfWriteImage()-> Error occurred when writing the image file '"<< ap_pathToImgs.at(d1) << "' !" << endl);
          return 1;
        }
    }
    else
    {
      Cerr("***** IntfWriteImage()-> Error occurred while trying to write the image file at the path:  '"<< ap_pathToImgs.at(0) << "' !" << endl);
      return 1;
    }
  }
  else // We have a number of data file equal to the number of image to load
  {
    if(ap_pathToImgs.size()!=ap_dim[0]) // Check we have a number of file corresponding to the number of images to load
    {
      Cerr("***** IntfWriteImage()-> Error : number of interfile images ("<< ap_pathToImgs.size() << ") not consistent with the number of images to load (" << ap_dim[0] << ") !" << endl);
      return 1;
    }
    
    for(uint32_t d1=0 ; d1<ap_dim[0] ; d1++)
    {
      ofstream img_file(ap_pathToImgs.at(d1).append(".img").c_str(), ios::binary | ios::out); // Should be one different path by file
      
      if(img_file)
      {
        // Call to writing function.
        if(IntfWriteData(&img_file, a2p_outImgMtx[d1], ap_dim[1], vb) )
        {
          Cerr("***** IntfWriteImage()-> Error occurred when writing the image file '"<< ap_pathToImgs.at(d1)  << "' !" << endl);
          return 1;
        }
      }
      else
      {
        Cerr("***** IntfWriteImage()-> Error occurred while trying to write the image file at the path:  '"<< ap_pathToImgs.at(d1) << "' !" << endl);
        return 1;
      }
    }
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteImage()
  \param vector<string> ap_pathToImgs: List of string containing the paths to the image to write
  \param FLTNB**** a4p_outImgMtx : 4 dimensional image matrix (3D temporal + 3D)
  \param int ap_imgDim[4] : Dimensions of ap_outImgMtx
  \param int vb : Verbosity level
  \brief Write Interfile raw data whose path is provided in parameter, using the image matrix provided in parameter.
  \brief For 4 dimensional image matrices
  \return 0 if success, positive value otherwise.
*/
int IntfWriteImage(vector<string> ap_pathToImgs, FLTNB**** a4p_outImgMtx, uint32_t ap_dim[4], int vb)
{
  if(vb >= 5) Cout("IntfWriteImage()" << endl);
  
  if(ap_pathToImgs.size() == 1) // We have one single data file containing all the data
  {
    // As the content will be appended, make sure there is no existing file with such name
    // remove it otherwise
    ifstream fcheck(ap_pathToImgs.at(0).append(".img").c_str());
    if(fcheck.good())
    {
      fcheck.close();
      #ifdef _WIN32
      string dos_instruction = "del " + ap_pathToImgs.at(0);
      system(dos_instruction.c_str());
      #else
      remove(ap_pathToImgs.at(0).c_str());
      #endif
    }

    // Add app to the options as we will write the image in several step.
    ofstream img_file(ap_pathToImgs.at(0).c_str(), ios::binary | ios::out | ios::app);

    if(img_file)
    {
      // Call to writing function for each image matrix
      for(uint32_t d1=0 ; d1<ap_dim[0] ; d1++)
        for(uint32_t d2=0 ; d2<ap_dim[1] ; d2++)
          for(uint32_t d3=0 ; d3<ap_dim[2] ; d3++)
          if(IntfWriteData(&img_file, a4p_outImgMtx[d1][d2][d3], ap_dim[3], vb) )
          {
            int idx_img = d1*ap_dim[1]*ap_dim[2] + d2*ap_dim[2] + d3;
            Cerr("***** IntfWriteImage()-> Error occurred when writing the image file '"<< ap_pathToImgs.at(idx_img) << "' !" << endl);
            return 1;
          }
    }
    else
    {
      Cerr("***** IntfWriteImage()-> Error occurred while trying to write the image file '"<< ap_pathToImgs.at(0) << "' !" << endl);
      return 1;
    }
    
  }
  else // We have a number of data file equal to the number of image to load
  {
    // Check we have a number of file corresponding to the number of images to load
    if(ap_pathToImgs.size() != ap_dim[0]*ap_dim[1]*ap_dim[2]) 
    {
      Cerr("***** IntfWriteImage()-> Error : number of interfile images (="<< ap_pathToImgs.size() 
                                       << ") not consistent with the number of images to load (=" 
                                       << ap_dim[0]*ap_dim[1]*ap_dim[2] << ") !" << endl);
      return 1;
    }

    for(uint32_t d1=0 ; d1<ap_dim[0] ; d1++)
      for(uint32_t d2=0 ; d2<ap_dim[1] ; d2++)
        for(uint32_t d3=0 ; d3<ap_dim[2] ; d3++)
        {
          int idx_img = d1*ap_dim[1]*ap_dim[2] + d2*ap_dim[2] + d3;
          ofstream img_file(ap_pathToImgs.at(idx_img).append(".img").c_str(), ios::binary | ios::out); // Should be one different path by file
          
          if(img_file)
          {
            // Call to writing function.
            if(IntfWriteData(&img_file, a4p_outImgMtx[d1][d2][d3], ap_dim[3], vb) )
            {
              Cerr("***** IntfWriteImage()-> Error occurred when writing the image file '"<< ap_pathToImgs.at(idx_img)  << "' !" << endl);
              return 1;
            }
          }
          else
          {
            Cerr("***** IntfWriteImage()-> Error occurred while trying to write the image file at the path:  '"<< ap_pathToImgs.at(idx_img) << "' !" << endl);
            return 1;
          }
        }
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfWriteData
  \param ap_oFile : Ofstream pointing to an image file
  \param ap_outImgMtx : 3D image matrix with reconstruction dimensions/voxel size containing the image data
  \param a_nbVox : A number of voxels in the 3D image matrix with reconstruction dimensions/voxel size
  \param vb : Verbosity level
  \brief Write the content of the image matrix in the file pointed by ofstream
  \todo  keep original orientations ? Would require a loop on voxels and a reindexing of voxel indexes
  \todo  check writing error
  \return 0 if success, positive value otherwise.
*/
int IntfWriteData(ofstream* ap_oFile, FLTNB* ap_outImgMatrix, int a_nbVox, int vb)
{
  if(vb >= 5) Cout("IntfWriteData() " << endl);

  // TODO : Keep original orientations ? Would require a loop on voxels and a reindexing of voxel indexes
  // TODO : Perhaps check here if nothing wrong happened during the writing (deletion of objects or writing to full disk). Not sure about how this is done for ofstream
  if (ap_oFile->write(reinterpret_cast<char*>(ap_outImgMatrix), a_nbVox*sizeof(FLTNB)) ) // implement error check here
  {
    //Cerr("***** IntfWriteData()-> Error occurred when trying to write the image file !" << endl);
    //return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyInitFields
  \param ap_IF : Structure containing Interfile keys
  \brief Init the keys of an the Interfile structure passed in parameter to default values
*/
void IntfKeyInitFields(Intf_fields* ap_IF)
{
  ap_IF->nb_bed_positions = 0;
  // The default bed relative position is intentionaly set to FLT_MIN to guess if the corresponding field has been supplied when reading an interfile image.
  // In any unexpected case, at least the FLT_MIN value is very close to 0. (e-38), so as it is interpreted in mm, it should not deviate the results in a significative extent.
  ap_IF->bed_relative_position = FLT_MIN;
  ap_IF->bed_position_provided = false;
  ap_IF->originating_system = "";
  ap_IF->path_to_image = "";
  ap_IF->endianness = User_Endianness; 
  ap_IF->data_offset = 0;
  ap_IF->nb_format = IntfKeyGetPixTypeStr();
  ap_IF->nb_dims = 0;
  for(int d=0 ; d<7 ; d++)
    ap_IF->mtx_size[d] = 0;
  for(int d=0 ; d<3 ; d++)
    ap_IF->vox_offset[d] = 0.;
  for(int d=0 ; d<3 ; d++)
    ap_IF->vox_size[d] = -1.;
    
  ap_IF->slice_thickness_mm = -1.;
  ap_IF->ctr_to_ctr_separation = 0;
  ap_IF->nb_time_frames = 1;
  ap_IF->nb_resp_gates = 1;
  ap_IF->nb_card_gates = 1;
  ap_IF->nb_total_imgs = 1;
  ap_IF->nb_bytes_pixel = 0;
  ap_IF->slice_orientation = 0;
  ap_IF->pat_rotation = 0;
  ap_IF->pat_orientation = 0;
  ap_IF->rescale_slope = 1.;
  ap_IF->rescale_intercept = 0.;
  ap_IF->quant_units = 1.;
             
  for(int d=0 ; d<3 ; d++)
  {
    ap_IF->cmtx_size[d] = 0;
    ap_IF->cvox_size[d] = -1.;
    ap_IF->cvox_offset[d] = 0.;
  }
  ap_IF->is_mtx_size_different = false;
  ap_IF->data_type = -1;
  ap_IF->study_duration = 1.;
  // ap_IF->image_duration = vector<float>;
  // ap_IF->image_start_time = vector<float>;
  // ap_IF->frame_group_pause = vector<float>;
  ap_IF->nb_time_windows = 0;
  ap_IF->process_status = "reconstructed";

  ap_IF->nb_detector_heads = 0;
  ap_IF->nb_energy_windows = 0;
  ap_IF->nb_projections = 1;
  ap_IF->extent_rotation = -1.;
  ap_IF->direction_rotation = "CW";
  ap_IF->first_angle = -1.;
  ap_IF->projection_angles = "";
  ap_IF->radius = "";  
};




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeySetFieldsOutput
  \param ap_IF : Structure containing Interfile keys
  \param ap_ID : oImageDimensionsAndQuantification object containing additional infos about reconstruction
  \brief Init the keys of the Interfile header of an image to be written on disk
  \details Init the keys of the Interfile structure passed in parameter for output writing \n
           using the ImageDimensions object containing information about the reconstruction
*/
void IntfKeySetFieldsOutput(Intf_fields* ap_IF, oImageDimensionsAndQuantification* ap_ID)
{
  // Set header metadata using Image Dimensions object
  ap_IF->endianness = User_Endianness;
  ap_IF->nb_format = IntfKeyGetPixTypeStr(); 
  ap_IF->nb_dims = 3;
  if(ap_ID->GetNbTimeFrames()>1) ap_IF->nb_dims++;
  if(ap_ID->GetNbRespGates()>1) ap_IF->nb_dims++;
  if(ap_ID->GetNbCardGates()>1) ap_IF->nb_dims++;
  ap_IF->mtx_size[0] = ap_ID->GetNbVoxX();
  ap_IF->mtx_size[1] = ap_ID->GetNbVoxY();
  ap_IF->mtx_size[2] = ap_ID->GetNbVoxZ();
  ap_IF->vox_size[0] = ap_ID->GetVoxSizeX();
  ap_IF->vox_size[1] = ap_ID->GetVoxSizeY();
  ap_IF->vox_size[2] = ap_ID->GetVoxSizeZ();
  ap_IF->vox_offset[0] = ap_ID->GetOffsetX();
  ap_IF->vox_offset[1] = ap_ID->GetOffsetY();
  ap_IF->vox_offset[2] = ap_ID->GetOffsetZ();
    
  ap_IF->slice_thickness_mm = ap_ID->GetVoxSizeZ();
  ap_IF->ctr_to_ctr_separation = 0;
  ap_IF->nb_time_frames = ap_ID->GetNbTimeFrames();
  ap_IF->nb_resp_gates = ap_ID->GetNbRespGates();
  ap_IF->nb_card_gates = ap_ID->GetNbCardGates();
  ap_IF->nb_total_imgs = ap_IF->nb_time_frames *
                         ap_IF->nb_resp_gates *
                         ap_IF->nb_card_gates*
                         ap_IF->mtx_size[2];
                         
  ap_IF->nb_bytes_pixel = sizeof(FLTNB);
  //ap_IF->slice_orientation = 0; // slice orientation : (=0, default)
  //ap_IF->pat_rotation = 0;      // patient rotation : supine (=0, default)
  //ap_IF->pat_orientation = 0;   // slice orientation : head-in (=0, default)
  ap_IF->rescale_slope=1.;     // multiplicative calibration values.
  ap_IF->rescale_intercept=0.; // additive calibration values.
  ap_IF->quant_units=1.;     // multiplicative calibration values.
  
  // Just required for reading, not for writing
  //ap_IF->cmtx_size[0] = ap_ID->GetNbVoxX();
  //ap_IF->cmtx_size[1] = ap_ID->GetNbVoxY();
  //ap_IF->cmtx_size[2] = ap_ID->GetNbVoxZ();
  //ap_IF->cvox_size[0] = ap_ID->GetVoxSizeX();
  //ap_IF->cvox_size[1] = ap_ID->GetVoxSizeY();
  //ap_IF->cvox_size[2] = ap_ID->GetVoxSizeZ();
  //ap_IF->cvox_offset[0] = ap_ID->GetOffsetX();
  //ap_IF->cvox_offset[1] = ap_ID->GetOffsetY();
  //ap_IF->cvox_offset[2] = ap_ID->GetOffsetZ();
  
  
  //ap_IF->is_mtx_size_different = false;
  ap_IF->data_type = IntfKeyGetOutputImgDataType(ap_ID);
  ap_IF->study_duration = ap_ID->GetFinalTimeStopInSec(0) -
                          ap_ID->GetFrameTimeStartInSec(0,0);
  for(int fr=0 ; fr<ap_ID->GetNbTimeFrames() ; fr++)
  {
    ap_IF->image_duration.push_back(ap_ID->GetFrameDurationInSec(0, fr));
    ap_IF->frame_group_pause.push_back((fr == 0) ? 0 
                                                 : ap_ID->GetFrameTimeStartInSec(0,fr) - ap_ID->GetFrameTimeStopInSec(0,fr-1));
  }
  ap_IF->nb_time_windows = ap_ID->GetNbRespGates() *
                           ap_ID->GetNbCardGates();
  //ap_IF->process_status;

  //ap_IF->detector_heads = 0;
  ap_IF->nb_energy_windows = ap_IF->nb_resp_gates *
                          ap_IF->nb_card_gates;
  //ap_IF->nb_projections
  //ap_IF->extent_rotation;
  //ap_IF->extent_rotation;
  //ap_IF->first_angle;
  //ap_IF->projection_angles;
  //ap_IF->radius;
  // Set the number of bed positions
  ap_IF->nb_bed_positions = ap_ID->GetNbBeds();
  // Set the bed position to 0 if multiple bed positions
  ap_IF->bed_relative_position = 0.;
  // Set it to the actual value if only one bed
  if (ap_ID->GetNbBeds()==1) ap_IF->bed_relative_position = ap_ID->GetBedPosition(0);
  // Set the flag that say if the bed relative position was provided from the datafile or not; in order
  // to know if this information should be written in the image header or not
  ap_IF->bed_position_provided = ap_ID->GetProvidedBedPositionFlag();
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyPrintFields
  \param ap_IF
  \brief Print all the keys of the Intf_fields structure passed in parameter, as well as their values for debugging purposes
*/
void IntfKeyPrintFields(Intf_fields a_IF)
{
  Cout("// ------ IntfKeyPrintFields ------ // " << endl << endl);
  Cout("path_to_image         : " << a_IF.path_to_image << endl);
  Cout("endianness            : " << IntfKeyGetEndianStr(a_IF.endianness) << endl);
  Cout("data_offset           : " << a_IF.data_offset << endl);
  Cout("nb_format             : " << a_IF.nb_format << endl);
  Cout("nb_dims               : " << unsigned(a_IF.nb_dims) << endl); // uint_8t doesn't output with cout, hence the (us) cast
  for(int i=0 ; i<7 ; i++)
    Cout("mtx_size["<<i<<"]           : " << a_IF.mtx_size[i] << endl);
  for(int i=0 ; i<3 ; i++)
    Cout("vox_offset["<<i<<"]         : " << a_IF.vox_offset[i] << endl);
  for(int i=0 ; i<3 ; i++)
    Cout("vox_size["<<i<<"]           : " << a_IF.vox_size[i] << endl);
    
  Cout("slice_thickness_mm    : " << a_IF.slice_thickness_mm << endl);
  Cout("ctr_to_ctr_separation : " << a_IF.ctr_to_ctr_separation << endl);
  Cout("nb_time_frames        : " << a_IF.nb_time_frames << endl);
  Cout("nb_resp_gates         : " << a_IF.nb_resp_gates << endl);
  Cout("nb_card_gates         : " << a_IF.nb_card_gates << endl);
  Cout("nb_total_imgs         : " << a_IF.nb_total_imgs << endl);
  Cout("nb_bytes_pixel        : " << unsigned(a_IF.nb_bytes_pixel) << endl); // uint_8t doesn't output with cout, hence the (us) cast
  Cout("slice_orientation     : " << a_IF.slice_orientation << endl);
  Cout("pat_rotation          : " << a_IF.pat_rotation << endl);
  Cout("pat_orientation       : " << a_IF.pat_orientation << endl);
  Cout("rescale_slope         : " << a_IF.rescale_slope << endl);
  Cout("rescale_intercept     : " << a_IF.rescale_intercept << endl);
  Cout("quantification_units  : " << a_IF.quant_units << endl);
  for(int i=0 ; i<3 ; i++)
    Cout("cmtx_size["<<i<<"]          : " << a_IF.cmtx_size[i] << endl);
  for(int i=0 ; i<3 ; i++)
    Cout("cvox_size["<<i<<"]          : " << a_IF.cvox_size[i] << endl);
  for(int i=0 ; i<3 ; i++)
    Cout("cvox_offset["<<i<<"]          : " << a_IF.cvox_offset[i] << endl);
    
  Cout("is_mtx_size_different : " << a_IF.is_mtx_size_different << endl);
  Cout("data_type             : " << a_IF.data_type << endl);
  Cout("study_duration        : " << a_IF.study_duration << endl);
  Cout("image_duration(fr)    : " << endl);
  for(uint32_t fr=0 ; fr<a_IF.image_duration.size() ; fr++)
    Cout("image_duration("<<fr+1<<")     : " << a_IF.image_duration.at(fr) << endl);
  Cout("image_start_time(fr)    : " << endl);
  for(uint32_t fr=0 ; fr<a_IF.image_start_time.size() ; fr++)
    Cout("image_start_time("<<fr+1<<")     : " << a_IF.image_start_time.at(fr) << endl);
  Cout("pause_duration(fr)    : " << endl);
  for(uint32_t fr=0 ; fr<a_IF.frame_group_pause.size() ; fr++)
    Cout("pause_duration("<<fr+1<<")     : " << a_IF.frame_group_pause.at(fr) << endl);
  
    //ap_IF->image_duration = vector<float>;
  Cout("nb_time_windows       : " << a_IF.nb_time_windows << endl);
  Cout("process_status        : " << a_IF.process_status << endl);

  // SPECT and projection related data
  Cout("nb_detector_heads     : " << a_IF.nb_detector_heads << endl);
  Cout("nb_energy_windows     : " << a_IF.nb_energy_windows << endl);
  Cout("nb_projections        : " << a_IF.nb_projections << endl);
  Cout("extent_rotation       : " << a_IF.extent_rotation << endl);
  Cout("direction_rotation    : " << a_IF.direction_rotation << endl);
  Cout("first_angle           : " << a_IF.first_angle << endl);
  Cout("projection_angles     : " << a_IF.projection_angles << endl);
  Cout("radius                : " << a_IF.radius << endl);
  Cout("// ------ ------------------ ------ // " << endl << endl);
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfRecoverKey
  \param ap_Key : Structure to recover the parsed key components (key, value,..)
  \param a_line : String to process
  \brief Process the line passed in parameter and write the key information in the ap_Key Intf_key member structure
  \details .korig : Get original line without comments
           .kcase : Get key without spaces and without comments
           .klcase: Same as kcase, in lower case
           .kvalue: Value of the key, without spaces
  \todo Check that IntfToLowerCase() doesn't cause issue with some characters or some ASCII file format (unicode, etc..)
  \return 0 if success, positive value otherwise.
*/
int IntfRecoverKey(Intf_key* ap_Key, const string& a_line)
{
  string intf_sep = ":=";
  
  // Remove any comment from the original key line
  int pos = a_line.find_first_of(';',0);
  ap_Key->korig = a_line.substr(0, pos);
  
  // check for interfile separator.
  pos = ap_Key->korig.find_first_of(intf_sep);
  
  // If interfile key is not found (not an interfile key or wrong), just add it at the end of the key and proceed
  if (ap_Key->korig.find(intf_sep) == string::npos)
    ap_Key->korig.append(":=");
  
  // Get key title
  ap_Key->kcase = ap_Key->korig.substr(0,pos);
  // Get key value
  ap_Key->kvalue = ap_Key->korig.substr(pos+2);

  // Get case insensitive key, and without space
  ap_Key->klcase = ap_Key->kcase;
  IntfToLowerCase(&ap_Key->klcase); // TODO Safe ?
  
  // Get key value
  ap_Key->kvalue = ap_Key->korig.substr(pos+2);

  // Clear the space before the first and after the last element in the key and value;
  IntfEraseSpaces(&ap_Key->klcase);
  IntfEraseSpaces(&ap_Key->kvalue);
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfCheckKeyMatch(Intf_key ap_Key, const string& a_field)
  \param ap_Key : Structure containing the parsed key components (key, value,..)
  \param a_line : String containing an interfile key
  \brief Check if the key matches the string passed in parameter
  \todo Be sure it is appropriate to use Intf_key.klcase for each key comparison
  \return 1 if success, 0 otherwise (not found).
*/
int IntfCheckKeyMatch(Intf_key ap_Key, const string& a_field)
{
/*
  // TODO Check if appropriate to use .klcase in any situation for key comparison
  // SS: but the klcase has already been applied by the IntfRecoverKey called before calling this IntfCheckKeyMatch function ?
  if(ap_Key.klcase == a_field)
    return 1;
  else
    return 0;
*/
  // ----------------
  // SS: we now remove all blanks, tabulations, carriage returns, end of line, and compare the strings without any of these characters
  // ----------------
  string a_copy_of_the_key = ap_Key.klcase;
  string a_copy_of_the_field = a_field;
  size_t found_char_at_pos = string::npos;
  // Process the key
  while ( (found_char_at_pos=a_copy_of_the_key.find_first_of(" \t\r\n")) != string::npos) a_copy_of_the_key.replace(found_char_at_pos,1,"");
  // Process the field
  while ( (found_char_at_pos=a_copy_of_the_field.find_first_of(" \t\r\n")) != string::npos) a_copy_of_the_field.replace(found_char_at_pos,1,"");
  // Now compare them
//cout << "Copy of the key: " << a_copy_of_the_key << endl;
//cout << "Copy of the field: " << a_copy_of_the_field << endl;
  if (a_copy_of_the_key==a_copy_of_the_field) return 1;
  else return 0;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyIsArray()
  \param ap_Key
  \brief Check if the key passed in parameter is an array (contains brackets '{' and '}' )
  \return 1 if success, 0 otherwise (not array).
*/
int IntfKeyIsArray(Intf_key ap_Key)
{
  if(ap_Key.kvalue.find("{") != string::npos && 
     ap_Key.kvalue.find("}") != string::npos)
    return 1;
    
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyGetArrayNbElts
  \param a_Key
  \brief Return the number of elts in an Interfile array Key
  \return the number of elements in the array key, or negative value if error
*/
int IntfKeyGetArrayNbElts(Intf_key a_Key)
{
  int nb_elts = 0;
  string val_str = a_Key.kvalue;
  
  size_t pos = val_str.find_first_of('{')+1;

  while (pos < val_str.length()) 
  {
    size_t pos_c = val_str.find_first_of(",", pos);

    // no comma found, looking for closing bracket
    if(pos_c == string::npos)
    {
      pos_c = val_str.find_first_of("}", pos);
      if(! (IntfKeyIsArray(a_Key)) )
      {
        Cerr("***** IntfKeyGetArrayNbElts-> Error : closing bracket not found in interfile array key : "<< a_Key.korig << " !" << endl);
        return -1;
      }
      else
        return nb_elts+1; // add last elt before the end bracket
    }
    pos = pos_c+1;
    nb_elts++;
  }

  // return error if we end of key if reached without closing bracket
  Cerr("***** IntfKeyGetArrayNbElts-> Error : closing bracket not found in interfile array key : "<< a_Key.korig << " !" << endl);
  return -1;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyGetMaxArrayKey
  \param ap_Key
  \brief Return the maximum value in an array (key value contains brackets '{,,}' )
  \return the max value in the array key.
*/
int IntfKeyGetMaxArrayKey(Intf_key ap_Key)
{
  int value, max=0;
  string val_str = ap_Key.kvalue;
  if (val_str == "") return(max);

  size_t pos = val_str.find_first_of('{')+1;
  
  while (pos < val_str.length()) 
  {
    size_t  pos_c = val_str.find_first_of(",", pos);
    
    // no comma found, looking for closing bracket
    if(pos_c == string::npos) pos_c = val_str.find_first_of("}", pos);
        
    if(ConvertFromString(val_str.substr(pos,pos_c-pos), &value) )
    {
      Cerr("***** IntfKeyGetMaxArrayKey()-> An error occured when trying to recover the following value from the array key : "<< val_str.substr(pos,pos_c-pos) << " !" << endl);
      return 1;
    }
    
    if (value > max) max = value;
    
    pos = pos_c+1;
  }

  return max;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn int IntfKeyGetArrayElts()
  \param ap_Key
  \param T* ap_return : Templated parameter in which the elts will be returned
  \brief Get all the elements in an array key in a templated array passed in parameter. 
         It assumes the return variable has been instanciated with a correct number of elements.
  \return 0 if success, positive value otherwise
*/
template<class T>
int IntfKeyGetArrayElts(Intf_key a_Key, T* ap_return)
{
  string val_str = a_Key.kvalue;
  
  // Check if interfile key
  if(! (IntfKeyIsArray(a_Key)) )
  {
    Cerr("***** IntfKeyGetArrayElts-> Error : Problem reading the following interfile array key : "<< a_Key.korig << " !" << endl);
    return 1;
  }
  
  if (val_str == "")
  {
    Cerr("***** IntfKeyGetArrayElts-> Error : no elements in the array key : "<< a_Key.korig << " !" << endl);
    return 1;
  }
  
  size_t pos = val_str.find_first_of('{')+1;

  int elt = 0;
  
  while (pos < val_str.length()) 
  {
    size_t pos_c = val_str.find_first_of(",", pos);

    // no comma found, looking for closing bracket
    if(pos_c == string::npos) pos_c =  val_str.find_first_of("}", pos);
        
    if(ConvertFromString(val_str.substr(pos,pos_c-pos), &ap_return[elt]) )
    {
      Cerr("***** IntfKeyGetMaxArrayKey()-> An error occured when trying to recover the following value from the array key : "<< val_str.substr(pos,pos_c-pos) << " !" << endl);
      return 1;
    }
    
    pos = pos_c+1;
    elt++;
  }

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfGetVoxIdxSHTOrientation(Intf_fields a_IF, int a_voxId)
  \param ap_IF
  \param a_voxId : index of the voxel in a 1-D image vector
  \brief Compute a voxel index corresponding to the default orientation (Sup/Hin/Trans) \n
         from the orientation informations contained in the Intf_fields passed in parameter
  \todo NOT CURRENTLY IMPLEMENTED
  \return new voxel index
*/  
int IntfGetVoxIdxSHTOrientation(Intf_fields a_IF, int a_voxId)
{
  int dimX=a_IF.vox_size[0], dimY=a_IF.vox_size[1], dimZ=a_IF.vox_size[2];
  int dimXY=dimX*dimY;
  
  // Get x,y,z original image indexes
  int z = a_voxId/dimXY;
  int y = (a_voxId - z*dimXY) / dimX;
  int x = a_voxId - z*dimXY - y*dimX;
  
  // new image indexes
  int X = x; 
  int Y = y;
  int Z = z; 
  
  // Convert to default orientations (TRANSVERSE / SUPINE / HEADIN)
  if(a_IF.slice_orientation == INTF_TRANSVERSE) // X=x, Y=y, Z=z
  {
    if(a_IF.pat_orientation == INTF_HEADIN)
    {
      if(a_IF.pat_rotation == INTF_SUPINE)
      {
        // nothing to change
        return a_voxId;
      }
      else // a_IF.pat_rotation == INTF_PRONE
      {
        X = dimX-x;
        Y = dimY-y;
      }
    }
    else // a_IF.pat_orientation == INTF_FEETIN
    {
      if(a_IF.pat_rotation == INTF_SUPINE)
      {
        X = dimX-x;
        Z = dimZ-z;
      }
      else // a_IF.pat_rotation == INTF_PRONE
      {
        Y = dimY-y;
        Z = dimZ-z;
      }
    }
  }
  else if(a_IF.slice_orientation == INTF_CORONAL) // X=x, Y=z, Z=y
  {
    if(a_IF.pat_orientation == INTF_HEADIN)
    {
      if(a_IF.pat_rotation == INTF_SUPINE)
      {
        Y = z;
        Z = y;
      }
      else // a_IF.pat_rotation == INTF_PRONE
      {
        X = dimX-x;
        Y = dimZ-z;
        Z = y;
      }
    }
    else // a_IF.pat_orientation == INTF_FEETIN
    {
      if(a_IF.pat_rotation == INTF_SUPINE)
      {
        X = dimX-x;
        Y = z;
        Z = dimY-y;
      }
      else // a_IF.pat_rotation == INTF_PRONE
      {
        X = dimX-x;
        Y = dimZ-z;
        Z = dimY-y;
      }
    }
  }
  else // a_IF.slice_orientation == INTF_SAGITTAL // X=z, Y=x, Z=y
  {
    if(a_IF.pat_orientation == INTF_HEADIN)
    {
      if(a_IF.pat_rotation == INTF_SUPINE)
      {
        X = dimZ-z;
        Y = dimX-x;
        Z = y;
      }
      else // a_IF.pat_rotation == INTF_PRONE
      {
        X = z;
        Y = x;
        Z = y;
      }
    }
    else // a_IF.pat_orientation == INTF_FEETIN
    {
      if(a_IF.pat_rotation == INTF_SUPINE)
      {
        X = z;
        Y = dimX-x;
        Z = dimY-y;
      }
      else // a_IF.pat_rotation == INTF_PRONE
      {
        X = dimZ-z;
        Y = x;
        Z = dimY-y;
      }
    }
  }
  
  a_voxId = X + Y*dimX + Z*dimX*dimY;
  return a_voxId;
}    




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyGetEndianStr()
  \param a_val : 
  \brief return the endian string corresponding to the value passed in parameter (see module INTF_ENDIANNESS).
  \return "BIG_ENDIAN" if 0, \n
          "LITTLE_ENDIAN" if 1, \n
          "UNKNOWN" otherwise
*/
string IntfKeyGetEndianStr(int a_val)
{
  if(a_val == 0) return "BIGENDIAN";
  if(a_val == 1) return "LITTLEENDIAN";
  return "UNKNOWN";
}
    

/*
  \fn IntfKeyGetModalityStr
  \param a_modalityIdx
  \brief Convert the integer provided in parameter to the string related
         to the corresponding modality as defined by the scanner objects
  \todo Add other modalities as we implement them
  \return string corresponding to the modality
*/
string IntfKeyGetModalityStr(int a_modalityIdx)
{
  if(a_modalityIdx == 0)
    return "PET";
  else if(a_modalityIdx == 1)
    return "SPECT";
  /*
  else if(a_modalityIdx == 2)
    return "CT";
  else if(a_modalityIdx == 3) // Pinhole
    return "SPECT";
  */
  else
    return "UNKNOWN";
}
    



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn int IntfKeyGetInputImgDataType()
  \param a_str : original string
  \brief Get the image data type corresponding to the the input string
  \return int value corresponding to the image data type (see module INTF_IMG_TYPE).
*/
int IntfKeyGetInputImgDataType(const string& a_str)
{
  if (a_str == "static")
    return INTF_IMG_STATIC;
  if (a_str == "dynamic")
    return INTF_IMG_DYNAMIC;
  if (a_str == "gated")
    return INTF_IMG_GATED;
  if (a_str == "tomographic")
    return INTF_IMG_SPECT;
  if (a_str == "gspect")
    return INTF_IMG_GSPECT;
  if (a_str == "pet")
    return INTF_IMG_PET;
    
  return INTF_IMG_UNKNOWN;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyGetOutImgDataType
  \param ap_ID
  \brief Get the image data type corresponding to the image metadata passed in parameter
  \return int value corresponding to the image data type (see module INTF_IMG_TYPE).
*/
int IntfKeyGetOutputImgDataType(oImageDimensionsAndQuantification* ap_ID)
{
  sScannerManager* p_scanMgr = sScannerManager::GetInstance(); 
  
  int data_type = 0;
  
  // Update data type key
  if(ap_ID->GetNbTimeFrames() > 1)
    data_type = INTF_IMG_DYNAMIC;
  else if(p_scanMgr->GetScannerType() == 2 &&
          (ap_ID->GetNbRespGates() > 1 ||
           ap_ID->GetNbCardGates() > 1 ) )
    data_type = INTF_IMG_GSPECT;
  else if(ap_ID->GetNbRespGates() > 1 ||
          ap_ID->GetNbCardGates() > 1 )
    data_type = INTF_IMG_GATED;
  else if(p_scanMgr->GetScannerType() == 2)
    data_type = INTF_IMG_SPECT;
  else
    data_type = INTF_IMG_STATIC;

  return data_type;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyGetPixTypeStr()
  \brief Return the string corresponding to the nb of bytes in the type FLTNB
  \return string corresponding to the pixel data type
*/
string IntfKeyGetPixTypeStr()
{
  // Return either "long long float" (long double type), "long float" (doubletype ) or "short float" (float type)
  if (sizeof(FLTNB) == 4) return "short float";
  else if (sizeof(FLTNB) == 8) return "long float";
  else if (sizeof(FLTNB) == 16) return "long long float";
  else
  {
    Cerr("***** oInterfileIO::IntfKeyGetPixTypeStr() -> Size of current float type (" << sizeof(FLTNB) << ") does not correspond to a known type !" << endl);
    Exit(EXIT_FAILURE);
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfKeyGetPatOrientation()
  \param ap_IF
  \brief Get the complete patient orientation from an Intf_fields structure according to 
         the values of keys 'slice orientation', 'patient rotation', and 'patient orientation'
  \return int value corresponding to the patient orientation (see module PATIENT_ORIENTATION).
*/
int IntfKeyGetPatOrientation(Intf_fields ap_IF)
{
  switch (ap_IF.pat_rotation) 
  {
    case INTF_SUPINE: 
      switch (ap_IF.pat_orientation) 
      {
        case INTF_HEADIN: 
          switch (ap_IF.slice_orientation) 
          {
            case INTF_TRANSVERSE:
              return(INTF_SUPINE_HEADIN_TRANSAXIAL); break;
            case INTF_SAGITTAL   :
              return(INTF_SUPINE_HEADIN_SAGITTAL);   break;
            case INTF_CORONAL    :
              return(INTF_SUPINE_HEADIN_CORONAL);    break; 
          }
          break;
        
        case INTF_FEETIN: 
          switch(ap_IF.slice_orientation) 
          {
            case INTF_TRANSVERSE: 
              return(INTF_SUPINE_FEETIN_TRANSAXIAL); break;
            case INTF_SAGITTAL   :
              return(INTF_SUPINE_FEETIN_SAGITTAL);   break;
            case INTF_CORONAL    :
              return(INTF_SUPINE_FEETIN_CORONAL);    break;
          }
          break;
      }
      break;
      
    case INTF_PRONE : 
      switch (ap_IF.pat_orientation) 
      {
        case INTF_HEADIN:  
          switch (ap_IF.slice_orientation) 
          {
            case INTF_TRANSVERSE: 
              return(INTF_PRONE_HEADIN_TRANSAXIAL);  break;
            case INTF_SAGITTAL   :
              return(INTF_PRONE_HEADIN_SAGITTAL);    break;
            case INTF_CORONAL    :
              return(INTF_PRONE_HEADIN_CORONAL);     break;
          }
          break;
        case INTF_FEETIN:
          switch (ap_IF.slice_orientation) 
          {
            case INTF_TRANSVERSE : 
              return(INTF_PRONE_FEETIN_TRANSAXIAL);  break;
            case INTF_SAGITTAL   :
              return(INTF_PRONE_FEETIN_SAGITTAL);    break;
            case INTF_CORONAL    :
              return(INTF_PRONE_FEETIN_CORONAL);     break;
          }
          break;
      }
      break;
  }

  return(INTF_SUPINE_HEADIN_TRANSAXIAL); // default
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetUserEndianness()
  \brief Check user/host computer endianness and write it to the global variable User_Endianness
  \details This function should be called once during the initialization step of the algorithm (currently, the singleton initialization)
  \todo Maybe better compute it in a preprocessor macro
*/
void GetUserEndianness()
{
  int n = 1;
  User_Endianness = (*(char *)&n) == 1 ? 1 : 0 ;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfEraseSpaces()
  \param string* input_str
  \brief Erase space, blank characters (\t\r\n), and '!' before and after the characters in the string passed in parameter
*/
void IntfEraseSpaces(string* input_str)
{
  input_str->erase(0, input_str->find_first_not_of(" !\t\r\n")); // Erase all blank stuff before first character
  input_str->erase(input_str->find_last_not_of(" \t\r\n")+1 , input_str->length()); // Erase all blank stuff after last character
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IntfToLowerCase()
  \param string* ap_str : original string
  \brief Set all characters of the string passed in parameter to lower case
  \todo May have issue with non ASCII characters, file decoding, etc..
*/
void IntfToLowerCase(string* ap_str)
{
  std::transform(ap_str->begin(), ap_str->end(), ap_str->begin(), ::tolower);
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn void IntfKeyGetPatientNameTag()
  \brief Recover datafile name(s) stored in sOutputManager in one string
*/
string IntfKeyGetPatientNameTag()
{
  // Recover here the name of the datafile used to reconstruct the images [ bed ]
  string patient_tag_str = "";
  sOutputManager* p_outputMgr = sOutputManager::GetInstance();

  vector<string> dfName = p_outputMgr->GetDataFileName();

  if(dfName.size() > 1)
  {
    patient_tag_str +=  "{ " + dfName[0];
    for(size_t n=1 ; n<dfName.size() ; n++ )
      patient_tag_str +=  ", " + dfName[n];
    patient_tag_str +=  "} ";
  }
  else if (dfName.size() == 1)
    patient_tag_str += dfName[0];
  else
    patient_tag_str += "unknown data";

  return patient_tag_str;
}
      
         
      
      
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SwapBytes()
  \param T *ap_type : Variable of type T
  \brief Use std::reverse to swap the bits of a variable of any type
  \details Used for Little/Big Endian conversion
*/
template <class T>
void SwapBytes(T *ap_type)
{
  char *buffer = reinterpret_cast<char*>(ap_type);
  std::reverse(buffer, buffer + sizeof(T));

  //int i, j;
  //for (i=0,j=bytes-1;i < (bytes/2); i++, j--) 
  //{
  //   ptr[i]^=ptr[j]; ptr[j]^=ptr[i]; ptr[i]^=ptr[j];
  //}
  
  //for (int i = 0; i < size; ++i) 
  //    buffer[i] = ((buffer[i]>> 8) | (buffer[i] << 8));
}
