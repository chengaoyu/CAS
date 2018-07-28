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
  \ingroup management
  \brief   This file is used for all kind of different functions designed for options parsing and ASCII file reading. \n
           All functions included here are not class members.
*/

#ifndef GOPTIONS_HH
#define GOPTIONS_HH 1

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/**
 * @defgroup READ_CHECK_KEYWORD Read check keyword
 *
 *    \brief Keywords to check whether reading in text file was successful or not \n
 *           Defined in gOptions.hh
 * @{
 */

/** Constant corresponding to a mandatory field, which has to be recovered from a file (=true) */
#define KEYWORD_MANDATORY true
/** Constant corresponding to an optional field, not necessary to be recovered from a file (=false) */
#define KEYWORD_OPTIONAL false
/** Constant corresponding to the case of a mandatory keyword not found (=1) */
#define KEYWORD_MANDATORY_NOT_FOUND 1
/** Constant corresponding to the case of an optional keyword read with success (=0) */
#define KEYWORD_OPTIONAL_SUCCESS 0
/** Constant corresponding to the case of an optional keyword rising an error (=1) */
#define KEYWORD_OPTIONAL_ERROR 1
/** Constant corresponding to the case of an optional keyword not found (=2) */
#define KEYWORD_OPTIONAL_NOT_FOUND 2
/** @} */

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
// Functions for command line parameters reading
// Functions parsing a string composed of several numbers separated by commas.
/*!
  \fn      int ReadStringOption(const string& a_input, T* ap_return, int a_nbElts, const string& sep, const string& a_option)
  \param   a_input : string to parse
  \param   ap_return : (templated) array in which the parsed elements will be returned
  \param   a_nbElts : a number of elements to parse
  \param   sep : the separator (usually comma)
  \param   a_option : string indicating from where the function has been called
  \brief   Parse the 'a_input' string corresponding to the 'a_option' into 'a_nbElts' elements, using the 'sep' separator. \n
           The results are returned in the templated 'ap_return' dynamic templated array. \n
           Call "ConvertFromString()" to perform the correct conversion depending on the type of the data to convert
  \return  0 if success, and positive value otherwise.
*/
template<typename T> int ReadStringOption(const string& a_input, T* ap_return, int a_nbElts, const string& sep, const string& a_option);

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
// Functions for ASCII file input reading
/*!
  \fn      int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, bool a_mandatoryFlag)
  \param   a_file : string containing the path to the ASCII file
  \param   a_keyword : key related to the data we want to recover
  \param   ap_return : templated array in which the data will be returned
                       the recovered data will be converted to its type
                       its size should be equal to the number of elts to recover
  \param   a_nbElts : number of elements to recover
  \param   a_mandatoryFlag : boolean indicating the data to recover is  \n
                             mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found) \n
                             or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \brief   Look for "a_nbElts" elts in the "a_file" file matching the "a_keyword" string passed as parameter  \n
           and return the corresponding value(s) in the "ap_return" templated array. \n
           This function expect the data to recover to be written on one line after the key
  \details This function assumes the following separators : \n
           ":" is used as separator between the keyworld and the value \n
           "#" is used for comment. Every following characters will be discarded \n
           "," is used to separate elements if more than one are required
  \return  0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<typename T> int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, bool a_mandatoryFlag);

/*!
  \fn      int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag)
  \param   a_file : string containing the path to the ASCII file
  \param   a_keyword : key related to the data we want to recover
  \param   ap_return : templated array in which the data will be returned \n
                       the recovered data will be converted to its type \n
                       its size should be equal to the nb of lines*nb of elts to recover
  \param   a_nbElts : number of elements to recover (on each line)
  \param   a_nbLines : number of lines to recover
  \param   a_mandatoryFlag : boolean indicating the data to recover is  \n
                             mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found) \n
                             or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \brief   Look for "a_nbLines" lines of "a_nbElts" elts in the 'a_file' file matching the "a_keyword" string  \n
           passed as parameter and return the corresponding value(s) in the "ap_return" templated 1D array. \n
           This function expects the following parsing : \n
           KEY :  \n
           line1 : elt1,elt2,...,eltN \n
           line2 : elt1,elt2,...,eltN \n
           (...) \n
           lineN : elt1,elt2,...,eltN
  \details This function assumes the following separators : \n
           ":" is used as separator between the keyworld and the value \n
           "#" is used for comment. Every following characters will be discarded \n
           "," is used to separate elements if more than one are required
  \return  0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<typename T> int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);

/*!
  \fn      int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag)
  \param   a_file : string containing the path to the ASCII file
  \param   a_keyword : key related to the data we want to recover
  \param   a2p_return : 2D templated array in which the data will be returned \n
                       the recovered data will be converted to its type \n
                       its size should be equal [nb of lines][nb of elts to recover]
  \param   a_nbElts : number of elements to recover (on each line)
  \param   a_nbLines : number of lines to recover
  \param   a_mandatoryFlag : boolean indicating the data to recover is  \n
                             mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found) \n
                             or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \brief   Look for "a_nbLines" lines of "a_nbElts" elts in the 'a_file' file matching the "a_keyword" string  \n
           passed as parameter and return the corresponding value(s) in the "a2p_return" 2D array. \n
           The first and second dimensions of "a2p_return" correspond to the line and elements respectively \n
           This function expects the following parsing : \n
           KEY :  \n
           line1 : elt1,elt2,...,eltN \n
           line2 : elt1,elt2,...,eltN \n
           (...) \n
           lineN : elt1,elt2,...,eltN
  \details This function assumes the following separators : \n
           ":" is used as separator between the keyworld and the value \n
           "#" is used for comment. Every following characters will be discarded \n
           "," is used to separate elements if more than one are required
  \return  0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<typename T> int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);

/*!
  \fn      int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag)
  \param   a_file : string containing the path to the ASCII file
  \param   a_keyword : key related to the data we want to recover
  \param   ap_return : templated array in which the data will be returned \n
                       the recovered data will be converted to its type \n
                       its size should be equal to the number of elts to recover
  \param   a_nbElts : number of elements to recover
  \param   a_mandatoryFlag : boolean indicating the data to recover is  \n
                             mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found) \n
                             or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \param   a_firstTag : the requested key will be looked for after the first occurence of this string
  \param   a_lastTag : the requested key will be looked for before the first occurence of this string
  \brief   Look for "a_nbElts" elts in the "a_file" file matching the "a_keyword" string passed as parameter \n
           and return the corresponding value(s) in the "ap_return" templated array. \n
           Additionnal two parameters allow to search within two specific tags \n
           This function expects the following parsing : \n
           FIRST_TAG \n
           (...) \n
           KEY : elt1,elt2,...,eltN  \n
           (...) \n
           LAST_TAG
  \details This function assumes the following separators : \n
           ":" is used as separator between the keyworld and the value \n
           "#" is used for comment. Every following characters will be discarded \n
           "," is used to separate elements if more than one are required
  \return  0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<typename T> int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);

/*!
  \fn      int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag)
  \param   a_file : string containing the path to the ASCII file
  \param   a_keyword : key related to the data we want to recover
  \param   a2p_return : 2D templated array in which the data will be returned \n
                       the recovered data will be converted to its type \n
                       its size should be equal [nb of lines][nb of elts to recover]
  \param   a_nbElts : number of elements to recover (on each line)
  \param   a_nbLines : number of lines to recover
  \param   a_mandatoryFlag : boolean indicating the data to recover is  \n
                             mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found) \n
                             or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \param   a_firstTag : the requested key will be looked for after the first occurence of this string
  \param   a_lastTag : the requested key will be looked for before the first occurence of this string
  \brief   Look for "a_nbLines" lines of "a_nbElts" elts in the "a_file" file matching the "a_keyword" string  \n
           passed as parameter and return the corresponding value(s) in the "a2p_return" templated 2D array. \n
           The first and second dimensions of "a2p_return" correspond to the line and elements respectively \n
           Additionnal two parameters allow to search within two specific tags \n
           This function expects the following parsing : \n
           KEY : \n
           FIRST_TAG \n
           (...) \n
           line1 : elt1,elt2,...,eltN \n
           line2 : elt1,elt2,...,eltN \n
           (...) \n
           lineN : elt1,elt2,...,eltN \n
           (...) \n
           LAST_TAG
  \details This function assumes the following separators : \n
           ":" is used as separator between the keyworld and the value \n
           "#" is used for comment. Every following characters will be discarded \n
           "," is used to separate elements if more than one are required
  \return  0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<typename T> int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
// Functions for data conversion
/*!
  \fn      int ConvertFromString(const string& a_str, string* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Copy the 'a_str' string in the position pointed by 'a_result'
  \details The only purposes of this function is to have an            unified templated conversion function for each type
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, string* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, float* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in float and copy the result in the variable pointed by 'a_result'
  \details Uses strtod to check errors with str to float conversion \n
           (implementation similar to c++11 std::stof() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, float* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, double* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in double and copy the result in the position pointed by 'a_result'
  \details Uses strtod to check errors with str to double conversion \n
           (implementation similar to c++11 std::stod() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, double* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, long double* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in long double and copy the result in the position pointed by 'a_result'
  \details Uses strtod to check errors with str to ldouble conversion \n
           (implementation similar to c++11 std::stod() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, long double* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, int* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in int and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to int conversion \n
           (implementation similar to c++11 std::stoi() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, int* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, int64_t* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in long int and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to int64_t (8 bytes int) conversion \n
           (implementation similar to c++11 std::stol() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, int64_t* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, uint8_t* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in uint16_t and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to uint8 conversion \n
           (implementation similar to c++11 std::stoi() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, uint8_t* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, uint16_t* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in uint16_t and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to uint16 conversion \n
           (implementation similar to c++11 std::stoi() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, uint16_t* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, uint32_t* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in uint32_t and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to uint32 conversion \n
           (implementation similar to c++11 std::stoi() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, uint32_t* a_result);

/*!
  \fn      int ConvertFromString(const string& a_str, bool* a_result)
  \param   a_str : string to convert
  \param   a_result : variable which will recover the result
  \brief   Convert the 'a_str' string in bool and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to bool conversion \n
           (implementation similar to c++11 std::stoi() ).
  \return  0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, bool* a_result);


/*!
  \fn      bool FLTNBIsEqual(FLTNB a, FLTNB b, FLTNB a_eps)
  \param   a : 1st FLTNB nb to compare
  \param   b : 2nd FLTNB nb to compare
  \param   a_eps : epsilon for comparison
  \brief   Comparison of FLTNB numbers
  \return  true if equal according to the provided epsilon, false otherwise
*/
bool FLTNBIsEqual(FLTNB a, FLTNB b, FLTNB a_eps);


/*!
  \fn      string GetFileFromPath(const string& a_pathToFile)
  \param   a_pathToFile
  \brief   Simply return the file from a path string passed in parameter
  \return  The path.
*/
string GetFileFromPath(const string& a_pathToFile);

/*!
  \fn      string GetPathOfFile(const string& a_pathToFile)
  \param   a_pathToFile
  \brief   Simply return the path to the directory of a file path string passed in parameter
  \return  The path.
*/
string GetPathOfFile(const string& a_pathToFile);

/*!
  \fn      string ConvertAllSlashOcurrencesToBackSlash()
  \param   const string& a_path
  \brief   Simply convert all occurrences of "/" to "\"
  \return  A new converted string
*/
string ConvertAllSlashOccurrencesToBackSlash(const string& a_path);

#endif
