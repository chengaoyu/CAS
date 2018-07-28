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

  \brief Implementation of class gOptions functions
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "sOutputManager.hh"
// For error handling
#include <errno.h> 
#include <limits.h>
#ifdef _WIN32
// Avoid compilation errors due to mix up between std::min()/max() and
// min max macros
#undef min
#undef max
#endif

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadStringOption
  \param a_input : string to parse
  \param ap_return : (templated) array in which the parsed elements will be returned
  \param a_nbElts : a number of elements to parse
  \param sep : the separator (usually comma)
  \param a_option : string indicating from where the function has been called
  \brief Parse the 'a_input' string corresponding to the 'a_option' into 'a_nbElts' elements, using the 'sep' separator. 
         The results are returned in the templated 'ap_return' dynamic templated array. 
         Call "ConvertFromString()" to perform the correct conversion depending on the type of the data to convert
  \return 0 if success, and positive value otherwise.
*/
template<class T> 
int ReadStringOption(const string& a_input, T* ap_return, int a_nbElts, const string& sep, const string& a_option)
{

  size_t pos = 0;
  size_t pos2 = a_input.find_first_of(sep, 0);

  if (a_nbElts>1 && pos2 == string::npos)
  {
    Cerr("***** gOptions::ReadStringOption() -> Error : '" << sep << "' not found in option: " << a_option << endl);
    return 1;
  }
  else
  {
    for (int i=0 ; i<a_nbElts ; i++)
    {
      string elt = a_input.substr(pos, pos2-pos); //substr(position, nbEltToRead)
      
      if(ConvertFromString(elt, &ap_return[i]))
      {
        Cerr("***** gOptions::ReadStringOption() -> Error when trying to read input data for option: " << a_option <<  endl);
        return 1;
      } 

      pos = pos2+1;
      pos2 = a_input.find_first_of(sep, pos);

      // TODO Errors in case of too many/few elements in the string
    }
  }
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadDataASCIIFile
  \param a_file : string containing the path to the ASCII file
  \param a_keyword : key related to the data we want to recover
  \param ap_return : templated array in which the data will be returned
                     the recovered data will be converted to its type
                     its size should be equal to the number of elts to recover
  \param a_nbElts : number of elements to recover
  \param a_mandatoryFlag : boolean indicating the data to recover is 
                           mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found)
                        or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \brief Look for "a_nbElts" elts in the "a_file" file matching the "a_keyword" string passed as parameter 
         and return the corresponding value(s) in the "ap_return" templated array.
         This function expects the following parsing :
         KEY : elt1,elt2,...,eltN
  \details This function assumes the following separators :
           ":" is used as separator between the keyworld and the value
           "#" is used for comment. Every following characters will be discarded
           "," is used to separate elements if more than one are required
  \return 0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<class T> 
int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, bool a_mandatoryFlag)
{
  ifstream input_file(a_file.c_str(), ios::in);
  string line;
  string sep = ":";
  string sep_comment = "#";
  string sep_elt = ",";
  
  // Check file
  if (input_file)
  {
    while(!input_file.eof())
    {
      getline(input_file, line);
  
      //remove comment
      if (line.find(sep_comment) != string::npos) line = line.substr(0, line.find_first_of(sep_comment)) ;
      
      if (line.find(a_keyword) != string::npos)
      //if (line.compare(a_keyword) == 0)
      {
        //remove field in string
        line = line.substr(line.find_first_of(sep)+1);

        //clear every spaces in the line string;
        //line.erase(remove_if(line.begin(), line.end(), (int(*)(int))isspace), line.end()); // Explicit type (int(*)(int)) is required to tell compiler which function to take the address of.

        // Erase all blank stuff before first character
        line.erase(0, line.find_first_not_of(" !\t\r\n")); // Erase all blank stuff before first character
        line.erase(line.find_last_not_of(" \t\r\n")+1 , line.length());
        
        size_t pos = 0;
        size_t pos2 = line.find_first_of(sep_elt, pos);

        // Check if separators were found
        if (a_nbElts>1 && pos2 == string::npos)
        {
          Cerr("***** gOptions::ReadDataASCIIFile() -> The required separator : '" << sep_elt << "' not found for tag: " << a_keyword << endl);
          return 1;
        }
        else // Read one element, or several (a_nbElts) elements separated with ','
        {
          for (int i=0 ; i<a_nbElts ; i++)
          {
            // Check if we reach the end of the line before initializing all elements
            // SS: corrected this incorrect line: if(pos<0)
            if (pos==string::npos)
            {
              Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file '" << a_file << "'." << endl);
              Cerr("*****                                  Expected to read " << a_nbElts << " elements, but only " << i+1 << " were found." << endl);
              return 1;
            }
            
            // Parse the line only if more than 1 elt are requested. Just pick the line otherwise
            string elt = (a_nbElts>1) ? line.substr(pos, pos2-pos) : line ;
            
            if (ConvertFromString(elt, &ap_return[i]))
            {
              Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file '" << a_file << "'." << endl);
              return 1;
            }
            
            
            //pos = pos2+1;
            // return string::npos if pos2 not found, meaning we reach the end of the line
            // SS: corrected this incorrect line pos = (pos2<0) ? -1 : pos2+1 ;
            pos = (pos2==string::npos) ? string::npos : pos2+1;
            pos2 = line.find_first_of(sep_elt, pos);
            
          }
        }
        return 0;
      }
  
    }
    // Throw an error message if the tag is mandatory
    if (a_mandatoryFlag == true) 
    {
      Cerr("***** gOptions::ReadDataASCIIFile() -> Error when reading file '" << a_file << "'. Tag '" << a_keyword << "' was not found." << endl);
      return KEYWORD_MANDATORY_NOT_FOUND;
    }
    else
    {
      return KEYWORD_OPTIONAL_NOT_FOUND;
    }
  }
  else
  {
    Cerr("***** gOptions::ReadDataASCIIFile() -> Couldn't find or read data-file '"<< a_file << "' !" << endl);
    return 1;
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadDataASCIIFile
  \param a_file : string containing the path to the ASCII file
  \param a_keyword : key related to the data we want to recover
  \param ap_return : templated array in which the data will be returned
                     the recovered data will be converted to its type
                     its size should be equal to the nb of lines*nb of elts to recover
  \param a_nbElts : number of elements to recover (on each line)
  \param a_nbLines : number of lines to recover
  \param a_mandatoryFlag : boolean indicating the data to recover is 
                           mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found)
                        or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \brief Look for "a_nbLines" lines of "a_nbElts" elts in the 'a_file' file matching the "a_keyword" string 
         passed as parameter and return the corresponding value(s) in the "ap_return" templated 1D array.
         This function expects the following parsing :
         KEY : 
         line1 : elt1,elt2,...,eltN
         line2 : elt1,elt2,...,eltN
         (...)
         lineN : elt1,elt2,...,eltN
  \details This function assumes the following separators :
           ":" is used as separator between the keyworld and the value
           "#" is used for comment. Every following characters will be discarded
           "," is used to separate elements if more than one are required
  \return 0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<class T> 
int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag)
{  
  ifstream input_file(a_file.c_str(), ios::in);
  string line;
  string sep = ":";
  string sep_comment = "#";
  string sep_elt = ",";

  // Check file
  if (input_file)
  {
    while (!input_file.eof())
    {
      getline(input_file, line);
  
      //remove comment
      if (line.find(sep_comment)) line = line.substr(0, line.find_first_of(sep_comment));
  
      if (line.find(a_keyword) != string::npos)
      //if (line.compare(a_keyword) == 0)
      {
        for (int l=0 ; l<a_nbLines ; l++)
        {
          getline(input_file, line);
          //remove field in string
          line = line.substr(line.find_first_of(sep)+1);
  
          //clear the space before the first element in the line string;
          //line.erase(remove_if(line.begin(), line.end(), (int(*)(int))isspace), line.end()); // Explicit type (int(*)(int)) is required to tell compiler which function to take the address of.

          // Erase all blank stuff before first character
          line.erase(0, line.find_first_not_of(" !\t\r\n")); // Erase all blank stuff before first character
          line.erase(line.find_last_not_of(" \t\r\n")+1 , line.length());
        
          size_t pos = 0;
          size_t pos2 = line.find_first_of(sep_elt, pos);
              
          // Check if separators were found
          if (a_nbElts>1 && pos2 == string::npos)
          {
            Cerr("***** gOptions::ReadDataASCIIFile() -> The required separator : '" << sep_elt << "' not found for tag: " << a_keyword << endl);
            return 1;
          }
          else // Read one element, or several (a_nbElts) elements separated with ','
          {
            for (int i=0 ; i<a_nbElts ; i++)
            {
              // Check if we reach the end of the line before initializing all elements
              // SS: corrected this incorrect line: if(pos<0)
              if (pos==string::npos)
              {
                Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file '" << a_file << "'." << endl);
                Cerr("*****                                  Expected to read " << a_nbElts << " elements, but only " << i+1 << " were found." << endl);
                return 1;
              }
              
              // Parse the line only if more than 1 elt are requested. Just pick the line otherwise
              string elt = (a_nbElts>1) ? line.substr(pos, pos2-pos) : line ;
    
              if (ConvertFromString(elt, &ap_return[l*a_nbElts+i]))
              {
                Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file " << a_file <<  endl);
                return 1;
              } 
    
              //pos = pos2+1;
              // return string::npos if pos2 not found, meaning we reach the end of the line
              // SS: corrected this incorrect line pos = (pos2<0) ? -1 : pos2+1 ;
              pos = (pos2==string::npos) ? string::npos : pos2+1;
              pos2 = line.find_first_of(sep_elt, pos);
            }
          }
        }  
        return 0; 
      }
    }
    
    // Throw an error message if the tag is mandatory
    if (a_mandatoryFlag == true) 
    {
      Cerr("***** gOptions::ReadDataASCIIFile() -> Error when reading file '" << a_file << "'. Tag '" << a_keyword << "' was not found." << endl);
      return KEYWORD_MANDATORY_NOT_FOUND;
    }
    else
    {
      return KEYWORD_OPTIONAL_NOT_FOUND;
    }
  }
  else
  {
    Cerr("***** gOptions::ReadDataASCIIFile() -> Couldn't find or read data-file '" << a_file << "' !" << endl);
    return 1;
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadDataASCIIFile
  \param a_file : string containing the path to the ASCII file
  \param a_keyword : key related to the data we want to recover
  \param a2p_return : 2D templated array in which the data will be returned
                     the recovered data will be converted to its type
                     its size should be equal [nb of lines][nb of elts to recover]
  \param a_nbElts : number of elements to recover (on each line)
  \param a_nbLines : number of lines to recover
  \param a_mandatoryFlag : boolean indicating the data to recover is 
                           mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found)
                        or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \brief Look for "a_nbLines" lines of "a_nbElts" elts in the 'a_file' file matching the "a_keyword" string 
         passed as parameter and return the corresponding value(s) in the "a2p_return" 2D array.
         The first and second dimensions of "a2p_return" correspond to the line and elements respectively
         This function expects the following parsing :
         KEY : 
         line1 : elt1,elt2,...,eltN
         line2 : elt1,elt2,...,eltN
         (...)
         lineN : elt1,elt2,...,eltN
  \details This function assumes the following separators :
           ":" is used as separator between the keyworld and the value
           "#" is used for comment. Every following characters will be discarded
           "," is used to separate elements if more than one are required
  \return 0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<class T> 
int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag)
{
  ifstream input_file(a_file.c_str(), ios::in);
  string line;
  string sep = ":";
  string sep_comment = "#";
  string sep_elt = ",";

  // Check file
  if (input_file)
  {
    while (!input_file.eof())
    {
      getline(input_file, line);
  
      //remove comment
      if (line.find(sep_comment)) line = line.substr(0, line.find_first_of(sep_comment));
  
      if (line.find(a_keyword) != string::npos)
      //if (line.compare(a_keyword) == 0)
      {
        for (int l=0 ; l<a_nbLines ; l++)
        {
          getline(input_file, line);
          //remove field in string
          line = line.substr(line.find_first_of(sep)+1);
  
          //clear the space before the first element in the line string;
          //line.erase(remove_if(line.begin(), line.end(), (int(*)(int))isspace), line.end()); // Explicit type (int(*)(int)) is required to tell compiler which function to take the address of.

          // Erase all blank stuff before first character
          line.erase(0, line.find_first_not_of(" !\t\r\n")); // Erase all blank stuff before first character
          line.erase(line.find_last_not_of(" \t\r\n")+1 , line.length());
          
          size_t pos = 0;
          size_t pos2 = line.find_first_of(sep_elt, pos);
              
          // Check if separators were found
          if (a_nbElts>1 && pos2 == string::npos)
          {
            Cerr("***** gOptions::ReadDataASCIIFile() -> The required separator : '" << sep_elt << "' not found for tag: " << a_keyword << endl);
            return 1;
          }
          else // Read one element, or several (a_nbElts) elements separated with ','
          {
            for (int i=0 ; i<a_nbElts ; i++)
            {
              // Check if we reach the end of the line before initializing all elements
              // SS: corrected this incorrect line: if(pos<0)
              if (pos==string::npos)
              {
                Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file '" << a_file << "'." << endl);
                Cerr("*****                                  Expected to read " << a_nbElts << " elements, but only " << i+1 << " were found." << endl);
                return 1;
              }
              
              // Parse the line only if more than 1 elt are requested. Just pick the line otherwise
              string elt = (a_nbElts>1) ? line.substr(pos, pos2-pos) : line ;
    
              if (ConvertFromString(elt, &a2p_return[l][i]))
              {
                Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file " << a_file <<  endl);
                return 1;
              } 
    
              //pos = pos2+1;
              // return string::npos if pos2 not found, meaning we reach the end of the line
              // SS: corrected this incorrect line pos = (pos2<0) ? -1 : pos2+1 ;
              pos = (pos2==string::npos) ? string::npos : pos2+1;
              pos2 = line.find_first_of(sep_elt, pos);
            }
          }
        }  
        return 0; 
      }
    }
    
    // Throw an error message if the tag is mandatory
    if (a_mandatoryFlag == true ) 
    {
      Cerr("***** gOptions::ReadDataASCIIFile() -> Error when reading file '" << a_file << "'. Tag '" << a_keyword << "' was not found." << endl);
      return KEYWORD_MANDATORY_NOT_FOUND;
    }
    else
    {
      return KEYWORD_OPTIONAL_NOT_FOUND;
    }
  }
  else
  {
    Cerr("***** gOptions::ReadDataASCIIFile() -> Couldn't find or read data-file '" << a_file << "' !" << endl);
    return 1;
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/*
  \fn ReadDataASCIIFile
  \param a_file : string containing the path to the ASCII file
  \param a_keyword : key related to the data we want to recover
  \param ap_return : templated array in which the data will be returned
                     the recovered data will be converted to its type
                     its size should be equal to the number of elts to recover
  \param a_nbElts : number of elements to recover
  \param a_mandatoryFlag : boolean indicating the data to recover is 
                           mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found)
                        or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \param a_firstTag : the requested key will be looked for after the first occurence of this string
  \param a_lastTag : the requested key will be looked for before the first occurence of this string
  \brief Look for "a_nbElts" elts in the "a_file" file matching the "a_keyword" string passed as parameter
         and return the corresponding value(s) in the "ap_return" templated array.
         Additionnal two parameters allow to search within two specific tags
         This function expects the following parsing :
         FIRST_TAG
         (...)
         KEY : elt1,elt2,...,eltN 
         (...)
         LAST_TAG
  \details This function assumes the following separators :
           ":" is used as separator between the keyworld and the value
           "#" is used for comment. Every following characters will be discarded
           "," is used to separate elements if more than one are required
  \return 0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<class T> 
int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag)
{
  ifstream input_file(a_file.c_str(), ios::in);
  string line;
  string sep = ":";
  string sep_comment = "#";
  string sep_elt = ",";
  bool search_on = false;
  
  // Check file
  if(input_file)
  {
    while(!input_file.eof())
    {
      getline(input_file, line);
  
      //remove comment
      if (line.find(sep_comment) != string::npos) line = line.substr(0, line.find_first_of(sep_comment)) ;
      
      if( line.find(a_firstTag) != string::npos) search_on = true;
      if( line.find(a_lastTag) != string::npos || line.find("eof") != string::npos ) search_on = false;
      
      if( search_on )
      {      
        if (line.find(a_keyword) != string::npos)
        //if (line.compare(a_keyword) == 0)
        {
          //remove field in string
          line = line.substr(line.find_first_of(sep)+1);

          //clear the space before the first element in the line string;
          //line.erase(remove_if(line.begin(), line.end(), (int(*)(int))isspace), line.end()); // Explicit type (int(*)(int)) is required to tell compiler which function to take the address of.

          // Erase all blank stuff before first character
          line.erase(0, line.find_first_not_of(" !\t\r\n")); // Erase all blank stuff before first character
          line.erase(line.find_last_not_of(" \t\r\n")+1 , line.length());
          
          size_t pos = 0;
          size_t pos2 = line.find_first_of(sep_elt, pos);

          // Check if separators were found
          if (a_nbElts>1 && pos2 == string::npos)
          {
            Cerr("***** gOptions::ReadDataASCIIFile() -> The required separator : '" << sep_elt << "' not found for tag: " << a_keyword << endl);
            return 1;
          }
          else // Read one element, or several (a_nbElts) elements separated with ','
          {
            for (int i=0 ; i<a_nbElts ; i++)
            {
              
              // Check if we reach the end of the line before initializing all elements
              // SS: corrected this incorrect line: if(pos<0)
              if (pos==string::npos)
              {
                Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file '" << a_file << "'." << endl);
                Cerr("*****                                  Expected to read " << a_nbElts << " elements, but only " << i+1 << " were found." << endl);
                return 1;
              }
            
              // Parse the line only if more than 1 elt are requested. Just pick the line otherwise
              string elt = (a_nbElts>1) ? line.substr(pos, pos2-pos) : line ;
            
              if(ConvertFromString(elt, &ap_return[i]))
              {
                Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file '" << a_file << "'." << endl);
                return 1;
              }
            
              //pos = pos2+1;
              // return string::npos if pos2 not found, meaning we reach the end of the line
              // SS: corrected this incorrect line pos = (pos2<0) ? -1 : pos2+1 ;
              pos = (pos2==string::npos) ? string::npos : pos2+1;
              pos2 = line.find_first_of(sep_elt, pos);
            }
          }
          return 0;
        }
      }
      
      
      
      
  
    }
    // Throw an error message if the tag is mandatory
    if(a_mandatoryFlag == true) 
    {
      Cerr("***** gOptions::ReadDataASCIIFile() -> Error when reading file '" << a_file << "'. Tag '" << a_keyword << "' was not found." << endl);
      return KEYWORD_MANDATORY_NOT_FOUND;
    }
    else
    {
      return KEYWORD_OPTIONAL_NOT_FOUND;
    }
  }
  else
  {
    Cerr("***** gOptions::ReadDataASCIIFile() -> Couldn't find or read data-file '"<< a_file << "' !" << endl);
    return 1;
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadDataASCIIFile
  \param a_file : string containing the path to the ASCII file
  \param a_keyword : key related to the data we want to recover
  \param a2p_return : 2D templated array in which the data will be returned
                     the recovered data will be converted to its type
                     its size should be equal [nb of lines][nb of elts to recover]
  \param a_nbElts : number of elements to recover (on each line)
  \param a_nbLines : number of lines to recover
  \param a_mandatoryFlag : boolean indicating the data to recover is 
                           mandatory (KEYWORD_MANDATORY) (error value (=1) will be returned if not found)
                        or optional (KEYWORD_OPTIONAL) (warning value (=2) will be returned if not found)
  \param a_firstTag : the requested key will be looked for after the first occurence of this string
  \param a_lastTag : the requested key will be looked for before the first occurence of this string
  \brief Look for "a_nbLines" lines of "a_nbElts" elts in the "a_file" file matching the "a_keyword" string 
         passed as parameter and return the corresponding value(s) in the "a2p_return" templated 2D array.
         The first and second dimensions of "a2p_return" correspond to the line and elements respectively
         Additionnal two parameters allow to search within two specific tags
         This function expects the following parsing :
         KEY :
          FIRST_TAG
         (...)
         line1 : elt1,elt2,...,eltN
         line2 : elt1,elt2,...,eltN
         (...)
         lineN : elt1,elt2,...,eltN
         (...)
         LAST_TAG
  \details This function assumes the following separators :
           ":" is used as separator between the keyworld and the value
           "#" is used for comment. Every following characters will be discarded
           "," is used to separate elements if more than one are required
  \return 0 if success, and positive value otherwise (1 if error, 2 if tag not found).
*/
template<class T> 
int ReadDataASCIIFile(const string& a_file, const string& a_keyword, T** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag)
{
  ifstream input_file(a_file.c_str(), ios::in);
  string line;
  string sep = ":";
  string sep_comment = "#";
  string sep_elt = ",";
  bool search_on = false;
  
  // Check file
  if(input_file)
  {
    while(!input_file.eof())
    {
      getline(input_file, line);
  
      //remove comment
      if (line.find(sep_comment) != string::npos) line = line.substr(0, line.find_first_of(sep_comment)) ;
      
      if( line.find(a_firstTag) != string::npos) search_on = true;
      if( line.find(a_lastTag) != string::npos || line.find("eof") != string::npos ) search_on = false;

      if( search_on )
      {      
          for (int l=0 ; l<a_nbLines ; l++)
          {
            getline(input_file, line);
            //remove field in string
            line = line.substr(line.find_first_of(sep)+1);
    
            //clear the space before the first element in the line string;
            //line.erase(remove_if(line.begin(), line.end(), (int(*)(int))isspace), line.end()); // Explicit type (int(*)(int)) is required to tell compiler which function to take the address of.

            // Erase all blank stuff before and after first character
            line.erase(0, line.find_first_not_of(" !\t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n")+1 , line.length());
  
            size_t pos = 0;
            size_t pos2 = line.find_first_of(sep_elt, pos);
                
            // Check if separators were found
            if (a_nbElts>1 && pos2 == string::npos)
            {
              Cerr("***** gOptions::ReadDataASCIIFile() -> The required separator : '" << sep_elt << "' not found for tag: " << a_keyword << endl);
              return 1;
            }
            else // Read one element, or several (a_nbElts) elements separated with ','
            {
              for (int i=0 ; i<a_nbElts ; i++)
              {
                // Check if we reach the end of the line before initializing all elements
                // SS: corrected this incorrect line: if(pos<0)
                if (pos==string::npos)
                {
                  Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file '" << a_file << "'." << endl);
                  Cerr("*****                                  Expected to read " << a_nbElts << " elements, but only " << i+1 << " were found." << endl);
                  return 1;
                }
                
                // Parse the line only if more than 1 elt are requested. Just pick the line otherwise
                string elt = (a_nbElts>1) ? line.substr(pos, pos2-pos) : line ;

                if (ConvertFromString(elt, &a2p_return[l][i]))
                {
                  Cerr("***** gOptions::ReadDataASCIIFile() -> Exception when trying to read tag '" << a_keyword << "' in file " << a_file <<  endl);
                  return 1;
                } 
      
                //pos = pos2+1;
                // return string::npos if pos2 not found, meaning we reach the end of the line
                // SS: corrected this incorrect line pos = (pos2<0) ? -1 : pos2+1 ;
                pos = (pos2==string::npos) ? string::npos : pos2+1;
                pos2 = line.find_first_of(sep_elt, pos);
              }
            }
            
          }
          return 0;
      }  
    }
    // Throw an error message if the tag is mandatory
    if(a_mandatoryFlag == true) 
    {
      Cerr("***** gOptions::ReadDataASCIIFile() -> Error when reading file '" << a_file << "'. Tag '" << a_keyword << "' was not found." << endl);
      return KEYWORD_MANDATORY_NOT_FOUND;
    }
    else
    {
      return KEYWORD_OPTIONAL_NOT_FOUND;
    }
  }
  else
  {
    Cerr("***** gOptions::ReadDataASCIIFile() -> Couldn't find or read data-file '"<< a_file << "' !" << endl);
    return 1;
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Copy the 'a_str' string in the position pointed by 'a_result'
  \details The only purposes of this function is to have an 
           unified templated conversion function for each type
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, string* a_result)
{
  *a_result = a_str.c_str();
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in float and copy the result in the variable pointed by 'a_result'
  \details Uses strtod to check errors with str to float conversion
          (implementation similar to c++11 std::stof() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, float* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  double val = strtod(p, &end);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into float " << endl);
    return 1;
  }
  if (errno == ERANGE) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into float " << endl);
    return 1;
  }
  
  *a_result = static_cast<float>(val);
  
  return 0;  
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in double and copy the result in the position pointed by 'a_result'
  \details Uses strtod to check errors with str to double conversion
          (implementation similar to c++11 std::stod() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, double* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  double val = strtod(p, &end);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into double " << endl);
    return 1;
  }
  if (errno == ERANGE) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into double " << endl);
    return 1;
  }
  
  *a_result = val;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in long double and copy the result in the position pointed by 'a_result'
  \details Uses strtod to check errors with str to ldouble conversion
          (implementation similar to c++11 std::stod() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, long double* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  long double val = strtold(p, &end);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into double " << endl);
    return 1;
  }
  if (errno == ERANGE) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into double " << endl);
    return 1;
  }
  
  *a_result = val;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in int and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to int conversion
          (implementation similar to c++11 std::stoi() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, int* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  int64_t val = strtol(p, &end, 10);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into int " << endl);
    return 1;
  }
  if (errno==ERANGE || val<INT_MIN || val>INT_MAX) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into int " << endl);
    return 1;
  }
  
  *a_result = static_cast<int>(val);
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in int64_t and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to lint conversion
          (implementation similar to c++11 std::stol() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, int64_t* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  int64_t val = strtol(p, &end, 10);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into int64_t " << endl);
    return 1;
  }
  if (errno == ERANGE) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into int64_t " << endl);
    return 1;
  }
  
  *a_result = val;

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in uint16_t and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to uint8 conversion
          (implementation similar to c++11 std::stoi() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, uint8_t* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  int64_t val = strtol(p, &end, 10);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into uint16 " << endl);
    return 1;
  }
  if (errno == ERANGE || val<0 || val>UCHAR_MAX) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into uint16 " << endl);
    return 1;
  }
  
  *a_result = val;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in uint16_t and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to uint16 conversion
          (implementation similar to c++11 std::stoi() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, uint16_t* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  int64_t val = strtol(p, &end, 10);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into uint16 " << endl);
    return 1;
  }
  if (errno == ERANGE || val<0 || val>USHRT_MAX) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into uint16 " << endl);
    return 1;
  }
  
  *a_result = val;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in uint32_t and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to uint32 conversion
          (implementation similar to c++11 std::stoi() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, uint32_t* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  int64_t val = strtol(p, &end, 10);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into uint32 " << endl);
    return 1;
  }
  if (errno == ERANGE || val<0 || val>UINT_MAX) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into uint32 " << endl);
    return 1;
  }
  
  *a_result = val;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertFromString
  \param a_str : string to convert
  \param a_result : variable which will recover the result
  \brief Convert the 'a_str' string in bool and copy the result in the position pointed by 'a_result'
  \details Uses strtol to check errors with str to bool conversion
        (implementation similar to c++11 std::stoi() ).
  \return 0 if success, and positive value otherwise
*/
int ConvertFromString(const string& a_str, bool* a_result)
{
  const char* p = a_str.c_str();
  char* end;
  errno = 0;
  int64_t val = strtol(p, &end, 10);
  
  if (p == end) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Invalid argument exception while trying to convert '" << a_str << "' into bool " << endl);
    return 1;
  }
  if (errno == ERANGE || val<0 || val>1) 
  {
    Cerr("***** gOptions::ConvertFromString() -> Out of range exception while trying to convert '" << a_str << "' into bool " << endl);
    return 1;
  }
  
  *a_result = val;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetFileFromPath()
  \param a_pathToFile
  \brief Simply return the file from a path string passed in parameter
  \return The path.
*/
string GetFileFromPath(const string& a_pathToFile)
{
  string path = a_pathToFile;

  int pos = path.find_last_of(OS_SEP);
  if (path.find_last_of(OS_SEP) == string::npos)
    return path;
    
  path = path.substr(pos+1);
  return path;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetPathOfFile()
  \param a_pathToFile
  \brief Simply return the path to the directory of a file path string passed in parameter
  \return The path.
*/
string GetPathOfFile(const string& a_pathToFile)
{
  string path = a_pathToFile;
  
  int pos = path.find_last_of(OS_SEP);
  if (path.find_last_of(OS_SEP) == string::npos)
    return "";
  
  path = path.substr(0,pos+1);

  return path;
}


/*
  \fn      string ConvertAllSlashOcurrencesToBackSlash()
  \param   const string& a_path
  \brief   Simply convert all occurrences of "/" to "\"
  \return  A new converted string
*/
string ConvertAllSlashOccurrencesToBackSlash(const string& a_path)
{
  string result = a_path;
  size_t position;
  while ( (position  = result.find_first_of("/")) != string::npos )
  {
    result.replace(position,1,"\\");
  }
  return result;
}


/*
  \fn      bool FLTNBIsEqual(FLTNB a, FLTNB b, FLTNB a_eps)
  \param   a : 1st FLTNB nb to compare
  \param   b : 2nd FLTNB nb to compare
  \param   a_eps : epsilon for comparison
  \brief   Comparison of FLTNB numbers
  \return  true if equal according to the provided epsilon, false otherwise
*/
bool FLTNBIsEqual(FLTNB a, FLTNB b, FLTNB a_eps)
{
  FLTNB absA = abs(a);
  FLTNB absB = abs(b);
  FLTNB diff = abs(a - b);
  
  if (a == b)
  {
    return true;
  } 
  else if (a == 0 || 
           b == 0 || 
           diff < std::numeric_limits<FLTNB>::min()) 
  { // a or b is 0, or extremely close to it
    return diff < (a_eps * std::numeric_limits<FLTNB>::min() );
  } 
  else 
  { // relative error
    return diff / min((absA + absB), std::numeric_limits<FLTNB>::max()) < a_eps;
  }
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
// Explicit template instantiation 
// (that way no one would use the templated versions of these functions for not yet implemented type (long double, unsigned, ect..) )
template int ReadStringOption<string>(const string& a_input, string* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<int>(const string& a_input, int* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<int64_t>(const string& a_input, int64_t* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<float>(const string& a_input, float* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<double>(const string& a_input, double* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<long double>(const string& a_input, long double* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<uint8_t>(const string& a_input, uint8_t* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<uint16_t>(const string& a_input, uint16_t* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<uint32_t>(const string& a_input, uint32_t* ap_return, int a_nbElts, const string& sep, const string& a_option);
template int ReadStringOption<bool>(const string& a_input, bool* ap_return, int a_nbElts, const string& sep, const string& a_option);

template int ReadDataASCIIFile<string>(const string& a_file, const string& a_keyword, string* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<int>(const string& a_file, const string& a_keyword, int* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<int64_t>(const string& a_file, const string& a_keyword, int64_t* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<float>(const string& a_file, const string& a_keyword, float* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<double>(const string& a_file, const string& a_keyword, double* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<long double>(const string& a_file, const string& a_keyword, long double* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint8_t>(const string& a_file, const string& a_keyword, uint8_t* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint16_t>(const string& a_file, const string& a_keyword, uint16_t* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint32_t>(const string& a_file, const string& a_keyword, uint32_t* ap_return, int a_nbElts, bool a_mandatoryFlag);
template int ReadDataASCIIFile<bool>(const string& a_file, const string& a_keyword, bool* ap_return, int a_nbElts, bool a_mandatoryFlag);

template int ReadDataASCIIFile<string>(const string& a_file, const string& a_keyword, string* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<int>(const string& a_file, const string& a_keyword, int* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<int64_t>(const string& a_file, const string& a_keyword, int64_t* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<float>(const string& a_file, const string& a_keyword, float* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<double>(const string& a_file, const string& a_keyword, double* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<long double>(const string& a_file, const string& a_keyword, long double* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint8_t>(const string& a_file, const string& a_keyword, uint8_t* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint16_t>(const string& a_file, const string& a_keyword, uint16_t* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint32_t>(const string& a_file, const string& a_keyword, uint32_t* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<bool>(const string& a_file, const string& a_keyword, bool* ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);

template int ReadDataASCIIFile<string>(const string& a_file, const string& a_keyword, string** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<int>(const string& a_file, const string& a_keyword, int** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<int64_t>(const string& a_file, const string& a_keyword, int64_t** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<float>(const string& a_file, const string& a_keyword, float** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<double>(const string& a_file, const string& a_keyword, double** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<long double>(const string& a_file, const string& a_keyword, long double** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint8_t>(const string& a_file, const string& a_keyword, uint8_t** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint16_t>(const string& a_file, const string& a_keyword, uint16_t** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<uint32_t>(const string& a_file, const string& a_keyword, uint32_t** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);
template int ReadDataASCIIFile<bool>(const string& a_file, const string& a_keyword, bool** ap_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag);

template int ReadDataASCIIFile<string>(const string& a_file, const string& a_keyword, string* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<int>(const string& a_file, const string& a_keyword, int* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<int64_t>(const string& a_file, const string& a_keyword, int64_t* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<float>(const string& a_file, const string& a_keyword, float* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<double>(const string& a_file, const string& a_keyword, double* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<long double>(const string& a_file, const string& a_keyword, long double* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<uint8_t>(const string& a_file, const string& a_keyword, uint8_t* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<uint16_t>(const string& a_file, const string& a_keyword, uint16_t* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<uint32_t>(const string& a_file, const string& a_keyword, uint32_t* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<bool>(const string& a_file, const string& a_keyword, bool* ap_return, int a_nbElts, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);

template int ReadDataASCIIFile<string>(const string& a_file, const string& a_keyword, string** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<int>(const string& a_file, const string& a_keyword, int** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<int64_t>(const string& a_file, const string& a_keyword, int64_t** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<float>(const string& a_file, const string& a_keyword, float** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<double>(const string& a_file, const string& a_keyword, double** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<long double>(const string& a_file, const string& a_keyword, long double** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<uint8_t>(const string& a_file, const string& a_keyword, uint8_t** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<uint16_t>(const string& a_file, const string& a_keyword, uint16_t** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<uint32_t>(const string& a_file, const string& a_keyword, uint32_t** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
template int ReadDataASCIIFile<bool>(const string& a_file, const string& a_keyword, bool** a2p_return, int a_nbElts, int a_nbLines, bool a_mandatoryFlag, string a_firstTag, string a_lastTag);
