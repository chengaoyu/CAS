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

  \brief Implementation of class gDataConversionUtilities functions
*/

#include "gVariables.hh"
#include "gDataConversionUtilities.hh"
#include <iomanip>

/*
   Miscelleanous functionalities
 */

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckGATECommand
  \param a_key: string containing a GATE command to recover
  \param a_line : string containing the line to check
  \brief Check if the line contains the provided GATE command. In this case, \n
         parse the line and returns the values in a vector of strings
  \details Values are converted in mm if 'cm' is found in the line 
  \return the vector of strings containing the elements of the line.
*/
vector<string> CheckGATECommand(const string& a_key, const string& a_line)
{
  vector<string> values;

  // cut any part after comment symbol
  string line = a_line;

  line = (line.find("#") != string::npos)        ? 
         line.substr(0, line.find_first_of("#")) : 
         line;

  size_t foundAdress = line.find(a_key);

  if (foundAdress != string::npos) 
  {
    values = Split(a_line);
    values.erase (values.begin());
  }
  
  ConvertValuesTomm(values);
  return values;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn split
  \param a_line : string to split
  \brief Split the line provided in parameter into a vector of strings (separator is blankspace)
  \return vector of string containing the splitted elements of the line
*/
vector<string> Split(string a_line)
{
  string buf; 
  stringstream ss(a_line); 

  vector<string> tokens; 

  while (ss >> buf)
    tokens.push_back(buf);

  return tokens;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertValuesTomm
  \param ap_v : vector of strings
  \brief Check if the vector of strings passed in parameter contains the 'cm' unit \n
         In this case, convert all its values in mm
*/
void ConvertValuesTomm(vector<string>& ap_v)
{
  // loop on values
  for(uint16_t i=0; i < ap_v.size(); i++)
    // check if values were provided in cm
    if (ap_v[i] == "cm")
      // Convert all values to mm (skip command key)
      for (int j=0; j<i; j++)
      {
        float val;
        ConvertFromString(ap_v[j], &val);
        ap_v[j] = toString(10 * val);
      }
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn WriteVector
  \param file : the output file
  \param a_key : key to write
  \param a_vals : vector containing the key values
  \brief Write the key and its values in the file provided in parameter
  \return 0 if success, positive value otherwise
*/
template <class T>
int WriteVector(ofstream& file, const string& a_key, vector <T> a_vals)
{  
  int n = a_vals.size();
  
  if (n > 0)
  {
    file << a_key ;
    for (int i=0; i < n; i++)
    {
      stringstream ss;
      ss << a_vals[i];
      if (i == n-1)
        file << ss.str() << endl;
      else
        file << ss.str() << ",";
    }
  }
  else
    file << a_key << "0" <<endl;

  return 0;
}

template int WriteVector(ofstream& file, const string& a_key, vector <double> a_vals);



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn WriteVector
  \param file : the output file
  \param a_key : key to write
  \param a_vals : vector containing the key values
  \brief Write the key and its values in the file provided in parameter
  \return 0 if success, positive value otherwise
*/
int WriteVector(ofstream& file, const string& a_key, vector <string> a_vals)
{  
  int n = a_vals.size();
  
  if (n > 0)
  {
    file << a_key ;
    for (int i=0; i < n; i++)
    {
      if (i == n-1)
        file << a_vals[i] << endl;
      else
        file << a_vals[i] << ",";
    }
  }
  else
    file << a_key << "0" <<endl;

  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn WriteVector
  \param file : the output file
  \param a_key : key to write
  \param a_vals : vector containing the key values
  \brief Write the key and its values contained in a 2 level vector of strings in the file provided in parameter
  \return 0 if success, positive value otherwise
*/
int WriteVector(ofstream& file, const string& a_key, vector <vector<string> > a_vals)
{
  int n = a_vals.size();
  file << a_key;
  for (int i=0; i < n; i++)
    for (int j=0; j < 3; j++)
      if (i == n-1 && j == 2)
        file << a_vals[i][j] << endl;
      else
        file << a_vals[i][j] << ",";
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetGATEMacFiles(const string& a_pathMac, vector<string> ap_pathToMacFiles)
  \param a_pathMac : path to a GATE main macro file
  \param ap_pathToMacFiles : array containing the paths of macro files
  \brief Extract the paths to each macro file contained in the main macro file.     
  \return 0 if success, positive value otherwise
*/
int GetGATEMacFiles(const string& a_pathMac, vector<string> &ap_pathToMacFiles)
{  
  ifstream mac_file(a_pathMac, ios::in);
  
  if(mac_file)
  {
    string quickLine;

    while(getline(mac_file, quickLine))
    {
      vector <string> values;

      values = CheckGATECommand("/control/execute", quickLine);

      if (values.size()>0)
        ap_pathToMacFiles.push_back(GetPathOfFile(a_pathMac)+values[0]);
    }
  }
  else
  {
    Cerr("***** GetGATEMacFiles() -> Couldn't open mac file "<< a_pathMac << " !" << endl);
    return -1;
  }
  
  mac_file.close();
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetGATESystemType
  \param a_pathMac : path to a GATE macro file
  \brief Read a GATE macro file and identify the system type from the 'gate/systems/' command lines
  \return system type, as described by the macro GATE_SYS_TYPE, -1 if error
*/
int GetGATESystemType(const string& a_pathMac)
{
  int system_type = GATE_SYS_UNKNOWN;
  
  // Recover path to all macro files from the main mac file
  vector<string> path_mac_files;
  path_mac_files.push_back(a_pathMac);
  
  if(GetGATEMacFiles(a_pathMac , path_mac_files))
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover paths to GATE macro files !" << endl);
    return -1;
  }

  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
    cout << f << " : " << path_mac_files[f] << endl;
    
  // Loop on all macro files
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream mac_file(path_mac_files[f], ios::in);
    
    if(mac_file)
    {  
      string line;
      while(getline(mac_file, line))
      {
        // Check a 'gate/systems' command line
        if(line.find("/gate/systems/") != string::npos )
        {
          if(line.find("/gate/systems/ecat") != string::npos )
            system_type = GATE_SYS_ECAT;
          
          if(line.find("/gate/systems/cylindricalPET") != string::npos )
            system_type = GATE_SYS_CYLINDRICAL;

          if(line.find("/gate/systems/SPECThead") != string::npos )
            system_type = GATE_SYS_SPECT;
            
          if(line.find("/gate/systems/OPET")          != string::npos ||
             line.find("/gate/systems/CTSCANNER")     != string::npos ||
             line.find("/gate/systems/CPET")          != string::npos ||
             line.find("/gate/systems/ecatAccel")     != string::npos ||
             line.find("/gate/systems/OpticalSystem") != string::npos )
            {
              Cerr("unsupported system detected (line = " << line <<") ! "<< endl);
              Cerr("supported systems for this script are cylindricalPET, SPECThead, and ecat" << endl);
            }
        }
      }
    }
    else
    {
      Cerr("***** GetGATESystemType() -> Error : Couldn't open mac file "<< path_mac_files[f] << " !" << endl);
      return -1;
    }
    
    mac_file.close();
  }
  
  return system_type;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetGATEAliasesCylindrical(vector<string>  path_mac_files
                                string&         rsector_name,
                                string&         module_name,
                                string&         submodule_name,
                                string&         crystal_name,
                                vector<string>& layers_name,
                                int             vb )
  \param path_mac_files
  \param rsector_name
  \param module_name
  \param submodule_name
  \param crystal_name
  \param layers_name
  \param vb
  \brief Loop over a list of path to GATE macro files passed in parameter 
         to recover aliases of the different parts of a cylindricalPET
  \return 0 if success, positive value otherwise
*/
int GetGATEAliasesCylindrical(vector<string>  path_mac_files,
                              string&         rsector_name,
                              string&         module_name,
                              string&         submodule_name,
                              string&         crystal_name,
                              vector<string>& layers_name,
                              int             vb )
{    
  int depth = 0;
  
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream mac_file(path_mac_files[f].c_str(), ios::in);
    
    if(mac_file)
    {
      string quickLine;
      while(getline(mac_file, quickLine))
      {
        vector <string> values;
  
        values = CheckGATECommand("/gate/systems/cylindricalPET/rsector/attach", quickLine);
  
        if (values.size()>0)
        {
          rsector_name = values[0];
          depth++;
        }
  
        values = CheckGATECommand("/gate/systems/cylindricalPET/module/attach", quickLine);
        if (values.size()>0)
        {
          module_name = values[0];
          depth++;
        }
        
        values = CheckGATECommand("/gate/systems/cylindricalPET/submodule/attach", quickLine);
        if (values.size()>0)
        {
          submodule_name = values[0];
          depth++;
        }
        
        values = CheckGATECommand("/gate/systems/cylindricalPET/crystal/attach", quickLine);
        if (values.size()>0)
        {
          crystal_name = values[0];
          depth++;
        } 
  
        // layers
        for(int l=0 ; l<GATE_NB_MAX_LAYERS ; l++)
        {
          stringstream ss;
          ss << l;
          values = CheckGATECommand("/gate/systems/cylindricalPET/layer"+ss.str()+"/attach", quickLine);
          if (values.size()>0)
          {
            layers_name.push_back(values[0]);
            depth++;
          }
        }
          
      }
    }
    else
    {
      Cerr("***** GetGATEAliasesCylindrical()->Couldn't open mac file "<< path_mac_files[f] << " !" << endl);
      return 1;
    }
    
    mac_file.close();
  }


  // Check we have the required number of elements for the system
  if (depth < 2)
  {
    Cout("***** GetGATEAliasesCylindrical() :: Error : Missing elements in the system architecture" << endl <<
         "      At least two of the following lines are required :" << endl <<
         "         - /gate/systems/cylindricalPET/rsector/attach" << endl <<
         "         - /gate/systems/cylindricalPET/module/attach" << endl <<
         "         - /gate/systems/cylindricalPET/submodule/attach" << endl <<
         "         - /gate/systems/cylindricalPET/crystal/attach" << endl <<
         "         - /gate/systems/cylindricalPET/layeri[i=0..3]/attach" << endl);
    return 1;
  }
  else
  {
    // Interpret the first element as the rsector
    if(rsector_name.empty())
    {
      if(module_name.empty())
      {
        //submodule as rsector
        rsector_name = submodule_name;
        submodule_name = "";
      }
      else
      {
        //module as rsector
        rsector_name = module_name;
        module_name = "";
      }
    }
    
    if(vb >= 2)
    {
      if(!rsector_name.empty())   Cout("Detected rsector container's name : " << rsector_name << endl);
      if(!module_name.empty())    Cout("Detected module container's name : " << module_name << endl);
      if(!submodule_name.empty()) Cout("Detected submodule container's name : " << submodule_name << endl);
      if(!crystal_name.empty())   Cout("Detected crystal container's name : " << crystal_name << endl);
      for(size_t l=0 ; l<layers_name.size() ; l++)
        if(!layers_name[l].empty())  Cout("Detected layer #"<< l << " container's name : " << layers_name[l] << endl);
    }
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetGATEAliasesEcat(vector<string> path_mac_files,
                         string&        block_name,
                         string&        crystal_name,
                         int            vb )
  \param path_mac_files
  \param block_name
  \param crystal_name
  \param vb
  \brief Loop over a list of path to GATE macro files passed in parameter 
         to recover aliases of the different parts of an ecat system
  \return 0 if success, positive value otherwise
*/
int GetGATEAliasesEcat(vector<string> path_mac_files,
                       string&        block_name,
                       string&        crystal_name,
                       int            vb )
{
  int depth = 0;
  
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream mac_file(path_mac_files[f].c_str(), ios::in);
        
    if(mac_file)
    {
      // Reading the .mac file line by line and finding the names of the different parts of the architecture
      string quickLine;
      while(getline(mac_file, quickLine))
      {
        vector <string> values;
             
        values = CheckGATECommand("/gate/systems/ecat/block/attach", quickLine);
        if (values.size()>0)
        {
          block_name = values[0];
          depth++;
        }
       
        values = CheckGATECommand("/gate/systems/ecat/crystal/attach", quickLine);
        if (values.size()>0)
        {
          crystal_name = values[0];
          depth++;
        }
      }
    }
    else
    {
      Cerr("***** GetGATEAliasesEcat()-> Couldn't open mac file "<< path_mac_files[f].c_str()<< " !" << endl);
      return 1;
    }
  }
    
    
  // Check we have the required number of elements for the system
  if (depth < 2)
  {
    Cerr("***** GetGATEAliasesEcat() :: Error : Missing elements in the system architecture" << endl
      << "                   The following lines are required :" << endl
      << "                   - /gate/systems/ecat/block/attach" << endl
      << "                   - /gate/systems/ecat/crystal/attach" << endl);
    return 1;
  }
  else 
  {
    if(vb >= 2)
      Cout("First container's name (usually block) is : " << block_name << endl
        << "Second container's name (usually crystal) is : " << crystal_name << endl);
  }
  
  return 0;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetGATEAliasesSPECT(vector<string> path_mac_files,
                          string&        base_name,
                          string&        crystal_name,
                          string&        pixel_name,
                           int            vb )
  \param path_mac_files
  \param base_name
  \param crystal_name
  \param pixel_name
  \param vb
  \brief Loop over a list of path to GATE macro files passed in parameter 
         to recover aliases of the different parts of an ecat system
  \return 0 if success, positive value otherwise
*/
int GetGATEAliasesSPECT(vector<string> path_mac_files,
                        string&        base_name,
                        string&        crystal_name,
                        string&        pixel_name,
                        int            vb )
{
  int depth = 0;
  
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream mac_file(path_mac_files[f].c_str(), ios::in);
        
    if(mac_file)
    {
      // Reading the .mac file line by line and finding the names of the different parts of the architecture
      string quickLine;

      while(getline(mac_file, quickLine))
      {
        vector <string> values;

        values = CheckGATECommand("/gate/systems/SPECThead/base/attach", quickLine);
        if (values.size()>0)
        {
          base_name = values[0];
          depth++;
        }
        
        values = CheckGATECommand("/gate/systems/SPECThead/crystal/attach", quickLine);
        if (values.size()>0)
        {
          crystal_name = values[0];
          depth++;
        }
       
        values = CheckGATECommand("/gate/systems/SPECThead/pixel/attach", quickLine);
        if (values.size()>0)
        {
          pixel_name = values[0];
          depth++;
        }
      }

    }
    else
    {
      Cerr("***** GetGATEAliasesSPECT()-> Couldn't open mac file "<< path_mac_files[f].c_str()<< " !" << endl);
      return 1;
    }
  }
    
  // Check we have the required number of elements for the system
  if (depth < 1)
  {
    Cerr("***** GetGATEAliasesSPECT() :: Error : Missing elements in the system architecture" << endl
      << "                              The following line is required :" << endl
      << "                              - /gate/systems/SPECThead/crystal/attach" << endl);
    return 1;
  }
  else 
  {
    if(vb >= 2)
      Cout("Crystal container's name is : " << crystal_name << endl);
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertIDecat
  \param nBlocksPerRing
  \param nBlocksLine
  \param nCrystalsTransaxial
  \param nCrystalsAxial
  \param crystalID
  \param blockID
  \brief Compute a CASToR crystal index of a GATE ecat system from its indexes (block/crystal) and the system layout
  \return CASToR crystal index
*/
uint32_t ConvertIDecat( int32_t nBlocksPerRing, 
                        int32_t nBlocksLine, 
                        int32_t nCrystalsTransaxial, 
                        int32_t nCrystalsAxial, 
                        int32_t crystalID, 
                        int32_t blockID)
{    
  int32_t nCrystalsPerRing = nBlocksPerRing * nCrystalsTransaxial;

  int32_t ringID = (int32_t)( blockID/nBlocksPerRing ) * nCrystalsAxial 
                  + (int32_t)( crystalID/nCrystalsTransaxial ); 

  int32_t castorID = nCrystalsPerRing * ringID 
                    + nCrystalsTransaxial*( blockID % nBlocksPerRing ) 
                    + crystalID % nCrystalsTransaxial;

  return castorID;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertIDSPECTRoot1
  \param a_headID : head index as provided in the root file
  \param a_rotAngle : rotation angle (deg) as provided in the root file
  \param a_angStep : angular step (deg) computed from the macro file
  \param a_nProjectionsByHead : total number of projections for each head
  \brief Compute a CASToR projection index of a GATE SPECThead system
  \return CASToR crystal index
*/
uint32_t ConvertIDSPECTRoot1( int32_t a_headID,
                              float_t a_rotAngle,
                              float_t a_angStep,
                              uint32_t a_nProjectionsByHead)
{
  // Compute angular index from the angle position of the head and the angular step
  int32_t angID = round(a_rotAngle/a_angStep);

  // Get final index for this head
  int32_t castorID = a_headID*a_nProjectionsByHead + angID;

  return castorID;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertIDSPECTRoot2
  \param a_nbSimulatedPixels
  \param a_nPixTrs
  \param a_nPixAxl
  \param a_headID
  \param a_crystalID
  \param a_pixelID
  \param a_rotAngle
  \param a_headAngPitch
  \param a_crystalSizeAxl
  \param a_crystalSizeTrs
  \param a_gPosX
  \param a_gPosY
  \param a_gPosZ
  \brief Compute a CASToR crystal index of a GATE SPECThead system
  \return CASToR crystal index. Return a higher number than the max index if error
*/
uint32_t ConvertIDSPECTRoot2( uint32_t a_nbSimulatedPixels,
                              uint32_t a_nPixTrs, 
                              uint32_t a_nPixAxl, 
                              int32_t a_headID,
                              int32_t a_crystalID,
                              int32_t a_pixelID,
                              float_t a_rotAngle,
                              float_t a_headAngPitch,
                              float_t a_crystalSizeAxl,
                              float_t a_crystalSizeTrs,
                              float_t a_gPosX,
                              float_t a_gPosY,
                              float_t a_gPosZ)
{
  int32_t castorID = 0;

  // we have a pixel matrix, just return the pixelID
  if (a_nbSimulatedPixels > 1)
  {
    castorID = a_pixelID;
  }
  // Compute the pixelID according to the global XYZ positions in the crystal
  // Cheap implementation. Does not take the DOI into account (would depend on collimator)
  else
  {
    // Compute pixel sizes
    FLTNB sizePixTrs = a_crystalSizeTrs/a_nPixTrs;
    FLTNB sizePixAxl = a_crystalSizeAxl/a_nPixAxl;

    // Compute axial ID
    uint32_t axialID = (uint32_t)(( a_gPosZ + a_nPixAxl/2*sizePixAxl) / sizePixAxl);
  
    // Compute transaxial ID
    // Compute head position angle ( in GATE referential) 
    float_t ang = a_headID*a_headAngPitch + a_rotAngle;
    
    // Get transaxial position
    float_t sin_a = sin(-ang*M_PI/180);
    float_t cos_a = cos(-ang*M_PI/180);
    float_t trs_pos = a_gPosX*sin_a + a_gPosY*cos_a ;

    // Compute transaxial ID
    uint32_t transID = (uint32_t)(( trs_pos + a_nPixTrs/2*sizePixTrs) / sizePixTrs);

    if ( axialID < a_nPixAxl && transID < a_nPixTrs )
    {
      castorID = axialID*a_nPixTrs + transID;
    }
    else // return more than the max crystal ID to check for errors
    {
      castorID = a_nPixAxl*a_nPixTrs; 
    }
    
  }

  return castorID;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ConvertIDcylindrical
  \param nRsectorsAngPos
  \param nRsectorsAxial
  \param a_rsectorIdOrder
  \param nModulesTransaxial
  \param nModulesAxial
  \param nSubmodulesTransaxial
  \param nSubmodulesAxial
  \param nCrystalsTransaxial
  \param nCrystalsAxial
  \param nLayers
  \param nb_crystal_per_layer
  \param nLayersRptTransaxial
  \param nLayersRptAxial
  \param layerID
  \param crystalID
  \param submoduleID
  \param moduleID
  \param rsectorID
  \brief Compute a CASToR crystal index of a GATE cylindricalPET system from its indexes (rsector/module/submodule/crystal) and the system layout
  \return CASToR crystal index
*/
uint32_t ConvertIDcylindrical(uint32_t  nRsectorsAngPos, 
                              uint32_t  nRsectorsAxial, 
                              int       a_rsectorIdOrder,
                              uint32_t  nModulesTransaxial,
                              uint32_t  nModulesAxial,
                              uint32_t  nSubmodulesTransaxial,
                              uint32_t  nSubmodulesAxial,
                              uint32_t  nCrystalsTransaxial, 
                              uint32_t  nCrystalsAxial, 
                              uint8_t   nLayers,
                              uint32_t* nCrystalPerLayer, 
                       vector<uint32_t> nLayersRptTransaxial,
                       vector<uint32_t> nLayersRptAxial,
                               int32_t  layerID,
                               int32_t  crystalID, 
                               int32_t  submoduleID, 
                               int32_t  moduleID, 
                               int32_t  rsectorID)
{ 
  // Castor ID definition
  uint32_t castorID = 0;
  uint8_t  layer = 0;


  if(nLayersRptTransaxial.size()==0 && nLayersRptAxial.size()==0)
  {
    // layerID represents the actual layer level
    
    // add the number of crystals contained in previous layers as
    // CASToR indexes all crystals of a layer ring before the next layer
    layer = layerID;
    
    for(int l=0 ; l<layer; l++)
      castorID += nCrystalPerLayer[l];

    int32_t nTrsCrystalsPerSubmodule = nCrystalsTransaxial;
    int32_t nTrsCrystalsPerModule = nTrsCrystalsPerSubmodule * nSubmodulesTransaxial;
    int32_t nTrsCrystalsPerRsector = nTrsCrystalsPerModule * nModulesTransaxial;
    int32_t nCrystalsPerRing = nTrsCrystalsPerRsector * nRsectorsAngPos;
    
    // Rsector axial(=ring) and transaxial(=angular) ID
    // Fastest ordering orientation (axial or transaxial) depends on the repeaters
    int32_t rsectorAxlID = 0 ; 
    int32_t rsectorTrsID = 0 ;
    
    // standard, transaxial first
    if(a_rsectorIdOrder == 0)
    {
      rsectorAxlID = rsectorID/nRsectorsAngPos ; 
      rsectorTrsID = (int32_t)(rsectorID%nRsectorsAngPos) ;
    }
    else // using cubic array, axial first
    {
      rsectorAxlID = rsectorID%nRsectorsAxial ; 
      rsectorTrsID = (int32_t)(rsectorID/nRsectorsAxial) ;
    }
    
    
    int32_t ringID = rsectorAxlID * nModulesAxial * nSubmodulesAxial * nCrystalsAxial
                   + (int32_t)(moduleID/nModulesTransaxial) * nSubmodulesAxial * nCrystalsAxial
                   + (int32_t)(submoduleID/nSubmodulesTransaxial) * nCrystalsAxial
                   + (int32_t)(crystalID/nCrystalsTransaxial);
    
    castorID += nCrystalsPerRing * ringID 
             +  nTrsCrystalsPerRsector * rsectorTrsID
             +  nTrsCrystalsPerModule * ( moduleID % nModulesTransaxial ) 
             +  nTrsCrystalsPerSubmodule * ( submoduleID % nSubmodulesTransaxial ) 
             +  crystalID % nCrystalsTransaxial;
  }
  
  else
  {
    // layerID represents a crystal layer element

    // Get the total number of crystals in the first layer
    uint32_t sum_crystals = nLayersRptTransaxial[layer] 
                          * nLayersRptAxial[layer];

    // Get the layer which the crystal belongs to
    while (layerID >= (int32_t)sum_crystals)
    {
      layer++;
      sum_crystals += nLayersRptTransaxial[layer] 
                    * nLayersRptAxial[layer];
    }

    // add the number of crystals contained in previous layers as
    // CASToR indexes all crystals of a layer ring before the next layer
    for(int l=0 ; l<layer ; l++)
      castorID += nCrystalPerLayer[l];
      
    int32_t nTrsCrystalsPerSubmodule = nCrystalsTransaxial * nLayersRptTransaxial[layer];
    int32_t nTrsCrystalsPerModule = nTrsCrystalsPerSubmodule * nSubmodulesTransaxial;
    int32_t nTrsCrystalsPerRsector = nTrsCrystalsPerModule * nModulesTransaxial;
    int32_t nCrystalsPerRing = nTrsCrystalsPerRsector * nRsectorsAngPos;
    
    // Rsector axial(=ring) and transaxial(=angular) ID
    // Fastest ordering orientation (axial or transaxial) depends on the repeaters
    int32_t rsectorAxlID = 0 ; 
    int32_t rsectorTrsID = 0 ;
    
    // standard, transaxial first
    if(a_rsectorIdOrder == 0)
    {
      rsectorAxlID = rsectorID/nRsectorsAngPos ; 
      rsectorTrsID = (int32_t)(rsectorID%nRsectorsAngPos) ;
    }
    else // using cubic array, axial first
    {
      rsectorAxlID = rsectorID%nRsectorsAxial ; 
      rsectorTrsID = (int32_t)(rsectorID/nRsectorsAxial) ;
    }
    
    
    int32_t ringID = rsectorAxlID * nModulesAxial * nSubmodulesAxial * nCrystalsAxial * nLayersRptAxial[layer]
                   + (int32_t)(moduleID/nModulesTransaxial) * nSubmodulesAxial * nCrystalsAxial * nLayersRptAxial[layer]
                   + (int32_t)(submoduleID/nSubmodulesTransaxial) * nCrystalsAxial * nLayersRptAxial[layer]
                   + (int32_t)(crystalID/nCrystalsTransaxial) * nLayersRptAxial[layer];
    
    if(nLayersRptTransaxial.empty() || nLayersRptAxial.empty() )
      ringID += (int32_t)(layerID/nLayersRptTransaxial[layer]);
    
    castorID += nCrystalsPerRing * ringID 
             +  nTrsCrystalsPerRsector * rsectorTrsID 
             +  nTrsCrystalsPerModule * ( moduleID % nModulesTransaxial ) 
             +  nTrsCrystalsPerSubmodule * ( submoduleID % nSubmodulesTransaxial ) 
             +  crystalID % nCrystalsTransaxial
             +  layerID % nLayersRptTransaxial[layer];
  }
  
  return castorID;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ComputeKindGATEEvent
  \param eventID1
  \param eventID2
  \param comptonPhantom1
  \param comptonPhantom2
  \param rayleighPhantom1
  \param rayleighPhantom2
  \brief Determine kind of a given coincidence event, from its attributes 
  \return kind of the coincidence : unknown (=0, default), true(=1), single scat(=2), multiple scat(=3), random(=4)) (default value =0)
*/
int ComputeKindGATEEvent(uint32_t eventID1, uint32_t eventID2, 
                         int comptonPhantom1, int comptonPhantom2, 
                         int rayleighPhantom1, int rayleighPhantom2)
{
  if (eventID1 != eventID2) 
    //random
    return KIND_RDM;  
  else
  {
    if (comptonPhantom1 == 0 && comptonPhantom2 == 0 &&
        rayleighPhantom1 == 0 && rayleighPhantom2 == 0)
        //true
        return KIND_TRUE;
    else
    {
      if (comptonPhantom1 == 1 || comptonPhantom2 == 1 ||
         rayleighPhantom1 == 1 || rayleighPhantom2 == 1)
          //single scat
          return KIND_SCAT;
      if (comptonPhantom1 > 1 || comptonPhantom2 > 1 ||
         rayleighPhantom1 > 1 || rayleighPhantom2 > 1)
           //multiple scat
          return KIND_MSCAT;
    }
  }
  // unknown
  return KIND_UNKNOWN; 
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadMacCylindrical
  \param a_pathMac : path to a macro file
  \param nLayers : nb of crystal layers
  \param nCrystalsAxial : nb of axial crystals (in a submodule)
  \param nCrystalsTransaxial : nb of transaxial crystals (in a submodule)
  \param nLayersRptAxial : Array containing the number of axial crystals in each layer
  \param nLayersRptTransaxial : Array containing the number of transaxial crystals in each layer
  \param nSubmodulesAxial : nb of axial submodules (in a module)
  \param nSubmodulesTransaxial : nb of transaxial submodules (in a module)
  \param nModulesAxial : nb of axial modules (in a rsector)
  \param nModulesTransaxial : nb of transaxial modules (in a rsector)
  \param nRsectorsAxial : nb of axial rsectors
  \param nRsectorsAngPos : nb of rsectors per ring
  \param rsector_id_order : ordering of rsector id
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param vb : verbosity
  \brief Recover informations about the scanner element of a cylindricalPET system and acquisition duration, from a GATE macro file
  \return 0 if success, positive value otherwise
*/
int ReadMacCylindrical(string a_pathMac,
                      uint8_t &nLayers,
                     uint32_t *nb_crystal_per_layer,
                     uint32_t &nCrystalsTot,
                     uint32_t &nCrystalsAxial,
                     uint32_t &nCrystalsTransaxial, 
             vector<uint32_t> &nLayersRptAxial,
             vector<uint32_t> &nLayersRptTransaxial,
                     uint32_t &nSubmodulesAxial, 
                     uint32_t &nSubmodulesTransaxial,
                     uint32_t &nModulesAxial, 
                     uint32_t &nModulesTransaxial,
                     uint32_t &nRsectorsAxial,
                     uint32_t &nRsectorsAngPos, 
                     int      &rsector_id_order,
                     uint32_t &start_time_ms, 
                     uint32_t &duration_ms,
                          int vb)
{
  vector<string> path_mac_files;
  path_mac_files.push_back(a_pathMac);

  // Recover path to all macro files from the main mac file
  if(GetGATEMacFiles(a_pathMac , path_mac_files))
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover paths to GATE macro files !" << endl);
    return 1;
  }

  string rsector_name = "";
  string module_name = "";
  string submodule_name = "";
  string crystal_name = "";
  string mod_rptr_type = "cubicArray";
  string smod_rptr_type = "cubicArray";
  string cry_rptr_type = "cubicArray";
  string lay_rptr_type = "cubicArray";
  
  vector <string> layers_name;
  bool is_rsector_Y_axis = false;
  
  // Recover aliases of the different parts of the architecture
  if(GetGATEAliasesCylindrical(path_mac_files, rsector_name, module_name, submodule_name, crystal_name, layers_name, vb) )
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover aliases for the elements of the cylindricalPET !" << endl);
    return 1;
  }
  
  // Recover nb of detected layers
  nLayers = layers_name.size();

  // Loop to recover all other informations
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream mac_file(path_mac_files[f].c_str(), ios::in);
    
    string line;
    double time_start =-1., 
           time_stop =-1.,
           time_slices =-1.;
    
    while(getline(mac_file, line))
    {
      vector <string> values;
      string kword ="";


      // RSECTORS

      kword = "/gate/"+rsector_name+"/placement/setTranslation";
      values = CheckGATECommand(kword, line);
      
      // Check where the first rsector has been created
      if (values.size()>0)
      {
        FLTNB rsector_pos_X =0.,
              rsector_pos_Y =0.;
              
        if(ConvertFromString(values[0], &rsector_pos_X) ||
           ConvertFromString(values[1], &rsector_pos_Y) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        // Get the axis on which the first rsector is positionned
        if(rsector_pos_Y!=0) is_rsector_Y_axis = true;
      }
      
      
      kword = "/gate/"+rsector_name+"/ring/setRepeatNumber";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nRsectorsAngPos) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }

  
      // Check if (axial) linear repeaters for rsectors
      kword = "/gate/"+rsector_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &nRsectorsAxial) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }

      
      // Check if cubic array repeaters for rsectors
      kword = "/gate/"+rsector_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(kword, line);
      
      if (values.size()>0)
      {
        // Using a cubic array after a ring repeater (nRsectorsAngPos>1) leads to rsector being ordered axially first (and not transaxially)
        // Track this using a flag
        rsector_id_order = nRsectorsAngPos>1 ? 1 : 0 ; 
        
        if(ConvertFromString(values[0], &nRsectorsAxial) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
  
  
  
      // MODULES 
      
      // Box/cubicArray repeaters
      kword = is_rsector_Y_axis ?
              "/gate/"+module_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+module_name+"/cubicArray/setRepeatNumberY";
              
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nModulesTransaxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }
  
      kword = "/gate/"+module_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nModulesAxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }

      // linear repeaters instead of box
      kword = "/gate/"+module_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nModulesAxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }


      // SUBMODULES 
      
      kword = is_rsector_Y_axis ?
              "/gate/"+submodule_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+submodule_name+"/cubicArray/setRepeatNumberY";
              
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nSubmodulesTransaxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }
  
      kword = "/gate/"+submodule_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nSubmodulesAxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }

      // linear repeaters instead of box
      kword = "/gate/"+submodule_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nSubmodulesAxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }
      


      // CRYSTALS

      kword = is_rsector_Y_axis ?
              "/gate/"+crystal_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+crystal_name+"/cubicArray/setRepeatNumberY";
              
      kword = "/gate/"+crystal_name+"/cubicArray/setRepeatNumberY";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nCrystalsTransaxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }
  
      kword = "/gate/"+crystal_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nCrystalsAxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      } 
  
  
      // linear repeaters instead of box
      kword = "/gate/"+crystal_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nCrystalsAxial) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }
  
  
      // LAYERS
        
      // Check if there are any repeaters on crystal layers
      for(int l=0 ; l<nLayers ; l++)
      {
        kword = is_rsector_Y_axis ?
              "/gate/"+layers_name[l]+"/cubicArray/setRepeatNumberX":
              "/gate/"+layers_name[l]+"/cubicArray/setRepeatNumberY";
              
        values = CheckGATECommand(kword, line);
        if (values.size()>0)
        {
          int32_t val;
          if(ConvertFromString( values[0] , &val) )
          {
            Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
            return 1;
          }
          nLayersRptTransaxial.push_back(val);
        }
        
        kword = "/gate/"+layers_name[l]+"/cubicArray/setRepeatNumberZ";
        values = CheckGATECommand(kword, line);
        if (values.size()>0)
        {
          int32_t val;
          if(ConvertFromString( values[0] , &val) )
          {
            Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
            return 1;
          }
          nLayersRptAxial.push_back(val);
        }
        
        // linear repeaters instead of box
        kword = "/gate/"+layers_name[l]+"/linear/setRepeatNumber";
        values = CheckGATECommand(kword, line);
        if (values.size()>0)
        {
          int32_t val;
          if(ConvertFromString( values[0] , &val) )
          {
            Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
            return 1;
          }
          nLayersRptAxial.push_back(val);
        }
      
      
      }
      
      
      kword = "/gate/application/setTimeStart";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &time_start) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_start *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
  
        start_time_ms = time_start;
      }
      
      kword = "/gate/application/setTimeStop";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &time_stop) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_stop *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
      }
  
      kword = "/gate/application/addSlice";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        double time_slice_tmp=0;
        
        if(ConvertFromString( values[0] , &time_slice_tmp) )
        {
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_slice_tmp *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacCylindrical()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
        
        time_slices += time_slice_tmp;
      }
    }
  
    // Compute duration in ms
    duration_ms = (time_slices>0) ? 
                  (uint32_t)(time_slices-time_start) : 
                  duration_ms;
                  
    duration_ms = (time_start>=0 && time_stop>=0)  ?
                  (uint32_t)(time_stop-time_start) :
                  duration_ms ;

    mac_file.close();
  }

  // Computing the total number of crystals in the scanner
  if(nLayers == 0)
  {
    nCrystalsTot = nRsectorsAngPos * nRsectorsAxial 
                 * nModulesTransaxial * nModulesAxial 
                 * nSubmodulesTransaxial * nSubmodulesAxial 
                 * nCrystalsTransaxial * nCrystalsAxial;
  }
  else                           
    for(int l=0 ; l<nLayers ; l++)
    {
      uint32_t nb_crystals_layer = nRsectorsAngPos * nRsectorsAxial 
                                 * nModulesTransaxial * nModulesAxial 
                                 * nSubmodulesTransaxial * nSubmodulesAxial 
                                 * nCrystalsTransaxial * nCrystalsAxial;
        
      // Add layer elements if repeaters have been used on layers
      if(nLayersRptTransaxial.size()>0 || nLayersRptAxial.size()>0    )
         nb_crystals_layer *= nLayersRptTransaxial[l] * nLayersRptAxial[l];
      
      nb_crystal_per_layer[l] = nb_crystals_layer;
      
      nCrystalsTot += nb_crystals_layer; 
    }
    

  if(vb >= 2)
  {
    Cout(endl);
    Cout("-----------------------------------------------------------" << endl);
    Cout("ReadMacCylindrical()-> Information recovered from mac file:" << endl);
    Cout("-----------------------------------------------------------" << endl);
    Cout("Number of rsectors angular position: " << nRsectorsAngPos << endl);
    Cout("Number of axial rsectors: " << nRsectorsAxial << endl);
    Cout("Number of axial modules: " << nModulesAxial << endl);
    Cout("Number of transaxial modules: " << nModulesTransaxial << endl);
    Cout("Number of axial submodules: " << nSubmodulesAxial << endl);
    Cout("Number of transaxial submodules: " << nSubmodulesTransaxial << endl);
    Cout("Number of axial crystals: " << nCrystalsAxial << endl);
    Cout("Number of transaxial crystals: " << nCrystalsTransaxial << endl);
    if (nLayers>=1)
    {
      Cout("Number of layers: " << (uint16_t)nLayers << endl);  // cast to uint16_t for output purposes
      for(int l=0 ; l<nLayers ; l++)
        Cout("Layer "<< l <<" : Number of crystals: " << nb_crystal_per_layer[l] << endl);
    }
    Cout("Total number of crystals (including layers): " << nCrystalsTot << endl);
    Cout("Acquisition start time (ms): " << start_time_ms << endl);
    Cout("Acquisition duration (ms): " << duration_ms << endl);
    Cout("-----------------------------------------------------------" << endl << endl);
  }

  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadMacECAT
  \param a_pathMac : path to a macro file
  \param nCrystalsAxial : nb of axial crystals
  \param nCrystalsTransaxial : nb of transaxial crystals
  \param nBlocksLine : nb of axial blocks
  \param nBlocksPerRing : nb of blocks per ring
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param vb : verbosity
  \brief Recover informations about the scanner element of an ECAT system, and acquisition duration, from a GATE macro file
  \return 0 if success, positive value otherwise
*/
int ReadMacECAT(  string a_pathMac,
                uint32_t &nCrystalsTot,
                uint32_t &nCrystalsAxial,
                uint32_t &nCrystalsTransaxial,
                uint32_t &nBlocksLine,
                uint32_t &nBlocksPerRing,
                uint32_t &start_time_ms,
                uint32_t &duration_ms,
                     int vb)
{
  // Recover path to all macro files from the main mac file
  vector<string> path_mac_files;
  path_mac_files.push_back(a_pathMac);
  if(GetGATEMacFiles(a_pathMac , path_mac_files))
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover paths to GATE macro files !" << endl);
    return 1;
  }
  
  string block_name = "block";
  string crystal_name = "crystal";
  bool is_block_Y_axis = false;
  
  // Recover aliases of the different parts of the architecture
  if(GetGATEAliasesEcat(path_mac_files, block_name, crystal_name, vb) )
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover aliases for the elements of the ecat !" << endl);
    return 1;
  }
  
    
  // 2nd loop to recover all other informations
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream mac_file(path_mac_files[f].c_str(), ios::in);
    
    string line;
    double time_start=-1., 
           time_stop=-1.,
           time_slices=-1.;
    
    while(getline(mac_file, line))
    {
      vector <string> values;
      string kword ="";
  
      kword = "/gate/"+block_name+"/placement/setTranslation";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        FLTNB block_pos_X=0.,
              block_pos_Y=0.;
              
        if(ConvertFromString(values[0], &block_pos_X) ||
           ConvertFromString(values[1], &block_pos_Y) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        // Get the axis on which the first rsector is positionned
        if(block_pos_Y!=0) is_block_Y_axis = true;
      }
  
      kword = "/gate/"+block_name+"/ring/setRepeatNumber";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nBlocksPerRing) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }
  
      kword = "/gate/"+block_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nBlocksLine) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }
  
      kword = is_block_Y_axis ?
            "/gate/"+crystal_name+"/cubicArray/setRepeatNumberX":
            "/gate/"+crystal_name+"/cubicArray/setRepeatNumberY";
              
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nCrystalsTransaxial) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      }
  
      kword = "/gate/"+crystal_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &nCrystalsAxial) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
      } 

      
      kword = "/gate/application/setTimeStart";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &time_start) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_start *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacECAT()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
      
        start_time_ms = time_start;
      }
      
      kword = "/gate/application/setTimeStop";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &time_stop) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_stop *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacECAT()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
      }
  
      kword = "/gate/application/addSlice";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        double time_slice_tmp=0;
        
        if(ConvertFromString( values[0] , &time_slice_tmp) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_slice_tmp *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacECAT()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
        
        time_slices = time_slice_tmp;
      }
    }

    // Compute duration in ms
    duration_ms = (time_slices>0) ? 
                  (uint32_t)(time_slices-time_start) : 
                  duration_ms;
                  
    duration_ms = (time_start>=0 && time_stop>=0)  ?
                  (uint32_t)(time_stop-time_start) :
                  duration_ms ;

    mac_file.close();
  }
  
  // Computing the total number of crystals in the scanner
  nCrystalsTot = nCrystalsTransaxial * nCrystalsAxial 
               * nBlocksLine * nBlocksPerRing;


  if(vb >= 2)
  {
    Cout(endl);
    Cout("-----------------------------------------------------" << endl);
    Cout("ReadMacECAT()-> Information recovered from mac file:" << endl);
    Cout("-----------------------------------------------------" << endl);
    Cout("Number of blocks per ring: " << nBlocksPerRing << endl);
    Cout("Number of axial blocks: " << nBlocksLine << endl);
    Cout("Number of axial crystals: " << nCrystalsAxial << endl);
    Cout("Number of transaxial crystals: " << nCrystalsTransaxial << endl);
    Cout("Total number of crystals: " << nCrystalsTot << endl);
    Cout("Acquisition start time (ms): " << start_time_ms << endl);
    Cout("Acquisition duration (ms): " << duration_ms << endl);
    Cout("-----------------------------------------------------" << endl << endl);
  }  
  
  return 0; 
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadMacSPECT
  \param a_pathMac : path to a macro file
  \param distToDetector : distance between center of rotation and detector surface
  \param nHeads : nb of cameras
  \param nPixAxl : nb of transaxial pixels
  \param nPixTrs : nb of axial pixels
  \param crystalSizeAxl : crystal axial dimension
  \param crystalSizeTrs : crystal transaxial dimension
  \param head1stAngle : head first angle
  \param headAngPitch : angular pitch between heads
  \param headRotSpeed : angle between projection
  \param headRotDirection : rotation direction of the head (CW or CCW)
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param vb : verbosity
  \brief Recover informations about the scanner element of an ECAT system, and acquisition duration, from a GATE macro file
  \return 0 if success, positive value otherwise
*/
int ReadMacSPECT( string   a_pathMac,
                  float_t  &a_distToDetector,
                  uint32_t &a_nHeads,
                  uint32_t &a_nPixAxl,
                  uint32_t &a_nPixTrs,
                  float_t  &a_crystalSizeAxl,
                  float_t  &a_crystalSizeTrs,
                  uint32_t &a_nProjectionsTot,
                  uint32_t &a_nProjectionsByHead,
                  float_t  &a_head1stAngle,
                  float_t  &a_headAngPitch,
                  float_t  &a_headAngStepDeg,
                  int      &a_headRotDirection,
                  uint32_t &a_start_time_ms,
                  uint32_t &a_duration_ms,
                       int vb)
{
  // Recover path to all macro files from the main mac file
  vector<string> path_mac_files;
  path_mac_files.push_back(a_pathMac);
  
  if(GetGATEMacFiles(a_pathMac , path_mac_files))
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover paths to GATE macro files !" << endl);
    return 1;
  }
  
  string head_name = "SPECThead";
  string crystal_name = "crystal";
  string pixel_name = "pixel";
  string head_orbit_name = "";
  bool is_head_Y_axis = false;
  
  // Recover aliases of the different parts of the architecture
  if(GetGATEAliasesSPECT(path_mac_files, head_name, crystal_name, pixel_name, vb) )
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover aliases for the elements of the ecat !" << endl);
    return 1;
  }
  
  // Init variables to recover some data
  double time_start=-1., 
         time_stop=-1.,
         time_slice=-1,
         time_slices=-1.,
         time_slice_ms=-1.,
         head_rot_speed =-1.;
  
  // 2nd loop to recover all other informations
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream mac_file(path_mac_files[f].c_str(), ios::in);
    
    string line;
    
    while(getline(mac_file, line))
    {
      vector <string> values;
      string kword ="";
      kword = "/gate/"+head_name+"/placement/setTranslation";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        FLTNB head_pos_X=0.,
              head_pos_Y=0.;
              
        if(ConvertFromString(values[0], &head_pos_X) ||
           ConvertFromString(values[1], &head_pos_Y) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      
        // Get the axis on which the first rsector is positionned
        if(head_pos_Y!=0) is_head_Y_axis = true;
        
        a_distToDetector = is_head_Y_axis ? abs(head_pos_Y) : abs(head_pos_X) ;
      }

  

      kword = "/gate/"+head_name+"/ring/setRepeatNumber";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &a_nHeads))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }

      kword = "/gate/"+head_name+"/ring/setFirstAngle";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &a_head1stAngle))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      kword = "/gate/"+head_name+"/ring/setAngularPitch";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &a_headAngPitch))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      kword = "/gate/"+head_name+"/moves/insert";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &head_orbit_name))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }

      kword = "/gate/"+head_name+"/"+head_orbit_name+"/setSpeed";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &head_rot_speed))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      kword = "/gate/"+head_name+"/"+head_orbit_name+"/setPoint2";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        int rot_direction;
        if(ConvertFromString(values[2], &rot_direction))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        a_headRotDirection = rot_direction>0 ? GEO_ROT_CW : GEO_ROT_CCW;
      }
      
      // --- Crystals ---           
      kword = is_head_Y_axis ?
        "/gate/"+crystal_name+"/geometry/setXLength":
        "/gate/"+crystal_name+"/geometry/setYLength";
              
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &a_crystalSizeTrs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }

      }
  
      kword = "/gate/"+crystal_name+"/geometry/setZLength";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &a_crystalSizeAxl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      // --- Pixels ---     
      kword = is_head_Y_axis ?
              "/gate/"+pixel_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+pixel_name+"/cubicArray/setRepeatNumberY";
              
              
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &a_nPixTrs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      kword = "/gate/"+pixel_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &a_nPixAxl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      kword = "/gate/application/setTimeStart";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &time_start) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_start *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacECAT()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
      
        a_start_time_ms = time_start;
      }

      kword = "/gate/application/setTimeSlice";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &time_slice) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_slice *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacECAT()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
          
        time_slice_ms = time_slice;
      }
      
      kword = "/gate/application/setTimeStop";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &time_stop) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_stop *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacECAT()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
      }
  
      kword = "/gate/application/addSlice";
      values = CheckGATECommand(kword, line);
      if (values.size()>0)
      {
        double time_slice_tmp=0;
        
        if(ConvertFromString( values[0] , &time_slice_tmp) )
        {
          Cerr("***** dataConversionUtilities::readMacECAT()-> Error occured while trying to get value from entry '"<< kword << "' in file " << path_mac_files[f] << "!");
          return 1;
        }
        
        // Convert in ms
        if (values.size()>1)
        {
          if(values[1] == "s") time_slice_tmp *= 1000; 
        }
        else
          Cerr("***** dataConversionUtilities::readMacECAT()-> WARNING : can't read unit of '"<< kword << ". Assuming time in seconds");
        
        time_slices = time_slice_tmp;
      }
    }

    // Compute duration in ms
    
    // Time slices were provided
    a_duration_ms = (time_slices>0) ? 
                    (uint32_t)(time_slices-time_start) : 
                    a_duration_ms;

    // Time stop/start is provided
    a_duration_ms = (time_start>=0 && time_stop>=0)  ?
                    (uint32_t)(time_stop-time_start) :
                    a_duration_ms;


    // Compute the nb of projections and total number of crystals in the system
    a_nProjectionsByHead = a_duration_ms / time_slice_ms;
    a_nProjectionsTot = a_nHeads*a_nProjectionsByHead;

    // Compute angular pitch if not provided in the mac files (=-1)
    a_headAngPitch = (a_headAngPitch<0) ?
                      360./a_nHeads  :
                      a_headAngPitch ;
                    
    a_headAngStepDeg = head_rot_speed*time_slice_ms/1000.;
    
    if(head_rot_speed<0)
    {
      Cerr("***** GetGATESystemType() -> Error couldn't find line '/gate/"+head_name+"/"+head_orbit_name+"/setSpeed' !" << endl);
      Cerr("                             This information is mandatory to compute the projection angle step." << endl);
      return 1;
    }
    
    if(time_slice_ms<0)
    {
      Cerr("***** GetGATESystemType() -> Error couldn't find line '/gate/application/setTimeSlice'  !" << endl);
      Cerr("                             This information is mandatory to compute the projection angle step." << endl);
      return 1;
    }

    if(a_duration_ms == 0)
    {
      if(time_stop <0)
      {
        Cerr("***** GetGATESystemType() -> Error couldn't compute acquisition find line '/gate/application/setTimeStop'  !" << endl);
        Cerr("                             This information is mandatory to compute the acquisition duration." << endl);
        return 1;
      }
      if(time_start <0)
      {
        Cerr("***** GetGATESystemType() -> Error couldn't compute acquisition find line '/gate/application/setTimeStart'  !" << endl);
        Cerr("                             This information is mandatory to compute the acquisition duration." << endl);
        return 1;
      }
    }
    
    mac_file.close();
  }  // end loop on .mac files
  
  if(vb >= 2)
  {
    Cout(endl);
    Cout("-----------------------------------------------------" << endl);
    Cout("ReadMacSPECT()-> Information recovered from mac file:" << endl);
    Cout("-----------------------------------------------------" << endl);
    Cout("Distance to detector: " << a_distToDetector << endl);
    Cout("Number of heads: " << a_nHeads << endl);
    Cout("Number of axial pixels: " << a_nPixAxl << endl);
    Cout("Number of transaxial pixels: " << a_nPixTrs << endl);
    Cout("Crystal axial size: " << a_crystalSizeAxl << endl);
    Cout("Crystal transaxial size: " << a_crystalSizeTrs << endl);
    Cout("Number of projections per head: " << a_nProjectionsByHead << endl);
    Cout("Total number of projections: " << a_nProjectionsTot << endl);
    Cout("Head(s) first transaxial angle: " << a_head1stAngle << endl);
    Cout("Head(s) angular pitch: " << a_headAngPitch << endl);
    Cout("Angular step between projections (deg): " << a_headAngStepDeg << endl);
    Cout("Rotation direction (0=CW, 1=CCW): " << a_headRotDirection << endl);
    Cout("Acquisition start time (ms): " << a_start_time_ms << endl);
    Cout("Acquisition duration (ms): " << a_duration_ms << endl);
    Cout("-----------------------------------------------------" << endl << endl);
  }

  return 0; 
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadIntfSPECT
  \param a_pathIntf : path to the interfile header
  \param distToDetector : distance between center of rotation and detector surface
  \param nHeads : nb of cameras
  \param nPixAxl : nb of transaxial pixels
  \param nPixTrs : nb of axial pixels
  \param crystalSizeAxl : crystal axial dimension
  \param crystalSizeTrs : crystal transaxial dimension
  \param head1stAngle : head first angle
  \param headAngPitch : angular pitch between heads
  \param headRotSpeed : angle between projection
  \param headRotDirection : rotation direction of the head (CW or CCW)
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param vb : verbosity
  \brief Recover informations about the scanner element of an ECAT system, and acquisition duration, from a GATE macro file
  \return 0 if success, positive value otherwise
*/
int ReadIntfSPECT(string      a_pathIntf,
                  float_t     &a_distToDetector,
                  uint32_t    &a_nHeads,
                  uint32_t    &a_nPixAxl,
                  uint32_t    &a_nPixTrs,
                  float_t     &a_crystalSizeAxl,
                  float_t     &a_crystalSizeTrs,
                  uint32_t    &a_nProjectionsTot,
                  uint32_t    &a_nProjectionsByHead,
                  float_t     &a_head1stAngle,
                  float_t     &a_headAngPitchDeg,
                  float_t     &a_headAngStepDeg,
                  int         &a_headRotDirection,
                  uint32_t    &a_start_time_ms,
                  uint32_t    &a_duration_ms,
                  int         vb)
{
  // Read all required information in Interfile header

  string key;

  key = "matrix size [1]";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &a_nPixTrs, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  key = "matrix size [2]";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &a_nPixAxl, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  FLTNB size_pix_trs = 1.,
        size_pix_axl = 1.;
        
  key = "scaling factor (mm/pixel) [1]";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &size_pix_trs, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  key = "scaling factor (mm/pixel) [2]";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &size_pix_axl, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  
  key = "total number of images";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &a_nProjectionsTot, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  
  key = "number of projections";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &a_nProjectionsByHead, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  key = "number of detector heads";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &a_nHeads, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  key = "study duration (sec)";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &a_duration_ms, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  // convert duration in ms
  a_duration_ms *= 1000;
  
  
  FLTNB size_crystal_X =0.,
        size_crystal_Y =0.;
  
  key = "crystal x dimension (cm)";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &size_crystal_X, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  
  key = "crystal y dimension (cm)";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &size_crystal_Y, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  
  key = "crystal z dimension (cm)";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &a_crystalSizeAxl, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  
  // Just pick the larger value to get the transaxial crystal size
  // convert in mm (values in cm in Interfile)
  if(size_crystal_X>0 && size_crystal_Y>0)
    a_crystalSizeTrs = (size_crystal_X>size_crystal_Y) ? size_crystal_X*10. : size_crystal_Y*10.;
  a_crystalSizeAxl = (a_crystalSizeAxl>0) ? a_crystalSizeAxl*10. : a_crystalSizeAxl ;
  

  FLTNB head_pos_X =-1.,
        head_pos_Y =-1.,
        head_pos_Z =-1.;
        
  key = "head x translation (cm)";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &head_pos_X, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  key = "head y translation (cm)";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &head_pos_Y, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  key = "head z translation (cm)";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &head_pos_Z, 1, KEYWORD_OPTIONAL) == 1)
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  
  // Pick the largest value to initialize distance betweeen COR and detector, converted in mm (cm by default)
  if(head_pos_X > head_pos_Y)
    a_distToDetector = head_pos_X > head_pos_Z ? head_pos_X*10 : head_pos_Z*10;
  else
    a_distToDetector = head_pos_Y > head_pos_Z ? head_pos_Y*10 : head_pos_Z*10;


  key = "direction of rotation";
  string rot_dir;
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &rot_dir, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  a_headRotDirection = (rot_dir == "CCW") ? GEO_ROT_CCW : GEO_ROT_CW ;
  
  // Could have several keys with this name (each for each head)
  // It will recover the value corresponding to the first match
  key = "start angle";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &a_head1stAngle, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  
  key = "extent of rotation";
  uint32_t extent_rotation =360;
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &extent_rotation, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  a_headAngStepDeg = (FLTNB)extent_rotation / a_nProjectionsByHead;
  
  float_t first_angle = 0;
  float_t second_angle = extent_rotation;
  
  key = "start angle";
  if(IntfKeyGetValueFromFile(a_pathIntf.c_str(), key.c_str(), &first_angle, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }

  key = "start angle";
  if(IntfKeyGetRecurringValueFromFile(a_pathIntf.c_str(), key.c_str(), &second_angle, 1, KEYWORD_MANDATORY, 1))
  {
    Cerr("***** castor-GATERootToCastor :: Error when trying to read key: '"<< key << "' in interfile : " << a_pathIntf << "!" << endl);
    Cerr("                                 Either key not found or conversion error" << endl);
    return 1;
  }
  
  a_headAngPitchDeg = second_angle - first_angle;  
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      CreateGeomWithCylindrical()
  \param   a_pathMac : string containing the path to a GATE macro file
  \param   a_pathGeom : string containing the path to a CASToR output geom file
  \brief   Read a GATE macro file containing the description of a cylindricalPET system, and convert it to a geom file
  \return  0 if success, positive value otherwise
*/
int CreateGeomWithCylindrical(string a_pathMac, string a_pathGeom)
{
  /* Declaring variables */
  
  // Scanner
  string modality;
  string scanner_name = GetFileFromPath(a_pathGeom.substr(0,a_pathGeom.find(".geom")));
  double scanner_radius = 0.;
  uint32_t number_of_elements = 0;
  string description = "PET system extracted from GATE macro: " + a_pathMac;
  
  // Rsectors
  string rsector_name = "";
  uint32_t number_of_rsectors_ang = 1; // repeater on a ring
  uint32_t number_of_rsectors_axl = 1; // linear repeater
  double rsector_step_axl = 0.;
  double rsector_gap_axl = 0.;
  double rsector_first_angle = 0.;
  double rsector_angular_span = 360.;
  double rsector_pos_X = 0.;
  double rsector_pos_Y = 0.;
  uint32_t rsector_nb_zshifts = 0;
  vector <double> vec_rsector_Z_Shift;
  bool is_rsector_Y_axis = false;
  double rsector_size_trs;
  double rsector_size_axl;
  
  // Modules
  string module_name = "";
  uint32_t number_of_modules_trs = 1;
  uint32_t number_of_modules_axl = 1;
  double module_step_trs = 0.;
  double module_step_axl = 0.;
  double module_gap_trs = 0.;
  double module_gap_axl = 0.;
  double module_size_trs;
  double module_size_axl;
  
  // Submodules
  string submodule_name = "";
  uint32_t number_of_submodules_trs = 1;
  uint32_t number_of_submodules_axl = 1;
  double submodule_step_trs = 0.;
  double submodule_step_axl = 0.;
  double submodule_gap_trs = 0.;
  double submodule_gap_axl = 0.;
  double submodule_size_trs;
  double submodule_size_axl;
  
  // Crystals 
  string crystal_name = "";
  uint32_t number_of_crystals_trs = 1;
  uint32_t number_of_crystals_axl = 1;
  double crystal_step_trs = 0.;
  double crystal_step_axl = 0.;
  double crystal_gap_trs = 0.;
  double crystal_gap_axl = 0.;
  double crystal_size_depth = 0.;
  double crystal_size_trs = 0.;
  double crystal_size_axl = 0.;
  
  
  // Layers
  int number_of_layers = 0;
  string n_layers = "1"; //default value for number_of_layers
  vector <uint32_t> number_of_lyr_elts_trs;
  vector <uint32_t> number_of_lyr_elts_axl;
  
  vector <string> layers_names;
  vector <vector <double> > layers_positions;
  vector <double> layers_size_depth;
  vector <double> layers_size_trs;
  vector <double> layers_size_axl;
  
  vector <double> layers_step_trs;
  vector <double> layers_step_axl;
  
  // Optional
  uint32_t voxels_number_trs;
  uint32_t voxels_number_axl;
  double fov_trs, 
         fov_axl;
  double mean_depth_of_interaction = -1.;
  double min_angle_diff = 0.;
  
  
  
  // If we have multiple layers, we need to enter multiple values in the .geom.
  vector <string>  vec_scanner_radius;
  
  // Rsector vectors
  vector <string>  vec_number_of_rsectors_ang;
  vector <string>  vec_number_of_rsectors_axl;
  vector <string>  vec_rsector_gap_trs;
  vector <string>  vec_rsector_gap_axl;
  vector <string>  vec_rsector_first_angle;
  
  
  // Module vectors
  vector <string>  vec_number_of_modules_trs;
  vector <string>  vec_number_of_modules_axl;
  vector <string>  vec_module_gap_trs;
  vector <string>  vec_module_gap_axl;
  
  // Submodule vectors
  vector <string>  vec_number_of_submodules_trs;
  vector <string>  vec_number_of_submodules_axl;
  vector <string>  vec_submodule_gap_trs;
  vector <string>  vec_submodule_gap_axl;
  
  // Crystal vectors
  vector <string>  vec_number_of_crystals_trs;
  vector <string>  vec_number_of_crystals_axl;
  vector <string>  vec_crystal_gap_trs;
  vector <string>  vec_crystal_gap_axl;
  
  // Optionnal
  vector <string>  vec_mean_depth_of_interaction;
  vector <string>  vec_min_angle_diff;
  
  vector<string> path_mac_files;
  path_mac_files.push_back(a_pathMac);

  // Recover path to all macro files from the main mac file
  if(GetGATEMacFiles(a_pathMac , path_mac_files))
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover paths to GATE macro files !" << endl);
    return 1;
  }
  
  // Recover aliases of the different parts of the architecture
  if(GetGATEAliasesCylindrical(path_mac_files, rsector_name, module_name, submodule_name, crystal_name, layers_names, 2) )
  {
    Cerr("***** GetGATESystemType() -> Error occured when trying to recover aliases for the elements of the cylindricalPET !" << endl);
    return 1;
  }
    
  // Recover nb of detected layers
  n_layers = layers_names.size();
  

  // Loop to recover all other informations
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream systemMac(path_mac_files[f].c_str(), ios::in);
  
    // Reading the .mac file line by line and update the variables if a matching adress is found
    string line;
    while(getline(systemMac, line))
    {
      vector <string> values;
        
      // Scanner
      modality = "PET";
      
      string entry = "";
      
      entry = "/gate/"+rsector_name+"/placement/setTranslation";
      values = CheckGATECommand(entry, line);
      
      // Assumes that first rsector located on the X or Y axis
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &rsector_pos_X) ||
           ConvertFromString(values[1], &rsector_pos_Y) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
    
        // Check first rsector cartesian coordinates is 0 on X or Y axis
        // Throw error otherwise
        if(rsector_pos_X!=0 && rsector_pos_Y !=0)
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          Cerr("                              Rsector cartesian coordinates on either the X or Y axis expected to be equal to 0 "<< endl);
          return 1;
        }
        
        // Get the axis on which the first rsector is positionned
        if(rsector_pos_Y!=0) is_rsector_Y_axis = true;
          
        scanner_radius += is_rsector_Y_axis ? abs(rsector_pos_Y) : abs(rsector_pos_X) ;
      }

      // Adjust scanner radius (from center to rsector surface) according to block size
      entry = is_rsector_Y_axis ? 
              "/gate/"+rsector_name+"/geometry/setYLength" :
              "/gate/"+rsector_name+"/geometry/setXLength";
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double rsector_size_depth = 0.;
        if(ConvertFromString(values[0], &rsector_size_depth) ) 
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }

        scanner_radius -= rsector_size_depth/2;
      }
      

      // rsector transaxial size
      entry = is_rsector_Y_axis ? 
              "/gate/"+rsector_name+"/geometry/setXLength":
              "/gate/"+rsector_name+"/geometry/setYLength";
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &rsector_size_trs) ) 
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      
      // --- Rsectors ---
      
      entry = "/gate/"+rsector_name+"/ring/setModuloNumber";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &rsector_nb_zshifts) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      entry = "/gate/"+rsector_name+"/ring/setRepeatNumber";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_rsectors_ang) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }



      // Check if (axial) linear repeaters for rsectors
      entry = "/gate/"+rsector_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_rsectors_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }




      entry = "/gate/"+rsector_name+"/linear/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[2],  &rsector_step_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      // Check if cubic array repeaters for rsectors
      entry = is_rsector_Y_axis ?
              "/gate/"+rsector_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+rsector_name+"/cubicArray/setRepeatNumberY";
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        uint32_t number_of_rsectors_trs;
        
        if(ConvertFromString(values[0], &number_of_rsectors_trs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        // Check if more than one transaxial rsector is present.
        // Not currently implemented as the geometry model can vary according to whether a ring repeater is inserted before or after the cubicarray repeater
        // Throw error in this situation
        if(number_of_rsectors_trs>1)
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Error while trying to parse line: " << line<< endl);
          Cerr("                                                             The GATE system contains more than one 'transaxial' rsector " << endl);
          Cerr("                                                             The current implementation does not support such cylindricalPET model" << endl);
          Cerr("                                                             Manual implementation of the system is required (ex: model the transaxial rsectors as modules)" << endl);
          return 1;
        }
      }
  
      
      entry = "/gate/"+rsector_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(entry, line);
      
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_rsectors_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      
      entry = "/gate/"+rsector_name+"/cubicArray/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[2], &rsector_step_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      
      // axial FOV computed from rsector z-length 
      entry = "/gate/"+rsector_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &rsector_size_axl) ) 
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
      }


      entry = "/gate/"+rsector_name+"/ring/setFirstAngle";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &rsector_first_angle) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      entry = "/gate/"+rsector_name+"/ring/setAngularSpan";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &rsector_angular_span) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      // Up to 8 zshifts in Gate
      for (int i=1; i <8; i++)
      {
        entry = "/gate/"+rsector_name+"/ring/setZShift"+toString(i);
        values = CheckGATECommand(entry, line);
        if (values.size()>0)
        {
          double zshift = 0.;
          if(ConvertFromString(values[0], &zshift) )
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
        
          vec_rsector_Z_Shift.push_back(zshift);
        }
      }

  
  
      // --- Modules ---
      
      entry = is_rsector_Y_axis ?
              "/gate/"+module_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+module_name+"/cubicArray/setRepeatNumberY";
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_modules_trs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
  

      
      entry = "/gate/"+module_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(entry, line);
      
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_modules_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      
      entry = is_rsector_Y_axis ?
              "/gate/"+module_name+"/geometry/setXLength":
              "/gate/"+module_name+"/geometry/setYLength";
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &module_size_trs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
  
      entry = "/gate/"+module_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &module_size_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      entry = "/gate/"+module_name+"/cubicArray/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        string trs_step = is_rsector_Y_axis ?
                          values[0] :
                          values[1] ;
        
        if(ConvertFromString(trs_step,  &module_step_trs) ||
           ConvertFromString(values[2], &module_step_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  


  
      // linear repeaters instead of box
      entry = "/gate/"+module_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(entry, line);
      
      if (values.size()>0)
      {
        if(ConvertFromString(values[0] , &number_of_modules_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      entry = "/gate/"+module_name+"/linear/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {        
        if(ConvertFromString(values[2], &module_step_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      // --- Submodules ---
      entry = is_rsector_Y_axis ?
              "/gate/"+submodule_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+submodule_name+"/cubicArray/setRepeatNumberY";
              
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_submodules_trs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      entry = "/gate/"+submodule_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_submodules_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  

      
      
      entry = is_rsector_Y_axis ?
              "/gate/"+submodule_name+"/geometry/setXLength":
              "/gate/"+submodule_name+"/geometry/setYLength";
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &submodule_size_trs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      entry = "/gate/"+submodule_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &submodule_size_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      entry = "/gate/"+submodule_name+"/cubicArray/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        string trs_step = is_rsector_Y_axis ?
                          values[0] :
                          values[1] ;
        
        if(ConvertFromString(trs_step,  &submodule_step_trs) ||
           ConvertFromString(values[2], &submodule_step_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
  
  
      // linear repeaters instead of box
      entry = "/gate/"+submodule_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(entry, line);
      
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &number_of_submodules_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      

      entry = "/gate/"+submodule_name+"/linear/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {        
        if(ConvertFromString(values[2], &submodule_step_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      
      
  
      // --- Crystals ---     
      entry = is_rsector_Y_axis ?
              "/gate/"+crystal_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+crystal_name+"/cubicArray/setRepeatNumberY";
              
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_crystals_trs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      entry = "/gate/"+crystal_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_crystals_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      entry = "/gate/"+crystal_name+"/geometry/setXLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double x_length;
        if(ConvertFromString(values[0], &x_length) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        if (is_rsector_Y_axis)
          crystal_size_trs = x_length;
        else
          crystal_size_depth = x_length;
      }
  
      entry = "/gate/"+crystal_name+"/geometry/setYLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double y_length;
        if(ConvertFromString(values[0], &y_length) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        if (is_rsector_Y_axis)
          crystal_size_depth = y_length;
        else
          crystal_size_trs = y_length;
          
      }
  
      entry = "/gate/"+crystal_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &crystal_size_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      entry = "/gate/"+crystal_name+"/cubicArray/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        string trs_step = is_rsector_Y_axis ?
                          values[0] :
                          values[1] ;
                          
        if(ConvertFromString(trs_step,  &crystal_step_trs) ||
           ConvertFromString(values[2], &crystal_step_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      // linear repeaters instead of box
      entry = "/gate/"+crystal_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(entry, line);
      
      if (values.size()>0)
      {
        if(ConvertFromString( values[0] , &number_of_crystals_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      entry = "/gate/"+crystal_name+"/linear/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {        
        if(ConvertFromString(values[2], &crystal_step_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      


      // --- Layers ---
      entry = "/gate/"+crystal_name+"/daughters/name";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {   
        layers_names.push_back(values[0]);     
        number_of_layers++;   
      }

      // loop on the layers to get their specific information as they are integrated in layers_names.
      for (int i=0; i < number_of_layers; i++)
      {
        // assumes that first block is positionned on the x-axis
        entry = "/gate/"+layers_names[i]+"/placement/setTranslation";
        values = CheckGATECommand(entry, line);
  
        if (values.size()>0)
        {
          vector<double> layer_pos;
          for(int d=0 ; d<3 ; d++)
          {
            double pos;
            if(ConvertFromString(values[d], &pos) )
            {
              Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
              return 1;
            }
            layer_pos.push_back(pos);
          }
          
          if(is_rsector_Y_axis)
          {
            double pos = layer_pos[1];
            layer_pos[1] = layer_pos[0];
            layer_pos[0] = pos;
          }
                          
          layers_positions.push_back(layer_pos);
        }
  
  
        // assumes that first block is positionned on the x-axis
        entry = "/gate/"+layers_names[i]+"/geometry/setXLength";
        values = CheckGATECommand(entry, line);
  
        if (values.size()>0)
        {
          double xlength;
          if(ConvertFromString(values[0], &xlength) )
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          
          if (is_rsector_Y_axis)
            layers_size_trs.push_back(xlength);
          else
            layers_size_depth.push_back(xlength);
        }
  
        // assumes that first block is positionned on the x-axis
        entry = "/gate/"+layers_names[i]+"/geometry/setYLength";
        values = CheckGATECommand(entry, line);
        if (values.size()>0)
        {
          double ylength;
          if(ConvertFromString(values[0], &ylength) )
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          if (is_rsector_Y_axis)
            layers_size_depth.push_back(ylength);
          else
            layers_size_trs.push_back(ylength);
        }
  
        entry = "/gate/"+layers_names[i]+"/geometry/setZLength";
        values = CheckGATECommand(entry, line);
        if (values.size()>0)
        {
          double zlength;
          if(ConvertFromString(values[0], &zlength) )
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          layers_size_axl.push_back(zlength);
        }
  
        // Check if a repeater if applied to the layers elements
        // In this case, the total number of crystals for a layer[i]
        // will be crystal elts*layers elts[i] 
        entry = is_rsector_Y_axis ?
                "/gate/"+layers_names[i]+"/cubicArray/setRepeatNumberX":
                "/gate/"+layers_names[i]+"/cubicArray/setRepeatNumberY";
              
        values = CheckGATECommand(entry, line);
        if (values.size()>0)
        {
          uint32_t step_trs;
          if(ConvertFromString(values[0], &step_trs))
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          number_of_lyr_elts_trs.push_back(step_trs);
        }
        
        entry = "/gate/"+layers_names[i]+"/cubicArray/setRepeatNumberZ";
        values = CheckGATECommand(entry, line);
        if (values.size()>0)
        {
          uint32_t step_axl;
          if(ConvertFromString(values[0], &step_axl))
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          number_of_lyr_elts_axl.push_back(step_axl);
        }
      
        entry = "/gate/"+layers_names[i]+"/cubicArray/setRepeatVector";
        values = CheckGATECommand(entry, line);
        if (values.size()>0)
        {
          string trs_step = is_rsector_Y_axis ?
                            values[0] :
                            values[1] ;
                          
          double step_trs, step_axl;
          if(ConvertFromString(trs_step,  &step_trs) ||
             ConvertFromString(values[2], &step_axl) )
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          
          
          layers_step_trs.push_back(step_trs);
          layers_step_axl.push_back(step_axl);
        }


        // linear repeaters instead of box
        entry = "/gate/"+layers_names[i]+"/linear/setRepeatNumber";
        values = CheckGATECommand(entry, line);
        
        if (values.size()>0)
        {
          uint32_t step_axl;
          if(ConvertFromString( values[0] , &step_axl) )
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          number_of_lyr_elts_axl.push_back(step_axl);
        }


        entry = "/gate/"+layers_names[i]+"/linear/setRepeatVector";
        values = CheckGATECommand(entry, line);
        if (values.size()>0)
        {
          double step_axl;
          if(ConvertFromString(values[2], &step_axl))
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          
          layers_step_axl.push_back(step_axl);
        }
      
      }
      
      
      entry = "/gate/digitizer/Coincidences/minSectorDifference";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        FLTNB min_sector_diff= 0.;
        
        if(ConvertFromString(values[0], &min_sector_diff))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        min_angle_diff = 360./number_of_rsectors_ang*min_sector_diff;
      }
      
    }
    systemMac.close();

  }

  // Compute total number of crystal elements
  if(number_of_layers == 0)
  {
    number_of_elements = number_of_rsectors_ang
                       * number_of_rsectors_axl 
                       * number_of_modules_trs 
                       * number_of_modules_axl 
                       * number_of_submodules_trs 
                       * number_of_submodules_axl 
                       * number_of_crystals_trs 
                       * number_of_crystals_axl;
  }
  else  
    for(int l=0 ; l<number_of_layers ; l++)
    {
      int32_t nb_crystals_layer = number_of_rsectors_ang
                                * number_of_rsectors_axl 
                                * number_of_modules_trs 
                                * number_of_modules_axl 
                                * number_of_submodules_trs 
                                * number_of_submodules_axl 
                                * number_of_crystals_trs 
                                * number_of_crystals_axl;
                                
      // Add layer elements if repeaters have been used on layers
      if(number_of_lyr_elts_trs.size()>0 || number_of_lyr_elts_axl.size()>0 )
         nb_crystals_layer *= number_of_lyr_elts_trs[l] * number_of_lyr_elts_axl[l];
      
      number_of_elements += nb_crystals_layer; 
    }

   
  // Compute sizes and gaps of each cylindrical PET structures. 
  
  // Compute crystal gaps
  if (crystal_step_axl - crystal_size_axl >= 0)
    crystal_gap_axl = crystal_step_axl - crystal_size_axl;
        
  if (crystal_step_trs - crystal_size_trs >= 0)    
    crystal_gap_trs = crystal_step_trs - crystal_size_trs;
  
  // Compute submodule size from crystals size, gaps & number
  submodule_size_axl = crystal_size_axl*number_of_crystals_axl + crystal_gap_axl*(number_of_crystals_axl-1);
  submodule_size_trs = crystal_size_trs*number_of_crystals_trs + crystal_gap_trs*(number_of_crystals_trs-1);
  
  // Compute submodule gaps
  if (submodule_step_axl - submodule_size_axl >= 0)
    submodule_gap_axl = submodule_step_axl - submodule_size_axl;
  
  if (submodule_step_trs - submodule_size_trs >= 0)
    submodule_gap_trs = submodule_step_trs - submodule_size_trs;  
    
  // Compute module size from submodule size, gaps & number
  module_size_axl = submodule_size_axl*number_of_submodules_axl + submodule_gap_axl*(number_of_submodules_axl-1);
  module_size_trs = submodule_size_trs*number_of_submodules_trs + submodule_gap_trs*(number_of_submodules_trs-1);

  // Compute module gaps
  if (module_step_axl - module_size_axl >= 0)
    module_gap_axl = module_step_axl - module_size_axl;
  
  if (module_step_trs - module_size_trs >= 0)
    module_gap_trs = module_step_trs - module_size_trs;  
    
  // Compute rsector size from crystals size, gaps & number
  rsector_size_axl = module_size_axl*number_of_modules_axl + module_gap_axl*(number_of_modules_axl-1);
  rsector_size_trs = module_size_trs*number_of_modules_trs + module_gap_trs*(number_of_modules_trs-1);


  // Compute rsector gaps. 
  if (rsector_step_axl - rsector_size_axl >= 0)
    rsector_gap_axl = rsector_step_axl - rsector_size_axl;
    

  // Axial FOV
  fov_axl = rsector_size_axl * number_of_rsectors_axl // rsectors
          + (number_of_rsectors_axl-1)*rsector_gap_axl; // gaps
  
  // Just use twice the total number of crystals to arbitrary define a number of voxels
  //voxels_number_axl = fov_axl/4 + 1;
  voxels_number_axl = ( number_of_crystals_axl
                      * number_of_submodules_axl
                      * number_of_modules_axl
                      * number_of_rsectors_axl )
                    * 2;
  
  // Transaxial FOV defined as half diameter/1.5
  fov_trs = (2*scanner_radius ) / 1.5;
  
  // Arbitrary define the transaxial voxel size as axial_voxel_size / 1.5
  //voxels_number_trs = fov_trs/4 + 1;
  voxels_number_trs = ( number_of_crystals_trs
                      * number_of_submodules_trs
                      * number_of_modules_trs
                      * number_of_rsectors_ang )
                    / 3;
  
  
  // Compute first angle
  // CASToR x and y axis for rotation inverted in comparison with GATE
  rsector_first_angle -= round(atan2f(rsector_pos_X , rsector_pos_Y) * 180. / M_PI);

  // Computing the other parts of variables, if their is more than one layer
  if (number_of_layers > 0)
  {
    for (int l=0; l < number_of_layers ; l++)
    {
      // Scanner vectors
      vec_scanner_radius.push_back( toString(scanner_radius+layers_positions[l][0]) );

      // Rsector vectors
      vec_number_of_rsectors_ang.push_back( toString(number_of_rsectors_ang) );
      vec_number_of_rsectors_axl.push_back( toString(number_of_rsectors_axl) );
      vec_rsector_gap_axl.push_back( toString(rsector_gap_axl) );
      vec_rsector_first_angle.push_back( toString(rsector_first_angle) );
      
      // Module vectors
      vec_number_of_modules_trs.push_back( toString(number_of_modules_trs) );
      vec_number_of_modules_axl.push_back( toString(number_of_modules_axl) );
      vec_module_gap_trs.push_back( toString(module_gap_trs) );
      vec_module_gap_axl.push_back( toString(module_gap_axl) );
      
      // Submodule vectors
      vec_number_of_submodules_trs.push_back( toString(number_of_submodules_trs) );
      vec_number_of_submodules_axl.push_back( toString(number_of_submodules_axl) );
      vec_submodule_gap_trs.push_back( toString(submodule_gap_trs) );
      vec_submodule_gap_axl.push_back( toString(submodule_gap_axl) );
          
      // Crystal vectors
      uint32_t nb_tot_trs_cry = (number_of_lyr_elts_trs.size()>0) ? 
                                number_of_lyr_elts_trs[l]*number_of_crystals_trs : 
                                number_of_crystals_trs ;

      uint32_t nb_tot_axl_cry = (number_of_lyr_elts_axl.size()>0) ? 
                                number_of_lyr_elts_axl[l]*number_of_crystals_axl : 
                                number_of_crystals_axl ;
                                
      vec_number_of_crystals_trs.push_back( toString(nb_tot_trs_cry) );
      vec_number_of_crystals_axl.push_back( toString(nb_tot_axl_cry) );
      

      // Compute the gap from the layer variables
      // If those variables have not been set, use the crystals gap previously computed
      if(layers_step_trs.size()>0 ||
         layers_step_axl.size()>0)
      {
        if (layers_step_trs[l] - layers_size_trs[l] >= 0)
          crystal_gap_trs = layers_step_trs[l] - layers_size_trs[l];
    
        if (layers_step_axl[l] - layers_size_axl[l] >= 0)
          crystal_gap_axl = layers_step_axl[l] - layers_size_axl[l];
      }
      
      vec_crystal_gap_trs.push_back( toString(crystal_gap_trs) );
      vec_crystal_gap_axl.push_back( toString(crystal_gap_axl) );
      
      // Optionnal
      vec_mean_depth_of_interaction.push_back( toString(mean_depth_of_interaction) );
      vec_min_angle_diff.push_back( toString(min_angle_diff) );
    }
  }
  else
  {
    // Layer dimensions = crystal dimensions
    layers_size_depth.push_back( crystal_size_depth );
    layers_size_trs.push_back( crystal_size_trs );
    layers_size_axl.push_back( crystal_size_axl );

    
    // Scanner vectors
    vec_scanner_radius.push_back( toString(scanner_radius) );
        
    // Rsector vectors
    vec_number_of_rsectors_ang.push_back( toString(number_of_rsectors_ang) );
    vec_number_of_rsectors_axl.push_back( toString(number_of_rsectors_axl) );
    vec_rsector_gap_axl.push_back( toString(rsector_gap_axl) );
    vec_rsector_first_angle.push_back( toString(rsector_first_angle) );
        
    // Module vectors
    vec_number_of_modules_trs.push_back( toString(number_of_modules_trs) );
    vec_number_of_modules_axl.push_back( toString(number_of_modules_axl) );
    vec_module_gap_trs.push_back( toString(module_gap_trs) );
    vec_module_gap_axl.push_back( toString(module_gap_axl) );
        
    // Submodule vectors
    vec_number_of_submodules_trs.push_back( toString(number_of_submodules_trs) );
    vec_number_of_submodules_axl.push_back( toString(number_of_submodules_axl) );
    vec_submodule_gap_trs.push_back( toString(submodule_gap_trs) );
    vec_submodule_gap_axl.push_back( toString(submodule_gap_axl) );
        
    // Crystal vectors
    vec_number_of_crystals_trs.push_back( toString(number_of_crystals_trs) );
    vec_number_of_crystals_axl.push_back( toString(number_of_crystals_axl) );
    vec_crystal_gap_trs.push_back( toString(crystal_gap_trs) );
    vec_crystal_gap_axl.push_back( toString(crystal_gap_axl) );
        
    // Optionnal
    vec_mean_depth_of_interaction.push_back( toString(mean_depth_of_interaction) );
    vec_min_angle_diff.push_back( toString(min_angle_diff) );  
    
    // Update the real number of layers
    number_of_layers = 1;
  }
    
    
  // Write the .geom file
  ofstream fileGeom(a_pathGeom.c_str(), ios::out | ios::trunc);  
  if(fileGeom)
  {
    fileGeom <<"# comments" << endl;
    fileGeom <<"#       Y                                        _________  "<< endl;
    fileGeom <<"#       |                                       / _ \\     \\ "<< endl;
    fileGeom <<"#       |                                      | / \\ |     |"<< endl;    
    fileGeom <<"#       |_____ Z                               | | | |     |"<< endl;        
    fileGeom <<"#        \\                                     | | | |     |" << endl;                  
    fileGeom <<"#         \\                                    | \\_/ |     |"  << endl;                   
    fileGeom <<"#          X                                    \\___/_____/"   << endl;      
    fileGeom <<"# Left-handed axis orientation"<< endl;
    fileGeom <<"# scanner axis is z" << endl;
    fileGeom <<"# positions in millimeters"<< endl;
    fileGeom <<"# Use comma without space as separator in the tables." << endl;

    fileGeom << ""<< endl;

    // Mandatory fields
    fileGeom << "# MANDATORY FIELDS"<< endl;
    fileGeom << "modality : " << modality << endl;
    fileGeom << "scanner name : " << scanner_name << endl;
    fileGeom << "number of elements              : " << number_of_elements << endl; 
    fileGeom << "number of layers : " << number_of_layers << endl;
    fileGeom << "" << endl;
    fileGeom << "voxels number transaxial        : " << voxels_number_trs << endl;
    fileGeom << "voxels number axial                : " << voxels_number_axl << endl;

    fileGeom << "field of view transaxial        : " << fov_trs << endl;
    fileGeom << "field of view axial                : " << fov_axl << endl << endl;
    fileGeom << "description        : " << description << endl;
    fileGeom << "" << endl;
    WriteVector(fileGeom,"scanner radius : ",vec_scanner_radius);
    WriteVector(fileGeom,"number of rsectors              : ",vec_number_of_rsectors_ang);
    WriteVector(fileGeom,"number of crystals transaxial    : ",vec_number_of_crystals_trs);
    WriteVector(fileGeom, "number of crystals axial            : ",vec_number_of_crystals_axl);
    fileGeom << ""<< endl;
    WriteVector(fileGeom, "crystals size depth                : ", layers_size_depth);
    WriteVector(fileGeom, "crystals size transaxial          : ", layers_size_trs);
    WriteVector(fileGeom, "crystals size axial                 : ", layers_size_axl);
    fileGeom << ""<< endl;
    fileGeom << ""<< endl;
    
    // Optional fields
    fileGeom << "# OPTIONAL FIELDS"<< endl;
    WriteVector(fileGeom,"rsectors first angle              : ",vec_rsector_first_angle);
    WriteVector(fileGeom, "number of rsectors axial            : ",vec_number_of_rsectors_axl);
    WriteVector(fileGeom,"rsector gap transaxial                : ",vec_rsector_gap_trs);
    WriteVector(fileGeom,"rsector gap axial                        : ",vec_rsector_gap_axl);
    WriteVector(fileGeom,"number of modules transaxial    : ",vec_number_of_modules_trs);
    WriteVector(fileGeom, "number of modules axial            : ",vec_number_of_modules_axl);
    WriteVector(fileGeom,"module gap transaxial                : ",vec_module_gap_trs);
    WriteVector(fileGeom,"module gap axial                        : ",vec_module_gap_axl);
    WriteVector(fileGeom,"number of submodules transaxial    : ",vec_number_of_submodules_trs);
    WriteVector(fileGeom, "number of submodules axial            : ",vec_number_of_submodules_axl);
    WriteVector(fileGeom,"submodule gap transaxial                : ",vec_submodule_gap_trs);
    WriteVector(fileGeom,"submodule gap axial                        : ",vec_submodule_gap_axl);
    WriteVector(fileGeom,"crystal gap transaxial                : ",vec_crystal_gap_trs);
    WriteVector(fileGeom,"crystal gap axial                        : ",vec_crystal_gap_axl);
    WriteVector(fileGeom, "mean depth of interaction       :  ", vec_mean_depth_of_interaction); 
    fileGeom << "rotation direction       : CCW " << endl; // default for GATE
    fileGeom << ""<< endl;

    // Write only if different from default values
    if(min_angle_diff > 0.)           fileGeom << "min angle difference        : " << min_angle_diff << endl; 
      
    // Convert angular span to the CASToR convention 
    if(rsector_angular_span >= 360.0005 ||
       rsector_angular_span <= 359.9995 ) // GATE bounds
    {
      rsector_angular_span *= (double)number_of_rsectors_ang/(double)(number_of_rsectors_ang-1);
      fileGeom << "rsectors angular span        : " << rsector_angular_span << endl;
    }
      
    if(rsector_nb_zshifts > 0)       fileGeom << "rsectors nbZShift                    :" << rsector_nb_zshifts << endl;
    if(!vec_rsector_Z_Shift.empty()) WriteVector(fileGeom, "rsectors ZShift                    : ", vec_rsector_Z_Shift);
    
    fileGeom.close();
    
    cout << "Output geom file written at :" << a_pathGeom << endl;
  }
  else
  {
    Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Couldn't open geom file for writing "<< a_pathGeom << " !" << endl);
    return 1;
  }
  
  
  
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      CreateGeomWithECAT()
  \param   a_pathMac : string containing the path to a GATE macro file
  \param   a_pathGeom : string containing the path to a CASToR output geom file
  \brief   Read a GATE macro file containing the description of an ecat system, and convert it to a geom file
  \return  0 if success, positive value otherwise
*/
int CreateGeomWithECAT(string a_pathMac, string a_pathGeom)
{
  /* Declaring variables */
  
  // Scanner
  string modality;
  string scanner_name = GetFileFromPath(a_pathGeom.substr(0,a_pathGeom.find_first_of(".geom")));
  double scanner_radius = 0.;
  uint32_t number_of_elements = 0;
  string description = "ECAT system extracted from GATE macro: " + a_pathMac;
  
  // blocks
  string block_name = "block";
  uint32_t number_of_blocks = 1;
  uint32_t number_of_blocks_trs = 1;
  uint32_t number_of_blocks_axl = 1;
  double block_step_trs = 0.;
  double block_step_axl = 0.;
  double block_gap_trs = 0.;
  double block_gap_axl = 0.;
  double block_size_Y;
  double block_size_Z;
  double block_pos_X = 0.;
  double block_pos_Y = 0.;
  double block_first_angle = 0.;
  double block_angular_span = 360.;
  uint32_t block_nb_zshifts = 0;
  vector <double> vec_block_Z_Shift;
  bool is_block_Y_axis = false;


  
  // Crystals 
  string crystal_name = "crystal";
  uint32_t number_of_crystals_trs = 1;
  uint32_t number_of_crystals_axl = 1;
  double crystal_step_trs = 0.;
  double crystal_step_axl = 0.;
  double crystal_gap_trs = 0.;
  double crystal_gap_axl = 0.;
  double crystal_size_depth = 0.;
  double crystal_size_trs = 0.;
  double crystal_size_axl = 0.;
  
  // Optional
  uint32_t voxels_number_trs;
  uint32_t voxels_number_axl;
  double fov_trs;
  double fov_axl;
  double min_angle_diff = 0.;

  vector<string> path_mac_files;
  path_mac_files.push_back(a_pathMac);

  // Recover path to all macro files from the main mac file
  if(GetGATEMacFiles(a_pathMac , path_mac_files))
  {
    Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Error occured when trying to recover paths to GATE macro files !" << endl);
    return 1;
  }
  
  
  // Recover aliases of the different parts of the architecture
  if(GetGATEAliasesEcat(path_mac_files, block_name, crystal_name, 2) )
  {
    Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Error occured when trying to recover aliases for the elements of the ecat !" << endl);
    return 1;
  }


  // Loop to recover all other informations
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream systemMac(path_mac_files[f].c_str(), ios::in);
      
    string line;
    while(getline(systemMac, line))
    {
      // Scanner
      modality = "PET";
      vector <string> values;
      string entry = "";
  
  
      entry = "/gate/"+block_name+"/placement/setTranslation";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &block_pos_X) ||
           ConvertFromString(values[1], &block_pos_Y) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        // Check first block cartesian coordinates is 0 on X or Y axis
        // Throw error otherwise
        if(block_pos_X!=0 && block_pos_Y !=0)
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          Cerr("                                                     Block cartesian coordinates on either the X or Y axis expected to be equal to 0 "<< endl);
          return 1;
        }
        
        // Get the axis on which the first rsector is positionned
        if(block_pos_Y!=0) is_block_Y_axis = true;
          
        scanner_radius += is_block_Y_axis ? abs(block_pos_Y) : abs(block_pos_X) ;
      }


      // Adjust scanner radius (from center to block surface) according to block size
      entry = is_block_Y_axis ? 
              "/gate/"+block_name+"/geometry/setYLength" :
              "/gate/"+block_name+"/geometry/setXLength";
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double block_size = 0.;
        if(ConvertFromString(values[0], &block_size) ) 
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }

        scanner_radius -= block_size/2;
      }
      
      
      
      entry = "/gate/"+block_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &fov_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        // 4mm voxels by default
        voxels_number_axl = fov_axl/4 + 1;
      }
      
      
      // blocks
      entry = "/gate/"+block_name+"/ring/setRepeatNumber";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_blocks))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
  
      entry = "/gate/"+block_name+"/ring/setFirstAngle";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &block_first_angle))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      entry = "/gate/"+block_name+"/ring/setAngularSpan";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &block_angular_span))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      entry = "/gate/"+block_name+"/ring/setModuloNumber";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &block_nb_zshifts) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      for (int i=1; i<8; i++)
      {
        entry = "/gate/"+block_name+"/ring/setZShift"+toString(i);
        values = CheckGATECommand(entry, line);
        if (values.size()>0)
        {
          double zshift;
          if(ConvertFromString(values[0], &zshift))
          {
            Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
            return 1;
          }
          vec_block_Z_Shift.push_back(zshift);
        }
      } 
      
      
      
      
      // Modules
      entry = "/gate/"+block_name+"/linear/setRepeatNumber";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_blocks_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }

  
      entry = "/gate/"+block_name+"/linear/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[2], &block_step_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      
      // Crystals 
      
      // assumes that first block is positionned on the x-axis
      entry = is_block_Y_axis ?
              "/gate/"+crystal_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+crystal_name+"/cubicArray/setRepeatNumberY";
                
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_crystals_trs))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      entry = "/gate/"+crystal_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_crystals_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
  
      entry = "/gate/"+crystal_name+"/cubicArray/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        string trs_step = is_block_Y_axis ?
                          values[0] :
                          values[1] ;
  
        if(ConvertFromString(values[1], &crystal_step_trs) ||
           ConvertFromString(values[2], &crystal_step_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
  
      // assumes that first block is positionned on the x-axis
      entry = "/gate/"+crystal_name+"/geometry/setXLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double x_length;
        if(ConvertFromString(values[0], &x_length))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        if (is_block_Y_axis)
          crystal_size_trs = x_length;
        else
          crystal_size_depth = x_length;
      }
      
      
      // assumes that first block is positionned on the x-axis
      entry = "/gate/"+crystal_name+"/geometry/setYLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double y_length;
        if(ConvertFromString(values[0], &y_length))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        if (is_block_Y_axis)
          crystal_size_depth = y_length;
        else
          crystal_size_trs = y_length;
      }
      
      entry = "/gate/"+crystal_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &crystal_size_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      
  
      
      entry = "/gate/digitizer/Coincidences/minSectorDifference";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {  
        double min_sector_diff;
        if(ConvertFromString(values[0], &min_sector_diff))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        min_angle_diff = 360/number_of_blocks*min_sector_diff;
      }

    }
    systemMac.close();
  }
  
  number_of_elements = number_of_blocks * 
                       number_of_blocks_axl *
                       number_of_crystals_trs * 
                       number_of_crystals_axl;
  


  // Transaxial FOV defined as scanner radius / 2
  fov_trs = (2*scanner_radius ) / 1.5;
  
  // 4mm voxels by default
  voxels_number_trs = ( number_of_crystals_trs
                      * number_of_blocks_trs)
                    / 3;


  if (crystal_step_axl - crystal_size_axl >= 0)
    crystal_gap_axl = crystal_step_axl - crystal_size_axl;
        
  if (crystal_step_trs - crystal_size_trs >= 0)
    crystal_gap_trs = crystal_step_trs - crystal_size_trs;


  // Compute submodule size from crystals size, gaps & number
  block_size_Z = crystal_size_axl*number_of_crystals_axl + crystal_gap_axl*(number_of_crystals_axl-1);
  block_size_Y = crystal_size_trs*number_of_crystals_trs + crystal_gap_trs*(number_of_crystals_trs-1);
  
  // Compute submodule gaps
  if (block_step_axl - block_size_Z >= 0)
    block_gap_axl = block_step_axl - block_size_Z;
  
  if (block_step_trs - block_size_Y >= 0)
    block_gap_trs = block_step_trs - block_size_Y;  

  // CASToR x and y axis for rotation inverted in comparison with GATE
  block_first_angle -= round(atan2f(block_pos_X , block_pos_Y) * 180. / M_PI);

    
  // Writing the .geom file
  ofstream fileGeom(a_pathGeom.c_str(), ios::out | ios::trunc);  

  if(fileGeom)
  {
    fileGeom <<"# comments" << endl;
    fileGeom <<"#       Y                                        _________  "<< endl;
    fileGeom <<"#       |                                       / _ \\     \\ "<< endl;
    fileGeom <<"#       |                                      | / \\ |     |"<< endl;
    fileGeom <<"#       |_____ Z                               | | | |     |"<< endl;
    fileGeom <<"#        \\                                     | | | |     |" << endl;
    fileGeom <<"#         \\                                    | \\_/ |     |"  << endl;
    fileGeom <<"#          X                                    \\___/_____/"   << endl;
    fileGeom <<"# Left-handed axis orientation"<< endl;
    fileGeom <<"# scanner axis is z" << endl;
    fileGeom <<"# positions in millimeters"<< endl;
    fileGeom <<"# Use comma without space as separator in the tables." << endl;


    fileGeom << "modality : " << modality << endl;
    fileGeom << "scanner name : " << scanner_name << endl;
    fileGeom << "number of elements              : " << number_of_elements << endl; 
    fileGeom << "number of layers : " << "1" << endl;
    fileGeom << "" << endl;
    fileGeom << "# default reconstruction parameters" << endl;
    fileGeom << "voxels number transaxial        : " << voxels_number_trs << endl;
    fileGeom << "voxels number axial                : " << voxels_number_axl << endl;

    fileGeom << "field of view transaxial        : " << fov_trs << endl;
    fileGeom << "field of view axial                : " << fov_axl << endl;
    fileGeom << "" << endl;
    fileGeom << "description        : " << description << endl;
    fileGeom << "" << endl;
    fileGeom << "scanner radius : " << scanner_radius << endl;
    fileGeom << "number of rsectors              : " << number_of_blocks << endl;
    fileGeom << "number of crystals transaxial    : " << number_of_crystals_trs << endl;
    fileGeom << "number of crystals axial            : " << number_of_crystals_axl << endl;

    fileGeom << ""<< endl;
    fileGeom << "# The 4 following parameters could be defined in arrays (SizeLayer1,SizeLayer2,SizeLayer3,etc..) if their is more than one layer"<< endl;
    fileGeom << "crystals size depth                : " << crystal_size_depth << endl;
    fileGeom << "crystals size transaxial          : " << crystal_size_trs << endl;
    fileGeom << "crystals size axial                 : " << crystal_size_axl << endl;
    fileGeom << ""<< endl;
    
    // Optional fields
    fileGeom << "rsectors first angle              : " << block_first_angle << endl;
    fileGeom << "number of modules transaxial    : " << number_of_blocks_trs << endl;
    fileGeom << "number of modules axial            : " << number_of_blocks_axl << endl;
    fileGeom << "module gap transaxial                : " << block_gap_trs << endl;
    fileGeom << "module gap axial                        : " << block_gap_axl << endl;
    fileGeom << "crystal gap transaxial                : " << crystal_gap_trs << endl;
    fileGeom << "crystal gap axial                        : " << crystal_gap_axl << endl;
    fileGeom << "rotation direction       : CCW " << endl; // default for GATE
    fileGeom << ""<< endl;

    // Write only if different from default values
    if(min_angle_diff > 0.)           fileGeom << "min angle difference        : " << min_angle_diff << endl; 

    // Convert angular span to the CASToR convention 
    if(block_angular_span >= 360.0005 ||
       block_angular_span <= 359.9995 ) // GATE bounds
    {
      block_angular_span *= (double)(number_of_blocks)/(double)(number_of_blocks-1);
      fileGeom << "rsectors angular span        : " << block_angular_span << endl;
    }

    if(block_nb_zshifts > 0)       fileGeom << "rsectors nbZShift                    :" << block_nb_zshifts << endl;
    if(!vec_block_Z_Shift.empty()) WriteVector(fileGeom, "rsectors ZShift                    : ", vec_block_Z_Shift);

    fileGeom.close();
    
    Cout("Output geom file written at :" << a_pathGeom << endl);
  }
  else
  {
    Cerr("***** dataConversionUtilities::CreateGeomWithECAT()-> Couldn't open geom file for writing "<< a_pathGeom << " !" << endl);
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      CreateGeomWithSPECT()
  \param   a_pathMac : string containing the path to a GATE macro file
  \param   a_pathGeom : string containing the path to a CASToR output geom file
  \brief   Read a GATE macro file containing the description of a SPECThead system, and convert it to a geom file
  \return  0 if success, positive value otherwise
*/
int CreateGeomWithSPECT(string a_pathMac, string a_pathGeom)
{
  /* Declaring variables */
  string modality;
  string scanner_name = GetFileFromPath(a_pathGeom.substr(0,a_pathGeom.find_first_of(".geom")));
  FLTNB scanner_radius = -1.;
  string description = "SPECT camera extracted from GATE macro: " + a_pathMac;
  
  // head(s)
  string head_name = "SPECThead"; // not used. Correspond to SPECThead
  uint32_t number_of_heads = 0;
  FLTNB head_pos_X = 0.;
  FLTNB head_pos_Y = 0.;
  FLTNB head_first_angle = 0.;
  FLTNB head_angular_pitch = -1.;
  string head_orbit_name = "";
  FLTNB head_rotation_speed = 0.;
  bool is_head_Y_axis = false;
  
  // crystal
  string crystal_name = "crystal";
  FLTNB crystal_size_trs = 0.;
  FLTNB crystal_size_axl = 0.;
  FLTNB crystal_depth = 0;
  
  // pixel
  string pixel_name = "pixel";
  uint32_t number_of_pixels_trs = 1;
  uint32_t number_of_pixels_axl = 1;
  FLTNB pix_size_trs = 0.;
  FLTNB pix_size_axl = 0.;
  FLTNB pix_step_trs = 0.;
  FLTNB pix_step_axl = 0.;
  FLTNB pix_gap_trs = 0.;
  FLTNB pix_gap_axl = 0.;

  // collimator parameters
  string focal_model_trs = "constant";
  uint16_t nb_coeff_model_trs = 1;
  FLTNB coeff_model_trs = 0.;
  string focal_model_axl = "constant";
  uint16_t nb_coeff_model_axl = 1;
  FLTNB coeff_model_axl = 0.;

  // others
  uint32_t voxels_number_trs;
  uint32_t voxels_number_axl;
  double fov_trs;
  double fov_axl;
  
  vector<string> path_mac_files;
  path_mac_files.push_back(a_pathMac);

  // Recover path to all macro files from the main mac file
  if(GetGATEMacFiles(a_pathMac , path_mac_files))
  {
    Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Error occured when trying to recover paths to GATE macro files !" << endl);
    return 1;
  }
  
  
  // Recover aliases of the different parts of the architecture
  if(GetGATEAliasesSPECT(path_mac_files, head_name, crystal_name, pixel_name, 2) )
  {
    Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Error occured when trying to recover aliases for the elements of the SPECThead !" << endl);
    return 1;
  }

  // Loop to recover all other informations
  for(uint16_t f=0 ; f<path_mac_files.size() ; f++)
  {
    ifstream systemMac(path_mac_files[f].c_str(), ios::in);
      
    string line;
    while(getline(systemMac, line))
    {
      // Scanner
      modality = "SPECT_CONVERGENT";
      vector <string> values;
      string entry = "";
  
      // assumes that first block is positionned on the x-axis
      entry = "/gate/"+head_name+"/placement/setTranslation";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &head_pos_X) ||
           ConvertFromString(values[1], &head_pos_Y) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        // Check first block cartesian coordinates is 0 on X or Y axis
        // Throw error otherwise
        if(head_pos_X!=0 && head_pos_Y !=0)
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          Cerr("                                                     Block cartesian coordinates on either the X or Y axis expected to be equal to 0 "<< endl);
          return 1;
        }
        
        // Get the axis on which the first rsector is positionned
        if(head_pos_Y!=0) is_head_Y_axis = true;
          
        scanner_radius = is_head_Y_axis ? abs(head_pos_Y) : abs(head_pos_X) ;
      }
  
  
      entry = "/gate/"+head_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &fov_axl))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        // 4mm voxels by default
        voxels_number_axl = fov_axl/4 + 1;
      }
      
      
      entry = "/gate/"+head_name+"/ring/setRepeatNumber";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_heads))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
  
      entry = "/gate/"+head_name+"/ring/setFirstAngle";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &head_first_angle))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      entry = "/gate/"+head_name+"/ring/setAngularPitch";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &head_angular_pitch))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
  
      entry = "/gate/"+head_name+"/moves/insert";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &head_orbit_name))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }

      entry = "/gate/"+head_orbit_name+"/setSpeed";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &head_rotation_speed))
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      

      
      // --- Crystals ---           
      entry = "/gate/"+crystal_name+"/geometry/setXLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double x_length;
        if(ConvertFromString(values[0], &x_length) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        if (is_head_Y_axis)
          crystal_size_trs = x_length;
        else
          crystal_depth = x_length;
      }
  
      entry = "/gate/"+crystal_name+"/geometry/setYLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double y_length;
        if(ConvertFromString(values[0], &y_length) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        if (is_head_Y_axis)
          crystal_depth = y_length;
        else
          crystal_size_trs = y_length;
          
      }
  
      entry = "/gate/"+crystal_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &crystal_size_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      
      
  
      // --- Pixels ---     
      entry = is_head_Y_axis ?
              "/gate/"+pixel_name+"/cubicArray/setRepeatNumberX":
              "/gate/"+pixel_name+"/cubicArray/setRepeatNumberY";
              
              
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_pixels_trs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      entry = "/gate/"+pixel_name+"/cubicArray/setRepeatNumberZ";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &number_of_pixels_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      entry = "/gate/"+pixel_name+"/geometry/setXLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double x_length;
        if(ConvertFromString(values[0], &x_length) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        if (is_head_Y_axis)
          pix_size_trs = x_length;
      }
  
      entry = "/gate/"+pixel_name+"/geometry/setYLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        double y_length;
        if(ConvertFromString(values[0], &y_length) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        if (!is_head_Y_axis)
          pix_size_trs = y_length;
          
      }
  
      entry = "/gate/"+pixel_name+"/geometry/setZLength";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0], &pix_size_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }
      
      
      entry = "/gate/"+pixel_name+"/cubicArray/setRepeatVector";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        string trs_step = is_head_Y_axis ?
                          values[0] :
                          values[1] ;
                          
        if(ConvertFromString(trs_step,  &pix_step_trs) ||
           ConvertFromString(values[2], &pix_step_axl) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
      }


      // collimator parameter
      entry = is_head_Y_axis ?
              "/gate/fanbeam/geometry/setFocalDistanceY":
              "/gate/fanbeam/geometry/setFocalDistanceX";
              
      entry = "/gate/fanbeam/geometry/setFocalDistanceX";
      values = CheckGATECommand(entry, line);
      if (values.size()>0)
      {
        if(ConvertFromString(values[0],  &coeff_model_trs) )
        {
          Cerr("***** dataConversionUtilities::CreateGeomWithCylindrical()-> Conversion error occured while trying to parse line: " << line<< endl);
          return 1;
        }
        
        focal_model_trs = "polynomial";
      }
      
    }
    systemMac.close();
  }
  
  // Transaxial FOV defined as scanner radius / 2
  fov_trs = scanner_radius/2;
  // 4mm voxels by default
  voxels_number_trs = fov_trs/2 + 1;
  
  uint32_t nb_pixels = number_of_pixels_axl * number_of_pixels_trs;
  
  pix_size_axl = nb_pixels>1 ? pix_size_axl : crystal_size_axl;
  pix_size_trs = nb_pixels>1 ? pix_size_trs : crystal_size_trs;
  
  // Compute gaps
  if (pix_step_axl - pix_size_axl >= 0)
      pix_gap_axl = pix_step_axl - pix_size_axl;
  
  if (pix_step_trs - pix_size_trs >= 0)
      pix_gap_trs = pix_step_trs - pix_size_trs;


  // CASToR x and y axis for rotation inverted in comparison with GATE
  head_first_angle = round(atan2f(head_pos_X , head_pos_Y) * 180. / M_PI)
                   - head_first_angle; 


  // Writing the .geom file
  ofstream fileGeom(a_pathGeom.c_str(), ios::out | ios::trunc);  

  if(fileGeom)
  {
    fileGeom << "modality : " << modality << endl;
    fileGeom << "scanner name : " << scanner_name << endl;
    fileGeom << "number of detector heads : " << number_of_heads << endl; 
    fileGeom << "trans number of pixels : " << number_of_pixels_trs << endl;
    fileGeom << "trans pixel size : " << pix_size_trs << endl;
    fileGeom << "trans gap size : " << pix_gap_trs << endl;
    fileGeom << "axial number of pixels : " << number_of_pixels_axl << endl;
    fileGeom << "axial pixel size : " << pix_size_axl << endl;
    fileGeom << "axial gap size : " << pix_gap_axl << endl;
    
    fileGeom << "detector depth : " << crystal_depth << endl;
    
    fileGeom << "scanner radius : " << scanner_radius;
    for(size_t h=1 ; h<number_of_heads ; h++)
      fileGeom << "," << scanner_radius;
    fileGeom <<  endl;
    
    fileGeom << "# Collimator configuration : "<< endl << endl;
    for(size_t h=0 ; h<number_of_heads ; h++)
    {
      fileGeom << "head" << h+1 << ":" << endl; 
      fileGeom << "trans focal model: " << focal_model_trs << endl;
      fileGeom << "trans number of coef model: " << nb_coeff_model_trs << endl;
      fileGeom << "trans parameters: " << coeff_model_trs << endl;
      fileGeom << "axial focal model: " << focal_model_axl << endl;
      fileGeom << "axial number of coef model: " << nb_coeff_model_axl << endl;
      fileGeom << "axial parameters: " << coeff_model_axl << endl;
      fileGeom << endl;
    }
    
    fileGeom << "" << endl;
    fileGeom << "# default reconstruction parameters" << endl;
    fileGeom << "voxels number transaxial        : " << voxels_number_trs << endl;
    fileGeom << "voxels number axial                : " << voxels_number_axl << endl;

    fileGeom << "field of view transaxial        : " << fov_trs << endl;
    fileGeom << "field of view axial                : " << fov_axl << endl << endl ;
    fileGeom << ""<< endl;
    
    fileGeom << "# description" << endl;
    fileGeom << "description        : " << description << endl;

    fileGeom.close();
    
    Cout("Output geom file written at :" << a_pathGeom << endl);
  }
  else
  {
    Cerr("***** dataConversionUtilities::CreateGeomWithSPECT()-> Couldn't open geom file for writing "<< a_pathGeom << " !" << endl);
    return 1;
  }
  
  return 0;
}
