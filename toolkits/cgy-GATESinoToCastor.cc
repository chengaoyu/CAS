
/*!
  \file
  \ingroup utils
  \brief This program convert a GATE datafile in sinogram format to a hist-mode datafile in CASToR format
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "iDataFilePET.hh"
#include "iDataFileSPECT.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include "gDataConversionUtilities.hh"
#include <iomanip> //std::setprecision

#include <random> //std::shuffle
#include <chrono> // std::chrono::system_clock


#ifdef CASTOR_ROOT
  #include "TROOT.h"
  #include "TApplication.h"
  #include "TGClient.h"
  #include "TCanvas.h"
  #include "TSystem.h"
  #include "TTree.h"
  #include "TBranch.h"
  #include "TFile.h"
#endif


#include "oImageSpace.hh"
#include "oProjectorManager.hh"

struct LOR{
    int crystID1;
    int crystID2;
};


/*!
  \fn      ShowHelp()
  \param a_returnCode
  \brief   Display main command line options for cgy-GATESinoToCastor
*/
void ShowHelp(int a_returnCode)
{
  cout << endl;
  cout << "Usage:  cgy-GATESinoToCastor      -i   path_to_ifile.ima -OR- -il path_to_ifile.txt "<< endl;
  cout << "                                     -s   scanner_alias "<< endl;
  cout << "                                     -o   path_to_out_file "<< endl;
  cout << "                                     -m   path_to_macro.mac " << endl;
  cout << "                                    [-t   only convert true prompts]" << endl;
  cout << "                                    [-oh  histogram datafile output]" << endl;
  cout << "                                    [-atn path_to_atn_image(cm-1)]" << endl;
  cout << "                                    [-k   recover coincidence kind]" << endl;
  cout << "                                    [-ist isotope_alias" << endl;
  cout << "                                    [-ot  time offsets in s]" << endl;
  cout << "                                    [-vb   verbosity]" << endl;
  cout << endl;
  cout << endl;
  cout << "[Main settings]:" << endl;
  cout << "  -i  path_to_file.ima   : give an input sinogram datafile to convert" << endl;
  cout << "  -il path_to_file.txt   : give an input text file containing a list of sinogram files to convert" << endl;
  cout << "                         : they will be concatenated into 1 CASToR datafile" << endl;
  cout << "                         : the path to each datafile to convert must be entered on a newline" << endl;
  cout << "  -m  path_to_macro.mac  : gives the input GATE macro file used for the GATE simulation" << endl;
  cout << "  -o  path_to_out_file   : give the path to the output file which will be created inside this folder (no default)" << endl;
  cout << "  -s  scanner_alias      : provide the name of the scanner used for to acquire the original data" << endl;
  cout << "                         : Must correspond to a .geom or .hscan file in the config/scanner repository." << endl;
  cout << "                         : A geom file can be created from the macro files using the facultative option -geo (see below)" << endl;
  cout << endl;
  cout << "[Optional settings]:" << endl;
  cout << "  -t                     : only the true prompts will be converted" << endl;
  cout << "  -oh                    : Indicate if the output datafile should be written in histogram format (default : List-mode)" << endl;
  cout << "  -atn path_image:proj   : For PET histogram output, provide an attenuation image (cm-1) related to the acquisition." << endl;
  cout << "                           Analytic projection will be performed during the data conversion in order to estimate PET attenuation correction factors (acf) for each histogram event" << endl;
  cout << "                           path_image : path to the attenuation image" << endl;
  cout << "                           proj       : (optional) projector to use for analytic projection. Defaut projector = Incremental siddon" << endl;
  cout << "  -k                     : For List-mode output, write kind of coincidence (true/scatter/rdm/...) in the output datafile (disabled by default)" << endl;
  cout << "  -ist isotope_alias     : provide alias of the isotope used during the simulation"<< endl;
  cout << "                           Aliases for supported PET isotopes and their parameters are listed in config/misc/isotopes_pet"<< endl;
  cout << "                           Other isotopes can be added in the same file"<< endl;
  cout << "  -isrc path_to_img:dims : Provide name and dimensions (separated by a colon) of an image generated with the source (annihilation) XYZ position"<< endl;
  cout << "                           The option must start with the path to the output image which will be generated." << endl;
  cout << "                           Dimensions and voxel sizes of the image must be provided using commas as separators, as in the following template:" << endl;
  cout << "                           path/to/image:dimX,dimY,dimZ,voxSizeX,voxSizeY,voxSizeZ"<< endl;
  cout << "  -geo                   : Generate a CASToR geom file from the provided GATE macro file(s)"<< endl;
  cout << "                           A geom file with the 'scanner_alias' (from the -s option) basename will be generated in the scanner repository (default location : /config/scanner)" << endl;
  cout << "  -sp_bins nbinsT,nbinsA : Option specific to simulation using SPECThead systems, with root output."<< endl;
  cout << "                           Give transaxial and axial number of bins for projections, separated by a comma."<< endl;
  cout << "                           Pixel sizes will be computed from the crystal surface and the transaxial/axial number of bins" << endl;
  cout << "  -ot list               : Provide a serie of time offset in seconds to apply before each input file"<< endl;
  cout << "                           (In the case of converting several datafiles of a dynamic acquisition, timestamps of events may be resetted for each file" << endl;
  cout << "                           This variable allows to manually increment the time between each datafile(s) if required" << endl;
  cout << "                           The number of time offsets must be equal to the number of input files, provided by -i or -if options." << endl;
  cout << "                          'list' is a list of time offsets, separated by ','" << endl;
  cout << endl;
  cout << "[Miscellaneous settings]:" << endl;
  cout << "  -vb             : give the verbosity level, from 0 (no verbose) to above 5 (at the event level) (default: 1)." << endl;
  cout << endl;
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
  Exit(a_returnCode);
}


/*!
 * \brief SinoIDToRingID
 * \param sinoID
 * \param num_Ring
 * \param ring1ID
 * \param ring2ID
 * get the ID of two rings which form the input sinoID in GATE ECAT system.
 */
int SinoIDToRingID(int sinoID, int num_Ring ,int& ringID1,int& ringID2){
    if(sinoID<num_Ring){ // the ring difference is zero.
        ringID1 = sinoID;
        ringID2 = sinoID;
        return 1;
    }
    else{ //the ring difference is bigger than zero.
        sinoID -= num_Ring;
        int ringdiff = 0;
        while( sinoID >= 0 ){
            ringdiff++;
            sinoID -= 2*(num_Ring-ringdiff);
        }
        if(ringdiff>=num_Ring){ //the ring difference is out of range.
            return -1;
        }
        sinoID += 2*(num_Ring-ringdiff);
        if(sinoID >= num_Ring-ringdiff){ //negative ringdiff ( ringID1 > ringID2)
            ringID2 = sinoID - (num_Ring-ringdiff);
            ringID1 = ringID2 + ringdiff;
        }
        else{ // positive ringdiff (ringID2 > ringID1)
            ringID1 = sinoID;
            ringID2 = ringdiff+ringID1;
        }
        return 1;
    }
}

/*!
 * \brief BinIDToCrystID
 * \param binID
 * \param crystID1
 * \param crystID2
 * get the crystIDs of LOR which form the input binID of a 2D sinogram in GATE ECAT system.
 */
void BinIDToCrystID(int binID, int num_RadialElem, int num_CrystPerRing, int& crystID1, int& crystID2){
    int binViewID = binID/num_RadialElem;
    int binElemID = binID%num_RadialElem;
    int sigma = binElemID + num_CrystPerRing/2 - num_RadialElem/2;
    if(sigma >= num_RadialElem){
        sigma -= num_CrystPerRing;
        sigma =std::abs(sigma);
        crystID1 = binViewID + (num_RadialElem-sigma)/2;
        crystID2 = crystID1 + sigma;
    }
    else{
        crystID1 = binViewID - (num_RadialElem-sigma)/2;
        if(crystID1 < 0){
            crystID1 += num_RadialElem;
        }
        crystID2 = crystID1 - sigma;
        if (crystID2 < 0){
            crystID2 += num_CrystPerRing;
        }
    }
}


void ConvertToGATEEcatID(int num_BlockPerRing, int num_CrystPerBlockY,int& ringID1, int& ringID2,  int& crystID1, int& crystID2){
    int num_CrystPerRing = num_BlockPerRing*num_CrystPerBlockY;
    // swap the ringID to keep pace with GATE.
    int tempRingID = ringID1;
    ringID1 = ringID2;
    ringID2 = tempRingID;

    //
    crystID1 -= num_CrystPerRing/4;
    crystID2 -= num_CrystPerRing/4;
    if (crystID1 < 0)  crystID1 += num_CrystPerRing;
    if (crystID2 <0)   crystID2 += num_CrystPerRing;
    //
    crystID1 += num_CrystPerBlockY/2;
    crystID2 += num_CrystPerBlockY/2;
    if(crystID1 >= num_CrystPerRing) crystID1 -= num_CrystPerRing;
    if(crystID2 >= num_CrystPerRing) crystID2 -= num_CrystPerRing;

}

/*
  Main program
*/

int main(int argc, char** argv)
{
  //debug

  // No argument, then show help
  if (argc==1) ShowHelp(0);

  // ============================================================================================================
  // Parameterized variables with their default values
  // ============================================================================================================

  // String which gathers the path to the input data filename provided by the user. no default
  string input_file = "";
  vector<string> path_to_input_file;

  // Path to user provided data and output files/images
  string path_to_out_file = "";
  string path_to_data_filename = "";
  string path_to_header_filename = "";
  string path_to_mac_file = "";
  string path_to_atn_image = "";
  string path_to_src_image = "";

  // Scanner name provided by the user
  string scanner_name = "";

  // GATE system type
  int GATE_system_type = GATE_SYS_UNKNOWN;

  // Only recover true GEvents
  bool true_only_flag = false;

  // Verbosity
  int vb = 0;

  // Histogram output
  bool histo_out_flag = false;

  // Estimate acf (histogram output)
  bool estimate_acf_flag = false;
  // Projector for analytic projection
  string options_projector  = "incrementalSiddon";

  // Recover kind of coincidence (list-mode output)
  bool kind_flag = false;

  // Isotope alias
  string istp_alias = "";

  // Variables related to the source image
  bool src_img_flag = false;
  INTNB dim_src_img[3];
  FLTNB vox_src_img[3];
  FLTNB* p_src_img = NULL;

  // Generate geom file
  bool geom_flag = false;

  // Input is interfile (for SPECT simulation) -> not supported (def = false)
  bool input_is_intf = false;

  // SPECT bins
  uint16_t spect_bin_axl = 0,
           spect_bin_trs = 0;

  // Time offsets
  string offset_time_str = "";
  uint32_t* offset_time_list = NULL;

  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================
  for (int i=1; i<argc; i++)
  {
    string option = (string)argv[i];

    if (option=="-h" || option=="--help" || option=="-help") ShowHelp(0);

    // Just one file is provided
    if (option=="-i")
    {
      if (path_to_input_file.size() > 0)
      {
        Cerr("***** cgy-GATESinoToCastor :: the following file names have already been provided (-i/-il options must be used ONCE): " << endl);
        for (size_t i=0 ; i<path_to_input_file.size() ; i++)
          Cout(path_to_input_file[i] << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        if (argv[i+1] == NULL)
        {
          Cerr("***** cgy-GATESinoToCastor :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
        {
          input_file = argv[i+1];
          path_to_input_file.push_back(input_file);
        }
        i++;
      }
    }

    // Read a text file containing a list of root datafiles to read
    else if (option=="-il")
    {
      if (path_to_input_file.size() > 0)
      {
        Cerr("***** cgy-GATESinoToCastor :: the following file names have already been provided (-i/-il options must be used ONCE) " << endl);
        for (size_t i=0 ; i<path_to_input_file.size() ; i++)
          Cout(path_to_input_file[i] << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        if (argv[i+1] == NULL)
        {
          Cerr("***** cgy-GATESinoToCastor :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
        {
          ifstream ifile(argv[i+1], ios::in);
          string line;

          // Read list of rootfgiles
          if(ifile)
          {
            // Get the position of the list_file, then append the name of the content datafiles to this position.
            input_file = argv[i+1];

            // Recover the files
            while (getline(ifile, line))
            {
              string path_to_datafile = GetPathOfFile(input_file) + OS_SEP + line;
              path_to_input_file.push_back(path_to_datafile);
            }
          }
          else
          {
            Cerr("***** cgy-GATESinoToCastor :: Error, can't read txt file: " << argv[i+1] << endl);
            Exit(EXIT_FAILURE);
          }

          ifile.close();
        }
        i++;
      }
    }
    // Mac file
    else if (option=="-m")
    {
      path_to_mac_file = (string)argv[i+1];
      i++;
    }
    // Scanner alias
    else if (option=="-s")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** cgy-GATESinoToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        scanner_name = argv[i+1];
      i++;
    }
    // Output CASToR datafile
    else if (option=="-o")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** cgy-GATESinoToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        path_to_out_file = argv[i+1];
      i++;
    }
    // Recover only trues
    else if (option=="-t")
    {
      #ifdef CASTOR_ROOT
        true_only_flag = true;
      #else
        Cerr("***** cgy-GATESinoToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** cgy-GATESinoToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }

    // Time offsets
    else if (option=="-ot")
    {
      if (i>=argc-1)
      {
        Cerr("***** cgy-GATESinoToCastor :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      offset_time_str = (string)argv[i+1];
      i++;
    }

    // Output datafile in histogram mode
    else if (option=="-oh")
    {
      histo_out_flag = true;
    }

    else if (option=="-atn")
    {
      if (i>=argc-1)
      {
        Cerr("***** cgy-GATESinoToCastor :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }

      path_to_atn_image = argv[i+1];

      // check first if a projector has been provided (':' included)
      size_t pos = path_to_atn_image.find_first_of(":");

      if (pos != string::npos)
      {
        options_projector = path_to_atn_image.substr(pos+1);
        path_to_atn_image = path_to_atn_image.substr(0,pos);
      }

      estimate_acf_flag = true;
      i++;
    }

    // Output datafile in histogram mode
    else if (option=="-k")
    {
      kind_flag = true;
    }

    // Provide an isotope alias
    else if (option=="-ist")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** cgy-GATESinoToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        istp_alias = argv[i+1];

      i++;
    }

    // Generate image from the sourceID
    else if (option=="-isrc")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** cgy-GATESinoToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        // Recover the path to the output image
        string input = argv[i+1];
        int pos_colon = input.find_first_of(":");
        path_to_src_image = input.substr(0,pos_colon);
        input = input.substr(pos_colon + 1);

        // Get string section related to dimensions
        string p_dims_str[6];
        if (ReadStringOption(input.c_str(), p_dims_str, 6, ",", option))
        {
          Cerr("***** cgy-GATESinoToCastor :: Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
          Exit(EXIT_FAILURE);
        }

        // Recover dimensions
        for(int i=0 ; i<3 ; i++)
          if (ConvertFromString(p_dims_str[i], &dim_src_img[i]))
          {
            Cerr("***** cgy-GATESinoToCastor :: Conversion error for elt " << p_dims_str[i] << " for option " << option << " !" << endl);
            Exit(EXIT_FAILURE);
          }

        // Recover dimensions
        for(int i=0 ; i<3 ; i++)
          if (ConvertFromString(p_dims_str[i+3], &vox_src_img[i]))
          {
            Cerr("***** cgy-GATESinoToCastor :: Conversion error for elt " << p_dims_str[i+3] << " for option " << option << " !" << endl);
            Exit(EXIT_FAILURE);
          }

        // Initilize source image
        p_src_img = new FLTNB[dim_src_img[0]*dim_src_img[1]*dim_src_img[2]];

        for(int v=0 ; v<dim_src_img[0]*dim_src_img[1]*dim_src_img[2] ; v++)
          p_src_img[v] = 0;

        src_img_flag = true;
      }
      i++;
    }

    // Generate geom file
    else if (option=="-geo")
    {
      geom_flag = true;
    }

    // SPECT bins
    else if (option=="-sp_bins")
    {
      string input = argv[i+1];
      uint16_t s_bins[2];
      if (ReadStringOption(input.c_str(), s_bins, 2, ",", option))
      {
        Cerr("***** cgy-GATESinoToCastor :: Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
        Exit(EXIT_FAILURE);
      }

      spect_bin_trs = s_bins[0];
      spect_bin_axl = s_bins[1];

      i++;
    }


    // Verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** cgy-GATESinoToCastor :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      vb = atoi(argv[i+1]);
      i++;
    }

    else
    {
      Cerr("***** cgy-GATESinoToCastor :: Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }


  // ============================================================================================================
  // Mandatory checks:
  // ============================================================================================================

  // Basic initialization checks (minimal initializations mandatory for the next steps)

  // data files
  if (path_to_input_file.empty() )
  {
    Cerr("***** cgy-GATESinoToCastor :: Please provide at least one data filename (-i or -if)" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
  {
    /*
    // Check if we have interfile
    if(path_to_input_file.size() == 1)
    {
      string check;
      int rvalue = 0;
      rvalue = IntfKeyGetValueFromFile(path_to_input_file[0], "interfile", &check, 1, KEYWORD_OPTIONAL);

      if(rvalue == 1)
      {
        // error
        Cerr("***** cgy-GATESinoToCastor :: Error when trying to read file: " << path_to_input_file[0] << "!" << endl);
        Exit(EXIT_FAILURE);
      }
      else if(rvalue == 0)
        // key found, we have an interfile as input
        input_is_intf = true;
    }
    */
    if(vb >= 2)
    {
      Cout(" Selected root data-file(s) to convert: " << endl);
      for (size_t i=0 ; i<path_to_input_file.size() ; i++)
        Cout(path_to_input_file[i] << endl);
    }
  }



  if(!input_is_intf)
  {
    // Check ROOT is enabled if the input file is not interfile (SPECT)
    #ifndef CASTOR_ROOT
    Cerr("***** cgy-GATESinoToCastor :: CASToR must be compiled with ROOT to read input root files (check installation instructions)" << endl);
    Exit(EXIT_FAILURE);
    #endif
  }
  else if(src_img_flag) // intf input + image of the source -> error
  {
    Cerr("***** cgy-GATESinoToCastor :: Can't use -isrc with interfile dataset ");
    Cerr(" (image of the source can only be generated from root file) !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Check file(s) existence
  for (size_t i=0 ; i<path_to_input_file.size() ; i++)
  {
    ifstream input_file(path_to_input_file[i].c_str(), ios::in);
    if (!input_file)
    {
      Cerr("***** gOptions::ReadDataASCIIFile() -> Couldn't find or read data-file '"<< path_to_input_file[i] << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // output directory
  if (path_to_out_file.empty() )
  {
    Cerr("***** cgy-GATESinoToCastor :: Please provide the output file name (-o)" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
    if(vb >= 2) Cout(" selected output file:" << path_to_out_file << endl);

  // macro
  if (path_to_mac_file.empty())
  {
    Cerr("***** cgy-GATESinoToCastor :: Please provide the macro file associated to the GATE root datafile (-m) :" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
    if(vb >= 2) Cout(" selected macro file: " << path_to_mac_file << endl);

  // scanner
  if (scanner_name.empty())
  {
    Cerr("***** cgy-GATESinoToCastor :: Please provide a scanner alias (-s) :" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
    if(vb >= 2) Cout(" selected scanner alias: " << scanner_name << endl);


  // (Required for options using interfile I/O)
  GetUserEndianness();


  // ============================================================================================================
  // SOutputManager object initialisation:
  // ============================================================================================================

  sOutputManager* p_outputManager = sOutputManager::GetInstance();
  // Set verbose level
  p_outputManager->SetVerbose(vb);
  // Set MPI rank
  p_outputManager->SetMPIRank(0);

  // Set path to the config directory
  if (p_outputManager->CheckConfigDir(""))
  {
    Cerr("***** cgy-GATESinoToCastor :: A problem occured while checking for the config directory path !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize output directory and base name
  if (p_outputManager->InitOutputDirectory(path_to_out_file, ""))
  {
    Cerr("*****cgy-GATESinoToCastor ::  A problem occured while initializing output directory !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Log command line
  if (p_outputManager->LogCommandLine(argc,argv))
  {
    Cerr("***** cgy-GATESinoToCastor :: A problem occured while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Output progression options
  cout << std::fixed;
  cout << std::setprecision(2);


  // ============================================================================================================
  // Check system type from the macro file
  // ============================================================================================================
  GATE_system_type = GetGATESystemType(path_to_mac_file);

  if(GATE_system_type<0)
  {
    // Unknown system
    Cerr("***** cgy-GATESinoToCastor :: GATE system type not found : " << endl);
    Cerr("                           This script only supports conversion for cylindricalPET, ecat and SPECThead systems"  << endl);
    Cerr("                           The system type is recovered from the lines '/gate/systems/...'"  << endl);
    Exit(EXIT_FAILURE);
  }


  // ============================================================================================================
  // Generate the geom file from the mac file(s) is required
  // ============================================================================================================
  if(geom_flag)
  {
    string scanner_repository = sOutputManager::GetInstance()->GetPathToConfigDir() + "scanner" + OS_SEP;
    string path_to_geom = scanner_repository + scanner_name + ".geom";

    // Check system type
    switch ( GATE_system_type )
    {
      case GATE_SYS_CYLINDRICAL:
        if( vb>=2 )Cout(endl << " --- CylindricalPET system detected. Proceeding to conversion... --- " << endl << endl);

        if(CreateGeomWithCylindrical(path_to_mac_file , path_to_geom) )
        {
          Cerr("***** cgy-GATESinoToCastor :: An error occured while trying to process mac file for cylindrical system: " << path_to_mac_file << endl);
          Exit(EXIT_FAILURE);
        }
        break;

      case GATE_SYS_ECAT:
        if( vb>=2 )Cout(endl << " --- ECAT system detected. Proceeding to conversion... --- " << endl << endl);

        if(CreateGeomWithECAT(path_to_mac_file , path_to_geom) )
        {
          Cerr("***** castor-GATEMacToGeom :: An error occured while trying to process mac file for ecat system: " << path_to_mac_file  << endl);
          Exit(EXIT_FAILURE);
        }
        break;
      // TODO
      case GATE_SYS_SPECT:
        if( vb>=2 )Cout(endl << " --- SPECThead system detected. Proceeding to conversion... --- " << endl << endl);

        if(CreateGeomWithSPECT(path_to_mac_file , path_to_geom) )
        {
          Cerr("***** castor-GATEMacToGeom :: An error occured while trying to process mac file for SPECT system: " << path_to_mac_file  << endl);
          Exit(EXIT_FAILURE);
        }
        break;

      default: // Unknown system
        Cerr("***** cgy-GATESinoToCastor :: System type not found : " << endl);
        Cerr("                   This script only supports conversion for cylindricalPET ecat and SPECThead systems"  << endl);
        Cerr("                   The system type is recovered from the lines '/gate/systems/...'"  << endl);
        Exit(EXIT_FAILURE);
        break;
    }

    if( vb>=2 )Cout(endl << " --- Conversion completed --- " << endl << endl);
  }


  // ============================================================================================================
  // ScannerManager object initialisation:
  // ============================================================================================================

  sScannerManager* p_scannermanager = sScannerManager::GetInstance();
  p_scannermanager->SetVerbose(vb);

  // Check if the scanner exists and get the name from ScannerManager
  scanner_name = (scanner_name.find(OS_SEP)) ?
                  scanner_name.substr(scanner_name.find_last_of(OS_SEP)+1) :
                  scanner_name;

  if(p_scannermanager->FindScannerSystem(scanner_name) )
  {
    Cerr("**** cgy-GATESinoToCastor :: A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Build output file names
  path_to_data_filename = path_to_out_file + ".Cdf";
  path_to_header_filename = path_to_out_file + ".Cdh";


  // ============================================================================================================
  // Instanciate/Initialize CASToR DataFile
  // ============================================================================================================

  // Instantiate & Initialize oImageDimensionsAndQuantification object, required for datafile generation (number of threads)
  oImageDimensionsAndQuantification* p_ID = new oImageDimensionsAndQuantification();

  // Set isotope if provided
  if(!istp_alias.empty())
  {
    if(p_ID->SetPETIsotope(0, istp_alias) )
    {
      Cerr("**** cgy-GATESinoToCastor :: A problem occurred while checking isotope name !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // Instanciate & Initialize iDataFilePET and Event objects
  vDataFile* Out_data_file = NULL;
  vEvent* Event = NULL;

  uint16_t max_nb_lines_per_event = 1; // No compression for root files

  if(GATE_system_type == GATE_SYS_SPECT)
  {
    Out_data_file = new iDataFileSPECT();
    iDataFileSPECT* p_datafile = (dynamic_cast<iDataFileSPECT*>(Out_data_file));
    p_datafile->SetDataType(TYPE_SPECT);
    p_datafile->SetIsotope(istp_alias);
    histo_out_flag = true; // force histogram output for SPECT
  }
  else
  {
    Out_data_file = new iDataFilePET();
    iDataFilePET* p_datafile = (dynamic_cast<iDataFilePET*>(Out_data_file));
    p_datafile->SetDataType(TYPE_PET);
    p_datafile->SetIsotope(istp_alias);

    // ACF computed
    if(estimate_acf_flag)
      p_datafile->SetAtnCorrectionFlagOn();

    p_datafile->SetMaxNumberOfLinesPerEvent(max_nb_lines_per_event);
  }

  Out_data_file->SetImageDimensionsAndQuantification(p_ID);
  Out_data_file->SetHeaderDataFileName(path_to_header_filename);
  // Out_data_file->SetPercentageLoad(0); // 0 (default)
  Out_data_file->SetVerbose(0);


  // Init Histogram-mode Event
  if(histo_out_flag)
  {
    Out_data_file->SetDataMode(MODE_HISTOGRAM);

    // Instanciate histogram Event depending on modality
    if(GATE_system_type == GATE_SYS_SPECT)
      Event = new iEventHistoSPECT();
    else
      Event = new iEventHistoPET();

    Event->SetNbLines(max_nb_lines_per_event);
    Event->AllocateID();
  }

  // Init List-mode Event
  else
  {
    Out_data_file->SetDataMode(MODE_LIST);

    // Instanciate histogram Event depending on modality
    if(GATE_system_type == GATE_SYS_SPECT)
    {
      Event = new iEventListSPECT();
      // record coincidence kind or not
      if(kind_flag)
        ((iDataFileSPECT*)Out_data_file)->SetEventKindFlagOn();
    }
    else
    {
      Event = new iEventListPET();
      // record coincidence kind or not
      if(kind_flag)
        ((iDataFilePET*)Out_data_file)->SetEventKindFlagOn();
    }

    Event->SetNbLines(max_nb_lines_per_event);
    Event->AllocateID();
  }

  Out_data_file->PROJ_InitFile();
  Out_data_file->ComputeSizeEvent();
  Out_data_file->PrepareDataFile();


  // ============================================================================================================
  // If acf estimation is enabled for histogram output, initialize all required objects for analytic projection
  // (Image space, projector and scanner geometry)
  // ============================================================================================================
  oProjectorManager* p_ProjectorManager = new oProjectorManager();
  oImageSpace* p_ImageSpace = new oImageSpace();


  if(estimate_acf_flag)
  {
    // Check if histogram output has been enabled, throw error otherwise
    if(!histo_out_flag)
    {
      Cerr("***** cgy-GATESinoToCastor :: Estimation of acf from an attenuation image (-atn option) is only available for histogram datafile output !" << endl
        << "                                 Add the (-oh) option  to the command line to enable this option." << endl);
      Exit(1);
    }

    // Check if system is SPECT, throw error in this case
    if(GATE_system_type == GATE_SYS_SPECT)
    {
      Cerr("***** cgy-GATESinoToCastor :: Estimation of acf from an attenuation image (-atn option) only available for PET systems ! (detected system is SPECT)" << endl);
      Exit(1);
    }

    Intf_fields IF;
    IntfKeyInitFields(&IF);
    if(IntfReadHeader(path_to_atn_image, &IF, vb) )
    {
      Cerr("***** cgy-GATESinoToCastor :: An error occurred while trying to read the interfile header of attenuation file " << path_to_atn_image << " !" << endl);
      Exit(1);
    }

    // --- oImageDimensionsAndQuantification initialization ---
    p_ID->SetNbVoxX(IF.mtx_size[0]);
    p_ID->SetNbVoxY(IF.mtx_size[1]);
    p_ID->SetNbVoxZ(IF.mtx_size[2]);
    p_ID->SetNbThreads("1");
    p_ID->SetNbBeds(1);
    p_ID->SetVoxSizeX(IF.vox_size[0]);
    p_ID->SetVoxSizeY(IF.vox_size[1]);
    p_ID->SetVoxSizeZ(IF.vox_size[2]);
    p_ID->SetFOVOutMasking(0., 0);
    p_ID->SetFOVSizeX(-1.);
    p_ID->SetFOVSizeY(-1.);
    p_ID->SetFOVSizeZ(-1.);
    p_ID->SetOffsetX(0);
    p_ID->SetOffsetY(0);
    p_ID->SetOffsetZ(0);
    p_ID->SetVerbose(vb);
    p_ID->SetNbRespGates(1);
    p_ID->SetNbCardGates(1);
    p_ID->SetFrames("");

    if (p_ID->CheckParameters())
    {
      Cerr("***** cgy-GATESinoToCastor :: A problem occured while checking image dimensions parameters !" << endl);
      Exit(1);
    }
    if (p_ID->Initialize())
    {
      Cerr("***** cgy-GATESinoToCastor :: A problem occured while initializing image dimensions !" << endl);
      Exit(1);
    }

    // Initialization of DynamicDataManager class, related 4D data splitting management
    if (p_ID->InitDynamicData( "", 0, 0, 0, 0, 1, 1 ) )
    {
      Cerr("***** cgy-GATESinoToCastor :: A problem occured while initializing Dynamic data manager's class !" << endl);
      Exit(EXIT_FAILURE);
    }


    // --- Image space initialization ---
    p_ImageSpace->SetImageDimensionsAndQuantification(p_ID);
    p_ImageSpace->SetVerbose(vb);

    // Allocate memory for main image
    p_ImageSpace->LMS_InstantiateImage();

    // Read attenuation image
    if(p_ImageSpace->PROJ_LoadInitialImage(path_to_atn_image) )
    {
      Cerr("***** cgy-GATESinoToCastor :: Error during image initialization !" << endl);
      Exit(EXIT_FAILURE);
    }


    // --- Build Scanner geometry ---
    if(p_scannermanager->BuildScannerObject() )
    {
      Cerr("**** cgy-GATESinoToCastor :: A problem occurred during scanner object construction !" << endl);
      Exit(EXIT_FAILURE);
    }

    if(p_scannermanager->InstantiateScanner() )
    {
      Cerr("**** cgy-GATESinoToCastor :: A problem occurred while creating Scanner object !" << endl);
      Exit(EXIT_FAILURE);
    }

    if(p_scannermanager->BuildLUT() )
    {
      Cerr("***** cgy-GATESinoToCastor :: A problem occurred while generating/reading the LUT !" << endl);
      Exit(EXIT_FAILURE);
    }

    // Check the scanner manager parameters and initialize the scanner
    if (p_scannermanager->CheckParameters())
    {
      Cerr("***** cgy-GATESinoToCastor :: A problem occured while checking scanner manager parameters !" << endl);
      Exit(1);
    }

    if (p_scannermanager->Initialize())
    {
      Cerr("***** cgy-GATESinoToCastor :: A problem occured while initializing scanner !" << endl);
      Exit(1);
    }


    // --- Projector Manager initialization---
    p_ProjectorManager->SetScanner(p_scannermanager->GetScannerObject());
    p_ProjectorManager->SetImageDimensionsAndQuantification(p_ID);
    p_ProjectorManager->SetDataFile(Out_data_file);
    p_ProjectorManager->SetComputationStrategy(FIXED_LIST_COMPUTATION_STRATEGY);
    p_ProjectorManager->SetOptionsForward(options_projector);
    p_ProjectorManager->SetOptionsBackward(options_projector);
    p_ProjectorManager->SetOptionsCommon("");
    p_ProjectorManager->SetVerbose(vb);

    // Check parameters
    if (p_ProjectorManager->CheckParameters())
    {
      Cerr("***** cgy-GATESinoToCastor :: A problem occured while checking projector manager's parameters !" << endl);
      Exit(EXIT_FAILURE);
    }

    // Initialize projector manager
    if (p_ProjectorManager->Initialize())
    {
      Cerr("***** cgy-GATESinoToCastor :: A problem occured while initializing projector manager !" << endl);
      Exit(EXIT_FAILURE);
    }
  }


  // ============================================================================================================
  // Parse Time offsets
  // ============================================================================================================
  offset_time_list = new uint32_t[path_to_input_file.size()];
  for(size_t f=0 ; f<path_to_input_file.size() ; f++)
    offset_time_list[f] = 0;

  // Parse offset_time_str, if it contains any data
  if(offset_time_str != "")
  {
    vector<string> offsets;
    size_t comma_pos = 0;
    while ((comma_pos=offset_time_str.find_first_of(",")) != string::npos)
    {
      string offset = offset_time_str.substr(0,comma_pos);
      offsets.push_back(offset);
      offset_time_str = offset_time_str.substr(comma_pos+1);
    }

    // Check we have a correct number of offsets
    if(offsets.size() != path_to_input_file.size())
    {
      Cerr("**** cgy-GATESinoToCastor :: Unmatching number of offsets with -ot option ! "
           << offsets.size() << " have been found, while "<< path_to_input_file.size() <<"input file have been provided !" << endl);
      Exit(EXIT_FAILURE);
    }

    for(size_t o=0 ; o<offsets.size() ; o++)
    {
      double offset;
      if(ConvertFromString(offsets[o] , &offset) )
      {
        Cerr("**** cgy-GATESinoToCastor :: Error while trying to convert offset : "<< offsets[o] <<" in ms ! " << endl);
        Exit(EXIT_FAILURE);
      }

      // convert in ms
      offset_time_list[o] = (uint32_t)offset*1000;
    }
  }


  // ============================================================================================================
  // *********************************************** CONVERSION *************************************************
  // ============================================================================================================
  Cout(" --- Start conversion of datafile(s) : " << input_file <<  endl
    << "     using mac file: " << path_to_mac_file << endl
    << "     CASToR output header datafile: " << path_to_header_filename << endl
    << "     CASToR output binary datafile: " << path_to_data_filename << endl << endl);

  // ============================================================================================================
  // Variables initialization
  // ============================================================================================================

  // Counter for the number of events (and the number of trues, scatters and randoms for simulated data)
  uint64_t nLORs_tot = 0,
           nLORs_trues = 0,
           nLORs_rdms = 0,
           nLORs_scatters =0,
           nLORs_unknown =0,
           nBINs = 0;

  // Scanner variables (PET)
  uint32_t nCrystalsTot = 0;
  uint32_t nRsectorsAngPos = 1, nRsectorsAxial = 1;
  // index order of rsectors (transaxial/axial) for cylindricalPET depending on the GATE macro
  int      rsector_id_order = 0;
  uint32_t nModulesTransaxial = 1, nModulesAxial = 1;
  uint32_t nSubmodulesTransaxial = 1, nSubmodulesAxial = 1;
  uint32_t nCrystalsTransaxial = 1, nCrystalsAxial = 1;
  uint32_t nBlocksLine = 1, nBlocksPerRing = 1;
  uint8_t  nLayers = 1;
  uint32_t nb_crystal_per_layer[4] = {0,0,0,0};

  // Scanner variables (SPECT)
  uint32_t nProjectionsByHead =1;
  uint32_t nProjectionsTot = 1;
  uint32_t nHeads =1;
  uint32_t nPixelsAxial = 1;
  uint32_t nPixelsTransaxial = 1;
  uint32_t nbSimulatedPixels = 1;
  float_t  distToDetector = 0.;
  int      headRotDirection = GEO_ROT_CW;
  float_t  head1stAngleDeg = 0;
  float_t  headAngPitchDeg = -1;
  float_t  headAngStepDeg = -1;
  float_t  crystalSizeAxl=-1.,
           crystalSizeTrs=-1.;
  FLTNB*   p_proj_spect_image=NULL;

  // layer repeaters
  vector<uint32_t> nLayersRptTransaxial, nLayersRptAxial;

  // Castor variables
  uint8_t** p_bins = NULL;
  uint32_t bins_elts = 0;
  uint32_t start_time_ms = 0; // 0 as default initialization
  uint32_t duration_ms = 1000; // 1 second as default initialization


  #ifdef CASTOR_ROOT
  uint32_t castorID1=0;
  uint32_t castorID2=0;
  uint8_t kind;

  // ROOT data variables
  int32_t crystalID1=0, crystalID2=0;

  // cylindricalPET specific
  int32_t rsectorID1=0  , rsectorID2=0;
  int32_t moduleID1=0   , moduleID2=0;
  int32_t submoduleID1=0, submoduleID2=0;
  int32_t layerID1=0    , layerID2=0;

  // ecat specific
  int32_t blockID1=0, blockID2=0;

  // SPECThead specific
  float_t rotAngle;
  float_t gPosX, gPosY, gPosZ;
  int32_t headID=0;
  int32_t pixelID=0;


  // others
  int32_t eventID1= 0, eventID2= 0;
  int32_t comptonPhantomID1 = 0, comptonPhantomID2= 0;
  int32_t rayleighPhantomID1 = 0,rayleighPhantomID2 = 0;
  double_t time1= 0, time2= 0;
  int32_t sourceID1=0, sourceID2=0;
  float_t sourcePosX1=0., sourcePosY1=0., sourcePosZ1=0.;


  // ============================================================================================================
  // ROOT objects declaration
  // ============================================================================================================

//  TApplication *Tapp = new TApplication("tapp",0,0);
//  TTree** GEvents = new TTree *[path_to_input_file.size()];
//  TFile** Tfile_root = new TFile*[path_to_input_file.size()];

//  if(!input_is_intf)
//  {
//    // Compute the total number of LORs in the dataset
//    for (size_t iFic=0 ; iFic<path_to_input_file.size() ; iFic++)
//    {
//      Tfile_root[iFic] = new TFile(path_to_input_file[iFic].c_str(),"READ","ROOT file with histograms");
//      if(GATE_system_type == GATE_SYS_SPECT)
//        GEvents[iFic] = (TTree*)Tfile_root[iFic]->Get("Singles");
//      else
//        GEvents[iFic] = (TTree*)Tfile_root[iFic]->Get("Coincidences");

//      nLORs_tot += GEvents[iFic]->GetEntries();
//    }
//  }
  #endif


  // Indexes for progression output
  uint64_t printing_index = 0;
  uint64_t printing_ratio = (nLORs_tot>10000) ? 10000 : nLORs_tot/10;

  // ============================================================================================================
  // Recover system geometric elements values
  // ============================================================================================================
  // PET systems
  {
    // cylindricalPET system
    if(GATE_system_type == GATE_SYS_CYLINDRICAL)
    {
      if(ReadMacCylindrical(path_to_mac_file,
                                     nLayers,
                        nb_crystal_per_layer,
                                nCrystalsTot,
                              nCrystalsAxial,
                         nCrystalsTransaxial,
                             nLayersRptAxial,
                        nLayersRptTransaxial,
                            nSubmodulesAxial,
                       nSubmodulesTransaxial,
                               nModulesAxial,
                          nModulesTransaxial,
                              nRsectorsAxial,
                             nRsectorsAngPos,
                            rsector_id_order,
                               start_time_ms,
                                 duration_ms,
                                          vb) )
      {
        Cerr("**** cgy-GATESinoToCastor :: Error when reading mac file : "<< path_to_mac_file <<" !" << endl);
        Exit(EXIT_FAILURE);
      }
    }

    else // ECAT
    {
      // Reading the macro file
      if(ReadMacECAT(path_to_mac_file,
                         nCrystalsTot,
                       nCrystalsAxial,
                  nCrystalsTransaxial,
                          nBlocksLine,
                       nBlocksPerRing,
                        start_time_ms,
                          duration_ms,
                                   vb) )
      {
        Cerr("**** cgy-GATESinoToCastor :: Error when reading mac file : "<< path_to_mac_file <<" !" << endl);
        Exit(EXIT_FAILURE);
      }
    }

    if(vb>=2) Cout("Detected number of crystals in the system : " << nCrystalsTot << endl);

      // Histogram bin vector initialization using the total number of crystals
    if(histo_out_flag)
    {
      Cout(" Allocating memory for bins... " << endl <<
           " Warning : this step can require huge amount of memory if the system contains a high number of crystals !" << endl);

      bins_elts = nCrystalsTot;

      p_bins = new uint8_t*[bins_elts];
      for (size_t c=0; c<bins_elts ; c++)
      {
        p_bins[c] = new uint8_t[nCrystalsTot-c];

        for (size_t c2=0; c2<nCrystalsTot-c ; c2++)
          p_bins[c][c2] = 0;
      }

      Cout(" Memory allocation for bins completed " << endl << endl);
    }

  }  // end of PET section

/*!
  to be modified
  */
  // ============================================================================================================
  // Loop on the input datafile(s) and process event by event
  // ============================================================================================================
  if(!input_is_intf) // SPECT projection interfile : data already processed
    for (size_t iFic=0 ; iFic<path_to_input_file.size() ; iFic++)
    {
      if(vb>=2) Cout(endl << "Processing datafile " << iFic << " : " << path_to_input_file[iFic] << "..." << endl);

      // Set variables of the root tree
      // If we got a cylindricalPET system
      #ifdef CASTOR_ROOT
      int num_CrystPerRing = nBlocksPerRing*nCrystalsTransaxial;
      int num_Ring = nBlocksLine*nCrystalsAxial;
      int nBinView = num_CrystPerRing/2;
      int nBinElem = num_CrystPerRing/2;
      //read the sinogram;
      int sinoSize = nBinView*nBinElem*num_Ring*num_Ring;
      vector<uint16_t> sinoValue(sinoSize);
      std::ifstream sinoInfile(path_to_input_file[iFic], std::ios_base::binary|std::ios_base::in);
      for (int i = 0; i < sinoSize;i++ ){
          sinoInfile.read((char*)&sinoValue[i],sizeof(uint16_t));
      }
      sinoInfile.close();

      int size2Dsino = nBinView * nBinElem;
      //randomly use the sinogram
      std::vector<int> sinoRandomIndex;
      for (int i = 0 ; i < size2Dsino; i++){
          sinoRandomIndex.push_back(i);
      }
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::shuffle (sinoRandomIndex.begin (), sinoRandomIndex.end (), std::default_random_engine (seed));

      std::vector<LOR> event_list;

      int sinoID = 0;
      int crystID1, crystID2;
      int ringID1,ringID2;
      for (sinoID = 0; sinoID < num_Ring; sinoID++){
        // read a 2D sinogram of GATE and convert to castor list mode data.
        SinoIDToRingID(sinoID,num_Ring,ringID1,ringID2);
        for(int iBin = 0; iBin < size2Dsino ;iBin++){
                if(sinoValue[sinoID*size2Dsino + sinoRandomIndex[iBin] ] > 0){
                    BinIDToCrystID(sinoRandomIndex[iBin],nBinElem,num_CrystPerRing,crystID1,crystID2);
                    ConvertToGATEEcatID(nBlocksPerRing,nCrystalsTransaxial,ringID1,ringID2,crystID1,crystID2);
                    castorID1 = crystID1 + ringID1 * num_CrystPerRing;
                    castorID2 = crystID2 + ringID2 * num_CrystPerRing;
                    for( int i = 0; i < sinoValue[sinoID*size2Dsino + sinoRandomIndex[iBin] ] ;i++){

                        LOR lor;
                        lor.crystID1 = castorID1;
                        lor.crystID2 = castorID2;
                        event_list.push_back(lor);

                        kind = KIND_TRUE;
                        // Count nb LORs according to kind
                        switch (kind)
                        {
                          case 1:
                            nLORs_trues++;
                            break;

                          case 2:
                            nLORs_scatters++;
                            break;

                          case 3:
                            nLORs_scatters++;
                            break;

                          case 4:
                            nLORs_rdms++;
                            break;

                          default:
                            nLORs_unknown++;
                        }

                        // Skip next step if event is not true if we only recover true
                        if (true_only_flag==true && kind != KIND_TRUE)
                          continue;

                        // --- Write Event ---

                        // SPECT event

                    }

                }
        }
      }

      nLORs_tot = nLORs_trues;
      unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
      std::shuffle (event_list.begin (), event_list.end (), std::default_random_engine (seed1));
      for(int iEvent = 0; iEvent < event_list.size();iEvent++){
          LOR lor = event_list[iEvent];
          castorID1 = lor.crystID1;
          castorID2 = lor.crystID2;
          if(GATE_system_type == GATE_SYS_SPECT)
          {
            // HISTOGRAM mode
            if(histo_out_flag)
              p_bins[castorID1][castorID2]++;
            else
            // LIST-mode
            {
              // Write event in the datafile
              uint32_t time_in_ms = time1*1000;
              Event->SetNbLines(1);
              Event->SetID1(0, castorID1); // 1st argument : line, 2nd argument :index
              Event->SetID2(0, castorID2); // 1st argument : line, 2nd argument :index
              Event->SetTimeInMs(time_in_ms);

             ((iEventListSPECT*)Event)->SetKind(kind);

              Out_data_file->PROJ_WriteEvent(Event, 0);
            }
          }
          else // PET event
          {
            // HISTOGRAM mode
            if(histo_out_flag)
            {
              if (castorID1 < castorID2)
                p_bins[castorID1][castorID2-castorID1-1]++;
              else
                p_bins[castorID2][castorID1-castorID2-1]++;
            }
            // LIST-mode
            else
            {
              // Write event in the datafile
              uint32_t time_in_ms = time1*1000;
              int nb_lines_in_event = 1; // 1 by default for GATE root files

              Event->SetNbLines(nb_lines_in_event);
              Event->SetID1(0, castorID1); // 1st argument : line, 2nd argument :index
              Event->SetID2(0, castorID2); // 1st argument : line, 2nd argument :index
              Event->SetTimeInMs(time_in_ms);

             ((iEventListPET*)Event)->SetKind(kind);

             Out_data_file->PROJ_WriteEvent(Event, 0);
            }
          }
      }

      #endif

      if(vb>=2) Cout(endl << "DataFile " << iFic << " : " << path_to_input_file[iFic] << " process OK" << endl);

    } // end of loop on input datafiles

  // Free Root objects
  #ifdef CASTOR_ROOT
//  delete[] Tfile_root;
//  delete[] GEvents;
//  delete Tapp;
  #endif



  // ============================================================================================================
  // Write the HISTOGRAM datafile and header
  // ============================================================================================================
  {
    if (histo_out_flag)
    {
      printing_index=0;

      // Writing the datafile
      Cout(endl << "Generate the histogram datafile" << endl);
      uint32_t nb_bins = (nCrystalsTot*nCrystalsTot - nCrystalsTot)/2;
      printing_ratio = (nb_bins>1000) ? 1000 : nb_bins/10;

      // Loop on the crystal IDs
      for (size_t id1=0 ; id1 <bins_elts ; id1++)
        for (size_t id2 = id1+1; id2 < nCrystalsTot;id2++)
        {
          uint32_t nb_events_in_bin = p_bins[id1][id2-id1-1];

          int nb_lines = 1; // 1 by default for GATE root file;
          Event->SetNbLines(nb_lines);
          Event->SetID1(0, id1);  // 1st argument : line, 2nd argument :index
          Event->SetID2(0, id2);  // 1st argument : line, 2nd argument :index
          Event->SetEventValue(0, nb_events_in_bin);

          // Estimate acf for the bin
          if(estimate_acf_flag)
          {
            // Compute the system matrix elements for the two crystals
            oProjectionLine* line = p_ProjectorManager->ComputeProjectionLine(Event, 0);

            // Compute forward projection of attenuation image
            FLTNB fp = 0.;
            if (line->NotEmptyLine())
              fp = line->ForwardProject(p_ImageSpace->m4p_image[0][0][0]) * 0.1; // 0.1 -> convert in mm-1

            // Write atn correction factor in Event
            ((iEventPET*)Event)->SetAttenuationCorrectionFactor(1/std::exp(-fp));
          }

          // Write event in DataFile
          Out_data_file->WriteEvent(Event, 0);
          nBINs++;

          // Progression
          if (printing_index%((nb_bins)/printing_ratio) == 0)
          {
            FLTNB percent = ( ((FLTNB)(printing_index+1))/((FLTNB)nb_bins) ) * ((FLTNB)100);
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
              << percent << " \%                    ";
          }
          printing_index++;
        }

      Out_data_file->SetStartTime((FLTNB)(start_time_ms)/1000); // start time in s
      Out_data_file->SetDuration((FLTNB)(duration_ms)/1000); // duration in s
      Out_data_file->SetNbEvents( nBINs );
      Out_data_file->WriteHeader();

      Cout(endl << "The output histogram contains " << nBINs  << " events." << endl);
    }
    // Write the List-mode datafile header
    else
    {
      Out_data_file->SetStartTime((FLTNB)(start_time_ms)/1000); // start time in s
      Out_data_file->SetDuration((FLTNB)(duration_ms)/1000); // duration in s
      //uint32_t nLORs = (true_only_flag) ? nLORs_trues : nLORs_tot;
      Out_data_file->SetNbEvents(nLORs_trues);
      Out_data_file->WriteHeader();
    }
  }  // end of PET section


  // Write the source image
  if(src_img_flag)
  {
    if (vb>=2) Cout("Writing source image ..." << endl);

    // Initialize Intf_fields object with the source image dimensions
    Intf_fields IF;
    IntfKeyInitFields(&IF);

    IF.mtx_size[0] = dim_src_img[0];
    IF.mtx_size[1] = dim_src_img[1];
    IF.mtx_size[2] = dim_src_img[2];
    IF.vox_size[0] = vox_src_img[0];
    IF.vox_size[1] = vox_src_img[1];
    IF.vox_size[2] = vox_src_img[2];
    IF.nb_total_imgs = dim_src_img[2];

    if (IntfWriteImgFile(path_to_src_image.append("_src"), p_src_img, IF, vb) )
    {
      Cerr("***** cgy-GATESinoToCastor :: Error writing Interfile of output image !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  Cout(endl << "The simulated dataset contained " << nLORs_tot      << " LORs. " << endl);

  if (vb>=2) Cout("Writing raw datafile ..." << endl); // todo just change name if only one datafile


  Out_data_file->PROJ_WriteData();
  Out_data_file->PROJ_DeleteTmpDataFile();

  Cout(endl << " --- Conversion completed --- " << endl);

  // ============================================================================================================
  // End
  // ============================================================================================================

  if (vb>=2) Cout(" Deallocating objects ..." << endl);

  // Delete objects
  if(histo_out_flag)
  {
    for(size_t b=0; b<bins_elts ; b++)
      delete[] p_bins[b];

    delete[]p_bins;

    if(Event)
    {
      if(GATE_system_type == GATE_SYS_SPECT)
        delete (iEventHistoSPECT*)Event;
      else
        delete (iEventHistoPET*)Event;
    }
  }
  else
    if(Event)
    {
      if(GATE_system_type == GATE_SYS_SPECT)
        delete (iEventListSPECT*)Event;
      else
        delete (iEventListPET*)Event;
    }

  delete[]offset_time_list;
  delete Out_data_file;

  if(src_img_flag && p_src_img)
    delete[] p_src_img;

  // Free objects created for analytic projection (acf estimation)
  if(estimate_acf_flag)
    p_ImageSpace->LMS_DeallocateImage();
  if(p_ImageSpace)       delete p_ImageSpace;
  if(p_ProjectorManager) delete p_ProjectorManager;
  if(p_ID)               delete p_ID;

  Cout(" ---          END         --- " << endl << endl);

  return 0;
}
