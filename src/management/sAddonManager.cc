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

  \brief Implementation of class sAddonManager
*/

#include "sAddonManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

sAddonManager* sAddonManager::mp_Instance = NULL;
sAddonManager::sAddonManager() {;}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sAddonManager::ShowHelpProjector()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Get the list of projectors
  std::map<string,maker_projector> list = sAddonManager::GetInstance()->mp_listOfProjectors;
  std::map<string,maker_projector>::iterator iter;
  cout << endl << "Here is the list of all implemented projectors along with their options:" << endl << endl;
  for (iter = list.begin(); iter!=list.end(); iter++)
  {
    // Print out the name of this projector
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << iter->first << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    if (iter->second)
    {
      // Create it
      vProjector* projector = iter->second();
      // Print specific help
      projector->ShowHelp();
      // Delete it
      delete projector;
    }
    cout << endl;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sAddonManager::ShowHelpOptimizer()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Get the list of optimizers
  std::map<string,maker_optimizer> list = sAddonManager::GetInstance()->mp_listOfOptimizers;
  std::map<string,maker_optimizer>::iterator iter;
  cout << endl << "Here is the list of all implemented optimizers along with their options:" << endl << endl;
  for (iter = list.begin(); iter!=list.end(); iter++)
  {
    // Print out the name of this optimizer
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << iter->first << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    if (iter->second)
    {
      // Create it
      vOptimizer* optimizer = iter->second();
      // Print specific help
      optimizer->ShowHelp();
      // Delete it
      delete optimizer;
    }
    cout << endl;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sAddonManager::ShowHelpScanner()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Get the list of scanners
  std::map<string,maker_scanner> list = sAddonManager::GetInstance()->mp_listOfScannerTypes;
  std::map<string,maker_scanner>::iterator iter;
  cout << endl << "Here is the list of all implemented systems along with their options:" << endl << endl;
  for (iter = list.begin(); iter!=list.end(); iter++)
  {
    // Print out the name of this scanner
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << iter->first << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    if (iter->second)
    {
      // Create it
      vScanner* scanner = iter->second();
      // Print specific help
      scanner->ShowHelp();
      // Delete it
      delete scanner;
    }
    cout << endl;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sAddonManager::ShowHelpImageConvolver()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Get the list of image convolvers
  std::map<string,maker_image_convolver> list = sAddonManager::GetInstance()->mp_listOfImageConvolvers;
  std::map<string,maker_image_convolver>::iterator iter;
  cout << endl << "Here is the list of all implemented image convolvers along with their options:" << endl << endl;
  for (iter = list.begin(); iter!=list.end(); iter++)
  {
    // Print out the name of this convolver
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << iter->first << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    if (iter->second)
    {
      // Create it
      vImageConvolver* convolver = iter->second();
      // Print specific help
      convolver->ShowHelp();
      // Delete it
      delete convolver;
    }
    cout << endl;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sAddonManager::ShowHelpImageProcessingModule()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Get the list of image processing modules
  std::map<string,maker_image_processing_module> list = sAddonManager::GetInstance()->mp_listOfImageProcessingModules;
  std::map<string,maker_image_processing_module>::iterator iter;
  cout << endl << "Here is the list of all implemented image processing modules along with their options:" << endl << endl;
  for (iter = list.begin(); iter!=list.end(); iter++)
  {
    // Print out the name of this processing module
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << iter->first << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    if (iter->second)
    {
      // Create it
      vImageProcessingModule* module = iter->second();
      // Print specific help
      module->ShowHelp();
      // Delete it
      delete module;
    }
    cout << endl;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sAddonManager::ShowHelpPenalty()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Get the list of penalties
  std::map<string,maker_penalty> list = sAddonManager::GetInstance()->mp_listOfPenalties;
  std::map<string,maker_penalty>::iterator iter;
  cout << endl << "Here is the list of all implemented penalties along with their options:" << endl << endl;
  for (iter = list.begin(); iter!=list.end(); iter++)
  {
    // Print out the name of this penalty
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << iter->first << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    if (iter->second)
    {
      // Create it
      vPenalty* penalty = iter->second();
      // Print specific help
      //penalty->ShowHelp();
      // Delete it
      delete penalty;
    }
    cout << endl;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sAddonManager::ShowHelpDynamicModel()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Get the list of dynamic models
  std::map<string,maker_dynamic_model> list = sAddonManager::GetInstance()->mp_listOfDynamicModels;
  std::map<string,maker_dynamic_model>::iterator iter;
  cout << endl << "Here is the list of all implemented models along with their options:" << endl << endl;
  for (iter = list.begin(); iter!=list.end(); iter++)
  {
    // Print out the name of this dynamic model
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << iter->first << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    if (iter->second)
    {
      // Create it
      vDynamicModel* dynamic_model = iter->second();
      // Print specific help
      dynamic_model->ShowHelp();
      // Delete it
      delete dynamic_model;
    }
    cout << endl;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void sAddonManager::ShowHelpDeformation()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Get the list of deformations
  std::map<string,maker_deformation> list = sAddonManager::GetInstance()->mp_listOfDeformations;
  std::map<string,maker_deformation>::iterator iter;
  cout << endl << "Here is the list of all implemented image deformation algorithms along with their options:" << endl << endl;
  for (iter = list.begin(); iter!=list.end(); iter++)
  {
    // Print out the name of this deformation
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << iter->first << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    if (iter->second)
    {
      // Create it
      vDeformation* deformation = iter->second();
      // Print specific help
      deformation->ShowHelp();
      // Delete it
      delete deformation;
    }
    cout << endl;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
