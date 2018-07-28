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

  \brief Implementation of class sRandomNumberGenerator
*/


#include "sRandomNumberGenerator.hh"

sRandomNumberGenerator *sRandomNumberGenerator::mp_Instance = NULL;

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \brief Constructor of sRandomNumberGenerator. Do nothing by default as it is a singleton clasee
*/
sRandomNumberGenerator::sRandomNumberGenerator()
{
  mp_Instance = NULL;
  m_verbose = 5;
  mp_Engines.clear();
  mp_extraEngines.clear();
  m_seed = -1;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief Destructor of sRandomNumberGenerator. Do nothing by default
*/
sRandomNumberGenerator::~sRandomNumberGenerator()
{
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int sRandomNumberGenerator::Initialize(int a_nbThreads, int a_nbExtra)
{
  if(m_verbose >=3) Cout("sRandomNumberGenerator::Initialize ..."<< endl); 

  int64_t seed = -1;

  // call a real random generator, only once because this operating system
  // random source might be limited in size and time
  std::random_device rd;

  // save the initial seed
  m_seed = rd();

  #ifdef CASTOR_MPI
  int mpi_rank_temp = 0;
  int mpi_size_temp = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_temp);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_temp);

  // helper variable
  int64_t temp_seed = -1;

  // if more than one MPI instance, generate seeds for each instance 
  // and dispatch them, otherwise keep the initial seed
  if (mpi_size_temp>1)
  {
    if (mpi_rank_temp==0)
    {
      Engine mpi_generator(m_seed);
      for (int p=0; p<mpi_size_temp; p++) 
      {
        m_seed = mpi_generator();

         Cout("sRandomNumberGenerator::Seed for rank " << p << " is " << m_seed << endl); 

        if (p==0)
          temp_seed = m_seed;
        else
          MPI_Send(&m_seed, 1, MPI_LONG, p, 0, MPI_COMM_WORLD);
      }
      m_seed = temp_seed;
    }
    else
    {
      MPI_Status* status = NULL;
      MPI_Recv(&m_seed, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, status);
    }
    // wait for all the processes to have their seeds
    MPI_Barrier(MPI_COMM_WORLD);
  }
  #endif

  // seed the initial temporary PRNG, which will generate seeds for actual PRNGs (threaded and non-threaded)
  Engine initialGenerator(m_seed);

  for (int th=0; th<a_nbThreads; th++)
  {
    seed = initialGenerator();
    if (m_verbose>=3) Cout("sRandomNumberGenerator::Initialize() -> Seed for thread " << th << " : " << seed << endl);
    mp_Engines.push_back(Engine(seed));
  }

  for(int ex=0; ex<a_nbExtra; ex++)
  {
    seed = initialGenerator();
    if(m_verbose>=3) Cout("sRandomNumberGenerator::Initialize() - >Seed for additional nonthreaded generator "<<ex<<" : " << seed << endl);
    mp_extraEngines.push_back(Engine(seed));
  }

  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int sRandomNumberGenerator::Initialize(int64_t a_seed, int a_nbThreads, int a_nbExtra)
{
  if(m_verbose >=3) Cout("sRandomNumberGenerator::Initialize with provided seed "<<a_seed<< endl); 

  // temporary variable
  int64_t seed = -1;

  if (a_seed<0)
  {
    Cout("***** sRandomNumberGenerator::Initialize()-> Error : seed for RNG should be >=0 !" << endl);
    return 1;
  }

  // save the initial seed
  m_seed = a_seed;

  #ifdef CASTOR_MPI
  int mpi_rank_temp = 0;
  int mpi_size_temp = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_temp);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_temp);

  // helper variable
  int64_t temp_seed = -1;

  // if more than one MPI instance, generate seeds for each instance 
  // and dispatch them, otherwise keep the initial seed
  if (mpi_size_temp>1)
  {
    if (mpi_rank_temp==0)
    {
      Engine mpi_generator(a_seed);
      for (int p=0; p<mpi_size_temp; p++) 
      {
        m_seed = mpi_generator();
     
         Cout("sRandomNumberGenerator::Seed for rank "<<p<<" is "<<m_seed<<""<< endl); 

        if (p==0)
          temp_seed = m_seed;
        else
          MPI_Send(&m_seed, 1, MPI_LONG, p, 0, MPI_COMM_WORLD);
      }
      m_seed = temp_seed;
    }
    else
    {
      MPI_Status* status = NULL;
      MPI_Recv(&m_seed, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  #endif

  // seed the initial temporary PRNG, which will generate seeds for actual PRNGs (threaded and non-threaded)
  Engine initialGenerator(m_seed);

  for(int th=0; th<a_nbThreads; th++)
  {
    seed = initialGenerator();
    mp_Engines.push_back(Engine(seed));
    if (m_verbose>=3) Cout("sRandomNumberGenerator::Initialize() -> Seed for thread " << th << " : " << seed << endl);
  }

  for(int ex=0; ex<a_nbExtra; ex++)
  {
    seed = initialGenerator();
    mp_extraEngines.push_back(Engine(seed));
    if (m_verbose>=3) Cout("sRandomNumberGenerator::Initialize() -> Seed for additional nonthreaded random generator "<<ex<< " : " << seed << endl);
  }

  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

HPFLTNB sRandomNumberGenerator::GenerateRdmNber()
{
  #ifdef CASTOR_VERBOSE
  if (m_verbose >=4) Cout("sRandomNumberGenerator::GenerateRdmNber..." << endl); 
  #endif

  int id=0;
  #ifdef CASTOR_OMP
  id = omp_get_thread_num();
  #endif
  return mp_Distribution(mp_Engines[id]);
} 



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

HPFLTNB sRandomNumberGenerator::GenerateExtraRdmNber(int a_nb)
{
  #ifdef CASTOR_VERBOSE
  if(m_verbose >=4) Cout("sRandomNumberGenerator::GenerateExtraRdmNber..."<< endl); 
  #endif
  
  return mp_Distribution(mp_extraEngines[a_nb]);
} 


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

sRandomNumberGenerator::Engine& sRandomNumberGenerator::GetExtraGenerator(int a_nb)
{
   #ifdef CASTOR_VERBOSE
  if(m_verbose >=4) Cout("sRandomNumberGenerator::GetExtraGenerator..."<< endl); 
  #endif
  
  return mp_extraEngines[a_nb];
} 
