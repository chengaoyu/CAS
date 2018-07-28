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
  \ingroup  management
  \brief    This header file is mainly used to declare some macro definitions and all includes needed from the standard c library
  \details  It is not associated with any implementation file as it does not define any function. \n
            This file should be included at the beginning of any files in this project.
*/

#ifndef GVARIABLES_HH
#define GVARIABLES_HH 1

// All includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <string>
#include <cmath>
#include <ctime>
#include <map>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <random>
#include <cfloat>
#include <chrono>
#include <iomanip>
#ifdef CASTOR_OMP
#include <omp.h>
#endif
#ifdef CASTOR_MPI
#include <mpi.h>
#endif
#ifdef _WIN32
#include <windows.h>
#endif

// The std namespace is used throughout the whole project
using namespace std;

// The current CASToR version
#define CASTOR_VERSION "2.0.1"

/**
 * @defgroup PROGRAM_PRECISION The precision of the computation and I/O operations
 *
 *    \brief Precision of the computation and I/O operations \n
 *           Defined in gVariables.hh
 * @{
 */

/** Define the precision for main program implementation (matrices, most computations, output images) */
#define FLTNB     float
/** Define high precision for computations where at least the double precision is required */
#define HPFLTNB   double
/** Define the precision of the program implementation (for matrices and computations), specific to MPI operations */
#define FLTNBMPI  MPI_FLOAT
/** Define the precision of the input/output datafile read/write */
#define FLTNBDATA float
/** Define the precision of the input/output scanner LUT read/write */
#define FLTNBLUT  float
/** Define the integer type for matrices dimensions */
//#define INTNB     int64_t
#define INTNB     int

/** @} */

/** Define the exit value for DEBUG use */
#define EXIT_DEBUG 10

/** Define the value of the standard 2*sqrt(2*ln(2)) */
#define TWO_SQRT_TWO_LN_2 2.354820045

/** Define the value of 1/sqrt(2*PI) */
#define INV_SQRT_2_PI 0.398942280

/** Define the value of the speed of light in mm/ps */
#define SPEED_OF_LIGHT 0.299792458

#ifdef _WIN32
#define M_PI 3.141592654
#endif

#endif
