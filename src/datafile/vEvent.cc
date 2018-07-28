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
  \ingroup datafile

  \brief Implementation of class vEvent
*/

#include "vEvent.hh"
#include "vDataFile.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vEvent::vEvent()
{
  m_timeInMs = 0;
  m_nbLines = 0;
  mp_ID1 = NULL;
  mp_ID2 = NULL;
  m_dataType = TYPE_UNKNOWN;
  m_dataMode = MODE_UNKNOWN;
  m_verbose = -1;
  m_eventValue = 0.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vEvent::~vEvent() 
{
  if (mp_ID1 != NULL) free(mp_ID1);
  if (mp_ID1 != NULL) free(mp_ID2);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vEvent::AllocateID() 
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("vEvent::AllocateID() -> Allocate buffers for indices" << endl);
  // Check number of lines
  if (m_nbLines<1)
  {
    Cerr("***** vEvent::AllocateID() -> Error, number of lines has not been initialized (<1) !" << endl);
    return 1;
  }
  else
  {
    // Allocate buffer indices
    mp_ID1 = (uint32_t*)malloc(m_nbLines*sizeof(uint32_t));
    mp_ID2 = (uint32_t*)malloc(m_nbLines*sizeof(uint32_t));
  }
  // Call the pure virtual function implemented in child classes for the allocation of data specific to child classes
  if (AllocateSpecificData())
  {
    Cerr("***** vEvent::AllocateID() -> Error when trying to allocated specific data for the Event !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
