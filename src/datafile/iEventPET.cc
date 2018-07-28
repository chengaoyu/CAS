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

  \brief Implementation of class iEventPET
*/

#include "iEventPET.hh"
#include "vDataFile.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventPET::iEventPET() : vEvent() 
{
  m_dataType = TYPE_PET;
  m_dataSpec = SPEC_EMISSION;
  m_eventRdmRate = 0.;
  m_atnCorrFactor = 1.;
  m_eventNormFactor = 1.;
  m_nbLines = 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventPET::~iEventPET() {}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventPET::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("iEventPET::Describe() -> Display contents" << endl);
  Cout("Time: " << m_timeInMs << " ms" << endl);
  Cout("Number of lines: " << m_nbLines << endl);
  for (uint16_t l=0; l<m_nbLines; l++) Cout("  --> ID1: " << mp_ID1[l] << " | ID2: " << mp_ID2[l] << endl);
  Cout("Random rate: " << m_eventRdmRate << endl);
  Cout("Normalization factor: " << m_eventNormFactor << endl);
  Cout("ACF: " << m_atnCorrFactor << endl);
  Cout(flush);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
