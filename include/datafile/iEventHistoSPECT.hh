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
  \ingroup  datafile
  \brief    Declaration of class iEventHistoSPECT
*/

#ifndef IEVENTHISTOSPECT_HH
#define IEVENTHISTOSPECT_HH 1

#include "iEventSPECT.hh"

/*!
  \class   iEventHistoSPECT
  \brief   Inherit from iEventSPECT. Class for SPECT histogram mode events 
  \details It manages data and functions specific to histo mode SPECT.
*/
class iEventHistoSPECT : public iEventSPECT
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventHistoSPECT constructor. 
               Initialize the member variables to their default values.
    */
    iEventHistoSPECT();
    /*!
      \brief   iEventHistoSPECT destructor. 
    */
    ~iEventHistoSPECT();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      inline int iEventHistoSPECT::AllocateSpecificData()
      \brief   Function allowing the allocation of specific data. Return 0 by default for iEventHistoSPECT
      \return  0 is success, positive value otherwise
    */
    inline int AllocateSpecificData()
           {return 0;}
    /*!
      \fn      void iEventHistoSPECT::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline FLTNB iEventHistoSPECT::GetEventValue()
      \param   a_bin (0 if noTOF)
      \return  the event value (a_bin is dedicated to PET and ignored for SPECT)
    */
    inline FLTNB GetEventValue(int a_bin)
           {return m_eventValue;}
    /*!
      \fn      inline void iEventHistoSPECT::SetEventValue()
      \param   a_bin (0 if noTOF)
      \param   a_value
      \brief   Cast the FLTNBDATA value passed in parameters in FLTNB, and use it to set the event value (a_bin is dedicated to PET and ignored for SPECT)
    */
    inline void SetEventValue(int a_bin, FLTNBDATA a_value)
           {m_eventValue = (FLTNB)a_value;}
    /*!
      \fn      virtual INTNB vEvent::GetNbValueBins()
      \brief   Get the number of event value bins
      \return  Single value so 1 bin
    */
    inline INTNB GetNbValueBins()
           {return 1;}

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  private:
};

#endif
