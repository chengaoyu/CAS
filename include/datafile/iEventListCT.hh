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
  \brief    Declaration of class iEventListCT
*/

#ifndef IEVENTLISTCT_HH
#define IEVENTLISTCT_HH 1

#include "iEventCT.hh"

/*!
  \class   iEventListCT
  \brief   Inherit from iEventCT. Class for CT list-mode events 
  \details It manages data and functions specific to list mode CT.
*/
class iEventListCT : public iEventCT
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventListCT constructor. 
               Initialize the member variables to their default values.
    */
    iEventListCT();
    /*!
      \brief   iEventListCT destructor. 
    */
    ~iEventListCT();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      inline int iEventListCT::AllocateSpecificData()
      \brief   Function allowing the allocation of specific data. Return 0 by default for iEventListCT
      \return  0 is success, positive value otherwise
    */
    inline int AllocateSpecificData()
           {return 0;}
    /*!
      \fn      void iEventListCT::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline uint8_t iEventListCT::GetKind() 
      \return  the kind of coincidence
    */
    inline uint8_t GetKind()
           {return m_kind;}
    /*!
      \fn      inline FLTNB iEventListCT::GetEventValue()
      \return  1 as default for a list-mode Event (a_bin is dedicated to histogram mode, and ignored for list-mode)
    */
    inline FLTNB GetEventValue(int a_bin)
           {return 1.;}
    /*!
      \fn      inline void iEventListCT::SetKind()
      \param   a_value
      \brief   Set the kind of event
    */
    inline void SetKind(uint8_t a_value)
           {m_kind = a_value;}
    /*!
      \fn      void iEventListCT::SetEventValue()
      \param   a_bin
      \param   a_value
      \brief   Throw a warning (depending of verbosity) as the event value of a list-mode Event should be equal to 1 and not modified
    */
    void SetEventValue(int a_bin, FLTNBDATA a_value);
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
    uint8_t m_kind;  /*!< Event type : unknown (=0), true(=1), single scat(=2), multiple scat(=3)) Default value =0  */
};

#endif
