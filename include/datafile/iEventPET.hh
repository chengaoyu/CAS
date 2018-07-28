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
  \brief    Declaration of class iEventPET
*/

#ifndef IEVENTPET_HH
#define IEVENTPET_HH 1

#include "vEvent.hh"

/*!
  \class   iEventPET
  \brief   Inherit from vEvent. Main PET class for the Event objects
  \details This class is designed to be an abstract class that should not be used on its own; only its children are used. \n
           It manages data and functions common to both class decicated to histogram and list mode PET.
*/
class iEventPET : public vEvent
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventPET constructor. 
               Initialize the member variables to their default values.
    */
    iEventPET();
    /*!
      \brief   iEventPET destructor
    */
    virtual ~iEventPET();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      virtual int iEventPET::AllocateSpecificData() = 0
      \brief   Pure virtual function implemented in the child classes, dedicated to the allocation of specific data in the child classes
      \return  0 is success, positive value otherwise
    */
    virtual int AllocateSpecificData() = 0;
    /*!
      \fn      virtual void iEventPET::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    virtual void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline FLTNB iEventPET::GetEventRdmRate()
      \return  the correction term for randoms given as a rate in 1/s
    */
    inline FLTNB GetEventRdmRate()
           {return m_eventRdmRate;}
    /*!
      \fn      inline FLTNB iEventPET::GetNormFactor()
      \return  the normalization term
    */
    inline FLTNB GetNormFactor()
           {return m_eventNormFactor;}
    /*!
      \fn      inline FLTNB iEventPET::GetAtnCorrFactor()
      \return  the correction term for attenuation
    */
    inline FLTNB GetAtnCorrFactor()
           {return m_atnCorrFactor;}
    /*!
      \fn      inline void iEventPET::SetRandomRate()
      \param   a random rate value in 1/s
      \brief   Cast the FLTNBDATA value passed as parameter in FLTNB, and set it to the random correction term
    */
    inline void SetRandomRate(FLTNBDATA a_value)
           {m_eventRdmRate = (FLTNB)a_value;}
    /*!
      \fn      inline void iEventPET::SetNormalizationFactor()
      \param   a normalization term
      \brief   Cast the FLTNBDATA value passed as parameter in FLTNB, and set it to the normalization term
    */
    inline void SetNormalizationFactor(FLTNBDATA a_value)
           {m_eventNormFactor = (FLTNB)a_value;}
    /*!
      \fn      inline void iEventPET::SetAttenuationCorrectionFactor()
      \param   an attenuation correction factor
      \brief   Cast the FLTNBDATA value passed as parameter in FLTNB, and set it to the attenuation correction term
    */
    inline void SetAttenuationCorrectionFactor(FLTNBDATA a_value)
           {m_atnCorrFactor = (FLTNB)a_value;}
    /*!
      \fn      inline FLTNB iEventPET::GetMultiplicativeCorrections()
      \return  the product of the multiplicative corrections terms for this event
    */
    inline FLTNB GetMultiplicativeCorrections()
           {return m_atnCorrFactor * m_eventNormFactor;}
    /*!
      \fn      virtual FLTNB iEventPET::GetEventValue() = 0
      \param   a bin (0 if only one line)
      \brief   Pure virtual function implemented in the child classes
      \return  the value of the event
    */
    virtual FLTNB GetEventValue(int a_bin) = 0;
    /*!
      \fn      virtual FLTNB iEventPET::GetAdditiveCorrections() = 0
      \param   a bin (0 if only one line)
      \brief   Pure virtual function implemented in the child classes
      \return  the sum of the additive corrections terms for this event
    */
    virtual FLTNB GetAdditiveCorrections(int a_bin) = 0;
    /*!
      \fn      virtual void iEventPET::SetEventValue() = 0
      \param   a bin (0 if only one line)
      \param   a_value
      \brief   Set the event value, this is a pure virtual function implemented in the child classes
    */
    virtual void SetEventValue(int a_bin, FLTNBDATA a_value) = 0;
    /*!
      \fn      virtual void iEventPET::SetScatterRate() = 0
      \param   a bin (0 if only one line)
      \param   a_value
      \brief   Set the scatter rate correction term in 1/s, this is a pure virtual function implemented in the child classes
    */
    virtual void SetScatterRate(int a_bin, FLTNBDATA a_value) = 0;
    /*!
      \fn      virtual INTNB vEvent::GetNbValueBins() = 0
      \brief   Get the number of event value bins
      \return  the number of value bins
    */
    virtual INTNB GetNbValueBins() = 0;

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  protected:
    FLTNB m_eventRdmRate;    /*!< Correction term for randoms rate (unit: s-1). Default value =0.0 */
    FLTNB m_eventNormFactor; /*!< Normalization term. Default value =1.0 */
    FLTNB m_atnCorrFactor;   /*!< Correction term for attenuation. Default value =1.0 */
};

#endif
