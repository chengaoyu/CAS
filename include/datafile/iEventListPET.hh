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
  \brief    Declaration of class iEventListPET
*/

#ifndef IEVENTLISTPET_HH
#define IEVENTLISTPET_HH 1

#include "iEventPET.hh"

/**
 * @defgroup EVENT_KIND Event kind
 *
 *    \brief Type of event in list-mode PET (true, scatter, random, ...)  \n
 *           Defined in iEventPET.hh
 * 
 * @{
 */
/** Constant corresponding to an event of unknown kind (=0) */
#define KIND_UNKNOWN 0
/** Constant corresponding to a true event (=1) */
#define KIND_TRUE 1
/** Constant corresponding to a scatter event (=2) */
#define KIND_SCAT 2
/** Constant corresponding to a multiple scatter event (=3) */
#define KIND_MSCAT 3
/** Constant corresponding to a random event (=4) */
#define KIND_RDM 4
/** @} */



/*!
  \class   iEventListPET
  \brief   Inherit from iEventPET. Class for PET list-mode events
  \details It manages data and functions specific to list mode PET.
*/
class iEventListPET : public iEventPET
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventListPET constructor. 
               Initialize the member variables to their default values.
    */
    iEventListPET();    
    /*!
      \brief   iEventListPET destructor. 
    */
    ~iEventListPET();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      inline int iEventListPET::AllocateSpecificData()
      \brief   Function allowing the allocation of specific data. Return 0 by default for iEventListPET
      \return  0 is success, positive value otherwise
    */
    inline int AllocateSpecificData() {return 0;}
    /*!
      \fn      void iEventListPET::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline uint8_t iEventListPET::GetKind() 
      \return  the kind of coincidence
    */
    inline uint8_t GetKind()
           {return m_kind;}
    /*!
      \fn      inline FLTNB iEventListPET::GetEventScatRate()
      \param   a_bin (0 if noTOF)
      \return  the scatter correction rate in 1/s for this list-mode event (a_bin is dedicated to histogram mode, and ignored for list-mode)
    */
    inline FLTNB GetEventScatRate(int a_bin)
           {return m_eventScatRate;}
    /*!
      \fn      inline uint32_t iEventListPET::GetTOFMeasurement()
      \return  the TOF measurement value
    */
    inline FLTNB GetTOFMeasurement()
           {return m_TOFMeasurement;}
    /*!
      \fn      inline FLTNB* iEventListPET::GetPOI1()
      \return  the pointer containing the 3 indices for the Point Of Interaction of the 1st crystal
    */
    inline FLTNB* GetPOI1()
           {return mp_POI1;}
    /*!
      \fn      inline FLTNB* iEventListPET::GetPOI2()
      \return  the pointer containing the 3 indices for the Point Of Interaction of the 2nd crystal
    */
    inline FLTNB* GetPOI2()
           {return mp_POI2;}
    /*!
      \fn      inline FLTNB iEventListPET::GetPOI1()
      \param   axis
      \return  the Point of interaction of the 1st crystal along the axis specified in parameter
    */
    inline FLTNB GetPOI1(uint8_t axis)
           {return mp_POI1[axis];}
    /*!
      \fn      inline FLTNB iEventListPET::GetPOI2()
      \param   axis
      \return  the Point of interaction of the 2nd crystal along the axis specified in parameter
    */
    inline FLTNB GetPOI2(uint8_t axis)
           {return mp_POI2[axis];}
    /*!
      \fn      inline FLTNB iEventListPET::GetEventValue()
      \param   a bin (0 if noTOF)
      \return  1 as default for a list-mode Event (a_bin is dedicated to histogram mode, and ignored for list-mode)
    */
    inline FLTNB GetEventValue(int a_bin)
           {return 1.;}
    /*!
      \fn      inline void iEventListPET::SetKind()
      \param   a_value
      \brief   Set the kind of coincidence
    */
    inline void SetKind(uint8_t a_value)
           {m_kind = a_value;}
    /*!
      \fn      inline void iEventListPET::SetScatterRate()
      \param   a bin (no use for list-mode)
      \param   a_value
      \brief   Cast the FLTNBDATA value passed as parameter in FLTNB, and set it to the scatter correction rate in 1/s
               (a_bin is dedicated to histogram mode, and ignored for list-mode)
    */
    inline void SetScatterRate(int a_bin, FLTNBDATA a_value)
           {m_eventScatRate = (FLTNB)a_value;}
    /*!
      \fn      inline void iEventListPET::SetTOFMeasurement()
      \param   a_value
      \brief   Initialize the TOFmeasurement with a value passed in parameters
    */
    inline void SetTOFMeasurement(FLTNB a_value)
           {m_TOFMeasurement = a_value;}
    /*!
      \fn      inline void iEventListPET::SetSpatialTOFMeasurementRange()
      \param   a_value
      \brief   Initialize the spatialTOFMeasurementRange with a value passed in parameters
    */
    inline void SetSpatialTOFMeasurementRange(FLTNB a_value)
           {m_spatialTOFMeasurementRange = a_value;}
    /*!
      \fn      inline void iEventListPET::SetPOI1()
      \param   a_axis
      \param   a_value
      \brief   Initialize the POI of the crystal #1 with a value for the specific axis.
    */
    inline void SetPOI1(int a_axis, FLTNBDATA a_value)
           {mp_POI1[a_axis] = (FLTNB)a_value;}
    /*!
      \fn      inline void iEventListPET::SetPOI2()
      \param   a_axis
      \param   a_value
      \brief   Initialize the POI of the crystal #2 with a value for the specific axis.
    */
    inline void SetPOI2(int a_axis, FLTNBDATA a_value)
           {mp_POI2[a_axis] = (FLTNB)a_value;}
    /*!
      \fn      void iEventListPET::SetEventValue()
      \param   a_bin (0 if noTOF)
      \param   a_value
      \brief   Throw a warning (depending of verbosity) as the event value of a list-mode Event should be equal to 1 and not modified ?
    */
    void SetEventValue(int a_bin, FLTNBDATA a_value); 
    /*!
      \fn      FLTNB iEventListPET::GetAdditiveCorrections()
      \param   a bin (0 if noTOF)
      \return  The sum of additive correction terms (the TOF bin parameter is ignored for list-mode PET)
    */
    FLTNB GetAdditiveCorrections(int a_bin);
    /*!
      \fn      virtual INTNB vEvent::GetNbValueBins()
      \brief   Get the number of event value bins
      \return  Single value so 1 bin
    */
    inline INTNB GetNbValueBins()
           {return 1;}
    /*!
      \fn      void iEventListPET::SetHasTOFinfo()
      \brief   Set whether this event contains TOF information or not, 
               awareness of TOF info existence is important for some operations
      \return  None
    */
    void SetHasTOFinfo(bool a_hasTOFinfo)
           {m_hasTOFinfo = a_hasTOFinfo;}
    /*!
      \fn      inline void iEventListPET::MultiplyAdditiveCorrections()
      \param   FLTNB a_factor
      \brief   Divide additive corrections by the provided factor (scatters and randoms)
    */
    inline void MultiplyAdditiveCorrections(FLTNB a_factor)
           {m_eventRdmRate *= a_factor; m_eventScatRate *= a_factor;}

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  private:
    uint8_t m_kind;            /*!< Coincidence type : unknown (=0), true(=1), single scat(=2), multiple scat(=3), random(=4)) Default value =KIND_UNKNOWN  */
    FLTNB m_eventScatRate;     /*!< Correction term for scatter rate (unit: s-1). Default value =0.0 */
    FLTNB m_TOFMeasurement;    /*!< TOF time in ps, corresponding to the difference in arrival time between the two photons. Default value =0.0 */
    FLTNB mp_POI1[3];          /*!< Position of interaction in the crystal #1 along each axis (mm). Default value =0.0;0.0;-1.0 */
    FLTNB mp_POI2[3];          /*!< Position of interaction in the crystal #2 along each axis (mm). Default value =0.0;0.0;-1.0 */
    FLTNB m_spatialTOFMeasurementRange;  /*!< Spatial equivalent of the maximum range of values for TOF delta time measurements, in mm */
    bool m_hasTOFinfo;         /*!< Flag indicating whether TOF information exists or not */
};

#endif
