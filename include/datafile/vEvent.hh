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
  \brief    Declaration of class vEvent
*/

#ifndef VEVENT_HH
#define VEVENT_HH 1

#include "gVariables.hh"

/*!
  \class   vEvent
  \brief   Mother class for the Event objects
  \details This class is designed to be a mother virtual class that should not be used on its own; only its children are used. \n
           The pure virtual GetEventIndices method is implemented in each children in order to get the number of lines included \n
           in the event and the associated indices (two crystal indices for PET, one crystal index and one view index for SPECT, etc).
*/
class vEvent
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   vEvent constructor. 
               Initialize the member variables to their default values.
    */
    vEvent();
    /*!
      \brief   vEvent destructor
    */
    virtual ~vEvent();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      int vEvent::AllocateID()
      \brief   Instantiate the mp_ID1 and mp_ID2 indices arrays
      \details This function instantiate the mp_ID1 and mp_ID2 indices arrays using the m_nbLines filled (assuming it has been initialized before),  \n
               and call the AllocateSpecificData() function implemented in child classes
      \return  0 is success, positive value otherwise
    */
    int AllocateID();
    /*!
      \fn      virtual int vEvent::AllocateSpecificData() = 0
      \brief   Pure virtual function implemented in the child classes, dedicated to the allocation of specific data in the child classes
      \return  0 is success, positive value otherwise
    */
    virtual int AllocateSpecificData() = 0;
    /*!
      \fn      virtual void vEvent::Describe() = 0
      \brief   This function can be used to get a description of the event printed out
    */
    virtual void Describe() = 0;

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline uint32_t vEvent::GetTimeInMs()
      \return  the timestamp of the Event
    */
    inline uint32_t GetTimeInMs()
           {return m_timeInMs;}
    /*!
      \fn      inline uint16_t vEvent::GetNbLines()
      \return  the number of lines in the Event
    */
    inline uint16_t GetNbLines()
           {return m_nbLines;};
    /*!
      \fn      inline int vEvent::GetID1()
      \param   a line index (0 if the event contains one line, any number if the events contains several lines as in PET compression)
      \return  the indice of the 1st ID of the Event corresponding to the line index passed as parameter
    */
    inline uint32_t GetID1(int a_line)
           {return mp_ID1[a_line];}
    /*!
      \fn      inline int vEvent::GetID2()
      \param   a line index (0 if the event contains one line, any number if the events contains several lines as in PET compression)
      \return  the indice of the 2nd ID of the Event corresponding to the line index passed as parameter
    */
    inline uint32_t GetID2(int a_line)
           {return mp_ID2[a_line];}
    /*!
      \fn      inline int* vEvent::GetEventID1()
      \return  the pointer containing the indices of the 1st ID of the event
    */
    inline uint32_t* GetEventID1()
           {return mp_ID1;}
    /*!
      \fn      inline int* vEvent::GetEventID2()
      \return  the pointer containing the indices of the 2nd ID of the event
    */
    inline uint32_t* GetEventID2()
           {return mp_ID2;}
    /*!
      \fn      inline int vEvent::GetDataType()
      \return  the type of the event (PET, SPECT, Transmission)
    */
    inline int GetDataType()
           {return m_dataType;}
    /*!
      \fn      inline int vEvent::GetDataMode()
      \return  the mode of the event (PET, SPECT, Transmission)
    */
    inline int GetDataMode()
           {return m_dataMode;}
    /*!
      \fn      inline void vEvent::SetTimeInMs()
      \param   a time value in ms
      \brief   Set the timestamp of the Event.
    */
    inline void SetTimeInMs(uint32_t a_value)
           {m_timeInMs = a_value;};
    /*!
      \fn      virtual void vEvent::SetEventValue() = 0
      \param   a_bin
      \param   a_value
      \brief   This function is implemented by child classes 
      \details Cast the FLTNBDATA value passed in parameters in FLTNB, and use it to set the event value of the specific TOF bin
    */
    virtual void SetEventValue(int a_bin, FLTNBDATA a_value) = 0; 
    /*!
      \fn      inline void vEvent::SetNbLines()
      \param   a number of lines
      \brief   Set the number of lines of the Event.
    */
    inline void SetNbLines(uint16_t a_value)
           {m_nbLines = a_value;};
    /*!
      \fn      inline void vEvent::SetID1()
      \param   a line index
      \param   a value for the ID
      \brief   Set the indice associated with the line index for the 1st ID of the Event
    */
    inline void SetID1(int a_line, uint32_t a_value)
           {mp_ID1[a_line] = a_value;}
    /*!
      \fn      inline void vEvent::SetID2()
      \param   a line index
      \param   a value for the ID
      \brief   Set the indice associated with the line index for the 2nd ID of the Event
    */
    inline void SetID2(int a_line, uint32_t a_value)
           {mp_ID2[a_line] = a_value;}
    /*!
      \fn      inline void vEvent::SetVerbose()
      \param   a verbose level
      \brief   Set verbosity
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      virtual FLTNB vEvent::GetEventValue() = 0
      \param   a bin
      \brief   Pure virtual function implemented in the child classes
      \return  the value of the event
    */
    virtual FLTNB GetEventValue(int a_bin) = 0;
    /*!
      \fn      virtual FLTNB vEvent::GetAdditiveCorrections() = 0
      \brief   Pure virtual function implemented in the child classes
      \return  the sum of the additive corrections terms for this event
    */
    virtual FLTNB GetAdditiveCorrections(int a_bin) = 0;
    /*!
      \fn      virtual FLTNB vEvent::GetBlankValue()
      \brief   This is a pure virtual function implemented in the child classes
      \return  the blank measurement if relevant, 0. otherwise
    */
    inline virtual FLTNB GetBlankValue()
           {return 0.;}
    /*!
      \fn      virtual FLTNB vEvent::GetMultiplicativeCorrections() = 0
      \brief   This is a pure virtual function implemented in the child classes
      \return  the product of the multiplicative corrections terms for this event
    */
    virtual FLTNB GetMultiplicativeCorrections() = 0;
    /*!
      \fn      virtual INTNB vEvent::GetNbValueBins() = 0
      \brief   Get the number of event value bins
      \return  the number of value bins
    */
    virtual INTNB GetNbValueBins() = 0;
    /*!
      \fn      virtual void vEvent::MultiplyAdditiveCorrections() = 0
      \param   FLTNB a_factor
      \brief   Divide additive corrections by the provided factor
    */
    virtual void MultiplyAdditiveCorrections(FLTNB a_factor) = 0;

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  protected:
    uint32_t m_timeInMs; /*!< Timestamp of the event in ms. Default value =0 */
    uint16_t m_nbLines;  /*!< Number of lines in the event. Default value =0 */
    uint32_t* mp_ID1;    /*!< Pointer containing the indice(s) of the 1st ID of the Event. Default value =NULL  */
    uint32_t* mp_ID2;    /*!< Pointer containing the indice(s) of the 2nd ID of the Event. Default value =0 */
    FLTNB m_eventValue;  /*!< Amount of data in the bin. Default value =0.0 (or 1. for list-mode event) */
    int m_dataType;      /*!< This integer is used to specify the type of the event as in the vDataFile; PET, SPECT or CT. Default value =Unknown */
    int m_dataMode;      /*!< This integer is used to specify the mode of the event as in the vDataFile; LIST or HISTO. Default value =Unknown */
    int m_dataSpec;                           /*!< Flag indicating the physical specificity of the data: SPEC_EMISSION or SPEC_TRANSMISSION */
    int m_verbose;       /*!< Verbosity. Default value =-1 */
};

#endif
