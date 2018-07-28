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
  \brief Declaration of class sChronoManager
*/

#ifndef SCHRONOMANAGER_HH
#define SCHRONOMANAGER_HH 1

#include "gVariables.hh"
#include "sOutputManager.hh"

// Typedef for times and durations
typedef std::chrono::time_point<std::chrono::system_clock> ChronoTime;
typedef std::chrono::duration<int64_t,std::nano> DurationNano;
typedef std::chrono::milliseconds Ms;
typedef std::chrono::seconds Secs;
typedef std::chrono::minutes Mins;
typedef std::chrono::hours Hs;

/*!
  \class   sChronoManager
  \brief   This class is designed to manage some profiling of the code
  \details We'll see for details...
*/
class sChronoManager
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public :
    /*!
      \fn      public static sChronoManager* sChronoManager::GetInstance()
      \brief   Instantiate the singleton if not already done, then return the pointer to its instance
      \return  mp_Instance
    */
    static sChronoManager *GetInstance()
    {
      if (mp_Instance == NULL) mp_Instance = new sChronoManager;
      return mp_Instance;
    };
    /*!
      \fn      public sChronoManager::~sChronoManager()
      \brief   The destructor of sChronoManager
      \details This is the default and unique destructor. It simply replace the instance pointer by NULL.
    */
    ~sChronoManager() {mp_Instance = NULL;}

  // -----------------------------------------------------------------------------------------
  // Public functions
  public:
    /*!
      \fn      public int sChronoManager::CheckParameters()
      \brief   Check validity of all parameters
      \return  0 if success, another value otherwise
    */
    int CheckParameters();
    /*!
      \fn      public int sChronoManager::Initialize()
      \brief   Initialize all thread-safe buffers for profiling
      \return  0 if success, another value otherwise
    */
    int Initialize();
    /*!
      \fn      public void sChronoManager::Display()
      \brief   Display the results of the duration buffers
    */
    void Display();

  // -----------------------------------------------------------------------------------------
  // Public inline functions for updating timers
  public:
    /*!
      \fn     inline public void sChronoManager::StartIterativeDataUpdateStep1()
      \param  int a_thread
      \brief  Start the timer for duration of iterative data update step 1
    */
    inline void StartIterativeDataUpdateStep1(int a_thread)
           {mp_startIterativeDataUpdateStep1[a_thread] = std::chrono::system_clock::now();}
    /*!
      \fn     inline public void sChronoManager::StopIterativeDataUpdateStep1()
      \param  int a_thread
      \brief  Stop the timer for duration of iterative data update step 1
    */
    inline void StopIterativeDataUpdateStep1(int a_thread)
           {mp_durationIterativeDataUpdateStep1[a_thread] += std::chrono::system_clock::now() - mp_startIterativeDataUpdateStep1[a_thread];}
    /*!
      \fn     inline public void sChronoManager::StartIterativeDataUpdateStep2()
      \param  int a_thread
      \brief  Start the timer for duration of iterative data update step 2
    */
    inline void StartIterativeDataUpdateStep2(int a_thread)
           {mp_startIterativeDataUpdateStep2[a_thread] = std::chrono::system_clock::now();}
    /*!
      \fn     inline public void sChronoManager::StopIterativeDataUpdateStep2()
      \param  int a_thread
      \brief  Stop the timer for duration of iterative data update step 2
    */
    inline void StopIterativeDataUpdateStep2(int a_thread)
           {mp_durationIterativeDataUpdateStep2[a_thread] += std::chrono::system_clock::now() - mp_startIterativeDataUpdateStep2[a_thread];}
    /*!
      \fn     inline public void sChronoManager::StartIterativeDataUpdateStep3()
      \param  int a_thread
      \brief  Start the timer for duration of iterative data update step 3
    */
    inline void StartIterativeDataUpdateStep3(int a_thread)
           {mp_startIterativeDataUpdateStep3[a_thread] = std::chrono::system_clock::now();}
    /*!
      \fn     inline public void sChronoManager::StopIterativeDataUpdateStep3()
      \param  int a_thread
      \brief  Stop the timer for duration of iterative data update step 3
    */
    inline void StopIterativeDataUpdateStep3(int a_thread)
           {mp_durationIterativeDataUpdateStep3[a_thread] += std::chrono::system_clock::now() - mp_startIterativeDataUpdateStep3[a_thread];}
    /*!
      \fn     inline public void sChronoManager::StartIterativeDataUpdateStep4()
      \param  int a_thread
      \brief  Start the timer for duration of iterative data update step 4
    */
    inline void StartIterativeDataUpdateStep4(int a_thread)
           {mp_startIterativeDataUpdateStep4[a_thread] = std::chrono::system_clock::now();}
    /*!
      \fn     inline public void sChronoManager::StopIterativeDataUpdateStep4()
      \param  int a_thread
      \brief  Stop the timer for duration of iterative data update step 4
    */
    inline void StopIterativeDataUpdateStep4(int a_thread)
           {mp_durationIterativeDataUpdateStep4[a_thread] += std::chrono::system_clock::now() - mp_startIterativeDataUpdateStep4[a_thread];}
    /*!
      \fn     inline public void sChronoManager::StartConvolution()
      \brief  Start the timer for duration of convolution
    */
    inline void StartConvolution()
           {m_startConvolution = std::chrono::system_clock::now();}
    /*!
      \fn     inline public void sChronoManager::StopConvolution()
      \brief  Stop the timer for duration of convolution
    */
    inline void StopConvolution()
           {m_durationConvolution += std::chrono::system_clock::now() - m_startConvolution;}
    /*!
      \fn     inline public void sChronoManager::StartCustomStep()
      \param  int a_thread
      \param  int a_step
      \brief  Start the timer for duration of custom step of the given index for the given thread
    */
    inline void StartCustomStep(int a_thread, int a_step)
           {mpp_startCustomSteps[a_step][a_thread] = std::chrono::system_clock::now();}
    /*!
      \fn     inline public void sChronoManager::StopCustomStep()
      \brief  Stop the timer for duration of custom step of the given index for the given thread
    */
    inline void StopCustomStep(int a_thread, int a_step)
           {mpp_durationCustomSteps[a_step][a_thread] += std::chrono::system_clock::now() - mpp_startCustomSteps[a_step][a_thread];}

  // -----------------------------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline public sChronoManager::SetNbThreads()
      \param   int a_nbThreadsForProjection
      \param   int a_nbThreadsForImageComputation
      \brief   Set the number of threads for both projection and image computations
    */
    inline void SetNbThreads(int a_nbThreadsForProjection, int a_nbThreadsForImageComputation)
           {m_nbThreadsForProjection = a_nbThreadsForProjection; m_nbThreadsForImageComputation = a_nbThreadsForImageComputation;
            m_nbThreadsMax = m_nbThreadsForProjection; if (m_nbThreadsForImageComputation>m_nbThreadsForProjection) m_nbThreadsMax = m_nbThreadsForImageComputation;}
    /*!
      \fn      inline public sChronoManager::SetNbCustomSteps()
      \param   int a_nbCustomSteps
      \brief   Set the number of custom steps for profiling
    */
    inline void SetNbCustomSteps(int a_nbCustomSteps)
           {m_nbCustomSteps = a_nbCustomSteps;}
    /*!
      \fn      inline public sChronoManager::SetVerbose()
      \param   int a_verbose
      \brief   Set the verbose level
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}

  // -------------------------------------------------------------------
  // Private constructor
  private:
    /*!
      \fn      private sChronoManager::sChronoManager()
      \brief   The constructor of sChronoManager
      \details This is the default and unique constructor. It does nothing.
    */
    sChronoManager();
    // Prevent the compiler to generate methods to copy the object
    sChronoManager(sChronoManager const&){};     
    void operator=(sChronoManager const&){};

  // -------------------------------------------------------------------
  // Private data members
  private:
    static sChronoManager *mp_Instance; /*!< Pointer to this singleton object */
    // Number of threads
    int m_nbThreadsForProjection;
    int m_nbThreadsForImageComputation;
    int m_nbThreadsMax;
    // Permanent profiled steps
    DurationNano* mp_durationIterativeDataUpdateStep1;
    DurationNano* mp_durationIterativeDataUpdateStep2;
    DurationNano* mp_durationIterativeDataUpdateStep3;
    DurationNano* mp_durationIterativeDataUpdateStep4;
    ChronoTime* mp_startIterativeDataUpdateStep1;
    ChronoTime* mp_startIterativeDataUpdateStep2;
    ChronoTime* mp_startIterativeDataUpdateStep3;
    ChronoTime* mp_startIterativeDataUpdateStep4;
    DurationNano m_durationConvolution;
    ChronoTime m_startConvolution;
    // Custom profiled steps
    int m_nbCustomSteps;
    DurationNano** mpp_durationCustomSteps;
    ChronoTime** mpp_startCustomSteps;
    // Verbose
    int m_verbose;
};

#endif
