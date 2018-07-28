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
  \brief   Implementation of file to memory mapping
  \details This class was taken from Stephan Brumme website and
           adapted to the CASToR project with almost no modifications.
           Many thanks to him!
           Copyright (c) 2013 Stephan Brumme. All rights reserved.
           see http://create.stephan-brumme.com/disclaimer.html
*/

#include "oMemoryMapped.hh"

#include <stdexcept>
#include <cstdio>
#include <iostream>

// OS-specific
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
// Windows
#include <windows.h>
#else
// Linux
// enable large file support on 32 bit systems
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#ifdef  _FILE_OFFSET_BITS
#undef  _FILE_OFFSET_BITS
#endif
#define _FILE_OFFSET_BITS 64
// and include needed headers
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#endif

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// do nothing, must use open()
oMemoryMapped::oMemoryMapped()
: _filename   (),
  _filesize   (0),
  _hint       (Normal),
  _mappedBytes(0),
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
  _mappedFile (NULL),
#endif
  _file       (0),
  _mappedView (NULL)
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// open file, mappedBytes = 0 maps the whole file
oMemoryMapped::oMemoryMapped(const std::string& filename, size_t mappedBytes, CacheHint hint)
: _filename   (filename),
  _filesize   (0),
  _hint       (hint),
  _mappedBytes(mappedBytes),
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
  _mappedFile (NULL),
#endif
  _file       (0),
  _mappedView (NULL)
{
  Open(filename, mappedBytes, hint);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// close file (see close() )
oMemoryMapped::~oMemoryMapped()
{
  Close();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// open file
int oMemoryMapped::Open(const std::string& filename, size_t mappedBytes, CacheHint hint)
{
  // already open ?
  if (IsValid())
  {
    Cerr("***** oMemoryMapped::Open() -> File is already open !" << endl);
    return 1;
  }

  _file       = 0;
  _filesize   = 0;
  _hint       = hint;
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
  _mappedFile = NULL;
#endif
  _mappedView = NULL;

#if defined(_WIN32) || defined(CASTOR_USE_MINGW)

  // ===================================================================
  // Windows
  // ===================================================================

  DWORD winHint = 0;
  switch (_hint)
  {
    case Normal:         winHint = FILE_ATTRIBUTE_NORMAL;     break;
    case SequentialScan: winHint = FILE_FLAG_SEQUENTIAL_SCAN; break;
    case RandomAccess:   winHint = FILE_FLAG_RANDOM_ACCESS;   break;
    default: break;
  }

  // open file
  _file = ::CreateFileA(filename.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, winHint, NULL);
  if (!_file)
  {
    Cerr("***** oMemoryMapped::Open() -> Failed to create windows file from function CreateFileA() !" << endl);
    return 1;
  }

  // file size
  LARGE_INTEGER result;
  if (!GetFileSizeEx(_file, &result))
  {
    Cerr("***** oMemoryMapped::Open() -> Failed to get file size from windows function GetFileSizeEx() !" << endl);
    return 1;
  }
  _filesize = static_cast<uint64_t>(result.QuadPart);

  // convert to mapped mode
  _mappedFile = ::CreateFileMapping(_file, NULL, PAGE_READONLY, 0, 0, NULL);
  if (!_mappedFile)
  {
    Cerr("***** oMemoryMapped::Open() -> Failed to convert file to mapped mode from windows function CreateFileMapping() !" << endl);
    return 1;
  }

#else

  // ===================================================================
  // Linux
  // ===================================================================

  // open file
  //_file = ::open(filename.c_str(), O_RDONLY | O_LARGEFILE);
  _file = ::open(filename.c_str(), O_RDONLY);
  
  if (_file == -1)
  {
    _file = 0;
    Cerr("***** oMemoryMapped::Open() -> Failed to open file from unix function open() !" << endl);
    return 1;
  }

  // file size
  struct stat statInfo;
  if (fstat(_file, &statInfo) < 0)
  {
    Cerr("***** oMemoryMapped::Open() -> Failed to get correct file size from unix function fstat() !" << endl);
    return 1;
  }

  _filesize = statInfo.st_size;
#endif

  // initial mapping
  if (Remap(0, mappedBytes))
  {
    Cerr("***** oMemoryMapped::Open() -> A problem occured while calling the oMemoryMapped::Remap() function !" << endl);
    return 1;
  }

  // check
  if (!_mappedView)
  {
    Cerr("***** oMemoryMapped::Open() -> Failed to get a correct mapped view after calling the oMemoryMapped::Remap() function !" << endl);
    return 1;
  }

  // everything's fine
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// close file
void oMemoryMapped::Close()
{
  // kill pointer
  if (_mappedView)
  {
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
    ::UnmapViewOfFile(_mappedView);
#else
    ::munmap(_mappedView, _filesize);
#endif
    _mappedView = NULL;
  }

#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
  if (_mappedFile)
  {
    ::CloseHandle(_mappedFile);
    _mappedFile = NULL;
  }
#endif

  // close underlying file
  if (_file)
  {
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
    ::CloseHandle(_file);
#else
    ::close(_file);
#endif
    _file = 0;
  }

  _filesize = 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// access position, no range checking (faster)
unsigned char oMemoryMapped::operator[](size_t offset) const
{
  return ((unsigned char*)_mappedView)[offset];
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// access position, including range checking
unsigned char oMemoryMapped::at(size_t offset) const
{
  // checks
  if (!_mappedView)
    throw std::invalid_argument("No view mapped");
  if (offset >= _filesize)
    throw std::out_of_range("View is not large enough");

  return operator[](offset);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// raw access
const unsigned char* oMemoryMapped::GetData() const
{
  return (const unsigned char*)_mappedView;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// true, if file successfully opened
bool oMemoryMapped::IsValid() const
{
  return _mappedView != NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// get file size
uint64_t oMemoryMapped::size() const
{
  return _filesize;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// get number of actually mapped bytes
size_t oMemoryMapped::mappedSize() const
{
  return _mappedBytes;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// replace mapping by a new one of the same file, offset MUST be a multiple of the page size
int oMemoryMapped::Remap(uint64_t offset, size_t mappedBytes)
{
  if (!_file)
  {
    Cerr("***** oMemoryMapped::Remap() -> Cannot remap a file that has not been created !" << endl);
    return 1;
  }

  if (mappedBytes == WholeFile)
    mappedBytes = _filesize;

  // close old mapping
  if (_mappedView)
  {
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
    ::UnmapViewOfFile(_mappedView);
#else
    ::munmap(_mappedView, _mappedBytes);
#endif
    _mappedView = NULL;
  }

  // don't go further than end of file
  if (offset > _filesize)
  {
    Cerr("***** oMemoryMapped::Remap() -> Provided offset is after the end of file !" << endl);
    return 1;
  }
  if (offset + mappedBytes > _filesize)
    mappedBytes = size_t(_filesize - offset);

#if defined(_WIN32) || defined(CASTOR_USE_MINGW)

  // ===================================================================
  // Windows
  // ===================================================================

  DWORD offsetLow  = DWORD(offset & 0xFFFFFFFF);
  DWORD offsetHigh = DWORD(offset >> 32);
  _mappedBytes = mappedBytes;

  // get memory address
  _mappedView = ::MapViewOfFile(_mappedFile, FILE_MAP_READ, offsetHigh, offsetLow, mappedBytes);

  if (_mappedView == NULL)
  {
    _mappedBytes = 0;
    _mappedView  = NULL;
    Cerr("***** oMemoryMapped::Remap() -> Mapped view is null after calling windows function MapViewOfFile() !" << endl);
    return 1;
  }

#else

  // ===================================================================
  // Linux
  // ===================================================================

  // new mapping
  //_mappedView = ::mmap64(NULL, mappedBytes, PROT_READ, MAP_SHARED, _file, offset);
  _mappedView = ::mmap(NULL, mappedBytes, PROT_READ, MAP_SHARED, _file, offset);
  
  if (_mappedView == MAP_FAILED)
  {
    _mappedBytes = 0;
    _mappedView  = NULL;
    Cerr("***** oMemoryMapped::Remap() -> Mapping failed after calling the unix function mmap64() !" << endl);
    return 1;
  }

  _mappedBytes = mappedBytes;

  // tweak performance
  int linuxHint = 0;
  switch (_hint)
  {
    case Normal:         linuxHint = MADV_NORMAL;     break;
    case SequentialScan: linuxHint = MADV_SEQUENTIAL; break;
    case RandomAccess:   linuxHint = MADV_RANDOM;     break;
    default: break;
  }
  // assume that file will be accessed soon
  //linuxHint |= MADV_WILLNEED;
  // assume that file will be large
  //linuxHint |= MADV_HUGEPAGE;
//  linuxHint [= MADV_NOHUGEPAGE;

  ::madvise(_mappedView, _mappedBytes, linuxHint);

#endif

  // end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/// get OS page size (for remap)
int oMemoryMapped::GetPageSize()
{
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
  SYSTEM_INFO sysInfo;
  GetSystemInfo(&sysInfo);
  return sysInfo.dwAllocationGranularity;
#else
  return sysconf(_SC_PAGESIZE); //::getpagesize();
#endif
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
