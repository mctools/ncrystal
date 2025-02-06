#ifndef NCrystal_Exception_hh
#define NCrystal_Exception_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ncrystal_api_h
#  include "NCrystal/ncapi.h"
#endif

#include <stdexcept>
#include <sstream>
#include <string>

/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// All Exceptions from NCrystal code will be of type NCrystal::Error::ErrType, where       //
// ErrType will be one of:                                                                 //
//                                                                                         //
//   FileNotFound  : Could not find specified input file                                   //
//   DataLoadError : Troubles while loading data from input file                           //
//   MissingInfo   : Code does not have required info for requested operation              //
//   CalcError     : Error encountered during calculations (inconsistent result, NaN, ...) //
//   LogicError    : Implementation errors.                                                //
//   BadInput      : Invalid input to function or method                                   //
//                                                                                         //
// All error types inherit from NCrystal::Error::Exception, which again inherits from      //
// std::runtime_exception. std::runtime_exception::what() provides the error message, and  //
// the file and line number of where the exception was generated can be queried with       //
// getFile() and getLineNo() respectively.                                                 //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace Error {

    class NCRYSTAL_API Exception : public std::runtime_error {
    public:
      explicit Exception(const std::string& msg, const char * f, unsigned l) noexcept;
      explicit Exception(const char * msg,  const char * f, unsigned l) noexcept;
      virtual const char * getTypeName() const noexcept = 0;
      const char * getFile() const noexcept { return m_file; }
      unsigned getLineNo() const noexcept { return m_lineno; }
      Exception( const Exception & o ) noexcept;
      Exception & operator= ( const Exception & o ) noexcept;
      virtual ~Exception() noexcept;
    private:
      const char * m_file;
      unsigned m_lineno;
    };

#define NCRYSTAL_ADD_ERROR_TYPE(ErrType)                                                               \
    struct NCRYSTAL_API ErrType : public Exception {                                                                \
      explicit ErrType(const std::string& msg, const char * f, unsigned l) noexcept : Exception(msg,f,l) {} \
      explicit ErrType(const char * msg,  const char * f, unsigned l) noexcept : Exception(msg,f,l) {}      \
      virtual const char * getTypeName() const noexcept { return #ErrType; }                            \
      virtual ~ErrType() noexcept;                                                                      \
    }
    //List of error types (destructors implemented in .cc):
    NCRYSTAL_ADD_ERROR_TYPE(FileNotFound);
    NCRYSTAL_ADD_ERROR_TYPE(DataLoadError);
    NCRYSTAL_ADD_ERROR_TYPE(MissingInfo);
    NCRYSTAL_ADD_ERROR_TYPE(CalcError);
    NCRYSTAL_ADD_ERROR_TYPE(LogicError);
    NCRYSTAL_ADD_ERROR_TYPE(BadInput);
#undef NCRYSTAL_ADD_ERROR_TYPE
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// Macro's for easy and consistent throwing from within NCrystal code (don't use THROW2    //
// in tight CPU-critical code):                                                            //
//                                                                                         //
//   NNCRYSTAL_THROW(ErrType,"some hardcoded message")                                     //
//   NNCRYSTAL_THROW2(ErrType,"some "<<flexible<<" message")                               //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

#ifdef NCRYSTAL_THROW
#  undef NCRYSTAL_THROW
#endif
#ifdef NCRYSTAL_THROW2
#  undef NCRYSTAL_THROW2
#endif

#define NCRYSTAL_THROW(ErrType, msg)                            \
  {                                                             \
    throw ::NCrystal::Error::ErrType(msg,__FILE__,__LINE__);    \
  }

#define NCRYSTAL_THROW2(ErrType, msg)     \
  {                                       \
    std::ostringstream err_oss;           \
    err_oss << msg;                       \
    NCRYSTAL_THROW(ErrType,err_oss.str()) \
  }

//Custom nc_assert which throws LogicErrors in dbg builds; nc_assert2
//avoids unused variable warnings in opt builds; nc_assert_always is
//enabled in all builds. Note that since these assert's throw, they
//should not be used in destructors.

#ifdef ncrystal_str
#  undef ncrystal_str
#endif
#ifdef ncrystal_xstr
#  undef ncrystal_xstr
#endif
#ifdef nc_assert_always
#  undef nc_assert_always
#endif
#ifdef nc_assert
#  undef nc_assert
#endif
#ifdef nc_assert2
#  undef nc_assert2
#endif
#ifdef nc_assert_rv
#  undef nc_assert_rv
#endif
#ifdef nc_assert_always_rv
#  undef nc_assert_always_rv
#endif

#define ncrystal_str(s) #s
#define ncrystal_xstr(s) ncrystal_str(s)
#define nc_assert_always(x) do { if (!(x)) { NCRYSTAL_THROW(LogicError,\
                              "Assertion failure: " ncrystal_xstr(x)); } } while(0)

#ifndef NDEBUG
#  define nc_assert(x) nc_assert_always(x)
#  define nc_assert2(x) nc_assert_always(x)
#else
#  define nc_assert(x) do {} while(0)
#  define nc_assert2(x) do { (void)sizeof(x); } while(0)//use but dont evaluate x
#endif

namespace NCRYSTAL_NAMESPACE {
  namespace assert_details {
    //For C++11 constexpr function a separate nc_assert(..) statement won't
    //work. Instead one can use "return nc_assert_rv(foo), bar;"
    inline constexpr int assert_always_fct_impl( bool x, const char * msg,  const char * file, unsigned line )
    {
      return  x ? 1 : throw ::NCrystal::Error::LogicError(msg,file,line);
    }
#ifndef NDEBUG
    inline constexpr int assert_fct_impl( bool x, const char * msg,  const char * file, unsigned line )
    {
      return  x ? 1 : throw ::NCrystal::Error::LogicError(msg,file,line);
    }
#  define nc_assert_rv(x) ::NCrystal::assert_details::assert_fct_impl( !!(x), "Assertion failure: " ncrystal_xstr(x), __FILE__,__LINE__)
#else
    inline constexpr bool assert_fct_impl() noexcept { return true; }
#  define nc_assert_rv(x) ::NCrystal::assert_details::assert_fct_impl()
#endif
#define nc_assert_always_rv(x) ::NCrystal::assert_details::assert_always_fct_impl( !!(x), "Assertion failure: " ncrystal_xstr(x), __FILE__,__LINE__)
  }
}

#ifdef nc_not_implemented
#  undef nc_not_implemented
#endif
#define nc_not_implemented do { NCRYSTAL_THROW(LogicError, "NotImplemented") } while(0)
#endif
