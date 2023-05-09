////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/NCException.hh"

namespace NCrystal {

  namespace Error {

    Exception::Exception(const std::string& msg, const char * f, unsigned l) throw()
      : std::runtime_error(msg),
        m_file(f),
        m_lineno(l)
    {
    }

    Exception::Exception(const char * msg,  const char * f, unsigned l) throw()
      : std::runtime_error(msg),
        m_file(f),
        m_lineno(l)
    {
    }

    Exception::Exception( const Exception & o ) throw()
      : std::runtime_error(o),
        m_file(o.m_file),
        m_lineno(o.m_lineno)
    {
    }

    Exception & Exception::operator= ( const Exception & o ) throw()
    {
      std::runtime_error::operator=(o);
      m_file = o.m_file;
      m_lineno = o.m_lineno;
      return *this;
    }

    Exception::~Exception() throw()
    {
    }

    FileNotFound::~FileNotFound() throw() {}
    MissingInfo::~MissingInfo() throw() {}
    DataLoadError::~DataLoadError() throw() {}
    LogicError::~LogicError() throw() {}
    CalcError::~CalcError() throw() {}
    BadInput::~BadInput() throw() {}
  }

}
