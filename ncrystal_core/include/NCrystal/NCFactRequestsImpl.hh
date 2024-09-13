#ifndef NCrystal_FactRequestsImpl_hh
#define NCrystal_FactRequestsImpl_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/NCMatCfg.hh"
#include "NCrystal/NCInfo.hh"

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// Implementation of the ProcessRequestsBase class needed to implement      //
// AbsorptionRequest and ScatterRequest NCFactRequests.hh. This needs to be //
// defined in this separate file, to support VisualStudio compilation.      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace FactImpl {

    template <class TRequest>
    class ProcessRequestBase {
    public:
      //Common interface for Scatter and Absorption requests.

      //Info object:
      UniqueIDValue infoUID() const;
      const Info& info() const;
      InfoPtr infoPtr() const;
      bool isMultiPhase() const { return info().isMultiPhase(); }

      //if isMultiPhase, requests for child phases can be generated with:
      std::size_t nPhases() const;
      TRequest createChildRequest( unsigned ichild ) const;

      //Easily create modified request with parameters of string applied.
      TRequest modified( const std::string& ) const;
      TRequest modified( const char* ) const;
      //For instance, if req is a ScatterRequest we might do
      //  auto req_noelas = req.modified("elas=0");

      //Miscellaneous:
      void stream( std::ostream & ) const;
      void streamParamsOnly( std::ostream & ) const;
      bool operator<(const ProcessRequestBase&) const;
      bool operator==(const ProcessRequestBase&) const;
      TRequest cloneThinned() const;//only cmp operators should be used after thinning
      bool isThinned() const;
      const Cfg::CfgData& rawCfgData() const;
      const DataSourceName& dataSourceName() const;//for err messages only

      //Can only construct from "trivial" MatCfg objects (those with
      //isTrivial()==true), or alternatively by providing an InfoPtr and
      //possibly CfgData (only relevant data items will be used)
      ProcessRequestBase( const MatCfg& );
      ProcessRequestBase( InfoPtr );
      ProcessRequestBase( InfoPtr, const Cfg::CfgData& );
    private:
      ProcessRequestBase( no_init_t ) {}
      struct internal_t {};
      ProcessRequestBase( internal_t, InfoPtr, const Cfg::CfgData* );
      TRequest modified( internal_t, const char*, std::size_t ) const;
      Cfg::CfgData m_data;
      OptionalInfoPtr m_infoPtr;
      UniqueIDValue m_infoUID;
      DataSourceName m_dataSourceName;
      bool cmpDataLT(const ProcessRequestBase&) const;
      bool cmpDataEQ(const ProcessRequestBase&) const;
    };

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {

    template<typename TRequest>
    TRequest ProcessRequestBase<TRequest>::cloneThinned() const
    {
      TRequest res{ no_init };
      res.m_data = m_data;
      res.m_infoUID = m_infoUID;
      res.m_dataSourceName = m_dataSourceName;
      return res;
    }

    template <class TR>
    inline UniqueIDValue ProcessRequestBase<TR>::infoUID() const { return m_infoUID; }

    template <class TR>
    inline const Info& ProcessRequestBase<TR>::info() const { nc_assert( m_infoPtr!=nullptr ); return *m_infoPtr; }

    template <class TR>
    inline InfoPtr ProcessRequestBase<TR>::infoPtr() const {  nc_assert( m_infoPtr!=nullptr ); return m_infoPtr; }

    template <class TR>
    inline bool ProcessRequestBase<TR>::isThinned() const { return m_infoPtr == nullptr;  }

    template <class TR>
    inline const DataSourceName& ProcessRequestBase<TR>::dataSourceName() const { return m_dataSourceName; }

    template <class TR>
    inline const Cfg::CfgData& ProcessRequestBase<TR>::rawCfgData() const { return m_data; }

    template <class TR>
    inline bool ProcessRequestBase<TR>::operator<( const ProcessRequestBase& o ) const
    {
      if ( m_infoUID != o.m_infoUID )
        return m_infoUID < o.m_infoUID;
      return cmpDataLT( o );
    }

    template <class TR>
    inline bool ProcessRequestBase<TR>::operator==( const ProcessRequestBase& o ) const
    {
      if ( m_infoUID != o.m_infoUID )
        return false;
      return cmpDataEQ( o );
    }

    template <class TR>
    ProcessRequestBase<TR>::ProcessRequestBase( InfoPtr infoptr )
      : ProcessRequestBase(internal_t(), std::move(infoptr), nullptr )
    {
    }

    template <class TR>
    ProcessRequestBase<TR>::ProcessRequestBase( InfoPtr infoptr, const Cfg::CfgData& data )
      : ProcessRequestBase(internal_t(), std::move(infoptr), &data )
    {
    }

    template <class TR>
    inline std::size_t ProcessRequestBase<TR>::nPhases() const
    {
      return info().isMultiPhase() ? info().getPhases().size() : 0;
    }

  }
}

#endif
