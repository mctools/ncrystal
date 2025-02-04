#ifndef NCrystal_FactRequestsImpl_hh
#define NCrystal_FactRequestsImpl_hh

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

// This file contains the implementation of the ProcessRequestData class, which
// is an implementation detail used in NCFactRequests.hh to implement the
// ScatterRequest and AbsorptionRequest classes.

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {
    namespace detail {

      //Validate MatCfg object ok for XXXRequest constructor:
      void validateMatCfgState( const MatCfg& );

      class NCRYSTAL_API ProcessRequestData {
      public:

        struct NCRYSTAL_API ParamDefs {
          //We use this light-weight struct to capture the difference between
          //ScatterRequest and AbsorptionRequest (to avoid any sort of public
          //inheritance tree).
          typedef bool(*VarFilter)(Cfg::detail::VarId);
          typedef void(*CheckParamConsistency)(const Cfg::CfgData&);
          VarFilter varFilter = nullptr;
          CheckParamConsistency checkParamConsistency = nullptr;
        };

        //Info object:
        UniqueIDValue infoUID() const;
        const Info& info() const;
        InfoPtr infoPtr() const;
        bool isMultiPhase() const { return info().isMultiPhase(); }

        //if isMultiPhase, requests for child phases can be generated with:
        std::size_t nPhases() const;
        ProcessRequestData createChildRequest( unsigned ichild ) const;

        //Easily create modified request with parameters of string applied.
        ProcessRequestData modified( const std::string& ) const;
        ProcessRequestData modified( const char* ) const;
        //For instance, if req is a ScatterRequest we might do
        //  auto req_noelas = req.modified("elas=0");

        //Miscellaneous:
        void stream( std::ostream & ) const;
        void streamParamsOnly( std::ostream & ) const;
        bool operator<(const ProcessRequestData&) const;
        bool operator==(const ProcessRequestData&) const;
        ProcessRequestData cloneThinned() const;//only cmp operators should be used after thinning
        bool isThinned() const;
        const Cfg::CfgData& rawCfgData() const;
        const DataSourceName& dataSourceName() const;//for err messages only

        //Can only construct from "trivial" MatCfg objects (those with
        //isTrivial()==true), or alternatively by providing an InfoPtr and
        //possibly CfgData (only relevant data items will be used)
        ProcessRequestData( const MatCfg&, ParamDefs );
        ProcessRequestData( InfoPtr, ParamDefs );
        ProcessRequestData( InfoPtr, const Cfg::CfgData&, ParamDefs );

      private:
        ProcessRequestData( no_init_t ) {}
        struct internal_t {};
        ProcessRequestData( internal_t,
                            InfoPtr,
                            const Cfg::CfgData*,
                            ParamDefs );
        ProcessRequestData modified( internal_t,
                                     const char*,
                                     std::size_t ) const;
        Cfg::CfgData m_data;
        OptionalInfoPtr m_infoPtr;
        UniqueIDValue m_infoUID;
        DataSourceName m_dataSourceName;
        ParamDefs m_paramDefs;
        bool cmpDataLT(const ProcessRequestData&) const;
        bool cmpDataEQ(const ProcessRequestData&) const;
      };
    }
  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {
    namespace detail {
      inline bool ProcessRequestData::isThinned() const
      {
        return m_infoPtr == nullptr;
      }

      inline const DataSourceName& ProcessRequestData::dataSourceName() const
      {
        return m_dataSourceName;
      }

      inline const Cfg::CfgData& ProcessRequestData::rawCfgData() const
      {
        return m_data;
      }

      inline bool ProcessRequestData::operator<( const ProcessRequestData& o ) const
      {
        if ( m_infoUID != o.m_infoUID )
          return m_infoUID < o.m_infoUID;
        return cmpDataLT( o );
      }

      inline bool ProcessRequestData::operator==( const ProcessRequestData& o ) const
      {
        if ( m_infoUID != o.m_infoUID )
          return false;
        return cmpDataEQ( o );
      }

      inline ProcessRequestData::ProcessRequestData( InfoPtr infoptr,
                                              ParamDefs pd )
        : ProcessRequestData(internal_t(), std::move(infoptr), nullptr, pd )
      {
      }

      inline ProcessRequestData::ProcessRequestData( InfoPtr infoptr,
                                              const Cfg::CfgData& data,
                                              ParamDefs pd )
        : ProcessRequestData(internal_t(), std::move(infoptr), &data, pd )
      {
      }

      inline std::size_t ProcessRequestData::nPhases() const
      {
        return info().isMultiPhase() ? info().getPhases().size() : 0;
      }

      inline UniqueIDValue ProcessRequestData::infoUID() const
      {
        return m_infoUID;
      }

    }
  }
}

#endif
