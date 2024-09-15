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

// This file contains the implementation of the ProcessRequestData class, which
// is an implementation detail used in NCFactRequests.hh to implement the
// ScatterRequest and AbsorptionRequest classes.

#include "NCrystal/NCMatCfg.hh"
#include "NCrystal/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {
    namespace detail {

      //Validate MatCfg object ok for XXXRequest constructor:
      void validateMatCfgState( const MatCfg& );

      class NCRYSTAL_API ProcessRequestData {
      public:

        using VarIdFilter = std::function<bool(Cfg::detail::VarId)>;

        //Info object:
        UniqueIDValue infoUID() const;
        const Info& info() const;
        InfoPtr infoPtr() const;
        bool isMultiPhase() const { return info().isMultiPhase(); }

        //if isMultiPhase, requests for child phases can be generated with:
        std::size_t nPhases() const;
        ProcessRequestData createChildRequest( unsigned ichild,
                                               const VarIdFilter& ) const;

        //Easily create modified request with parameters of string applied.
        ProcessRequestData modified( const std::string&, const VarIdFilter& ) const;
        ProcessRequestData modified( const char*, const VarIdFilter& ) const;
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
        ProcessRequestData( const MatCfg&, const VarIdFilter& );
        ProcessRequestData( InfoPtr, const VarIdFilter& );
        ProcessRequestData( InfoPtr, const Cfg::CfgData&, const VarIdFilter& );

      private:

        ProcessRequestData( no_init_t ) {}
        struct internal_t {};
        ProcessRequestData( internal_t,
                            InfoPtr,
                            const Cfg::CfgData*,
                            const VarIdFilter& );
        ProcessRequestData modified( internal_t,
                                     const char*,
                                     std::size_t,
                                     const VarIdFilter& ) const;
        Cfg::CfgData m_data;
        OptionalInfoPtr m_infoPtr;
        UniqueIDValue m_infoUID;
        DataSourceName m_dataSourceName;
        bool cmpDataLT(const ProcessRequestData&) const;
        bool cmpDataEQ(const ProcessRequestData&) const;
      };

    }
  }
}

#endif
