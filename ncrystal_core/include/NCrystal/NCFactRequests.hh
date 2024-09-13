#ifndef NCrystal_FactRequests_hh
#define NCrystal_FactRequests_hh

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
#include "NCrystal/NCFactRequestsImpl.hh"

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Requests (aka "keys") for the NCrystal Info, Scatter, and Absorption        //
// factories. They can all be created directly from a (trivial) MatCfg object, //
// although in the case of Scatter and Absorption requests, creation of the    //
// request objects then implies a call to createInfo(..) behind the scenes.    //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace FactImpl {

    class NCRYSTAL_API InfoRequest final {
    public:

      //Parameters:
      Temperature get_temp() const;
      double get_dcutoff() const;
      double get_dcutoffup() const;
      std::string get_atomdb() const;
      std::vector<VectS> get_atomdb_parsed() const;
      std::string get_infofactory() const;

      //Text data:
      TextDataUID textDataUID() const;
      const std::string& getDataType() const;
      const TextData& textData() const;
      TextDataSP textDataSP() const;

      //Miscellaneous:
      void stream( std::ostream & ) const;
      void streamParamsOnly( std::ostream & ) const;
      bool operator<(const InfoRequest&) const;
      bool operator==(const InfoRequest&) const;
      const DataSourceName& dataSourceName() const;//for err messages only
      InfoRequest cloneThinned() const;//only cmp operators should be used after thinning
      bool isThinned() const;
      const Cfg::CfgData& rawCfgData() const;

      //Can only construct from "trivial" MatCfg objects (those with
      //isTrivial()==true):
      InfoRequest( const MatCfg& );

      static constexpr const char * requestTypeName() noexcept { return "InfoRequest"; }
      void checkParamConsistency() const;
    private:
      InfoRequest( no_init_t ) {}
      Cfg::CfgData m_data;
      OptionalTextDataSP m_textDataSP;
      TextDataUID m_textDataUID;
      DataSourceName m_dataSourceName;
      bool cmpDataLT(const InfoRequest&) const;
      bool cmpDataEQ(const InfoRequest&) const;
    };

    class NCRYSTAL_API ScatterRequest final : public ProcessRequestBase<ScatterRequest> {
      friend class ProcessRequestBase<ScatterRequest>;//For VSCode
    public:

      //Parameters (basic):
      int get_vdoslux() const;
      bool get_coh_elas() const;
      bool get_incoh_elas() const;
      bool get_sans() const;
      std::string get_inelas() const;
      std::string get_scatfactory() const;

      //Parameters (advanced: single crystal, ucn, ...):
      MosaicityFWHM get_mos() const;
      OrientDir get_dir1() const;
      OrientDir get_dir2() const;
      double get_mosprec() const;
      double get_sccutoff() const;
      double get_dirtol() const;
      const LCAxis& get_lcaxis() const;
      std::int_least32_t get_lcmode() const;

      bool isSingleCrystal() const;
      bool isLayeredCrystal() const;
      SCOrientation createSCOrientation() const;

      StrView get_ucnmode_str() const;
      Optional<UCNMode> get_ucnmode() const;

      //Enable constructors:
      using ProcessRequestBase<ScatterRequest>::ProcessRequestBase;
      static bool varIsApplicable(Cfg::detail::VarId);

      static constexpr const char * requestTypeName() noexcept { return "ScatterRequest"; }
      void checkParamConsistency() const;
    };

    class NCRYSTAL_API AbsorptionRequest final : public ProcessRequestBase<AbsorptionRequest> {
      friend class ProcessRequestBase<AbsorptionRequest>;//For VSCode
    public:

      //Parameters:
      std::string get_absnfactory() const;

      //Enable constructors:
      using ProcessRequestBase<AbsorptionRequest>::ProcessRequestBase;
      static bool varIsApplicable(Cfg::detail::VarId);

      static constexpr const char * requestTypeName() noexcept { return "AbsorptionRequest"; }
      void checkParamConsistency() const;
    };

    //Stream support:
    NCRYSTAL_API std::ostream& operator<<(std::ostream&, const InfoRequest&);
    NCRYSTAL_API std::ostream& operator<<(std::ostream&, const ScatterRequest&);
    NCRYSTAL_API std::ostream& operator<<(std::ostream&, const AbsorptionRequest&);
  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {

    inline TextDataUID InfoRequest::textDataUID() const { return m_textDataUID; }
    inline const std::string& InfoRequest::getDataType() const { nc_assert( m_textDataSP!=nullptr ); return m_textDataSP->dataType(); }
    inline const TextData& InfoRequest::textData() const { nc_assert( m_textDataSP!=nullptr ); return *m_textDataSP; }
    inline TextDataSP InfoRequest::textDataSP() const { nc_assert( m_textDataSP!=nullptr ); return m_textDataSP; }

    inline InfoRequest InfoRequest::cloneThinned() const
    {
      InfoRequest res{ no_init };
      res.m_data = m_data;
      res.m_textDataUID = m_textDataUID;
      res.m_dataSourceName = m_dataSourceName;
      return res;
    }

    inline bool InfoRequest::isThinned() const { return m_textDataSP == nullptr;  }
    inline const DataSourceName& InfoRequest::dataSourceName() const { return m_dataSourceName; }
    inline const Cfg::CfgData& InfoRequest::rawCfgData() const { return m_data; }

    inline bool InfoRequest::operator<( const InfoRequest& o ) const
    {
      if ( m_textDataUID != o.m_textDataUID )
        return m_textDataUID < o.m_textDataUID;
      return cmpDataLT( o );
    }
    inline bool InfoRequest::operator==( const InfoRequest& o ) const
    {
      if ( m_textDataUID != o.m_textDataUID )
        return false;
      return cmpDataEQ( o );
    }

    inline std::ostream& operator<<(std::ostream& os, const InfoRequest& req ) { req.stream(os); return os; }
    inline std::ostream& operator<<(std::ostream& os, const ScatterRequest& req ) { req.stream(os); return os; }
    inline std::ostream& operator<<(std::ostream& os, const AbsorptionRequest& req ) { req.stream(os); return os; }

  }
}

#endif
