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

    // template <class TRequest>
    // class NCRYSTAL_API ProcessRequestBase {
    // public:
    //   //Common interface for Scatter and Absorption requests.

    //   //Info object:
    //   UniqueIDValue infoUID() const;
    //   const Info& info() const;
    //   InfoPtr infoPtr() const;
    //   bool isMultiPhase() const { return info().isMultiPhase(); }

    //   //if isMultiPhase, requests for child phases can be generated with:
    //   std::size_t nPhases() const;
    //   TRequest createChildRequest( unsigned ichild ) const;

    //   //Easily create modified request with parameters of string applied.
    //   TRequest modified( const std::string& ) const;
    //   TRequest modified( const char* ) const;
    //   //For instance, if req is a ScatterRequest we might do
    //   //  auto req_noelas = req.modified("elas=0");

    //   //Miscellaneous:
    //   void stream( std::ostream & ) const;
    //   void streamParamsOnly( std::ostream & ) const;
    //   bool operator<(const ProcessRequestBase&) const;
    //   bool operator==(const ProcessRequestBase&) const;
    //   TRequest cloneThinned() const;//only cmp operators should be used after thinning
    //   bool isThinned() const;
    //   const Cfg::CfgData& rawCfgData() const;
    //   const DataSourceName& dataSourceName() const;//for err messages only

    //   //Can only construct from "trivial" MatCfg objects (those with
    //   //isTrivial()==true), or alternatively by providing an InfoPtr and
    //   //possibly CfgData (only relevant data items will be used)
    //   ProcessRequestBase( const MatCfg& );
    //   ProcessRequestBase( InfoPtr );
    //   ProcessRequestBase( InfoPtr, const Cfg::CfgData& );
    // private:
    //   ProcessRequestBase( no_init_t ) {}
    //   struct internal_t {};
    //   ProcessRequestBase( internal_t, InfoPtr, const Cfg::CfgData* );
    //   TRequest modified( internal_t, const char*, std::size_t ) const;
    //   Cfg::CfgData m_data;
    //   OptionalInfoPtr m_infoPtr;
    //   UniqueIDValue m_infoUID;
    //   DataSourceName m_dataSourceName;
    //   bool cmpDataLT(const ProcessRequestBase&) const;
    //   bool cmpDataEQ(const ProcessRequestBase&) const;
    // };

    class NCRYSTAL_API ScatterRequest final {
      detail::ProcessRequestData m_data;
      struct internal_t{};
      ScatterRequest( internal_t, detail::ProcessRequestData&& data )
        : m_data(std::move(data)) {}
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

      static bool varIsApplicable(Cfg::detail::VarId);

      static constexpr const char * requestTypeName() noexcept { return "ScatterRequest"; }
      void checkParamConsistency() const;


      //Info object:
      UniqueIDValue infoUID() const { return m_data.infoUID(); }
      const Info& info() const { return m_data.info(); }
      InfoPtr infoPtr() const { return m_data.infoPtr(); }
      bool isMultiPhase() const { return m_data.isMultiPhase(); }

      //if isMultiPhase, requests for child phases can be generated with:
      std::size_t nPhases() const { return m_data.nPhases(); }
      ScatterRequest createChildRequest( unsigned ichild ) const
      {
        return { internal_t(),
                 m_data.createChildRequest( ichild, varIsApplicable ) };
      }

      //Easily create modified request with parameters of string applied.
      ScatterRequest modified( const std::string& str) const
      {
        return { internal_t(), m_data.modified( str, varIsApplicable ) };
      }

      ScatterRequest modified( const char* cstr ) const
      {
        return { internal_t(), m_data.modified( cstr, varIsApplicable ) };
      }

      //For instance, if req is a ScatterRequest we might do
      //  auto req_noelas = req.modified("elas=0");

      //Miscellaneous:
      void stream( std::ostream& os ) const { m_data.stream(os); }
      void streamParamsOnly( std::ostream& os ) const { m_data.streamParamsOnly(os); }
      bool operator<(const ScatterRequest& o) const { return m_data < o.m_data; }
      bool operator==(const ScatterRequest& o) const { return m_data == o.m_data; }

      ScatterRequest cloneThinned() const//only cmp operators should be used after thinning
      {
        return { internal_t(), m_data.cloneThinned() };
      }

      bool isThinned() const { return m_data.isThinned(); }
      const Cfg::CfgData& rawCfgData() const { return m_data.rawCfgData(); }
      const DataSourceName& dataSourceName() const { return m_data.dataSourceName(); }

        //Can only construct from "trivial" MatCfg objects (those with
        //isTrivial()==true), or alternatively by providing an InfoPtr and
        //possibly CfgData (only relevant data items will be used)
      ScatterRequest( const MatCfg& cfg ) : m_data( cfg, varIsApplicable ) {}
      ScatterRequest( InfoPtr info ) : m_data( std::move(info), varIsApplicable ) {}
      ScatterRequest( InfoPtr info, const Cfg::CfgData& cfgdata ) : m_data( std::move(info), cfgdata, varIsApplicable ) {}

    };

    class NCRYSTAL_API AbsorptionRequest final {
      detail::ProcessRequestData m_data;
      struct internal_t{};
      AbsorptionRequest( internal_t, detail::ProcessRequestData&& data )
        : m_data(std::move(data)) {}
    public:

      //Parameters:
      std::string get_absnfactory() const;

      static bool varIsApplicable(Cfg::detail::VarId);

      static constexpr const char * requestTypeName() noexcept { return "AbsorptionRequest"; }
      void checkParamConsistency() const;


      //Info object:
      UniqueIDValue infoUID() const { return m_data.infoUID(); }
      const Info& info() const { return m_data.info(); }
      InfoPtr infoPtr() const { return m_data.infoPtr(); }
      bool isMultiPhase() const { return m_data.isMultiPhase(); }

      //if isMultiPhase, requests for child phases can be generated with:
      std::size_t nPhases() const { return m_data.nPhases(); }
      AbsorptionRequest createChildRequest( unsigned ichild ) const
      {
        return { internal_t(),
                 m_data.createChildRequest( ichild, varIsApplicable ) };
      }

      //Easily create modified request with parameters of string applied.
      AbsorptionRequest modified( const std::string& str) const
      {
        return { internal_t(), m_data.modified( str, varIsApplicable ) };
      }

      AbsorptionRequest modified( const char* cstr ) const
      {
        return { internal_t(), m_data.modified( cstr, varIsApplicable ) };
      }

      //Miscellaneous:
      void stream( std::ostream& os ) const { m_data.stream(os); }
      void streamParamsOnly( std::ostream& os ) const { m_data.streamParamsOnly(os); }
      bool operator<(const AbsorptionRequest& o) const { return m_data < o.m_data; }
      bool operator==(const AbsorptionRequest& o) const { return m_data == o.m_data; }

      AbsorptionRequest cloneThinned() const//only cmp operators should be used after thinning
      {
        return { internal_t(), m_data.cloneThinned() };
      }

      bool isThinned() const { return m_data.isThinned(); }
      const Cfg::CfgData& rawCfgData() const { return m_data.rawCfgData(); }
      const DataSourceName& dataSourceName() const { return m_data.dataSourceName(); }

        //Can only construct from "trivial" MatCfg objects (those with
        //isTrivial()==true), or alternatively by providing an InfoPtr and
        //possibly CfgData (only relevant data items will be used)
      AbsorptionRequest( const MatCfg& cfg ) : m_data( cfg, varIsApplicable ) {}
      AbsorptionRequest( InfoPtr info ) : m_data( std::move(info), varIsApplicable ) {}
      AbsorptionRequest( InfoPtr info, const Cfg::CfgData& cfgdata ) : m_data( std::move(info), cfgdata, varIsApplicable ) {}

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
