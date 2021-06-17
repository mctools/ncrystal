////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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
#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/NCSCOrientation.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCAtomUtils.hh"
#include "NCrystal/internal/NCFileUtils.hh"
#include "NCrystal/internal/NCString.hh"
#include <sstream>
#include <iomanip>

namespace NC = NCrystal;

struct NC::MatCfg::Impl {
  Impl()
  {
#ifndef NDEBUG
    //Verify that parnames is sorted correctly (NB: when changing order, also
    //update partypes and PARAMETERS enum)!
    //From C++20 std::string's are constexpr-compatible, and we can do the check
    //with static_asserts instead:
    static std::atomic<unsigned> first = {0};
    if (first.fetch_add(1)==1) {
      for (int i = 1; i<PAR_NMAX; ++i) {
        nc_assert(parnames[i-1]<parnames[i]);
      }
    }
#endif
  }

  //clone:
  Impl( const Impl& o )
    : m_textDataUID( o.m_textDataUID ),
      m_textDataType( o.m_textDataType ),
      m_datafile_orig( o.m_datafile_orig ),
      m_ignoredfilecfg( o.m_ignoredfilecfg )
  {
    //clone parlist:
    for (int i = PAR_FIRST; i < PAR_NMAX; ++i) {
      if (o.m_parlist[i])
        m_parlist[i] = std::unique_ptr<ValBase>(o.m_parlist[i]->clone());
    }
  }

  void setOrientation( const SCOrientation& sco );
  std::string extractFileCfgStr( const TextData& input ) const;

  static std::pair<bool,std::string> parseIgnoreFileCfg( std::string s )
  {
    //String like "<whitespace>ignorefilecfg<whitespace>;<whatever>" should
    //return (true,"<whatever>"), otherwise return (false,s) (possibly trimmed).
    trim(s);
    constexpr const char s_ignorefilecfg[] = "ignorefilecfg";
    constexpr const unsigned n_ignorefilecfg = 13;
    static_assert(sizeof(s_ignorefilecfg)==n_ignorefilecfg+1,"");//+1 for the null char
    auto isws = [](char c) { return c==' ' || c=='\t'|| c=='\r'|| c=='\n' ; };
    if ( startswith( s, s_ignorefilecfg ) ) {
      std::string::size_type i = n_ignorefilecfg;
      while ( i < s.size() && isws( s[i] ) )
        ++i;//skip whitespace
      if ( i == s.size() )
        return { true, {} };
      if ( s[i] == ';' )
        return { true, s.substr(i+1) };
    }
    return { false, s };
  }

  void applyStrCfg( const std::string& str );

  TextDataUID m_textDataUID;
  std::string m_textDataType;
  std::string m_datafile_orig;//as passed to MatCfg constructor
  bool m_ignoredfilecfg;

  //Important!: Keep the following list in alphabetical order and synchronised
  //with parnames and partypes further down the file!
  enum PARAMETERS { PAR_FIRST = 0,
                    PAR_absnfactory = 0,
                    PAR_atomdb,
                    PAR_coh_elas,
                    PAR_dcutoff,
                    PAR_dcutoffup,
                    PAR_dir1,
                    PAR_dir2,
                    PAR_dirtol,
                    PAR_incoh_elas,
                    PAR_inelas,
                    PAR_infofactory,
                    PAR_lcaxis,
                    PAR_lcmode,
                    PAR_mos,
                    PAR_mosprec,
                    PAR_packfact,
                    PAR_scatfactory,
                    PAR_sccutoff,
                    PAR_temp,
                    PAR_vdoslux,
                    PAR_NMAX };
  using ParametersSet = std::set<PARAMETERS>;

  enum VALTYPE { VALTYPE_DBL, VALTYPE_BOOL, VALTYPE_INT, VALTYPE_STR, VALTYPE_ORIENTDIR, VALTYPE_VECTOR, VALTYPE_ATOMDB };
  static std::array<std::string,PAR_NMAX> parnames;
  static std::array<VALTYPE,PAR_NMAX> partypes;

  struct ValBase {
    virtual ~ValBase() = default;
    virtual std::unique_ptr<ValBase> clone() const = 0;
    virtual void set_from_strrep(const std::string& s) = 0;
    virtual std::string to_strrep(bool forcache) const = 0;
  };

  //Array where we keep the actual configuration:
  std::array<std::unique_ptr<ValBase>,PAR_NMAX> m_parlist;

  bool hasPar(PARAMETERS par) const { return m_parlist[par]!=nullptr; }

  struct ValDbl : public ValBase {
    enum UnitType { UnitNone, UnitAngle, UnitTemp, UnitLength };
    typedef double value_type;
    static const VALTYPE value_type_enum = VALTYPE_DBL;
    ValDbl() : ValBase(), unittype(UnitNone) {};
    virtual ~ValDbl() = default;
    std::unique_ptr<ValBase> clone() const final { return std::make_unique<ValDbl>(*this); }
    void set_from_strrep(const std::string& s) final
    {
      static std::string alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
      std::string tmporig = s;
      std::string tmp = s;
      trim(tmp);
      double unitfact = 1.0;
      double unitoffset = 0.0;
      if (unittype!=UnitNone&&tmp.size()>1&&contains(alpha,*(tmp.rbegin()))) {//*.rbegin since .back is c++11 only
        size_t iunit = tmp.size();
        while (iunit>0 && contains(alpha,tmp.at(iunit-1)))
          --iunit;
        std::string unit = tmp.substr(iunit);
        tmp.resize(iunit);
        trim(tmp);
        tmporig = tmp + unit;
        double u = -1.0;
        switch(unittype) {
        case UnitAngle:
          if (unit=="rad") { u = 1.0; }
          else if (unit=="deg") { u = kDeg; }
          else if (unit=="arcmin") { u = kArcMin; }
          else if (unit=="arcsec") { u = kArcSec; }
          break;
        case UnitLength:
          if (unit=="Aa") { u = 1.0; }
          else if (unit=="nm") { u = 10.0; }
          else if (unit=="mm") { u = 1e7;}
          else if (unit=="cm") { u = 1e8;}
          else if (unit=="m") { u = 1e10;}
          break;
        case UnitTemp:
          if (unit=="K") { u = 1.0; }
          else if (unit=="C") { u = 1.0; unitoffset = 273.15; }
          else if (unit=="F") { u = 1/1.8; unitoffset = 273.15-(32/1.8); }
          break;
        case UnitNone:
          break;
        default:
          break;
        }
        if (u<=0.0)
          NCRYSTAL_THROW2(BadInput,"Invalid unit: "<<unit);
        unitfact = u;
      }
      set(unitoffset + unitfact * str2dbl(tmp));//checks nan
      origstrrep = tmporig;
      trim(origstrrep);
    }
    void set(double v) {
      if (ncisnan(v))
        NCRYSTAL_THROW(BadInput,"Attempting to set number to NaN");
      value = v;
      origstrrep.clear();
    }
    std::string to_strrep(bool forcache) const {

      if (!forcache && !origstrrep.empty())
        return origstrrep;
      std::stringstream s;
      if (forcache)
        s << std::setprecision(16);
      s << value;
      return s.str();
    }

    //data:
    double value;
    void setUnitType(UnitType ut) { unittype = ut; }
  private:
    UnitType unittype;
    std::string origstrrep;//original input (if available), for lossless reproduction.
  };

  struct ValInt : public ValBase {
    typedef int value_type;
    static const VALTYPE value_type_enum = VALTYPE_INT;
    ValInt() : ValBase(){};
    virtual ~ValInt(){}
    std::unique_ptr<ValBase> clone() const final { return std::make_unique<ValInt>(*this); }
    void set_from_strrep(const std::string& s) final
    {
      set(str2int(s));//checks nan
    }
    void set(int v) {
      value = v;
    }
    std::string to_strrep(bool) const {
      std::stringstream s;
      s << value;
      return s.str();
    }
    //data:
    int value;
  };

  struct ValBool : public ValBase {
    typedef bool value_type;
    static const VALTYPE value_type_enum = VALTYPE_BOOL;
    ValBool() : ValBase(){}
    virtual ~ValBool(){}
    std::unique_ptr<ValBase> clone() const final { return std::make_unique<ValBool>(*this); }
    void set_from_strrep(const std::string& s) final
    {
      if (s=="true"||s=="1") { value = true; }
      else if (s=="false"||s=="0") { value = false; }
      else { NCRYSTAL_THROW2(BadInput,"Could not convert \""<<s
                             <<"\" to boolean value (should be \"true\", \"1\", \"false\" or \"0\")") }
    }
    std::string to_strrep(bool) const { return value?"true":"false"; }
    void set(bool v) { value = v; }
    bool value;
  };

  //reduce potential escaping worries in various contexts by making sure we
  //never use these special characters (in addition to a SimpleAscii check):
#define NCMATCFG_FORBIDDEN_CHARS "\"'|><(){}[]"

  struct ValStr : public ValBase {
    typedef std::string value_type;
    static const VALTYPE value_type_enum = VALTYPE_STR;
    ValStr() : ValBase() {}
    virtual ~ValStr(){}
    std::unique_ptr<ValBase> clone() const final { return std::make_unique<ValStr>(*this); }
    void set_from_strrep(const std::string& s) final { set(s); }
    void set(const std::string& s) {
      if (!isSimpleASCII(s,false,false))
        NCRYSTAL_THROW(BadInput,"Non-ASCII characters or tab/newlines in string value!");
      if (contains_any(s,NCMATCFG_FORBIDDEN_CHARS)||contains_any(s,"=;"))
        NCRYSTAL_THROW(BadInput,"Forbidden characters in string value!");
      value = s;
    }
    std::string to_strrep(bool) const { return value; }
    std::string value;
  };

  struct ValAtomDB : public ValBase {
    typedef std::vector<VectS> value_type;
    static const VALTYPE value_type_enum = VALTYPE_ATOMDB;
    ValAtomDB() : ValBase() {}
    virtual ~ValAtomDB(){}
    std::unique_ptr<ValBase> clone() const final { return std::make_unique<ValAtomDB>(*this); }
    void set_from_strrep(const std::string& s) final {
      value_type v;
      for (auto& line: split2(s,0,'@') ) {
        strreplace(line,":"," ");
        v.emplace_back( split2(line) );
      }
      set(v);
    }
    void set(const value_type& s) {
      value.clear();
      value.reserve(s.size());
      unsigned iline = 0;
      for (auto& line : s) {
        if (line.empty())
          continue;
        for (auto& word: line) {
          if (!isSimpleASCII(word,false,false))
            NCRYSTAL_THROW(BadInput,"Non-ASCII characters or tab/newlines in atomdb parameter!");
          if (contains_any(word,NCMATCFG_FORBIDDEN_CHARS)||contains_any(word,"=;"))
            NCRYSTAL_THROW(BadInput,"Forbidden characters in atomdb parameter!");
        }
        try {
          validateAtomDBLine( line );
        } catch (Error::BadInput&e) {
          NCRYSTAL_THROW2(BadInput,"Invalid entry in atomdb cfg parameter in the line: \""<<joinstr(line)<<"\". Error is: "<<e.what());
        }

        //Check for position of "nodefaults" keyword:
        if (line.size()==1&&line.at(0)=="nodefaults") {
          if (iline>0)
            NCRYSTAL_THROW2(BadInput,"Invalid entry in atomdb cfg parameter (\"nodefaults\" must be the first line)");
        }
        ++iline;
        value.push_back(std::move(line));
      }
      value_as_string = to_strrep_impl();
    }
    std::string to_strrep(bool) const { return value_as_string; }
    std::string to_strrep_impl() const {
      std::string res;
      if (value.empty())
        return res;
      auto n = value.size();
      for (decltype(n) i = 0; i < n; ++i) {
        res += joinstr( value.at(i), ":" );
        if ( i+1 < n )
          res += "@";
      }
      return res;
    }
    value_type value;
    std::string value_as_string;
  };

  struct ValOrientDir : public ValBase {
    ValOrientDir() : ValBase(){}
    virtual ~ValOrientDir(){}
    std::unique_ptr<ValBase> clone() const final { return std::make_unique<ValOrientDir>(*this); }
    void set_from_strrep(const std::string& s) final
    {
      std::string st = s; trim(st);
      VectS parts = split2(st,0,'@');
      if (parts.size()!=3||!parts.at(0).empty())
        NCRYSTAL_THROW2(BadInput,"Bad syntax for orientation: \""<<s<<"\"");
      std::string& c = parts.at(1);
      std::string& l = parts.at(2);
      int c_is_hkl(-1);
      if (startswith(c,"crys:")) { c = c.substr(5); c_is_hkl = 0; }
      else if (startswith(c,"crys_hkl:")) { c = c.substr(9); c_is_hkl = 1; }
      if (c_is_hkl==-1||!startswith(l,"lab:"))
        NCRYSTAL_THROW2(BadInput,"Bad syntax for orientation: \""<<s<<"\"");
      l = l.substr(4);
      trim(c);
      trim(l);
      VectS partsc = split2(c,0,',');
      VectS partsl = split2(l,0,',');
      if (partsc.size()!=3||partsl.size()!=3)
        NCRYSTAL_THROW2(BadInput,"Bad syntax for orientation: \""<<s<<"\"");
      set(c_is_hkl,
          str2dbl(partsc.at(0)),str2dbl(partsc.at(1)),str2dbl(partsc.at(2)),
          str2dbl(partsl.at(0)),str2dbl(partsl.at(1)),str2dbl(partsl.at(2)));
      origstrrep = s;
      trim(origstrrep);
    }
    std::string to_strrep(bool) const {
      if (!origstrrep.empty())
        return origstrrep;
      std::stringstream s;
      s.precision(17);//TODO: should we not only have high res if argument bool
                      //is true? (here and elsewhere). In fact, we should
                      //probably make sure that the argument bool can be changed
                      //in toStrCfg, and get high-res cache keys (e.g. in Geant4
                      //hooks).
      s << (crystal_is_hkl?"@crys_hkl:":"@crys:")
        << crystal[0] << "," << crystal[1] << ","
        << crystal[2] << "@lab:" << lab[0] << ","
        << lab[1] << "," << lab[2];
      return s.str();
    }
    void set(bool cishkl, double c1, double c2, double c3, double l1, double l2, double l3)
    {
      if (ncisnan(c1)||ncisnan(c2)||ncisnan(c3)||ncisnan(l1)||ncisnan(l2)||ncisnan(l3))
        NCRYSTAL_THROW(BadInput,"Attempting to set number to NaN");
      crystal_is_hkl = cishkl;
      crystal[0] = c1; crystal[1] = c2; crystal[2] = c3;
      lab[0] = l1; lab[1] = l2; lab[2] = l3;
      origstrrep.clear();
    }
    bool crystal_is_hkl;
    double crystal[3];
    double lab[3];
  private:
    std::string origstrrep;//original input (if available), for lossless reproduction.
  };

  struct ValVector : public ValBase {
    ValVector() : ValBase(){}
    virtual ~ValVector(){}
    std::unique_ptr<ValBase> clone() const final { return std::make_unique<ValVector>(*this); }
    void set_from_strrep(const std::string& s) final
    {
      std::string st = s;
      trim(st);
      VectS parts = split2(st,0,',');
      if (parts.size()!=3)
        NCRYSTAL_THROW2(BadInput,"Bad syntax for vector value: \""<<s<<"\"");
      trim(parts.at(0));
      trim(parts.at(1));
      trim(parts.at(2));
      this->set( Vector{ str2dbl(parts.at(0)),
                        str2dbl(parts.at(1)),
                        str2dbl(parts.at(2)) } );
      origstrrep = s;
      trim(origstrrep);
    }
    std::string to_strrep(bool) const {
      if (!origstrrep.empty())
        return origstrrep;
      std::stringstream s;
      s.precision(17);
      s << val[0] << "," << val[1] << "," << val[2];
      return s.str();
    }
    void set( const Vector& v )
    {
      if (ncisnan(v[0])||ncisnan(v[1])||ncisnan(v[2]))
        NCRYSTAL_THROW(BadInput,"Attempting to set number to NaN");
      val = v;
      origstrrep.clear();
    }
    Vector val;
  private:
    std::string origstrrep;//original input (if available), for lossless reproduction.
  };

  template <class ValType>
  const ValType* getValType(PARAMETERS par) const {
    const ValBase * vb = m_parlist[par].get();
    nc_assert( vb==nullptr || dynamic_cast<const ValType*>(vb) );
    return static_cast<const ValType*>(vb);
  }

  template <class ValType>
  const ValType* getValTypeThrowIfNotAvail(PARAMETERS par) const {
    const ValType * vt = getValType<ValType>(par);
    if (!vt)
      NCRYSTAL_THROW2(MissingInfo,"Value for parameter "<<parnames[par]<<" not available");
    return vt;
  }

  template <class ValType>
  void addUnitsForValType(ValType*, PARAMETERS) {}

  template <class ValType>
  ValType* getValTypeForSet(PARAMETERS par) {
    ValBase * vb = m_parlist[par].get();
    if (vb) {
      nc_assert( dynamic_cast<ValType*>(vb) );
      return static_cast<ValType*>(vb);
    } else {
      auto vt = std::make_unique<ValType>();
      addUnitsForValType<ValType>(vt.get(),par);//allow certain units for certain parameters
      auto vt_rawptr = vt.get();
      m_parlist[par] = std::move(vt);
      return vt_rawptr;
    }
  }

  template <class ValType>
  const typename ValType::value_type& getVal(PARAMETERS par, const typename ValType::value_type & code_default_val  ) const
  {
    nc_assert( ValType::value_type_enum == partypes[par] );
    const ValType * vt = getValType<ValType>(par);
    return vt ? vt->value : code_default_val;
  }

  template <class ValType>
  const typename ValType::value_type& getValNoFallback(PARAMETERS par) const
  {
    nc_assert( ValType::value_type_enum == partypes[par] );
    return getValTypeThrowIfNotAvail<ValType>(par)->value;
  }

  template <class ValType>
  void setVal(PARAMETERS par, const typename ValType::value_type& val )
  {
    getValTypeForSet<ValType>(par)->set(val);
  }

  int strNameToParIdx(const std::string& name) const {
    const std::string * itB = &parnames[0];
    const std::string * itE = itB+PAR_NMAX;
    const std::string* it = std::lower_bound(itB,itE,name);
    if ( it == itE || *it != name )
      NCRYSTAL_THROW2(BadInput,"Unknown parameter: \""<<name<<"\"");
    nc_assert( it-itB >= 0 && it-itB < PAR_NMAX );
    return it-itB;
  }

  void setValByStr( std::string name, const std::string& value )
  {
#ifndef NDEBUG
    static bool first = true;
    if (first) {
      first = false;
      nc_assert(strNameToParIdx(parnames[0])==0);
      nc_assert(strNameToParIdx(parnames[PAR_NMAX-1])==PAR_NMAX-1);
    }
#endif
    //Handle pseudo-parameters (special aliases and backwards compat.):
    if (name=="bragg") {
      name="coh_elas";
    } else if (name=="elas") {
      ValBool tmp;
      tmp.set_from_strrep(value);
      setVal<Impl::ValBool>(Impl::PAR_coh_elas,tmp.value);
      setVal<Impl::ValBool>(Impl::PAR_incoh_elas,tmp.value);
      return;
    } else if (name=="bkgd") {
      if ( value=="none" || value == "0" ) {
        setVal<Impl::ValBool>(Impl::PAR_incoh_elas,false);
        setVal<Impl::ValStr>(Impl::PAR_inelas,"none");
        return;
      } else {
        NCRYSTAL_THROW(BadInput,"The \"bkgd\" parameter is obsolete and is available for backwards compatibility "
                       "only with the values \"0\" or \"none\". For control of inelastic or incoherent-elastic "
                       "scattering, one must now instead use the parameters \"incoh_elas\" and \"inelas\".");
      }
    }

    int paridx = strNameToParIdx(name);

    if ( value.empty() && partypes[paridx] != VALTYPE_STR )//only string parameters can construct from empty strings
      NCRYSTAL_THROW2( BadInput, "Missing parameter value for parameter \""<<name<<"\"" );

    switch(partypes[paridx]) {
    case VALTYPE_DBL: getValTypeForSet<ValDbl>((PARAMETERS)paridx)->set_from_strrep(value); return;
    case VALTYPE_INT: getValTypeForSet<ValInt>((PARAMETERS)paridx)->set_from_strrep(value); return;
    case VALTYPE_BOOL: getValTypeForSet<ValBool>((PARAMETERS)paridx)->set_from_strrep(value); return;
    case VALTYPE_STR: getValTypeForSet<ValStr>((PARAMETERS)paridx)->set_from_strrep(value); return;
    case VALTYPE_ORIENTDIR: getValTypeForSet<ValOrientDir>((PARAMETERS)paridx)->set_from_strrep(value); return;
    case VALTYPE_VECTOR: getValTypeForSet<ValVector>((PARAMETERS)paridx)->set_from_strrep(value); return;
    case VALTYPE_ATOMDB: getValTypeForSet<ValAtomDB>((PARAMETERS)paridx)->set_from_strrep(value); return;
    default:
      nc_assert_always(false);
    }
  }

  static void decodeopts(const std::string& optstr, std::map<std::string,std::string>& opts2val, bool skipname = true )
  {
    opts2val.clear();
    VectS parts = split2(optstr,0,':');
    VectS::iterator it(parts.begin()), itE(parts.end());
    nc_assert_always(it!=itE);
    if (skipname)
      ++it;//skip main opt name
    VectS subparts;
    subparts.reserve(2);
    static std::string alphalowercase = "abcdefghijklmnopqrstuvwxyz";
    static std::string alphalowercasenumunderscore  = "abcdefghijklmnopqrstuvwxyz0123456789_";
    for (;it!=itE;++it) {
      trim(*it);
      if (it->empty())
        continue;
      subparts.clear();
      if (contains(*it,'@')) {
        split(subparts,*it,0,'@');
        for (std::size_t i = 0; i<subparts.size();++i)
          trim(subparts.at(i));
        if ( subparts.size()!=2||subparts.at(0).empty()||subparts.at(1).empty()||contains_any(subparts.at(1),"<>:=") ) {
          NCRYSTAL_THROW2(BadInput,"Syntax error in options: \""<<optstr<<"\"");
        }
        if ( !contains_only(subparts.at(0),alphalowercasenumunderscore)||!contains(alphalowercase,subparts.at(0)[0]) ) {
          NCRYSTAL_THROW2(BadInput,"Syntax error in options. Invalid option name: \""<<subparts.at(0)<<"\"");
        }
      } else {
        subparts.push_back(*it);
        subparts.push_back("<flag>");
      }
      if ( opts2val.find(subparts.at(0))!=opts2val.end() ) {
        NCRYSTAL_THROW2(BadInput,"Syntax error in options. Option specified multiple times: \""<<subparts.at(0)<<"\"");
      }
      opts2val[subparts.at(0)]=subparts.at(1);
    }
  }

  static std::string decodeopt_name(const std::string& optstr)
  {
    std::string name;
    if (!contains(optstr,':')) {
      name = optstr;
    } else {
      VectS parts = split2(optstr,1,':');
      nc_assert_always(!parts.empty());
      trim(parts[0]);
      name=parts[0];
    }
    return name;
  }

  static bool decodeopt_flag(const std::string& optstr, const std::string& flagname)
  {
    if (!contains(optstr,':'))
      return false;
    std::map<std::string,std::string> opts2val;
    decodeopts(optstr, opts2val );
    std::map<std::string,std::string>::const_iterator it=opts2val.find(flagname);
    if (it==opts2val.end())
      return false;
    if (it->second!="<flag>")
      NCRYSTAL_THROW2(BadInput,"Syntax error in flag: \""<<flagname<<"\" (takes no value)");
    return true;
  }

  static double decodeopt_dbl(const std::string& optstr, const std::string& parname, double defval)
  {
    if (!contains(optstr,':'))
      return defval;
    std::map<std::string,std::string> opts2val;
    decodeopts(optstr, opts2val);
    std::map<std::string,std::string>::const_iterator it=opts2val.find(parname);
    if (it==opts2val.end())
      return defval;
    return str2dbl(it->second);
  }

  static int decodeopt_int(const std::string& optstr, const std::string& parname, int defval)
  {
    if (!contains(optstr,':'))
      return defval;
    std::map<std::string,std::string> opts2val;
    Impl::decodeopts(optstr, opts2val);
    std::map<std::string,std::string>::const_iterator it=opts2val.find(parname);
    if (it==opts2val.end())
      return defval;
    return str2int(it->second);
  }

  static void decodedopt_validate(const std::string& optstr,
                                  const std::set<std::string>& recognised_opts)
  {
    if (!contains(optstr,':'))
      return;
    std::string name = decodeopt_name(optstr);
    std::map<std::string,std::string> opts2val;
    decodeopts(optstr, opts2val);
    std::map<std::string,std::string>::const_iterator it(opts2val.begin()),itE(opts2val.end());
    for (;it!=itE;++it)
      if (!recognised_opts.count(it->first)) {
        NCRYSTAL_THROW2(BadInput,"The flag \""<<it->first<<"\" is not supported by the chosen"
                        " factory for a mode of \""<<name<<"\"");
      }
  }

  void dump( const MatCfg *, std::ostream& out, bool add_endl, const MatCfg::Impl::ParametersSet * only_pars ) const;
  bool compareIgnoringTextDataUID(const MatCfg&, const MatCfg::Impl::ParametersSet * only_pars ) const;
  std::string toStrCfg( bool include_datafile, const MatCfg::Impl::ParametersSet * only_pars ) const;

  static const MatCfg::Impl::ParametersSet * onlyInfoPars() {
    static MatCfg::Impl::ParametersSet info_pars = {
      MatCfg::Impl::PAR_atomdb,
      MatCfg::Impl::PAR_dcutoff,
      MatCfg::Impl::PAR_dcutoffup,
      MatCfg::Impl::PAR_infofactory,
      MatCfg::Impl::PAR_temp
    };
    return &info_pars;
  }

  SCOrientation createSCOrientation( const MatCfg& cfg ) const
  {
    auto dir1 = getValTypeThrowIfNotAvail<ValOrientDir>(PAR_dir1);
    auto dir2 = getValTypeThrowIfNotAvail<ValOrientDir>(PAR_dir2);
    SCOrientation orient;
    double tol = cfg.get_dirtol();
    if ( dir1->crystal_is_hkl )
      orient.setPrimaryDirection( HKLPoint{dir1->crystal[0],dir1->crystal[1],dir1->crystal[2]},
                                  LabAxis{dir1->lab[0],dir1->lab[1],dir1->lab[2]} );
    else
      orient.setPrimaryDirection( CrystalAxis{dir1->crystal[0],dir1->crystal[1],dir1->crystal[2]},
                                  LabAxis{dir1->lab[0],dir1->lab[1],dir1->lab[2]} );
    if ( dir2->crystal_is_hkl )
      orient.setSecondaryDirection( HKLPoint{dir2->crystal[0],dir2->crystal[1],dir2->crystal[2]},
                                    LabAxis{dir2->lab[0],dir2->lab[1],dir2->lab[2]}, tol );
    else
      orient.setSecondaryDirection( CrystalAxis{dir2->crystal[0],dir2->crystal[1],dir2->crystal[2]},
                                    LabAxis{dir2->lab[0],dir2->lab[1],dir2->lab[2]}, tol );
    nc_assert_always(orient.isComplete());//should not be called otherwise
    return orient;
  }

};

namespace NCrystal {
  //Need fallback string values here so we can return references to them:
  static const std::string s_matcfg_str_empty = std::string();
  static const std::string s_matcfg_str_auto = std::string("auto");
  static const std::string s_matcfg_str_none = std::string("none");

  //Important!: Keep the following two lists ordered (parnames sorted
  //alphabetically) and synchronised between themselves as well as the
  //PARAMETERS enum earlier in the file!
  std::array<std::string,MatCfg::Impl::PAR_NMAX> MatCfg::Impl::parnames = { "absnfactory",
                                                   "atomdb",
                                                   "coh_elas",
                                                   "dcutoff",
                                                   "dcutoffup",
                                                   "dir1",
                                                   "dir2",
                                                   "dirtol",
                                                   "incoh_elas",
                                                   "inelas",
                                                   "infofactory",
                                                   "lcaxis",
                                                   "lcmode",
                                                   "mos",
                                                   "mosprec",
                                                   "packfact",
                                                   "scatfactory",
                                                   "sccutoff",
                                                   "temp",
                                                   "vdoslux" };
  std::array<MatCfg::Impl::VALTYPE,MatCfg::Impl::PAR_NMAX> MatCfg::Impl::partypes = { VALTYPE_STR,
                                                             VALTYPE_ATOMDB,
                                                             VALTYPE_BOOL,
                                                             VALTYPE_DBL,
                                                             VALTYPE_DBL,
                                                             VALTYPE_ORIENTDIR,
                                                             VALTYPE_ORIENTDIR,
                                                             VALTYPE_DBL,
                                                             VALTYPE_BOOL,
                                                             VALTYPE_STR,
                                                             VALTYPE_STR,
                                                             VALTYPE_VECTOR,
                                                             VALTYPE_INT,
                                                             VALTYPE_DBL,
                                                             VALTYPE_DBL,
                                                             VALTYPE_DBL,
                                                             VALTYPE_STR,
                                                             VALTYPE_DBL,
                                                             VALTYPE_DBL,
                                                             VALTYPE_INT };
  template<>
  void MatCfg::Impl::addUnitsForValType(ValDbl* vt, PARAMETERS par) {
    switch(par) {
    case PAR_mos:
    case PAR_dirtol:
      vt->setUnitType(ValDbl::UnitAngle);
      return;
    case PAR_temp:
      vt->setUnitType(ValDbl::UnitTemp);
      return;
    case PAR_dcutoff:
    case PAR_dcutoffup:
      vt->setUnitType(ValDbl::UnitLength);
      return;
    default:
      break;
    }
  }
}


bool NC::MatCfg::ignoredEmbeddedConfig() const
{
  return m_impl->m_ignoredfilecfg;
}

std::string NC::MatCfg::toEmbeddableCfg() const
{
  std::stringstream out;
  out << "NCRYSTALMATCFG[" << m_impl->toStrCfg(false,nullptr) << ']';
  return out.str();
}

std::string NC::MatCfg::Impl::toStrCfg( bool include_datafile, const MatCfg::Impl::ParametersSet * only_pars ) const
{
  std::stringstream out;
  if (include_datafile) {
    out << m_datafile_orig;
    if ( m_ignoredfilecfg )
      out << ";ignorefilecfg";
  }
  ValBase* vb;
  bool empty(out.str().empty());
  for (int i = PAR_FIRST; i<PAR_NMAX; ++i) {
    if ( ( vb = this->m_parlist[i].get() ) ) {
      if ( only_pars && only_pars->count(static_cast<decltype(PAR_FIRST)>(i))==0 )
        continue;
      if (!empty)
        out<<';';
      out << parnames[i]<<"="<<vb->to_strrep(false);
      empty = false;
    }
  }
  return out.str();
}

std::string NC::MatCfg::toStrCfg( bool include_datafile ) const
{
  return m_impl->toStrCfg( include_datafile, nullptr );
}

std::string NC::MatInfoCfg::toStrCfg( bool include_datafile ) const
{
  return m_cfg.m_impl->toStrCfg( include_datafile, MatCfg::Impl::onlyInfoPars() );
}

bool NC::MatCfg::isSingleCrystal() const
{
  return m_impl->hasPar(Impl::PAR_mos) || m_impl->hasPar(Impl::PAR_dir1) ||
    m_impl->hasPar(Impl::PAR_dir2) || m_impl->hasPar(Impl::PAR_dirtol);
}

bool NC::MatCfg::isLayeredCrystal() const
{
  return m_impl->hasPar(Impl::PAR_lcaxis);
}

void NC::MatCfg::checkConsistency() const
{
  const double parval_temp = get_temp().get();
  const double parval_dcutoff = get_dcutoff();
  const double parval_dcutoffup = get_dcutoffup();
  const double parval_packfact = get_packfact();
  const double parval_dirtol = get_dirtol();
  const double parval_sccutoff = get_sccutoff();
  if ( parval_temp!=-1.0 && (parval_temp<0.0||parval_temp>1e5) )
    NCRYSTAL_THROW(BadInput,"temp must be -1.0 or in the range (0.0,1e5]");
  if (parval_dcutoff!=-1) {
    if (parval_dcutoff<0.0)
      NCRYSTAL_THROW(BadInput,"dcutoff must be -1.0 or >=0.0");
    if (parval_dcutoff>=parval_dcutoffup)
      NCRYSTAL_THROW(BadInput,"dcutoff must be less than dcutoffup");
    if (!(parval_dcutoff>=1e-3&&parval_dcutoff<=1e5) && parval_dcutoff!=0 )
      NCRYSTAL_THROW(BadInput,"dcutoff must be -1 (hkl lists disabled), 0 (for automatic selection), or in range [1e-3,1e5]");
  }
  if (parval_packfact<=0.0||parval_packfact>1.0)
    NCRYSTAL_THROW(BadInput,"packfact must be in range (0.0,1.0]");
  if (parval_sccutoff<0.0)
    NCRYSTAL_THROW(BadInput,"sccutoff must be >=0.0");
  if (parval_dirtol<=0.0||parval_dirtol>kPi)
    NCRYSTAL_THROW(BadInput,"dirtol must be in range (0.0,pi]");
  const double parval_mosprec = get_mosprec();
  if ( ! (valueInInterval(0.9999e-7,0.10000001,parval_mosprec) ) )
    NCRYSTAL_THROW(BadInput,"mosprec must be in the range [1e-7,1e-1].");

  //inelas:
  std::string parval_inelas = get_inelas();
  if (parval_inelas.empty()||!contains_only(parval_inelas,"abcdefghijklmnopqrstuvwxyz_0123456789"))
    NCRYSTAL_THROW2(BadInput,"invalid inelas name specified: \""<<parval_inelas<<"\"");

  //infofactory:
  std::string parval_infofactory = get_infofactory();
  std::string parval_infofact_name = get_infofact_name();
  if (!contains_only(parval_infofact_name,"abcdefghijklmnopqrstuvwxyz_0123456789"))
    NCRYSTAL_THROW2(BadInput,"invalid infofactory name specified: \""<<parval_infofact_name<<"\"");
  if (parval_infofact_name.empty()&&contains(parval_infofactory,':'))
    NCRYSTAL_THROW2(BadInput,"infofactory options not allowed when not specifying specific factory");
  std::map<std::string,std::string> opts2val;
  Impl::decodeopts(parval_infofactory, opts2val);//decode to trigger any BadInput errors here

  //Now check the 4 SC parameters, only 1 of which has a code fallback value:
  int nOrient = (m_impl->hasPar(Impl::PAR_dir1)?1:0)
    + (m_impl->hasPar(Impl::PAR_dir2)?1:0)
    + (m_impl->hasPar(Impl::PAR_mos)?1:0);
  if (nOrient!=0 && nOrient<3)
    NCRYSTAL_THROW(BadInput,"Must set all or none of mos, dir1 and dir2 parameters");
  if (nOrient==0&&m_impl->hasPar(Impl::PAR_dirtol))
    NCRYSTAL_THROW(BadInput,"mos, dir1 and dir2 parameters must all be set when dirtol is set");

  if (nOrient) {
    //Check the validity of last SC parameters here!
    const double parval_mos = get_mos().dbl();

    if ( !( parval_mos > 0.0) || parval_mos > kPiHalf )
      NCRYSTAL_THROW(BadInput,"mos must be in range (0.0,pi/2]");
    //should be single crystal
    if (parval_packfact!=1.0)
      NCRYSTAL_THROW(BadInput,"Single crystal parameters are set, so packfact must be 1.0");

    //validate orientations by constructing SCOrientation object:
    (void)m_impl->createSCOrientation(*this);

  } else {
    //should be polycrystal. No extra validation needed for now, packfact was already validated above.
  }

  if (m_impl->hasPar(Impl::PAR_lcaxis)) {
    nc_assert(isLayeredCrystal());
    auto v = get_lcaxis();
    nc_assert_always(! ( ncisnan(v[0]) || ncisnan(v[1]) || ncisnan(v[2]) ) );//should have been caught
    double mag = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if ( ncisinf(mag) || ncisinf(v[0]) || ncisinf(v[1]) || ncisinf(v[2]) )
      NCRYSTAL_THROW(BadInput, "Infinities or too large values specified in lcaxis vector");
    if (!mag)
      NCRYSTAL_THROW(BadInput, "Null vector or too small values specified in lcaxis vector");
  }

  const int parval_vdoslux = get_vdoslux();
  if ( parval_vdoslux < 0 || parval_vdoslux > 5 ) {
    NCRYSTAL_THROW2(BadInput, "Specified invalid vdoslux value of "
                    <<parval_vdoslux<<" (must be integer from 0 to 5)");
  }

}

void NC::MatCfg::set_lcaxis( const LCAxis& axis )
{
  m_impl.modify()->getValTypeForSet<Impl::ValVector>(Impl::PAR_lcaxis)->set(axis.as<Vector>());
}

NC::LCAxis NC::MatCfg::get_lcaxis() const
{
  return m_impl->getValTypeThrowIfNotAvail<Impl::ValVector>(Impl::PAR_lcaxis)->val.as<LCAxis>();
}

void NC::MatCfg::set_dir1( const HKLPoint& c, const LabAxis& l )
{
  m_impl.modify()->getValTypeForSet<Impl::ValOrientDir>(Impl::PAR_dir1)->set(true,c[0],c[1],c[2],l[0],l[1],l[2]);
}

void NC::MatCfg::set_dir1( const CrystalAxis& c, const LabAxis& l )
{
  m_impl.modify()->getValTypeForSet<Impl::ValOrientDir>(Impl::PAR_dir1)->set(false,c[0],c[1],c[2],l[0],l[1],l[2]);
}

void NC::MatCfg::set_dir2( const HKLPoint& c, const LabAxis& l )
{
  m_impl.modify()->getValTypeForSet<Impl::ValOrientDir>(Impl::PAR_dir2)->set(true,c[0],c[1],c[2],l[0],l[1],l[2]);
}

void NC::MatCfg::set_dir2( const CrystalAxis& c, const LabAxis& l )
{
  m_impl.modify()->getValTypeForSet<Impl::ValOrientDir>(Impl::PAR_dir2)->set(false,c[0],c[1],c[2],l[0],l[1],l[2]);
}

std::pair<NC::Variant<NC::CrystalAxis,NC::HKLPoint>,NC::LabAxis> NC::MatCfg::get_dir1() const
{
  SCOrientation sco = createSCOrientation();
  nc_assert_always(sco.getLabDir(0).has_value());
  return { sco.getCrysDir(0), sco.getLabDir(0).value() };
}

std::pair<NC::Variant<NC::CrystalAxis,NC::HKLPoint>,NC::LabAxis> NC::MatCfg::get_dir2() const
{
  SCOrientation sco = createSCOrientation();
  nc_assert_always(sco.getLabDir(1).has_value());
  return { sco.getCrysDir(1), sco.getLabDir(1).value() };
}

void NC::MatCfg::setOrientation( const SCOrientation& sco )
{
  if (!sco.isComplete())
    NCRYSTAL_THROW(BadInput,"setOrientation called with incomplete SCOrientation object");
  auto mod = m_impl.modify();
  mod->setOrientation(sco);
  nc_assert(isSingleCrystal());
}

void NC::MatCfg::Impl::setOrientation( const SCOrientation& sco )
{
  if (!sco.isComplete())
    NCRYSTAL_THROW(BadInput,"Incomplete SCOrientation object - must set both primary and secondary directions.");
  ValOrientDir* p[2];
  p[0] = getValTypeForSet<ValOrientDir>(PAR_dir1);
  p[1] = getValTypeForSet<ValOrientDir>(PAR_dir2);
  nc_assert(p[0]&&p[1]);
  for ( int i = 0; i < 2; ++i ) {
    auto crysdir = sco.getCrysDir(i);
    nc_assert_always(!crysdir.empty());
    auto opt_labdir = sco.getLabDir(i);
    nc_assert_always(opt_labdir.has_value());
    auto labdir = opt_labdir.value();
    bool is_hkl(crysdir.has_value<HKLPoint>());
    Vector& cv( is_hkl? crysdir.get<HKLPoint>().as<Vector>() : crysdir.get<CrystalAxis>().as<Vector>() );
    p[i]->set( is_hkl, cv[0],cv[1],cv[2], labdir[0],labdir[1],labdir[2] );
  }
  setVal<ValDbl>(PAR_dirtol,sco.getTolerance());
}

NC::SCOrientation NC::MatCfg::createSCOrientation() const
{
  checkConsistency();//TODO: creates scorientation twice!
  if (!isSingleCrystal())
    NCRYSTAL_THROW(MissingInfo,"Can only create SCOrientation object for single crystals");
  if ( ! m_impl->hasPar(Impl::PAR_dir1) )
    NCRYSTAL_THROW(MissingInfo,"Can not create SCOrientation object without the dir1 parameter set");
  if ( ! m_impl->hasPar(Impl::PAR_dir2) )
    NCRYSTAL_THROW(MissingInfo,"Can not create SCOrientation object without the dir2 parameter set");
  return m_impl->createSCOrientation(*this);
}

void NC::MatCfg::applyStrCfg( const std::string& str )
{
  m_impl.modify()->applyStrCfg(str);
}

void NC::MatCfg::Impl::applyStrCfg( const std::string& str )
{
  if (!isSimpleASCII(str,true,true))
    NCRYSTAL_THROW(BadInput,"Non-ASCII characters in parameter specification!");

  if (contains_any(str,NCMATCFG_FORBIDDEN_CHARS))
    NCRYSTAL_THROW(BadInput,"Forbidden characters in parameter specification!");

  VectS par_and_val;
  VectS parts = split2(str,0,';');

  for (size_t i = 0; i<parts.size();++i) {
    trim(parts.at(i));
    if (parts.at(i).empty()) {
      //be flexible and simply ignore missing parts (so for instance
      //MatCfg("myfile.ncmat;") will still work).
      continue;
    }
    if (parts.at(i)=="ignorefilecfg") {
      NCRYSTAL_THROW2(BadInput,"The \"ignorefilecfg\" keyword can only be used in the MatCfg "
                      "constructor (and only directly after the filename)");
    }
    split(par_and_val,parts.at(i),0,'=');
    if (par_and_val.size()!=2) {
      NCRYSTAL_THROW2(BadInput,"Bad syntax in parameter specification: \""<<parts.at(i)<<"\"");
    }
    trim(par_and_val.at(0));
    trim(par_and_val.at(1));
    if (par_and_val.at(0).empty())
      NCRYSTAL_THROW(BadInput,"Missing parameter name");
    this->setValByStr(par_and_val.at(0),par_and_val.at(1));
  }
}

NC::MatCfg::MatCfg( const char* datafile_and_parameters )
  : MatCfg(std::string(datafile_and_parameters))
{
}

NC::MatCfg NC::MatCfg::createFromRawData( std::string&& data, std::string pars, std::string ext )
{
  return MatCfg( from_raw_t(), std::move(data), std::move(pars), std::move(ext) );
}

NC::MatCfg::MatCfg( const std::string& datafile_and_parameters )
  : MatCfg( [&datafile_and_parameters]() -> constructor_args
            {
              //Trim and split (once) on ';', throwing away empty parts:
              std::string input(datafile_and_parameters);
              trim( input );
              VectS parts = split2(input,1,';');
              for ( auto& e : parts )
                trim(e);
              //First and only required parameter is the datafile:
              if ( parts.empty() || parts.at(0).empty() )
                NCRYSTAL_THROW(MissingInfo,"Please supply name of data file");
              if (contains(parts.at(0),'=')) {
                //catch typical user error
                NCRYSTAL_THROW2(BadInput,"Filename contains a forbidden character ('='): "<<parts.at(0));
              }
              constructor_args args;
              if ( parts.size()>1 )
                args.pars = std::move(parts.at(1));
              args.origfn = std::move(parts.at(0));
              nc_assert( ! args.origfn.empty() );
              args.td = FactImpl::createTextData( args.origfn );
              return args;
            }())
{
}

NC::MatCfg::MatCfg( TextDataSP sp, std::string pars )
  : MatCfg( [&sp,&pars]() -> constructor_args
  {
    constructor_args args;
    args.td = std::move(sp);
    args.pars = std::move(pars);
    return args;
  }())
{
}

NC::MatCfg::MatCfg( from_raw_t, std::string&& data, std::string pars, std::string dataType )
  : MatCfg( [&data,&pars,&dataType]() -> constructor_args
  {
    constructor_args args;
    RawStrData rawdata(std::move(data));
    if ( dataType.empty() )
      dataType = FactImpl::guessDataType(rawdata);
    if ( dataType.empty() )
      NCRYSTAL_THROW2(BadInput,"Can not determine format of anonymous data (must be specified explicitly in this case):");
    args.td = makeSO<const TextData>( std::move(rawdata), TextData::DataType{std::move(dataType)} );
    args.pars = std::move(pars);
    return args;
  }())
{
}

NC::MatCfg::MatCfg( constructor_args&& args )
{
  auto mod = m_impl.modifyWithoutLocking();//we just constructed m_impl from
                                           //scratch, no need to lock as no one
                                           //else can refer to it.

  const TextData& textData = *args.td;
  m_textDataSP = std::move(args.td);

  //Cache TextData meta data which should remain even if thinned:
  mod->m_textDataUID = textData.dataUID();
  mod->m_textDataType = textData.dataType();

  //Original file name is needed for serialisation:
  mod->m_datafile_orig = args.origfn;

  //Apply parameter config "parname1=val1;...", combining the parts embedded in
  //the input data (if any) with the ones in args.pars, taking into account (and
  //noting down) the presence of an "ignorefilecfg" keyword.

  std::tie(mod->m_ignoredfilecfg,args.pars) = Impl::parseIgnoreFileCfg(args.pars);
  std::string cfgstr = ( mod->m_ignoredfilecfg
                         ? std::string()
                         : m_impl->extractFileCfgStr(textData) );
  if ( cfgstr.empty() )
    cfgstr = args.pars;
  else if ( !args.pars.empty() ) {
    cfgstr += ';';
    cfgstr += args.pars;
  }

  if (!cfgstr.empty())
    mod->applyStrCfg( cfgstr );
}

std::string NC::MatCfg::Impl::extractFileCfgStr( const TextData& input ) const
{
  std::string res;
  std::string pattern="NCRYSTALMATCFG";
  for ( const std::string& line : input ) {
    std::size_t pos = line.find(pattern);
    if ( pos == std::string::npos )
      continue;
    if (!res.empty())
      NCRYSTAL_THROW2(BadInput,"Input file contains more than one "<<pattern<<" specification: "<<m_datafile_orig);
    std::string s = line.substr(pos+pattern.size());
    if (s.empty()||s.at(0)!='[')
      NCRYSTAL_THROW2(BadInput,"Input file contains "<<pattern<<" which is not followed by a '[' character: "<<m_datafile_orig);
    if (s.find(pattern)!=std::string::npos)
      NCRYSTAL_THROW2(BadInput,"Input file contains more than one "<<pattern<<" specification on a single line: "<<m_datafile_orig);
    s = s.substr(1);
    pos = s.find(']');
    if ( pos == std::string::npos )
      NCRYSTAL_THROW2(BadInput,"Input file contains "<<pattern<<" without a closing ']' character: "<<m_datafile_orig);
    res = s.substr(0,pos);
    if (res.empty())
      res = " ";//for detection of multiple occurances
  }
  trim(res);
  return res;
}

//All this nice stuff is generated by the COWPimpl helper:
NC::MatCfg::MatCfg(const MatCfg&) = default;
NC::MatCfg& NC::MatCfg::operator=(const MatCfg&) = default;
NC::MatCfg::MatCfg( MatCfg&& ) = default;
NC::MatCfg& NC::MatCfg::operator=(MatCfg&&) = default;
NC::MatCfg::~MatCfg() = default;

NC::Temperature NC::MatCfg::get_temp() const { return Temperature{ m_impl->getVal<Impl::ValDbl>(Impl::PAR_temp,-1.0) }; }
double NC::MatCfg::get_dcutoff() const { return m_impl->getVal<Impl::ValDbl>(Impl::PAR_dcutoff,0.0); }
double NC::MatCfg::get_dcutoffup() const { return m_impl->getVal<Impl::ValDbl>(Impl::PAR_dcutoffup,kInfinity); }
double NC::MatCfg::get_packfact() const { return m_impl->getVal<Impl::ValDbl>(Impl::PAR_packfact,1.0); }
NC::MosaicityFWHM NC::MatCfg::get_mos() const { return MosaicityFWHM{ m_impl->getValNoFallback<Impl::ValDbl>(Impl::PAR_mos) }; }
double NC::MatCfg::get_mosprec() const { return m_impl->getVal<Impl::ValDbl>(Impl::PAR_mosprec,1e-3); }
double NC::MatCfg::get_sccutoff() const { return m_impl->getVal<Impl::ValDbl>(Impl::PAR_sccutoff,0.4); }
double NC::MatCfg::get_dirtol() const { return m_impl->getVal<Impl::ValDbl>(Impl::PAR_dirtol,1e-4); }
bool NC::MatCfg::get_coh_elas() const { return m_impl->getVal<Impl::ValBool>(Impl::PAR_coh_elas,true); }
bool NC::MatCfg::get_incoh_elas() const { return m_impl->getVal<Impl::ValBool>(Impl::PAR_incoh_elas,true); }
const std::string& NC::MatCfg::get_inelas() const {
  const std::string& ss = m_impl->getVal<Impl::ValStr>(Impl::PAR_inelas,s_matcfg_str_auto);
  if (isOneOf(ss,"none","0","sterile","false"))
    return s_matcfg_str_none;
  return ss;
 }
const std::string& NC::MatCfg::get_infofactory() const { return m_impl->getVal<Impl::ValStr>(Impl::PAR_infofactory,s_matcfg_str_empty); }
const std::string& NC::MatCfg::get_scatfactory() const { return m_impl->getVal<Impl::ValStr>(Impl::PAR_scatfactory,s_matcfg_str_empty); }
const std::string& NC::MatCfg::get_absnfactory() const { return m_impl->getVal<Impl::ValStr>(Impl::PAR_absnfactory,s_matcfg_str_empty); }
void NC::MatCfg::set_temp( Temperature v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValDbl>(Impl::PAR_temp,v.get()); }
void NC::MatCfg::set_dcutoff( double v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValDbl>(Impl::PAR_dcutoff,v); }
void NC::MatCfg::set_dcutoffup( double v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValDbl>(Impl::PAR_dcutoffup,v); }
void NC::MatCfg::set_packfact( double v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValDbl>(Impl::PAR_packfact,v); }
void NC::MatCfg::set_mos( MosaicityFWHM v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValDbl>(Impl::PAR_mos,v.dbl()); }
void NC::MatCfg::set_mosprec( double v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValDbl>(Impl::PAR_mosprec,v); }
void NC::MatCfg::set_sccutoff( double v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValDbl>(Impl::PAR_sccutoff,v); }
void NC::MatCfg::set_dirtol( double v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValDbl>(Impl::PAR_dirtol,v); }
void NC::MatCfg::set_coh_elas( bool v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValBool>(Impl::PAR_coh_elas,v); }
void NC::MatCfg::set_incoh_elas( bool v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValBool>(Impl::PAR_incoh_elas,v); }
void NC::MatCfg::set_inelas( const std::string& v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValStr>(Impl::PAR_inelas,v); }
void NC::MatCfg::set_infofactory( const std::string& v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValStr>(Impl::PAR_infofactory,v); }
void NC::MatCfg::set_scatfactory( const std::string& v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValStr>(Impl::PAR_scatfactory,v); }
void NC::MatCfg::set_absnfactory( const std::string& v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValStr>(Impl::PAR_absnfactory,v); }
void NC::MatCfg::set_lcmode( int v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValInt>(Impl::PAR_lcmode,v); }
int NC::MatCfg::get_lcmode() const { return m_impl->getVal<Impl::ValInt>(Impl::PAR_lcmode,0); }
void NC::MatCfg::set_vdoslux( int v ) { auto mod = m_impl.modify(); mod->setVal<Impl::ValInt>(Impl::PAR_vdoslux,v); }
int NC::MatCfg::get_vdoslux() const { return m_impl->getVal<Impl::ValInt>(Impl::PAR_vdoslux,3); }

const std::string& NC::MatCfg::get_atomdb() const {
  const Impl::ValAtomDB * vt = m_impl->getValType<Impl::ValAtomDB>(Impl::PAR_atomdb);
  return vt ? vt->value_as_string : s_matcfg_str_empty;
}

const std::vector<NC::VectS>& NC::MatCfg::get_atomdb_parsed() const
{
  const Impl::ValAtomDB * vt = m_impl->getValType<Impl::ValAtomDB>(Impl::PAR_atomdb);
  static decltype(vt->value) s_empty_atomdb;
  return vt ? vt->value : s_empty_atomdb;
}

void NC::MatCfg::set_atomdb( const std::string& v ) {
  m_impl.modify()->getValTypeForSet<Impl::ValAtomDB>(Impl::PAR_atomdb)->set_from_strrep(v);
}

bool NC::MatCfg::isPolyCrystal() const { return !isSingleCrystal(); }

std::string NC::MatCfg::get_infofact_name() const
{
  return Impl::decodeopt_name( get_infofactory() );
}

bool NC::MatCfg::get_infofactopt_flag(const std::string& flagname) const
{
  return Impl::decodeopt_flag(get_infofactory(),flagname);
}

double NC::MatCfg::get_infofactopt_dbl(const std::string& flagname, double defval) const
{
  return Impl::decodeopt_dbl(get_infofactory(),flagname,defval);
}

int NC::MatCfg::get_infofactopt_int(const std::string& flagname, int defval) const
{
  return Impl::decodeopt_int(get_infofactory(),flagname,defval);
}

void NC::MatCfg::infofactopt_validate(const std::set<std::string>& recognised_opts) const
{
  Impl::decodedopt_validate(get_infofactory(),recognised_opts);
}


namespace NCrystal {
  namespace {
    std::string factRequestsToString(const NC::MatCfg::FactRequested& req) {
      VectS parts;
      if (!req.specific.empty())
        parts.push_back(req.specific);
      for ( auto& e: req.excluded ) {
        parts.emplace_back("!");
        parts.back() += e;
      }
      return joinstr(parts,"@");
    }
    NC::MatCfg::FactRequested parseFactRequests( const std::string& requests, const char* parname ) {
      NC::MatCfg::FactRequested res;
      for (auto& e : split2(requests,0,'@') ) {
        trim(e);
        if (e.empty())
          continue;
        if (startswith(e,"!")) {
          auto exlname = e.substr(1);
          trim(exlname);
          if (!exlname.empty())
            res.excluded.insert(exlname);
          continue;
        }
        if (!res.specific.empty()) {
          NCRYSTAL_THROW2(BadInput,parname<<" parameter contains more than one entry (both \""
                          <<res.specific<<"\" and \""<<e<<"\") which is not supported (only negated "
                          "entries starting with \"!\" can appear in any number).");
        }
        res.specific = e;
      }
      if (!res.specific.empty()&&res.excluded.count(res.specific))
        NCRYSTAL_THROW2(BadInput,parname<<" parameter contains the name \""
                        <<res.specific<<"\" both as a required and as an excluded entry.");

      return res;
    }
  }
}

NC::MatCfg::FactRequested NC::MatCfg::get_scatfactory_parsed() const
{
  return parseFactRequests( get_scatfactory(), "scatfactory" );
}

NC::MatCfg::FactRequested NC::MatCfg::get_absnfactory_parsed() const
{
  return parseFactRequests( get_absnfactory(), "absnfactory");
}

void NC::MatCfg::set_scatfactory(const NC::MatCfg::FactRequested& req )
{
  set_scatfactory(factRequestsToString(req));
}

void NC::MatCfg::set_absnfactory(const NC::MatCfg::FactRequested& req )
{
  set_absnfactory(factRequestsToString(req));
}

NC::MatInfoCfg NC::MatCfg::createInfoCfg() const
{
  return { *this };
}

void NC::MatCfg::Impl::dump( const MatCfg * self,
                             std::ostream& out,
                             bool add_endl,
                             const MatCfg::Impl::ParametersSet * only_pars ) const
{
  std::string strcfg = self->m_impl->toStrCfg( false, only_pars );
  out << "MatCfg(\"";
  if ( !m_datafile_orig.empty() ) {
    out << m_datafile_orig;
  } else {
    auto dt = self->getDataType();
    if ( dt.empty() )
      out << "<anonymous-data>";
    else
      out << "<anonymous-"<<dt<<"-data>";
  }
  if (this->m_ignoredfilecfg)
    out << ";ignorefilecfg";
  if (!strcfg.empty())
    out << (strcfg[0]==';'?"":";") << strcfg;
  out<<"\")";
  if (add_endl)
    out<<std::endl;
}

void NC::MatInfoCfg::dump( std::ostream &out, bool add_endl ) const
{
  m_cfg.m_impl->dump(&m_cfg,out,add_endl,MatCfg::Impl::onlyInfoPars());
}

void NC::MatCfg::dump( std::ostream& out, bool add_endl ) const
{
  m_impl->dump(this,out,add_endl,nullptr);
}

const NC::TextDataUID NC::MatCfg::textDataUID() const
{
  return m_impl->m_textDataUID;
}

const std::string& NC::MatCfg::getDataType() const
{
  return m_impl->m_textDataType;
}

NC::MatCfg NC::MatCfg::cloneThinned() const
{
  MatCfg cfg(*this);
  cfg.m_textDataSP.reset();
  return cfg;
}

NC::TextDataSP NC::MatCfg::textDataSP() const
{
  if ( m_textDataSP == nullptr )
    NCRYSTAL_THROW(LogicError,"MatCfg::textDataSP/textData methods should not be"
                   " used in a MatCfg object which was thinned or moved-from.");
  return m_textDataSP;
}

bool NC::MatCfg::Impl::compareIgnoringTextDataUID(const MatCfg& o, const MatCfg::Impl::ParametersSet * only_pars ) const
{
  const Impl * oimpl = &*o.m_impl;
  if ( this == oimpl )
    return false;//same internal data instance, must be equal

  if ( this->m_datafile_orig != oimpl->m_datafile_orig )
    return this->m_datafile_orig < oimpl->m_datafile_orig;
  if ( this->m_textDataType != oimpl->m_textDataType )
    return this->m_textDataType < oimpl->m_textDataType;

  //To keep things simple we also compare "ignoredfilecfg", although for the
  //purposes of e.g. creating Info objects from MatCfg objects, it should not
  //play a role:
  if ( this->m_ignoredfilecfg != oimpl->m_ignoredfilecfg )
    return this->m_ignoredfilecfg;

  //Ok, just the actual parameters left:

  for (int i = Impl::PAR_FIRST; i<Impl::PAR_NMAX; ++i) {
    Impl::ValBase * v = this->m_parlist[i].get();
    Impl::ValBase * vo = oimpl->m_parlist[i].get();
    if ( !v && !vo )
      continue;//neither is set
    if ( only_pars && only_pars->count(static_cast<decltype(Impl::PAR_FIRST)>(i))==0 )
      continue;
    if ( bool(v) != bool(vo) )
      return int(bool(v)) < int(bool(vo));//one has set
    nc_assert( v && vo );
    //both set (rare, we only get here for the same input file after all).
    std::string vs = v->to_strrep(true);
    std::string vos = vo->to_strrep(true);
    if ( vs != vos )
      return vs < vos;
  }
  //identical:
  return false;
}

bool NC::MatCfg::operator<( const MatCfg& o ) const
{
  if ( m_impl->m_textDataUID != o.m_impl->m_textDataUID )
    return m_impl->m_textDataUID < o.m_impl->m_textDataUID;
  return m_impl->compareIgnoringTextDataUID(o,nullptr);
}

bool NC::MatInfoCfg::operator<( const MatInfoCfg& o ) const
{
  if ( m_cfg.m_impl->m_textDataUID != o.m_cfg.m_impl->m_textDataUID )
    return m_cfg.m_impl->m_textDataUID < o.m_cfg.m_impl->m_textDataUID;
  return m_cfg.m_impl->compareIgnoringTextDataUID( o.m_cfg, MatCfg::Impl::onlyInfoPars() );
}
