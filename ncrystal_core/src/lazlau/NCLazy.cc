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

#include "NCLazy.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/atomdb/NCAtomDBExtender.hh"
#include "NCrystal/internal/infobld/NCInfoBuilder.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace Lazy {
    namespace {
      struct CollectedData {
        enum class ColType { Dbl, PosInt, Str };
        static std::map<std::string,ColType>& colDefs() {
          static std::map<std::string,ColType> ss =
          {
            {"column_h",ColType::PosInt},
            {"column_k",ColType::PosInt},
            {"column_l",ColType::PosInt},
            {"column_F",ColType::PosInt},
            {"column_F2",ColType::PosInt},
            {"column_d",ColType::PosInt},
            {"column_j",ColType::PosInt},
            {"lattice_a",ColType::Dbl},
            {"lattice_b",ColType::Dbl},
            {"lattice_c",ColType::Dbl},
            {"lattice_aa",ColType::Dbl},
            {"lattice_bb",ColType::Dbl},
            {"lattice_cc",ColType::Dbl},
            {"SPCGRP",ColType::Str},//NB: Has special parsing, search for SPCGRP below!
            {"TITLE",ColType::Str},//NB: Has special parsing, search for TITLE below!
            //Our extensions:
            {"spacegroup",ColType::PosInt},
            {"temperature",ColType::Dbl},
            {"debye_temperature",ColType::Dbl},
            {"formula",ColType::Str},
            {"nformula_per_unitcell",ColType::PosInt},
            {"unit_F2",ColType::Str},
          };
          return ss;
        }

        class ParsedHdr {
        public:
          bool hasKey( const std::string& key )
          {
            return ( m_posints.find(key)!=m_posints.end()
                     || m_dbls.find(key)!=m_dbls.end()
                     || m_strs.find(key)!=m_strs.end() );
          }
          double getDbl(const std::string& key ) const { nc_assert(m_dbls.count(key)); return m_dbls.find(key)->second; }
          unsigned getPosInt(const std::string& key ) const { nc_assert(m_posints.count(key)); return m_posints.find(key)->second; }
          std::string getStr(const std::string& key ) const { nc_assert(m_strs.count(key)); return m_strs.find(key)->second; }

          double value_or( const std::string& key, double fallback )
          {
            auto it = m_dbls.find(key);
            return it == m_dbls.end() ? fallback : it->second;
          }

          unsigned value_or( const std::string& key, unsigned fallback )
          {
            auto it = m_posints.find(key);
            return it == m_posints.end() ? fallback : it->second;
          }

          template<class TMap, class TVal>
          void insertValue( const StrView& errprefix, const std::string& key, TMap& themap, const TVal& val ) {
            if ( themap.find(key)!=themap.end() && themap.find(key)->second != val )
              NCRYSTAL_THROW2(BadInput,errprefix<<"Key \""<<key<<"\" specified more than once (and with different values).");
            themap[key] = val;
          }
          void updateIfNonZero(const std::string& key, double value) {
            nc_assert(colDefs().count(key)&&colDefs().find(key)->second == ColType::Dbl);
            if ( ! ( value == 0.0 ) )
              m_dbls[key]=value;
          }

          void addLine( StrView errprefix, StrView l ) {
            nc_assert(l.startswith('#'));
            auto parts = l.substr(1).ltrimmed().split();
            if ( parts.size() < 1 )
              return;
            const auto& kw = parts.at(0).to_string();
            auto itType = colDefs().find(kw);
            if ( itType == colDefs().end() )
              return;
            if ( parts.size() < 2 )
              NCRYSTAL_THROW2(BadInput,errprefix<<"Missing value after keyword \""<<kw<<'"');
            auto valstr = parts.at(1);
            switch (itType->second) {
            case ColType::PosInt:
              {
                auto val = valstr.toInt64();
                if ( !val.has_value() || val.value()<1 || val.value() > std::numeric_limits<unsigned>::max() )
                  NCRYSTAL_THROW2(BadInput,errprefix<<"Key \""<<kw<<"\" has invalid value: \""<<valstr<<"\"");
                insertValue(errprefix,kw,m_posints,(unsigned)val.value());
              }
              break;
            case ColType::Dbl:
              {
                auto val = valstr.toDbl();
                if ( !val.has_value() || std::isnan(val.value()) || val.value() <= 0.0 || val.value() > 1e99 )
                  NCRYSTAL_THROW2(BadInput,errprefix<<"Key \""<<kw<<"\" has invalid value: \""<<valstr<<"\"");
                insertValue(errprefix,kw,m_dbls,val.value());
              }
              break;
            case ColType::Str:
              {
                if ( isOneOf(kw,"SPCGRP","TITLE") ) {
                  //need more than just the first "word" to be able to extract the spacegroup number.
                  const Span<const StrView> tmp( std::next(parts.begin()), parts.end() );
                  insertValue(errprefix,kw,m_strs,joinstr(tmp));
                } else {
                  insertValue(errprefix,kw,m_strs,valstr.to_string());
                }
              }
              break;
            default:
              nc_assert_always(false);
            };
          }

          void done(NC::Lazy::ParsedLazyData& res) {
            auto complainRequired = [](const char*key,const char*example)
            {
              NCRYSTAL_THROW2(BadInput,"Error in Lazy/Lau data: Missing required "<<key
                              <<" header field. Example of line to add to the header would be: "
                              <<"# "<<key<<" "<<example);
            };

            //Spacegroup info is special, it will preferably be taken from our
            //custom "spacegroup" field, but might optionally be extracted from a
            //header field "SPCGRP" where it is encoded (by e.g. cif2hkl?) like: "SPCGRP  P m -3 m [Number 221]"

            unsigned spacegroup(0);

            if ( this->hasKey("spacegroup") ) {
              spacegroup = this->getPosInt("spacegroup");
              if ( spacegroup < 1 || spacegroup > 230 )
                NCRYSTAL_THROW(BadInput,"Error in Lazy/Lau data: Invalid spacegroup value (must be number in 1..230).");
            } else {
              if ( !this->hasKey("SPCGRP") )
                complainRequired("spacegroup","225");
              int64_t SPCGRP_val ( 0 );
              std::string ss = lowerCase( this->getStr("SPCGRP") );
              strreplace( ss, "[", " [ " );
              strreplace( ss, "]", " ] " );
              auto sv = StrView(ss);
              auto i = sv.find( "number" );
              if ( i != StrView::npos ) {
                auto p = sv.substr( i + 6 ).split();
                if ( !p.empty() && p.at(0).toInt().has_value() )
                  SPCGRP_val = p.at(0).toInt().value();
              }
              if ( SPCGRP_val < 1 || SPCGRP_val > 230 )
                NCRYSTAL_THROW(BadInput,"Error in Lazy/Lau data: Could not extract valid spacegroup number from SPCGRP field"
                               " (must be number in 1..230, expected format like \"SPCGRP  P m -3 m [Number 221]\").");
              spacegroup = static_cast<unsigned>( SPCGRP_val );
            }
            nc_assert_always( spacegroup>=1 && spacegroup<=230 );

            double lattice_a = this->value_or("lattice_a",0.0);
            double lattice_b, lattice_c;
            if ( ! (lattice_a>0.0&&lattice_a<1e6) )
              NCRYSTAL_THROW(BadInput,"Error in Lazy/Lau data: Missing or invalid lattice_a value");
            lattice_b = this->value_or("lattice_b",0.0);
            lattice_c = this->value_or("lattice_c",0.0);
            checkAndCompleteLattice( spacegroup, lattice_a, lattice_b, lattice_c );

            double lattice_aa = this->value_or("lattice_aa",0.0);
            double lattice_bb = this->value_or("lattice_bb",0.0);
            double lattice_cc = this->value_or("lattice_cc",0.0);
            checkAndCompleteLatticeAngles( spacegroup, lattice_aa, lattice_bb, lattice_cc );

            this->updateIfNonZero("lattice_a",lattice_a);
            this->updateIfNonZero("lattice_b",lattice_b);
            this->updateIfNonZero("lattice_c",lattice_c);
            this->updateIfNonZero("lattice_aa",lattice_aa);
            this->updateIfNonZero("lattice_bb",lattice_bb);
            this->updateIfNonZero("lattice_cc",lattice_cc);

            std::vector<std::pair<const char*,const char*>> fields = {
              {"lattice_a","4.321"},
              {"lattice_b","4.321"},
              {"lattice_c","4.321"},
              {"lattice_aa","90"},
              {"lattice_bb","90"},
              {"lattice_cc","90"},
              {"column_h","1"},
              {"column_k","2"},
              {"column_l","3"},
            };

            NC::Optional<NC::DecodedChemForm> chemform;
            NC::Optional<unsigned> formulaPerUnitCell;
            if ( this->hasKey("nformula_per_unitcell") && !this->hasKey("formula") )
              NCRYSTAL_THROW(BadInput,"Error in Lazy/Lau data: nformula_per_unitcell should not be used without the formula field as well");

            if ( this->hasKey("formula") ) {
              chemform = decodeSimpleChemicalFormula( this->getStr("formula") );
              if ( !this->hasKey("nformula_per_unitcell") )
                complainRequired("nformula_per_unitcell","4");
              formulaPerUnitCell = this->getPosInt( "nformula_per_unitcell" );
            } else if ( this->hasKey("TITLE") ) {
              //Try to extract the formula from data like: "# TITLE   Al2 O3 [Cubic, Centric (-1 at origin)]"
              auto ss = this->getStr("TITLE");
              auto p1 = NC::StrView(ss).splitTrimmedNoEmpty('[');
              if ( !p1.empty() )
                chemform =  NC::tryDecodeSimpleChemicalFormula( p1.at(0).to_string() );
              if ( !chemform.has_value () )
                NCRYSTAL_THROW(BadInput,"Error in Lazy/Lau data: Could not extract valid chemical formula number from TITLE field"
                               " (expected format like \"TITLE Al2 O3 [whatever]\"). You should edit the file and insert dedicated \"formula\" field like \"#formula Al2O3\".");
              formulaPerUnitCell = 1;
            } else {
              complainRequired("formula","Al2O3");
            }

            for ( auto e : fields ) {
              if ( !this->hasKey(e.first) )
                complainRequired(e.first,e.second);
            }

            if ( !this->hasKey("column_F") && !this->hasKey("column_F2") )
              complainRequired("column_F2","4");

            std::string unitF2 = "barn";
            if ( this->hasKey("unit_F2") ) {
              if (!this->hasKey("column_F2"))
                NCRYSTAL_THROW2(BadInput,"Error in Lazy/Lau data: unit_F2 specified without column_F2 specified");
              if (this->hasKey("column_F"))
                NCRYSTAL_THROW2(BadInput,"Error in Lazy/Lau data: unit_F2 specified so column_F should not be specified");
              unitF2 = this->getStr("unit_F2");
              if ( !isOneOf(unitF2,"barn","fm^2"))
                NCRYSTAL_THROW2(BadInput,"Error in Lazy/Lau data: unit_F2 specified with forbidden"
                                " value. Allowed values are \"barn\" or \"fm^2\"");
              if ( unitF2 == "fm^2" )
                this->col_structfactsq_scalefactor = 0.01;
            }
            if ( this->hasKey("temperature") )
              res.temp = Temperature{ this->getDbl("temperature") };
            res.debye_temp = DebyeTemperature{ this->value_or("debye_temperature",300.0) };

            nc_assert_always( chemform.has_value() );
            res.chemform = chemform.value();
            std::uint64_t chemform_ntot = 0;
            for ( const auto& e : res.chemform )
              chemform_ntot += e.first;
            nc_assert( chemform_ntot >= 1 );
            auto& si = res.structInfo;
            si.lattice_a = lattice_a;
            si.lattice_b = lattice_b;
            si.lattice_c = lattice_c;
            si.alpha = lattice_aa;
            si.beta = lattice_bb;
            si.gamma = lattice_cc;
            si.spacegroup = spacegroup;

            nc_assert_always( formulaPerUnitCell.has_value() );
            auto natoms = chemform_ntot * formulaPerUnitCell.value();
            nc_assert_always( natoms < 60000 );
            si.n_atoms = static_cast<unsigned>(natoms);

            //Prepare stuff needed for processing the data lines:
            this->rec_lattice = getReciprocalLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                                         si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
            this->eqreflcalc.emplace(spacegroup);

            this->col_h = this->getPosInt("column_h") - 1;
            this->col_k = this->getPosInt("column_k") - 1;
            this->col_l = this->getPosInt("column_l") - 1;
            if ( this->hasKey("column_d") )
              this->col_d = this->getPosInt("column_d") - 1;
            if ( this->hasKey("column_j") )
              this->col_j = this->getPosInt("column_j") - 1;

            if ( this->hasKey("column_F2") ) {
              this->col_structfact_is_unsquared = false;
              this->col_structfact = this->getPosInt("column_F2") - 1;
            } else {
              nc_assert( this->hasKey("column_F") );
              this->col_structfact_is_unsquared = true;
              this->col_structfact = this->getPosInt("column_F") - 1;
            }
          }

          double dspacingFromHKL( int h, int k, int l ) const
          {
            nc_assert( rec_lattice.has_value() );
            return ::NCrystal::dspacingFromHKL( h, k, l, rec_lattice.value() );
          }
          std::size_t multFromHKL( int h, int k, int l ) const
          {
            nc_assert(eqreflcalc.has_value());
            return eqreflcalc.value().getEquivalentReflections(h,k,l).size() * 2;
          }

          unsigned col_h = 9999999;
          unsigned col_k = 9999999;
          unsigned col_l = 9999999;
          unsigned col_structfact = 9999999;
          Optional<unsigned> col_d, col_j;
          bool col_structfact_is_unsquared = false;
          double col_structfactsq_scalefactor = 1.0;
        private:
          std::map<std::string,unsigned> m_posints;
          std::map<std::string,double> m_dbls;
          std::map<std::string,std::string> m_strs;
          Optional<RotMatrix> rec_lattice;
          Optional<EqRefl> eqreflcalc;
        };

        ParsedHdr hdrData;
        unsigned ncols = 0;//set after first line - to check ncols does not change.
      };
    }
  }
}

NC::Lazy::ParsedLazyData NC::Lazy::parseLazyTextData( const TextData& td, const double& dcutoff )
{
  NC::Lazy::ParsedLazyData res;
  CollectedData cdat;
  auto parse_hdr_line = [&res,&cdat](const StrView& l)
  {
    constexpr auto errprefix = StrView::make("Error in Lazy/Lau data: ");
    cdat.hdrData.addLine(errprefix,l);
    nc_assert(l.startswith('#'));
    res.raw_header.push_back(l.to_string());
  };
  auto parse_line = [&res,&cdat,dcutoff](const StrView& line)
  {
    auto parts = line.split();
    auto& hdr = cdat.hdrData;
    if ( cdat.ncols==0 ) {
      //First data line:
      cdat.ncols = parts.size();
      auto maxcol = std::max(hdr.col_h,std::max(hdr.col_k,std::max(hdr.col_l,hdr.col_structfact)));
      if ( hdr.col_d.has_value() )
        maxcol = std::max(maxcol,hdr.col_d.value());
      if ( hdr.col_j.has_value() )
        maxcol = std::max(maxcol,hdr.col_j.value());
      if ( maxcol >= cdat.ncols )
        NCRYSTAL_THROW2(BadInput,"Invalid data (column_XX in header overflows actual number of data columns)");
    } else if ( cdat.ncols != parts.size() ) {
      NCRYSTAL_THROW2(BadInput,"Invalid data line (changed number of columns): "<<line);
    }

    auto hh = parts.at(hdr.col_h).toInt();
    auto kk = parts.at(hdr.col_k).toInt();
    auto ll = parts.at(hdr.col_l).toInt();
    auto fsq = parts.at(hdr.col_structfact).toDbl();
    Optional<double> dsp_listed;
    Optional<unsigned> mult_listed;
    if ( hdr.col_d.has_value() ) {
      auto tmp = parts.at(hdr.col_d.value()).toDbl();
      if ( !tmp.has_value() || !(tmp.value()>0.0 && tmp.value()<=1e6) )
        NCRYSTAL_THROW2(BadInput,"Bad dspacing value in data line: \""<<line<<"\"");
      dsp_listed = tmp.value();
    }
    if ( hdr.col_j.has_value() ) {
      auto tmp = parts.at(hdr.col_j.value()).toInt();
      if ( !tmp.has_value() || !(tmp.value()>=1 && tmp.value()<=10000000) )
        NCRYSTAL_THROW2(BadInput,"Bad multiplicity value in data line: \""<<line<<"\"");
      mult_listed = static_cast<unsigned>(tmp.value());
    }

    if ( hh.has_value() && kk.has_value() && ll.has_value() && fsq.has_value()
         && ncabs(hh.value()) <= 100000 && ncabs(kk.value()) <= 100000 && ncabs(ll.value()) <= 100000 ) {
      if ( hdr.col_structfact_is_unsquared )
        fsq.value() *= fsq.value();
      if (!(fsq.value()>=0) )
        NCRYSTAL_THROW2(BadInput,"Bad fsquared value in data line: \""<<line<<"\"");
      fsq.value() *= hdr.col_structfactsq_scalefactor;
      auto h = static_cast<int>(hh.value());
      auto k = static_cast<int>(kk.value());
      auto l = static_cast<int>(ll.value());
      const double dsp_calc = cdat.hdrData.dspacingFromHKL(h,k,l);
      if ( dsp_listed.has_value() ) {
        //Use the calculated value for higher precision, but check we are close
        //to the listed.
        constexpr double dsp_tol = 1e-2;
        if ( !floateq( dsp_listed.value(), dsp_calc, dsp_tol, dsp_tol ) )
          NCRYSTAL_THROW2(BadInput,"Listed dspacing value ("<<dsp_listed.value()
                          <<") does not match expected value ("<<dsp_calc
                          <<") in data line: \""<<line<<"\"");
      }

      if ( mult_listed.has_value() ) {
        //Verify that the value is not higher than that expected from the
        //symmetry group (which would be a sign of multiple sym. families in one
        //line => we can't use EqRefl to get the full list of normals, nor can
        //we verify the correct multiplicity).
        auto mult_calc = hdr.multFromHKL(h,k,l);
        auto mult = mult_listed.value();
        if ( mult > mult_calc )
          NCRYSTAL_THROW2(BadInput,"Listed multiplicity value ("<<mult
                          <<") is greater than expected value ("<<mult_calc
                          <<") in data line: \""<<line<<"\"");
        if ( mult != mult_calc && mult != 1 && mult != 2 && mult*2 != mult_calc )
          NCRYSTAL_THROW2(BadInput,"Listed multiplicity value "<<mult
                          <<" does not match expected value "<<mult_calc
                          <<" (nor is it 1, 2, or half of the expected value) in data line: \""<<line<<"\"");
      }
      if ( ( dcutoff<=0.0 || dsp_calc >= dcutoff ) && fsq.value() > 0.0 )
        res.hkllist.push_back( { fsq.value(), h, k, l } );
    } else {
      NCRYSTAL_THROW2(BadInput,"Invalid data line: "<<line);
    }
  };

  bool header_done(false);
  for ( auto& line_str : td ) {
    auto l = StrView(line_str).trimmed();
    if ( l.empty() )
      continue;//just an empty line (after trimming)
    if ( l.startswith('#') ) {
      if (header_done)
        continue;//just a comment, not part of the initial header
      parse_hdr_line(line_str);
    } else {
      header_done = true;
      auto n = l.find('#');
      if ( n!=StrView::npos)
        l = l.substr(0,n).trimmed();
      if ( l.empty() )
        continue;
      if ( cdat.ncols == 0 )
        cdat.hdrData.done(res);
      parse_line( l );
    }
  }

  if ( cdat.ncols == 0 )
    NCRYSTAL_THROW2(BadInput,"No data lines found in Lazy/Lau file");

  return res;
}

NC::InfoPtr NC::Lazy::buildInfo( const LazyCfgVars& cfg, const ParsedLazyData& data )
{
  auto hkllist = validateAndNormaliseHKLFsqList( data.structInfo.spacegroup, data.hkllist );
  if ( hkllist.size() != data.hkllist.size() ) {
    NCRYSTAL_MSG("Trimming and normalising provided HKL list from "
                 <<data.hkllist.size()<<" to "<<hkllist.size()<<" entries");
  } else {
    NCRYSTAL_MSG("Normalising provided HKL list with "
                 <<data.hkllist.size()<<" entries");
  }

  auto& si = data.structInfo;
  auto rec_lat = getReciprocalLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                          si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
  InfoBuilder::SinglePhaseBuilder builder;
  builder.dataSourceName = cfg.dataSourceName;

  HKLList out_hkllist;
  EqRefl sym(data.structInfo.spacegroup);
  //  out_hkllist.reserve(512);
  double dsp_min(kInfinity);
  for ( auto& e : hkllist ) {
    double dsp = dspacingFromHKL( e.h, e.k, e.l, rec_lat );
    if ( dsp < cfg.dcutoff || dsp > cfg.dcutoffup )
      continue;
    if ( out_hkllist.size() == decltype(out_hkllist)::nsmall+1 )
      out_hkllist.reserve_hint( 512 );
    out_hkllist.emplace_back();
    auto& hi = out_hkllist.back();
    hi.dspacing = dsp;
    dsp_min = ncmin( dsp_min, hi.dspacing );
    auto eqrefl = sym.getEquivalentReflections(e.h,e.k,e.l);
    hi.multiplicity = eqrefl.size() * 2;
    hi.hkl = eqrefl.front();
    hi.fsquared = e.fsquared;
  }
  out_hkllist.shrink_to_fit();
  PairDD dsp_range{ ( cfg.dcutoff <= 0.0 ? dsp_min : cfg.dcutoff ), cfg.dcutoffup };
  if ( dsp_range.first >= cfg.dcutoffup )
    dsp_range.first = 0.999 * dsp_range.second;

  builder.hklPlanes.emplace();
  builder.hklPlanes.value().dspacingRange = dsp_range;
  builder.hklPlanes.value().source = std::move(out_hkllist);

  builder.unitcell.emplace();
  builder.unitcell.value().structinfo = std::move(data.structInfo);

  if ( cfg.temp.has_value() && data.temp.has_value() && cfg.temp.value() != data.temp.value() )
    NCRYSTAL_THROW2(BadInput,"Requested T="<<cfg.temp.value()<<" in input \""
                    <<cfg.dataSourceName<<"\" which only supports T="<<data.temp.value());
  if ( data.temp.has_value() )
    builder.temperature = data.temp.value();
  else if ( cfg.temp.has_value() )
    builder.temperature = cfg.temp.value();
  else
    builder.temperature = Temperature{ 293.15 };

  builder.dynamics.emplace();
  auto& dyninfos = builder.dynamics.value();
  //dyninfos.reserve(data.chemform.size());

  //Setup atomdb (support atomdb cfg keyword):
  const bool cfg_atomdb_nodefs = ( !cfg.atomdb.empty() && cfg.atomdb.at(0).size()==1
                                   && cfg.atomdb.at(0).at(0)=="nodefaults");
  AtomDBExtender atomdb(!cfg_atomdb_nodefs);//<-- the database
  if ( !cfg.atomdb.empty() ) {
    for (unsigned i = (cfg_atomdb_nodefs?1:0); i<cfg.atomdb.size(); ++i)
      atomdb.addData(cfg.atomdb.at(i),AtomDBExtender::latest_version);
  }

  //Total number in chemform:
  std::uint64_t chemform_ntot = 0;
  for ( const auto& e : data.chemform )
    chemform_ntot += e.first;
  nc_assert( chemform_ntot >= 1 );

  //Setup dyninfolist:
  unsigned atom_idx_counter(0);
  for ( const auto& e : data.chemform ) {
    if (!e.second.isElement())
      NCRYSTAL_THROW2(BadInput,"Non-natural element found in chemical formula in input \""
                      <<cfg.dataSourceName<<"\"");
    const double fraction = ( double(e.first) / double(chemform_ntot) );
    auto atomdatasp = atomdb.lookupAtomData( elementZToName( e.second.Z() ) );
    IndexedAtomData iad{ atomdatasp, AtomIndex{atom_idx_counter++} };
    dyninfos.emplace_back( std::make_unique<DI_VDOSDebye>(fraction,
                                                          iad,
                                                          builder.temperature.value(),
                                                          data.debye_temp ) );
  }

  dyninfos.shrink_to_fit();
  return buildInfoPtr(std::move(builder));
}

NC::InfoPtr NC::Lazy::buildInfoFromLazyData( const FactImpl::InfoRequest& req )
{
  LazyCfgVars cfg;
  if ( req.get_temp().dbl() != -1.0 )
    cfg.temp = req.get_temp();
  cfg.dcutoff = req.get_dcutoff();
  cfg.dcutoffup = req.get_dcutoffup();
  cfg.dataSourceName = req.dataSourceName();
  cfg.atomdb = req.get_atomdb_parsed();
  auto parsedData = parseLazyTextData( req.textData() );
  return buildInfo( cfg, parsedData );
}

std::vector<NC::Lazy::HKLFsq>
NC::Lazy::validateAndNormaliseHKLFsqList( int spacegroup,
                                          const std::vector<HKLFsq>& orig )
{
  nc_assert_always(spacegroup>=1&&spacegroup<=230);
  std::vector<HKLFsq> l = orig;//NB: in principle this copying is a bit
                               //wasteful, but lazy files are not that important
                               //anyway for us.

  std::stable_sort(l.begin(), l.end(),
                   [](const HKLFsq& a, const HKLFsq & b) -> bool
                   {
                     return a.fsquared > b.fsquared;
                   });

  std::vector<HKLFsq> res;
  res.reserve(512);

  EqRefl sym(spacegroup);

  using HKLPt = std::tuple<int,int,int>;
  std::vector<HKLPt> families_seen;
  families_seen.reserve(1024);

  auto it = l.begin();
  auto itAllE = l.end();

  enum class ListType { All, OnePerGroup, HalfInEachGroup, Unknown };

  ListType listTypeSeen = ListType::Unknown;
  while ( it != itAllE ) {
    auto itE = std::find_if(it,itAllE,[&it](const HKLFsq& a){ return a.fsquared != (it->fsquared); });
    //All in [it,itE) have same fsquared. Sort so all that are sym equiv. with
    //the first entry goes first. If only one, that entry will be expanded
    //automatically to the full family. If multiple, we require the correct
    //family. In either case we move the it pointer up to the first unused entry
    //and proceed.
    auto eq_hkls = sym.getEquivalentReflections(it->h, it->k, it->l);
    auto inFamily = [&eq_hkls]( const HKLFsq& a ) { return eq_hkls.contains(HKL{ a.h, a.k, a.l }); };
    nc_assert(inFamily(*it));

    std::stable_sort( it, itE,
                     [&inFamily](const HKLFsq& a, const HKLFsq & b) -> bool
                     {
                       nc_assert( a.fsquared == b.fsquared );
                       bool a_in_family = inFamily(a);
                       if ( a_in_family == inFamily(b) )
                         return false;//"equal"
                       return a_in_family;
                     });
    auto itFamilyE = it;
    while ( itFamilyE != itE && inFamily(*itFamilyE) )
      ++itFamilyE;
    auto nfam = (std::size_t)std::distance(it,itFamilyE);

    ListType listType = ListType::Unknown;
    if ( nfam == 1 )
      listType = ListType::OnePerGroup;
    else if ( nfam == 2*eq_hkls.size() )
      listType = ListType::HalfInEachGroup;
    else if ( nfam == eq_hkls.size() )
      listType = ListType::All;
    if ( listType == ListType::Unknown )
      NCRYSTAL_THROW2(BadInput,"HKL list is not consistent with spacegroup. The hkl=("
                      <<it->h<<", "<<it->k<<", "<<it->l<<") value should have "
                      <<2*eq_hkls.size()<<" equivalent planes (or possibly just one or "
                      "half of these specified), but only found "
                      <<nfam<<" of these in the list");
    if ( listTypeSeen == ListType::Unknown )
      listTypeSeen = listType;
    else if ( listTypeSeen != listType )
      NCRYSTAL_THROW2(BadInput,"Inconsistent specification of HKL. The family with the hkl=("
                      <<it->h<<", "<<it->k<<", "<<it->l<<") use a different specification than"
                      " other entries in the file (entries should have NFAM, NFAM/2, or 1"
                      " equivalent planes listed - but this must be the same in the entire file");

    //Use the first entry in eq_hkls to identify the family, but be aware that
    //the sign of all entries in eq_hkls will be flipped if the HKL point passed
    //to getEquivalentReflections(..) was flipped. Thus, for consistency (and to
    //ensure the the check below involving "families_seen" works), we ensure
    //that the first non-zero hkl coordinate is not negative.:
    int fam_h(eq_hkls.begin()->h);
    int fam_k(eq_hkls.begin()->k);
    int fam_l(eq_hkls.begin()->l);

    //Just put the first entry into the list:
    families_seen.emplace_back( fam_h, fam_k, fam_l );
    res.push_back( HKLFsq{ it->fsquared, fam_h, fam_k, fam_l } );
    //Skip family:
    it = itFamilyE;
  }

  //Guard against half of the family members having different Fsq factor than
  //the others (which would fool the checks above):
  std::sort(families_seen.begin(),families_seen.end());
  auto itAdj = std::adjacent_find(families_seen.begin(),families_seen.end());
  if ( itAdj != families_seen.end() ) {
    auto eq_hkls = sym.getEquivalentReflections( std::get<0>(*itAdj), std::get<1>(*itAdj), std::get<2>(*itAdj) );
    std::ostringstream ss;
    for ( auto& ee : eq_hkls ) {
      ss<<" ("<<ee.h<<','<<ee.k<<','<<ee.l<<")"
        <<" ("<<-ee.h<<','<<-ee.k<<','<<-ee.l<<")";
    }
    NCRYSTAL_THROW2(BadInput,"HKL list is not consistent with spacegroup."
                    " Members of the following symmetry family of HKL planes have"
                    " been found with differing F^2 values:"<<ss.str());
  }
  res.shrink_to_fit();
  return res;

}
