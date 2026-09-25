////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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


#include "NCBrowseQuery.hh"
#include "NCrystal/factories/NCDataSources.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/misc/NCCompositionUtils.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <map>
#include <set>
#include <sstream>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace BrowseQuery {
    namespace {

      using FileListEntry = DataSources::FileListEntry;

      //Rank for sorting (higher first), consistent with the Python
      //NCrystal.browseFiles function:
      std::int64_t priorityRank( const Priority& p )
      {
        if ( !p.canServiceRequest() )
          return -2;
        if ( p.needsExplicitRequest() )
          return -1;
        return static_cast<std::int64_t>( p.priority() );
      }

      struct Entry {
        FileListEntry fle;
        bool hidden;
      };

      //All browsable entries, sorted by priority, factory, source and
      //name. Entries are hidden if an earlier entry has the same name.
      std::vector<Entry> allEntries()
      {
        auto v = DataSources::listAvailableFiles();
        std::stable_sort( v.begin(), v.end(),
                          []( const FileListEntry& a,
                              const FileListEntry& b )
        {
          auto ra = priorityRank(a.priority);
          auto rb = priorityRank(b.priority);
          if ( ra != rb )
            return ra > rb;
          if ( a.factName != b.factName )
            return a.factName < b.factName;
          if ( a.source != b.source )
            return a.source < b.source;
          return a.name < b.name;
        } );
        std::vector<Entry> res;
        res.reserve( v.size() );
        std::set<std::string> seen;
        for ( auto& e : v ) {
          bool hidden = !seen.insert( e.name ).second;
          res.push_back( Entry{ std::move(e), hidden } );
        }
        return res;
      }

      std::string fullKey( const FileListEntry& e )
      {
        return e.factName + "::" + e.name;
      }

      void streamPriority( std::ostream& os, const Priority& p )
      {
        if ( !p.canServiceRequest() )
          streamJSON( os, "Unable" );
        else if ( p.needsExplicitRequest() )
          streamJSON( os, "OnlyOnExplicitRequest" );
        else
          streamJSON( os, static_cast<std::uint64_t>( p.priority() ) );
      }

      //Initial comment lines of NCMAT data (those before the first
      //section), without the '#' and dedented. Mirrors the Python
      //function _extractInitialHeaderCommentsFromNCMATData:
      VectS headerComments( const TextData& td )
      {
        VectS res;
        bool first = true;
        for ( auto& line : td ) {
          if ( first ) {
            first = false;
            continue;
          }
          StrView sv( line );
          auto t = sv.ltrimmed();
          if ( t.startswith('@') )
            break;
          if ( t.startswith('#') )
            res.push_back( t.substr(1).rtrimmed().to_string() );
        }
        auto canDedent = [&res]()
        {
          bool anyNonEmpty = false;
          for ( auto& e : res ) {
            if ( e.empty() )
              continue;
            anyNonEmpty = true;
            if ( e.front() != ' ' )
              return false;
          }
          return anyNonEmpty;
        };
        while ( canDedent() ) {
          for ( auto& e : res )
            if ( !e.empty() )
              e.erase( 0, 1 );
        }
        return res;
      }

      void collectSinglePhases( const Info& info,
                                std::vector<const Info*>& out )
      {
        if ( !info.isMultiPhase() ) {
          out.push_back( &info );
          return;
        }
        for ( auto& ph : info.getPhases() )
          collectSinglePhases( ph.second, out );
      }

      const char * dynInfoTypeName( const DynamicInfo& di )
      {
        if ( dynamic_cast<const DI_Sterile*>(&di) )
          return "sterile";
        if ( dynamic_cast<const DI_FreeGas*>(&di) )
          return "freegas";
        if ( dynamic_cast<const DI_ScatKnlDirect*>(&di) )
          return "scatknl";
        if ( dynamic_cast<const DI_VDOS*>(&di) )
          return "vdos";
        if ( dynamic_cast<const DI_VDOSDebye*>(&di) )
          return "vdosdebye";
        return "other";
      }

      //Physics properties of loaded material as a JSON dict:
      std::string infoPropsJSON( const Info& info )
      {
        std::vector<const Info*> sps;
        collectSinglePhases( info, sps );
        bool crystalline = false;
        std::set<std::string> ditypes;
        for ( auto sp : sps ) {
          if ( sp->isCrystalline() )
            crystalline = true;
          for ( auto& di : sp->getDynamicInfoList() )
            ditypes.insert( dynInfoTypeName( *di ) );
        }
        Optional<unsigned> sg, natoms;
        if ( !info.isMultiPhase() && info.hasStructureInfo() ) {
          auto& si = info.getStructureInfo();
          if ( si.spacegroup )
            sg = si.spacegroup;
          natoms = si.n_atoms;
        }
        auto bd = CompositionUtils::createFullBreakdown(
          info.getComposition(), nullptr,
          CompositionUtils::PreferNaturalElements );
        std::ostringstream os;
        os << "{\"composition\":"
           << CompositionUtils::fullBreakdownToJSON( bd );
        streamJSONDictEntry( os, "xsect_absorption",
                             info.getXSectAbsorption().dbl() );
        streamJSONDictEntry( os, "xsect_free",
                             info.getXSectFree().dbl() );
        streamJSONDictEntry( os, "density", info.getDensity().dbl() );
        streamJSONDictEntry( os, "numberdensity",
                             info.getNumberDensity().dbl() );
        Optional<double> temp;
        if ( info.hasTemperature() )
          temp = info.getTemperature().dbl();
        streamJSONDictEntry( os, "temperature", temp );
        streamJSONDictEntry( os, "stateofmatter",
                             Info::toString( info.stateOfMatter() ) );
        streamJSONDictEntry( os, "crystalline", crystalline );
        streamJSONDictEntry( os, "spacegroup", sg );
        streamJSONDictEntry( os, "natoms_unitcell", natoms );
        streamJSONDictEntry( os, "dyninfo_types",
                             VectS( ditypes.begin(), ditypes.end() ) );
        const auto nphases = static_cast<unsigned>(
          info.isMultiPhase() ? info.getPhases().size() : 1 );
        streamJSONDictEntry( os, "nphases", nphases,
                             JSONDictPos::LAST );
        return os.str();
      }

      std::string errMsg( const std::exception& e )
      {
        std::string s = e.what();
        return s.empty() ? std::string("<unknown error>") : s;
      }

      //JSON dict for a single entry:
      std::string entryJSON( const Entry& entry, bool cheap )
      {
        const auto& e = entry.fle;
        const std::string key = fullKey( e );
        std::ostringstream os;
        streamJSONDictEntry( os, "name", e.name, JSONDictPos::FIRST );
        streamJSONDictEntry( os, "fullkey", key );
        streamJSONDictEntry( os, "factory", e.factName );
        streamJSONDictEntry( os, "source", e.source );
        os << ",\"priority\":";
        streamPriority( os, e.priority );
        streamJSONDictEntry( os, "hidden", entry.hidden );
        Optional<std::string> error;
        try {
          auto td = FactImpl::createTextData( TextDataPath( key ) );
          streamJSONDictEntry( os, "datatype", td->dataType() );
          os << ",\"comments\":";
          if ( td->dataType() == "ncmat" )
            streamJSON( os, headerComments( *td ) );
          else
            streamJSON( os, json_null_t{} );
          if ( !cheap ) {
            auto info = FactImpl::createInfo( MatCfg( key ) );
            os << ",\"info\":" << infoPropsJSON( *info );
          }
        } catch ( std::exception& ex ) {
          error = errMsg( ex );
        } catch ( ... ) {
          error = std::string("<unknown error>");
        }
        if ( error.has_value() )
          streamJSONDictEntry( os, "error", error.value() );
        os << '}';
        return os.str();
      }

      template<class TFactList>
      VectS factNames( const TFactList& factlist )
      {
        VectS res;
        for ( auto& f : factlist )
          res.emplace_back( f->name() );
        return res;
      }

      unsigned parseUInt( StrView sv, const char * what )
      {
        auto v = sv.trimmed().toInt();
        if ( !v.has_value() || v.value() < 0 )
          NCRYSTAL_THROW2( BadInput, "Invalid " << what
                           << " in browsedb query: \"" << sv << '"' );
        return static_cast<unsigned>( v.value() );
      }
    }
  }
}

void NC::BrowseQuery::browseDB( std::ostream& os,
                                const std::vector<StrView>& args_orig )
{
  auto args = args_orig;
  if ( args.empty() ) {
    //Dict of TextData factories and number of entries:
    std::map<std::string,unsigned> counts;
    for ( auto& e : DataSources::listAvailableFiles() )
      ++counts[ e.factName ];
    os << '{';
    bool first = true;
    for ( auto& f : FactImpl::getTextDataFactoryList() ) {
      if ( !first )
        os << ',';
      first = false;
      streamJSON( os, f->name() );
      os << ':';
      auto it = counts.find( f->name() );
      streamJSON( os, it == counts.end() ? 0u : it->second );
    }
    os << '}';
    return;
  }

  const std::string factname = args.front().trimmed().to_string();
  args.erase( args.begin() );
  bool cheap = false;
  if ( !args.empty() && args.back().trimmed() == "cheap" ) {
    cheap = true;
    args.pop_back();
  }
  unsigned ichunk = 0;
  unsigned nchunks = 1;
  if ( args.size() == 2 ) {
    ichunk = parseUInt( args.at(0), "chunk index I" );
    nchunks = parseUInt( args.at(1), "number of chunks N" );
  } else if ( !args.empty() ) {
    NCRYSTAL_THROW( BadInput, "Invalid browsedb query (usage:"
                    " [\"util\",\"browsedb\",FACTNAME,"
                    "(I,N,)(\"cheap\")])" );
  }
  if ( !( nchunks >= 1 && ichunk < nchunks ) )
    NCRYSTAL_THROW2( BadInput, "Invalid chunk specification in browsedb"
                     " query (must have 0<=I<N): I=" << ichunk
                     << ", N=" << nchunks );
  if ( !FactImpl::hasTextDataFactory( factname ) )
    NCRYSTAL_THROW2( BadInput, "Unknown TextData factory in browsedb"
                     " query: \"" << factname << '"' );

  std::vector<Entry> entries;
  for ( auto& e : allEntries() )
    if ( e.fle.factName == factname )
      entries.push_back( std::move(e) );

  //The I'th of N contiguous chunks:
  const std::size_t n = entries.size();
  const std::size_t ibegin = ( n * ichunk ) / nchunks;
  const std::size_t iend = ( n * ( ichunk + 1 ) ) / nchunks;

  VectS results( iend - ibegin );
  {
    //NB: Jobs run in parallel only if factory threads are enabled:
    FactoryJobs jobs;
    for ( auto i : ncrange( ibegin, iend ) ) {
      jobs.queue( [i,ibegin,cheap,&entries,&results]()
      {
        vectAt( results, i - ibegin )
          = entryJSON( vectAt( entries, i ), cheap );
      } );
    }
    jobs.waitAll();
  }
  os << '[';
  for ( auto i : ncrange( results.size() ) ) {
    if ( i )
      os << ',';
    os << vectAt( results, i );
  }
  os << ']';
}

void NC::BrowseQuery::browseFactories( std::ostream& os )
{
  streamJSONDictEntry( os, "textdata",
                       factNames( FactImpl::getTextDataFactoryList() ),
                       JSONDictPos::FIRST );
  streamJSONDictEntry( os, "info",
                       factNames( FactImpl::getInfoFactoryList() ) );
  streamJSONDictEntry( os, "scatter",
                       factNames( FactImpl::getScatterFactoryList() ) );
  auto absfacts = factNames( FactImpl::getAbsorptionFactoryList() );
  streamJSONDictEntry( os, "absorption", absfacts, JSONDictPos::LAST );
}
