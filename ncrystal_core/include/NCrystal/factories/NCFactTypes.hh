#ifndef NCrystal_FactTypes_hh
#define NCrystal_FactTypes_hh

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

#include "NCrystal/text/NCTextData.hh"

//////////////////////////////////////////////
// Various types associated with factories. //
//////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  class NCRYSTAL_API TextDataPath {
  public:

    /////////////////////////////////////////////////////////////////////////////////////
    // Generalised path initialised with the same syntax as filenames in MatCfg
    // strings. This is an absolute, relative, or virtual path to a file,
    // optionally prefixed with the name of a TextData factory separated by
    // "::". The factory name is intended for usage in the global TextData
    // factories (cf. NCFactImpl.hh) allowing a convenient way for users to
    // request specific factories (e.g. "relpath::Al_sg225.ncmat" will select
    // the "relpath" factory which only looks in the current directory for
    // files, and "stdlib::Al_sg225.ncmat" will only provide files from the
    // standard NCrystal library).
    //
    // NB: Relative paths will be interpreted relatively to the current working
    //     directory at the point of usage (which is usually right away, when
    //     passing strings to the MatCfg constructor or the global createTextData
    //     factory function.
    //
    // A few special rules are imposed by the constructor (will throw if invalid input):
    //
    //  1) Path must be non empty
    //  2) Paths starting with "~/" will have the "~" replaced with the value of
    //     the users home directory (std::getenv("HOME")). Other occurances of
    //     '~' are not allowed.
    //  3) Factory names (if provided) must be alpha-numeric+"_-" without spaces.

    TextDataPath( const std::string& );
    TextDataPath( const char * );
    ncconstexpr17 const std::string& path() const noexcept { return m_path; }
    ncconstexpr17 const std::string& fact() const noexcept { return m_fact; }
    std::string toString() const;//as "fact::path" if fact set, else "path"
    bool operator<( const TextDataPath& ) const noexcept;
    bool operator==( const TextDataPath& ) const noexcept;
    bool operator!=( const TextDataPath& ) const noexcept;
  private:
    std::string m_path, m_fact;
  };
  std::ostream& operator<< ( std::ostream& , const TextDataPath& );

  class NCRYSTAL_API Priority final {

    /////////////////////////////////////////////////////////////////////////////
    // Priority value returned by FactImpl query(..) methods. Will typically   //
    // be returned like:                                                       //
    //                                                                         //
    //   "return Priority::Unable;"  or "return Priority{999};"                //
    //                                                                         //
    // Where 999 is the chosen priority (a higher number means a higher        //
    // priority}. Note that core NCrystal factories always return a priority   //
    // value of 100 if applicable (so they specifically return either          //
    // Priority::Unable or Priority{100}), so externally developed plugins     //
    // should return a higher value than that in order to beat those.          //
    //                                                                         //
    // Finally, a factory can also return Priority::OnlyOnExplicitRequest. In  //
    // that case, it is considered able to service the request, but will never //
    // be asked to do so unless the factory is explicitly by the user (by      //
    // specifying the factory name to the infofactory, scatfactory, or         //
    // absnfactory cfg parameter).                                             //
    /////////////////////////////////////////////////////////////////////////////

  public:
    enum SpecialCase { Unable, OnlyOnExplicitRequest };
    constexpr Priority( SpecialCase ) noexcept;
    ncconstexpr17 Priority( const uint_fast32_t priority_value );
    constexpr bool canServiceRequest() const;
    constexpr bool needsExplicitRequest() const;
    constexpr uint_fast32_t priority() const;
    constexpr Priority( const Priority& ) noexcept = default;
    ncconstexpr17 Priority& operator=( const Priority& ) = default;
    constexpr bool operator==( const Priority& ) const noexcept;
    constexpr bool operator!=( const Priority& ) const noexcept;
  private:
    uint_fast32_t m_value;
  };
  NCRYSTAL_API std::ostream& operator<< ( std::ostream& , Priority );

  class NCRYSTAL_API TextDataSource final {
  public:
    //Text data source, either in-memory or as a path to an on-disk file.  A
    //dataType can optionally be specified (if not, an error will ensue later if
    //NCrystal is not able to determine the format though other means, such as
    //looking at the actual data -- which always works for NCMAT data --, or the
    //file-name extension (e.g. a file named Foo.bar is assumed to have data
    //type "bar").
    static TextDataSource createFromOnDiskPath( std::string, std::string dataType = {},
                                                std::string suggestedDataSourceName = {} );
    static TextDataSource createFromInMemData( RawStrData, std::string dataType = {},
                                               std::string suggestedDataSourceName = {} );
    const Variant<std::string,RawStrData>&& data() &&;
    const std::string& dataType() const noexcept;

    //Suggested data source name (just leave empty if basename(..) on the
    //associated textdatapath.path() is enough, this method is mostly intended
    //for "quick-factories" where the names are actually filenames (virtual or
    //on-disk).
    const std::string& suggestedDataSourceName() const noexcept { return m_sdsn; }
  private:
    TextDataSource() = default;
    Variant<std::string,RawStrData> m_data;
    std::string m_dataType;
    std::string m_sdsn;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline constexpr Priority::Priority( SpecialCase sc ) noexcept : m_value( sc==Unable ? 0 : 0x3 ) {}
  inline constexpr bool Priority::canServiceRequest() const { return m_value & 0x1; }
  inline constexpr bool Priority::needsExplicitRequest() const { return m_value & 0x2; }
  inline constexpr uint_fast32_t Priority::priority() const { return m_value / 4; }
  inline constexpr bool Priority::operator==( const Priority& o ) const noexcept { return m_value == o.m_value; }
  inline constexpr bool Priority::operator!=( const Priority& o ) const noexcept { return m_value != o.m_value; }

  inline ncconstexpr17 Priority::Priority( const uint_fast32_t ppp ) : m_value{ 4*ppp + 0x1 }
  {
    if (!(ppp>=1&&ppp<=1000000000))
      NCRYSTAL_THROW(BadInput,"Priority must be in range 1-1000000000");
  }

  inline std::ostream& operator<< ( std::ostream& os , Priority p )
  {
    if ( !p.canServiceRequest() )
      return os << "Priority{Unable}";
    if ( p.needsExplicitRequest() )
      return os << "Priority{OnlyOnExplicitRequest}";
    return os << "Priority{"<< p.priority()<<"}";
  }

  inline TextDataSource TextDataSource::createFromInMemData( RawStrData data, std::string dt, std::string sdsn )
  {
    TextDataSource t;
    t.m_data = std::move(data);
    t.m_dataType = std::move(dt);
    t.m_sdsn = std::move(sdsn);
    return t;
  }

  inline TextDataSource TextDataSource::createFromOnDiskPath( std::string path, std::string dt, std::string sdsn )
  {
    TextDataSource t;
    t.m_data = std::move(path);
    t.m_dataType = std::move(dt);
    t.m_sdsn = std::move(sdsn);
    return t;
  }

  inline const Variant<std::string,RawStrData>&& TextDataSource::data() &&
  {
    return std::move(m_data);
  }

  inline const std::string& TextDataSource::dataType() const noexcept
  {
    return m_dataType;
  }

  inline bool TextDataPath::operator<( const TextDataPath& o ) const noexcept
  {
    if ( m_path != o.m_path )
      return m_path < o.m_path;
    return m_fact < o.m_fact;
  }

  inline bool TextDataPath::operator==( const TextDataPath& o ) const noexcept
  {
    return m_path == o.m_path && m_fact == o.m_fact;
  }

  inline bool TextDataPath::operator!=( const TextDataPath& o ) const noexcept
  {
    return m_path != o.m_path || m_fact != o.m_fact;
  }

  inline std::string TextDataPath::toString() const
  {
    std::string s;
    s.reserve( 2 + m_path.size() + m_fact.size() );
    if (!m_fact.empty()) {
      s += m_fact;
      s += "::";
    }
    s += m_path;
    return s;
  }

  inline TextDataPath::TextDataPath( const char * c ) : TextDataPath( std::string(c) ) {}

  inline std::ostream& operator<< ( std::ostream& os, const TextDataPath& fp )
  {
    if ( !fp.fact().empty() )
      os <<  fp.fact() << "::";
    os << fp.path();
    return os;
  }
}

#endif
