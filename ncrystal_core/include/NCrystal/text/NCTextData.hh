#ifndef NCrystal_TextData_hh
#define NCrystal_TextData_hh

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

#include "NCrystal/core/NCVariant.hh"
#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  ////////////////////////////////////////////
  // TextData objects and associated types. //
  ////////////////////////////////////////////

  class TextData;
  using TextDataSP = shared_obj<const TextData>;
  using OptionalTextDataSP = std::shared_ptr<const TextData>;

  class NCRYSTAL_API TextDataUID {
  public:
    constexpr TextDataUID();//unset
    static TextDataUID createNewUniqueIDValue();
    constexpr bool isUnset() const noexcept;
    constexpr bool operator<( const TextDataUID& ) const noexcept;
    constexpr bool operator==( const TextDataUID& ) const noexcept;
    constexpr bool operator!=( const TextDataUID& ) const noexcept;
    uint64_t value() const noexcept;//==0 when unset
    void set( TextDataUID );
  private:
    UniqueIDValue m_value;
    TextDataUID(UniqueIDValue);
  };

  class RawStrData {
  public:
    //Low-level string data object which does not have any sort of meta-data. Is
    //cheap to copy around and should not have any life-time issues, since it
    //internally wraps either a long-lived const char pointer or a
    //shared_obj<std::string> object. The data is always contiguous and
    //null-terminated (thus can only store data with null-free encodings like
    //ASCII or UTF-8). Constructors taking std::string also optionally takes a
    //srcdescr parameter which is used for more meaningful error messages in
    //case of unsupported encodings.
    struct static_data_ptr_t {};
    RawStrData( static_data_ptr_t, const char * ) noexcept;
    RawStrData( std::string&&, const char * srcdescr = nullptr );
    RawStrData( shared_obj<std::string>, const char * srcdescr = nullptr );
    RawStrData( std::string&&, const DataSourceName& );
    RawStrData( shared_obj<std::string>, const DataSourceName& );

    RawStrData( const RawStrData& ) = default;
    RawStrData( RawStrData&& ) = default;
    RawStrData& operator=( const RawStrData& ) = default;
    RawStrData& operator=( RawStrData&& ) = default;

    ncconstexpr17 const char * begin() const noexcept { return m_b; }
    ncconstexpr17 const char * end() const noexcept { return m_e; }

    //Calculate simple checksum of the data:
    uint64_t calcCheckSum() const;

    //Check if kept data is byte-wise identical to data in argument:
    bool hasSameContent( const char* dataBegin, const char* dataEnd ) const;
    bool hasSameContent( const std::string& ) const;
    bool hasSameContent( const RawStrData& ) const;

    //Expose checksum algorithm for usage without RawStrData objects:
    static uint64_t checkSumFromRawStringData(const char*begin, const char*end);

  private:
    const char *m_b, *m_e;
    optional_shared_obj<std::string> m_s;
  };

  class NCRYSTAL_API TextData : private MoveOnly {
  public:

    // Text data accessible line by line, with associated meta-data. This always
    // include a UID (useful for comparison and downstream caching purposes) and
    // the data type (e.g. "ncmat"). Optionally available is the last known
    // on-disk path to a file with the same content, which might be useful in
    // case the data needs to be passed to 3rd party software which can only
    // work with physical files.
    //
    // Text data objects are easily line-iterable, easily providing lines
    // (without newline characters): for( auto& line : mytextdata ) {...}.
    // Of course, the raw underlying data buffer can also be accessed if needed.
    //
    // The raw data must be ASCII or UTF-8 text, with line endings \n=CR=0x0A
    // (Unix) or \r\n=LF+CR=0x0D0A (Windows/dos). Other encodings might work
    // only if 0x00, 0x0A, 0x0D bytes do not occur in them outside of line
    // endings.
    //
    // Notice that ancient pre-OSX Mac line-endings \r=LF=0x0D are not
    // supported, and iterators will actually throw an error upon encountering
    // them. This is done on purpose, since files with \r on unix might hide
    // content when inspected in a terminal can be either confusing, a potential
    // security issue, or both.

    //When dereferenced, iterators provide each line as a std::string (with
    //newline characters removed and properly null terminated content):
    class Iterator;
    Iterator begin() const;
    Iterator end() const;
    Iterator cbegin() const { return begin(); }
    Iterator cend() const { return end(); }

    //Raw access to underlying data:
    const RawStrData& rawData() const noexcept;
    std::string rawDataCopy() const;

    //Construct by combining raw string data with optional meta-data, which are:
    //   1) A data type (e.g. "ncmat", "lau", ...) describing data format.
    //   2) A data source name (only to be used for e.g. error messages). This
    //      flexible description can for instance be a filename.
    //   3) A last known physical location (as resolved absolute on-disk
    //      path). The indicated file should have content which is byte-to-byte
    //      identical to the raw string data provided in a separate argument,
    //      and is intended to be used solely for software which needs to deal
    //      in on-disk files.

    struct DataType { std::string value; };
    struct LastKnownOnDiskAbsPath { std::string value; };

    TextData( RawStrData, DataType,
              Optional<DataSourceName> = NullOpt,
              Optional<LastKnownOnDiskAbsPath> = NullOpt );

    ncconstexpr17 const TextDataUID& dataUID() const noexcept;

    //Data type ("ncmat", "lau", ...):
    const std::string& dataType() const noexcept;

    //Description (for debug output and error messages). If not provided in the
    //constructor, an alternative will be constructed. Should ideally be just
    //the filename, presented in essentially the same form originally entered by
    //users.
    const DataSourceName& dataSourceName() const noexcept { return m_dsn; }

    //Some consumers might only be able to deal with physical files. For such
    //files, we keep the last known on-disk path around. Code reading from such
    //an on-disk file directly should call verifyOnDiskFileUnchanged() just
    //before and after the code reading is carried out. This verification will
    //throw an error in case the file content are not as expected:
    const Optional<std::string>& getLastKnownOnDiskLocation() const noexcept;
    void verifyOnDiskFileUnchanged() const;

    //Check if meta-data (everything except the RawStrData and UID) is
    //identical:
    bool hasIdenticalMetaData( const TextData& ) const noexcept;

    //Full iterator declaration and internals:
    class Iterator {
    public:
      //Readonly iterator keep std::string buffer with current line. Only
      //support postfix increment, as prefix increment would cost an extra
      //string allocation.
      using value_type = std::string;
      Iterator& operator++();
      ncconstexpr17 const value_type* operator->() const noexcept;
      ncconstexpr17 const value_type& operator*() const noexcept;
      ncconstexpr17 bool operator==(const Iterator&) const noexcept;
      ncconstexpr17 bool operator!=(const Iterator&) const noexcept;
      ncconstexpr17 bool operator<(const Iterator&) const noexcept;
      Iterator( const Iterator& );
      Iterator& operator=( const Iterator& );
      Iterator( Iterator&& );
      Iterator& operator=( Iterator&& );
    private:
      friend class TextData;
      Iterator(const char *);
      struct is_end_t{};
      Iterator(const char *, is_end_t );
      void setup();
      std::string m_buf;
      const char* m_data;
      const char* m_nextData;
    };

  private:
    RawStrData m_data;
    Optional<std::string> m_optOnDisk;
    DataSourceName m_dsn;
    std::string m_dt;
    TextDataUID m_uid;
  public:
    //For factory infrastructure:
    struct internal_with_unset_textdatauid_t {};
    TextData( internal_with_unset_textdatauid_t, RawStrData, DataType,
              Optional<DataSourceName> = NullOpt,
              Optional<LastKnownOnDiskAbsPath> = NullOpt );
    static TextData internal_consumeAndSetNewUID( TextData&& td_with_no_uid );

  };

  std::ostream& operator<< ( std::ostream& , const TextData& );
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  inline TextDataUID::TextDataUID( UniqueIDValue v ) : m_value(v) {}
  inline constexpr TextDataUID::TextDataUID() : m_value(UniqueIDValue{0}) {}
  inline constexpr bool TextDataUID::isUnset() const noexcept { return m_value==UniqueIDValue{0}; }
  inline constexpr bool TextDataUID::operator<(const TextDataUID&o) const noexcept { return m_value<o.m_value; }
  inline constexpr bool TextDataUID::operator==(const TextDataUID&o) const noexcept { return m_value==o.m_value; }
  inline constexpr bool TextDataUID::operator!=(const TextDataUID&o) const noexcept { return m_value!=o.m_value; }
  inline uint64_t TextDataUID::value() const noexcept { return m_value.value; }
  inline void TextDataUID::set( TextDataUID o ) { m_value = o.m_value; }
  inline TextDataUID TextDataUID::createNewUniqueIDValue() { return TextDataUID{ UniqueID().getUniqueID() }; }

  inline RawStrData::RawStrData( static_data_ptr_t, const char * ptr ) noexcept
    : m_b(ptr), m_e( std::next( ptr, std::strlen(ptr) ) )
  {
  }

  inline RawStrData::RawStrData( std::string&& ss, const char * srcdescr )
    : RawStrData( makeSO<std::string>(std::move(ss)), srcdescr )
  {
  }

  inline RawStrData::RawStrData( std::string&& ss, const DataSourceName&  dsn)
    : RawStrData( makeSO<std::string>(std::move(ss)), dsn.str().c_str() )
  {
  }

  inline RawStrData::RawStrData( shared_obj<std::string> ss, const DataSourceName&  dsn)
    : RawStrData( std::move(ss), dsn.str().c_str() )
  {
  }

  inline bool RawStrData::hasSameContent( const std::string& ss ) const
  {
    return hasSameContent( ss.c_str(), std::next(ss.c_str(),ss.size()) );
  }

  inline bool RawStrData::hasSameContent( const RawStrData& o ) const
  {
    return ( ( m_b == o.m_b && m_e == o.m_e )
             ? true//trivially referring to exact same data
             : hasSameContent( o.m_b, o.m_e ) );//more careful analysis needed
  }

  inline const std::string& TextData::dataType() const noexcept { return m_dt; }
  inline const Optional<std::string>& TextData::getLastKnownOnDiskLocation() const noexcept { return m_optOnDisk; }
  inline const RawStrData& TextData::rawData() const noexcept { return m_data; }
  inline std::string TextData::rawDataCopy() const
  {
    return std::string( m_data.begin(),
                        std::distance(m_data.begin(),m_data.end()) );
  }


  inline ncconstexpr17 const TextDataUID& TextData::dataUID() const noexcept { return m_uid; }

  inline TextData::Iterator::Iterator(const char * data, is_end_t)
    : m_data(data), m_nextData(data)
  {
    nc_assert(*m_data=='\0');
  }

  inline TextData::Iterator::Iterator(const char * data)
    : m_data(data)
  {
    setup();
  }

  inline TextData::Iterator TextData::begin() const
  {
    return Iterator( m_data.begin() );
  }

  inline TextData::Iterator TextData::end() const
  {
    return Iterator( m_data.end(), Iterator::is_end_t() );
  }

  inline TextData::Iterator& TextData::Iterator::operator++()
  {
    m_data = m_nextData;
    setup();
    return *this;
  }

  inline ncconstexpr17 const TextData::Iterator::value_type* TextData::Iterator::operator->() const noexcept { return &m_buf; }
  inline ncconstexpr17 const TextData::Iterator::value_type& TextData::Iterator::operator*() const noexcept { return m_buf; }
  inline ncconstexpr17 bool TextData::Iterator::operator==(const Iterator& o) const noexcept { return m_data == o.m_data; }
  inline ncconstexpr17 bool TextData::Iterator::operator!=(const Iterator& o) const noexcept { return m_data != o.m_data; }
  inline ncconstexpr17 bool TextData::Iterator::operator<(const Iterator& o) const noexcept { return m_data < o.m_data; }

  inline TextData::Iterator::Iterator( const Iterator& o )
  {
    *this = o;
  }

  inline TextData::Iterator& TextData::Iterator::operator=( const Iterator& o )
  {
    m_buf.clear();
    m_buf += o.m_buf;//possibly preserve our existing allocation
    m_data = o.m_data;
    m_nextData = o.m_nextData;
    return *this;
  }

  inline TextData::Iterator::Iterator( Iterator&& o )
  {
    *this = o;
  }

  inline TextData::Iterator& TextData::Iterator::operator=( Iterator&& o )
  {
    return *this = o;
  }

  inline std::ostream& operator<< ( std::ostream& os, const TextData& td)
  {
    //NB: Not just "os << m_data.begin()" since that would pass through \r\n
    //newlines:
    for ( auto& line : td )
      os << line << '\n';
    return os;
  }

  inline bool TextData::hasIdenticalMetaData( const TextData& o ) const noexcept
  {
    return m_dsn.str() == o.m_dsn.str() && m_dt == o.m_dt && m_optOnDisk == o.m_optOnDisk;
  }

  inline TextData::TextData( RawStrData rsd, DataType dt,
                             Optional<DataSourceName> dsn,
                             Optional<LastKnownOnDiskAbsPath> lp )
    : TextData( internal_with_unset_textdatauid_t{}, rsd, dt, std::move(dsn), std::move(lp) )
  {
    m_uid = TextDataUID::createNewUniqueIDValue();
  }

}

#endif
