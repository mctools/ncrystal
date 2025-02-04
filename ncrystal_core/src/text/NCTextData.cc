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
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include <sstream>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    const char * findNextNR0(const char * data) {

      //Search forward until finding the next \n, \r, \r\n or \0. Any \r\n
      //occurance are treated as a single newline separator so will return the
      //position of the \n in \r\n. Any \r not followed by \n will result in
      //BadInput as documented.

      //The search take advantage of the fact that in most input-files, almost
      //all non-newline characters have one of the 4 upper bits set, while the
      //newline chars \r=LF=0x0D=0b00001101 and \n=CR=0x0A=0b00001010 do
      //not. Thus, any character with one of the upper four bits set is neither
      //a newline or null character. And when none of the four upper bits are
      //set, it is almost always because we are dealing with a newline or null
      //character. In fact, the lowest-value non-newline non-null character
      //found in most texts is a space with decimal value 32, 0b00100000. The
      //most notable exception is that some files might contain tab (\t)
      //characters, with value 9=0b00001001.

      //NB: We express the search pattern like '\xF0' to avoid char signedness
      //platform issues.
      while ( '\xF0' & *data )
        ++data;
      if ( *data == '\0' || *data == '\n' )
        return data;
      if ( *data == '\r' ) {
        if ( *std::next(data) == '\n' )
          return std::next(data); //dos new-line
        else
          NCRYSTAL_THROW(BadInput,"Data with ancient pre-OSX Mac line-endings is explicitly not allowed!");
      }
      //For the rare case of files with non-newline, non-null control chars, the
      //search continues here:
      return findNextNR0(++data);
    }
  }
}

void NC::TextData::Iterator::setup() {
  m_buf.clear();
  if ( *m_data == '\0' ) {
    //already at end:
    m_nextData = m_data;
    return;
  }
  const char * e = m_nextData = findNextNR0(m_data);
  //Move e so it points at the last char of the data we want to present to
  //the user (we don't show \r, \n or \r\n chars):
  if ( *m_nextData == '\n' && *(std::prev(m_nextData)) == '\r' )
    --e;
  //Copy over the line data and ensure null-termination;
  auto newn = static_cast<std::string::size_type>(std::distance(m_data,e)+1);
  if (  newn > m_buf.capacity() )
    m_buf.reserve( std::max<std::string::size_type>( newn, static_cast<std::string::size_type>(256) ) );
  m_buf.append(m_data, newn );
  m_buf.back() = '\0';
  m_buf.pop_back();
  //Make sure m_nextData points at the beginning of the next line (or the final
  //\0 char):
  if ( *m_nextData != '\0' )
    ++m_nextData;
}

void NC::TextData::verifyOnDiskFileUnchanged() const
{
  if ( !m_optOnDisk.has_value() )
    NCRYSTAL_THROW(BadInput,"TextData::verifyOnDiskFileUnchanged called for object without on-disk location");

  const std::string& path = m_optOnDisk.value();
  Optional<std::string> current_content = readEntireFileToString( path );
  if ( !current_content.has_value() )
    NCRYSTAL_THROW2(BadInput,"File disappeared or became unreadable: "<<path);
  if ( !m_data.hasSameContent( current_content.value() ) )
    NCRYSTAL_THROW2(BadInput,"File unexpectedly changed content while being used: "<<path);
}

uint64_t NC::RawStrData::checkSumFromRawStringData(const char* c_begin, const char* c_end)
{
  //Very simply checksum alg, using randomly generated initial value (which will
  //be the checksum of an empty range):
  uint_fast64_t checksum(static_cast<uint64_t>(0x2254a62a1af0a16bull));
  unsigned shift = 0;
  for ( ;c_begin!=c_end; ++c_begin ) {
    checksum += ( static_cast<unsigned char>(*c_begin) << shift);
    shift = ( shift + 8) % 64;
  }
  return static_cast<uint64_t>(checksum);
}

uint64_t NC::RawStrData::calcCheckSum() const
{
  return checkSumFromRawStringData(m_b,m_e);
}

bool NC::RawStrData::hasSameContent( const char* dataBegin, const char* dataEnd ) const
{
  auto it = m_b;
  auto itE = m_e;
  auto it2 = dataBegin;
  auto it2E = dataEnd;
  const auto n = (std::size_t)(std::distance(it,itE));
  const auto n2 = (std::size_t)(std::distance(it2,it2E));
  if ( n != n2 )
    return false;
  if ( n == 0 )
    return true;
  if ( it == it2 )
    return true;//pointing at the same data
  return std::memcmp( it, it2, n ) == 0;
}

NC::TextData NC::TextData::internal_consumeAndSetNewUID( TextData&& td_with_no_uid )
{
  nc_assert_always(td_with_no_uid.dataUID().isUnset());
  NC::TextData res = std::move(td_with_no_uid);
  res.m_uid = TextDataUID::createNewUniqueIDValue();
  return res;
}

NC::TextData::TextData( internal_with_unset_textdatauid_t, RawStrData data,
                        DataType dt,
                        Optional<DataSourceName> dsn,
                        Optional<LastKnownOnDiskAbsPath> ondisk )
  : m_data(std::move(data)),
    m_dt(dt.value),
    m_uid()
{
  nc_assert(m_uid.isUnset());
  if ( m_dt.empty() || !isAlphaNumeric(m_dt)  )
    NCRYSTAL_THROW(BadInput,"Error: Data type must be alpha numeric and non-empty.");
  if ( ondisk.has_value() )
    m_optOnDisk = ondisk.value().value;
  if ( dsn.has_value() && !dsn.value().str().empty() ) {
    m_dsn = std::move(dsn.value());
  } else {
    std::ostringstream s;
    s << "(anonymous TextData, "<<std::distance( m_data.begin(), m_data.end() )<<"bytes"
      <<", type="<<m_dt;
    s<<")";
    //Explicitly do not add m_optOnDisk, it is only for from-file-reading.
    m_dsn = s.str();
  }
}

NC::RawStrData::RawStrData( shared_obj<std::string> d, const char * srcdescr )
  : m_s( std::move(d) )
{
  //Verify input data does not contain unexpected null chars:
  const std::string& s = *m_s;
  auto n = s.size();
  m_b = s.c_str();
  m_e = m_b + n;
  if ( std::strlen( m_b ) != (std::size_t)n ) {
    //Some extraneous null character must have spoiled it!
    NCRYSTAL_THROW2(BadInput,"Invalid text data"
                    <<(srcdescr?" in ":"")
                    <<(srcdescr?srcdescr:"")
                    <<": Data is not in UTF-8 or ASCII format.");
  }
}
