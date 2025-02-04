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

#include "NCrystal/internal/utils/NCStrView.hh"

namespace NC = NCrystal;

NC::StrView NC::StrView::trimmed() const
{
  nc_assert(has_value());
  const char * it = m_data;
  const char * itE = it + m_size;
  while ( it!=itE && isWhiteSpace(*it) )
    ++it;
  while ( it < itE && isWhiteSpace(*std::prev(itE)) )
    --itE;
  return StrView( it, std::distance(it,itE) );
}

NC::StrView NC::StrView::ltrimmed() const
{
  nc_assert(has_value());
  const char * it = m_data;
  const char * itE = it + m_size;
  while ( it!=itE && isWhiteSpace(*it) )
    ++it;
  return StrView( it, std::distance(it,itE) );
}

NC::StrView NC::StrView::rtrimmed() const
{
  nc_assert(has_value());
  const char * it = m_data;
  const char * itE = it + m_size;
  while ( it < itE && isWhiteSpace(*std::prev(itE)) )
    --itE;
  return StrView( it, std::distance(it,itE) );
}

NC::StrView::size_type NC::StrView::find_first_of( const char * ss ) const {
  //Unfortunately we can not just call the otherwise efficient std::strpbrk,
  //since that expects a null terminated string.
  size_type res = npos;
  auto it = ss;
  while (*it)
    res = std::min<size_type>(res, find(*it++));
  return res;
}

NC::StrView::size_type NC::StrView::find_first_of( StrView ss ) const {
  //Unfortunately we can not just call the otherwise efficient std::strpbrk,
  //since that expects a null terminated string.
  if ( m_size > 1024 && ss.size() > 1 ) {
    //Avoid searching very large strings ss.size() times, if we are lucky enough
    //to get a hit in the beginning.
    auto i = this->substr(0,128).find_first_of(ss);
    if ( i != StrView::npos )
      return i;
  }
  size_type res = npos;
  for ( char ch : ss )
    res = std::min<size_type>(res, find(ch));
  return res;
}

bool NC::StrView::contains_only( StrView sv ) const
{
  nc_assert(has_value());
  nc_assert(sv.has_value());
  for ( auto& e : *this )
    if ( !sv.contains(e) )
      return false;
  return true;
}

const char * NC::detail::strstr_nonullterm( const char * cstr, std::size_t nc,
                                            const char * pattern, std::size_t np )
{
  nc_assert( cstr && pattern && np );
  char p0 = pattern[0];
  while ( true ) {
    if ( nc < np )
      return nullptr;
    auto it = (const char*)std::memchr( cstr, p0, nc );
    if (!it)
      return nullptr;
    nc -= std::distance( cstr, it );
    cstr = it;
    if ( nc >= np && std::memcmp( cstr, pattern, np ) == 0 )
      return cstr;
    //Full pattern not matched here, must start next search one additional
    //step in.
    nc -= 1;
    ++cstr;
  }
}

NC::StrView NC::WordIterator::next()
{
  const char * it = m_text.data();
  const char * itE = it + m_text.size();
  //skip leading WS:
  while ( it!=itE && m_ws.contains(*it) )
    ++it;
  //find next WS:
  const char * it2 = it;
  while ( it2!=itE && !m_ws.contains(*it2) )
    ++it2;
  StrView next_word(it,std::distance(it,it2));
  m_text = StrView(it2,std::distance(it2,itE));
  return next_word;
}

void NC::streamWrappedText( std::ostream& os, StrView text, const WordWrapCfg& cfg )
{
  if ( !text.isSimpleASCII() )
    NCRYSTAL_THROW(BadInput,"Text passed to streamWrappedText must be pure ASCII");
  if ( !(cfg.colwidth > cfg.prefix.size() + 1) && cfg.overflow_is_error )
    NCRYSTAL_THROW(BadInput,"Too small colwidth for given prefix");
  auto word_iter = WordIterator(text,cfg.whitespace);
  const auto actual_width = cfg.colwidth - cfg.prefix.size();
  bool at_line_start = true;
  std::int64_t nleft = actual_width;
  if ( cfg.initial_offset.has_value() && cfg.initial_offset.value() < cfg.colwidth )
    nleft = cfg.colwidth - cfg.initial_offset.value();

  if ( !cfg.initial_offset.has_value() && !cfg.prefix.empty() )
    os << cfg.prefix;

  auto addNewLine = [&os,&nleft,&cfg,&at_line_start,actual_width](){
    os << '\n' << cfg.prefix;
    nleft = actual_width;
    at_line_start = true;
  };

  while (true) {
    auto nextword = word_iter.next();
    if ( nextword.empty()) {
      if ( cfg.ensure_final_newline && !at_line_start )
        os << '\n';
      return;
    }
    auto needs = ( at_line_start ? nextword.size() : nextword.size() + 1 );//+1 for a space in front
    if ( static_cast<std::int64_t>(needs) > nleft ) {
      if ( cfg.overflow_is_error && nextword.size() > actual_width )
        NCRYSTAL_THROW2(BadInput,"Overflow error - word too long to wrap: \""<<nextword<<"\"");
      if ( at_line_start ) {
        //must put it here and live with the overflow
        os << nextword;
        addNewLine();
        continue;
      } else {
        //go to next line before outputting:
        addNewLine();
        os << nextword;
        nleft -= nextword.size();
        at_line_start = false;
        if ( nleft<=0 ) {
          addNewLine();
          at_line_start = true;
        }

      }
    } else {
      if ( !at_line_start)
        os << ' ';
      os << nextword;
      nleft -= needs;
      at_line_start = false;
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  static_assert( StrView::make("foo").size() == 3,"");
  static_assert( StrView::make("").size() == 0,"");
  static_assert( StrView::make("b").size() == 1,"");
  static_assert( StrView::make("").has_value(),"" );
  static_assert( !NC::StrView().has_value(),"" );
  static_assert(constexpr_strcmp(StrView::make(""),StrView::make(""))==0,"");
  static_assert(constexpr_strcmp(StrView::make("a"),StrView::make("a"))==0,"");
  static_assert(constexpr_strcmp(StrView::make("abCDEF"),StrView::make("abCDEF"))==0,"");
  static_assert(constexpr_strcmp(StrView::make("a"),StrView::make("b"))<0,"");
  static_assert(constexpr_strcmp(StrView::make("b"),StrView::make("a"))>0,"");
  static_assert(constexpr_strcmp(StrView::make("aa"),StrView::make("bb"))<0,"");
  static_assert(constexpr_strcmp(StrView::make("bb"),StrView::make("aa"))>0,"");
  static_assert(constexpr_strcmp(StrView::make("aa"),StrView::make("aaa"))<0,"");
  static_assert(constexpr_strcmp(StrView::make("aaa"),StrView::make("aa"))>0,"");
  static_assert(StrView::make("bla bla").constexprLessThan(StrView::make("bla bla yo")),"");
}
