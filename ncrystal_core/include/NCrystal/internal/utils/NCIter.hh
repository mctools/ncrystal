#ifndef NCrystal_Iter_hh
#define NCrystal_Iter_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  //A bit of template magic to enable loops with functionality similar to
  //Python's enumerate, accessing both values and index numbers of
  //containers.
  //
  //Examples:
  //
  //for ( auto&& e : enumerate(mycontainer) )
  //   othercontainer[e.idx] = e.val;
  //for ( auto&& e : enumerate(mycontainer) )
  //   e.val += 0.1*e.idx;
  //
  //The by-value form is also OK, but it might still be possible to modify the
  //elements:
  //
  //for ( auto e : enumerate(mycontainer) )
  //   e.val += 0.1*e.idx;
  //
  //It is also possible to explicitly impose read-only access (for code clarity perhaps):
  //
  //  for (const auto& e : enumerate(mycontainer))
  //     othercontainer[e.idx] = e.val;
  //
  // The form (auto& e: enumerate(mycontainer)) does not always compile
  // (depending on constness), but is anyway not needed given the options above.

  namespace detail {
    template <typename TIter>
    struct EnumIter
    {
      std::size_t idx;
      TIter it;
      bool operator != (const EnumIter & other) const { return it != other.it; }
      void operator ++() { ++idx; ++it; }
      struct Entry {
        std::size_t idx;
        decltype( *TIter()) val;
      };
      Entry operator*() const { return Entry{idx, *it}; }
    };

    template <typename T, typename TIter>
    struct EnumIterWrapper
    {
      static_assert(std::is_reference<T>::value,"Inefficient enumerate(..) usage");
      T container;
      EnumIter<TIter> begin() { return EnumIter<TIter>{ 0, std::begin(container) }; }
      EnumIter<TIter> end() { return EnumIter<TIter>{ 0, std::end(container) }; }
    };

  }
  template <typename T,
            typename TIter = decltype(std::begin(std::declval<T>())),
            typename TIterE = decltype(std::end(std::declval<T>()))>
  inline constexpr detail::EnumIterWrapper<T,TIter> enumerate( T && container )
  {
    static_assert(std::is_same<TIter,TIterE>::value,"");
    return detail::EnumIterWrapper<T,TIter>{ std::forward<T>(container) };
  }
}


#endif
