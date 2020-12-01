#ifndef NCrystal_Iter_hh
#define NCrystal_Iter_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include <iterator>
#include <utility>
#include <initializer_list>

namespace NCrystal {

  //A bit of template magic to enable loops with functionality similar to
  //Python's enumerate, accessing both values and index numbers of
  //containers.
  //
  //Examples:
  //
  //for (auto&& e : enumerate(mycontainer))
  //   othercontainer[e.idx] = e.val;
  //for (auto&& e : enumerate(mycontainer))
  //   e.val += 0.1*e.idx;
  //
  //It is also possible to explicitly impose read-only access (for code clarity perhaps):
  //
  //  for (const auto& e : enumerate(mycontainer))
  //     othercontainer[e.idx] = e.val;
  //
  //The form (auto& e: enumerate(mycontainer)) does not always compile
  //(depending on constness), but is anyway not needed given the options above.

  template <typename TIter>
  class enumerate_iter_t : std::iterator<std::forward_iterator_tag, typename std::iterator_traits<TIter>::value_type> {
  public:
    using reference = typename std::iterator_traits<TIter>::reference;
  private:
    std::size_t m_idx = 0;
    TIter m_it;
  public:
    explicit enumerate_iter_t(TIter&& it): m_it{ std::forward<TIter>(it) } {}
    enumerate_iter_t& operator++() { ++m_it; ++m_idx; return *this; }
    bool operator!=(const enumerate_iter_t& o) const { return m_it != o.m_it; }
    struct entry_t { const std::size_t idx; reference val; };
    struct const_entry_t { const std::size_t idx; const reference val; };
    entry_t operator*() { return { m_idx, *m_it }; }
    const_entry_t operator*() const { return { m_idx, *m_it }; }
  };
  template <typename T> struct enumerate_wrapper_t { T& container; };
  template <typename T>
  auto begin(enumerate_wrapper_t<T>& w) -> enumerate_iter_t<decltype(std::begin(w.container))>
  {
    return enumerate_iter_t<decltype(std::begin(w.container))>(std::begin(w.container));
  }
  template <typename T>
  auto end(enumerate_wrapper_t<T>& w) -> enumerate_iter_t<decltype(std::begin(w.container))>
  {
    return enumerate_iter_t<decltype(std::begin(w.container))>(std::end(w.container));
  }
  template <typename T>
  enumerate_wrapper_t<T> enumerate(T&& container) { return {container}; }
  template <typename T>
  enumerate_wrapper_t<std::initializer_list<T>> enumerate(std::initializer_list<T>&& container) { return {container}; }

}


#endif
