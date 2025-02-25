
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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// As an addition to the usual NCrystal license pasted above, note that THIS  //
// PARTICULAR FILE is also placed into the Public Domain, so that it might    //
// be easily adopted into the code-base of any project wishing to use it:     //
//                                                                            //
// This is free and unencumbered software released into the public domain.    //
//                                                                            //
// Anyone is free to copy, modify, publish, use, compile, sell, or            //
// distribute this software, either in source code form or as a compiled      //
// binary, for any purpose, commercial or non-commercial, and by any          //
// means.                                                                     //
//                                                                            //
// In jurisdictions that recognize copyright laws, the author or authors      //
// of this software dedicate any and all copyright interest in the            //
// software to the public domain. We make this dedication for the benefit     //
// of the public at large and to the detriment of our heirs and               //
// successors. We intend this dedication to be an overt act of                //
// relinquishment in perpetuity of all present and future rights to this      //
// software under copyright law.                                              //
//                                                                            //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,            //
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF         //
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.     //
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR          //
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,      //
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR      //
// OTHER DEALINGS IN THE SOFTWARE.                                            //
//                                                                            //
// For more information, please refer to <https://unlicense.org/>             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ncrystal_load
#define ncrystal_load

#include <functional>
#include <memory>

namespace NCrystalVirtualAPI {

  //WARNING: Do NOT make ANY changes in the VirtAPI_Type1_v1 class, it is
  //required to stay exactly constant over time and compatible with the same
  //definition used to compile the NCrystal library! But changes to white space,
  //comments, and formatting is of course allowed.

  class VirtAPI_Type1_v1 {
  public:

    // This Virtual API was introduced in NCrystal 4.1.0.
    //
    // Notes about this Virtual API:
    //
    // * This methods available in this interface was chosen to be "just enough"
    //   for what is needed for usage in OpenMC. Since OpenMC extracts material
    //   composition and density via the Python API, and handles absorption
    //   elsewhere, such information was left out.
    // * Units in this interface are barn/atom for cross sections, and eV for
    //   neutron energy.
    // * The "neutron" parameter below is a pointer to an array of length 4,
    //   with values (ekin,ux,uy,uz) where (ux,uy,uz) is the direction of the
    //   neutron. When sampling, this array is modified directly.
    // * Due to concerns about multi-thread safety when used in OpenMC, this API
    //   does not allow for client-side caching control (for simplicity and
    //   MT safety).

    class ScatterProcess;
    virtual const ScatterProcess * createScatter( const char * cfgstr ) const = 0;
    virtual const ScatterProcess * cloneScatter( const ScatterProcess * ) const = 0;
    virtual void deallocateScatter( const ScatterProcess * ) const = 0;
    virtual double crossSectionUncached( const ScatterProcess&,
                                         const double* neutron ) const = 0;
    virtual void sampleScatterUncached( const ScatterProcess&,
                                        std::function<double()>& rng,
                                        double* neutron ) const = 0;
    //Plumbing:
    static constexpr unsigned interface_id = 1001;//1000*typenumber+version
    virtual ~VirtAPI_Type1_v1() = default;
    VirtAPI_Type1_v1() = default;
    VirtAPI_Type1_v1( const VirtAPI_Type1_v1& ) = delete;
    VirtAPI_Type1_v1& operator=( const VirtAPI_Type1_v1& ) = delete;
    VirtAPI_Type1_v1( VirtAPI_Type1_v1&& ) = delete;
    VirtAPI_Type1_v1& operator=( VirtAPI_Type1_v1&& ) = delete;
  };
}

using NCrystalAPI = NCrystalVirtualAPI::VirtAPI_Type1_v1;

std::shared_ptr<const NCrystalAPI> loadNCrystalAPI();

class NCrystalScatProc final {
public:
  NCrystalScatProc( const char * cfgstr )
    : m_api( loadNCrystalAPI() ),
      m_p( m_api->createScatter( cfgstr ) )
  {
  }
  ~NCrystalScatProc()
  {
    if ( m_p )
      m_api->deallocateScatter( m_p );
  }

  class NeutronState {
  public:
    double& ekin() noexcept { return data[0]; }
    double* u() noexcept { return &data[1]; }
    const double& ekin() const noexcept { return data[0]; }
    const double* u() const noexcept { return &data[1]; }
    double data[4];
  };

  double crossSection( const NeutronState& n ) const
  {
    return m_api->crossSectionUncached( *m_p, n.data );
  }
  void scatter( std::function<double()>& rng, NeutronState& n ) const
  {
    m_api->sampleScatterUncached( *m_p, rng, n.data );
  }

  NCrystalScatProc( const NCrystalScatProc& ) = delete;
  NCrystalScatProc& operator=( const NCrystalScatProc& ) = delete;
  NCrystalScatProc( NCrystalScatProc&& o  )
    : m_api(std::move(o.m_api)), m_p(nullptr)
  {
    std::swap( m_p, o.m_p );
  }
  NCrystalScatProc& operator=( NCrystalScatProc&& o )
  {
    if ( m_p ) {
      m_api->deallocateScatter( m_p );
      m_p = nullptr;
    }
    std::swap( m_api, o.m_api );
    std::swap( m_p, o.m_p );
    return *this;
  }

private:
  std::shared_ptr<const NCrystalAPI> m_api;
  const NCrystalAPI::ScatterProcess * m_p;
};

#endif
