#ifndef NCrystal_VDOSLux_hh
#define NCrystal_VDOSLux_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

    //////////////////////////////////////////////////////////////////////////
    // Strongly typed decoded VDOSLux, with easy conversion to/from raw int //
    //////////////////////////////////////////////////////////////////////////

  namespace VDOS {

    class VDOSLux final {
    public:
      explicit VDOSLux( int vdoslux_raw );
      VDOSLux();
      VDOSLux reduceLvl( unsigned ) const;
      int raw() const;
      unsigned lvl() const { return m_v & 127u; }
      bool isLegacy() const { return static_cast<bool>(m_v & 128u); }
      bool operator==(const VDOSLux& o) const { return m_v == o.m_v; }
      bool operator!=(const VDOSLux& o) const { return m_v != o.m_v; }
      bool operator<(const VDOSLux& o) const { return m_v < o.m_v; }
    private:
      unsigned m_v;
    };
  }
}

#endif
