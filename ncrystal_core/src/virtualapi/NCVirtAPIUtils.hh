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
#include <functional>

namespace NCRYSTAL_NAMESPACE {
  namespace VirtAPIUtils {

    class RNGWrapper : public RNGStream {
      //Intended for the case where a provider function clearly outlives the
      //RNGWrapper (i.e. during a scatter sampling).
    public:
      constexpr RNGWrapper(std::function<double()>* f) : m_f(f) {}

    protected:
      double actualGenerate() override
      {
        return std::max<double>(std::numeric_limits<double>::min(),(*m_f)());
      }
    private:
      std::function<double()> * m_f;
    };

  }
}
