////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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

#include "NCrystal/NCBkgdExtCurve.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCDefs.hh"

NCrystal::BkgdExtCurve::BkgdExtCurve(const Info* ci, bool thermalise)
  : ScatterXSCurve(ci,"BkgdExtCurve",thermalise),
    m_ci(0)
{
  nc_assert_always(ci);
  if (!ci->providesNonBraggXSects())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks NonBraggXSects needed for cross sections.");
  ci->ref();
  m_ci = ci;
  validate();
}

NCrystal::BkgdExtCurve::~BkgdExtCurve()
{
  m_ci->unref();
}

double NCrystal::BkgdExtCurve::crossSectionNonOriented(double ekin) const
{
  return m_ci->xsectScatNonBragg(ekin2wl(ekin));
}
