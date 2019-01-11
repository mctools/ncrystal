////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "NCNeutronSCL.hh"
#include "NCMath.hh"

NCrystal::NeutronSCL* NCrystal::NeutronSCL::m_instance = 0;

const NCrystal::NeutronSCL * NCrystal::NeutronSCL::instance()
{
  static NeutronSCL instance;
  return &instance;
}

NCrystal::NeutronSCL::NeutronSCL()
{
  addData("H", 1, 1.00794, -0.37409, 80.26, 0.3326);
  addData("He", 2, 4.002602, 0.326, 0.0, 0.00747);
  addData("Li", 3, 6.941, -0.19, 0.92, 70.5);
  addData("Be", 4, 9.012182, 0.779, 0.0018, 0.0076);
  addData("B", 5, 10.811, 0.53, 1.7, 767.0);
  addData("C", 6, 12.0107, 0.66484, 0.001, 0.0035);
  addData("N", 7, 14.0067, 0.936, 0.5, 1.9);
  addData("O", 8, 15.9994, 0.5805, 0.0, 0.00019);
  addData("F", 9, 18.9984032, 0.5654, 0.0008, 0.0096);
  addData("Ne", 10, 20.1797, 0.4566, 0.008, 0.039);
  addData("Na", 11, 22.98977, 0.363, 1.62, 0.53);
  addData("Mg", 12, 24.305, 0.5375, 0.08, 0.063);
  addData("Al", 13, 26.981538, 0.3449, 0.0082, 0.231);
  addData("Si", 14, 28.0855, 0.415071, 0.004, 0.171);
  addData("P", 15, 30.973761, 0.513, 0.005, 0.172);
  addData("S", 16, 32.065, 0.2847, 0.007, 0.53);
  addData("Cl", 17, 35.453, 0.95792, 5.3, 33.5);
  addData("Ar", 18, 39.948, 0.1909, 0.225, 0.675);
  addData("K", 19, 39.0983, 0.367, 0.27, 2.1);
  addData("Ca", 20, 40.078, 0.47, 0.05, 0.43);
  addData("Sc", 21, 44.95591, 1.21, 4.5, 27.5);
  addData("Ti", 22, 47.867, -0.337, 2.87, 6.09);
  addData("V", 23, 50.9415, -0.0443, 5.08, 5.08);
  addData("Cr", 24, 51.9961, 0.3635, 1.83, 3.05);
  addData("Mn", 25, 54.938049, -0.375, 0.4, 13.3);
  addData("Fe", 26, 55.845, 0.945, 0.4, 2.56);
  addData("Co", 27, 58.9332, 0.249, 4.8, 37.18);
  addData("Ni", 28, 58.6934, 1.03, 5.2, 4.49);
  addData("Cu", 29, 63.546, 0.7718, 0.55, 3.78);
  addData("Zn", 30, 65.409, 0.568, 0.077, 1.11);
  addData("Ga", 31, 69.723, 0.7288, 0.16, 2.75);
  addData("Ge", 32, 72.64, 0.8185, 0.18, 2.2);
  addData("As", 33, 74.9216, 0.658, 0.06, 4.5);
  addData("Se", 34, 78.96, 0.797, 0.32, 11.7);
  addData("Br", 35, 79.904, 0.679, 0.1, 6.9);
  addData("Kr", 36, 83.798, 0.781, 0.01, 25.0);
  addData("Rb", 37, 85.4678, 0.708, 0.5, 0.38);
  addData("Sr", 38, 87.62, 0.702, 0.06, 1.28);
  addData("Y", 39, 88.90585, 0.775, 0.15, 1.28);
  addData("Zr", 40, 91.224, 0.716, 0.02, 0.185);
  addData("Nb", 41, 92.90638, 0.7054, 0.0024, 1.15);
  addData("Mo", 42, 95.94, 0.6715, 0.04, 2.48);
  addData("Tc", 43, 98, 0.68, 0.5, 20.0);
  addData("Ru", 44, 101.07, 0.702, 0.4, 2.56);
  addData("Rh", 45, 102.9055, 0.59, 0.3, 144.8);
  addData("Pd", 46, 106.42, 0.591, 0.093, 6.9);
  addData("Ag", 47, 107.8682, 0.5922, 0.58, 63.3);
  addData("Cd", 48, 112.411, 0.483, 3.46, 2520.0);
  addData("In", 49, 114.818, 0.4065, 0.54, 193.8);
  addData("Sn", 50, 118.71, 0.6225, 0.022, 0.626);
  addData("Sb", 51, 121.76, 0.557, 0.0, 4.91);
  addData("Te", 52, 127.6, 0.568, 0.09, 4.7);
  addData("I", 53, 126.90447, 0.528, 0.31, 6.15);
  addData("Xe", 54, 131.293, 0.469, 0.0, 23.9);
  addData("Cs", 55, 132.90545, 0.542, 0.21, 29.0);
  addData("Ba", 56, 137.327, 0.507, 0.15, 1.1);
  addData("La", 57, 138.9055, 0.824, 1.13, 8.97);
  addData("Ce", 58, 140.116, 0.484, 0.0, 0.63);
  addData("Pr", 59, 140.90765, 0.458, 0.015, 11.5);
  addData("Nd", 60, 144.24, 0.769, 9.2, 50.5);
  addData("Pm", 61, 145, 1.26, 1.3, 168.4);
  addData("Eu", 63, 151.964, 0.53, 2.5, 4530.0);
  addData("Gd", 64, 157.25, 0.95, 151.0, 49700.0);
  addData("Tb", 65, 158.92534, 0.734, 0.004, 23.4);
  addData("Dy", 66, 162.5, 1.69, 54.4, 994.0);
  addData("Ho", 67, 164.93032, 0.844, 0.36, 64.7);
  addData("Er", 68, 167.259, 0.779, 1.1, 159.0);
  addData("Tm", 69, 168.93421, 0.707, 0.1, 100.0);
  addData("Yb", 70, 173.04, 1.241, 4.0, 34.8);
  addData("Lu", 71, 174.967, 0.721, 0.7, 74.0);
  addData("Hf", 72, 178.49, 0.777, 2.6, 104.1);
  addData("Ta", 73, 180.9479, 0.691, 0.01, 20.6);
  addData("W", 74, 183.84, 0.4755, 1.63, 18.3);
  addData("Re", 75, 186.207, 0.92, 0.9, 89.7);
  addData("Os", 76, 190.23, 1.07, 0.3, 16.0);
  addData("Ir", 77, 192.217, 1.06, 0.0, 425.0);
  addData("Pt", 78, 195.078, 0.96, 0.13, 10.3);
  addData("Au", 79, 196.96655, 0.79, 0.43, 98.65);
  addData("Hg", 80, 200.59, 1.2595, 6.6, 372.3);
  addData("Tl", 81, 204.3833, 0.8776, 0.21, 3.43);
  addData("Pb", 82, 207.2, 0.9401, 0.003, 0.171);
  addData("Bi", 83, 208.98038, 0.8532, 0.0084, 0.0338);
  addData("Ra", 88, 226, 1.0, 0.0, 12.8);
  addData("Th", 90, 232.0381, 1.031, 0.0, 7.37);
  addData("Pa", 91, 231.03588, 0.91, 0.1, 200.6);
  addData("U", 92, 238.02891, 0.8417, 0.005, 7.57);
}

void NCrystal::NeutronSCL::addData(std::string name, unsigned z, double m, double csl, double ixs, double cxs)
{
  nc_assert(m_natural_elements.find(name)==m_natural_elements.end());
  nc_assert(m_z_to_elemname.find(z)==m_z_to_elemname.end());
  m_natural_elements.insert(std::pair<std::string, IsotopeData>(name,IsotopeData(z,m,csl,ixs,cxs)));
  m_z_to_elemname[z] = name;
}

const std::string& NCrystal::NeutronSCL::getAtomName(unsigned z) const
{
  std::map<unsigned,std::string>::const_iterator it = m_z_to_elemname.find(z);
  static std::string emptystr;
  if (it==m_z_to_elemname.end())
    return emptystr;
  return it->second;
}

double NCrystal::NeutronSCL::getNeutronWeightedMass(const std::string& element) const
{
   return getAtomicMass(element)/const_neutron_atomic_mass;
}

double NCrystal::NeutronSCL::getAtomicMass(const std::string& element) const
{
  double mass = 0.;
  DataMap::const_iterator it = m_natural_elements.find(element);
  if (it != m_natural_elements.end())
    mass=it->second.mass;
  else
    NCRYSTAL_THROW2(BadInput,"element \"" << element << "\" is misspelled or not supported by current version of NCrystal");
  return mass;
}

unsigned NCrystal::NeutronSCL::getAtomicNumber(const std::string& element) const
{
  unsigned num = 0.;
  DataMap::const_iterator it = m_natural_elements.find(element);
  if (it != m_natural_elements.end())
    num=it->second.atomic_num;
  else
    NCRYSTAL_THROW2(BadInput,"element \"" << element << "\" is misspelled or not supported by current version of NCrystal");
  return num;
}


double NCrystal::NeutronSCL::getBoundXS(const std::string& element) const
{
  double tmp=getCoherentXS(element);
  tmp += getIncoherentXS(element);
  return tmp;
}


double NCrystal::NeutronSCL::getFreeXS(const std::string& element) const
{
  double bound_xs=getBoundXS(element);
  double awr=getNeutronWeightedMass(element);
  double tmp = (1+awr)/awr;
  return bound_xs/tmp/tmp;
}

double NCrystal::NeutronSCL::getCoherentXS(const std::string& element) const
{
  double tmp=getCoherentSL(element);
  return k4Pi*tmp*tmp;
}

double NCrystal::NeutronSCL::getCaptureXS(const std::string& element) const
{
  double tmp=0;
  DataMap::const_iterator it = m_natural_elements.find(element);
  if (it != m_natural_elements.end())
    tmp+=it->second.cap_xs;
  else
    NCRYSTAL_THROW2(BadInput,"element \"" << element << "\" is misspelled or not supported by current version of NCrystal");
  return tmp;

}

double NCrystal::NeutronSCL::getIncoherentXS(const std::string& element) const
{
  double tmp=0;
  DataMap::const_iterator it = m_natural_elements.find(element);
  if (it != m_natural_elements.end())
    tmp+=it->second.inc_xs;
  else
    NCRYSTAL_THROW2(BadInput,"element \"" << element << "\" is misspelled or not supported by current version of NCrystal");
  return tmp;
}

double NCrystal::NeutronSCL::getIncoherentSL(const std::string& element) const
{
  return std::sqrt( kInv4Pi * getIncoherentXS(element) );
}


double NCrystal::NeutronSCL::getCoherentSL(const std::string& element) const
{
  double tmp=0.;
  DataMap::const_iterator it = m_natural_elements.find(element);
  if (it != m_natural_elements.end())
    tmp=it->second.coh_sl;
  else
    NCRYSTAL_THROW2(BadInput,"element \"" << element << "\" is misspelled or not supported by current version of NCrystal");
  return tmp;


}
