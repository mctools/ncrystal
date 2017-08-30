////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCPhononDebye.hh"
#include "NCPhononData.hh"
#include "NCrystal/NCException.hh"
#include "NCMath.hh"
#include "NCNeutronSCL.hh"
#include "NCrystal/NCCalcBase.hh"

#include <algorithm>


NCrystal::PhononDebye::PhononDebye( double debye_energy, double kt,
                                    const std::string& ele_name,
                                    unsigned max_phonon_order)
  : m_debye(debye_energy),
    m_kt(kt),m_ele_name(ele_name),
    m_max_phononnum(max_phonon_order),
    m_max_wl2ekin(4.)
{
  if(max_phonon_order>1000)
    NCRYSTAL_THROW(BadInput,"PhononDebye model phonon expansion above 1000 orders is numerically unstable");

  m_dt = 0;
  m_cal = 0;

  unsigned num=100;
  double diff = m_debye/num;
  std::vector<double> xvec(10);
  std::vector<double> w(10);
  m_gamma0=0;  // ref. 1, eq II.3
  m_gamma0_bar =0; // ref. 1, eq II.52'
  m_delta0_bar=0; // ref. 1, eq II.55

  for(unsigned i=0;i<num;i++)
  {
    gauleg_10_ord(i*diff, (i+1)*diff, xvec, w);
    for(unsigned j=0;j<10;j++)
    {
      m_gamma0 += g0knl(xvec[j])*w[j];
      m_gamma0_bar += g0barknl(xvec[j])*w[j];
      m_delta0_bar += d0bar_power_knl(xvec[j])*w[j];
    }
  }

  m_delta0_bar = sqrt(m_delta0_bar/m_gamma0_bar);
  double A = NeutronSCL::instance()->getNeutronWeightedMass(m_ele_name);
  m_msd = const_hhm*m_gamma0/2./A;

  if(m_max_phononnum==0)
  {

    double x = (m_gamma0/A + 1/2./m_kt) * m_delta0_bar;
    double factor = sqrt(1+0.25*x*x)-0.5*x; // III.22
    double z = m_gamma0_bar/A * m_delta0_bar * factor * exp(0.5*factor*factor);  // III.22
    //printf("big three %f %f %f, x %f\n",m_gamma0,m_gamma0_bar,m_delta0_bar,x);

    para p;
    p.x=x;
    p.z=z;
    p.factor=factor;
    p.delta0_bar=m_delta0_bar;
    p.kt = m_kt;

    double sum_n1 = 0;
    for(unsigned i=1;i<100;i++)
    {
      double contribution = sigma_1(p,i) ;
      if(sum_n1!=0 && contribution/sum_n1<1e-5)
      {
        break;
      }
      else
      {
        sum_n1 +=contribution;
        //printf("%u, %e\n",i ,contribution);
      }
    }
    //printf ("sum for sig1 is %f\n",sum_n1);


    double sum_n2 = 0;
    for(unsigned i=1;i<100;i++)
    {
      double contribution = sigma_2(p,i) ;
      if(sum_n2!=0 && contribution/sum_n2<1e-5)
      {
        break;
      }
      else
      {
        sum_n2 +=contribution;
        //printf("%u, %e\n",i ,contribution);
      }
    }
    //printf ("sum for sig2 is %f\n",sum_n2);


    double tot_xs = 0;
    unsigned breaking=100;

    for(unsigned i=1;i<100;i++)
    {
      double contribution = pow(sum_n1,i) + pow(sum_n2,i) ;
      if(tot_xs!=0 && contribution/tot_xs<1e-5)
      {
        breaking=i;
        break;
      }
      else
      {
        tot_xs +=contribution;
        //printf("%u, %e\n",i ,contribution);
      }
    }
    m_max_phononnum =breaking*4;
    if(m_max_phononnum>100)
      m_max_phononnum=100;
  }

  //printf (" ord %u\n", m_max_phononnum);

}


void NCrystal::PhononDebye::doit(const std::vector<double> &ekin_vec, std::vector<double> &xs_vec, unsigned alpha_grid_size,  unsigned beta_sym_grid_size)
{
  //1. phonon expansion
  //1.1 make asym g1 array
  unsigned psize=51;
  m_dt = m_debye/(psize-1);


  std::vector<double> phonon_energy=linspace(0,m_debye,psize);
  std::vector<double> g1_arr(phonon_energy.size());
  g1_arr[0]=3/pow(m_debye,3)*(2*m_kt)/(2.*m_gamma0);
  for(unsigned i=1;i<g1_arr.size();i++)
  {
    g1_arr[i]=phonon_energy[i]*3/pow(m_debye,3)/sinh(phonon_energy[i]/(2*m_kt))/(2.*m_gamma0);
  }
  std::vector<double> g1_flip;
  flip(g1_arr,g1_flip);
  //concatenate(g1_flip,g1_arr,1);
  g1_flip.insert(g1_flip.end(), g1_arr.begin()+1, g1_arr.end());
  g1_arr.swap(g1_flip);
  double g1_start_t = -phonon_energy.back();
  for(unsigned i=0;i<g1_arr.size();i++)
  {
    g1_arr[i] *= exp(-(g1_start_t+i*m_dt)/2/m_kt);
  }
  // m_phonon_spec_arr.clear();
  // m_phonon_spec_arr[0]=std::make_pair(g1_arr,g1_start_t); 

  //1.2 phonon expansion
  m_cal = new PhononCalculator(g1_arr,g1_start_t,m_dt,m_max_phononnum,m_kt);
  //2. fill a 2d table
  //make sym kernel
  //upper limits suggested in F.G. Bischoff and M.L. Yeater, Nucl. Sci. Eng. 48, 266-280, 1972
  double lower_beta = 1e-3/m_kt;
  double lower_alpha =  4*lower_beta;
  double thoeredical_beta_uplim = m_max_wl2ekin/m_kt;
  //limited by the representation of the factor exp(beta/2) in double
  //to be safe, multiply by 1.9 instead of 2.0
  double numerical_uplim = log(std::numeric_limits<double>::max())*1.9;
  double upper_beta = ncmin(thoeredical_beta_uplim,numerical_uplim);
  double upper_alpha = 4*upper_beta;
  //hard coded size, should be convergent because of the well behaved input Debye phonon spectrum

  m_alpha = logspace(log10(lower_alpha),log10(upper_alpha),alpha_grid_size);
  m_beta = logspace(log10(lower_beta),log10(upper_beta),beta_sym_grid_size-1);
  m_beta.insert(m_beta.begin(),0.);
  calcSymSab(m_alpha,m_beta,m_sab);

  //kernel expansion
  std::vector<double> neg_beta;
  flip(m_beta,neg_beta,true);
  //concatenate(neg_beta,m_beta,1);
  neg_beta.insert(neg_beta.end(), m_beta.begin()+1, m_beta.end());
  std::swap(neg_beta,m_beta);
  std::vector<double> temp_sab(alpha_grid_size*(beta_sym_grid_size-1));

  for(unsigned i=0;i<beta_sym_grid_size-1;i++)  {
    std::copy(m_sab.begin() + (beta_sym_grid_size-1-i)*alpha_grid_size  ,
        m_sab.begin()+(beta_sym_grid_size-i)*alpha_grid_size  ,
        temp_sab.begin()+i*alpha_grid_size);
  }
  //concatenate(temp_sab,m_sab,0);
  temp_sab.insert(temp_sab.end(), m_sab.begin(), m_sab.end());
  std::swap(temp_sab,m_sab);


  if(m_sab.size()!=m_alpha.size()*m_beta.size())
    NCRYSTAL_THROW(BadInput,"ScatteringKernel::init just wrong ");

  //3. scale the m_sab to make it asym

  double sigma_b = NCrystal::NeutronSCL::instance()->getBoundXS(m_ele_name);
  double * sab_begin = &m_sab[0];
  size_t alpha_size = m_alpha.size();


  for(unsigned i=0;i<m_beta.size();i++)
  {
    nc_assert( alpha_size*(i+1) <= m_sab.size() );
    double * rowIt = sab_begin + alpha_size*i;
    double * rowItEnd = rowIt + alpha_size;
    nc_assert( rowItEnd>sab_begin && rowItEnd <= sab_begin+m_sab.size() );
    double c = exp(-0.5*m_beta.at(i)) * sigma_b * 0.25 * m_kt;
    for (;rowIt!=rowItEnd;++rowIt)
      *rowIt *= c;
  }

  //check for error
  for(std::vector<double>::const_iterator it = m_sab.begin();it!=m_sab.end();++it)
  {
    if(ncisnan(*it)||ncisinf(*it))
      NCRYSTAL_THROW(CalcError,"Calculated scattering kernel contains nan or inf. It is likely caused by that the factor exp(neutron_energy/Kt/2) or its reciprocal is beyond the valid range of a double.");
  }

  //4. 2D integration of the table for total xs
  unsigned en_mesh = ekin_vec.size();

  //first dimension integration
  xs_vec.resize(en_mesh,0.);
  std::vector<double> firstInt;
  firstInt.resize(en_mesh*m_beta.size(),0.);

  for(unsigned i=0;i<ekin_vec.size();i++)
  {
    getSecondarySpectrum( ekin_vec.at(i),&firstInt.at(m_beta.size()*i) );
  }

  //second dimension integration
  std::vector<double>::const_iterator secbegin = firstInt.begin();
  for(unsigned i=0;i<ekin_vec.size();i++)
  {
    std::vector<double> scatted_enp (secbegin + i*m_beta.size() , secbegin + (i+1) *m_beta.size() );
    std::vector<double> beta(m_beta);
    double beta_lower_limit = -ekin_vec.at(i)/m_kt;
    std::vector<double>::iterator second_bin = std::upper_bound(beta.begin(),beta.end(), beta_lower_limit);
    if (second_bin!=beta.begin())
    {
      unsigned distance = second_bin - beta.begin()-1;
      beta.at(distance)=beta_lower_limit;
    }
    xs_vec.at(i)=NCrystal::trapz(scatted_enp, beta )/ekin_vec.at(i);
  }

  //add incoherent zero order scattering cross section (zero order is a special
  //case, for the other orders we include both incoherent and coherent parts
  //above (this is the incoherent approximation), while the zero order coherent
  //contribution is provided by Bragg diffraction, treated elsewhere:

  double inco_sca_xs = NCrystal::NeutronSCL::instance()->getIncoherentXS(m_ele_name);
  double dwB = m_msd*8*M_PI*M_PI;

  for(unsigned i = 0;i<xs_vec.size();i++)
  {
    double wavelength_Aa = ekin2wl(ekin_vec.at(i));
    double temp=dwB/(wavelength_Aa*wavelength_Aa);
    xs_vec[i] += inco_sca_xs/2/temp*(1-exp(-2*temp));
  }

}

NCrystal::PhononDebye::~PhononDebye()
{
  delete m_cal;
}

double NCrystal::PhononDebye::getMSD() const
{
  return m_msd;
}
double NCrystal::PhononDebye::integrateAlphaInterval(double a1,double s1, double a2 , double s2 ) const
{
  if (s1==s2)
    return (a2-a1)*s1;
  if (s1==0. || s2==0.)
    return (a2-a1)*(s2+s1)/2.;
  return (a2-a1)/log(s2/s1)*(s2-s1);
}


void NCrystal::PhononDebye::getSecondarySpectrum(double kin, double* spec)
{
  if(kin<0.)
    NCRYSTAL_THROW(CalcError,"Incident neutron energy is negative.");

  //spec.resize(m_beta.size(), 0.);
  double lbeta=-kin/m_kt;
  double alpha_lower =0., alpha_upper=0.;

  for(unsigned b_idx=0;b_idx<m_beta.size();b_idx++)  {
    if (m_beta[b_idx] < lbeta) //lower than beta lower limit
      continue;

    getAlphaLimts(kin, m_beta[b_idx], alpha_lower, alpha_upper);
    if(alpha_lower >= alpha_upper )  //a singularity when kin==-kTB (scattered energy is zero). There is no scattered angle.
      continue;
    bool firstAlphaBin=true;
    const size_t alpha_size = m_alpha.size();
    double frontS = -1;
    double backS = -1;
    for(unsigned a_idx=0;a_idx<alpha_size;a_idx++)  {
      if(m_alpha[a_idx] <= alpha_lower )
        continue;
      if(m_alpha[a_idx] > alpha_lower )  {
        if(firstAlphaBin)  { // first alpha bin
          if(m_alpha[a_idx] > alpha_upper ) { //narrower than one interval
            frontS = (frontS==-1?getS(b_idx, alpha_lower):frontS);
            backS = (backS==-1?getS(b_idx, alpha_upper):frontS);
            spec[b_idx] = integrateAlphaInterval(alpha_lower, frontS,
                                                 alpha_upper, backS );
            break;
          }
          else  {
            frontS = (frontS==-1?getS(b_idx, alpha_lower):frontS);
            spec[b_idx] += integrateAlphaInterval(alpha_lower, frontS,
                                                  m_alpha[a_idx], m_sab[b_idx*alpha_size+a_idx]);
          }
          firstAlphaBin=false;
        }
        else if(m_alpha[a_idx] <alpha_upper) { //middle bins
          spec[b_idx]+= integrateAlphaInterval( m_alpha[a_idx-1], m_sab[b_idx*alpha_size+a_idx-1],
                                                m_alpha[a_idx],   m_sab[b_idx*alpha_size+a_idx]);
        }
        else if(m_alpha[a_idx] > alpha_upper) { // the last bin
          backS = (backS==-1?getS(b_idx, alpha_upper):frontS);

          spec[b_idx] += integrateAlphaInterval( m_alpha[a_idx-1], m_sab[b_idx*alpha_size+a_idx-1],
                                                 alpha_upper,   backS);
          break;

        }
      }
    }
  }


}

void NCrystal::PhononDebye::getAlphaLimts(double kin, double beta, double &lower, double& upper)
{
  double kTB=m_kt*beta;
  double dif = 2*sqrt(kin*(kin + kTB ));
  lower=(2*kin + kTB - dif)/m_kt;
  upper=(2*kin + kTB + dif)/m_kt;
}

double  NCrystal::PhononDebye::interpolate(double a, double fa, double b, double fb, double x) const
{
  nc_assert(a!=b);
  if( ! (fa*fb) ) //when loglin interpolation is invalid, use linlin in the place
    return  fa+(fb-fa)*(x-a)/(b-a);
  //Note from TK: The formula below was: "return exp(log(fa)+(log(fb/fa))*(x-a)/(b-a));",
  //but it was rewritten like this for efficiency:
  return fa*pow(fb/fa,(x-a)/(b-a));
}

double NCrystal::PhononDebye::getS(unsigned beta_index, double alpha)
{
  std::vector<double>::const_iterator it_up = std::upper_bound(m_alpha.begin(), m_alpha.end(), alpha);
  size_t dis = it_up - m_alpha.begin();
  const size_t bm = beta_index*m_alpha.size();
  if (dis) {
    nc_assert(it_up!=m_alpha.end() && bm+dis<m_sab.size());
    double * sab_bmdis_m1 = &m_sab[bm+(dis-1)];
    double * sab_bmdis = sab_bmdis_m1; ++sab_bmdis;
    return interpolate(*(it_up-1),  *sab_bmdis_m1,*it_up, *sab_bmdis, alpha);
  }
  //extrapolate
  nc_assert(m_alpha.size()>=2);
  return interpolate(m_alpha.front(),m_sab[bm],
                     m_alpha[1],m_sab[bm+1], alpha);
}


void NCrystal::PhononDebye::calcSymSab(const std::vector<double> &alpha,const std::vector<double> &beta, std::vector<double> & sab) const
{
  std::vector<double>::const_iterator itE = beta.end();
  for(std::vector<double>::const_iterator it = beta.begin();it!=itE;++it) {
    if((*it)<0)
      NCRYSTAL_THROW(CalcError,"member function calcSymSab only calculates neutron energy lost, so input beta"
                     " should always be positive. The cross section for energy gain is automatically calculated"
                     " in the theory of details balancing.");
  }

  sab.resize(alpha.size()*beta.size(),0.);

  size_t n_alpha = alpha.size();
  const double c1(2*m_kt/const_hhm);

  std::vector<double> ceo_arr;

  for(unsigned a_idx=0;a_idx<n_alpha;++a_idx)
  {

    double ksqr = alpha[a_idx]*c1;
    double doubleW = m_msd*ksqr;
    double exp_minus_doubleW = std::exp(-doubleW);
    ceo_arr.clear();
    ceo_arr.resize(m_max_phononnum+1,1.);

    for(unsigned i=1;i<m_max_phononnum+1;++i)
      ceo_arr[i] = doubleW*ceo_arr[i-1]/i;

    for(unsigned b_idx=0;b_idx<beta.size();++b_idx)
    {
      const double beta_b_idx = beta[b_idx];
      double eps= -beta_b_idx*m_kt; //Negative beta contains less noise

#if 1
      //Straight-forward sum
      double S=0.;
      for(unsigned n=1;n<=m_max_phononnum;++n)
      {
        double weight = m_cal->interpolate(n-1,eps);
        S += ceo_arr[n] * weight;
      }
#else
      //Sum using Neumaier's method, for increased numerical stability:
      //(disabled, the numerical issues are likely from NCFastConvolve::fftd, not here)
      double S = (m_max_phononnum ? ceo_arr[1] * m_cal->interpolate(0,eps) : 0.0);
      double correction = 0.0;
      for(unsigned n=2;n<=m_max_phononnum;++n)
      {
        double val = ceo_arr[n] * m_cal->interpolate(n-1,eps);
        double t = S + val;
        if (ncabs(S)>=ncabs(val)) {
          correction += (S - t) + val;
        } else {
          correction += (val - t) + S;
        }
        S = t;
      }
      S += correction;
#endif
      sab[b_idx*n_alpha+a_idx] = exp_minus_doubleW * exp(-0.5*beta_b_idx) * S * m_kt;
    }
  }
}

double NCrystal::PhononDebye::g0knl(double omega) const
{
  double d3=m_debye*m_debye*m_debye;
  if(omega==0.)
    return d3*(6*m_kt);
  return 3./d3*omega/tanh(omega/(2*m_kt));

}

double NCrystal::PhononDebye::g0barknl(double omega) const
{
  double d3=m_debye*m_debye*m_debye;
  if(omega==0.)
    return d3*(6*m_kt);
  return 3./d3*omega/sinh(omega/(2*m_kt));
}

double NCrystal::PhononDebye::d0bar_power_knl(double omega) const
{
  return g0barknl(omega)*omega*omega;
}

double NCrystal::PhononDebye::sigma_1 (const para& p, double n)
{
  double h1=1./8*p.x*(p.factor+p.x)/(1+0.25*p.x*p.x) - 1./8;
  return sqrt(p.delta0_bar/M_PI) /2 * p.factor /( sqrt(p.factor+.5*p.x)) * pow(p.z,n) * (1.+1./n*h1 );
}


double NCrystal::PhononDebye::sigma_2 (const para& p, double n)
{
  double p1 =  2./3 * p.factor*p.factor*p.factor  + 4/3.*p.factor*p.factor/2/p.kt*p.delta0_bar
      + 2./3/4/p.kt/p.kt*p.delta0_bar*p.delta0_bar*p.factor ;
  return 0.5/sqrt(p.delta0_bar*M_PI) *p.factor /sqrt(p.factor+0.5*p.x) * pow(p.z,n) * (n*p1 );
}
