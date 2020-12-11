
 ///////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                 //
// This source and associated header file contain both the NXS Crystallography library (Mirko      //
// Boin, HZB) [nxslib] and the Space Group Info library (Ralf W. Grosse-Kunstleve) [SgInfo], with  //
// modifications by NCrystal developers (clearly marked by comments in the code and reachable by   //
// searching for "NCrystal"). The modifications include a few bug-fixes, but were mainly needed    //
// in order to be able to contain all of both nxslib and SgInfo sources in a single .hh/.cc file   //
// pair, compilable as C++ and protected by a namespace. Both nxslib and SgInfo sources were       //
// lifted from the NXSG4 distribution (http://cern.ch/nxsg4/), since they were both reproduced     //
// there with clear license text, repeated in the following:                                       //
//                                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                 //
//  The nxs-Geant4 integration code (NXSG4) was written by:                                        //
//                                                                                                 //
//    Thomas Kittelmann, European Spallation Source ESS AB <thomas.kittelmann@esss.se>             //
//                                                                                                 //
//  The nxs library was written by:                                                                //
//                                                                                                 //
//    Mirko Boin, Helmholtz Centre Berlin for Materials and Energy <boin@helmholtz-berlin.de>      //
//                                                                                                 //
//  The nxs library includes the Space Group Info (SgInfo) library by Ralf W. Grosse-Kunstleve.    //
//                                                                                                 //
//  The NXSG4 code is not an official product of, nor endorsed by the Geant4 collaboration.        //
//                                                                                                 //
//  Copyright notices for nxs, NXSG4 and SgInfo which must be observed in order to                 //
//  guarantee their free usage follows here:                                                       //
//                                                                                                 //
//    nxs (c) 2010-14 Mirko Boin, Helmholtz Centre Berlin for Materials and Energy                 //
//    NXSG4 (c) 2013-14 Thomas Kittelmann, European Spallation Source ESS AB                       //
//    Permission to use and distribute this software and its documentation for noncommercial       //
//    use and without fee is hereby granted, provided that the two above copyright notices appears //
//    in all copies and that both that copyright notice and this permission notice appear in       //
//    the supporting documentation. It is not allowed to sell this software in any way. This       //
//    software is not in the public domain.                                                        //
//                                                                                                 //
//    Space Group Info (c) 1994-96 Ralf W. Grosse-Kunstleve                                        //
//    Permission to use and distribute this software and its documentation for noncommercial       //
//    use and without fee is hereby granted, provided that the above copyright notice appears      //
//    in all copies and that both that copyright notice and this permission notice appear in       //
//    the supporting documentation. It is not allowed to sell this software in any way. This       //
//    software is not in the public domain.                                                        //
//                                                                                                 //
 ///////////////////////////////////////////////////////////////////////////////////////////////////

#define SGCLIB_C__
#include "NCNXSLib.hh"

#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>

namespace NCrystal {
namespace nxs {

/*!
    \section abstract_sec Abstract
    nxs: a program library for neutron cross section calculations
    Purpose: Calculates neutron scattering and absorption cross sections for
    polycrystalline materials based on the composition of a crystallographic unit
    cell.

    @author Mirko Boin, Helmholtz-Zentrum Berlin f&uuml;r Materialien und Energy GmbH, <boin@helmholtz-berlin.de>
    @version 1.4

*/

/*! \mainpage nxs - a library for neutron cross section calculations

    \section intro_sec Introduction
    The nxs library for computing neutron scattering and absorption cross sections provides
    a number of C structs and functions to calculate wavelength-dependent cross section values for
    polycrystalline/powder-like materials. The definition of a material is represented by the composition
    of a unit cell (NXS_UnitCell). A unit cell is created from the specification of a space group and its
    unit cell parameters. The SgInfo routines from Ralf W. Grosse-Kunstleve is included here for such
    purposes. Monoatomic materials as well as multi-atomic compounds are created by adding NXS_AtomInfo
    atom information/properties. The library also provides a reading and saving routines to compose unit
    cells from nxs parameter files.


    \section usage_sec Usage

    Example:

    The below example shows howto quickly use the library routines to initialise a unit cell and
    calculate some cross sections.
    \code{.c}
      int numAtoms = nxs_readParameterFile( nxsFileName, &uc, &atomInfoList);
      if( numAtoms > 0 )
      {
        int i=0;
        nxs_initUnitCell(&uc);
        for( i=0; i<numAtoms; i++ )
          nxs_addAtomInfo( &uc, atomInfoList[i] );

        nxs_initHKL( &uc );

        double lambda=0.1;
        for( lambda=0.1; lambda<4.0; lambda+=0.1 )
        {
          printf("%f\n",nxs_Absorption(lambda, &uc ) );
        }
      }
    \endcode


    \section parameterfiles_sec Parameter files
    Currently, nxs provides routines to read and save nxs parameter files of one particular kind. It is a
    human-readable, INI-like, format to store the necessary information for the compositon of a crystallographic
    unit cell. An example of NaCl is given below:

    \verbatim
    #
    # This is an nxs parameter file
    #

    # define the unit cell parameters:
    #   space_group                      - the space group number or Hermann or Hall symbol [string]
    #   lattice_a, ...b, ...c            - the lattice spacings a,b,c [angstrom]
    #   lattice_alpha, ...beta, ...gamma - the lattice angles alpha,beta,gamma [degree]
    #   debye_temp                       - the Debye temperature [K]

    space_group=225
    lattice_a=5.64
    lattice_c=4.95
    lattice_alpha=90
    debye_temp=320

    # add atoms to the unit cell:
    # notation is "atom_number = name b_coh sigma_inc sigma_abs_2200 molar_mass x y z"
    #   name           - labels the current atom/isotope  [string]
    #   b_coh          - the coherent scattering length [fm]
    #   sigma_inc      - the incoherent scattering cross section [barns]
    #   sigma_abs_2200 - the absorption cross sect. at 2200 m/s [barns]
    #   molar_mass     - the Molar mass [g/mol]
    #   x y z          - the Wyckoff postion of the atom inside the unit cell
    #
    # e.g.: add_atom = Fe 9.45 0.4 2.56 55.85 0.0 0.0 0.0

    [atoms]
    add_atom=Na 3.63 1.62 0.53 22.99 0.0 0.0 0.0
    add_atom=Cl 9.577 5.3 33.5 35.45 0.5 0.5 0.5
    \endverbatim

    \section copyright_sec Copyright
    nxs - neutron cross sections (c) 2010-2014 Mirko Boin

    The nxs library includes the SgInfo library, whose free usage is granted by the following notice:

    Space Group Info (c) 1994-96 Ralf W. Grosse-Kunstleve
    Permission to use and distribute this software and its documentation for noncommercial
    use and without fee is hereby granted, provided that the above copyright notice appears
    in all copies and that both that copyright notice and this permission notice appear in
    the supporting documentation. It is not allowed to sell this software in any way. This
    software is not in the public domain.


    \section references_sec References

    \subsection mainreference_sec Main references

      \li Boin, M. (2012). J. Appl. Cryst. 45. 603-607. doi: 10.1107/S0021889812016056

    \subsection appreference_sec Applications of nxs

      \li Boin, M.; Hilger, A.; Kardjilov, N.; Zhang, S. Y.; Oliver, E. C.; James, J. A.; Randau, C. & Wimpory, R. C. (2011). J. Appl. Cryst. 44. 1040-1046. doi: 10.1107/S0021889811025970
      \li Strobl, M.; Hilger, A.; Boin, M.; Kardjilov, N.; Wimpory, R.; Clemens, D.; M&uuml;hlbauer, M.; Schillinger, B.; Wilpert, T.; Schulz, C.; Rolfs, K.; Davies, C. M.; O'Dowd, N.; Tiernan, P. & Manke, I. (2011). Nucl. Instrum. Methods Phys. Res., Sect. A 651. 149-155. doi: 10.1016/j.nima.2010.12.121
      \li Boin, M.; Wimpory, R. C.; Hilger, A.; Kardjilov, N.; Zhang, S. Y. & Strobl, M. (2012). J. Phys: Conf. Ser. 340. 012022. doi: 10.1088/1742-6596/340/1/012022
      \li Kittelmann, T.; Stefanescu, I.; Kanaki, K.; Boin, M.; Hall-Wilton, R. & Zeitelhack, K. (2014). J. Phys: Conf. Ser. 513. 022017. doi: 10.1088/1742-6596/513/2/022017
      \li Peetermans, S.; Tamaki, M.; Hartmann, S.; Kaestner, A.; Morgano, M. & Lehmann, E. H. (2014). Nucl. Instrum. Methods Phys. Res., Sect. A 757. 28-32. doi: 10.1016/j.nima.2014.04.033
 */





#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif


/* GLOBAL CONSTANTS */
/*
static const double BOLTZMANN_CONSTANT_evK = 8.617343E-5; // in [eV/K]
static const double BOLTZMANN_CONSTANT_JK = 1.3806504E-23; // in [J/K]
static const double PLANCK_CONSTANT_eVs = 4.13566733E-15; // in [eVs]
static const double PLANCK_CONSTANT_Js = 6.62606896E-34; // in [Js]
static const double PLANCK_BAR_eVs = 6.58211899E-16; // in [eVs] (h_ = h/2pi)
static const double PLANCK_BAR_Js = 1.054571628E-34; // in [Js]
static const double EV2J = 1.6021773E-19; // 1eV = 1.6021773E-19J
static const double J2EV = 6.2415063E18;
static const double MASS_NEUTRON_kg = 1.674927351E-27; // [kg]
static const double AVOGADRO = 6.0221417930E23; // [mol-1]
static const double ATOMIC_MASS_U_kg = 1.6605402E-27; // atomic_mass in kg
*/



/**
 * \fn static int _dhkl_compare( const void *par1, const void *par2 )
 * \brief Compares to dhkl values.
 *
 * This function is used within nxs_initHKL() to sort the hkl lattice planes.
 *
 * @param par1
 * @param par2
 * @return int
 */
static int _dhkl_compare( const void *par1, const void *par2 )
{
  NXS_HKL* p1 = (NXS_HKL*)par1;
  NXS_HKL* p2 = (NXS_HKL*)par2;

  /* compare dhkl first */
  if( fabs(p1->dhkl-p2->dhkl) > 1.0e-6 ) return p1->dhkl < p2->dhkl ? 1 : -1;
  /* compare FSquare*multiplicity */
  else if( fabs(p1->FSquare*p1->multiplicity - p2->FSquare*p2->multiplicity) > 1.0e-6 ) return p1->FSquare*p1->multiplicity < p2->FSquare*p2->multiplicity ? 1 : -1;
  /* compare h, then k, then l */
  else if( p1->h!=p2->h ) return p2->h < p1->h ? 1 : -1;
  else if( p1->k!=p2->k ) return p2->k < p1->k ? 1 : -1;
  else return p2->l < p1->l ? 1 : -1;
}



//static int _equivhkl_compare( const void *par1, const void *par2 )
//{
//  int hkl1, hkl2;
//  hkl1 = abs( (int)1E6*((NXS_EquivHKL*)par1)->h ) + abs( (int)1E3*((NXS_EquivHKL*)par1)->k ) + abs( ((NXS_EquivHKL*)par1)->l );
//  hkl2 = abs( (int)1E6*((NXS_EquivHKL*)par2)->h ) + abs( (int)1E3*((NXS_EquivHKL*)par2)->k ) + abs( ((NXS_EquivHKL*)par2)->l );
//  return ( (hkl1<hkl2)  ? 1 : 0 );
//}



/**
 * \fn static double _distance( double x1, double y1, double z1, double x2, double y2, double z2 )
 * \brief Calculates the distance between two atom positions (x1,y1,z1) and (x2,y2,z2).
 *
 * This function is used within _generateWyckoffPositions().
 *
 * @param x1
 * @param y1
 * @param z1
 * @param x2
 * @param y2
 * @param z2
 * @return int
 */
static double _distance( double x1, double y1, double z1, double x2, double y2, double z2 )
{
  double dx = x2 - x1;
  double dy = y2 - y1;
  double dz = z2 - z1;

  return sqrt( dx*dx + dy*dy + dz*dz );
}


/**
 * \fn static int _generateWyckoffPositions( NXS_UnitCell *uc )
 * \brief Generates Wyckoff position of a unit cell.
 *
 * Using SgInfo library functions this function generates Wyckoff atom position depending on the symmetry of the unit cell,
 * given via NXS_UnitCell and on the atom information added to it before, e.g. via nxs_addAtomInfo().
 *
 * @param uc NXS_UnitCell struct
 * @return nxs error code
 */
static int _generateWyckoffPositions( NXS_UnitCell *uc )
{
  unsigned int nTrV,  nLoopInv, i,j,iLoopInv,iAtomInfo;
  const T_RTMx *lsmx;
  T_RTMx SMx;
  int iList;
  T_SgInfo SgInfo;

  SgInfo = uc->sgInfo;

  nTrV = SgInfo.LatticeInfo->nTrVector;

  nLoopInv = Sg_nLoopInv(&SgInfo);

  uc->nAtoms = 0;

  for( iAtomInfo=0; iAtomInfo<uc->nAtomInfo; iAtomInfo++ )
  {
    NXS_AtomInfo ai;
    double x0, y0, z0;
    int *TrV;

    ai = uc->atomInfoList[iAtomInfo];
    ai.nAtoms = 0;
    x0 = ai.x[0];
    y0 = ai.y[0];
    z0 = ai.z[0];

    TrV = SgInfo.LatticeInfo->TrVector;
    for( j=0; j<nTrV; j++, TrV+=3 )
    {
      for( iLoopInv=0; iLoopInv<nLoopInv; iLoopInv++ )
      {
        int f = -1;
        if( iLoopInv==0 ) f = 1;

        lsmx = SgInfo.ListSeitzMx;
        for( iList=0; iList<SgInfo.nList; iList++, lsmx++ )
        {
          double x,y,z;
          unsigned int l;
          int isExclusive;
          for( i=0; i<9; i++ )
            SMx.s.R[i] = f * lsmx->s.R[i];
          for( i=0; i<3; i++ )
            SMx.s.T[i] = iModPositive(f * lsmx->s.T[i] + TrV[i], STBF);

          x = (double)SMx.s.R[0]*x0 + (double)SMx.s.R[1]*y0 + (double)SMx.s.R[2]*z0 + (double)SMx.s.T[0]/(double)STBF;
          y = (double)SMx.s.R[3]*x0 + (double)SMx.s.R[4]*y0 + (double)SMx.s.R[5]*z0 + (double)SMx.s.T[1]/(double)STBF;
          z = (double)SMx.s.R[6]*x0 + (double)SMx.s.R[7]*y0 + (double)SMx.s.R[8]*z0 + (double)SMx.s.T[2]/(double)STBF;

          x = x<0.0 ? fmod(1.0-fmod(-x,1.0),1.0) : fmod(x,1.0);
          y = y<0.0 ? fmod(1.0-fmod(-y,1.0),1.0) : fmod(y,1.0);
          z = z<0.0 ? fmod(1.0-fmod(-z,1.0),1.0) : fmod(z,1.0);

          l = 0;
          isExclusive = 1;
          while( l<ai.nAtoms && isExclusive )
          {
            if( _distance(x,y,z,ai.x[l],ai.y[l],ai.z[l]) < 1E-6 )
              isExclusive = 0;
            else
              l++;
          }

          if( isExclusive )
          {
            ai.x[ai.nAtoms] = x;
            ai.y[ai.nAtoms] = y;
            ai.z[ai.nAtoms] = z;
            ai.nAtoms++;
          }
        }
      }
    }

    uc->nAtoms += ai.nAtoms;
    uc->atomInfoList[iAtomInfo] = ai;
  }
  return NXS_ERROR_OK;
}


/**
 * \fn nxs_calcFSquare( NXS_HKL *hklReflex, NXS_UnitCell* uc )
 * \brief Calculates the structure factor squared |F<sub>hkl</sub>|<sup>2</sup>.
 *
 * This function calculates the structure factor squared depending on the hkl Miller indices, given
 * via NXS_HKL, and on the crystal system, given via the NXS_UnitCell.
 *
 * @param hklReflex NXS_HKL struct
 * @param uc NXS_UnitCell struct
 * @return |F<sub>hkl</sub>|<sup>2</sup>
 */
double nxs_calcFSquare( NXS_HKL *hklReflex, NXS_UnitCell *uc )
{
  unsigned int i, j;
  int h = hklReflex->h;
  int k = hklReflex->k;
  int l = hklReflex->l;
  double dhkl = hklReflex->dhkl;

  double structure_factor_square = 0.0;
  double real = 0.0;
  double imag = 0.0;

  for( i=0; i<uc->nAtomInfo; i++ )
  {
    double exp_B_iso = exp( -uc->atomInfoList[i].B_iso/4.0/dhkl/dhkl );
    double sin_exp = 0.0;
    double cos_exp = 0.0;

    for( j=0; j<uc->atomInfoList[i].nAtoms; j++ )
    {
      double x = uc->atomInfoList[i].x[j];
      double y = uc->atomInfoList[i].y[j];
      double z = uc->atomInfoList[i].z[j];

      if( fabs(x+y+z)<1E-6 )
        cos_exp += 1.0;
      else
      {
        double exponent = 2.0*M_PI*(x*h+y*k+z*l);
        sin_exp += sin( exponent );
        cos_exp += cos( exponent );
      }
    }
    sin_exp *= exp_B_iso * uc->atomInfoList[i].b_coherent;
    cos_exp *= exp_B_iso * uc->atomInfoList[i].b_coherent;
    real += cos_exp;
    imag += sin_exp;
  }
  structure_factor_square = (real*real + imag*imag);

  return structure_factor_square;
}




/**
 * \fn double nxs_calcDhkl( int h, int k, int l, NXS_UnitCell *uc )
 * \brief Calculates the lattice spacing.
 *
 * This function calculates the lattice spacing in &Aring; depending on the hkl Miller indices and the crystal system,
 * given via the NXS_UnitCell.
 *
 * @param h Miller index h
 * @param k Miller index k
 * @param l Miller index l
 * @param uc NXS_UnitCell struct
 * @return d<sub>hkl</sub> [&Aring;]
 */
double nxs_calcDhkl( int h, int k, int l, NXS_UnitCell *uc )
{
  double d_spacing = 0.0;
  double a = uc->a;
  double b = uc->b;
  double c = uc->c;
  /* Next 3 lines modified by NCrystal developers (original lines missed the factor of "*M_PI/180"): */
  double alpha = uc->alpha*M_PI/180;
  double beta = uc->beta*M_PI/180;
  double gamma = uc->gamma*M_PI/180;
  double t1, t2, t3, t4, t5, t6;

  //   h = abs(h);
  //   k = abs(k);
  //   l = abs(l);

  switch( uc->crystalSystem )
  {
    /* XS_Cubic */
    case 7:   d_spacing = a / sqrt( (double)((h*h) + (k*k) + (l*l)) );
              break;
    /* XS_Hexagonal */
    case 6:   d_spacing = a / sqrt( 4.0/3.0*(h*h + h*k + k*k) + (a*a/(c*c) * l*l) );
              break;
    /* XS_Trigonal */
    case 5:   d_spacing = sqrt(3.0)*a*c / sqrt( 4.0*(h*h+k*k+h*k)*c*c + 3.0*l*l*a*a );
              break;
    /* XS_Tetragonal */
    case 4:   d_spacing = a / sqrt( h*h + k*k + a*a/(c*c)*l*l );
              break;
    /* XS_Orthorhombic */
    case 3:   d_spacing = 1.0 / sqrt( h*h/a/a + k*k/b/b + l*l/c/c );
              break;
    /* XS_Monoclinic */
    case 2:   d_spacing =  a*b*c*sqrt( 1.0-cos(beta)*cos(beta) ) / sqrt( b*b*c*c*h*h + a*a*c*c*k*k*sin(beta)*sin(beta) + a*a*b*b*l*l - 2.0*a*b*b*c*h*l*cos(beta) );
              break;
    /* XS_Triclinic */
    case 1:   t1 = b*b*c*c * sin(alpha);
              t2 = a*a*c*c * sin(beta);
              t3 = a*a*b*b * sin(gamma);
              t4 = a*b*c*c * ( cos(alpha)*cos(beta)-cos(gamma) );
              t5 = a*a*b*c * ( cos(beta)*cos(gamma)-cos(alpha) );
              t6 = b*b*c*c * ( cos(gamma)*cos(alpha)-cos(beta) );

              d_spacing = 1.0/uc->volume/uc->volume * ( t1*h*h + t2*k*k + t3*l*l + 2.0*t4*h*k + 2.0*t5*k*l + 2.0*t6*h*l );
              d_spacing = sqrt( 1.0 / d_spacing );
              break;
    /* XS_Unknown */
    case 0:   break;
    default:  break;
  }

  return d_spacing;
}




/**
 * \fn static double _calcPhi_1( double theta )
 * \brief Calculates &phi;<sub>1</sub>(&theta;) as shown by Vogel (2000).
 *
 * References:
 *  - S. Vogel (2000) Thesis, Kiel University, Germany.
 *
 * @param theta = T/&theta;<sub>D</sub>
 * @return &phi;<sub>1</sub>(&theta;)
 */
static double _calcPhi_1( double theta )
{
  double phi_1, a_n, step1, step2, n, riemann_zeta_2, I_m;

  step1 = 1.0;
  step2 = 0.0;
  n = 1.0;
  a_n = 0.0;
  while ( fabs(step1-step2)>1E-6 )
  {
    step1 = step2;
    step2 = 1.0/( exp(n/theta)*n*n );

    a_n += step2;
    n = n + 1.0;
  }

  riemann_zeta_2 = M_PI*M_PI/6.0;
  I_m = theta * log( 1.0-exp(-1.0/theta) ) + theta*theta * (riemann_zeta_2 - a_n);

  phi_1 = 0.5 + 2.0*I_m;

  return phi_1;
}




/**
 * \fn static double _calcPhi_3( double theta )
 * \brief Calculates &phi;<sub>1</sub>(&theta;) as shown by Vogel (2000).
 *
 * References:
 *  - S. Vogel (2000) Thesis, Kiel University, Germany.
 *
 * @param theta = T/&theta;<sub>D</sub>
 * @return &phi;<sub>3</sub>(&theta;)
 */
static double _calcPhi_3( double theta )
{
  double phi_3;

  double a_n, step1, step2, n, riemann_zeta_4, I_m;

  step1 = 1.0;
  step2 = 0.0;
  n = 1.0;
  a_n = 0.0;
  while ( fabs(step1-step2)>1E-6 )
  {
    step1 = step2;
    step2 = 1.0/( exp(n/theta)*n*n ) * ( 0.5+theta/n+theta*theta/n/n );

    a_n += step2;
    n = n + 1.0;
  }

  riemann_zeta_4 = M_PI*M_PI*M_PI*M_PI/90.0;
  I_m = theta * log( 1.0-exp(-1.0/theta) ) + 6.0*theta*theta * (riemann_zeta_4*theta*theta - a_n);

  phi_3 = 0.25 + 2.0*I_m;

  return phi_3;
}


/**
 * \fn int nxs_initUnitCell( NXS_UnitCell *uc )
 * \brief Initializes a unit cell.
 *
 * Based on a space group (uc->spaceGroup) and on the lattice parameters (uc->a,..., uc->alpha,...), given by the user,
 * this function initializes a unit cell and calculates its volume. If atoms, i.e. NXS_AtomInfo, have been added before,
 * e.g. by nxs_addAtomInfo(), then they are lost as this routine resets the unit cell. The unit cell temperature (sample
 * temperature) is set to 293 K. Any wish to change the temperature must be applied after calling this function.
 * Internally, the routine uses SgInfo library function to parse the space group information and obtain the crystal system.
 *
 * @param uc NXS_UnitCell struct
 * @return nxs error code
 */
int nxs_initUnitCell( NXS_UnitCell *uc )
{
  /* at first some initialization for SgInfo */
  T_SgInfo SgInfo;
  const T_TabSgName *tsgn = NULL;
  double a,b,c,alpha,beta,gamma;

  SgInfo.MaxList = 192;
  SgInfo.ListSeitzMx = (T_RTMx*)malloc( SgInfo.MaxList * sizeof(*SgInfo.ListSeitzMx) );/*"(T_RTMx*)" added by NCrystal developers for C++ compilation*/
  /* no list info needed here */
  SgInfo.ListRotMxInfo = NULL;

  if( isdigit(uc->spaceGroup[0]) )
  {
    tsgn = FindTabSgNameEntry(uc->spaceGroup, 'A');
    if (tsgn == NULL)
      return free(SgInfo.ListSeitzMx),NXS_ERROR_NOMATCHINGSPACEGROUP; /* no matching table entry *//* "free(SgInfo.ListSeitzMx)," added by NCrystal developers to fix potential leak detected by static code analysis */
    strncpy(uc->spaceGroup,tsgn->HallSymbol,MAX_CHARS_SPACEGROUP-1);  /* "-1" added by NCrystal developers */
  }

  /* initialize SgInfo struct */
  InitSgInfo( &SgInfo );
  SgInfo.TabSgName = tsgn;
  if ( tsgn )
    SgInfo.GenOption = 1;

  ParseHallSymbol( uc->spaceGroup, &SgInfo );
  CompleteSgInfo( &SgInfo );
  Set_si( &SgInfo );

  uc->crystalSystem = SgInfo.XtalSystem;
  uc->sgInfo = SgInfo;

  /* get the unit cell volume depending on crystal system*/
  uc->volume = 0.0;
  a = uc->a;
  b = uc->b;
  c = uc->c;
  /* Next 3 lines modified by NCrystal developers (original lines missed the factor of "*M_PI/180"): */
  alpha = uc->alpha*M_PI/180;
  beta = uc->beta*M_PI/180;
  gamma = uc->gamma*M_PI/180;

  /* V = a * b * c * sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2
                            + 2 * cos(alpha) * cos(beta) * cos(gamma))
  */
  switch( uc->crystalSystem )
  {
    /* XS_Cubic */
    case 7:   uc->volume = a*a*a;
              break;
    /* XS_Hexagonal */
    case 6:   uc->volume = 0.866025*a*a*c; /* = a*a*c * sin(M_PI/3.0) */
              break;
    /* XS_Trigonal */
    case 5:   uc->volume = 0.866025*a*a*c; /* = a*a*c * sin(M_PI/3.0) */
              break;
    /* XS_Tetragonal */
    case 4:   uc->volume = a*a*c;
              break;
    /* XS_Orthorhombic */
    case 3:   uc->volume = a*b*c;
              break;
    /* XS_Monoclinic */
    case 2:   uc->volume = a*b*c * sin(beta);
              break;
    /* XS_Triclinic */
    case 1:   uc->volume = a*b*c * sqrt( 1.0 - cos(alpha)*cos(alpha) - cos(beta)*cos(beta) - cos(gamma)*cos(gamma) + 2.0*cos(alpha)*cos(beta)*cos(gamma) );
              break;
    /* XS_Unknown */
    case 0:   break;
    default:  break;
  }

  /*
  printf( "\n# --------------------\n# ");
  PrintFullHM_SgName(tsgn, ' ', stdout);
  printf( "\n# Crystal system: %s\n", XS_Name[uc->crystalSystem] );
  printf( "# TSGName: %s %i %s %s\n", SgInfo.TabSgName->HallSymbol, SgInfo.TabSgName->SgNumber,
          SgInfo.TabSgName->Extension, SgInfo.TabSgName->SgLabels );
  printf( "# Centric: %i\n", SgInfo.Centric );
  printf( "# LoopInv: %i\n", Sg_nLoopInv(&SgInfo) );
  printf( "# nList: %i\n", SgInfo.nList );
  printf( "# OriginShift: %i %i %i\n", SgInfo.OriginShift[0], SgInfo.OriginShift[1], SgInfo.OriginShift[2]);
  printf( "# OrderL: %i\n", SgInfo.OrderL );
  printf( "# OrderP: %i\n", SgInfo.OrderP );
  printf( "# Point group: %s\n", PG_Names[PG_Index(SgInfo.PointGroup)] );
  printf( "# nGenerator: %i\n", SgInfo.nGenerator );
  printf( "# GeneratorList: %i %i %i %i\n", SgInfo.Generator_iList[0], SgInfo.Generator_iList[1],
          SgInfo.Generator_iList[2], SgInfo.Generator_iList[3] );
  */

  // some default values to start with
  uc->temperature = 293.0;
  uc->mass = 0.0;
  uc->density = 0.0;
  uc->nAtomInfo = 0;
  uc->atomInfoList = NULL;

  if( uc->mph_c2 > 1E-6 )
    uc->__flag_mph_c2 = 0;

  return NXS_ERROR_OK;
}


/* implementation of R factors by Th. Kittelmann after A.K. Freund (1983) Nucl. Instr. Meth. 213, 495-501 */
static double rfacts[23] = {0};

static void initRfacts()
{
  // Bernoulli numbers B_n (only for n = 0...23)
  double bernoulli[23] = { 1, -1./2, 1./6, 0., -1./30., 0., 1./42, 0.,
                           -1./30., 0., 5./66, 0., -691./2730., 0.,
                           7./6, 0., -3617./510., 0., 43867./798, 0.,
                           -174611./330., 0., 854513./138 };
  double nfact = 1.0;
  unsigned int n;
  for( n=0; n <= 22; ++n )
  {
    if (n)
      nfact *= n;
    rfacts[n] = bernoulli[n] / (nfact*(n+5.0/2.0));
  }
}

static double calcR(double x)
{
  double sum = 0.0;
  double xx = 1.0/x;
  unsigned int n;

  // This is done only once at the beginning
  if( !rfacts[0] )
    initRfacts();

  for (n=0;n<22;++n)
  {
    sum += xx * rfacts[n];
    xx *= x;
  }
  sum += 0.5 * xx * rfacts[22];
  return sum;
}


/**
 * \fn int nxs_addAtomInfo( NXS_UnitCell *uc, NXS_AtomInfo ai )
 * \brief Adds an atom to a unit cell.
 *
 * This function adds an NXS_AtomInfo to a NXS_UnitCell. The Wyckoff atom positions, the unit cell mass and phi_1,
 * phi3 and B_iso are calculated implicitly.
 *
 * @param uc NXS_UnitCell struct
 * @param ai NXS_AtomInfo struct
 * @return nxs error code
 */
int nxs_addAtomInfo( NXS_UnitCell *uc, NXS_AtomInfo ai )
{
  double x;

  /* if debyeTemp was not given as a parameter */
  if( uc->debyeTemp<1E-6 )
  {
    // maybe a calculation/estimation of debyeTemp from atom-specific Debye temperatures
    // is possible in future; then such a calculation should be called here
    // Instead a standard value is set if no debyeTemp was given before
    uc->debyeTemp = 300.0;
  }

  uc->atomInfoList = (NXS_AtomInfo*)realloc( uc->atomInfoList, sizeof(NXS_AtomInfo)*(uc->nAtomInfo+1) );
  if( !uc->atomInfoList )
    return NXS_ERROR_MEMORYALLOCATIONFAILED;

  uc->nAtomInfo++;
  //ai.M_m = ai.molarMass*ATOMIC_MASS_U_kg/MASS_NEUTRON_kg;
  ai.M_m = ai.molarMass * 0.99140954426;

  // Single phonon part after A.K. Freund (1983) Nucl. Instr. Meth. 213, 495-501
  ai.sph =  sqrt( uc->debyeTemp ) / ai.M_m;

  x = uc->debyeTemp/uc->temperature;
  if( x <= 6 )
    ai.sph *= calcR(x);
  else
    ai.sph *= 3.29708964927644 * pow(x,-3.5);

  // from S. Vogel (2000) Thesis, Kiel University
  ai.phi_1 = _calcPhi_1( 1/x );
  ai.B_iso = 5.7451121E3 * ai.phi_1 / ai.molarMass / uc->debyeTemp;
  ai.phi_3 = _calcPhi_3( 1/x );

  uc->atomInfoList[uc->nAtomInfo-1] = ai;
  _generateWyckoffPositions( uc );

  uc->mass += (double)(uc->atomInfoList[uc->nAtomInfo-1].nAtoms) * uc->atomInfoList[uc->nAtomInfo-1].molarMass;

  //uc->density = uc->mass / uc->volume / AVOGADRO * 1E24; // [g/cm^3]
  uc->density = uc->mass / uc->volume / 0.60221417930; // [g/cm^3]

  /* if mph_c2 was not given as a parameter */
  if( uc->__flag_mph_c2 )
  {
    // from A.K. Freund (1983) Nucl. Instr. Meth. 213, 495-501
    // this needs further investigation as the fit in Freund's paper is not sufficient for all atoms,
    //uc->mph_c2 = 4.27 * exp( uc->mass*ATOMIC_MASS_U_kg/MASS_NEUTRON_kg/uc->nAtoms / 61.0 );
    uc->mph_c2 = 4.27 * exp( uc->mass*0.01625261547/uc->nAtoms );
  }

  return NXS_ERROR_OK;
}


/**
 * \fn int nxs_initHKL( NXS_UnitCell *uc )
 * \brief Initializes the hkl lattice planes.
 *
 * This function uses SgInfo library functions to initialise the hkl planes and implicitly calculates
 * the unit cell density, the C2 multi-phonon constant (if not given by the user before), average
 * coherent and incoherent cross section value and stores equivalent hkl lattice planes (if any) per hkl.
 *
 * @param uc UnitCell struct
 * @return nxs error code
 */
/* fix_incoh_xs parameter added in next line NCrystal developers (see header file for explanation): */
int nxs_initHKL( NXS_UnitCell *uc, int fix_incoh_xs )
{
  unsigned int ai, i,j;
  double tmp;
  T_SgInfo SgInfo;
  int minH, minK, minL, max_hkl, restriction, h,k,l;
  NXS_HKL *hkl;
  unsigned long index_count;

  /*
  int pos;
  printf( "# Generated Positions:\n" );
  for( ai=0; ai<uc->nAtomInfo; ai++ )
    for( pos=0; pos<uc->atomInfoList[ai].nAtoms; pos++ )
      printf( "#   %s   %.3f, %.3f, %.3f\n", uc->atomInfoList[ai].label,uc->atomInfoList[ai].x[pos],
              uc->atomInfoList[ai].y[pos], uc->atomInfoList[ai].z[pos] );
  printf( "Rel. Unit cell mass: %f [g/mol]\n", uc->mass );
  printf( "Unit cell volume: %f [AA^3]\n", uc->volume );
  printf( "Unit cell density: %f [g/cm^3]\n", uc->density );
  */

  /* at first calculate average sigmaCoherent and sigmaIncoherent for unit cell */

  tmp = 0.0;
  uc->avgSigmaIncoherent = 0.0;
  uc->avgSigmaCoherent = 0.0;

  for( ai=0; ai<uc->nAtomInfo; ai++ )
  {
    uc->avgSigmaCoherent += uc->atomInfoList[ai].b_coherent * uc->atomInfoList[ai].nAtoms;
    uc->avgSigmaIncoherent += uc->atomInfoList[ai].sigmaIncoherent * uc->atomInfoList[ai].nAtoms;
    tmp += uc->atomInfoList[ai].b_coherent * uc->atomInfoList[ai].b_coherent * uc->atomInfoList[ai].nAtoms;
  }
  tmp = tmp / uc->nAtoms;
  uc->avgSigmaIncoherent = uc->avgSigmaIncoherent / uc->nAtoms;
  uc->avgSigmaCoherent = uc->avgSigmaCoherent / uc->nAtoms;
  uc->avgSigmaCoherent = uc->avgSigmaCoherent * uc->avgSigmaCoherent;

  /* Check on fix_incoh_xs in next line added by NCrystal developers (see header file for explanation): */
  if ( ! fix_incoh_xs )
  uc->avgSigmaIncoherent += 0.04*M_PI * ( tmp - uc->avgSigmaCoherent );
  uc->avgSigmaCoherent *= 0.04*M_PI;

  /* some initialization for SgInfo */
  SgInfo = uc->sgInfo;

  /* start calculation of permitted reflections and multiplicities */
  max_hkl = uc->maxHKL_index;

  SetListMin_hkl( &SgInfo, max_hkl, max_hkl, &minH, &minK, &minL );

  /* how much hkl indices */
  /* Next line modified by NCrystal developers (original lines had minH instead of both minK and minL): */
  index_count = (max_hkl-minH+1)*(max_hkl-minK+1)*(max_hkl-minL+1);
  hkl = (NXS_HKL*)malloc(  sizeof(NXS_HKL)*index_count );
  if( !hkl )
    return NXS_ERROR_MEMORYALLOCATIONFAILED;

  i = 0;

  /* initialize all hkl indices */
  //   for (h=minH; h<=max_hkl; h++)
  //   for (k=minK; k<=max_hkl; k++)
  //   for (l=minL; l<=max_hkl; l++)
  for( h=max_hkl; h>=minH; h-- )
  for( k=max_hkl; k>=minK; k-- )
  for( l=max_hkl; l>=minL; l-- )
  {
    /* do not show hkls that are systematic absent for the space group */
    if( !IsSysAbsent_hkl( &SgInfo, h, k, l, &restriction ) )
    {
      char is_exclusive;
      /* exclude (hkl)=(000) */
      if( h==0 && k==0 && l==0 )
        continue;

      hkl[i].h = h;
      hkl[i].k = k;
      hkl[i].l = l;

      /* check if equivalent plane has been found before (and calculated) */
      is_exclusive = 1;
      for( j=0; j<i; j++ )
      {
        if( AreSymEquivalent_hkl(&SgInfo, h, k, l, hkl[j].h, hkl[j].k, hkl[j].l) )
        {
          is_exclusive = 0;
          /* sort hkl indices */
          if( 1E6*h+1E3*k+l > 1E6*hkl[j].h+1E3*hkl[j].k+hkl[j].l )
          {
            hkl[j].h = h;
            hkl[j].k = k;
            hkl[j].l = l;
          }
        }
        if( !is_exclusive )
          continue;
      }

      if( is_exclusive )
        i++;
    }
  }
  /* reduce the allocated memory */
  uc->nHKL = i;
  /* problems with realloc call detected by static code analysis was fixed by NCrystal developers, changing 3 lines...
  hkl = (NXS_HKL*)realloc( hkl, sizeof(NXS_HKL)*uc->nHKL );
  if( !hkl )
    return NXS_ERROR_MEMORYALLOCATIONFAILED;
  ... to: */
  NXS_HKL* realloc_hkl = (NXS_HKL*)realloc( hkl, (uc->nHKL?sizeof(NXS_HKL)*uc->nHKL:1) );
  if( !realloc_hkl ) {
    free(hkl);
    return NXS_ERROR_MEMORYALLOCATIONFAILED;
  }
  hkl = realloc_hkl;

  for( i=0; i<uc->nHKL; i++ )
  {
    /* store the equivalent lattice plane (hkl) */
    T_Eq_hkl eqHKL;
    NXS_EquivHKL *equivHKL;
    unsigned int nEqHKL;

    hkl[i].multiplicity = BuildEq_hkl( &SgInfo, &eqHKL, hkl[i].h, hkl[i].k, hkl[i].l );
    equivHKL = (NXS_EquivHKL*)malloc( sizeof(NXS_EquivHKL)*eqHKL.N );
    if( !equivHKL )
      return free(hkl),NXS_ERROR_MEMORYALLOCATIONFAILED;/* "free(hkl)," added by NCrystal developers to fix potential leak detected by static code analysis */


    nEqHKL = eqHKL.N;
    for( j=0; j<nEqHKL; j++ )
    {
      equivHKL[j].h = eqHKL.h[j];
      equivHKL[j].k = eqHKL.k[j];
      equivHKL[j].l = eqHKL.l[j];
    }
    hkl[i].equivHKL = equivHKL;

    /* get d-spacing and |F|^2 */
    hkl[i].dhkl = nxs_calcDhkl( hkl[i].h, hkl[i].k, hkl[i].l, uc );
    hkl[i].FSquare = nxs_calcFSquare( &(hkl[i]), uc );
  }
  /* end of initalizing */

  /* sort hkl lattice planes by d_hkl */
  qsort( hkl, uc->nHKL, sizeof(NXS_HKL), _dhkl_compare );

  uc->hklList = hkl;
  return NXS_ERROR_OK;
}





/**
 * \fn double nxs_CoherentElastic( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the coherent elastic scattering cross section.
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>coh_el</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_CoherentElastic( double lambda, NXS_UnitCell* uc )
{
  NXS_HKL *hkl = uc->hklList;

  double xsect_coh_el = 0.0;

  /* for all d-spacings... */
  unsigned int i;
  for( i=0; i<uc->nHKL; ++i )
  {
    double delta = lambda - 2.0*hkl[i].dhkl;
    if( delta<1E-6 )
    {
      /* calculate the elastic coherent cross section */
      xsect_coh_el += hkl[i].FSquare * hkl[i].multiplicity * hkl[i].dhkl;
    }
  }
  /* this is our final coherent elastic scattering cross section */
  xsect_coh_el = xsect_coh_el*1E-2 * lambda*lambda / (2.0*uc->volume);

  return xsect_coh_el;
}




/**
 * \fn double nxs_Absorption( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the absorbtion cross section.
 *
 * This function calculates the absorbtion cross section on the base of a given reference value at a
 * neutron velocity of \f$v = 2200 m/s\f$, i.e. 1.798 &Aring;
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>abs</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_Absorption( double lambda, NXS_UnitCell* uc )
{
  double xsect_abs, sigma = 0.0;
  unsigned int i;
  for( i=0; i<uc->nAtomInfo; i++ )
  {
    sigma += uc->atomInfoList[i].sigmaAbsorption * uc->atomInfoList[i].nAtoms;
  }
  xsect_abs = sigma / 1.798 * lambda;

  return xsect_abs;
}




/**
 * \fn double nxs_IncoherentElastic( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the incoherent elastic scattering cross section.
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>inc_el</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_IncoherentElastic( double lambda, NXS_UnitCell* uc )
{
  double xsect_inc_el, s_el_inc = 0.0;
  unsigned int i;
  for( i=0; i<uc->nAtomInfo; i++ )
  {
    double value = lambda*lambda / 2.0 / uc->atomInfoList[i].B_iso;
    s_el_inc += value * ( 1.0 - exp(-1.0/value) ) * uc->atomInfoList[i].nAtoms;
  }

  xsect_inc_el = s_el_inc * uc->avgSigmaIncoherent;
  return xsect_inc_el;
}




/**
 * \fn double nxs_IncoherentInelastic( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the incoherent inelastic scattering cross section.
 *
 * An approximation of the inelastic part by the incoherent cross section is used here.
 * This approach works sufficiently for thermal neutrons as shown by Binder (1970).
 *
 * References:
 *  - K. Binder (1970) Phys. Stat. Sol. 41, 767-779.
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>inc_inel</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_IncoherentInelastic( double lambda, NXS_UnitCell* uc )
{
  unsigned int i;
  double phi1_phi3, value, s_el_inc, A,  s_total_inc;
  double xsect_inc_inel = 0.0;

  for( i=0; i<uc->nAtomInfo; i++ )
  {
    phi1_phi3 = uc->atomInfoList[i].phi_1 * uc->atomInfoList[i].phi_3;
    value = lambda*lambda / 2.0 / uc->atomInfoList[i].B_iso;
    s_el_inc = value * ( 1.0 - exp(-1.0/value) );
    A = uc->atomInfoList[i].M_m;
    s_total_inc = A/(A+1.0)*A/(A+1.0) * ( 1.0 + 9.0 * phi1_phi3 * value/A/A );
    xsect_inc_inel += (s_total_inc - s_el_inc) * uc->atomInfoList[i].nAtoms;
  }

  xsect_inc_inel = xsect_inc_inel * uc->avgSigmaIncoherent;
  return xsect_inc_inel;
}




/**
 * \fn double nxs_CoherentInelastic( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the coherent inelastic scattering cross section.
 *
 * An approximation of the inelastic part by the incoherent cross section is used here.
 * This approach works sufficiently for thermal neutrons as shown by Binder (1970)
 *
 * References:
 *  - K. Binder (1970) Phys. Stat. Sol. 41, 767-779.
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>coh_inel</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_CoherentInelastic( double lambda, NXS_UnitCell* uc )
{
  unsigned int i;
  double phi1_phi3, value, s_el_inc, A,  s_total_inc;
  double xsect_coh_inel = 0.0;

  for( i=0; i<uc->nAtomInfo; i++ )
  {
    phi1_phi3 = uc->atomInfoList[i].phi_1 * uc->atomInfoList[i].phi_3;
    value = lambda*lambda / 2.0 / uc->atomInfoList[i].B_iso;
    s_el_inc = value * ( 1.0 - exp(-1.0/value) );
    A = uc->atomInfoList[i].M_m;
    s_total_inc = A/(A+1.0)*A/(A+1.0) * ( 1.0 + 9.0 * phi1_phi3 * value/A/A );
    xsect_coh_inel += (s_total_inc - s_el_inc) * uc->atomInfoList[i].nAtoms;
  }
  xsect_coh_inel = xsect_coh_inel * uc->avgSigmaCoherent;
  return xsect_coh_inel;
}



/**
 * \fn double nxs_TotalInelastic( double lambda, NXS_UnitCell* uc )
 * This function is provided for convenience and to be downward compatible to earlier versions of the nxs library.
 * It currently links to nxs_TotalInelastic_BINDER().
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>tot_inel</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_TotalInelastic( double lambda, NXS_UnitCell* uc )
{
  return nxs_TotalInelastic_BINDER( lambda, uc);
}


/**
 * \fn double nxs_TotalInelastic_BINDER( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the total inelastic scattering cross section.
 *
 * It is faster than calculating nxsCoherentInelastic() and nxsIncoherentInelastic() each, but
 * uses the same approach of approximating the individual inelastic parts by the incoherent cross
 * section.
 * This approach works sufficiently for thermal neutrons as shown by Binder (1970), but is actually
 * rather weak for the cold neutron region.
 *
 * References:
 *  - K. Binder (1970) Phys. Stat. Sol. 41, 767-779.
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>tot_inel_BINDER</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_TotalInelastic_BINDER( double lambda, NXS_UnitCell* uc )
{
  unsigned int i;
  double phi1_phi3, value, s_el_inc, A,  s_total_inc;
  double xsect_total_inel = 0.0;

  for( i=0; i<uc->nAtomInfo; i++ )
  {
    phi1_phi3 = uc->atomInfoList[i].phi_1 * uc->atomInfoList[i].phi_3;
    value = lambda*lambda / 2.0 / uc->atomInfoList[i].B_iso;
    s_el_inc = value * ( 1.0 - exp(-1.0/value) );
    A = uc->atomInfoList[i].M_m;
    s_total_inc = A/(A+1.0)*A/(A+1.0) * ( 1.0 + 9.0 * phi1_phi3 * value/A/A );
    xsect_total_inel += (s_total_inc - s_el_inc) * uc->atomInfoList[i].nAtoms;
  }

  return xsect_total_inel * (uc->avgSigmaCoherent+uc->avgSigmaIncoherent);
}


/**
 * \fn double nxs_TotalInelastic_COMBINED( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the total inelastic neutron cross section as a combination of two approaches.
 *
 * This function calculates the total inelastic scattering cross section
 * as a combination of two different approaches by returning nxs_SinglePhonon() + nxs_MultiPhonon_COMBINED().
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>tot_inel_COMBINED</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_TotalInelastic_COMBINED( double lambda, NXS_UnitCell* uc )
{
  return nxs_SinglePhonon( lambda, uc ) + nxs_MultiPhonon_COMBINED( lambda, uc );
}


/**
 * \fn double nxs_SinglePhonon( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the single phonon part of the total inelastic neutron cross section.
 *
 * This function calculates single phonon part of the total inelastic scattering cross section
 * per unit cell after Freund (1983).
 *
 * References:
 *  - A.K. Freund (1983) Nucl. Instr. Meth. 213, 495-501
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>sph</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_SinglePhonon( double lambda, NXS_UnitCell* uc )
{
  unsigned int i = 0;

  // Angstrom to eV
  // double neutronEnergy = PLANCK_CONSTANT_eVs*PLANCK_CONSTANT_eVs/(2*lambda*lambda*1E-20*MASS_NEUTRON_kg * J2EV);
  double neutronEnergy = 8.18042531017E-2 / lambda / lambda;

  double sph = 0.0;
  for( i=0; i<uc->nAtomInfo; i++ )
    sph += uc->atomInfoList[i].sph * uc->atomInfoList[i].nAtoms;

  return sph * (uc->avgSigmaCoherent+uc->avgSigmaIncoherent) / 35.90806936252971 / sqrt(neutronEnergy);
}


/**
 * \fn double nxs_MultiPhonon( double lambda, NXS_UnitCell* uc )
 * This function is provided for convenience and to be downward compatible to earlier versions of the nxs library.
 * It currently links to nxs_MultiPhonon_CASSELS().
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>mph</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_MultiPhonon( double lambda, NXS_UnitCell* uc )
{
  return nxs_MultiPhonon_CASSELS( lambda, uc );
}


/**
 * \fn double nxs_MultiPhonon_CASSELS( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the multiple phonon part of the total inelastic neutron cross section.
 *
 * This function calculates multiple phonon part of the total inelastic scattering cross section
 * per unit cell after Cassels (1950).
 *
 * References:
 *  - J.M. Cassels (1950) Prog. Nucl. Phys. 1, 185
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>mph_CASSELS</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_MultiPhonon_CASSELS( double lambda, NXS_UnitCell* uc )
{
  unsigned int i = 0;
  double mph = 0.0;
  for( i=0; i<uc->nAtomInfo; i++ )
  {
    double value = lambda*lambda / 2.0 / uc->atomInfoList[i].B_iso;
    double s_el_inc = value * ( 1.0 - exp(-1.0/value) );
    double A = uc->atomInfoList[i].M_m;
    mph += A/(A+1.0)*A/(A+1.0) * (1 - s_el_inc) * uc->atomInfoList[i].nAtoms;
  }
  return mph * (uc->avgSigmaCoherent+uc->avgSigmaIncoherent);//TK *= -> * (to avoid clang-tidy report)
}


// /* implementation of phi(x) by Th. Kittelmann after A.K. Freund (1983) Nucl. Instr. Meth. 213, 495-501 */
//static double phi_integrand(double s)
//{
//  //evaluate s/(exp(s)-1) with correct limit at s->0 to at least ?? significant figures
//  if (s>0.05||s<-0.05)
//    return s/(exp(s)-1);

//  return 1-0.5*s+s*s*(1.0/12.0)-s*s*s*s*(1.0/720)+s*s*s*s*s*s*(1.0/30240.0)-s*s*s*s*s*s*s*s*(1.0/1209600.);
//  //Refactored slightly faster:
//  //double s2=s*s;
//  //return 1.0+s*(-0.5+s*((1/12.0)+s2*((-1.0/720)+s2*((1.0/30240.0)-s2*(1.0/1209600.)))));
//}

//static double phi(double x)
//{
//  //Numerical integration of s/(exp(-s)-1) from 0 to x, divided by x. Does not
//  //need to be particularly fast as it is evaluated only once per material, so
//  //we use a simple application of integration via Simpson's rule (could easily
//  //be optimised further):
//  const unsigned n = 10000;
//  const double dx = x / n;
//  double sum = 0;
//  unsigned int i;

//  if (x==0)
//    return 1;//limiting value

//  for ( i=0; i<n; ++i)
//    sum += phi_integrand(i*dx)+4*phi_integrand((i+0.5)*dx)+phi_integrand((i+1)*dx);

//  return sum * dx/( 6.0 * x );//notice the minus to get the correct form
//}

/**
 * \fn double nxs_MultiPhonon_FREUND( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the multiple phonon part of the total inelastic neutron cross section.
 *
 * This function calculates multiple phonon part of the total inelastic scattering cross section
 * per unit cell after Freund (1983).
 *
 * References:
 *  - A.K. Freund (1983) Nucl. Instr. Meth. 213, 495-501
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>mph_FREUND</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_MultiPhonon_FREUND( double lambda, NXS_UnitCell* uc )
{
  // calculation of BT an B0 is not neccessary here, as BT + B0 = Biso
  // hence, the phi(x) calculation is obsolete

  unsigned int i = 0;
  double mph = 0.0;

  double C2_times_Energy = uc->mph_c2 *  8.18042531017E-2/lambda/lambda;
  for( i=0; i<uc->nAtomInfo; i++ )
  {
    double A = uc->atomInfoList[i].M_m;

    /* BT+B0=Biso, i.e. no need to calculate BT and B0 here */
    // double B0 = 2872.568576905348 / (A*uc->atomInfoList[i].debyeTemp);
    // double x = uc->atomInfoList[i].debyeTemp / uc->temperature;
    // double BT = 4.0 * B0 * phi(x) / x;

    double mphPerAtom = A*A / ( (1.0+A)*(1.0+A) ) * ( 1.0 - exp( -uc->atomInfoList[i].B_iso * C2_times_Energy) );
    mph += mphPerAtom * uc->atomInfoList[i].nAtoms;
  }
  return mph * (uc->avgSigmaCoherent+uc->avgSigmaIncoherent);
}


/**
 * \fn double nxs_MultiPhonon_COMBINED( double lambda, NXS_UnitCell* uc )
 * \brief Calculates the total inelastic neutron cross section as a combination of two approaches.
 *
 * This function calculates the multiple phone part of the total inelastic scattering cross section
 * as a combination of two different approaches. For neutron with $E > k_B \theta_D$ an incoherent
 * approximation for the inelastic scattering part (Cassels[1950]) is used. For $E < k_B \theta_D$ Freund's (1983)
 * multiple phonon cross section calculation is applied.
 * Here, a linear switchover in the region around $E = k_B \theta_D$ is computed between nxs_MultiPhonon_CASSELS() and
 * nxs_MultiPhonon_FREUND(). This approach was suggested by Adib (2007), for example.
 *
 * References:
 *  - A.K. Freund (1983) Nucl. Instr. Meth. 213, 495-501
 *  - J.M. Cassels (1950) Prog. Nucl. Phys. 1, 185
 *  - M. Adib & M. Fathalla (2008) Egyptian Nuclear Physics Association (ENPA), Egyptian Atomic Energy Authority (EAEA), 642 p, 39-54,
 *  NUPPAC 07; 6th Conference on Nuclear and Particle Physics; Luxor (Egypt); 17-21 Nov 2007
 *
 * @param lambda wavelength in &Aring;
 * @param uc NXS_UnitCell struct
 * @return &sigma;<sub>mph_COMBINED</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_MultiPhonon_COMBINED( double lambda, NXS_UnitCell* uc )
{
  double lambda_debye = 30.8106673293723 / sqrt( uc->debyeTemp );

  /* these limits are estimates from the graphite case in Adib (2007) to provide a universal approach (for all materials) */
  double lambda_Cassels = lambda_debye * 1.78789683887;
  double lambda_Freund = lambda_debye * 3.68096408002;

  if( lambda>=lambda_Freund )
    return nxs_MultiPhonon_FREUND( lambda, uc );
  else if( lambda<=lambda_Cassels )
    return nxs_MultiPhonon_CASSELS( lambda, uc );
  else
  {
    /* this is the linear switchover proposed by by Adib (2007) */
    return (lambda_Freund-lambda)/(lambda_Freund-lambda_Cassels)*nxs_MultiPhonon_CASSELS( lambda, uc ) + (lambda-lambda_Cassels)/(lambda_Freund-lambda_Cassels)*nxs_MultiPhonon_FREUND( lambda, uc );
  }
}


/**
 * \fn NXS_MarchDollase nxs_initMarchDollase( NXS_UnitCell* uc )
 * \brief Initialises a NXS_MarchDollase struct.
 *
 * @param uc NXS_UnitCell struct
 * @return NXS_MarchDollase struct
 */
NXS_MarchDollase nxs_initMarchDollase( NXS_UnitCell* uc )
{
  NXS_MarchDollase md;
  double n;
  unsigned int i;

  md.N = 101;
  md.M = 1001;
  md.nOrientations = 0;
  md.texture = NULL;
  md.unitcell = uc;

  /* calculate sin_phi and cos_phi */
  md.sin_phi = (double*)malloc( sizeof(double)*md.N );
  md.cos_phi = (double*)malloc( sizeof(double)*md.N );

  n = (double)md.N;
  for( i=0; i<n; i++ )
  {
    md.sin_phi[i] = sin( M_PI/n*(double)i - M_PI/2.0 );
    md.cos_phi[i] = cos( M_PI/n*(double)i - M_PI/2.0 );
  }
  return md;
}


/**
 * \fn void nxs_addTexture( NXS_MarchDollase* md, NXS_Texture texture )
 * \brief Adds texture to the unit cell via the NXS_MarchDollase struct.
 *
 * The function adds a preferred orientation (a,b,c) with weighting factors r and f, stored in texture,
 * to the unit cell via md and prepares for the March-Dollase texture calculation.
 * @param md NXS_MarchDollase struct
 * @param texture NXS_Texture struct
 */
void nxs_addTexture( NXS_MarchDollase* md, NXS_Texture texture )
{
  int a = texture.a;
  int b = texture.b;
  int c = texture.c;
  int h,k,l;
  unsigned int i,j;
  double cos_beta, sin_beta,  N_Nplus1, r, m;

  NXS_UnitCell *uc;
  NXS_HKL *hkl;

  if( a==0 && b==0 && c==0 )
    return;

  //   a = 1.0*texture.a + 0.0*texture.b + 0.0*texture.c;
  //   b = 0.0*texture.a + 0.0*texture.b - 1.0*texture.c;
  //   c = 0.0*texture.a + 1.0*texture.b + 0.0*texture.c;

  uc = md->unitcell;
  hkl = uc->hklList;

  /* calculate sin_beta and cos_beta for all planes */
  texture.sin_beta = (double**)malloc( sizeof(double*)*uc->nHKL );
  texture.cos_beta = (double**)malloc( sizeof(double*)*uc->nHKL );

  for( i=0; i<uc->nHKL; i++ )
  {
    unsigned int nEquivalent = hkl[i].multiplicity / 2;
    NXS_EquivHKL *equivHKL;

    if( nEquivalent==0 ) nEquivalent = 1;
    equivHKL = hkl[i].equivHKL;

    texture.sin_beta[i] = (double*)malloc( sizeof(double)*nEquivalent );
    texture.cos_beta[i] = (double*)malloc( sizeof(double)*nEquivalent );

    for( j=0; j<nEquivalent; j++ )
    {
      h = equivHKL[j].h;
      k = equivHKL[j].k;
      l = equivHKL[j].l;

      cos_beta = (double)(a*h+b*k+c*l) / sqrt((double)(a*a+b*b+c*c)) / sqrt((double)(h*h+k*k+l*l));
      if ( cos_beta > 1.0 )
        cos_beta = 1.0;
      else if ( cos_beta < -1.0 )
        cos_beta = -1.0;
      sin_beta = sin( acos(cos_beta) );

      texture.sin_beta[i][j] = sin_beta;
      texture.cos_beta[i][j] = cos_beta;
    }
  }

  /* calculate P_alpha_H in constant cos_alpha_H steps */
  N_Nplus1 = 1.0 /(double)(md->N+1.0);
  r = texture.r;
  texture.P_alpha_H = (double*)malloc( sizeof(double)*md->M );

  m = (double)md->M;
  for ( i=0; i<md->M; i++ )
  {
    double  cos_alpha_H, sin_alpha_H;
    cos_alpha_H = 2.0/m*i - 1.0;
    sin_alpha_H = sin( acos(cos_alpha_H) );
    texture.P_alpha_H[i] = pow( r*r*cos_alpha_H*cos_alpha_H + sin_alpha_H*sin_alpha_H/r, -1.5 ) * N_Nplus1;
  }

  md->nOrientations++;
  md->texture = (NXS_Texture*)realloc( md->texture, sizeof(NXS_Texture)*md->nOrientations );
  md->texture[md->nOrientations-1] = texture;
}


/**
 * \fn double nxs_CoherentElasticTexture( double lambda, NXS_MarchDollase* md )
 * \brief Calculates the coherent elastic scattering cross section with texture influence.
 *
 * The function computes the coherent elastic scattering cross section but also applies a
 * March-Dollase texture correction. The unit cell information must be stored to the
 * NXS_MarchDollase struct.
 * @param lambda wavelength in &Aring;
 * @param md NXS_MarchDollase struct
 * @return &sigma;<sub>coh_el_texture</sub> [barn = 10<sup>-24</sup> cm<sup>2</sup>]
 */
double nxs_CoherentElasticTexture( double lambda, NXS_MarchDollase* md )
{
  double xsect_coh_el = 0.0;
  // double energy = 8.18042531017E-2/(lambda*lambda);

  //  double x = 1.0;
  //  double y = 0.0;
  //  double z = 0.0;
  //  double length_xyz = sqrt( x*x + y*y + z*z );
  //   x = x/length_xyz;
  //   y = y/length_xyz;
  //   z = z/length_xyz;

  NXS_UnitCell *uc = md->unitcell;
  NXS_HKL *hkl = uc->hklList;
  NXS_Texture *texture = md->texture;
  unsigned int nEquivalent;

  /* for all d-spacings... */
  unsigned int i, j, k, l;

  for( i=0; i<uc->nHKL; i++ )
  {
    double delta = lambda - 2.0*hkl[i].dhkl;
    if( delta < -1E-6 )
    {
      double alpha_h, sin_alpha_h, cos_alpha_h,  cos_beta, sin_beta;
      double corr = 0.0;
      nEquivalent = hkl[i].multiplicity / 2;
      if( nEquivalent==0 ) nEquivalent = 1;

      /* calculate alpha_h */
      alpha_h = M_PI/2.0 - asin( lambda/2.0/hkl[i].dhkl );
      sin_alpha_h = sin( alpha_h );
      cos_alpha_h = cos( alpha_h );

      /* calculate the correction factor... */
      for( j=0; j<md->nOrientations; j++ )
      {
        for( k=0; k<nEquivalent; k++ )
        {
          cos_beta = texture[j].cos_beta[i][k];
          sin_beta = texture[j].sin_beta[i][k];
          for( l=0; l<md->N; l++ )
          {
            double sin_phi, a, cos_alpha_H;
            unsigned int index;
            sin_phi = md->sin_phi[l];
            // double cos_phi = md->cos_phi[l];

            a = cos_beta*cos_alpha_h - sin_beta*sin_alpha_h*sin_phi;
            //             double b = sin_beta*cos_phi;
            //             double c = cos_beta*sin_alpha_h + sin_beta*cos_alpha_h*sin_phi;
            cos_alpha_H = a/* + b + c*/;
            //             cos_alpha_H = a*x+b*y+c*z;

            index = (int)( (1.0+cos_alpha_H)/2.0*(double)md->M );
            corr += md->texture[j].P_alpha_H[index] * md->texture[j].f;
          }
        } /* end of nEquivalent */
      } /* end of nOrientations */

      /* calculate the elastic coherent cross section with texture influence */
      corr = corr / (double)(nEquivalent);
      xsect_coh_el += hkl[i].FSquare * hkl[i].multiplicity * hkl[i].dhkl * corr;
    }
  } /* end of hkls */

  /* this is our final coherent elastic scattering cross section with texture influence*/
  xsect_coh_el = xsect_coh_el*1E-2 * lambda*lambda / (2.0*uc->volume); // uc->nAtoms;

  return xsect_coh_el;
}


/* keys for the nxs parameter file */
const char *NXS_keys[] =
{
  "space_group",
  "lattice_a",
  "lattice_b",
  "lattice_c",
  "lattice_alpha",
  "lattice_beta",
  "lattice_gamma",
  "add_atom",
  "mph_c2",
  "debye_temp"
};


/**
 * \fn int nxs_readParameterFile( const char* fileName, NXS_UnitCell *uc, NXS_AtomInfo *atomInfoList[] )
 * \brief Reads a nxs parameter file.
 *
 * The function reads the nxs parameter file fileName and stores the data to NXS_UnitCell and the atom
 * information to atomInfoList.
 * In case of a successful read the funtion returns the number of atom info line read from the file.
 * Otherwise, an error code is returned. Then, the funtion nxs_initUnitCell() should be called after this.
 * @param fileName
 * @param uc NXS_UnitCell struct
 * @param atomInfoList NXS_AtomInfo struct
 * @return nxs error code
 */
int nxs_readParameterFile( const char* fileName, NXS_UnitCell *uc, NXS_AtomInfo *atomInfoList[] )
{
  unsigned int numAtoms = 0;
  NXS_AtomInfo ai;
  NXS_AtomInfo *aiList = NULL;

  char line[200];
  unsigned int max_keys = sizeof(NXS_keys)/sizeof(NXS_keys[0]);

  /* open the parameter file */
  FILE* file;
  file = fopen(fileName, "r");

  *uc = nxs_newUnitCell();

  if (!file)
  {
    return NXS_ERROR_READINGFILE;
  }

  /* to make sure parameters are initialized if not in the parameter file */

  while( fgets(line, sizeof(line), file) != NULL )
  {
    /* make sure 1st character isn't a space */
    char *ptr = line;
    char *par;
    while( *ptr && isspace(*ptr) )
      ptr++;

    /* find parameter and value pair */
    par = strtok( ptr, "=" );
    if( par != NULL )
    {
      unsigned int i=0;
      /* remove possible spaces from the end of the parameter term */
      char *key = par;
      while( isspace(key[strlen(key)-1]) )
        key[strlen(key)-1] = '\0';

      while( i < max_keys && strcmp(key, NXS_keys[i]) )
        i++;
      /* parameter found, now check the value */
      if( i < max_keys )
      {
        char *endptr;
        unsigned int nwords = 0;
        char isinword = 0;

        par = strtok( NULL, "=#!;" );

        /* remove possible spaces from begin and end of the value term */
        while( *par && isspace(*par) )
          par++;
        while( isspace(par[strlen(par)-1]) )
          par[strlen(par)-1] = '\0';

        switch( i )
        {
          case 0:                                    /* space_group */
            strncpy(uc->spaceGroup,par,MAX_CHARS_SPACEGROUP-1); break;  /* "-1" added by NCrystal developers */
          case 1:                                    /* lattice_a */
            uc->a = strtod(par, &endptr); break;
          case 2:                                    /* lattice_b */
            uc->b = strtod(par, &endptr); break;
          case 3:                                    /* lattice_c */
            uc->c = strtod(par, &endptr); break;
          case 4:                                    /* lattice_alpha */
            uc->alpha = strtod(par, &endptr); break;
          case 5:                                    /* lattice_beta */
            uc->beta = strtod(par, &endptr); break;
          case 6:                                    /* lattice_gamma */
            uc->gamma = strtod(par, &endptr); break;
          case 7:                                    /* add_atom */
            /* first count the words in that line */
            endptr = par;
            nwords = 0;
            for( ; *endptr; ++endptr)
            {
              if( !isspace(*endptr) )
              {
                if( isinword == 0 )
                {
                    isinword = 1;
                    ++nwords;
                }
              }
              else if( isinword )
                isinword = 0;
            }
            /* check number of words */
            if( nwords<8 )
              return NXS_ERROR_READINGFILE;
            else
            {
              strncpy(ai.label,strtok(par, " \t"),MAX_CHARS_ATOMLABEL-1);     /* atom name or label */  /* "-1" added by NCrystal developers */
              ai.b_coherent = strtod(strtok( NULL, " \t" ), &endptr);       /* b_coh */
              ai.sigmaIncoherent = strtod(strtok( NULL, " \t" ), &endptr);  /* sigma_inc */
              ai.sigmaAbsorption = strtod(strtok( NULL, " \t" ), &endptr);  /* sigma_abs_2200 */
              ai.molarMass = strtod(strtok( NULL, " \t" ), &endptr);        /* molar_mass */

              /* for compatibility to earlier versions of nxs the atom specific entry of a */
              /* Debye temperature is also allowed. If not given, -1.0 is set and continue with x y z */
              if( nwords>8 )
                ai.debyeTemp = strtod(strtok( NULL, " \t" ), &endptr);      /* debye_temp */
              else
                ai.debyeTemp = -1.0;

              /* the Wyckoff postion of the atom inside the unit cell */
              ai.x[0] = strtod(strtok( NULL, " \t" ), &endptr);             /* x */
              ai.y[0] = strtod(strtok( NULL, " \t" ), &endptr);             /* y */
              ai.z[0] = strtod(strtok( NULL, " \t" ), &endptr);             /* z */

              /* collect all "add_atom" entries from file to make sure     */
              /* uc is initialised first before nxs_addAtomInfo is applied */
              numAtoms++;
              aiList = (NXS_AtomInfo*)realloc( aiList, sizeof(NXS_AtomInfo)*numAtoms );
              aiList[numAtoms-1] = ai;
            }
            break;
          case 8:                                    /* mph_c2 */
            uc->mph_c2 = strtod(par, &endptr);
            break;
          case 9:                                    /* global debye_temp */
            uc->debyeTemp = strtod(par, &endptr); break;
          default:  break;

        } /* end of: switch( i ) */
      } /* end of: if( i < max_keys ) */
    } /* end of: if( par != NULL ) */
  } /* end of: while( fgets(line, sizeof(line), file) != NULL ) */

  /* finished with reading file */
  fclose(file);
  if( !numAtoms )
    return NXS_ERROR_NOATOMINFOFOUND;

  *atomInfoList = aiList;
  return numAtoms;
}


/**
 * \fn int nxs_saveParameterFile( const char* fileName, NXS_UnitCell *uc )
 * \brief Saves the current unit cell settings to a nxs parameter file.
 *
 * @param fileName
 * @param uc NXS_UnitCell struct
 * @return nxs error code
 */
int nxs_saveParameterFile( const char* fileName, NXS_UnitCell *uc )
{
  unsigned int i = 0;

  /* open the parameter file */
  FILE* file;
  file = fopen(fileName, "w");

  if (!file)
  {
    return NXS_ERROR_SAVINGFILE;
  }

  fprintf( file, "#\n# This is an nxs parameter file\n#\n%s = %s\n%s = %f\n%s = %f\n%s = %f\n%s = %f\n%s = %f\n%s = %f\n%s = %f\n%s = %f\n\n# label  b_coherent  sigma_inc  sigma_abs  molar_mass  debye_temp  x  y  z\n",
        NXS_keys[0], uc->spaceGroup, NXS_keys[1], uc->a, NXS_keys[2], uc->b, NXS_keys[3], uc->c, NXS_keys[4], uc->alpha, NXS_keys[5], uc->beta, NXS_keys[6], uc->gamma, NXS_keys[8], uc->mph_c2, NXS_keys[9], uc->debyeTemp);

  for(i=0; i<uc->nAtomInfo; i++)
  {
    fprintf( file, "%s = %s %f %f %f %f ",
         NXS_keys[7], uc->atomInfoList[i].label, uc->atomInfoList[i].b_coherent, uc->atomInfoList[i].sigmaIncoherent,
         uc->atomInfoList[i].sigmaAbsorption, uc->atomInfoList[i].molarMass );

    /* if debye temp value is stored to atomInfo store it */
    /* this is maybe removed in the future */
    if( uc->atomInfoList[i].debyeTemp<1E-6 )
      fprintf( file, "%f ", uc->atomInfoList[i].debyeTemp );

    fprintf( file, "%f %f %f\n", uc->atomInfoList[i].x[0], uc->atomInfoList[i].y[0], uc->atomInfoList[i].z[0] );
  }

  fclose(file);

  return NXS_ERROR_OK;
}


/**
 * \fn const char* nxs_version()
 * \brief Returns the library version number.
 *
 * @return const char*
 */
const char* nxs_version()
{
#ifdef NXSLIB_VERSION
  #define VERNUM2STR(x) VERSTR(x)
  #define VERSTR(x) #x
  #define RETURNVALUE VERNUM2STR( NXSLIB_VERSION )
#else
  #define RETURNVALUE (const char*)"unknown"
#endif

  return RETURNVALUE;
}


/**
 * \fn NXS_UnitCell nxs_newUnitCell()
 * \brief Returns an empty but correctly initialized unit cell.
 *
 * This funtion is called by nxs_readParameterFile(). So, there is only need to call this function if the individual
 * unit cell parameters will be assigned manually, i.e. without the file reading. Like for nxs_readParameterFile(), after
 * the assigning the unit cell parameters the funtion nxs_initUnitCell() should be called.
 *
 * @return NXS_UnitCell
 */
NXS_UnitCell nxs_newUnitCell()
{
  NXS_UnitCell uc;
  /* Next line modified by NCrystal developers to null out ALL fields: */
  memset(&uc,0,sizeof(uc));
  uc.spaceGroup[0] = '\0';
  uc.a = 0.0;
  uc.b = 0.0;
  uc.c = 0.0;
  uc.alpha = 0.0;
  uc.beta = 0.0;
  uc.gamma = 0.0;
  /* actual values should never be smaller than 0, hence a -1.0 is set as a start value */
  uc.mph_c2 = -1.0;  // if smaller 0 then mph_c2 is computed during nxs_addAtomInfo()
  uc.debyeTemp = -1.0; // if smaller 0 then debyeTemp is set during nxs_addAtomInfo()

  uc.__flag_mph_c2 = 1;

  return uc;
}
/*
  Space Group Info's (c) 1994-96 Ralf W. Grosse-Kunstleve
 */





static const char *Err_Ill_SMx_in_List =
  "Error: Illegal SeitzMx in list";


void SetSgError(const char *msg)
{
  if (SgError == NULL) SgError = msg;
}


int iModPositive(int ix, int iy)
{
  if (iy > 0)
  {
    ix %= iy;
    if (ix < 0) ix += iy;
  }

  return ix;
}


static void SwapRTMx(T_RTMx *Mx_a, T_RTMx *Mx_b)
{
  int     i;
  T_RTMx  BufMx;


  for (i = 0; i < 12; i++) BufMx.a[i] = Mx_a->a[i];
  for (i = 0; i < 12; i++) Mx_a->a[i] = Mx_b->a[i];
  for (i = 0; i < 12; i++) Mx_b->a[i] = BufMx.a[i];
}


static void CopyRotMxInfo(T_RotMxInfo *target, const T_RotMxInfo *source)
{
  memcpy(target, source, sizeof (*target));
}


static void SwapRotMxInfo(T_RotMxInfo *RMx_a, T_RotMxInfo *RMx_b)
{
  T_RotMxInfo  Buffer;

  memcpy(&Buffer, RMx_a,   sizeof (Buffer));
  memcpy(RMx_a,   RMx_b,   sizeof (Buffer));
  memcpy(RMx_b,   &Buffer, sizeof (Buffer));
}


int traceRotMx(const int *RotMx)
{
  return RotMx[0] + RotMx[4] + RotMx[8];
}


int deterRotMx(const int *RotMx)
{
  int     det;

  det =  RotMx[0] * (RotMx[4] * RotMx[8] - RotMx[5] * RotMx[7]);
  det -= RotMx[1] * (RotMx[3] * RotMx[8] - RotMx[5] * RotMx[6]);
  det += RotMx[2] * (RotMx[3] * RotMx[7] - RotMx[4] * RotMx[6]);

  return det;
}


void RotMx_t_Vector(int *R_t_V, const int *RotMx, const int *Vector, int FacTr)
{
  const int  *vec;


  if (FacTr > 0)
  {
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec;
    *R_t_V   %= FacTr; if (*R_t_V < 0) *R_t_V += FacTr;
     R_t_V++;
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec;
    *R_t_V   %= FacTr; if (*R_t_V < 0) *R_t_V += FacTr;
     R_t_V++;
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V   += *RotMx   * *vec;
    *R_t_V   %= FacTr; if (*R_t_V < 0) *R_t_V += FacTr;
  }
  else
  {
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V++ += *RotMx++ * *vec;
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V++ += *RotMx++ * *vec;
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V   += *RotMx   * *vec;
  }
}


void RotMxMultiply(int *rmxab, const int *rmxa, const int *rmxb)
{
  const int  *a, *b;

  /* no loops to be as fast as posslible */

  a = rmxa;
  b = rmxb;
  *rmxab  = *a++ * *b; b += 3; /* r11 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r12 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r13 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a++ * *b; b = rmxb;
   rmxab++;

  rmxa = a;
  *rmxab  = *a++ * *b; b += 3; /* r21 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r22 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r23 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a++ * *b; b = rmxb;
   rmxab++;

  rmxa = a;
  *rmxab  = *a++ * *b; b += 3; /* r31 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r32 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r33 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b;
}


void RotateRotMx(int *RotMx, const int *RMx, const int *InvRMx)
{
  int  BufMx[9];


  RotMxMultiply(BufMx, RotMx, InvRMx);
  RotMxMultiply(RotMx, RMx,   BufMx);
}


void SeitzMxMultiply(T_RTMx *smxab, const T_RTMx *smxa, const T_RTMx *smxb)
{
  const int  *ar, *a, *b, *bt;
  int        *ab;

  /* no loops to be as fast as posslible */

  ar = smxa->a;
  a  = smxa->a;
  b  = smxb->a;
  ab = smxab->a;

  *ab  = *a++ * *b; b += 3; /* r11 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r12 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r13 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b; b = smxb->a;
   ab++;

  ar = a;
  *ab  = *a++ * *b; b += 3; /* r21 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r22 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r23 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b; b = smxb->a;
   ab++;

  ar = a;
  *ab  = *a++ * *b; b += 3; /* r31 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r32 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r33 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b++; bt = b;
   ab++;

  ar = smxa->a;
  *ab  = *ar++ * *b++; /* t1 */
  *ab += *ar++ * *b++;
  *ab += *ar++ * *b; b = bt;
  *ab += *a++;
  *ab %= STBF; if (*ab < 0) *ab += STBF;
   ab++;

  *ab  = *ar++ * *b++; /* t2 */
  *ab += *ar++ * *b++;
  *ab += *ar++ * *b; b = bt;
  *ab += *a++;
  *ab %= STBF; if (*ab < 0) *ab += STBF;
   ab++;

  *ab  = *ar++ * *b++; /* t3 */
  *ab += *ar++ * *b++;
  *ab += *ar   * *b;
  *ab += *a;
  *ab %= STBF; if (*ab < 0) *ab += STBF;
}


void RTMxMultiply(T_RTMx *rtmxab, const T_RTMx *rtmxa, const T_RTMx *rtmxb,
                  int FacAug, int FacTr)
{
  const int  *ar, *a, *b, *bt;
  int        *ab;

  /* no loops to be as fast as posslible */

  ar = rtmxa->a;
  a  = rtmxa->a;
  b  = rtmxb->a;
  ab = rtmxab->a;

  *ab  = *a++ * *b; b += 3; /* r11 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r12 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r13 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b; b = rtmxb->a;
   ab++;

  ar = a;
  *ab  = *a++ * *b; b += 3; /* r21 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r22 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r23 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b; b = rtmxb->a;
   ab++;

  ar = a;
  *ab  = *a++ * *b; b += 3; /* r31 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r32 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r33 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b++; bt = b;
   ab++;

  if (FacTr > 0)
  {
    ar = rtmxa->a;
    *ab  = *ar++ * *b++; /* t1 */
    *ab += *ar++ * *b++;
    *ab += *ar++ * *b; b = bt;
    *ab += *a++ * FacAug;
    *ab %= FacTr; if (*ab < 0) *ab += FacTr;
     ab++;

    *ab  = *ar++ * *b++; /* t2 */
    *ab += *ar++ * *b++;
    *ab += *ar++ * *b; b = bt;
    *ab += *a++ * FacAug;
    *ab %= FacTr; if (*ab < 0) *ab += FacTr;
     ab++;

    *ab  = *ar++ * *b++; /* t3 */
    *ab += *ar++ * *b++;
    *ab += *ar   * *b;
    *ab += *a   * FacAug;
    *ab %= FacTr; if (*ab < 0) *ab += FacTr;
  }
  else
  {
    ar = rtmxa->a;
    *ab  = *ar++ * *b++; /* t1 */
    *ab += *ar++ * *b++;
    *ab += *ar++ * *b; b = bt;
    *ab += *a++ * FacAug;
     ab++;

    *ab  = *ar++ * *b++; /* t2 */
    *ab += *ar++ * *b++;
    *ab += *ar++ * *b; b = bt;
    *ab += *a++ * FacAug;
     ab++;

    *ab  = *ar++ * *b++; /* t3 */
    *ab += *ar++ * *b++;
    *ab += *ar   * *b;
    *ab += *a   * FacAug;
  }
}


void InverseRotMx(const int *RotMx, int *InvRotMx)
{
  InvRotMx[0] =   RotMx[4] * RotMx[8] - RotMx[5] * RotMx[7];
  InvRotMx[1] = - RotMx[1] * RotMx[8] + RotMx[2] * RotMx[7];
  InvRotMx[2] =   RotMx[1] * RotMx[5] - RotMx[2] * RotMx[4];
  InvRotMx[3] = - RotMx[3] * RotMx[8] + RotMx[5] * RotMx[6];
  InvRotMx[4] =   RotMx[0] * RotMx[8] - RotMx[2] * RotMx[6];
  InvRotMx[5] = - RotMx[0] * RotMx[5] + RotMx[2] * RotMx[3];
  InvRotMx[6] =   RotMx[3] * RotMx[7] - RotMx[4] * RotMx[6];
  InvRotMx[7] = - RotMx[0] * RotMx[7] + RotMx[1] * RotMx[6];
  InvRotMx[8] =   RotMx[0] * RotMx[4] - RotMx[1] * RotMx[3];
}


void InverseRTMx(const T_RTMx *RTMx, T_RTMx *InvRTMx)
{
  int        *iR;
  const int  *T;


  iR = InvRTMx->s.R;

  InverseRotMx(RTMx->s.R, iR);

  T = RTMx->s.T;

  InvRTMx->s.T[0] = - iR[0] * T[0] - iR[1] * T[1] - iR[2] * T[2];
  InvRTMx->s.T[1] = - iR[3] * T[0] - iR[4] * T[1] - iR[5] * T[2];
  InvRTMx->s.T[2] = - iR[6] * T[0] - iR[7] * T[1] - iR[8] * T[2];
}


int IsSMxTransl0(const T_LatticeInfo *LatticeInfo, const int *SeitzMxT)
{
  int        iTrV, nTrV, t;
  const int  *TrV;


  nTrV = LatticeInfo->nTrVector;
   TrV = LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++)
  {
        t =     (SeitzMxT[0] + TrV[0]) % STBF;
    if (t == 0) {
          t =   (SeitzMxT[1] + TrV[1]) % STBF;
      if (t == 0) {
            t = (SeitzMxT[2] + TrV[2]) % STBF;
        if (t == 0)
          return 1;
    }}

    TrV += 3;
  }

  return 0;
}


static int IsSpecialSeitzMx(T_SgInfo *SgInfo, const T_RTMx *SMx, int ExpandLT)
{
  int                  i, special, smx11;
  const T_LatticeInfo  *ExpLT;


  /* check rotation part for identity or inversion operation */

         smx11 = SMx->s.R[0];
  if (   smx11 !=  1
      && smx11 != -1) return 0;

  for (i = 1; i < 9; i++)
  {
    if (i % 4) {
      if (SMx->s.R[i] !=     0) return 0; }
    else {
      if (SMx->s.R[i] != smx11) return 0; }
  }

  if (smx11 == 1) special = SpecialSMx_Identity;
  else            special = SpecialSMx_Inversion;

  /* rotation part is identity or inversion
     check translation part now
   */

  if (IsSMxTransl0(SgInfo->LatticeInfo, SMx->s.T) == 1)
    return (special | SpecialSMx_Transl0);

  if (ExpandLT && (smx11 == 1 || SgInfo->Centric))
  {
    /* try to expand lattice type */

    ExpLT = NULL;

    switch (SgInfo->LatticeInfo->Code)
    {
      case 'P':
        if (IsSMxTransl0(LI_A, SMx->s.T) == 1)
        { ExpLT = LI_A; break; }
        if (IsSMxTransl0(LI_B, SMx->s.T) == 1)
        { ExpLT = LI_B; break; }
        if (IsSMxTransl0(LI_C, SMx->s.T) == 1)
        { ExpLT = LI_C; break; }
        if (IsSMxTransl0(LI_I, SMx->s.T) == 1)
        { ExpLT = LI_I; break; }
        if (IsSMxTransl0(LI_R, SMx->s.T) == 1)
        { ExpLT = LI_R; break; }
        if (IsSMxTransl0(LI_S, SMx->s.T) == 1)
        { ExpLT = LI_S; break; }
        if (IsSMxTransl0(LI_T, SMx->s.T) == 1)
        { ExpLT = LI_T; break; }/* FALLTHRU *//*FALLTHRU added by NCrystal developers for gcc 7.2.1 */
      case 'A':
      case 'B':
      case 'C':
        if (IsSMxTransl0(LI_F, SMx->s.T) == 1)
          ExpLT = LI_F;
    }

    if (ExpLT != NULL)
    {
      SgInfo->LatticeInfo = ExpLT;
      SgInfo->StatusLatticeTr = -1;
      return (special | SpecialSMx_Transl0);
    }
  }

  return special;
}


int CompareSeitzMx(const T_LatticeInfo *LatticeInfo,
                   const T_RTMx *SeitzMxA, const T_RTMx *SeitzMxB)
{
  int        i, iTrV, nTrV, t;
  const int  *TrV;


  /* compare rotation part */

  for (i = 0; i < 9; i++)
    if (SeitzMxA->s.R[i] != SeitzMxB->s.R[i]) return 1;

  /* rotation part is same
     check translation
   */

  nTrV = LatticeInfo->nTrVector;
   TrV = LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
        t =     (SeitzMxA->s.T[0] + TrV[0]) % STBF;
    if (t ==     SeitzMxB->s.T[0]) {
          t =   (SeitzMxA->s.T[1] + TrV[1]) % STBF;
      if (t ==   SeitzMxB->s.T[1]) {
            t = (SeitzMxA->s.T[2] + TrV[2]) % STBF;
        if (t == SeitzMxB->s.T[2])
          return 0;
    }}
  }

  return -1;
}


int GetRotMxOrder(const int *RotMx)
{
  int deter = deterRotMx(RotMx);

  if (deter == -1 || deter == 1)
  {
    switch (traceRotMx(RotMx))
    {
      case -3:                  return -1;
      case -2:                  return -6;
      case -1: if (deter == -1) return -4;
               else             return  2;
      case  0: if (deter == -1) return -3;
               else             return  3;
      case  1: if (deter == -1) return -2;
               else             return  4;
      case  2:                  return  6;
      case  3:                  return  1;
    }
  }

  return 0;
}


static int nNextBasis_of_DirCode(const int DirCode,
                                 const int **RMx, const int **InvRMx)
{
  switch (DirCode)
  {
    case  '.': *RMx =            *InvRMx = NULL;      return 1;
    case  '=':
    case  '"':
    case '\'':
    case  '|':
    case '\\': *RMx = RMx_3_111; *InvRMx = RMx_3i111; return 3;
    case  '*': *RMx = RMx_4_001; *InvRMx = RMx_4i001; return 4;
    default:
      break;
  }

  SetSgError("Internal Error: Corrupt DirCode");
  return -1;
}


int GetRotMxInfo(const int *RotMx, T_RotMxInfo *RotMxInfo)
{
  int        i;
  int        nNextBasis, iNextBasis;
  int        nLoopInv, iLoopInv;
  int        Order, AbsOrder;
  int        RMxCopy[9], MatchMx[9], InvMatchMx[9], REV[3];
  int        *mmx;
  const int  *NBRMx, *InvNBRMx;

  const T_TabXtalRotMx  *txrmx;


  Order = GetRotMxOrder(RotMx);

  if (RotMxInfo)
      RotMxInfo->Order = Order;

  if (Order)
  {
    AbsOrder = abs(Order);

    if (Order > 0)
      for (i = 0; i < 9; i++) RMxCopy[i] =  RotMx[i];
    else
      for (i = 0; i < 9; i++) RMxCopy[i] = -RotMx[i];

    for (txrmx = TabXtalRotMx; txrmx->Order; txrmx++)
      if (txrmx->Order == AbsOrder) break;

    while (txrmx->Order == AbsOrder)
    {
      nNextBasis = nNextBasis_of_DirCode(txrmx->DirCode, &NBRMx, &InvNBRMx);

      if (nNextBasis < 0)
        return 0;

      if (AbsOrder > 2) nLoopInv = 2;
      else              nLoopInv = 1;

      for (i = 0; i < 9; i++) MatchMx[i] = txrmx->RMx[i];

      for (iNextBasis = 0; iNextBasis < nNextBasis; iNextBasis++)
      {
        if (iNextBasis)
          RotateRotMx(MatchMx, NBRMx, InvNBRMx);

        mmx = MatchMx;

        for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
        {
          if (iLoopInv)
          {
            InverseRotMx(MatchMx, InvMatchMx);
            mmx = InvMatchMx;
          }

          for (i = 0; i < 9; i++)
            if (mmx[i] != RMxCopy[i]) break;

          if (i == 9) /* matrix has been found */
          {
            if (RotMxInfo)
            {
              RotMxInfo->Inverse = iLoopInv;

              if (nNextBasis == 3)
              {
                switch(iNextBasis)
                {
                  case 0: RotMxInfo->RefAxis = 'z'; break;
                  case 1: RotMxInfo->RefAxis = 'x'; break;
                  case 2: RotMxInfo->RefAxis = 'y'; break;
                }
              }
              else
                RotMxInfo->RefAxis = 'o';

              RotMxInfo->DirCode = txrmx->DirCode;

              for (i = 0; i < 3; i++)
                RotMxInfo->EigenVector[i] = txrmx->EigenVector[i];

              for (; iNextBasis--;)
              {
                RotMx_t_Vector(REV, NBRMx, RotMxInfo->EigenVector, 0);

                if (iNextBasis-- == 0)
                {
                  for (i = 0; i < 3; i++)
                    RotMxInfo->EigenVector[i] = REV[i];

                  break;
                }

                RotMx_t_Vector(RotMxInfo->EigenVector, NBRMx, REV, 0);
              }
            }

            return Order;
          }
        }
      }

      txrmx++;
    }
  }

  return 0;
}


const T_RotMxInfo *ListOrBufRotMxInfo(const T_SgInfo *SgInfo, int iList,
                                      T_RotMxInfo *BufRotMxInfo)
{
  T_RotMxInfo  *RMxI;


      RMxI = SgInfo->ListRotMxInfo;
  if (RMxI)
      RMxI += iList;
  else
  {
    RMxI = BufRotMxInfo;

    if (GetRotMxInfo(SgInfo->ListSeitzMx[iList].s.R, RMxI) == 0) {
      SetSgError(Err_Ill_SMx_in_List);
      return NULL;
    }
  }

  return RMxI;
}


static int CoreAdd2ListSeitzMx(T_SgInfo *SgInfo, const T_RTMx *NewSMx)
{
  int                  i, iList;
  T_RTMx               *lsmx;
  T_RotMxInfo          RotMxInfo;
  const T_LatticeInfo  *LI;

  static const char *Err_NonXtalOp =
    "Error: Generators produce non-crystallographic operation";


  if (SgInfo->GenOption) LI = SgInfo->LatticeInfo;
  else                   LI = LI_P;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
    if (CompareSeitzMx(LI, NewSMx, lsmx) == 0)
      return 0; /* matrix is not unique */

  if (GetRotMxInfo(NewSMx->s.R, &RotMxInfo) == 0) {
    SetSgError(Err_NonXtalOp);
    return -1;
  }

  if (SgInfo->nList >= SgInfo->MaxList)
  {
    if (SgInfo->nList >= 192)
      SetSgError(Err_NonXtalOp);
    else
      SetSgError("Internal Error: Allocated space for ListSeitzMx too small");

    return -1;
  }

  for (i = 0; i < 12; i++) lsmx->a[i] = NewSMx->a[i];

  if (SgInfo->ListRotMxInfo != NULL)
    CopyRotMxInfo(&SgInfo->ListRotMxInfo[SgInfo->nList], &RotMxInfo);

  SgInfo->nList++;

  return 1;
}


int Add2ListSeitzMx(T_SgInfo *SgInfo, const T_RTMx *NewSMx)
{
  int     i, special, Enter;
  int     iMult, jMult;
  T_RTMx  TrialSMx, *iMultSMx, *jMultSMx;

  static const char *Err_Ill_lattice_translation =
    "Error: Illegal lattice translation";


  if (SgInfo->nList == 0)
  {
    /* make sure identity matrix is first in list */

    InitSeitzMx(&TrialSMx, 1);

    if (CoreAdd2ListSeitzMx(SgInfo, &TrialSMx) < 0)
      return -1;;
  }

  for (i = 0; i < 9; i++)
        TrialSMx.s.R[i] = NewSMx->s.R[i];

  for (i = 0; i < 3; i++) {
        TrialSMx.s.T[i] = NewSMx->s.T[i] % STBF;
    if (TrialSMx.s.T[i] < 0)
        TrialSMx.s.T[i] += STBF;
  }

  iMult = SgInfo->nList;
  iMultSMx = &SgInfo->ListSeitzMx[iMult];

  jMult = 1;
  jMultSMx = &SgInfo->ListSeitzMx[1]; /* skip first = identity matrix */

  for (;;)
  {
    Enter = 1;

    special = IsSpecialSeitzMx(SgInfo, &TrialSMx, 1);

    if      (special & SpecialSMx_Identity)
    {
      if (! (special & SpecialSMx_Transl0)) {
        SetSgError(Err_Ill_lattice_translation);
        return -1;
      }

      if (SgInfo->GenOption)
        Enter = 0;
    }
    else if (special & SpecialSMx_Inversion)
    {
      if (! (special & SpecialSMx_Transl0))
      {
        if (SgInfo->Centric) {
          SetSgError(Err_Ill_lattice_translation);
          return -1;
        }

        SgInfo->InversionOffOrigin = 1;
      }
      else
      {
        if (SgInfo->InversionOffOrigin)
          SgInfo->Centric = 1;

        SgInfo->InversionOffOrigin = 0;

        if (SgInfo->GenOption)
        {
          if (SgInfo->Centric == 0)
              SgInfo->Centric = -1;

          Enter = 0;
        }
        else
          SgInfo->Centric = 1;
      }
    }

    if (Enter && CoreAdd2ListSeitzMx(SgInfo, &TrialSMx) < 0)
      return -1;

    if (SgInfo->GenOption < 0)
      break;

    if (jMult > iMult)
    {
      iMult++;
      iMultSMx++;

      jMult = 1;
      jMultSMx = &SgInfo->ListSeitzMx[1]; /* skip first = identity matrix */
    }

    if (iMult == SgInfo->nList)
      break;

    SeitzMxMultiply(&TrialSMx, jMultSMx, iMultSMx);

    jMult++;
    jMultSMx++;
  }

  return 0;
}


int AddInversion2ListSeitzMx(T_SgInfo *SgInfo)
{
  T_RTMx  SeitzMx;

  InitSeitzMx(&SeitzMx, -1);
  return Add2ListSeitzMx(SgInfo, &SeitzMx);
}


int AddLatticeTr2ListSeitzMx(T_SgInfo *SgInfo,
                             const T_LatticeInfo *LatticeInfo)
{
  int        iTrV, nTrV;
  const int  *TrV;
  T_RTMx     SeitzMx;


  InitRotMx(SeitzMx.s.R, 1);

  nTrV = LatticeInfo->nTrVector;
   TrV = &LatticeInfo->TrVector[3]; /* skip first vector which is always 000 */

  for (iTrV = 1; iTrV < nTrV; iTrV++)
  {
    SeitzMx.s.T[0] = *TrV++;
    SeitzMx.s.T[1] = *TrV++;
    SeitzMx.s.T[2] = *TrV++;

    if (Add2ListSeitzMx(SgInfo, &SeitzMx) < 0)
      return -1;
  }

  if (SgInfo->GenOption)
    SgInfo->StatusLatticeTr = 0;
  else
    SgInfo->StatusLatticeTr = 1;

  return 0;
}


static int RemoveLatticeTr(T_SgInfo *SgInfo)
{
  int          iList, jList, i;
  T_RTMx       *smxi, *smxj, *lastsmx;
  T_RotMxInfo  *lrmxi, *lastrmxi;


  if (SgInfo->LatticeInfo->Code == 'P')
    return 0;

  if (SgInfo->StatusLatticeTr == -1)
  {
    if (AddLatticeTr2ListSeitzMx(SgInfo, SgInfo->LatticeInfo) < 0)
      return -1;
  }

  iList = 0;

  while (iList < SgInfo->nList)
  {
    smxi = &SgInfo->ListSeitzMx[iList];
    jList = iList + 1;

    while (jList < SgInfo->nList)
    {
      smxj = &SgInfo->ListSeitzMx[jList];

      if (CompareSeitzMx(SgInfo->LatticeInfo, smxi, smxj) == 0)
      {
        /* copy last element to this place */

        SgInfo->nList--;
        lastsmx = &SgInfo->ListSeitzMx[SgInfo->nList];
        for (i = 0; i < 12; i++) smxj->a[i] = lastsmx->a[i];

        if (SgInfo->ListRotMxInfo != NULL)
        {
          lrmxi =    &SgInfo->ListRotMxInfo[jList];
          lastrmxi = &SgInfo->ListRotMxInfo[SgInfo->nList];
          CopyRotMxInfo(lrmxi, lastrmxi);
        }
      }
      else
        jList++;
    }

    iList++;
  }

  SgInfo->StatusLatticeTr = 0;

  return 0;
}


static int RemoveInversion(T_SgInfo *SgInfo)
{
  int          i, iList, special, deter, Centric, InversionOffOrigin;
  T_RTMx       *lsmx, *smx, ProperSMx;
  T_RotMxInfo  *lrmxi, *rmxi;


  Centric = 0;
  InversionOffOrigin = 0;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    special = IsSpecialSeitzMx(SgInfo, lsmx, 0);

    if (special & SpecialSMx_Inversion)
    {
      if (special & SpecialSMx_Transl0)
        Centric = 1;
      else
        InversionOffOrigin = 1;

      break;
    }
  }

  if (InversionOffOrigin && Centric) {
    SetSgError("Internal Error: Corrupt SgInfo");
    return -1;
  }

  if (Centric == 0)
  {
    if (InversionOffOrigin) {
      SgInfo->Centric = 0;
      SgInfo->InversionOffOrigin = 1;
    }
    else
    {
      if (SgInfo->Centric) SgInfo->Centric = -1;
      SgInfo->InversionOffOrigin = 0;
    }
  }
  else
  {
    SgInfo->InversionOffOrigin = 0;

    lsmx  = SgInfo->ListSeitzMx;
    lrmxi = SgInfo->ListRotMxInfo;
    iList = 0;

    while (iList < SgInfo->nList)
    {
      deter = deterRotMx(lsmx->s.R);

      if (deter == -1 && SgInfo->Centric == -1)
      {
        for (i = 0; i < 9; i++)
              ProperSMx.s.R[i] = -lsmx->s.R[i];

        for (i = 0; i < 3; i++) {
              ProperSMx.s.T[i] = -lsmx->s.T[i] % STBF;
          if (ProperSMx.s.T[i] < 0)
              ProperSMx.s.T[i] += STBF;
        }

        smx = SgInfo->ListSeitzMx;

        for (i = 0; i < SgInfo->nList; i++, smx++)
          if (CompareSeitzMx(LI_P, &ProperSMx, smx) == 0) break;

        if (i == SgInfo->nList)
        {
          for (i = 0; i < 12; i++) lsmx->a[i] = ProperSMx.a[i];

          deter = deterRotMx(lsmx->s.R);

          if (deter != 1 || (lrmxi && GetRotMxInfo(lsmx->s.R, lrmxi) == 0)) {
            SetSgError(Err_Ill_SMx_in_List);
            return -1;
          }
        }
      }

      if      (deter == -1)
      {
        /* copy last element to this place */

            SgInfo->nList--;
        if (SgInfo->nList != iList)
        {
          smx = &SgInfo->ListSeitzMx[SgInfo->nList];
          for (i = 0; i < 12; i++) lsmx->a[i] = smx->a[i];

          if (lrmxi)
          {
            rmxi = &SgInfo->ListRotMxInfo[SgInfo->nList];
            CopyRotMxInfo(lrmxi, rmxi);
          }
        }
      }
      else if (deter ==  1)
      {
        lsmx++;
        if (lrmxi != NULL) lrmxi++;
        iList++;
      }
      else {
        SetSgError(Err_Ill_SMx_in_List);
        return -1;
      }
    }

    SgInfo->Centric = -1;
  }

  return 0;
}


int ApplyOriginShift(T_SgInfo *SgInfo)
{
  int     HaveOrSh, iList, i;
  T_RTMx  *lsmx, SMx;
  int     OrSh[3], lo[3];


  HaveOrSh = 0;

  for (i = 0; i < 3; i++) {
        OrSh[i] = SgInfo->OriginShift[i] * (STBF / 12);
    if (OrSh[i]) HaveOrSh = 1;
  }

  if (HaveOrSh == 0)
    return 0;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    RotMx_t_Vector(lo, lsmx->s.R, OrSh, STBF);

    for (i = 0; i < 3; i++)
      lsmx->s.T[i] = iModPositive(lsmx->s.T[i] - lo[i] + OrSh[i], STBF);
  }

  if (SgInfo->Centric == -1)
  {
    lsmx = &SMx;

    InitSeitzMx(lsmx, -1);

    RotMx_t_Vector(lo, lsmx->s.R, OrSh, STBF);

    for (i = 0; i < 3; i++)
      lsmx->s.T[i] = iModPositive(lsmx->s.T[i] - lo[i] + OrSh[i], STBF);

    if (CoreAdd2ListSeitzMx(SgInfo, lsmx) < 0)
      return -1;

    SgInfo->Centric = 0;
    SgInfo->InversionOffOrigin = 1;
  }

  return 1;
}


static void TidyTranslation(T_SgInfo *SgInfo)
{
  int        iList;
  int        iTrV, nTrV;
  const int  *TrV;
  T_RTMx     *lsmx;
  int        t1, t2, t3, *lt1, *lt2, *lt3, mint1, mint2, mint3;
  int        n0t, n0mint;


  nTrV = SgInfo->LatticeInfo->nTrVector;

  lsmx  = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    mint1 = *(lt1 = &lsmx->s.T[0]);
    mint2 = *(lt2 = &lsmx->s.T[1]);
    mint3 = *(lt3 = &lsmx->s.T[2]);

    TrV = SgInfo->LatticeInfo->TrVector;

    for (iTrV = 0; iTrV < nTrV; iTrV++)
    {
      t1 = ((*lt1) + *TrV++) % STBF;
      t2 = ((*lt2) + *TrV++) % STBF;
      t3 = ((*lt3) + *TrV++) % STBF;

      n0t = (t1 == 0) + (t2 == 0) + (t3 == 0);
      n0mint = (mint1 == 0) + (mint2 == 0) + (mint3 == 0);

      if (    n0t > n0mint
          || (   n0t == n0mint
              && (   t1 + t2 + t3 <  mint1 + mint2 + mint3
                  || (   t1 + t2 + t3 == mint1 + mint2 + mint3
                      && (t1 < mint1 || (t1 == mint1 && t2 < mint2)))))) {
        mint1 = t1; mint2 = t2; mint3 = t3;
      }
    }

    *lt1 = mint1;
    *lt2 = mint2;
    *lt3 = mint3;
  }
}


static T_SgInfo *Pt_SgInfo_ListSortFunction = NULL;

static int SgInfoListSortFunction(const int *iList_a, const int *iList_b)
{
  int          val_a, val_b, n0_a, n0_b, i;
  T_RTMx       *smx_a, *smx_b;
  T_RotMxInfo  RotMxInfo_a, RotMxInfo_b, *rmxi_a, *rmxi_b;
  T_SgInfo     *SgInfo = Pt_SgInfo_ListSortFunction;


  if (SgError != NULL) return 0;

  if (SgInfo->ListRotMxInfo == NULL)
  {
    rmxi_a = &RotMxInfo_a;
    rmxi_b = &RotMxInfo_b;

    smx_a = &SgInfo->ListSeitzMx[*iList_a];
    smx_b = &SgInfo->ListSeitzMx[*iList_b];

    if (   GetRotMxInfo(smx_a->s.R, rmxi_a) == 0
        || GetRotMxInfo(smx_b->s.R, rmxi_b) == 0)
    {
      SetSgError(Err_Ill_SMx_in_List);
      return 0;
    }
  }
  else
  {
    rmxi_a = &SgInfo->ListRotMxInfo[*iList_a];
    rmxi_b = &SgInfo->ListRotMxInfo[*iList_b];
  }

  val_a = abs(rmxi_a->Order);
  val_b = abs(rmxi_b->Order);

  if (val_a == 1 && val_b != 1) return -1;
  if (val_a != 1 && val_b == 1) return  1;
  if (rmxi_a->Order == 1 && rmxi_b->Order != 1) return -1;
  if (rmxi_a->Order != 1 && rmxi_b->Order == 1) return  1;

  if (val_a != 1)
  {
    if (val_a > val_b) return -1;
    if (val_a < val_b) return  1;
    if (rmxi_a->Order > rmxi_b->Order) return -1;
    if (rmxi_a->Order < rmxi_b->Order) return  1;
  }

  n0_a = n0_b = 0;

  for (i = 0; i < 3; i++)
  {
    if (rmxi_a->EigenVector[i] == 0) n0_a++;
    if (rmxi_b->EigenVector[i] == 0) n0_b++;
  }
  if (n0_a > n0_b) return -1;
  if (n0_a < n0_b) return  1;

  val_a = val_b = 0;

  for (i = 0; i < 3; i++)
  {
    if (val_a < abs(rmxi_a->EigenVector[i]))
        val_a = abs(rmxi_a->EigenVector[i]);
    if (val_b < abs(rmxi_b->EigenVector[i]))
        val_b = abs(rmxi_b->EigenVector[i]);
  }
  if (val_a < val_b) return -1;
  if (val_a > val_b) return  1;

  val_a =  100 * abs(rmxi_a->EigenVector[2]);
  val_a +=  10 * abs(rmxi_a->EigenVector[0]);
  val_a +=       abs(rmxi_a->EigenVector[1]);
  val_b =  100 * abs(rmxi_b->EigenVector[2]);
  val_b +=  10 * abs(rmxi_b->EigenVector[0]);
  val_b +=       abs(rmxi_b->EigenVector[1]);

  if (n0_a < 2)
  {
    if (val_a < val_b) return -1;
    if (val_a > val_b) return  1;
  }
  else
  {
    if (val_a > val_b) return -1;
    if (val_a < val_b) return  1;
  }

  for (i = 0; i < 3; i++)
  {
    if (rmxi_a->EigenVector[i] > rmxi_b->EigenVector[i]) return -1;
    if (rmxi_a->EigenVector[i] < rmxi_b->EigenVector[i]) return  1;
  }

  if (rmxi_a->Inverse < rmxi_b->Inverse) return -1;
  if (rmxi_a->Inverse > rmxi_b->Inverse) return  1;

  smx_a = &SgInfo->ListSeitzMx[*iList_a];
  smx_b = &SgInfo->ListSeitzMx[*iList_b];

  for (i = 0; i < 3; i++)
  {
    if (smx_a->s.T[i] < smx_b->s.T[i]) return -1;
    if (smx_a->s.T[i] > smx_b->s.T[i]) return  1;
  }

  return 0;
}


static int PostSortSgInfoList(const T_SgInfo *SgInfo, int *List_iList)
{
  int                nList, iL_iL, jL_iL;
  T_RTMx             BufMx1, BufMx2, *smxab;
  const T_RTMx       *lsmx, *smxb;
  T_RotMxInfo        RotMxInfo;
  const T_RotMxInfo  *rmxi;
  int                nO_1, iO, save, i, i_;


  nList = SgInfo->nList;

  iL_iL = 0;

  while (iL_iL < nList)
  {
    lsmx = &SgInfo->ListSeitzMx[List_iList[iL_iL]];

        rmxi = ListOrBufRotMxInfo(SgInfo, List_iList[iL_iL], &RotMxInfo);
    if (rmxi == NULL)
      return -1;

    iL_iL++;

    iO = rmxi->Order;
    if (iO < 0 && iO % 2) iO *= 2;
    nO_1 = abs(iO) - 1;

    smxab = &BufMx2;
    smxb = lsmx;

    for (iO = 1; iO < nO_1; iO++)
    {
      SeitzMxMultiply(smxab, lsmx, smxb);

      for (jL_iL = iL_iL; jL_iL < nList; jL_iL++)
      {
        smxb = &SgInfo->ListSeitzMx[List_iList[jL_iL]];
        if (CompareSeitzMx(SgInfo->LatticeInfo, smxab, smxb) == 0)
          break;
      }

      if (jL_iL < nList)
      {
                              save = List_iList[jL_iL];

        for (i = i_ = jL_iL; i > iL_iL; i = i_)
          List_iList[i] = List_iList[--i_];

        List_iList[iL_iL++] = save;
      }

      smxb = smxab;
      if (iO % 2) smxab = &BufMx1;
      else        smxab = &BufMx2;
    }
  }

  return 0;
}


static void SortSgInfoList(T_SgInfo *SgInfo, int *List_iList)
{
  int          i, j, refi;
  int          nList;
  T_RTMx       *lsmx;
  T_RotMxInfo  *lrmxi;


  if (SgError != NULL) return;

  nList = SgInfo->nList;
  lsmx  = SgInfo->ListSeitzMx;
  lrmxi = SgInfo->ListRotMxInfo;
  Pt_SgInfo_ListSortFunction = SgInfo;

  for (i = 0; i < nList; i++) List_iList[i] = i;

  qsort((void *) List_iList, nList, sizeof (*List_iList),
        (int (*)(const void *, const void *)) SgInfoListSortFunction);

  Pt_SgInfo_ListSortFunction = NULL;
  if (SgError != NULL) return;

  if (PostSortSgInfoList(SgInfo, List_iList) != 0)
    return;

  for (i = 0; i < nList; i++)
  {
    j = List_iList[i];

    if (j != i)
    {
      for (refi = i + 1; refi < nList; refi++)
        if (List_iList[refi] == i) break;

      if (refi >= nList) {
        SetSgError("Internal Error: SortSgInfoList(): Corrupt List_iList");
        return;
      }

      SwapRTMx(&lsmx[i], &lsmx[j]);
      if (lrmxi != NULL)
        SwapRotMxInfo(&lrmxi[i], &lrmxi[j]);

      List_iList[refi] = j;
    }
  }
}


int FindSeitzMx(const T_SgInfo *SgInfo,
                int Order, int HonorSign, int RefAxis, int DirCode)
{
  int          iList;
  int          MatchingOrder;
  T_RTMx       *lsmx;
  T_RotMxInfo  *lrmxi, RotMxInfo;


  lsmx  = SgInfo->ListSeitzMx;
  lrmxi = SgInfo->ListRotMxInfo;

  if (lrmxi == NULL) lrmxi = &RotMxInfo;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    if (lrmxi == &RotMxInfo)
    {
      if (GetRotMxInfo(lsmx->s.R, lrmxi) == 0) {
        SetSgError(Err_Ill_SMx_in_List);
        return -1;
      }
    }

    if (HonorSign == 0)
      MatchingOrder = (abs(Order) == abs(lrmxi->Order));
    else
      MatchingOrder = (Order == lrmxi->Order);

    if (   MatchingOrder
        && lrmxi->Inverse == 0
        && (RefAxis == 0 || RefAxis == lrmxi->RefAxis)
        && (DirCode == 0 || DirCode == lrmxi->DirCode))
    {
      if (DirCode != '*') return iList;

      if (   lrmxi->EigenVector[0] == 1
          && lrmxi->EigenVector[1] == 1
          && lrmxi->EigenVector[2] == 1) return iList;
    }

    if (lrmxi != &RotMxInfo) lrmxi++;
  }

  return -1;
}


static int FindXtalSystem(T_SgInfo *SgInfo)
{
  int                iList, i;
  int                HonorSign = 0, CheckEnantiomorph;
  const T_RTMx       *lsmx;
  T_RotMxInfo        RotMxInfo;
  const T_RotMxInfo  *lrmxi;

  enum Axis { N1 = 0, N2, N3, N4, N6, EndOfAxis };
  int                         N_count[EndOfAxis];


  SgInfo->XtalSystem    = XS_Unknown;
  SgInfo->UniqueRefAxis = 0;
  SgInfo->UniqueDirCode = 0;
  SgInfo->ExtraInfo     = EI_Unknown;

  CheckEnantiomorph = 0;

  for (i = 0; i < EndOfAxis; i++) N_count[i] = 0;

  lsmx  = SgInfo->ListSeitzMx;
  lrmxi = SgInfo->ListRotMxInfo;

  if (lrmxi == NULL) lrmxi = &RotMxInfo;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    if (lrmxi == &RotMxInfo)
    {
      if (GetRotMxInfo(lsmx->s.R, &RotMxInfo) == 0) {
        SetSgError(Err_Ill_SMx_in_List);
        return XS_Unknown;
      }
    }

    switch(abs(lrmxi->Order))
    {
      case 1:  i = N1; break;
      case 2:  i = N2; break;
      case 3:  i = N3; break;
      case 4:  i = N4; break;
      case 6:  i = N6; break;
      default:
        SetSgError("Internal Error: FindXtalSystem(): Corrupt ListRotMxInfo");
        return XS_Unknown;
    }

    if (lrmxi->Inverse == 0) /* skip N^-1 */
      N_count[i]++;

    if (lrmxi != &RotMxInfo) lrmxi++;
  }

  i = EndOfAxis;

  if (SgInfo->InversionOffOrigin == 1)
  {
    for (i = 0; i < EndOfAxis; i++)
    {
      if (N_count[i] % 2) break;
      N_count[i] /= 2;
    }
  }

  if (i == EndOfAxis)
  {
    if      (N_count[N3] == 4) SgInfo->XtalSystem = XS_Cubic;
    else if (N_count[N3] >  1) SgInfo->XtalSystem = XS_Unknown;
    else if (N_count[N6] == 1) SgInfo->XtalSystem = XS_Hexagonal;
    else if (N_count[N6] >  0) SgInfo->XtalSystem = XS_Unknown;
    else if (N_count[N3] == 1) SgInfo->XtalSystem = XS_Trigonal;
    else if (N_count[N4] == 1) SgInfo->XtalSystem = XS_Tetragonal;
    else if (N_count[N4] >  0) SgInfo->XtalSystem = XS_Unknown;
    else if (N_count[N2] >  2) SgInfo->XtalSystem = XS_Orthorhombic;
    else if (N_count[N2] >  0) SgInfo->XtalSystem = XS_Monoclinic;
    else if (N_count[N1] >  0) SgInfo->XtalSystem = XS_Triclinic;

    HonorSign = 1;
    iList = FindSeitzMx(SgInfo, -1, HonorSign, 'o', '.');
    if (iList < 0) HonorSign = 0;

    switch(SgInfo->XtalSystem)
    {
      case XS_Monoclinic:
        iList = FindSeitzMx(SgInfo, 2, HonorSign, 0, '=');
        if (iList < 0) SgInfo->XtalSystem = XS_Unknown;
        break;
      case XS_Tetragonal:
        CheckEnantiomorph = 1;
        iList = FindSeitzMx(SgInfo, 4, HonorSign, 0, '=');
        if (iList < 0) SgInfo->XtalSystem = XS_Unknown;
        break;
      case XS_Trigonal:
        CheckEnantiomorph = 1;
          iList = FindSeitzMx(SgInfo, 3, HonorSign, 0, '=');
        if (iList < 0)
          iList = FindSeitzMx(SgInfo, 3, HonorSign, 0, '*');
        if (iList < 0) SgInfo->XtalSystem = XS_Unknown;
        break;
      case XS_Hexagonal:
        CheckEnantiomorph = 1;
        iList = FindSeitzMx(SgInfo, 6, HonorSign, 0, '=');
        if (iList < 0) SgInfo->XtalSystem = XS_Unknown;
        break;
      case XS_Cubic:
        iList = FindSeitzMx(SgInfo, 4, HonorSign, 0, '=');
        if (iList >= 0) CheckEnantiomorph = 1;
        break;
      default:
        iList = -1;
        break;
    }
  }

  if (SgInfo->XtalSystem == XS_Unknown) {
    SetSgError("Error: Can't determine crystal system");
  }
  else if (iList >= 0)
  {
        lrmxi = ListOrBufRotMxInfo(SgInfo, iList, &RotMxInfo);
    if (lrmxi == NULL) {
      SgInfo->XtalSystem = XS_Unknown;
      return XS_Unknown;
    }

    if (SgInfo->XtalSystem != XS_Cubic)
    {
      SgInfo->UniqueRefAxis = lrmxi->RefAxis;
      SgInfo->UniqueDirCode = lrmxi->DirCode;

      if (SgInfo->XtalSystem == XS_Trigonal && lrmxi->DirCode == '=')
      {
        switch (lrmxi->RefAxis)
        {
          case 'z':
            switch (SgInfo->LatticeInfo->Code)
            {
              case 'R': SgInfo->ExtraInfo = EI_Obverse; break;
              case 'T': SgInfo->ExtraInfo = EI_Reverse; break;
            }
            break;
          case 'y':
            switch (SgInfo->LatticeInfo->Code)
            {
              case 'S': SgInfo->ExtraInfo = EI_Obverse; break;
              case 'R': SgInfo->ExtraInfo = EI_Reverse; break;
            }
            break;
          case 'x':
            switch (SgInfo->LatticeInfo->Code)
            {
              case 'T': SgInfo->ExtraInfo = EI_Obverse; break;
              case 'S': SgInfo->ExtraInfo = EI_Reverse; break;
            }
            break;
        }
      }
    }

    if (   HonorSign == 0 /* no inversion matrix */
        && SgInfo->LatticeInfo->Code == 'P'
        && CheckEnantiomorph == 1)
    {
      lsmx = &SgInfo->ListSeitzMx[iList];

      if (GetRotMxOrder(lsmx->s.R) > 1)
      {
        i =  lsmx->s.T[0] * lrmxi->EigenVector[0];
        i += lsmx->s.T[1] * lrmxi->EigenVector[1];
        i += lsmx->s.T[2] * lrmxi->EigenVector[2];

        if (i % (STBF / 2)) SgInfo->ExtraInfo = EI_Enantiomorphic;
      }
    }
  }

  return SgInfo->XtalSystem;
}


static int BuildGenerator_iList(T_SgInfo *SgInfo)
{
  int  iList, iList_1, nG;
  int  SgInfo_CI, HonorSign, Flag3asterisk, FlagPG;


  SgInfo_CI = (SgInfo->Centric || SgInfo->InversionOffOrigin);

  SgInfo->PointGroup = PG_Unknown;

        nG   = SgInfo->nGenerator = 0;
#define G_iL   SgInfo->Generator_iList

  HonorSign = 1;
  iList_1 = FindSeitzMx(SgInfo, -1, HonorSign, 'o', '.');
  if (iList_1 < 0) HonorSign = 0;

  if (SgInfo->XtalSystem == XS_Unknown)
    FindXtalSystem(SgInfo);

  switch(SgInfo->XtalSystem)
  {
    case XS_Triclinic:
      if (iList_1 < 0)
        iList_1 = FindSeitzMx(SgInfo, 1, HonorSign, 'o', '.');
      if (iList_1 >= 0) G_iL[nG++] = iList_1;

      if (SgInfo_CI)
        SgInfo->PointGroup = PG_1b;
      else
        SgInfo->PointGroup = PG_1;

      SgInfo->nGenerator = nG;
      return 0;

    case XS_Monoclinic:
      iList = FindSeitzMx(SgInfo, 2, HonorSign, 0, '=');
      if (iList < 0) break;
      G_iL[nG++] = iList;

      if (SgInfo_CI)
        SgInfo->PointGroup = PG_2_m;
      else if (deterRotMx(SgInfo->ListSeitzMx[iList].s.R) == -1)
        SgInfo->PointGroup = PG_m;
      else
        SgInfo->PointGroup = PG_2;

      if (iList_1 >= 0) G_iL[nG++] = iList_1;

      SgInfo->nGenerator = nG;
      return 0;

    case XS_Orthorhombic:
      iList = FindSeitzMx(SgInfo, 2, HonorSign, 'z', '=');
      if (iList >= 0) G_iL[nG++] = iList;

      iList = FindSeitzMx(SgInfo, 2, HonorSign, 'x', '=');
      if (iList >= 0) G_iL[nG++] = iList;

      if (nG < 2)
      {
        iList = FindSeitzMx(SgInfo, 2, HonorSign, 'y', '=');
        if (iList >= 0) G_iL[nG++] = iList;
      }

      if (nG != 2) break;

      if (SgInfo_CI)
        SgInfo->PointGroup = PG_mmm;
      else if (   deterRotMx(SgInfo->ListSeitzMx[G_iL[0]].s.R) == -1
               || deterRotMx(SgInfo->ListSeitzMx[G_iL[1]].s.R) == -1)
        SgInfo->PointGroup = PG_mm2;
      else
        SgInfo->PointGroup = PG_222;

      if (iList_1 >= 0) G_iL[nG++] = iList_1;

      SgInfo->nGenerator = nG;
      return 0;

    case XS_Tetragonal:
      iList = FindSeitzMx(SgInfo, 4, HonorSign, 0, '=');
      if (iList < 0) break;
      G_iL[nG++] = iList;

      if (          SgInfo->UniqueRefAxis != 'x')
      {
        iList = FindSeitzMx(SgInfo, 2, HonorSign, 'x', '=');
        if (iList >= 0) G_iL[nG++] = iList;
      }
      if (nG < 2 && SgInfo->UniqueRefAxis != 'z')
      {
        iList = FindSeitzMx(SgInfo, 2, HonorSign, 'z', '=');
        if (iList >= 0) G_iL[nG++] = iList;
      }
      if (nG < 2 && SgInfo->UniqueRefAxis != 'y')
      {
        iList = FindSeitzMx(SgInfo, 2, HonorSign, 'y', '=');
        if (iList >= 0) G_iL[nG++] = iList;
      }

      if (nG < 2)
      {
        if (SgInfo_CI)
          SgInfo->PointGroup = PG_4_m;
        else if (deterRotMx(SgInfo->ListSeitzMx[G_iL[0]].s.R) == -1)
          SgInfo->PointGroup = PG_4b;
        else
          SgInfo->PointGroup = PG_4;
      }
      else
      {
        if (SgInfo_CI)
          SgInfo->PointGroup = PG_4_mmm;
        else if (deterRotMx(SgInfo->ListSeitzMx[G_iL[0]].s.R) == -1)
        {
          if (deterRotMx(SgInfo->ListSeitzMx[G_iL[1]].s.R) == -1)
            SgInfo->PointGroup = PG_4bm2;
          else
            SgInfo->PointGroup = PG_4b2m;
        }
        else
        {
          if (deterRotMx(SgInfo->ListSeitzMx[G_iL[1]].s.R) == -1)
            SgInfo->PointGroup = PG_4mm;
          else
            SgInfo->PointGroup = PG_422;
        }
      }

      if (iList_1 >= 0) G_iL[nG++] = iList_1;

      SgInfo->nGenerator = nG;
      return 0;

    case XS_Trigonal:
    case XS_Hexagonal:
      Flag3asterisk = 0;

      if (SgInfo->XtalSystem == XS_Trigonal)
      {
          iList = FindSeitzMx(SgInfo, 3, HonorSign, 0, '=');
        if (iList < 0)
        {
          iList = FindSeitzMx(SgInfo, 3, HonorSign, 0, '*');
          Flag3asterisk = 1;
        }
      }
      else
        iList = FindSeitzMx(SgInfo, 6, HonorSign, 0, '=');

      if (iList < 0) break;
      G_iL[nG++] = iList;

      iList = FindSeitzMx(SgInfo, 2, HonorSign, 0, '\'');
      if (iList >= 0) G_iL[nG++] = iList;

      FlagPG = -1;

      if (nG < 2)
      {
        iList = FindSeitzMx(SgInfo, 2, HonorSign, 0, '"');
        if (iList >= 0) G_iL[nG++] = iList;
        FlagPG = 1;
      }

      if (SgInfo->XtalSystem == XS_Trigonal)
      {
        if (nG < 2)
        {
          if (SgInfo_CI) SgInfo->PointGroup = PG_3b;
          else           SgInfo->PointGroup = PG_3;
        }
        else
        {
          if (Flag3asterisk == 1)
          {
            if (SgInfo_CI)
              SgInfo->PointGroup = PG_3bm;
            else
            {     FlagPG = deterRotMx(SgInfo->ListSeitzMx[G_iL[1]].s.R);
              if (FlagPG == -1) SgInfo->PointGroup = PG_3m;
              else              SgInfo->PointGroup = PG_32;
            }
          }
          else if (FlagPG == -1)
          {
            if (SgInfo_CI)
              SgInfo->PointGroup = PG_3b1m;
            else
            {     FlagPG = deterRotMx(SgInfo->ListSeitzMx[G_iL[1]].s.R);
              if (FlagPG == -1) SgInfo->PointGroup = PG_31m;
              else              SgInfo->PointGroup = PG_312;
            }
          }
          else
          {
            if (SgInfo_CI)
              SgInfo->PointGroup = PG_3bm1;
            else
            {     FlagPG = deterRotMx(SgInfo->ListSeitzMx[G_iL[1]].s.R);
              if (FlagPG == -1) SgInfo->PointGroup = PG_3m1;
              else              SgInfo->PointGroup = PG_321;
            }
          }
        }
      }
      else
      {
        if (nG < 2)
        {
          if (SgInfo_CI)
            SgInfo->PointGroup = PG_6_m;
          else if (deterRotMx(SgInfo->ListSeitzMx[G_iL[0]].s.R) == -1)
            SgInfo->PointGroup = PG_6b;
          else
            SgInfo->PointGroup = PG_6;
        }
        else
        {
          if (SgInfo_CI)
            SgInfo->PointGroup = PG_6_mmm;
          else if (deterRotMx(SgInfo->ListSeitzMx[G_iL[0]].s.R) == -1)
          {
            if (deterRotMx(SgInfo->ListSeitzMx[G_iL[1]].s.R) == FlagPG)
              SgInfo->PointGroup = PG_6b2m;
            else
              SgInfo->PointGroup = PG_6bm2;
          }
          else if (deterRotMx(SgInfo->ListSeitzMx[G_iL[1]].s.R) == -1)
            SgInfo->PointGroup = PG_6mm;
          else
            SgInfo->PointGroup = PG_622;
        }
      }

      if (iList_1 >= 0) G_iL[nG++] = iList_1;

      SgInfo->nGenerator = nG;
      return 0;

    case XS_Cubic:
      FlagPG = 0;

        iList = FindSeitzMx(SgInfo, 4, HonorSign, 'z', '=');
      if (iList < 0) {
        iList = FindSeitzMx(SgInfo, 2, HonorSign, 'z', '=');
        FlagPG = 1;
      }
      if (iList < 0) break;
      G_iL[nG++] = iList;

      iList = FindSeitzMx(SgInfo, 2, HonorSign, 'x', '=');
      if (iList < 0) break;
      G_iL[nG++] = iList;

      iList = FindSeitzMx(SgInfo, 3, HonorSign, 'o', '*');
      if (iList < 0) break;
      G_iL[nG++] = iList;

      if (FlagPG)
      {
        if (SgInfo_CI) SgInfo->PointGroup = PG_m3b;
        else           SgInfo->PointGroup = PG_23;
      }
      else
      {
        if (SgInfo_CI)
          SgInfo->PointGroup = PG_m3bm;
        else if (deterRotMx(SgInfo->ListSeitzMx[G_iL[0]].s.R) == -1)
          SgInfo->PointGroup = PG_4b3m;
        else
          SgInfo->PointGroup = PG_432;
      }

      if (iList_1 >= 0) G_iL[nG++] = iList_1;

      SgInfo->nGenerator = nG;
      return 0;

    default:
      break;
  }

#undef G_iL

  return -1;
}


static int BuildHSym(T_SgInfo *SgInfo)
{
  int                NeedDash, HaveOrSh, nGRT, iList, iG, ip, os, i;
  int                AbsOrder, RefAxis, DirCode, ScrewPrimitive, Screw;
  int                PreviousRotation, PreviousRefAxis, PreviousDirCode;
  int                nTrV, iTrV, OrSh[3], RO[3], Transl[3];
  const int          *TrV, *ht;
  T_RTMx             SMx_1;
  const T_RTMx       *SMx;
  const T_RotMxInfo  *rmxi;
  char               *hsym, *hsym_mark;

  struct
  {
    //TK reordered fields (to avoid clang-tidy report about padding)
    const T_RotMxInfo  *RMxI;
    int                Transl[3];
    T_RotMxInfo         RMxI_Buf;
  }
  GRT[sizeof SgInfo->Generator_iList / sizeof (*SgInfo->Generator_iList) + 1];

  const char *Digits = "0123456";


  if (SgInfo->nGenerator == 0) {
    SetSgError("Internal Error: BuildHSym(): Empty generator list");
    return -1;
  }

  HaveOrSh = 0;

  for (i = 0; i < 3; i++) {
        OrSh[i] = SgInfo->OriginShift[i] * (STBF / 12);
    if (OrSh[i]) HaveOrSh = 1;
  }

  NeedDash = 0;
  nGRT = 0;

  for (iG = 0; iG < SgInfo->nGenerator; iG++)
  {
    iList = SgInfo->Generator_iList[iG];

        GRT[nGRT].RMxI = ListOrBufRotMxInfo(SgInfo, iList, &GRT[nGRT].RMxI_Buf);
    if (GRT[nGRT].RMxI == NULL)
      return -1;

    SMx = &SgInfo->ListSeitzMx[iList];

    RotMx_t_Vector(RO, SMx->s.R, OrSh, STBF);

    for (i = 0; i < 3; i++)
      GRT[nGRT].Transl[i] = iModPositive(SMx->s.T[i] + RO[i] - OrSh[i], STBF);

    if (GRT[nGRT].RMxI->Order == -1)
    {
      for (i = 0; i < 3; i++)
        if (GRT[nGRT].Transl[i] != 0) break;

      if (i == 3) NeedDash = 1;
      else        nGRT++;
    }
    else
      nGRT++;
  }

  if (SgInfo->Centric)
  {
    if (HaveOrSh == 0)
      NeedDash = 1;
    else
    {
      for (iG = 0; iG < nGRT; iG++)
        if (GRT[iG].RMxI->Order == 1) break;

      InitSeitzMx(&SMx_1, -1);

      if (GetRotMxInfo(SMx_1.s.R, &GRT[iG].RMxI_Buf) != -1) {
        SetSgError("Internal Error: BuildHSym(): Corrupt GetRotMxInfo()");
        return -1;
      }

      GRT[iG].RMxI = &GRT[iG].RMxI_Buf;

      for (i = 0; i < 3; i++)
        GRT[iG].Transl[i] = iModPositive(-2 * OrSh[i], STBF);

      if (iG == nGRT)
        nGRT++;
    }
  }

  hsym = SgInfo->HallSymbol;

  for (i = 0; i <= MaxLenHallSymbol; i++)
    *hsym++ = '\0';

  PreviousRotation = 0;
  PreviousRefAxis = 0;
  PreviousDirCode = 0;

  hsym = SgInfo->HallSymbol;

  if (NeedDash)
    *hsym++ = '-';
  else
    *hsym++ = ' ';

  *hsym++ = SgInfo->LatticeInfo->Code;

  nTrV = SgInfo->LatticeInfo->nTrVector;

  for (iG = 0; iG < nGRT; iG++)
  {
    rmxi = GRT[iG].RMxI;

    AbsOrder = abs(rmxi->Order);
    RefAxis = rmxi->RefAxis;
    DirCode = rmxi->DirCode;
    if (RefAxis == 'o') RefAxis = 0;
    if (DirCode == '=' || DirCode == '.') DirCode = 0;

    if (iG == 0)
    {
      if (RefAxis == 'z') RefAxis = 0;
    }
    else
    {
      if      (AbsOrder == 2)
      {
        if      (PreviousRotation == 2 || PreviousRotation == 4)
        {
          if (RefAxis == 'x') RefAxis = 0;
        }
        else if (PreviousRotation == 3 || PreviousRotation == 6)
        {
          if (   PreviousDirCode == '*'
              || RefAxis == PreviousRefAxis) RefAxis = 0;
          if (DirCode == '\'') DirCode = 0;
        }
      }
      else if (AbsOrder == 3)
      {
        if (DirCode == '*') DirCode = 0;
      }
    }

    PreviousRotation = AbsOrder;
    PreviousRefAxis = rmxi->RefAxis;
    PreviousDirCode = rmxi->DirCode;

    *hsym++ = ' ';
    if (rmxi->Order < 0) *hsym++ = '-';
    *hsym++ = Digits[AbsOrder];
    if (RefAxis) *hsym++ = RefAxis;
    if (DirCode) *hsym++ = DirCode;

    TrV = SgInfo->LatticeInfo->TrVector;

    for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
    {
      for (i = 0; i < 3; i++)
        if ((GRT[iG].Transl[i] + TrV[i]) % STBF != 0)
          break;

      if (i == 3)
        break;
    }

    if (iTrV < nTrV)
      continue; /* next iG */

    hsym_mark = hsym;

    TrV = SgInfo->LatticeInfo->TrVector;

    for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3, hsym = hsym_mark)
    {
      for (i = 0; i < 3; i++)
        Transl[i] = iModPositive(GRT[iG].Transl[i] + TrV[i], STBF);

      /* Next line commented by NCrystal developers to avoid static analysis error (value assigned to Screw never used): */
      /*Screw = 0;*/

      for (ip = 0; ip < 3; ip++)
        if (rmxi->EigenVector[ip] != 0) break;

      if (ip < 3 && rmxi->EigenVector[ip] == 1)
      {
        for (i = ip + 1; i < 3; i++)
          if (rmxi->EigenVector[i] != 0) break;

        if (i == 3)
        {
          ScrewPrimitive = STBF / AbsOrder;
          Screw = Transl[ip] / ScrewPrimitive;
              i = Screw * ScrewPrimitive;
          if (i % 3)
          {
            *hsym++ = Digits[Screw];
            Transl[ip] -= i;
          }
        }
      }

      ht = HallTranslations;

      while (*ht)
      {
        for (i = 0; i < 3; i++)
          if (Transl[i] < ht[i + 1]) break;

        if (i == 3)
        {
          for (i = 0; i < 3; i++)
            Transl[i] -= ht[i + 1];

          *hsym++ = (char) *ht;
        }

        ht += 4;
      }

      for (i = 0; i < 3; i++)
        if (Transl[i] != 0)
          break;

      if (i == 3)
        break;
    }

    if (iTrV == nTrV)
      return 0;
  }

  if (nGRT == 0)
  {
    *hsym++ = ' ';
    *hsym++ = '1';
  }

  if (HaveOrSh)
  {
    *hsym++ = ' ';
    *hsym++ = '(';

    for (i = 0; i < 3; i++)
    {
      if (i) *hsym++ = ' ';

          os = iModPositive(SgInfo->OriginShift[i], 12);
      if (os > 6)
      {
        *hsym++ = '-';
        os = 12 - os;
      }

      *hsym++ = Digits[os];
    }

    *hsym++ = ')';
  }

  *hsym = '\0';

  if (SgInfo->HallSymbol[MaxLenHallSymbol] != '\0') {
    SetSgError("Internal Error: BuildHSym(): MaxLenHallSymbol too small");
    return -1;
  }

  return 1;
}


static int BuildHallSymbol(T_SgInfo *SgInfo, int FixedOriginShift)
{
  int     ix, iy, iz;
  int     status;

  static const int ShiftTable[] = { 0, 1, -1, 2, -2, 3 };


  if (SgError != NULL) return -1;

  if (SgInfo->nGenerator == 0)
  {
    if (BuildGenerator_iList(SgInfo) != 0)
    {
      SetSgError("Error: Can't build generator list");
      return -1;
    }
  }

  if (FixedOriginShift)
  {
        status = BuildHSym(SgInfo);
    if (status == 1)
      return 0;
  }
  else
  {
    for (ix = 0; ix < 6; ix++)
    {
      SgInfo->OriginShift[0] = ShiftTable[ix];

      for (iy = 0; iy < 6; iy++)
      {
        SgInfo->OriginShift[1] = ShiftTable[iy];

        for (iz = 0; iz < 6; iz++)
        {
          SgInfo->OriginShift[2] = ShiftTable[iz];

              status = BuildHSym(SgInfo);
          if (status < 0)
            return -1;

          if (status == 1)
            return 0;
        }
      }
    }
  }

  SetSgError("Error: Can't build Hall Symbol");
  return -1;
}


void InitSgInfo(T_SgInfo *SgInfo)
{
  int  i;


  SgInfo->GenOption = 0;
  SgInfo->Centric = 0;
  SgInfo->InversionOffOrigin = 0;
  SgInfo->LatticeInfo = LI_P;
  SgInfo->StatusLatticeTr = 0;
  for (i = 0; i < 3; i++)
    SgInfo->OriginShift[i] = 0;
  SgInfo->nList = 0;

  SgInfo->OrderL = 0;
  SgInfo->OrderP = 0;
  SgInfo->XtalSystem = XS_Unknown;
  SgInfo->UniqueRefAxis = 0;
  SgInfo->UniqueDirCode = 0;
  SgInfo->ExtraInfo = EI_Unknown;
  SgInfo->PointGroup = PG_Unknown;
  SgInfo->nGenerator = 0;
  SgInfo->HallSymbol[0] = '\0';
  SgInfo->TabSgName = NULL;
  SgInfo->CCMx_LP = NULL;
  SgInfo->n_si_Vector = -1;
}


int CompleteSgInfo(T_SgInfo *SgInfo)
{
  int                List_iList[192];
  const T_TabSgName  *tsgn;


  if (SgInfo->StatusLatticeTr == -1)
  {
    if (AddLatticeTr2ListSeitzMx(SgInfo, SgInfo->LatticeInfo) < 0)
      return -1;
  }

  if (ApplyOriginShift(SgInfo) < 0)
    return -1;

  if (SgInfo->nList >  ((int)sizeof List_iList) / ((int)sizeof(*List_iList))) {
    SetSgError("Internal Error: CompleteSgInfo()");
    return -1;
  }

  if (SgInfo->nList > 1)
  {
    SortSgInfoList(SgInfo, List_iList);
    if (SgError != NULL) return -1;
  }

  if (RemoveLatticeTr(SgInfo) != 0)
    return -1;

  if (RemoveInversion(SgInfo) != 0)
    return -1;

  TidyTranslation(SgInfo);

  if (SgInfo->nList > 1)
  {
    SortSgInfoList(SgInfo, List_iList);
    if (SgError != NULL) return -1;
  }
                             SgInfo->OrderP = SgInfo->nList;
  if (SgInfo->Centric == -1) SgInfo->OrderP *= 2;

  SgInfo->OrderL = SgInfo->OrderP * SgInfo->LatticeInfo->nTrVector;

  if (BuildHallSymbol(SgInfo, 0) != 0)
    return -1;

  for (tsgn = TabSgName; tsgn->HallSymbol; tsgn++)
    if (   strcmp(tsgn->HallSymbol, SgInfo->HallSymbol) == 0
        && (   SgInfo->TabSgName == NULL
            || SgInfo->TabSgName == tsgn))
      break;

  if (SgInfo->TabSgName != NULL && tsgn->HallSymbol == NULL)
  {
    if (SgError) return -1;

    sprintf(SgErrorBuffer,
      "Internal Error: Input/Output HallSymbol mismatch: %s <> %s",
      SgInfo->TabSgName->HallSymbol, SgInfo->HallSymbol);

    SetSgError(SgErrorBuffer);
    return -1;
  }

  if (tsgn->HallSymbol)
    SgInfo->TabSgName = tsgn;

  SgInfo->CCMx_LP = NULL;

  switch (SgInfo->LatticeInfo->Code)
  {
    case 'P': SgInfo->CCMx_LP = CCMx_PP; break;
    case 'A': SgInfo->CCMx_LP = CCMx_AP; break;
    case 'B': SgInfo->CCMx_LP = CCMx_BP; break;
    case 'C': SgInfo->CCMx_LP = CCMx_CP; break;
    case 'I': SgInfo->CCMx_LP = CCMx_IP; break;
    case 'R':
      switch (SgInfo->UniqueRefAxis) {
        case   0:
        case 'z': SgInfo->CCMx_LP = CCMx_RP_z; break;
        case 'y': SgInfo->CCMx_LP = CCMx_RP_y; break;
        default: break;
      }
      break;
    case 'S':
      switch (SgInfo->UniqueRefAxis) {
        case   0:
        case 'y': SgInfo->CCMx_LP = CCMx_SP_y; break;
        case 'x': SgInfo->CCMx_LP = CCMx_SP_x; break;
        default: break;
      }
      break;
    case 'T':
      switch (SgInfo->UniqueRefAxis) {
        case   0:
        case 'x': SgInfo->CCMx_LP = CCMx_TP_x; break;
        case 'z': SgInfo->CCMx_LP = CCMx_TP_z; break;
        default: break;
      }
      break;
    case 'F': SgInfo->CCMx_LP = CCMx_FP; break;
    default:
      break;
  }

  if (SgInfo->CCMx_LP == NULL) {
    SetSgError("Internal Error: Illegal lattice code");
    return -1;
  }

  return 0;
}


int CB_SMx(T_RTMx *CSiC,
           const T_RTMx *CBMx, const T_RTMx *SMx, const T_RTMx *InvCBMx)
{
  int     i;
  T_RTMx  BufMx;


  RTMxMultiply(&BufMx, SMx,  InvCBMx, CTBF / STBF, CTBF);
  RTMxMultiply(CSiC,   CBMx, &BufMx,  CRBF,        CRBF * CTBF);

  for (i = 0; i < 9; i++)
  {
    if (CSiC->s.R[i] % (CRBF * CRBF)) {
      SetSgError("Internal Error: Corrupt CBMx/SMx/InvCBMx");
      return -1;
    }

    CSiC->s.R[i] /= (CRBF * CRBF);
  }

  for (i = 0; i < 3; i++)
  {
    if (CSiC->s.T[i] % (CRBF * (CTBF / STBF))) {
      SetSgError("Internal Error: Out of STBF range");
      return -1;
    }

    CSiC->s.T[i] /= (CRBF * (CTBF / STBF));
  }

  return 0;
}


int TransformSgInfo(const T_SgInfo *SgInfo,
                    const T_RTMx *CBMx, const T_RTMx *InvCBMx,
                    T_SgInfo *BC_SgInfo)
{
  int           iList, f, i;
  int           nTrV, iTrV, nLoopInv, iLoopInv;
  const int     *TrV;
  T_RTMx        SMx, BC_SMx;
  const T_RTMx  *lsmx;


  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
      {
        for (i = 0; i < 9; i++)
          SMx.s.R[i] = f * lsmx->s.R[i];

        for (i = 0; i < 3; i++)
          SMx.s.T[i] = f * lsmx->s.T[i] + TrV[i];

        if (CB_SMx(&BC_SMx, CBMx, &SMx, InvCBMx) != 0)
          return -1;

        if (Add2ListSeitzMx(BC_SgInfo, &BC_SMx) < 0)
          return -1;
      }
    }
  }

  return 0;
}


#undef SGCLIB_C__
/*
  Space Group Info's (c) 1994-96 Ralf W. Grosse-Kunstleve
 */



#define SGCOREDEF__


static int InitialCBMxR(T_SgInfo *SgInfo,
                        const T_LatticeInfo **NewLatticeInfo,
                        int *NewPointGroup,
                        int *IniCBMxR, int *IniInvCBMxR)
{
  int                  Code, NewPG, deterCCMx, fac, i;
  const T_LatticeInfo  *NewLI;
  const int            *CCMx;


  Code  = SgInfo->LatticeInfo->Code;
  NewLI = SgInfo->LatticeInfo;
  NewPG = SgInfo->PointGroup;
  CCMx  = CCMx_PP;

  switch(SgInfo->XtalSystem)
  {
    case XS_Triclinic:
      NewLI = LI_P;
      CCMx  = SgInfo->CCMx_LP;
      break;

    case XS_Monoclinic:
    case XS_Tetragonal:
      switch (SgInfo->UniqueRefAxis)
      {
        case 'z': if      (Code == 'C') {
                    NewLI = LI_P; CCMx = SgInfo->CCMx_LP; }
                  else if (Code == 'F') {
                    NewLI = LI_I; CCMx = CCMx_FI_z; }
                  break;
        case 'y': if      (Code == 'B') {
                    NewLI = LI_P; CCMx = SgInfo->CCMx_LP; }
                  else if (Code == 'F') {
                    NewLI = LI_I; CCMx = CCMx_FI_y; }
                  break;
        case 'x': if      (Code == 'A') {
                    NewLI = LI_P; CCMx = SgInfo->CCMx_LP; }
                  else if (Code == 'F') {
                    NewLI = LI_I; CCMx = CCMx_FI_x; }
                  break;
        default:
          goto ReturnError;
      }

      if (   SgInfo->XtalSystem == XS_Tetragonal
          && SgInfo->LatticeInfo != NewLI)
      {
        if      (NewPG == PG_4b2m) NewPG = PG_4bm2;
        else if (NewPG == PG_4bm2) NewPG = PG_4b2m;
      }

      break;

    case XS_Orthorhombic:
      break;

    case XS_Trigonal:
      NewLI = LI_P;
      CCMx  = SgInfo->CCMx_LP;

      if (Code == 'R' || Code == 'S' || Code == 'T')
      {
        if      (NewPG == PG_321)  NewPG = PG_32;
        else if (NewPG == PG_3m1)  NewPG = PG_3m;
        else if (NewPG == PG_3bm1) NewPG = PG_3bm;
      }

      break;

    case XS_Hexagonal:
      break;

    case XS_Cubic:
      break;

    default:
      goto ReturnError;
  }

      deterCCMx = deterRotMx(CCMx);
  if (deterCCMx < 1 || CRBF % deterCCMx)
    goto ReturnError;

  fac = CRBF / deterCCMx;

  InverseRotMx(CCMx, IniInvCBMxR);

  for (i = 0; i < 9; i++) {
       IniCBMxR[i] = CRBF * CCMx[i];
    IniInvCBMxR[i] *= fac;
  }

  *NewLatticeInfo = NewLI;
  *NewPointGroup  = NewPG;

  return deterCCMx;

  ReturnError:

  SetSgError("Internal Error: InitialCBMxR()");
  return -1;
}


static int CBR_RMx(int *RRMx,
                   const int *CBMxR, const int *RMx, const int *InvCBMxR)
{
  int  i, BufMx[9];


  RotMxMultiply(BufMx, RMx,   InvCBMxR);
  RotMxMultiply(RRMx,  CBMxR, BufMx);

  for (i = 0; i < 9; i++)
  {
    if (RRMx[i] % (CRBF * CRBF)) {
      SetSgError("Internal Error: CBR_SMx()");
      return -1;
    }

    RRMx[i] /= (CRBF * CRBF);
  }

  return 0;
}


static void RotateCBMxR(const int *RMx, const int *InvRMx,
                        int *CBMxR, int *InvCBMxR)
{
  int  i, BufMx[9];


  RotMxMultiply(BufMx, RMx, CBMxR);
  for (i = 0; i < 9; i++)    CBMxR[i] = BufMx[i];

  /* matrix algebra: (A * B)^-1 = B^-1 * A^-1 */

  RotMxMultiply(BufMx, InvCBMxR, InvRMx);
  for (i = 0; i < 9; i++) InvCBMxR[i] = BufMx[i];
}


static int AlignUniqueAxis(const T_SgInfo *SgInfo,
                           const T_SgInfo *GenSgI,
                           int *CBMxR, int *InvCBMxR,
                           const int **AlignRMx)
{
  int          i, iListS, DirCode;
  int          UAMx[9], RotEV[3];
  const int    *RMx, *InvRMx, *lsmxR;
  T_RotMxInfo  RMxI_S, *RMxI_G;


  if (GenSgI->nList < 2)
    goto ReturnError;

  RMxI_G = &GenSgI->ListRotMxInfo[1];

  if (abs(RMxI_G->Order) == 3) DirCode = 0;
  else                         DirCode = '=';

      iListS =   FindSeitzMx(SgInfo,  RMxI_G->Order, 1, 0, DirCode);
  if (iListS < 0)
  {
    if (SgInfo->Centric == 0)
      return 0;

        iListS = FindSeitzMx(SgInfo, -RMxI_G->Order, 1, 0, DirCode);
    if (iListS < 0)
      goto ReturnError;

    for (i = 0; i < 9; i++)
      UAMx[i] = -SgInfo->ListSeitzMx[iListS].s.R[i];

    lsmxR = UAMx;
  }
  else
    lsmxR = SgInfo->ListSeitzMx[iListS].s.R;

  if (CBR_RMx(UAMx, CBMxR, lsmxR, InvCBMxR) != 0)
    goto ReturnError;

  if (GetRotMxInfo(UAMx, &RMxI_S) != RMxI_G->Order)
    goto ReturnError;

  if (RMxI_S.DirCode != RMxI_G->DirCode)
    return 0;

  RMx = InvRMx = RMx_1_000;

  for (;;)
  {
    RotMx_t_Vector(RotEV, RMx, RMxI_S.EigenVector, 0);

    for (i = 0; i < 3; i++)
      if (RotEV[i] != RMxI_G->EigenVector[i]) break;

    if (i == 3)
      break;

    if      (RMxI_S.DirCode == '=')
    {
      if      (RMx == RMx_1_000) {
               RMx =  RMx_3_111; InvRMx = RMx_3i111; }
      else if (RMx == RMx_3_111) {
               RMx =  RMx_3i111; InvRMx = RMx_3_111; }
      else
        goto ReturnError;
    }
    else if (RMxI_S.DirCode == '*')
    {
      if      (RMx == RMx_1_000) {
               RMx =  RMx_4_001; InvRMx = RMx_4i001; }
      else if (RMx == RMx_4_001) {
               RMx =  RMx_2_001; InvRMx = RMx_2_001; }
      else if (RMx == RMx_2_001) {
               RMx =  RMx_4i001; InvRMx = RMx_4_001; }
      else
        goto ReturnError;
    }
    else
      goto ReturnError;
  }

  if (RMx != RMx_1_000)
    RotateCBMxR(RMx, InvRMx, CBMxR, InvCBMxR);

  if (AlignRMx)
     *AlignRMx = RMx;

  return 1;

  ReturnError:

  SetSgError("Internal Error: AlignUniqueAxis()");
  return -1;
}


static const T_RTMx *GetSMxWithSameRot(const int *WtdRotMx,
                                       const T_SgInfo *SgInfo, T_RTMx *BufMx)
{
  int           iList, i;
  const T_RTMx  *lsmx;


  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    for (i = 0; i < 9; i++)
      if (WtdRotMx[i] !=  lsmx->s.R[i])
        break;

    if (i == 9)
      return lsmx;

    if (SgInfo->Centric != -1)
      continue;

    for (i = 0; i < 9; i++)
      if (WtdRotMx[i] != -lsmx->s.R[i])
        break;

    if (i == 9)
    {
      for (i = 0; i < 9; i++)
            BufMx->s.R[i] = -lsmx->s.R[i];

      for (i = 0; i < 3; i++) {
            BufMx->s.T[i] = -lsmx->s.T[i] % STBF;
        if (BufMx->s.T[i] < 0)
            BufMx->s.T[i] += STBF;
      }

      return BufMx;
    }
  }

  return NULL;
}


static int BuildFreeMx(const int *EigenVector, int Order, int *FreeMx)
{
  static const int GeneratorEigenVectors[] =
    {
       001,   0,  0,  1,
       010,   0,  1,  0,
       100,   1,  0,  0,
       110,   1,  1,  0,
      -110,   1, -1,  0,
       111,   1,  1,  1,
         0
    };

  int        i;
  const int  *gev;


  for (i = 0; i < 9; i++)
    FreeMx[i] = 0;

  if (Order == -1 || Order == -3 || Order == -4 || Order == -6)
    return 0;

  for (gev = GeneratorEigenVectors; *gev++ != 0; gev += 3)
  {
    for (i = 0; i < 3; i++)
      if (EigenVector[i] != gev[i])
        break;

    if (i == 3)
      break;
  }

  gev--;

  if      (Order == -2)
  {
    switch (*gev)
    {
      case  001: FreeMx[0] =  1; FreeMx[4] =  1;                 return 0;
      case  010: FreeMx[8] =  1; FreeMx[0] =  1;                 return 0;
      case  100: FreeMx[4] =  1; FreeMx[8] =  1;                 return 0;
      case  110: FreeMx[1] =  1; FreeMx[4] = -1; FreeMx[8] =  1; return 1;
      case -110: FreeMx[1] =  1; FreeMx[4] =  1; FreeMx[8] =  1; return 1;
      default:
        break;
    }
  }
  else if (Order > 1)
  {
    switch (*gev)
    {
      case  001: FreeMx[8] =  1;                                 return 0;
      case  010: FreeMx[4] =  1;                                 return 0;
      case  100: FreeMx[0] =  1;                                 return 0;
      case  110: FreeMx[0] =  1; FreeMx[3] =  1;                 return 1;
      case -110: FreeMx[0] =  1; FreeMx[3] = -1;                 return 1;
      case  111: FreeMx[0] =  1; FreeMx[3] =  1; FreeMx[6] =  1; return 1;
      default:
        break;
    }
  }

  SetSgError("Internal Error: BuildFreeMx()");
  return -1;
}


static int StartFixAxes(const T_SgInfo *SgInfo,
                        const T_SgInfo *GenSgI, const int *iGen,
                        T_RTMx *CBMx, T_RTMx *InvCBMx,
                        T_RTMx *SMxG, T_RTMx *SMxS_G,
                        int *FreeMx, int *TryAgain)
{
  int                iG, Order, i;
  const int          *EV;
  T_RTMx             SMxG_S, BufMx;
  const T_RTMx       *SMx;
  const T_RotMxInfo  *RMxI_G;


  if (*iGen == -3)
    iG = 1;
  else
    iG = *iGen;

  if (iG == -1)
  {
    Order = -1;
    EV    = NULL;
  }
  else
  {
    if (iG < 1 || iG >= GenSgI->nList)
      goto ReturnError;

            RMxI_G = &GenSgI->ListRotMxInfo[iG];
    Order = RMxI_G->Order;
    EV    = RMxI_G->EigenVector;

    if (iG != *iGen)
    {
          Order *= -1;
      if (Order != *iGen)
        goto ReturnError;
    }
  }

  if (Order == -1)
  {
    if (GenSgI->Centric == -1) {
      InitSeitzMx(SMxG, -1);
    }
    else
    {
      for (iG = 1; iG < GenSgI->nList; iG++)
        if (GenSgI->ListRotMxInfo[iG].Order == -1)
          break;

      if (iG == GenSgI->nList)
        goto ReturnError;

      SMx = &GenSgI->ListSeitzMx[iG];

      for (i = 0; i < 12; i++) SMxG->a[i] = SMx->a[i];
    }
  }
  else
  {
    SMx = &GenSgI->ListSeitzMx[iG];

    if (iG == *iGen)
      for (i = 0; i < 12; i++) SMxG->a[i] = SMx->a[i];
    else
    {
      for (i = 0; i < 9; i++)
            SMxG->s.R[i] = -SMx->s.R[i];

      for (i = 0; i < 3; i++) {
            SMxG->s.T[i] = -SMx->s.T[i] % STBF;
        if (SMxG->s.T[i] < 0)
            SMxG->s.T[i] += STBF;
      }
    }
  }

  if (CB_SMx(&SMxG_S, InvCBMx, SMxG, CBMx) != 0)
    return -1;

      SMx = GetSMxWithSameRot(SMxG_S.s.R, SgInfo, &BufMx);
  if (SMx == NULL)
    return 0;

  if (CB_SMx(SMxS_G, CBMx, SMx, InvCBMx) != 0)
    return -1;

  for (i = 0; i < 9; i++)
    if (SMxS_G->s.R[i] != SMxG->s.R[i])
      goto ReturnError;

  *TryAgain = BuildFreeMx(EV, Order, FreeMx);/*indentation of this line fixed by NCrystal developers*/
  if (*TryAgain < 0)
    return -1;

  return 1;

  ReturnError:

  SetSgError("Internal Error: StartFixAxes()");
  return -1;
}


static int FindInvertableMx(const int *Mx, int *InvMx,
                            int *nActive, int *irActive, int *icActive)
{
  int  Init, deterMx, i;


  if (*nActive == 0 || *nActive == 3)
    return 0;

  if (*nActive == -1)
  {
    Init = 1;

        deterMx = deterRotMx(Mx);
    if (deterMx)
    {
      InverseRotMx(Mx, InvMx);

      *nActive = 3;
      return deterMx;
    }
  }
  else
    Init = 0;

  if (Init || *nActive == 2)
  {
    for (;;)
    {
      if (Init)
      {
        irActive[0] = 0;
        irActive[1] = 1;
        icActive[0] = 0;
        icActive[1] = 1;
        Init = 0;
      }
      else
      {
        if (++icActive[1] == 3) {
          if (++icActive[0] == 2) {
            if (++irActive[1] == 3) {
              if (++irActive[0] == 2) {
                Init = 1;
                break;
              }
              else {
                irActive[1] = irActive[0] + 1;
                icActive[0] = 0;
                icActive[1] = 1;
              }
            }
            else {
              icActive[0] = 0;
              icActive[1] = 1;
            }
          }
          else {
            icActive[1] = icActive[0] + 1;
          }
        }
      }

      InvMx[0] =   Mx[irActive[1] * 3 + icActive[1]];
      InvMx[1] = - Mx[irActive[0] * 3 + icActive[1]];
      InvMx[2] = - Mx[irActive[1] * 3 + icActive[0]];
      InvMx[3] =   Mx[irActive[0] * 3 + icActive[0]];

          deterMx = InvMx[3] * InvMx[0] - InvMx[1] * InvMx[2];
      if (deterMx) {
        *nActive = 2;
        return deterMx;
      }
    }
  }

  if (*nActive == 2)
    return 0;

  if (Init) i = 0;
  else      i = irActive[0] * 3 + icActive[0] + 1;

  for ( ; i < 9; i++)
  {
    if (Mx[i]) {
      irActive[0] = i / 3;
      icActive[0] = i % 3;
      *nActive = 1;
      return Mx[i];
    }
  }

  if (*nActive == 1)
    return 0;

  *nActive = 0;
  return 1;
}


static int SetInvCBMxT(const int *CBMxT, const int *InvCBMxR, int *InvCBMxT)
{
  int  i;


  RotMx_t_Vector(InvCBMxT, InvCBMxR, CBMxT, CRBF * CTBF);

  for (i = 0; i < 3; i++)
  {
    if (InvCBMxT[i] % CRBF) {
      SetSgError("Internal Error: SetInvCBMxT()");
      return -1;
    }

    if (InvCBMxT[i])
        InvCBMxT[i] = CTBF - InvCBMxT[i] / CRBF;
  }

  return 0;
}


static int FixAxes(const T_SgInfo *SgInfo,
                   const T_SgInfo *GenSgI, const int *iGen,
                   T_RTMx *CBMx, T_RTMx *InvCBMx,
                   int *FreeMx, int TryAgain)
{
  int        i, NextTryAgain;
  int        IniCBMxT[3], SingleFreeMx[9];
  T_RTMx     SMxG, SMxS_G;
  int        NextFreeMxBuf[9], R_I_FMxBuf[9];
  int        R_I[9], *R_I_FMx, InvR_I_FMx[9], deterR_I_FMx;
  int        S_G[3], CmpS_G[3], RedSh[3], Sh[3], *NextFreeMx;
  int        nActive, irActive[3], icActive[3];
  int        nTrV, iTrV;
  const int  *TrV;

  icActive[0]=icActive[1]=icActive[2]=0;
  irActive[0]=irActive[1]=irActive[2]=0;

  if (FreeMx == NULL) {
    for (i = 0; i < 3; i++) {
         CBMx->s.T[i] = 0;
      InvCBMx->s.T[i] = 0;
    }
  }

  i = StartFixAxes(SgInfo, GenSgI, iGen, CBMx, InvCBMx,
                   &SMxG, &SMxS_G, SingleFreeMx, &NextTryAgain);
  if (i != 1)
    return i;

  if (FreeMx) {
    RotMxMultiply(NextFreeMxBuf, SingleFreeMx, FreeMx);
    NextFreeMx =  NextFreeMxBuf;
  }
  else
    NextFreeMx = SingleFreeMx;

  for (i = 0; i < 9; i++)
    R_I[i] = SMxG.s.R[i];

  for (i = 0; i < 9; i += 4)
    R_I[i] -= 1;

  if (FreeMx) {
    RotMxMultiply(R_I_FMxBuf, R_I, FreeMx);
    R_I_FMx =     R_I_FMxBuf;
  }
  else
    R_I_FMx = R_I;

  for (i = 0; i < 3; i++)
    IniCBMxT[i] = CBMx->s.T[i];

  nActive = -1;

  for (;;)
  {
    deterR_I_FMx = FindInvertableMx(R_I_FMx, InvR_I_FMx,
                                    &nActive, irActive, icActive);
    if (deterR_I_FMx == 0)
      break;

    nTrV = GenSgI->LatticeInfo->nTrVector;
     TrV = GenSgI->LatticeInfo->TrVector;

    for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
    {
      for (i = 0; i < 3; i++) {
        S_G[i] =   (CTBF / STBF)
                 * ((SMxS_G.s.T[i] - SMxG.s.T[i] - TrV[i]) % STBF);
        RedSh[i] = 0;
      }

      switch(nActive)
      {
        case 1:
          RedSh[icActive[0]] = S_G[irActive[0]];
          break;
        case 2:
          RedSh[icActive[0]] =   InvR_I_FMx[0] * S_G[irActive[0]]
                               + InvR_I_FMx[1] * S_G[irActive[1]];
          RedSh[icActive[1]] =   InvR_I_FMx[2] * S_G[irActive[0]]
                               + InvR_I_FMx[3] * S_G[irActive[1]];
          break;
        case 3:
          RotMx_t_Vector(RedSh, InvR_I_FMx, S_G, 0);
          break;
        default:
          break;
      }

      if (FreeMx)
      {
        RotMx_t_Vector(Sh, FreeMx, RedSh, 0);

        for (i = 0; i < 3; i++)
          Sh[i] %= (CTBF * abs(deterR_I_FMx));
      }
      else
      {
        for (i = 0; i < 3; i++)
          Sh[i] = RedSh[i] % (CTBF * abs(deterR_I_FMx));
      }

      RotMx_t_Vector(CmpS_G, R_I, Sh, 0);

      for (i = 0; i < 3; i++)
        if ((CmpS_G[i] - S_G[i] * deterR_I_FMx) % (CTBF * abs(deterR_I_FMx)))
          break;

      if (i < 3)
        continue;

      if (deterR_I_FMx != 1)
      {
        for (i = 0; i < 3; i++)
        {
          if (Sh[i] % abs(deterR_I_FMx))
            goto ReturnError;

          Sh[i] /= deterR_I_FMx;
        }
      }

      for (i = 0; i < 3; i++) {
            CBMx->s.T[i] = IniCBMxT[i] + Sh[i] % CTBF;
        if (CBMx->s.T[i] < 0)
            CBMx->s.T[i] += CTBF;
      }

      if (SetInvCBMxT(CBMx->s.T, InvCBMx->s.R, InvCBMx->s.T) != 0)
        return -1;

      if (iGen[1] == 0)
        return 1;

          i = FixAxes(SgInfo, GenSgI, &iGen[1], CBMx, InvCBMx,
                      NextFreeMx, NextTryAgain);
      if (i != 0)
        return i;
    }

    if (TryAgain == 0)
      break;
  }

  return 0;

  ReturnError:

  SetSgError("Internal Error: FixAxes()");
  return -1;
}


static int CompleteCBMx(const T_SgInfo *SgInfo, const T_LatticeInfo *NewLI,
                        const T_SgInfo *GenSgI,
                        const int *IniCBMxR, const int *IniInvCBMxR,
                        T_RTMx       *CBMx,  T_RTMx       *InvCBMx)
{
  int  iGen[5], i;


  if (SgInfo->XtalSystem == XS_Triclinic)
  {
    for (i = 0; i < 9; i++) {
         CBMx->s.R[i] =    IniCBMxR[i];
      InvCBMx->s.R[i] = IniInvCBMxR[i];
    }

    if (GenSgI->PointGroup == PG_1)
    {
      for (i = 0; i < 3; i++) {
           CBMx->s.T[i] = 0;
        InvCBMx->s.T[i] = 0;
      }
      return 1;
    }

    iGen[0] = -1;
    iGen[1] =  0;

    return FixAxes(SgInfo, GenSgI, iGen, CBMx, InvCBMx, NULL, 0);
  }

  if (SgInfo->XtalSystem == XS_Monoclinic)
  {
    int        iCCs, BufRMx[9];
    int        RMxCCs_Buf[9], RMxCCn_Buf[9], InvRMxCCn_Buf[9], RotLTrV[3];
    const int  *RMxAA, *RMxCCs, *RMxCCn, *InvRMxCCn, *TrV;
    T_RTMx     BufCBMx, BufInvCBMx;


    if (NewLI->nTrVector != 1 && NewLI->nTrVector != 2)
      goto ReturnError;

    for (i = 0; i < 9; i++) {
         BufCBMx.s.R[i] =    IniCBMxR[i];
      BufInvCBMx.s.R[i] = IniInvCBMxR[i];
    }

        i = AlignUniqueAxis(SgInfo, GenSgI,
                            BufCBMx.s.R, BufInvCBMx.s.R, &RMxAA);
    if (i != 1)
      return i;

    if (GenSgI->nList < 2)
      goto ReturnError;

    for (i = 0; i < 9; i++) {
      RMxCCs_Buf[i] = RMx_2_110[i];
      RMxCCn_Buf[i] = RMx_3_001[i];
    }

    switch (GenSgI->ListRotMxInfo[1].RefAxis)
    {
      case 'z': break;
      case 'x': RotateRotMx(RMxCCs_Buf, RMx_3_111, RMx_3i111);
                RotateRotMx(RMxCCn_Buf, RMx_3_111, RMx_3i111);
                break;
      case 'y': RotateRotMx(RMxCCs_Buf, RMx_3i111, RMx_3_111);
                RotateRotMx(RMxCCn_Buf, RMx_3i111, RMx_3_111);
                break;
      default:
        goto ReturnError;
    }

    InverseRotMx(RMxCCn_Buf, InvRMxCCn_Buf);

                                           i = 0;
                                      iGen[i++] =  1;
    if (GenSgI->PointGroup == PG_2_m) iGen[i++] = -1;
    iGen[i  ] =  0;/*indentation of this line fixed by NCrystal developers*/

    RMxCCs = RMx_1_000;

    for (iCCs = 0; iCCs < 2; iCCs++, RMxCCs = RMxCCs_Buf)
    {
      RMxCCn = InvRMxCCn = RMx_1_000;

      for (;;)
      {
        if (NewLI->nTrVector == 2)
        {
          RotMx_t_Vector(RotLTrV, RMxAA,  &NewLI->TrVector[3], STBF);
          RotMx_t_Vector(BufRMx,  RMxCCn, RotLTrV,             STBF);
          RotMx_t_Vector(RotLTrV, RMxCCs, BufRMx,              STBF);

          TrV = &GenSgI->LatticeInfo->TrVector[3];

          for (i = 0; i < 3; i++)
            if (RotLTrV[i] != TrV[i])
              break;
        }

        if (NewLI->nTrVector == 1 || i == 3)
        {
          RotMxMultiply(BufRMx,    RMxCCn, BufCBMx.s.R);
          RotMxMultiply(CBMx->s.R, RMxCCs, BufRMx);

          RotMxMultiply(BufRMx,    BufInvCBMx.s.R, InvRMxCCn);
          RotMxMultiply(InvCBMx->s.R, BufRMx,         RMxCCs);

              i = FixAxes(SgInfo, GenSgI, iGen, CBMx, InvCBMx, NULL, 0);
          if (i != 0)
            return i;
        }

        if      (RMxCCn ==    RMx_1_000) {
                 RMxCCn =     RMxCCn_Buf; InvRMxCCn = InvRMxCCn_Buf; }
        else if (RMxCCn ==    RMxCCn_Buf) {
                 RMxCCn =  InvRMxCCn_Buf; InvRMxCCn =    RMxCCn_Buf; }
        else {
                 RMxCCn = NULL;
                 break;
        }
      }
    }

    return 0;
  }

  for (i = 0; i < 9; i++) {
       CBMx->s.R[i] =    IniCBMxR[i];
    InvCBMx->s.R[i] = IniInvCBMxR[i];
  }

  if (SgInfo->XtalSystem == XS_Orthorhombic)
  {
    int        iNextBasis;
    int        BufCBMxR[9], BufInvCBMxR[9];
    int         NLTrV_Buf1[3], NLTrV_Buf2[3];
    const int  *NLTrV, *GLTrV;


    if ((GenSgI->LatticeInfo->Code == 'I') != (NewLI->Code == 'I'))
      return 0;

    if (   NewLI->Code == 'A'
        || NewLI->Code == 'B'
        || NewLI->Code == 'C') {
      NLTrV =               &NewLI->TrVector[3];
      GLTrV = &GenSgI->LatticeInfo->TrVector[3]; }
    else {
      NLTrV = NULL;
      GLTrV = NULL; }
                                           i = 0;
                                      iGen[i++] =  1;
                                      iGen[i++] =  2;
    if (GenSgI->PointGroup == PG_mmm) iGen[i++] = -1;
    iGen[i  ] =  0;/*indentation of this line fixed by NCrystal developers*/

    for (iNextBasis = 0; iNextBasis < 6; iNextBasis++)
    {
      if      (iNextBasis % 2)
      {
        RotMxMultiply(   BufCBMxR, RMx_2_110,    CBMx->s.R);
        RotMxMultiply(BufInvCBMxR, InvCBMx->s.R, RMx_2_110);

        for (i = 0; i < 9; i++) {
             CBMx->s.R[i] =    BufCBMxR[i];
          InvCBMx->s.R[i] = BufInvCBMxR[i];
        }
      }
      else if (iNextBasis == 2) {
        RotMxMultiply(   CBMx->s.R, RMx_3_111,   IniCBMxR);
        RotMxMultiply(InvCBMx->s.R, IniInvCBMxR, RMx_3i111);
      }
      else if (iNextBasis) {
        RotMxMultiply(   CBMx->s.R, RMx_3i111,   IniCBMxR);
        RotMxMultiply(InvCBMx->s.R, IniInvCBMxR, RMx_3_111);
      }

      if (NLTrV)
      {
        if      (iNextBasis % 2) {
          RotMx_t_Vector(NLTrV_Buf1, RMx_2_110, NLTrV, STBF);
          NLTrV = NLTrV_Buf1;
        }
        else if (iNextBasis == 2) {
          RotMx_t_Vector(NLTrV_Buf2, RMx_3_111, &NewLI->TrVector[3], STBF);
          NLTrV = NLTrV_Buf2;
        }
        else if (iNextBasis) {
          RotMx_t_Vector(NLTrV_Buf2, RMx_3i111, &NewLI->TrVector[3], STBF);
          NLTrV = NLTrV_Buf2;
        }

        for (i = 0; i < 3; i++)
          if (NLTrV[i] != GLTrV[i])
            break;

        if (i < 3)
          continue;
      }

          i = FixAxes(SgInfo, GenSgI, iGen, CBMx, InvCBMx, NULL, 0);
      if (i != 0)
        return i;
    }
  }

  if (SgInfo->XtalSystem == XS_Tetragonal)
  {
        i = AlignUniqueAxis(SgInfo, GenSgI, CBMx->s.R, InvCBMx->s.R, NULL);
    if (i != 1)
      return -1;
                          i = 0;
                     iGen[i++] =  1;

    switch (GenSgI->PointGroup) {
      case PG_422:
      case PG_4mm:
      case PG_4b2m:
      case PG_4bm2:
      case PG_4_mmm: iGen[i++] =  2;
    }

    switch (GenSgI->PointGroup) {
      case PG_4_m:
      case PG_4_mmm: iGen[i++] = -1;
    }
                     iGen[i  ] =  0;

    return FixAxes(SgInfo, GenSgI, iGen, CBMx, InvCBMx, NULL, 0);
  }

  if (SgInfo->XtalSystem == XS_Trigonal)
  {
        i = AlignUniqueAxis(SgInfo, GenSgI, CBMx->s.R, InvCBMx->s.R, NULL);
    if (i != 1)
      return i;
                         i = 0;

    switch (GenSgI->PointGroup) {
      case PG_3:
      case PG_312:
      case PG_32:
      case PG_3m1:
      case PG_3m:   iGen[i++] =  1;
                    break;
      case PG_3b:
      case PG_3bm1:
      case PG_3b1m:
      case PG_3bm:  iGen[i++] = -3;
    }

    switch (GenSgI->PointGroup) {
      case PG_321:
      case PG_312:
      case PG_32:
      case PG_3m1:
      case PG_31m:
      case PG_3m:
      case PG_3bm1:
      case PG_3b1m:
      case PG_3bm:  iGen[i++] =  2;
    }

    switch (GenSgI->PointGroup) {
      case PG_321:
      case PG_31m:  iGen[i++] =  1;
    }
                    iGen[i  ] =  0;

    return FixAxes(SgInfo, GenSgI, iGen, CBMx, InvCBMx, NULL, 0);
  }

  if (SgInfo->XtalSystem == XS_Hexagonal)
  {
        i = AlignUniqueAxis(SgInfo, GenSgI, CBMx->s.R, InvCBMx->s.R, NULL);
    if (i != 1)
      return -1;
                         i = 0;

    switch (GenSgI->PointGroup) {
      case PG_6bm2:
      case PG_6b2m: iGen[i++] =  2;
    }
                    iGen[i++] =  1;

    switch (GenSgI->PointGroup) {
      case PG_622:
      case PG_6mm:
      case PG_6_mmm: iGen[i++] =  2;
    }

    switch (GenSgI->PointGroup) {
      case PG_6_m:
      case PG_6_mmm: iGen[i++] = -1;
    }
                     iGen[i  ] =  0;

    return FixAxes(SgInfo, GenSgI, iGen, CBMx, InvCBMx, NULL, 0);
  }

  if (SgInfo->XtalSystem == XS_Cubic)
  {                                            i = 0;
                                          iGen[i++] =  3;
                                          iGen[i++] =  1;
                                          iGen[i++] =  2;
    if (   GenSgI->PointGroup == PG_m3b
        || GenSgI->PointGroup == PG_m3bm) iGen[i++] = -1;
    iGen[i  ] =  0;/*indentation of this line fixed by NCrystal developers*/

    return FixAxes(SgInfo, GenSgI, iGen, CBMx, InvCBMx, NULL, 0);
  }

  return 0;

  ReturnError:

  SetSgError("Internal Error: CompleteCBMx()");
  return -1;
}


const T_TabSgName *FindReferenceSpaceGroup(T_SgInfo *SgInfo,
                                           T_RTMx *CBMx, T_RTMx *InvCBMx)
{
  int                  stat, NewPG, SgInfo_CI, OL_SgInfo, OL_GenSgI;
  const T_TabSgName    *tsgn;
  T_SgInfo             GenSgI;
  T_RTMx               GenSgI_ListSeitzMx[5];
  T_RotMxInfo          GenSgI_ListRotMxInfo[5];
  int                  iList, PrevSgNumber;
  int                  FacIniCBMxR;
  T_RotMxInfo          *lrmxi;
  const T_LatticeInfo  *NewLI;
  int                  IniCBMxR[9], IniInvCBMxR[9];


  GenSgI.MaxList       = 5;
  GenSgI.ListSeitzMx   = GenSgI_ListSeitzMx;
  GenSgI.ListRotMxInfo = GenSgI_ListRotMxInfo;

      FacIniCBMxR = InitialCBMxR(SgInfo, &NewLI, &NewPG, IniCBMxR, IniInvCBMxR);
  if (FacIniCBMxR < 0)
    return NULL;

      OL_SgInfo = SgInfo->OrderL;
  if (OL_SgInfo % FacIniCBMxR)
    goto ReturnError;

  OL_SgInfo /= FacIniCBMxR;

  SgInfo_CI = (SgInfo->Centric || SgInfo->InversionOffOrigin);

  PrevSgNumber = 0;

  for (tsgn = TabSgName; tsgn->HallSymbol; tsgn++)
  {
    if (tsgn->HallSymbol[1] == 'R')
      continue;

    if (VolAPointGroups[tsgn->SgNumber] != NewPG)
      continue;

    if (tsgn->SgNumber == PrevSgNumber)
      continue;

    PrevSgNumber = tsgn->SgNumber;

    InitSgInfo(&GenSgI);
    GenSgI.GenOption = -1;

    ParseHallSymbol(tsgn->HallSymbol, &GenSgI);

    if (SgError != NULL)
      return NULL;

    if (ApplyOriginShift(&GenSgI) < 0)
      return NULL;

    if (SgInfo_CI != (GenSgI.Centric || GenSgI.InversionOffOrigin))
      goto ReturnError;

    OL_GenSgI = GenSgI.LatticeInfo->nTrVector;

    if (SgInfo_CI)
      OL_GenSgI *= 2;

    lrmxi = &GenSgI.ListRotMxInfo[1];

    for (iList = 1; iList < GenSgI.nList; iList++, lrmxi++)
    {
      OL_GenSgI *= abs(lrmxi->Order);

      if (   (lrmxi->Order == -1 || lrmxi->Order == -3)
          && GenSgI.Centric == 0 && GenSgI.InversionOffOrigin == 0)
        goto ReturnError;
    }

    if (OL_GenSgI == OL_SgInfo)
    {
      if (NewLI->nTrVector != GenSgI.LatticeInfo->nTrVector)
        goto ReturnError;

      GenSgI.PointGroup = NewPG;

#if DEBUG_FindConventionalSetting
      fprintf(stdout, "%s ?= %s (%d)\n",
        SgInfo->HallSymbol, tsgn->HallSymbol, tsgn->SgNumber);
#endif

      stat = CompleteCBMx(SgInfo, NewLI, &GenSgI,
                          IniCBMxR, IniInvCBMxR,
                             CBMx,     InvCBMx);
      if (stat < 0)
        return NULL;

      if (stat)
        return tsgn;
    }
  }

  SetSgError("Internal Error: Space Group not found");
  return NULL;

  ReturnError:

  SetSgError("Internal Error: FindReferenceSpaceGroup()");
  return NULL;
}
/*
  Space Group Info's (c) 1994-96 Ralf W. Grosse-Kunstleve
 */





static const char *IErr_Inc_SymMx =
  "Internal Error: Inconsistent symmetry matrices";


int IsSysAbsent_hkl(const T_SgInfo *SgInfo,
                    int h, int k, int l, int *TH_Restriction)
{
  int           iTrV, nTrV;
  const int     *TrV;
  int           iList, mh, mk, ml, hm, km, lm;
  int           TH, THr, FlagMismatch;
  const T_RTMx  *lsmx;


  mh = -h;
  mk = -k;
  ml = -l;

  /* check list of symmetry operations
     take care of lattice type and "centric" flag */

  THr = -1;
  if (TH_Restriction != NULL) *TH_Restriction = THr;
  FlagMismatch = 0;

  nTrV = SgInfo->LatticeInfo->nTrVector;
  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    hm = lsmx->s.R[0] * h + lsmx->s.R[3] * k + lsmx->s.R[6] * l;
    km = lsmx->s.R[1] * h + lsmx->s.R[4] * k + lsmx->s.R[7] * l;
    lm = lsmx->s.R[2] * h + lsmx->s.R[5] * k + lsmx->s.R[8] * l;

    TrV = SgInfo->LatticeInfo->TrVector;

    for (iTrV = 0; iTrV < nTrV; iTrV++)
    {
      TH =  (lsmx->s.T[0] + *TrV++) * h;
      TH += (lsmx->s.T[1] + *TrV++) * k;
      TH += (lsmx->s.T[2] + *TrV++) * l;
      TH %= STBF; if (TH < 0) TH += STBF;

      if      (mh == hm && mk == km && ml == lm)
      {
        if (TH != 0 && SgInfo->Centric == -1)
          return -(iList + 1 + iTrV * SgInfo->nList);

        if (THr < 0)
          THr = TH;
        else if (THr != TH)
          FlagMismatch = 1; /* must be systematic absent */
                            /* will check later ...      */
      }
      else if ( h == hm &&  k == km &&  l == lm)
      {
        if (TH != 0)
          return  (iList + 1 + iTrV * SgInfo->nList);
      }
      else
        break;
    }
  }

  if (THr >= 0 && FlagMismatch) /* ... consistency check */
    SetSgError(IErr_Inc_SymMx);

  if (TH_Restriction != NULL)
  {
    if (SgInfo->Centric == -1) *TH_Restriction = 0;
    else                       *TH_Restriction = THr;
  }

  return 0;
}


int BuildEq_hkl(const T_SgInfo *SgInfo, T_Eq_hkl *Eq_hkl, int h, int k, int l)
{
  int       iList, hm, km, lm, i;
  T_RTMx    *lsmx;
  T_Eq_hkl  BufEq_hkl;


  if (Eq_hkl == NULL)
    Eq_hkl = &BufEq_hkl;

  Eq_hkl->M = 1;
  Eq_hkl->N = 1;
  Eq_hkl->h[0] = h;
  Eq_hkl->k[0] = k;
  Eq_hkl->l[0] = l;
  Eq_hkl->TH[0] = 0;

  if (! (h || k || l))
    return Eq_hkl->M; /* this is 000 */

  Eq_hkl->M++;

  /* check list of symmetry operations */

  lsmx = &SgInfo->ListSeitzMx[1]; /* skip first = identity matrix */

  for (iList = 1; iList < SgInfo->nList; iList++, lsmx++)
  {
    hm = lsmx->s.R[0] * h + lsmx->s.R[3] * k + lsmx->s.R[6] * l;
    km = lsmx->s.R[1] * h + lsmx->s.R[4] * k + lsmx->s.R[7] * l;
    lm = lsmx->s.R[2] * h + lsmx->s.R[5] * k + lsmx->s.R[8] * l;

    for (i = 0; i < Eq_hkl->N; i++)
    {
      if ( ( hm == Eq_hkl->h[i] &&  km == Eq_hkl->k[i] &&  lm == Eq_hkl->l[i])
        || (-hm == Eq_hkl->h[i] && -km == Eq_hkl->k[i] && -lm == Eq_hkl->l[i]))
        break;
    }

    if (i == Eq_hkl->N)
    {
      if (Eq_hkl->N >= 24) {
        SetSgError(IErr_Inc_SymMx);
        return 0;
      }

      Eq_hkl->h[i] = hm;
      Eq_hkl->k[i] = km;
      Eq_hkl->l[i] = lm;

          Eq_hkl->TH[i] = (  lsmx->s.T[0] * h
                           + lsmx->s.T[1] * k
                           + lsmx->s.T[2] * l) % STBF;
      if (Eq_hkl->TH[i] < 0)
          Eq_hkl->TH[i] += STBF;

      Eq_hkl->M += 2;
      Eq_hkl->N++;
    }
  }

  if (SgInfo->nList % Eq_hkl->N) /* another error trap */ {
    SetSgError(IErr_Inc_SymMx);
    return 0;
  }

  return Eq_hkl->M;
}


int AreSymEquivalent_hkl(const T_SgInfo *SgInfo, int h1, int k1, int l1,
                                                 int h2, int k2, int l2)
{
  int     iList, mh2, mk2, ml2, hm, km, lm;
  T_RTMx  *lsmx;


  mh2 = -h2;
  mk2 = -k2;
  ml2 = -l2;

  /* check list of symmetry operations */

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    hm = lsmx->s.R[0] * h1 + lsmx->s.R[3] * k1 + lsmx->s.R[6] * l1;
    km = lsmx->s.R[1] * h1 + lsmx->s.R[4] * k1 + lsmx->s.R[7] * l1;
    lm = lsmx->s.R[2] * h1 + lsmx->s.R[5] * k1 + lsmx->s.R[8] * l1;

    if      ( h2 == hm &&  k2 == km &&  l2 == lm)
      return  (iList + 1);

    else if (mh2 == hm && mk2 == km && ml2 == lm)
      return -(iList + 1);
  }

  return 0;
}


void SetListMin_hkl(const T_SgInfo *SgInfo,            int  Maxk, int  Maxl,
                                            int *Minh, int *Mink, int *Minl)
{
  *Minh = 0;

  switch(SgInfo->XtalSystem)
  {
    case XS_Triclinic:
      *Mink = -Maxk;
      *Minl = -Maxl;
      break;
    case XS_Monoclinic:
      if (SgInfo->UniqueRefAxis == 'z')
      {
        *Mink = -Maxk;
        *Minl = 0;
      }
      else
      {
        *Mink = 0;
        *Minl = -Maxl;
      }
      break;
    default:
      if (SgInfo->XtalSystem == XS_Trigonal && SgInfo->UniqueDirCode == '*')
        *Mink = -Maxk;
      else
        *Mink = 0;
      *Minl = 0;
      break;
  }
}


int IsSuppressed_hkl(const T_SgInfo *SgInfo, int Minh, int Mink, int Minl,
                                                       int Maxk, int Maxl,
                                             int    h, int    k, int    l)
{
  int     iList, mate, hm, km, lm;
  T_RTMx  *lsmx;


  /* check for Friedel mate first */

  hm = -h, km = -k, lm = -l;

  if (   (Minh <= hm && hm <=    h)
      && (Mink <= km && km <= Maxk)
      && (Minl <= lm && lm <= Maxl))
  {
    if (hm < h) return -1;
    else /* if (h == 0) */
      if (km < k) return -1;
      else if (k == 0)
        if (lm < l) return -1;
  }
  lsmx = &SgInfo->ListSeitzMx[1]; /* skip first = identity matrix */

  for (iList = 1; iList < SgInfo->nList; iList++, lsmx++)
  {
    /* check if equivalent hm, km, lm are inside loop range ... */

    hm = lsmx->s.R[0] * h + lsmx->s.R[3] * k + lsmx->s.R[6] * l;
    km = lsmx->s.R[1] * h + lsmx->s.R[4] * k + lsmx->s.R[7] * l;
    lm = lsmx->s.R[2] * h + lsmx->s.R[5] * k + lsmx->s.R[8] * l;

    for (mate = 0; mate < 2; mate++)
    {
      if (mate) hm = -hm, km = -km, lm = -lm; /* ... or friedel mate */

      if (   Minh <= hm && hm <=    h
          && Mink <= km && km <= Maxk
          && Minl <= lm && lm <= Maxl)
      {
        /* ... and were processed before */

        if (hm < h)
          return (mate ? -(iList + 1) : iList + 1);
        else /* if (hm == h) */
          if (km < k)
            return (mate ? -(iList + 1) : iList + 1);
          else if (km == k)
            if (lm < l)
              return (mate ? -(iList + 1) : iList + 1);
      }
    }
  }

  return 0;
}
/*
  Space Group Info's (c) 1994-96 Ralf W. Grosse-Kunstleve
 */



#define SGCOREDEF__


typedef struct
  {
    int         OriginChoice;
    int         CellChoice;
    int         BasisChoice;
    const char  *BT_or_UA;
  }
  T_ExtInfo;


static const char *Ext_BT_or_UA[] =
  {
    /*  0 */ "abc",
    /*  1 */ "ba-c",
    /*  2 */ "cab",
    /*  3 */ "-cba",
    /*  4 */ "bca",
    /*  5 */ "a-cb",
    /*  6 */         "bac", /* 6 -> 1 */
    /*  7 */         "cba", /* 7 -> 3 */
    /*  8 */         "acb", /* 8 -> 5 */
    /*  9 */ "-b", "b-", "bb", "bb", /* 10, 11, 12 ->  9 */
    /* 13 */ "-c", "c-", "bc", "cb", /* 14, 15, 16 -> 13 */
    /* 17 */ "-a", "a-", "ba", "ab", /* 18, 19, 20 -> 17 */
    /* 21 */ "b",
    /* 22 */ "c",
    /* 23 */ "a",
    NULL
  };


typedef struct
    {
      int     Improper, Rotation, RefAxis, DirCode, Screw;
      T_RTMx  SeitzMx;
    }
    T_HallGenerator;


#define SkipWhite(cp) while (*(cp) && (*(cp) == '_' || isspace(*(cp)))) (cp)++


static const char *IErr_Corrupt_TabSgName =
   "Internal Error: Corrupt TabSgName";


static int FindSchoenfliesSymbol(const char *SfSymbol)
{
  int         SgNumber;
  const char  **TabSymbol;
  const char  *s, *t;


  TabSymbol = SchoenfliesSymbols + 1;

  for (SgNumber = 1; SgNumber <= 230; SgNumber++)
  {
    t = *TabSymbol;
    s = SfSymbol;

    while (*t && *s)
    {
      if (   toupper(*t) != toupper(*s)
          && (*t != '^' || isalpha(*s) || isdigit(*s)))
        break;

      t++;
      s++;
    }

    if (*t == *s)
      return SgNumber;

    TabSymbol++;
  }

  return -1;
}


static int SgLabelCmp(const int SgNumber,
                      const char *SgLabel, const char *WtdLbl)
{
  const char  *sgl, *wl;


  /* first try: plain strcmp
   */
  sgl = SgLabel;

  for (wl = WtdLbl; ; wl++)
  {
    SkipWhite(wl);
    SkipWhite(sgl);

    if (*sgl == '\0' || *sgl == '=')
    {
      if (*wl == '\0') return 0;
      break;
    }

    if (*sgl == '-')
    {
      if (*wl != '-' && toupper(*wl) != 'B')
        break;
    }
    else if (toupper(*sgl) != toupper(*wl))
      break;

    sgl++;
  }

  /* second try: swap the dash (there should be only one)
   */
  sgl = SgLabel;

  for (wl = WtdLbl; ; wl++)
  {
    SkipWhite(wl);
    SkipWhite(sgl);

    if (*sgl == '-')
    {
      if (wl[1] != '-' && toupper(wl[1]) != 'B')
        break;
      if (toupper(sgl[1]) != toupper(*wl))
        break;

      sgl++;
      wl++;
    }
    else
    {
      if (*sgl == '\0' || *sgl == '=')
      {
        if (*wl == '\0') return 0;
        break;
      }

      if (toupper(*sgl) != toupper(*wl))
        break;
    }

    sgl++;
  }

  if (SgNumber >= 195) /* cubic space groups only */
  {
    /* third try: ignore the "-3" dash
     */
    sgl = SgLabel;

    for (wl = WtdLbl; ; wl++)
    {
      SkipWhite(wl);
      SkipWhite(sgl);

      if (*sgl == '-' && sgl[1] == '3')
        sgl++;

      if (*sgl == '\0' || *sgl == '=')
      {
        if (*wl == '\0') return 0;
        break;
      }

      if (toupper(*sgl) != toupper(*wl))
        break;

      sgl++;
    }
  }

  return -1;
}


static int ParseExtension(const char *Ext, T_ExtInfo *ExtInfo)
{
  int         i, mode;
  const char  *e, *t;


  ExtInfo->OriginChoice =
  ExtInfo->CellChoice =
  ExtInfo->BasisChoice = ' ';
  ExtInfo->BT_or_UA = "";

  mode = 0;

  while (*Ext)
  {
    if      (strchr("12",   *Ext) != NULL)
    {
      ExtInfo->CellChoice   =
      ExtInfo->OriginChoice = *Ext++;
    }
    else if (strchr("3",    *Ext) != NULL)
    {
      ExtInfo->CellChoice   = *Ext++;
    }
    else if (strchr("Ss",   *Ext) != NULL)
    {
      ExtInfo->OriginChoice = '1';
      Ext++;
    }
    else if (strchr("Zz",   *Ext) != NULL)
    {
      ExtInfo->OriginChoice = '2';
      Ext++;
    }
    else if (strchr("Hh",   *Ext) != NULL)
    {
      ExtInfo->BasisChoice = 'H';
      Ext++;
    }
    else if (strchr("Rr",   *Ext) != NULL)
    {
      ExtInfo->BasisChoice = 'R';
      Ext++;
    }
    else if (mode == 0)
      mode = 1;

    if (mode == 2)
      break;

    for (i = 0; Ext_BT_or_UA[i]; i++)
    {
      for (e = Ext, t = Ext_BT_or_UA[i]; *t; e++, t++)
        if (toupper(*e) != toupper(*t))
          break;

      if (*t == '\0')
      {
        if      (6 <= i && i <=  8)
          i = 2 * i - 11;
        else if (9 <= i && i <= 20)
          i = 9 + ((i - 9) / 4) * 4;

        ExtInfo->BT_or_UA = Ext_BT_or_UA[i];
        Ext = e;
        break;
      }
    }

    if (mode == 0)
      break;

    mode = 2;
  }

  if (*Ext)
    return -1;

  return 0;
}


static void ExpandMonoclinic(int unique_axis, const char *o, char *m)
{
  if (*o) *m++ = *o++;

  switch (tolower(unique_axis))
  {
    case 'a':
      while (*o) *m++ = *o++;
      *m++ = '1';
      *m++ = '1';
      break;
    case 'c':
      *m++ = '1';
      *m++ = '1';
      while (*o) *m++ = *o++;
      break;
    default:
      *m++ = '1';
      while (*o) *m++ = *o++;
      *m++ = '1';
      break;
  }

  *m = '\0';
}


const T_TabSgName *FindTabSgNameEntry(const char *UserSgName, int VolLetter)
{
#define                   MaxWtdLbl 20
  char     WtdLblOriginal[MaxWtdLbl    + 1];
  char     WtdLblModified[MaxWtdLbl    + 1];
  char    *WtdLbl;
#define          MaxWtdExt  5
  char    WtdExt[MaxWtdExt    + 1];
  int     WtdSgNumber;
  int     WtdLblOriginChoice;
  int     WtdLblBasisChoice;
  int     iwl, iwe;
  char    *wl, *we;

  int                i, IsExpanded, lbl_match;
  const char         *sgl;
  const T_TabSgName  *tsgn;
  int                 WtdCC;
  const char         *WtdUA;
  char                WtdUA_Buf[2];
  T_ExtInfo          ExtInfo, WtdExtInfo;


  if      (VolLetter == 0 || isspace(VolLetter))
    VolLetter = 'A';
  else if (VolLetter == '1')
    VolLetter = 'I';
  else
  {
           VolLetter = toupper(VolLetter);
    if (   VolLetter != 'I'
        && VolLetter != 'A')
      return NULL;
  }

  WtdLbl = WtdLblOriginal;

   wl = WtdLbl;
  iwl = 0;

  while (*UserSgName && *UserSgName != ':')
  {
    if (isspace(*UserSgName) == 0 && *UserSgName != '_')
    {
      if (iwl >= MaxWtdLbl)
        return NULL;

      *wl++ = *UserSgName;
      iwl++;
    }

    UserSgName++;
  }

  *wl = '\0';

  if (iwl == 0)
    return NULL;

   we = WtdExt;
  iwe = 0;
  *we = '\0';

  if (*UserSgName)
  {
    UserSgName++;

    while (*UserSgName)
    {
      if (isspace(*UserSgName) == 0 && *UserSgName != '_')
      {
        if (iwe >= MaxWtdExt)
          return NULL;

        *we++ = *UserSgName;
        iwe++;
      }

      UserSgName++;
    }
  }

  *we = '\0';

  WtdLblOriginChoice = ' ';
  WtdLblBasisChoice  = ' ';

  if (iwl > 1)
  {
    wl = &WtdLbl[iwl - 1];

    if      (*wl == 'S' || *wl == 's')
    { WtdLblOriginChoice = '1'; *wl = '\0'; iwl--; }
    else if (*wl == 'Z' || *wl == 'z')
    { WtdLblOriginChoice = '2'; *wl = '\0'; iwl--; }
    else if (*wl == 'H' || *wl == 'h')
    { WtdLblBasisChoice  = 'H'; *wl = '\0'; iwl--; }
    else if (*wl == 'R' || *wl == 'r')
    { WtdLblBasisChoice  = 'R'; *wl = '\0'; iwl--; }
  }

  if (isalpha(WtdLbl[0]))
    WtdSgNumber = FindSchoenfliesSymbol(WtdLbl);
  else
  {
    for (wl = WtdLbl; *wl; wl++)
      if (isdigit(*wl) == 0)
        return NULL;

    if (   sscanf(WtdLbl, "%d", &WtdSgNumber) != 1
        || WtdSgNumber <   1
        || WtdSgNumber > 230)
      return NULL;
  }

  if (ParseExtension(WtdExt, &WtdExtInfo) != 0)
    return NULL;

  if      (WtdExtInfo.OriginChoice == ' ')
           WtdExtInfo.OriginChoice =  WtdLblOriginChoice;
  else if (WtdExtInfo.OriginChoice != WtdLblOriginChoice
                                   && WtdLblOriginChoice != ' ')
    return NULL;

  if      (WtdExtInfo.BasisChoice == ' ')
           WtdExtInfo.BasisChoice =  WtdLblBasisChoice;
  else if (WtdExtInfo.BasisChoice != WtdLblBasisChoice
                                  && WtdLblBasisChoice != ' ')
    return NULL;

  if (   WtdExtInfo.OriginChoice != ' '
      && WtdExtInfo.BasisChoice  != ' ')
    return NULL;

  for (IsExpanded = 0; IsExpanded < 4; IsExpanded++)
  {
    for (tsgn = TabSgName; tsgn->HallSymbol; tsgn++)
    {
      if (   IsExpanded != 0
          && tsgn->SgNumber > 15)
        break;

      lbl_match = 0;

      if (WtdSgNumber == -1)
      {
        i = 1;
                sgl = tsgn->SgLabels;
        while (*sgl && i <= 2)
        {
          while (*sgl && strchr(" =\t", *sgl) != NULL) sgl++;

          if (SgLabelCmp(tsgn->SgNumber, sgl, WtdLbl) == 0)
          {
            lbl_match = i;
            break;
          }

          while (*sgl && strchr(" =\t", *sgl) == NULL) sgl++;

          i++;
        }
      }

      if (ParseExtension(tsgn->Extension, &ExtInfo) != 0) {
        SetSgError(IErr_Corrupt_TabSgName);
        return NULL;
      }

      if (WtdSgNumber == tsgn->SgNumber || lbl_match != 0)
      {
        if (   tsgn->SgNumber >=  3
            && tsgn->SgNumber <  16)
        {
          if (   WtdLblOriginChoice != ' '
              || WtdExtInfo.BasisChoice != ' '
              || (int) strlen(WtdExtInfo.BT_or_UA) > 2)
            continue; /* next tsgn */

          if (WtdSgNumber == tsgn->SgNumber)
          {
            if (WtdExtInfo.BT_or_UA[0])
              WtdUA = WtdExtInfo.BT_or_UA;
            else if (VolLetter == 'I')
            {
              if (   ExtInfo.BT_or_UA[0] != 'c'
                  && ExtInfo.BT_or_UA[1] != 'c')
                continue; /* next tsgn */

              if (   ExtInfo.CellChoice == ' '
                  && (   WtdExtInfo.CellChoice == ' '
                      || WtdExtInfo.CellChoice == '1'))
                return tsgn;

              i = 0;
              for (sgl = tsgn->SgLabels; *sgl; sgl++)
                if (*sgl == '=') i++;

              if (   i == 2
                  && (   WtdExtInfo.CellChoice == ' '
                      || WtdExtInfo.CellChoice == ExtInfo.CellChoice))
                return tsgn;

              continue; /* next tsgn */
            }
            else
              WtdUA = "b";
          }
          else /* if (lbl_match != 0) */
          {
            if (WtdExtInfo.BT_or_UA[0])
              WtdUA = WtdExtInfo.BT_or_UA;
            else if (lbl_match > 1)
              WtdUA = ExtInfo.BT_or_UA;
            else if (   VolLetter == 'I'
                     && ExtInfo.CellChoice == ' ')
              WtdUA = "c";
            else
              WtdUA = "b";
          }

          if (WtdExtInfo.CellChoice != ' ')
            WtdCC = WtdExtInfo.CellChoice;
          else if (ExtInfo.CellChoice == '1')
            WtdCC = ExtInfo.CellChoice;
          else
            WtdCC = ' ';

          if (strcmp(ExtInfo.BT_or_UA, WtdUA) == 0)
          {
            if (WtdCC == ' ' && lbl_match > 1)
              return tsgn;
            if (ExtInfo.CellChoice == WtdCC)
              return tsgn;
            if (ExtInfo.CellChoice == ' ' && WtdCC == '1')
              return tsgn;
            if (ExtInfo.CellChoice == '1' && WtdCC == ' ')
              return tsgn;
          }
        }
        else if (ExtInfo.BasisChoice != ' ')
        {
          if (   WtdExtInfo.OriginChoice != ' '
              || WtdExtInfo.CellChoice != ' '
              || WtdExtInfo.BT_or_UA[0] != '\0')
            continue; /* next tsgn */

          if (ExtInfo.BasisChoice == WtdExtInfo.BasisChoice)
            return tsgn;

          if (WtdExtInfo.BasisChoice == ' ')
          {
            if (ExtInfo.BasisChoice == 'R' && VolLetter == 'I')
              return tsgn;
            if (ExtInfo.BasisChoice == 'H' && VolLetter != 'I')
              return tsgn;
          }
        }
        else if (WtdExtInfo.BasisChoice == ' ')
        {
          if ( (WtdExtInfo.OriginChoice == ' ' && ExtInfo.OriginChoice == '1')
            || (WtdExtInfo.OriginChoice == '1' && ExtInfo.OriginChoice == ' ')
            ||  WtdExtInfo.OriginChoice ==        ExtInfo.OriginChoice)
          {
            if (WtdExtInfo.BT_or_UA[0])
            {
              if (WtdExtInfo.BT_or_UA == ExtInfo.BT_or_UA)
                return tsgn;
              if (   WtdExtInfo.BT_or_UA == Ext_BT_or_UA[0]
                  &&    ExtInfo.BT_or_UA[0] == '\0')
                return tsgn;
            }
            else
            {
              if (lbl_match != 0)
                return tsgn;
              if (ExtInfo.BT_or_UA[0] == '\0')
                return tsgn;
            }
          }
        }
      }
    }

    if (WtdSgNumber != -1)
      return NULL;

    if ((int) strlen(WtdExtInfo.BT_or_UA) > 2)
      return NULL;

    if (IsExpanded == 0)
    {
      iwl += 2;

      if (iwl > MaxWtdLbl)
        IsExpanded = 2;
      else
      {
        if (WtdExtInfo.BT_or_UA[0])
          WtdUA = WtdExtInfo.BT_or_UA;
        else
        {
          if (VolLetter == 'I')
            WtdUA = "c";
          else
            WtdUA = "b";
        }

        ExpandMonoclinic(WtdUA[0], WtdLblOriginal, WtdLblModified);

        WtdLbl = WtdLblModified;
      }
    }
    else if (IsExpanded == 1)
    {
      if (WtdExtInfo.BT_or_UA[0])
        return NULL;

      if (VolLetter == 'I')
        WtdUA = "b";
      else
        WtdUA = "c";

      ExpandMonoclinic(WtdUA[0], WtdLblOriginal, WtdLblModified);
    }

    if (IsExpanded == 2)
    {
      if (WtdExtInfo.BT_or_UA[0])
        return NULL;

      iwl -= 2;

      if (iwl < 2)
        return NULL;
                                            iwl--;
      WtdUA_Buf[0] = tolower(WtdLblOriginal[iwl]);
                             WtdLblOriginal[iwl] = '\0';
      WtdUA_Buf[1] = '\0';

      if (strchr("abc", WtdUA_Buf[0]) == NULL)
        return NULL;

      WtdUA = WtdUA_Buf;

      iwl += 2;

      if (iwl > MaxWtdLbl)
        return NULL;

      ExpandMonoclinic(WtdUA[0], WtdLblOriginal, WtdLblModified);

      WtdLbl = WtdLblModified;
    }
  }

  return NULL;
}


unsigned int SgID_Number(const T_TabSgName *tsgn)
{
  unsigned int  ID;
  int           iBT;
  const char    *UA;
  T_ExtInfo     ExtInfo;


  ID = tsgn->SgNumber;

  if (ParseExtension(tsgn->Extension, &ExtInfo) != 0)
    ID = 0;

  if (ID >= 3 && ID < 16)
  {
    UA = ExtInfo.BT_or_UA;

    if (   *UA != 'b'
        || (   ExtInfo.CellChoice != ' '
            && ExtInfo.CellChoice != '1'))
    {
      if (*UA == '-')
      {
        ID += 3000;
        UA++;
      }

      switch (*UA)
      {
        case 'b': ID += 10000; break;
        case 'c': ID += 20000; break;
        case 'a': ID += 30000; break;
        default:  ID = 0;      break;
      }

      if (ID != 0)
      {
        switch (ExtInfo.CellChoice)
        {
          case ' ':             break;
          case '1': ID += 1000; break;
          case '2': ID += 2000; break;
          case '3': ID += 3000; break;
          default:  ID = 0;     break;
        }
      }
    }
  }
  else
  {
    if (ExtInfo.BasisChoice == 'R')
      ID += 20000;
    else
    {
      if (ExtInfo.BT_or_UA[0])
      {
        for (iBT = 0; iBT < 6; iBT++)
          if (ExtInfo.BT_or_UA == Ext_BT_or_UA[iBT])
            break;
      }
      else
        iBT = 0;

      if (iBT < 6)
      {
        if (ExtInfo.OriginChoice == '2') ID += 20000;
        else if (iBT)                    ID += 10000;

        if (iBT)
          ID += (iBT + 1) * 1000;
      }
      else
        ID = 0;
    }
  }

  if (ID == 0)
    SetSgError(IErr_Corrupt_TabSgName);

  return ID;
}


int ParseSymXYZ(const char *SymXYZ, T_RTMx *SeitzMx, int FacTr)
{
  unsigned int  P_mode;
  int           Row, Column, Sign, GotXYZ, i;
  double        Value, Value1, Value2, Delta;


  for (i = 0; i < 12; i++) SeitzMx->a[i] = 0;

#define P_Blank   0x01u
#define P_Comma   0x02u
#define P_Plus    0x04u
#define P_Dash    0x08u
#define P_Slash   0x10u
#define P_Value1  0x20u
#define P_Value2  0x40u
#define P_XYZ     0x80u

  Value1 = 0.;

  Row    = 0;
  Sign   = 1;
  Value  = 0.;
  GotXYZ = 0;
  P_mode = P_Blank | P_Plus | P_Dash | P_Value1 | P_XYZ;

  do
  {
    switch (*SymXYZ)
    {
      case ' ':
      case '\t':
      case '_':
        if ((P_mode & P_Blank) == 0) return -1;
        break;
      case ',':
      case ';':
        if (Row == 2)                return -1;/* FALLTHRU *//*FALLTHRU added by NCrystal developers for gcc 7.2.1 */
      case '\0':
        if ((P_mode & P_Comma) == 0) return -1;
        if (GotXYZ == 0)             return -1;
        if (P_mode & P_Slash) Value += Value1;
        Value *= FacTr;
        if (Value < 0.) i = Value - .5;
        else            i = Value + .5;
        Delta = Value - i;
        if (Delta < 0.) Delta = -Delta;
        if (Delta > .01 * FacTr) return -1;
        i %= FacTr; if (i < 0) i += FacTr;
        SeitzMx->s.T[Row] = i;
        Row++;
        Sign   = 1;
        Value  = 0.;
        P_mode = P_Blank | P_Plus | P_Dash | P_Value1 | P_XYZ;
        GotXYZ = 0;
        break;
      case '+':
        if ((P_mode & P_Plus)  == 0) return -1;
        if (P_mode & P_Slash) Value += Value1;
        Sign =  1;
        if (P_mode & P_Value2)
          P_mode = P_Value2;
        else
          P_mode = P_Blank | P_Value1 | P_XYZ;
        break;
      case '-':
        if ((P_mode & P_Dash)  == 0) return -1;
        if (P_mode & P_Slash) Value += Value1;
        Sign = -1;
        if (P_mode & P_Value2)
          P_mode = P_Value2;
        else
          P_mode = P_Blank | P_Value1 | P_XYZ;
        break;
      case '/':
      case ':':
        if ((P_mode & P_Slash) == 0) return -1;
        Sign =  1;
        P_mode = P_Blank | P_Plus | P_Dash | P_Value2;
        break;
      case '.':
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        if      (P_mode & P_Value1)
        {
          if (sscanf(SymXYZ, "%lf%n", &Value1, &i) != 1) return -1;
          if (Sign == -1) Value1 = -Value1;
          P_mode = P_Blank | P_Comma | P_Plus | P_Dash | P_Slash;
        }
        else if (P_mode & P_Value2)
        {
          if (sscanf(SymXYZ, "%lf%n", &Value2, &i) != 1) return -1;
          if (Sign == -1) Value2 = -Value2;
          if (Value1 != 0.)
          {
            if (Value2 == 0.) return -1;
            Value += Value1 / Value2;
          }
          P_mode = P_Blank | P_Comma | P_Plus | P_Dash;
        }
        else
          return -1;
        SymXYZ += (i - 1);
        break;
      case 'X':
      case 'x': Column = 0; goto Process_XYZ;
      case 'Y':
      case 'y': Column = 1; goto Process_XYZ;
      case 'Z':
      case 'z': Column = 2;
       Process_XYZ:
        if ((P_mode & P_XYZ) == 0) return -1;
        i = Row * 3 + Column;
        if (SeitzMx->s.R[i] != 0) return -1;
        SeitzMx->s.R[i] = Sign;
        GotXYZ = 1;
        P_mode = P_Blank | P_Comma | P_Plus | P_Dash;
        break;
    }
  }
  while (*SymXYZ++);

  if (Row != 3) return -1;

  return 0;

#undef P_Blank
#undef P_Comma
#undef P_Plus
#undef P_Dash
#undef P_Slash
#undef P_Value1
#undef P_Value2
#undef P_XYZ
}


static int LookupRotMx(T_HallGenerator *HG)
{
  int                   i, f, refaxis, dircode;
  int                   iNextBasis, nNextBasis;
  const T_TabXtalRotMx  *txrmx;


  if (HG->Rotation <= 0) return 0;

  refaxis = HG->RefAxis;
  dircode = HG->DirCode;

  if (HG->Rotation == 1)
  {
    /* Next line commented by NCrystal developers to avoid static analysis error (value assigned to refaxis never used): */
    /*refaxis = 'o';*/
    dircode = '.';
    nNextBasis = 0;
  }
  else if (dircode == '*')
  {
    /* Next line commented by NCrystal developers to avoid static analysis error (value assigned to refaxis never used): */
    /*if (refaxis == 0) refaxis = 'o';*/
    nNextBasis = 0;
  }
  else
  {
    if (dircode == 0) dircode = '=';

    switch (refaxis)
    {
      case 'z': nNextBasis = 0; break;
      case 'x': nNextBasis = 1; break;
      case 'y': nNextBasis = 2; break;
      default:
        return 0;
    }
  }

  for (txrmx = TabXtalRotMx; txrmx->Order; txrmx++)
    if (txrmx->Order == HG->Rotation) break;

  while (txrmx->Order == HG->Rotation)
  {
    if (txrmx->DirCode == dircode)
    {
      if (HG->Improper == 0) f =  1;
      else                   f = -1;

      for (i = 0; i < 9; i++)
        HG->SeitzMx.s.R[i] = txrmx->RMx[i] * f;

      for (iNextBasis = 0; iNextBasis < nNextBasis; iNextBasis++)
        RotateRotMx(HG->SeitzMx.s.R, RMx_3_111, RMx_3i111);

      return 1;
    }

    txrmx++;
  }

  return 0;
}


int ParseHallSymbol(const char *hsym, T_SgInfo *SgInfo)
{
  int                  c, i, pos_hsym;
  const int            *ht;
  int                  Centric;
  const T_LatticeInfo  *LatticeInfo;
  int                  FieldType, PreviousFT;
  int                  iOriginShift, SignOriginShift;
  int                  digit, rotation, refaxis, dircode;
  const int            *translation;
  int                  PreviousRotation, PreviousRefAxis;
  int                  nHG, ClearHG;
  T_HallGenerator      HG;

  enum ListOfFieldTypes
    {
      FT_Delimiter,
      FT_Improper,
      FT_Digit,
      FT_Rotation,
      FT_RefAxis,
      FT_DirCode,
      FT_Translation,
      FT_OriginShift
    };

  static const char *Err_Ill_ori_shi_val =
    "Error: Illegal origin shift value";

  static const char *Err_Too_ori_shi_val =
    "Error: Too much origin shift values";


  Centric = 0;
  LatticeInfo = NULL;

  HG.Rotation = HG.RefAxis = HG.DirCode = HG.Screw = 0;

  nHG = 0;
  ClearHG = 1;
  FieldType = FT_Delimiter;
  PreviousRotation = 0;
  PreviousRefAxis = 0;
  iOriginShift = 0;
  SignOriginShift = 0;

  pos_hsym = 0;

  do
  {
    if (*hsym == '_' || *hsym == '.' || *hsym == '\t' || *hsym == '\0')
      c = ' ';
    else
      c = *hsym;

    pos_hsym++;

    if (LatticeInfo == NULL)
    {
      if (Centric == 0 && c == '-')
      {
        if (AddInversion2ListSeitzMx(SgInfo) < 0)
          return pos_hsym;
        Centric = 1;
      }
      else if (c != ' ')
      {
        c = toupper(c);

        switch (c)
        {
          case 'P': LatticeInfo = LI_P; break;
          case 'A': LatticeInfo = LI_A; break;
          case 'B': LatticeInfo = LI_B; break;
          case 'C': LatticeInfo = LI_C; break;
          case 'I': LatticeInfo = LI_I; break;
          case 'R': LatticeInfo = LI_R; break;
          case 'S': LatticeInfo = LI_S; break;
          case 'T': LatticeInfo = LI_T; break;
          case 'F': LatticeInfo = LI_F; break;
          default:
            SetSgError("Error: Illegal lattice code");
            return pos_hsym;
        }

        if (AddLatticeTr2ListSeitzMx(SgInfo, LatticeInfo) < 0)
          return pos_hsym;
      }
    }
    else if (FieldType != FT_OriginShift)
    {
      c = tolower(c);
      if      (c == 'q') c = '\'';
      else if (c == '+') c = '"';

      PreviousFT = FieldType;
      digit = rotation = refaxis = dircode = 0;
      translation = NULL;

      ht = HallTranslations;

      while (*ht)
      {
        if (c == *ht)
        {
          translation = ht;
          FieldType = FT_Translation;
          break;
        }
        ht += 4;
      }

      if (translation == NULL)
      {
        switch (c)
        {
          case  ' ': FieldType = FT_Delimiter; break;

          case  '-': FieldType = FT_Improper; break;

          case  '1': digit = 1; FieldType = FT_Digit; break;
          case  '2': digit = 2; FieldType = FT_Digit; break;
          case  '3': digit = 3; FieldType = FT_Digit; break;
          case  '4': digit = 4; FieldType = FT_Digit; break;
          case  '5': digit = 5; FieldType = FT_Digit; break;
          case  '6': digit = 6; FieldType = FT_Digit; break;

          case  'x':
          case  'y':
          case  'z': refaxis = c; FieldType = FT_RefAxis; break;

          case  '"':
          case '\'':
          case  '*': dircode = c; FieldType = FT_DirCode; break;

          case  '(': FieldType = FT_OriginShift; break;

          default:
            SetSgError("Error: Illegal character in Hall symbol");
            return pos_hsym;
        }

        if (FieldType == FT_Digit)
        {
          if (   ClearHG == 0
              && HG.Rotation > digit
              && HG.Screw == 0
              && HG.DirCode == 0)
          {
            HG.Screw = digit;
            FieldType = FT_Translation;
          }
          else if (digit == 5)
          {
            SetSgError("Error: Illegal 5-fold rotation");
            return pos_hsym;
          }
          else
          {
            rotation = digit;
            FieldType = FT_Rotation;
          }
        }
      }

      if (   ClearHG == 0
          && (    FieldType == FT_Delimiter
              ||  FieldType == FT_OriginShift
              ||  FieldType  < PreviousFT
              || (FieldType == PreviousFT && FieldType != FT_Translation))
          && ! (   FieldType == FT_RefAxis && HG.RefAxis == 0
                && PreviousFT == FT_DirCode))
      {
        if (HG.RefAxis == 0)
        {
          if (nHG == 0)
            HG.RefAxis = 'z';
          else
          {
            if (HG.Rotation == 2)
            {
              if      (PreviousRotation == 2 || PreviousRotation == 4)
                HG.RefAxis = 'x';
              else if (PreviousRotation == 3 || PreviousRotation == 6)
              {
                HG.RefAxis = PreviousRefAxis;
                if (HG.DirCode == 0) HG.DirCode = '\'';
              }
            }
            else if (HG.Rotation == 3)
            {
              if (HG.DirCode == 0) HG.DirCode = '*';
            }
          }
        }

        PreviousRefAxis = HG.RefAxis;
        PreviousRotation = HG.Rotation;

        if (LookupRotMx(&HG) == 0)
        {
          SetSgError("Error: Illegal generator or need explicit axis symbol");
          return pos_hsym - 1;
        }

        if (HG.Screw)
        {
          switch (HG.RefAxis)
          {
            case 'x': i =  0; break;
            case 'y': i =  1; break;
            case 'z': i =  2; break;
            default:  i = -1; break;
          }

          if (HG.DirCode != 0 || i < 0)
          {
            SetSgError("Error: Screw for non-principal direction");
            return pos_hsym - 1;
          }

          HG.SeitzMx.s.T[i] += STBF * HG.Screw / HG.Rotation;
        }

        for (i = 0; i < 3; i++)
          HG.SeitzMx.s.T[i] %= STBF;

        if (Add2ListSeitzMx(SgInfo, &HG.SeitzMx) < 0)
          return pos_hsym - 1;

        if (SgInfo->StatusLatticeTr == -1)
        {
          if (AddLatticeTr2ListSeitzMx(SgInfo, SgInfo->LatticeInfo) < 0)
            return pos_hsym - 1;
        }

        nHG++;
        ClearHG = 1;
      }

      if (FieldType != FT_Delimiter && FieldType != FT_OriginShift)
      {
        if (ClearHG)
        {
          HG.Improper = 0;
          HG.Rotation = 1;
          HG.RefAxis = 0;
          HG.DirCode = 0;
          HG.Screw = 0;
          for (i = 0; i < 12; i++) HG.SeitzMx.a[i] = 0;

          ClearHG = 0;
        }

        switch (FieldType)
        {
          case FT_Improper:    HG.Improper = 1;        break;
          case FT_Rotation:    HG.Rotation = rotation; break;
          case FT_RefAxis:     HG.RefAxis  = refaxis;  break;
          case FT_DirCode:     HG.DirCode  = dircode;  break;
          case FT_Translation:
            if (translation != NULL)
            {
              for (i = 0; i < 3; i++)
                HG.SeitzMx.s.T[i] += *(++translation);
            }
            break;
        }
      }
    }
    else /* FieldType == FT_OriginShift */
    {
      if (iOriginShift > 3) {
        SetSgError(Err_Too_ori_shi_val);
        return pos_hsym;
      }

      if (*hsym == '\0') c = ')';

      digit = -1;

      switch (c)
      {
        case ' ': break;

        case ')':
          if (iOriginShift != 3)
          {
            SetSgError("Error: Missing origin shift values");
            return pos_hsym;
          }
          iOriginShift++;
          FieldType = FT_Delimiter;
          break;

        case '-':
          if (SignOriginShift != 0) {
            SetSgError(Err_Ill_ori_shi_val);
            return pos_hsym;
          }
          SignOriginShift = 1;
          break;

        case '0': digit = 0; break;
        case '1': digit = 1; break;
        case '2': digit = 2; break;
        case '3': digit = 3; break;
        case '4': digit = 4; break;
        case '5': digit = 5; break;
        case '6': digit = 6; break;

        default:
          SetSgError(Err_Ill_ori_shi_val);
          return pos_hsym;
      }

      if (digit >= 0)
      {
        if (iOriginShift >= 3) {
          SetSgError(Err_Too_ori_shi_val);
          return pos_hsym;
        }
        if (SignOriginShift) digit *= -1;
        SignOriginShift = 0;
        SgInfo->OriginShift[iOriginShift++] = digit;
      }
    }
  }
  while (*hsym++ != '\0');

  if (LatticeInfo == NULL) {
    SetSgError("Error: Lattice type not specified");
    return pos_hsym;
  }

  return pos_hsym;
}


static const char *PrintSgLabel(const char *lbl, int space, int *n,
                                FILE *fpout)
{
  while (*lbl && *lbl != ' ')
  {
    if (*lbl == '_')
    {
      if (space)
      {
        putc(space, fpout);
        if (n) (*n)++;
      }
    }
    else
    {
      putc(*lbl, fpout);
      if (n) (*n)++;
    }

    lbl++;
  }

  return lbl;
}


int PrintFullHM_SgName(const T_TabSgName *tsgn, int space, FILE *fpout)
{
  int         n;
  const char  *lbl;


  lbl = tsgn->SgLabels;

  if (tsgn->SgNumber >= 3 && tsgn->SgNumber < 16)
    while (*lbl) if (*lbl++ == '=') break;

  SkipWhite(lbl);

  n = 0;

  PrintSgLabel(lbl, space, &n, fpout);

  lbl = tsgn->Extension;

  if (*lbl && strchr("12HhRr", *lbl))
  {
    putc(':', fpout);
    putc(*lbl, fpout);
    n += 2;
  }

  return n;
}


void PrintTabSgNameEntry(const T_TabSgName *tsgn, int Style, int space,
                         FILE *fpout)
{
  int         n;
  const char  *lbl, *SfSymbol;


  if (Style)
    n = fprintf(fpout, "%3d", tsgn->SgNumber);
  else
    n = fprintf(fpout,  "%d", tsgn->SgNumber);

  if (tsgn->Extension[0])
    n += fprintf(fpout, ":%s", tsgn->Extension);

  if (Style)
    while (n < 9) { putc(' ', fpout); n++; }

  putc(' ', fpout); n++;
  putc(' ', fpout); n++;

  if (tsgn->SgNumber >= 1 && tsgn->SgNumber <= 230)
    SfSymbol = SchoenfliesSymbols[tsgn->SgNumber];
  else
    SfSymbol = "";

  n += fprintf(fpout, "%s", SfSymbol);

  if (Style)
    while (n < 23) { putc(' ', fpout); n++; }

  putc(' ', fpout); n++;
  putc(' ', fpout); n++;

  if (tsgn->SgNumber >= 3 && tsgn->SgNumber < 16)
  {
    lbl = PrintSgLabel(tsgn->SgLabels, space, &n, fpout);

    if (tsgn->Extension[0])
      n += fprintf(fpout, ":%s", tsgn->Extension);

    putc(' ', fpout); putc('=', fpout); putc(' ', fpout); n += 3;

    n += PrintFullHM_SgName(tsgn, space, fpout);

    while (*lbl) if (*lbl++ == '=') break;
    while (*lbl) if (*lbl++ == '=') break;
    SkipWhite(lbl);

    if (*lbl)
    {
      putc(' ', fpout); putc('=', fpout); putc(' ', fpout); n += 3;

      PrintSgLabel(lbl, space, &n, fpout);
    }
  }
  else
    n += PrintFullHM_SgName(tsgn, space, fpout);

  if (Style)
    while (n < 51) { putc(' ', fpout); n++; }

  putc(' ', fpout);
  putc(' ', fpout);

  fprintf(fpout, "%s", tsgn->HallSymbol);
}


static int FindGCD2(int ri, int rj)
{
  int  rk;


  if (ri < 0) ri = -ri;

  if (rj)
  {
    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    if (ri < 0) ri = -ri;
  }

  return ri;
}


static void SimplifyFraction(int nume, int deno, int *o_nume, int *o_deno)
{
  int gcd = FindGCD2(nume, deno);
  if (gcd)
  {
    *o_nume = nume / gcd;
    *o_deno = deno / gcd;

    if (*o_deno < 0) {
      *o_nume *= -1;
      *o_deno *= -1;
    }
  }
}


const char *FormatFraction(int nume, int deno, int Decimal,
                           char *Buffer, int SizeBuffer)
{
  int          n=0;
  int          d=0;
  char         *cp, *cpp;
  static char  StaticBuffer[40];


  if (NULL == Buffer) {
              Buffer =        StaticBuffer;
          SizeBuffer = sizeof StaticBuffer / sizeof (*StaticBuffer);
  }

  Buffer[SizeBuffer - 1] = '\0';

  if (nume == 0)
  {
    Buffer[0] = '0';
    Buffer[1] = '\0';
  }
  if (Decimal)
  {
    sprintf(Buffer, "%.6g", (double) nume / deno);

         cp = Buffer;
    if (*cp == '-') cp++;
    if (*cp == '0') {
      cpp = cp + 1; while (*cp) *cp++ = *cpp++;
    }
  }
  else
  {
    SimplifyFraction(nume, deno, &n, &d);

    if (d == 1)
      sprintf(Buffer, "%d", n);
    else
      sprintf(Buffer, "%d/%d", n, d);
  }

  if (Buffer[SizeBuffer - 1] != '\0') {
      Buffer[SizeBuffer - 1] =  '\0';
    SetSgError("Internal Error: FormatFraction(): Buffer too small");
    return NULL;
  }

  return Buffer;
}


const char *RTMx2XYZ(const T_RTMx *RTMx, int FacRo, int FacTr,
                     int Decimal, int TrFirst, int Low,
                     const char *Seperator,
                     char *BufferXYZ, int SizeBufferXYZ)
{
  static const char *UpperXYZ = "XYZ";
  static const char *LowerXYZ = "xyz";

  int         i, j, p, iRo, iTr;
  char        *xyz, buf_tr[32];
  const char  *sep, *LetterXYZ, *ro, *tr;

  static char  StaticBufferXYZ[80];


  if (NULL == BufferXYZ) {
              BufferXYZ  =        StaticBufferXYZ;
          SizeBufferXYZ  = sizeof StaticBufferXYZ / sizeof (*StaticBufferXYZ);
  }

  BufferXYZ[SizeBufferXYZ - 1] = '\0';

  if (Low)
    LetterXYZ = LowerXYZ;
  else
    LetterXYZ = UpperXYZ;

  if (Seperator == NULL)
      Seperator = ",";

  xyz = BufferXYZ;

  for (i = 0; i < 3; i++)
  {
    if (i != 0)
      for (sep = Seperator; *sep; sep++) *xyz++ = *sep;

        iTr = iModPositive(RTMx->s.T[i], FacTr);
    if (iTr >  FacTr / 2)
        iTr -= FacTr;

    tr = FormatFraction(iTr, FacTr, Decimal,/*indentation of this line fixed by NCrystal developers*/
                            buf_tr, sizeof buf_tr / sizeof (*buf_tr));
    if (tr == NULL)
      return NULL;

    p = 0;

    if (  TrFirst && iTr) {
      if (*tr) p = 1;
      while (*tr) *xyz++ = *tr++;
    }

    for (j = 0; j < 3; j++)
    {
          iRo = RTMx->s.R[i * 3 + j];
      if (iRo)
      {
            ro = FormatFraction(iRo, FacRo, Decimal, NULL, 0);
        if (ro == NULL)
          return NULL;

        if      (*ro == '-')
          *xyz++ = *ro++;
        else if (*ro && p)
          *xyz++ = '+';

        if (ro[0] != '1' || ro[1] != '\0') {
          while (*ro) *xyz++ = *ro++;
          *xyz++ = '*';
        }

        *xyz++ = LetterXYZ[j];

        p = 1;
      }
    }

    if (! TrFirst && iTr)
    {
      if (*tr && *tr != '-' && p)
        *xyz++ = '+';

      while (*tr) *xyz++ = *tr++;
    }
  }

  *xyz = '\0';

  if (BufferXYZ[SizeBufferXYZ - 1] != '\0') {
      BufferXYZ[SizeBufferXYZ - 1] =  '\0';
    SetSgError("Internal Error: RTMx2XYZ(): BufferXYZ too small");
    return NULL;
  }

  return BufferXYZ;
}


void PrintMapleRTMx(const T_RTMx *RTMx, int FacRo, int FacTr,
                    const char *Label, FILE *fpout)
{
  int         i, j, nt;
  const int   *r, *t;
  const char  *ff;


  if (Label)
    fprintf(fpout, "%s", Label);

  fprintf(fpout, " := matrix(4,4, [");

  r = RTMx->s.R;
  t = RTMx->s.T;

  for (i = 0; i < 3; i++, t++)
  {
    putc(' ', fpout);

    for (j = 0; j < 3; j++, r++)
    {
          ff = FormatFraction(*r, FacRo, 0, NULL, 0);
      if (ff == NULL)
        return;

      fprintf(fpout, "%s,", ff);
    }

        nt = iModPositive(*t, FacTr);
    if (nt >  FacTr / 2)
        nt -= FacTr;

    ff = FormatFraction(nt, FacTr, 0, NULL, 0);/*indentation of this line fixed by NCrystal developers*/
    if (ff == NULL)
      return;

    fprintf(fpout, "%s,", ff);
  }

  fprintf(fpout, " 0,0,0,1]);\n");
}


static void PrintSeitzMx(const T_RTMx *SMx, FILE *fpout)
{
  int         i, nt;
  const char  *ff;
  const int   *r, *t;


  r = SMx->s.R;
  t = SMx->s.T;

  for (i = 0; i < 3; i++)
  {
    fprintf(fpout, " %2d", *r++);
    fprintf(fpout, " %2d", *r++);
    fprintf(fpout, " %2d", *r++);

        nt = iModPositive(*t++, STBF);
    if (nt >  STBF / 2)
        nt -= STBF;

    ff = FormatFraction(nt, STBF, 0, NULL, 0);/*indentation of this line fixed by NCrystal developers*/
    if (ff == NULL)
      return;

    fprintf(fpout, " %6s\n", ff);
  }

  putc('\n', fpout);
}


void ListSgInfo(const T_SgInfo *SgInfo, int F_XYZ, int F_Verbose, FILE *fpout)
{
  int           iList, i_si_v;
  char          buf[16];    /* Changed from buf[8] by NCrystal developers to avoid bad snprintf usage */
  const char    *xyz;
  const T_RTMx  *lsmx;
  T_RotMxInfo   *rmxi, RotMxInfo;


  iList = PG_Index(SgInfo->PointGroup);

  fprintf(fpout, "Point Group  %s\n", PG_Names[iList]);
  fprintf(fpout, "Laue  Group  %s\n",
    PG_Names[PG_Index(LG_Code_of_PG_Index[iList])]);

  fprintf(fpout, "%s\n", XS_Name[SgInfo->XtalSystem]);

  if (SgInfo->UniqueRefAxis != 0 || SgInfo->UniqueDirCode != 0)
  {
    fprintf(fpout, "Unique Axis  ");
    if (SgInfo->UniqueRefAxis != 0 && SgInfo->UniqueRefAxis != 'o')
      fprintf(fpout, "%c", SgInfo->UniqueRefAxis);
    if (SgInfo->UniqueDirCode != 0 && SgInfo->UniqueDirCode != '=')
      fprintf(fpout, "%c", SgInfo->UniqueDirCode);
    fprintf(fpout, "\n");
  }

  if (SgInfo->ExtraInfo != EI_Unknown)
    fprintf(fpout, "%s\n", EI_Name[SgInfo->ExtraInfo]);

  if (SgInfo->InversionOffOrigin)
    fprintf(fpout, "Note: Inversion operation off origin\n");

  putc('\n', fpout);

  fprintf(fpout, "Order   %3d\n", SgInfo->OrderL);
  fprintf(fpout, "Order P %3d\n", SgInfo->OrderP);
  putc('\n', fpout);

  if (SgInfo->n_si_Vector >= 0)
  {
    fprintf(fpout, "s.i.Vector  Modulus\n");
    for (i_si_v = 0; i_si_v < SgInfo->n_si_Vector; i_si_v++)
      fprintf(fpout, " %2d %2d %2d   %d\n",
        SgInfo->si_Vector[i_si_v * 3 + 0],
        SgInfo->si_Vector[i_si_v * 3 + 1],
        SgInfo->si_Vector[i_si_v * 3 + 2],
        SgInfo->si_Modulus[i_si_v]);
    putc('\n', fpout);
  }

  if (F_XYZ || F_Verbose)
  {
    fprintf(fpout, "#List   %3d\n", SgInfo->nList);
    putc('\n', fpout);

    lsmx = SgInfo->ListSeitzMx;
    rmxi = SgInfo->ListRotMxInfo;

    if (rmxi == NULL) rmxi = &RotMxInfo;

    for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
    {
      if (rmxi == &RotMxInfo)
      {
        if (GetRotMxInfo(lsmx->s.R, &RotMxInfo) == 0) {
          SetSgError("Error: Illegal SeitzMx in list");
          return;
        }
      }

      if (F_Verbose)
      {
        sprintf(buf, "(%d)", iList + 1);
        fprintf(fpout, "%-4s", buf);

        fprintf(fpout, "  %2d", rmxi->Order);
        if (rmxi->Inverse) fprintf(fpout, "^-1");
        else               fprintf(fpout, "   ");

        fprintf(fpout, " [%2d %2d %2d]",
                        rmxi->EigenVector[0],
                        rmxi->EigenVector[1],
                        rmxi->EigenVector[2]);

        if (rmxi->RefAxis) fprintf(fpout, " '%c'", rmxi->RefAxis);
        else               fprintf(fpout, "    ");
        if (rmxi->DirCode) fprintf(fpout, " '%c'", rmxi->DirCode);
        else               fprintf(fpout, "    ");

        fprintf(fpout, "    ");
      }

          xyz = RTMx2XYZ(lsmx, 1, STBF, 0, 0, 1, ", ", NULL, 0);
      if (xyz)
        fprintf(fpout, "%s", xyz);

      putc('\n', fpout);

      if (xyz == NULL)
        return;

      if (F_Verbose)
        PrintSeitzMx(lsmx, fpout);

      if (rmxi != &RotMxInfo) rmxi++;
    }

    if (iList && F_Verbose == 0)
      putc('\n', fpout);
  }
}
/*
  Space Group Info's (c) 1994-96 Ralf W. Grosse-Kunstleve
 */



#ifdef APP_INCLUDE
#endif

#ifndef AppMalloc
#define AppMalloc(ptr, n) (ptr) = malloc((n) * sizeof (*(ptr)))
#endif
#ifndef AppFree
#define AppFree(ptr, n) free(ptr)
#endif




/* Non elegant way to get s.i. vectors and moduli:
     1. Build field with legal reference points marked (TestField)
     2. Go through list of possible s.i. vects and mods:
        Verify with TestField
 */


void MarkLegalOrigins(const T_SgInfo *SgInfo, int *TestField)
{
  int           O[3], V[3], lx, ly, lz, mx, my, mz, i;
  int           IsFine, iList, iLoopInv, nLoopInv;
  int           BufMx[9];
  const T_RTMx  *lsmx;
  int           nTrV, iTrV;
  const int     *TrV;


  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;

  switch (SgInfo->LatticeInfo->Code)
  {
    default:
    case 'P': lx = ly = lz = 12;    break;
    case 'A': lx = ly = 12; lz = 6; break;
    case 'B': ly = lz = 12; lx = 6; break;
    case 'C': lz = lx = 12; ly = 6; break;
    case 'I': lx = ly = 12; lz = 6; break;
    case 'R': lx = ly = 12; lz = 4; break;
    case 'S': lz = lx = 12; ly = 4; break;
    case 'T': ly = lz = 12; lx = 4; break;
    case 'F': lx = 12; ly = lz = 6; break;
  }

  for (O[0] = 0; O[0] < 12; O[0]++)
  for (O[1] = 0; O[1] < 12; O[1]++)
  for (O[2] = 0; O[2] < 12; O[2]++)
  {
    IsFine = 1;

    for (iList = 0; IsFine && iList < SgInfo->nList; iList++)
    {
      lsmx = &SgInfo->ListSeitzMx[iList];

      for (iLoopInv = 0; IsFine && iLoopInv < nLoopInv; iLoopInv++)
      {
        if (iLoopInv == 0)
          for (i = 0; i < 9; i++)
          {
            if (i % 4) BufMx[i] =  lsmx->s.R[i];
            else       BufMx[i] =  lsmx->s.R[i] - 1;
          }
        else
          for (i = 0; i < 9; i++)
          {
            if (i % 4) BufMx[i] = -lsmx->s.R[i];
            else       BufMx[i] = -lsmx->s.R[i] - 1;
          }

        RotMx_t_Vector(V, BufMx, O, 12);

        TrV = SgInfo->LatticeInfo->TrVector;

        for (iTrV = 0; iTrV < nTrV; iTrV++)
        {
          mx = (V[0] * (STBF / 12) + *TrV++) % STBF;
          my = (V[1] * (STBF / 12) + *TrV++) % STBF;
          mz = (V[2] * (STBF / 12) + *TrV++) % STBF;

          if (mx == 0 && my == 0 && mz == 0)
            break;
        }

        if (iTrV == nTrV) IsFine = 0;
      }
    }

    if (! (O[0] < lx && O[1] < ly && O[2] < lz))
      IsFine = -IsFine;

    *TestField++ = IsFine;

#if DEBUG_MarkLegalOrigins
    if      (IsFine ==  1) putc(' ', stdout);
    else if (IsFine == -1) putc('#', stdout);
    if (IsFine != 0)
      fprintf(stdout, " %2d %2d %2d\n", O[0], O[1], O[2]);
#endif
  }
}


#define IsArbitraryShift(iShift) \
  (    (iShift) == 1 || (iShift) ==  5 \
    || (iShift) == 7 || (iShift) == 11)


int Verify_si(int h, int k, int l, const int *TestField)
{
  int    O[3], TH;


  for (O[0] = 0; O[0] < 12; O[0]++)
  for (O[1] = 0; O[1] < 12; O[1]++)
  for (O[2] = 0; O[2] < 12; O[2]++)
  {
    if (*TestField++)
    {
          TH = h * O[0] + k * O[1] + l * O[2];
          TH %= 12;
      if (TH) return 0;

      if (IsArbitraryShift(O[0])) TH += h;
      if (IsArbitraryShift(O[1])) TH += k;
      if (IsArbitraryShift(O[2])) TH += l;
      if (TH) return 0;
    }
  }

  return 1;
}


int Is_si(const T_SgInfo *SgInfo, int h, int k, int l)
{
  int        i_si_v, u;
  const int  *si_v, *si_m;


  si_v = SgInfo->si_Vector;
  si_m = SgInfo->si_Modulus;

  for (i_si_v = 0; i_si_v < SgInfo->n_si_Vector; i_si_v++)
  {
    u =  *si_v++ * h;
    u += *si_v++ * k;
    u += *si_v++ * l;

    if (*si_m) {
      if (u % (*si_m)) return 0; }
    else {
      if (u)           return 0; }

    si_m++;
  }

  return 1;
}


int Set_si(T_SgInfo *SgInfo)
{
  static const int TabTrial_si[] =
    {
      0,

      1,   0,  2, -1,  4,  /* I -4 */
      1,   2, -1,  0,  4,
      1,  -1,  0,  2,  4,

      1,   2,  4,  3,  6,  /* P 3 2 */
      1,   4,  3,  2,  6,
      1,   3,  2,  4,  6,

      1,   1,  1,  1,  4,
      1,   1,  1,  1,  2,
      1,   1,  1,  1,  0,

      1,   0,  0,  1,  2,
      1,   0,  1,  0,  2,
      1,   1,  0,  0,  2,

      1,   0,  0,  1,  0,
      1,   0,  1,  0,  0,
      1,   1,  0,  0,  0,

      2,   1, -1,  0,  3,
           0,  0,  1,  0,
      2,  -1,  0,  1,  3,
           0,  1,  0,  0,
      2,   0,  1, -1,  3,
           1,  0,  0,  0,

      2,   0,  1,  1,  4,  /* F 2x */
           1,  0,  0,  0,
      2,   1,  0,  1,  4,  /* F 2y */
           0,  1,  0,  0,
      2,   1,  1,  0,  4,  /* F 2z */
           0,  0,  1,  0,

      2,   1,  0,  0,  2,
           0,  0,  1,  2,
      2,   0,  1,  0,  2,
           0,  0,  1,  2,
      2,   1,  0,  0,  2,
           0,  1,  0,  2,

      2,   1,  1,  0,  2,
           0,  0,  1,  2,
      2,   1,  0,  1,  2,
           0,  1,  0,  2,
      2,   0,  1,  1,  2,
           1,  0,  0,  2,

      2,   1,  0,  0,  2,
           0,  0,  1,  0,
      2,   0,  1,  0,  2,
           0,  0,  1,  0,
      2,   1,  0,  0,  2,
           0,  1,  0,  0,

      2,   1,  0,  0,  0,
           0,  0,  1,  2,
      2,   0,  1,  0,  0,
           0,  0,  1,  2,
      2,   1,  0,  0,  0,
           0,  1,  0,  2,

      2,   1,  1,  0,  2,
           0,  0,  1,  0,
      2,   1,  0,  1,  2,
           0,  1,  0,  0,
      2,   0,  1,  1,  2,
           1,  0,  0,  0,

      2,   1,  0,  0,  0,
           0,  0,  1,  0,
      2,   0,  1,  0,  0,
           0,  0,  1,  0,
      2,   1,  0,  0,  0,
           0,  1,  0,  0,

      3,   1,  0,  0,  2,
           0,  1,  0,  2,
           0,  0,  1,  2,

      3,   1,  0,  0,  0,
           0,  1,  0,  2,
           0,  0,  1,  2,

      3,   1,  0,  0,  2,
           0,  1,  0,  0,
           0,  0,  1,  2,

      3,   1,  0,  0,  2,
           0,  1,  0,  2,
           0,  0,  1,  0,

      3,   1,  0,  0,  2,
           0,  1,  0,  0,
           0,  0,  1,  0,

      3,   1,  0,  0,  0,
           0,  1,  0,  2,
           0,  0,  1,  0,

      3,   1,  0,  0,  0,
           0,  1,  0,  0,
           0,  0,  1,  2,

      3,   1,  0,  0,  0,
           0,  1,  0,  0,
           0,  0,  1,  0,

      3,  -1,  0,  0,  2,  /* -A 1 */
           0, -1,  1,  4,
           0,  1,  1,  4,

      3,  -1,  0,  1,  4,  /* -B 1 */
           0, -1,  0,  2,
           1,  0,  1,  4,

      3,   1,  1,  0,  4,  /* -C 1 */
           1, -1,  0,  4,
           0,  0, -1,  2,

      3,  -1,  1,  1,  4,  /* -I 1 */
           1, -1,  1,  4,
           1,  1, -1,  4,

      3,   0,  1,  1,  4,  /* -F 1 */
           1,  0,  1,  4,
           1,  1,  0,  4,

      3,  -1,  0,  0,  0,  /* A 2x */
           0, -1,  1,  4,
           0,  1,  1,  4,

      3,  -1,  0,  1,  4,  /* B 2y */
           0, -1,  0,  0,
           1,  0,  1,  4,

      3,   1,  1,  0,  4,  /* C 2z */
           1, -1,  0,  4,
           0,  0, -1,  0,

      -1
    };

  int        h, k, l, iList;
  int        Maxh, Maxk, Maxl;
  int        Minh, Mink, Minl;
  int        nTestField, *TestField;
  int        nProperty, *Property, *pp;
  int        IsFine, would_be, is;
  int        i_si, *si_v;
  const int  *trial_si;


  SgInfo->n_si_Vector = -1;

                       nTestField = 12 * 12 * 12;
  AppMalloc(TestField, nTestField);
  if (TestField == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  MarkLegalOrigins(SgInfo, TestField);

  Maxh = Maxk = Maxl = 7;
  SetListMin_hkl(SgInfo, Maxk, Maxl, &Minh, &Mink, &Minl);

  nProperty =   (Maxh - Minh + 1)
              * (Maxk - Mink + 1)
              * (Maxl - Minl + 1);
  AppMalloc(Property, nProperty);
  if (Property == NULL) {
    SetSgError("Not enough core");
    AppFree(TestField, nTestField);
    return -1;
  }

  pp = Property;
  for (h = Minh; h <= Maxh; h++)
  for (k = Mink; k <= Maxk; k++)
  for (l = Minl; l <= Maxl; l++)
  {
    iList = IsSysAbsent_hkl(SgInfo, h, k, l, NULL);
    if (SgError != NULL)
    {
      AppFree(Property, nProperty);
      AppFree(TestField, nTestField);
      return -1;
    }

    if (iList == 0)
      *pp++ = Verify_si(h, k, l, TestField);
    else
      *pp++ = -1;
  }

  trial_si = TabTrial_si;
  while (*trial_si >= 0)
  {
    SgInfo->n_si_Vector = *trial_si++;
    si_v = SgInfo->si_Vector;
    for (i_si = 0; i_si < SgInfo->n_si_Vector; i_si++)
    {
      *si_v++ = *trial_si++;
      *si_v++ = *trial_si++;
      *si_v++ = *trial_si++;
      SgInfo->si_Modulus[i_si] = *trial_si++;
    }

    IsFine = 1;

    pp = Property;
    for (h = Minh; IsFine && h <= Maxh; h++)
    for (k = Mink; IsFine && k <= Maxk; k++)
    for (l = Minl; IsFine && l <= Maxl; l++)
    {
      is = *pp++;

      if (is >= 0)
      {
        would_be = Is_si(SgInfo, h, k, l);
        if (is != would_be)
          IsFine = 0;
      }
    }

    if (IsFine)
    {
      AppFree(Property, nProperty);
      AppFree(TestField, nTestField);
      return 0;
    }
  }

  SgInfo->n_si_Vector = -1;
  SetSgError("Internal Error: Can't determine s.i. vectors and moduli");

  AppFree(Property, nProperty);
  AppFree(TestField, nTestField);

  return -1;
}


void Set_uvw(const T_SgInfo *SgInfo, int h, int k, int l, int *uvw)
{
  int        i_si_v, u;
  const int  *si_v, *si_m;


  si_v = SgInfo->si_Vector;
  si_m = SgInfo->si_Modulus;

  for (i_si_v = 0; i_si_v < SgInfo->n_si_Vector; i_si_v++)
  {
    u =  *si_v++ * h;
    u += *si_v++ * k;
    u += *si_v++ * l;

    if (*si_m) u %= (*si_m);
    si_m++;

    uvw[i_si_v] = u;
  }
}

} // nxs namespace
} // NCrystal namespace
