
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#include "MBPolReferenceForce.h"
#include "MBPolReferenceDispersionForce.h"
#include "MBPolReferenceTwoBodyForce.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <string>
#include <iostream>

using std::pair;
using std::string;
using std::map;
using std::make_pair;

using std::vector;
using OpenMM::RealVec;
using namespace MBPolPlugin;


const double C6_HH = 2.009358600184719e+01; // kcal/mol * A^(-6)
const double C6_OH = 8.349556669872743e+01; // kcal/mol * A^(-6)
const double C6_OO = 2.373212214147944e+02; // kcal/mol * A^(-6)

const double d6_HH =  9.406475169954112e+00; // A^(-1)
const double d6_OH =  9.775202425217957e+00; // A^(-1)
const double d6_OO =  9.295485815062264e+00; // A^(-1)

MBPolReferenceDispersionForce::MBPolReferenceDispersionForce( ) : _nonbondedMethod(NoCutoff), _cutoff(1.0e+10) {

    _periodicBoxDimensions = RealVec( 0.0, 0.0, 0.0 );

    _c6d6Data[make_pair("H", "H")] = make_pair(C6_HH, d6_HH);
    _c6d6Data[make_pair("O", "H")] = make_pair(C6_OH, d6_OH);
    _c6d6Data[make_pair("O", "O")] = make_pair(C6_OO, d6_OO);

    c6d6Datatype :: const_iterator entry =
            _c6d6Data.find(make_pair("O", "H"));
//    if (entry == _c6d6Data.end())
//        entry = _c6d6Data.find(make_pair(allParticleElements[siteI], allParticleElements[siteJ]));
    pair<double, double> c6d6 = entry->second;
    std::cout << c6d6.first << std::endl;

}

MBPolReferenceDispersionForce::NonbondedMethod MBPolReferenceDispersionForce::getNonbondedMethod( void ) const {
    return _nonbondedMethod;
}

void MBPolReferenceDispersionForce::setNonbondedMethod( MBPolReferenceDispersionForce::NonbondedMethod nonbondedMethod ){
    _nonbondedMethod = nonbondedMethod;
}

void MBPolReferenceDispersionForce::setCutoff( double cutoff ){
    _cutoff  = cutoff;
}

double MBPolReferenceDispersionForce::getCutoff( void ) const {
    return _cutoff;
}

void MBPolReferenceDispersionForce::setPeriodicBox( const RealVec& box ){
    _periodicBoxDimensions = box;
}

RealVec MBPolReferenceDispersionForce::getPeriodicBox( void ) const {
    return _periodicBoxDimensions;
}

void MBPolReferenceDispersionForce::setC6d6Data( const c6d6Datatype& c6d6Data) {
    _c6d6Data = c6d6Data;
}



template <int N>
struct Factorial
{
    enum { value = N * Factorial<N - 1>::value };
};

template <>
struct Factorial<0>
{
    enum { value = 1 };
};

double tang_toennies(int n, const double& x)
{
    assert(n >= 0);

    int nn = n;

    double sum = 1.0 + x/nn;
    while (--nn != 0)
        sum = 1.0 + sum*x/nn;

    double tt = 1.0 - sum*std::exp(-x);

    if (std::fabs(tt) < 1.0e-8) {

        double term(1);
        for (nn = n; nn != 0; --nn)
            term *= x/nn;

        sum = 0.0;
        for (nn = n + 1; nn < 1000; ++nn) {
            term *= x/nn;
            sum += term;

            if (std::fabs(term/sum) < 1.0e-8)
                break;
        }

        tt = sum*std::exp(-x);
    }

    return tt;
}

inline double x6(const double& C6, const double& d6,
                 const OpenMM::RealVec& p1, const OpenMM::RealVec& p2,
                 OpenMM::RealVec& g1,       OpenMM::RealVec& g2)
{
    const double dx = (p1[0] - p2[0]);
    const double dy = (p1[1] - p2[1]);
    const double dz = (p1[2] - p2[2]);

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r;
    const double tt6 = tang_toennies(6, d6r);

    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;

    const double e6 = C6*tt6*inv_r6;

    const double if6 = 1.0/Factorial<6>::value;

    const double grd = 6*e6*inv_rsq - C6*std::pow(d6, 7)*if6*std::exp(-d6r)/r;

    g1[0] += dx*grd * cal2joule * -nm_to_A;
    g2[0] -= dx*grd * cal2joule * -nm_to_A;

    g1[1] += dy*grd * cal2joule * -nm_to_A;
    g2[1] -= dy*grd * cal2joule * -nm_to_A;

    g1[2] += dz*grd * cal2joule * -nm_to_A;
    g2[2] -= dz*grd * cal2joule * -nm_to_A;

    return - e6 * cal2joule;
}

RealOpenMM MBPolReferenceDispersionForce::calculatePairIxn( int siteI, int siteJ,
                                                      const std::vector<RealVec>& particlePositions,
                                                      const std::vector<string>& allParticleElements,
                                                      vector<RealVec>& forces ) const {

        if ((std::abs(siteI-siteJ) == 1) and (allParticleElements[siteI] != "O"))
            return 0;
        if ((std::abs(siteI-siteJ) == 2) and (allParticleElements[siteJ] == "O"))
            return 0;

        std::vector<RealVec> allPositions;

        allPositions.push_back(particlePositions[siteI] * nm_to_A);
        allPositions.push_back(particlePositions[siteJ] * nm_to_A);

        if( _nonbondedMethod == CutoffPeriodic )
            imageParticles(_periodicBoxDimensions, allPositions[0], allPositions[1]);

        c6d6Datatype :: const_iterator entry =
                _c6d6Data.find(make_pair(allParticleElements[siteI], allParticleElements[siteJ]));
        if (entry == _c6d6Data.end())
            entry = _c6d6Data.find(make_pair(allParticleElements[siteJ], allParticleElements[siteI]));
        pair<double, double> c6d6 = entry->second;

        std::cout << siteI << ", " << siteJ << " | " << allParticleElements[siteI] <<  allParticleElements[siteJ] << std::endl;
        std::cout << c6d6.first << ", " << c6d6.second << std::endl;
        RealOpenMM energy= x6(c6d6.first, c6d6.second, allPositions[0], allPositions[1], forces[siteI], forces[siteJ]);

        return energy;

}

RealOpenMM MBPolReferenceDispersionForce::calculateForceAndEnergy( int numParticles,
                                                             const vector<RealVec>& particlePositions,
                                                             const std::vector<string>& allParticleElements,
                                                             const NeighborList& neighborList,
                                                             vector<RealVec>& forces ) const {


    RealOpenMM energy = 0.;
    for( unsigned int ii = 0; ii < neighborList.size(); ii++ ){

        OpenMM::AtomPair pair       = neighborList[ii];
        int siteI                   = pair.first;
        int siteJ                   = pair.second;

        energy                     += calculatePairIxn( siteI, siteJ,
                particlePositions, allParticleElements, forces );

    }

    return energy;
}
