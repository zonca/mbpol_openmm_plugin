/* -------------------------------------------------------------------------- *
 *                               OpenMMMBPol                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/MBPolElectrostaticsForceImpl.h"
#include "openmm/mbpolKernels.h"
#include <stdio.h>
#include <cmath>

using namespace std;
using namespace  OpenMM;
using namespace MBPolPlugin;

using std::vector;

bool MBPolElectrostaticsForceImpl::initializedCovalentDegrees = false;
int MBPolElectrostaticsForceImpl::CovalentDegrees[]           = { 1,2,3,4,0,1,2,3};

MBPolElectrostaticsForceImpl::MBPolElectrostaticsForceImpl(const MBPolElectrostaticsForce& owner) : owner(owner) {
}

MBPolElectrostaticsForceImpl::~MBPolElectrostaticsForceImpl() {
}

void MBPolElectrostaticsForceImpl::initialize(ContextImpl& context) {

    const System& system = context.getSystem();
    if (owner.getNumElectrostatics() != system.getNumParticles())
        throw OpenMMException("MBPolElectrostaticsForce must have exactly as many particles as the System it belongs to.");

    // check cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == MBPolElectrostaticsForce::PME) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("MBPolElectrostaticsForce: The cutoff distance cannot be greater than half the periodic box size.");
    }   

    double quadrupoleValidationTolerance = 1.0e-05;
    for( int ii = 0; ii < system.getNumParticles(); ii++ ){

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, dampingFactor, polarity ;
        std::vector<double> molecularDipole;
        std::vector<double> molecularQuadrupole;
        std::vector<double> thole;

        owner.getElectrostaticsParameters( ii, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY,
                                      thole, dampingFactor, polarity );

       // check quadrupole is traceless and symmetric

       double trace = fabs( molecularQuadrupole[0] + molecularQuadrupole[4] + molecularQuadrupole[8] );
       if( trace > quadrupoleValidationTolerance ){
             std::stringstream buffer;
             buffer << "MBPolElectrostaticsForce: qudarupole for particle=" << ii;
             buffer << " has nonzero trace: " << trace << "; MBPOL plugin assumes traceless quadrupole.";
             throw OpenMMException(buffer.str());
       }
       if( fabs( molecularQuadrupole[1] - molecularQuadrupole[3] ) > quadrupoleValidationTolerance  ){
             std::stringstream buffer;
             buffer << "MBPolElectrostaticsForce: XY and YX components of quadrupole for particle=" << ii;
             buffer << "  are not equal: [" << molecularQuadrupole[1] << " " << molecularQuadrupole[3] << "];";
             buffer << " MBPOL plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }

       if( fabs( molecularQuadrupole[2] - molecularQuadrupole[6] ) > quadrupoleValidationTolerance  ){
             std::stringstream buffer;
             buffer << "MBPolElectrostaticsForce: XZ and ZX components of quadrupole for particle=" << ii;
             buffer << "  are not equal: [" << molecularQuadrupole[2] << " " << molecularQuadrupole[6] << "];";
             buffer << " MBPOL plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }

       if( fabs( molecularQuadrupole[5] - molecularQuadrupole[7] ) > quadrupoleValidationTolerance  ){
             std::stringstream buffer;
             buffer << "MBPolElectrostaticsForce: YZ and ZY components of quadrupole for particle=" << ii;
             buffer << "  are not equal: [" << molecularQuadrupole[5] << " " << molecularQuadrupole[7] << "];";
             buffer << " MBPOL plugin assumes symmetric quadrupole tensor.";
             throw OpenMMException(buffer.str());
       }

       // only 'Z-then-X', 'Bisector', Z-Bisect, ThreeFold  currently handled

        if( axisType != MBPolElectrostaticsForce::ZThenX     && axisType != MBPolElectrostaticsForce::Bisector &&
            axisType != MBPolElectrostaticsForce::ZBisect    && axisType != MBPolElectrostaticsForce::ThreeFold &&
            axisType != MBPolElectrostaticsForce::ZOnly      && axisType != MBPolElectrostaticsForce::NoAxisType ) {
             std::stringstream buffer;
             buffer << "MBPolElectrostaticsForce: axis type=" << axisType;
             buffer << " not currently handled - only axisTypes[ ";
             buffer << MBPolElectrostaticsForce::ZThenX   << ", " << MBPolElectrostaticsForce::Bisector  << ", ";
             buffer << MBPolElectrostaticsForce::ZBisect  << ", " << MBPolElectrostaticsForce::ThreeFold << ", ";
             buffer << MBPolElectrostaticsForce::NoAxisType;
             buffer << "] (ZThenX, Bisector, Z-Bisect, ThreeFold, NoAxisType) currently handled .";
             throw OpenMMException(buffer.str());
        }
    }
    kernel = context.getPlatform().createKernel(CalcMBPolElectrostaticsForceKernel::Name(), context);
    kernel.getAs<CalcMBPolElectrostaticsForceKernel>().initialize(context.getSystem(), owner);
}

double MBPolElectrostaticsForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcMBPolElectrostaticsForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> MBPolElectrostaticsForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcMBPolElectrostaticsForceKernel::Name());
    return names;
}

const int* MBPolElectrostaticsForceImpl::getCovalentDegrees( void ) {
    if( !initializedCovalentDegrees ){
        initializedCovalentDegrees                                      = true;
        CovalentDegrees[MBPolElectrostaticsForce::Covalent12]               = 1;
        CovalentDegrees[MBPolElectrostaticsForce::Covalent13]               = 2;
        CovalentDegrees[MBPolElectrostaticsForce::Covalent14]               = 3;
        CovalentDegrees[MBPolElectrostaticsForce::Covalent15]               = 4;
        CovalentDegrees[MBPolElectrostaticsForce::PolarizationCovalent11]   = 0;
        CovalentDegrees[MBPolElectrostaticsForce::PolarizationCovalent12]   = 1;
        CovalentDegrees[MBPolElectrostaticsForce::PolarizationCovalent13]   = 2;
        CovalentDegrees[MBPolElectrostaticsForce::PolarizationCovalent14]   = 3;
    }
    return CovalentDegrees;
}

void MBPolElectrostaticsForceImpl::getCovalentRange( const MBPolElectrostaticsForce& force, int atomIndex, const std::vector<MBPolElectrostaticsForce::CovalentType>& lists,
                                                 int* minCovalentIndex, int* maxCovalentIndex ){

    *minCovalentIndex =  999999999;
    *maxCovalentIndex = -999999999;
    for( unsigned int kk = 0; kk < lists.size(); kk++ ){
        MBPolElectrostaticsForce::CovalentType jj = lists[kk];
        std::vector<int> covalentList;
        force.getCovalentMap( atomIndex, jj, covalentList );
        for( unsigned int ii = 0; ii < covalentList.size(); ii++ ){
            if( *minCovalentIndex > covalentList[ii] ){
               *minCovalentIndex = covalentList[ii];
            }
            if( *maxCovalentIndex < covalentList[ii] ){
               *maxCovalentIndex = covalentList[ii];
            }
        }
    }   
    return;
}

void MBPolElectrostaticsForceImpl::getCovalentDegree( const MBPolElectrostaticsForce& force, std::vector<int>& covalentDegree ){
    covalentDegree.resize( MBPolElectrostaticsForce::CovalentEnd );
    const int* CovalentDegrees = MBPolElectrostaticsForceImpl::getCovalentDegrees();
    for( unsigned int kk = 0; kk < MBPolElectrostaticsForce::CovalentEnd; kk++ ){
        covalentDegree[kk] = CovalentDegrees[kk];
    }   
    return;
}

void MBPolElectrostaticsForceImpl::getElectrostaticPotential( ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                                          std::vector< double >& outputElectrostaticPotential ){
    kernel.getAs<CalcMBPolElectrostaticsForceKernel>().getElectrostaticPotential(context, inputGrid, outputElectrostaticPotential);
}

void MBPolElectrostaticsForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcMBPolElectrostaticsForceKernel>().copyParametersToContext(context, owner);
}
