/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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
/**
 * This tests the Cuda implementation of MBPolOneBodyForce.
 */

/* If you need to test the neighbor list you need to add
 * std::cout << "Number of Neighbor Pairs: " << *numPairs << std::endl;
 * at the end of double CudaCalcMBPolThreeBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy)
 * in CudaMBPolKernels.cpp around line 572
 * */


#include "OpenMMMBPol.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "openmm/reference/RealVec.h"

using namespace MBPolPlugin;
using namespace OpenMM;
using namespace std;

const double DegreesToRadians = 3.14159265/180.0;

extern "C" OPENMM_EXPORT void registerMBPolCudaKernelFactories();

void testNeighborList3( double boxDimension ) {

    std::string testName      = "testNeighborList3";

    System system;
    int numberOfParticles          = 9;   // only will run if there are more than 3 molecules
    MBPolThreeBodyForce* mbpolThreeBodyForce = new MBPolThreeBodyForce();
    double cutoff = 10;
    mbpolThreeBodyForce->setCutoff( cutoff );

    if( boxDimension > 0.0 ){
            Vec3 a( boxDimension, 0.0, 0.0 );
            Vec3 b( 0.0, boxDimension, 0.0 );
            Vec3 c( 0.0, 0.0, boxDimension );
            system.setDefaultPeriodicBoxVectors( a, b, c );
            mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffPeriodic);
        } else {
            mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffNonPeriodic);
        }

    unsigned int particlesPerMolecule = 3;

    std::vector<int> particleIndices(particlesPerMolecule);
    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        particleIndices[0] = jj;
        particleIndices[1] = jj+1;
        particleIndices[2] = jj+2;
        mbpolThreeBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
    }


    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);
    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy;

    positions[0]             = Vec3( 0, 0, 0  );
    positions[1]             = Vec3( 1, 0, 0  );
    positions[2]             = Vec3( 2, 0, 0  );

    positions[3]             = Vec3( 0, 0, 0  );
    positions[4]             = Vec3( 0, 1, 0  );
    positions[5]             = Vec3( 0, 2, 0  );

    positions[6]             = Vec3( 0, 0, 0 );
    positions[7]             = Vec3( 0, 0, 1 );
    positions[8]             = Vec3( 0, 0, 2  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    system.addForce(mbpolThreeBodyForce);
    std::string platformName;
    #define AngstromToNm 0.1
    #define CalToJoule   4.184

    platformName = "CUDA";
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();

    std::cout << "Expected number of Neighbor Pairs: " << 3 << std::endl;
}

void testNeighborList1( double boxDimension ) {

    std::string testName      = "testNeighborList1";

    System system;
    int numberOfParticles          = 9;   // only will run if there are more than 3 molecules
    MBPolThreeBodyForce* mbpolThreeBodyForce = new MBPolThreeBodyForce();
    double cutoff = 10;
    mbpolThreeBodyForce->setCutoff( cutoff );

    if( boxDimension > 0.0 ){
            Vec3 a( boxDimension, 0.0, 0.0 );
            Vec3 b( 0.0, boxDimension, 0.0 );
            Vec3 c( 0.0, 0.0, boxDimension );
            system.setDefaultPeriodicBoxVectors( a, b, c );
            mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffPeriodic);
        } else {
            mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffNonPeriodic);
        }

    unsigned int particlesPerMolecule = 3;

    std::vector<int> particleIndices(particlesPerMolecule);
    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        particleIndices[0] = jj;
        particleIndices[1] = jj+1;
        particleIndices[2] = jj+2;
        mbpolThreeBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
    }


    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);
    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy;

    positions[0]             = Vec3( 0, 0, 0  );
    positions[1]             = Vec3( 1, 0, 0  );
    positions[2]             = Vec3( 2, 0, 0  );

    positions[3]             = Vec3( 0, 0, 0  );
    positions[4]             = Vec3( 0, 1, 0  );
    positions[5]             = Vec3( 0, 2, 0  );

    positions[6]             = Vec3( 0, 0, 100 );
    positions[7]             = Vec3( 0, 0, 101 );
    positions[8]             = Vec3( 0, 0, 102  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    system.addForce(mbpolThreeBodyForce);
    std::string platformName;
    #define AngstromToNm 0.1
    #define CalToJoule   4.184

    platformName = "CUDA";
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();

    std::cout << "Expected number of Neighbor Pairs: " << 1 << std::endl;
}

void testNeighborList0( double boxDimension ) {

    std::string testName      = "testNeighborList0";

    System system;
    int numberOfParticles          = 9;   // only will run if there are more than 3 molecules
    MBPolThreeBodyForce* mbpolThreeBodyForce = new MBPolThreeBodyForce();
    double cutoff = 10;
    mbpolThreeBodyForce->setCutoff( cutoff );

    if( boxDimension > 0.0 ){
            Vec3 a( boxDimension, 0.0, 0.0 );
            Vec3 b( 0.0, boxDimension, 0.0 );
            Vec3 c( 0.0, 0.0, boxDimension );
            system.setDefaultPeriodicBoxVectors( a, b, c );
            mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffPeriodic);
        } else {
            mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffNonPeriodic);
        }

    unsigned int particlesPerMolecule = 3;

    std::vector<int> particleIndices(particlesPerMolecule);
    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        particleIndices[0] = jj;
        particleIndices[1] = jj+1;
        particleIndices[2] = jj+2;
        mbpolThreeBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
    }


    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);
    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy;

    positions[0]             = Vec3( 0, 0, 0  );
    positions[1]             = Vec3( 1, 0, 0  );
    positions[2]             = Vec3( 2, 0, 0  );

    positions[3]             = Vec3( 0, 100, 0  );
    positions[4]             = Vec3( 0, 101, 0  );
    positions[5]             = Vec3( 0, 102, 0  );

    positions[6]             = Vec3( 0, 0, 100 );
    positions[7]             = Vec3( 0, 0, 101 );
    positions[8]             = Vec3( 0, 0, 102  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    system.addForce(mbpolThreeBodyForce);
    std::string platformName;
    #define AngstromToNm 0.1
    #define CalToJoule   4.184

    platformName = "CUDA";
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();

    std::cout << "Expected number of Neighbor Pairs: " << 0 << std::endl;}

int main(int argc, char* argv[]) {
    try {
        registerMBPolCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("CudaPrecision", string(argv[1]));
        std::cout << "TestCudaMBPolNeighborList running test..." << std::endl;


        double boxDimension = 0;
        testNeighborList3( boxDimension );
        testNeighborList1( boxDimension );
        testNeighborList0( boxDimension );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
