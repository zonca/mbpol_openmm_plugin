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

void testNeighborListBody( double boxDimension, bool addPositionOffset ) {

    std::string testName      = "testNeighborListBody";

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

    positions[0]             = Vec3( 100, 0, 0  );
    positions[1]             = Vec3( 101, 0, 0  );
    positions[2]             = Vec3( 102, 0, 0  );

    positions[3]             = Vec3( 0, 100, 0  );
    positions[4]             = Vec3( 0, 101, 0  );
    positions[5]             = Vec3( 0, 102, 0  );

    positions[6]             = Vec3( 0, 0, 100 );
    positions[7]             = Vec3( 0, 0, 101 );
    positions[8]             = Vec3( 0, 0, 102  );

//    positions[9]             = Vec3( 5.588472140e-01,  -2.006699172e+00, 1.392786582e-01  );
//	positions[10]             = Vec3( 9.411558180e-01,  -1.541226676e+00,  -6.163293071e-01  );
//	positions[11]             = Vec3( 9.858551734e-01,  -1.567124294e+00, 8.830970941e-01  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    if (addPositionOffset) {
        // move second molecule 1 box dimension in Y direction
        positions[3][1] += boxDimension;
        positions[4][1] += boxDimension;
        positions[5][1] += boxDimension;
    }

    expectedForces[0]     = Vec3(  0.29919011, -0.34960381, -0.16238472 );
    expectedForces[1]     = Vec3(  0.34138467, -0.01255068, -0.00998383 );
    expectedForces[2]     = Vec3( -0.44376649,  0.03687577,  0.54604510 );
    expectedForces[3]     = Vec3( -0.01094164, -0.36171476, -0.05130395 );
    expectedForces[4]     = Vec3(  0.24939202,  1.29382952,  0.22930712 );
    expectedForces[5]     = Vec3( -0.13250943, -0.19313418, -0.34123592 );
    expectedForces[6]     = Vec3(  0.56722869,  0.46036139, -0.39999973 );
    expectedForces[7]     = Vec3( -0.75669111, -0.76132457, -0.29799486 );
    expectedForces[8]     = Vec3( -0.11328682, -0.11273867,  0.48755080 );


    // gradients => forces
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii] *= -1;
    }

    expectedEnergy        = 0.15586446;

    system.addForce(mbpolThreeBodyForce);
    std::string platformName;
    #define AngstromToNm 0.1
    #define CalToJoule   4.184

    platformName = "CUDA";
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        forces[ii][0] /= CalToJoule*10;
        forces[ii][1] /= CalToJoule*10;
        forces[ii][2] /= CalToJoule*10;
    }

    double tolerance = 1.0e-03;


    double energy = state.getPotentialEnergy() / CalToJoule;

    std::cout << "Energy: " << energy << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy << " Kcal/mol "<< std::endl;

    std::cout  << std::endl << "Forces:" << std::endl;

    for (int i=0; i<numberOfParticles; i++) {
           std::cout << "Force atom " << i << ": " << expectedForces[i] << " Kcal/mol/A <mbpol>" << std::endl;
           std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl << std::endl;
       }

       std::cout << "Comparison of energy and forces with tolerance: " << tolerance << std::endl << std::endl;



   ASSERT_EQUAL_TOL( expectedEnergy, energy, tolerance );

   for( unsigned int ii = 0; ii < forces.size(); ii++ ){
       ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
   }
   std::cout << "Test Successful: " << testName << std::endl << std::endl;
}

int main(int argc, char* argv[]) {
    try {
        registerMBPolCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("CudaPrecision", string(argv[1]));
        std::cout << "TestCudaMBPolNeighborList running test..." << std::endl;


        double boxDimension = 0;
        testNeighborListBody( boxDimension, false );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
