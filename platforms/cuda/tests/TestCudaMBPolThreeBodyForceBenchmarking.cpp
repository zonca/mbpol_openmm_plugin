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

#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */


using namespace MBPolPlugin;
using namespace OpenMM;
using namespace std;

const double DegreesToRadians = 3.14159265/180.0;

extern "C" OPENMM_EXPORT void registerMBPolCudaKernelFactories();

void testNMolecules( int N, double boxDimension) {

    System system;
    unsigned int particlesPerMolecule = 3;
    int numberOfParticles = N*particlesPerMolecule;
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

    std::vector<int> particleIndices(particlesPerMolecule);
    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        particleIndices[0] = jj;
        particleIndices[1] = jj+1;
        particleIndices[2] = jj+2;
        mbpolThreeBodyForce->addParticle( particleIndices);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);

    for (int i = 0; i < numberOfParticles; i++) {
        positions[i]            = Vec3( i, i, i);
    }

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
        		positions[i][j] *= 1e-1;
        }
    }
    system.addForce(mbpolThreeBodyForce);
    std::string platformName;
    platformName = "CUDA";
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();


   std::cout << "Test Successful"<< std::endl;

}

int main(int argc, char* argv[]) {
    try {
        registerMBPolCudaKernelFactories();

        std::cout << "TestCudaMBPolThreeBodyForce running test..." << std::endl;
        clock_t t;
        std::string testName;
        double boxDimension;
        int N;

        testName = "test3Molecules";
        boxDimension = 0;
        N = 3;
        std::cout << std::endl << testName << std::endl;
        t = clock();
        testNMolecules( N, boxDimension );
        t = clock() - t;
        std::cout << testName << " elapsed time = " << ((float)t)/CLOCKS_PER_SEC << "s" << std::endl;

//        testName = "test10Molecules";
//        boxDimension = 0;
//        N = 10;
//        std::cout << std::endl << testName << std::endl;
//        t = clock();
//        testNMolecules( N, boxDimension );
//        t = clock() - t;
//        std::cout << testName << " elapsed time = " << ((float)t)/CLOCKS_PER_SEC << "s" << std::endl;

//        testName = "test20Molecules";
//        boxDimension = 0;
//        N = 20;
//        std::cout << std::endl << testName << std::endl;
//        t = clock();
//        testNMolecules( N, boxDimension );
//        t = clock() - t;
//        std::cout << testName << " elapsed time = " << ((float)t)/CLOCKS_PER_SEC << "s" << std::endl;
//
//        testName = "test25Molecules";
//		boxDimension = 0;
//		N = 25;
//		std::cout << std::endl << testName << std::endl;
//		t = clock();
//		testNMolecules( N, boxDimension );
//		t = clock() - t;
//		std::cout << testName << " elapsed time = " << ((float)t)/CLOCKS_PER_SEC << "s" << std::endl;
//
        testName = "test100Molecules";
		boxDimension = 0;
		N = 100;
		std::cout << std::endl << testName << std::endl;
		t = clock();
		testNMolecules( N, boxDimension );
		t = clock() - t;
		std::cout << testName << " elapsed time = " << ((float)t)/CLOCKS_PER_SEC << "s" << std::endl;



    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << std::endl << "Done" << std::endl;
    return 0;
}
