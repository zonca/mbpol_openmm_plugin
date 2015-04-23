/* -------------------------------------------------------------------------- *
 *                                   OpenMMMBPolAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs                                                   *
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
 * This tests the CUDA implementation of MBPolElectrostaticsForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol, testname) {double _norm_ = std::sqrt(expected.dot(expected)); double _scale_ = _norm_ > 1.0 ? _norm_ : 1.0; if ((std::abs((expected[0])-(found[0]))/_scale_ > (tol)) || (std::abs((expected[1])-(found[1]))/_scale_ > (tol)) || (std::abs((expected[2])-(found[2]))/_scale_ > (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define CAL2JOULE 4.184

using namespace OpenMM;
using namespace std;
using namespace MBPolPlugin;

const double TOL = 1e-4;

extern "C" void registerMBPolCudaKernelFactories();

static void testWater3( ) {

    std::string testName      = "testWater3";


    int numberOfParticles     = 9;
    double cutoff             = 0.70;

    std::vector<double> outputElectrostaticsMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;


    // beginning of Electrostatics setup
    MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod = MBPolElectrostaticsForce::NoCutoff;

    System system;
    // box dimensions
    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();
    mbpolElectrostaticsForce->setNonbondedMethod( nonbondedMethod );
    mbpolElectrostaticsForce->setIncludeChargeRedistribution(false);

    unsigned int particlesPerMolecule = 3;

    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
    }

    std::vector<double> zeroDipole(3);
    std::vector<double> zeroQuadrupole(9);
    std::vector<double> thole(5);
    std::fill(zeroDipole.begin(), zeroDipole.end(), 0.);
    std::fill(zeroQuadrupole.begin(), zeroQuadrupole.end(), 0.);
    std::fill(thole.begin(), thole.end(), 0.4);

    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                            thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                            thole, 0.000294, 0.000294 );
    }

    system.addForce(mbpolElectrostaticsForce);

    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );

    positions[3]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[4]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[5]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );

    positions[6]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
    positions[7]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
    positions[8]             = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    std::string platformName;
    platformName = "CUDA";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    context.applyConstraints(1e-4); // update position of virtual site

    double tolerance          = 1.0e-04;

//    // test energy and forces
//
    State state                = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces   = state.getForces();
    double energy              = state.getPotentialEnergy();
//    std::cout << "Forces" << std::endl;


//    std::vector<Vec3> expectedForces(numberOfParticles);
//    expectedForces[0]         = Vec3( -1.029233628e-01,  1.752006876e-01, -2.394228296e-01  );
//    expectedForces[1]         = Vec3(  1.238286503e-01, -9.713944883e-02,  9.278441270e-02  );
//    expectedForces[2]         = Vec3( -1.992936921e-02, -8.084103617e-02,  1.660930712e-01  );
//    expectedForces[3]         = Vec3(  2.181116801e-01,  1.127169979e-01, -1.998507867e-01  );
//    expectedForces[4]         = Vec3( -1.021411513e-01, -6.244910893e-02,  1.595471969e-01  );
//    expectedForces[5]         = Vec3( -1.214347018e-01, -6.329887574e-02,  2.105405984e-02  );
//    expectedForces[6]         = Vec3(  1.708442625e-01,  1.860776100e-01,  2.249030303e-02  );
//    expectedForces[7]         = Vec3( -7.205290616e-02, -7.830256131e-02,  4.942309713e-02  );
//    expectedForces[8]         = Vec3( -9.430310162e-02, -9.196426456e-02, -7.211852443e-02  );
//    for (int i=0; i<numberOfParticles; i++) {
//        for (int j=0; j<3; j++) {
//            expectedForces[i][j] *= CAL2JOULE*10;
//        }
//    }


    //for( unsigned int ii = 0; ii < forces.size(); ii++ ){
    //    ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    //}

    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            forces[i][j] /= CAL2JOULE*10;
           }
       }

    std::cout << "Test start: " << testName << std::endl;

    std::cout  << std::endl << "Forces:" << std::endl;

    for (int i=0; i<numberOfParticles; i++) {
         std::cout << forces[i] << " Kcal/mol/A " << std::endl;
    }
    // Energy elec+ind(kcal/mol): -2.134083549e-02
    double expectedEnergy = -19.6545*CAL2JOULE;
    // ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );
    std::cout << "Energy: " << energy/CAL2JOULE << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy/CAL2JOULE << " Kcal/mol "<< std::endl;
    const double eps = 1.0e-4;

    double x_orig;

    std::vector<Vec3> finiteDifferenceForces(numberOfParticles);
    for (int i=0; i<numberOfParticles; i++) {
        finiteDifferenceForces.push_back(Vec3( 0.,  0., 0.  ));
    }
    for (int i=0; i<numberOfParticles; i++) {
        for (int xyz=0; xyz<3; xyz++) {
            x_orig = positions[i][xyz];

            positions[i][xyz] = x_orig + eps;
            context.setPositions(positions);
            state                = context.getState(State::Energy);
            const double Ep  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig + 2*eps;
            context.setPositions(positions);
            state                = context.getState(State::Energy);
            const double E2p  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - eps;
            context.setPositions(positions);
            state                = context.getState(State::Energy);
            const double Em   = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - 2*eps;
            context.setPositions(positions);
            state                = context.getState(State::Energy);
            const double E2m   = state.getPotentialEnergy();

            finiteDifferenceForces[i][xyz] = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
            positions[i][xyz] = x_orig;
        }

    }

    // Flip sign to convert gradient -> forces
    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            finiteDifferenceForces[i][j] /= -1*CAL2JOULE*10;
           }

       }
    std::cout  << std::endl << "Finite difference forces:" << std::endl;


    for (int i=0; i<numberOfParticles; i++) {
         std::cout << finiteDifferenceForces[i] << " Kcal/mol/A " << std::endl;
    }
    std::cout << "Test END: " << testName << std::endl << std::endl;

}


int main(int argc, char* argv[]) {
    try {
        std::cout << "TestCudaMBPolElectrostaticsForce running test..." << std::endl;
        registerMBPolCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("CudaPrecision", std::string(argv[1]));

        testWater3();

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
