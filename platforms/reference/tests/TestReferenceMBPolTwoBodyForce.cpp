/* -------------------------------------------------------------------------- *
 *                                   OpenMMMBPol                             *
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
 * This tests the Reference implementation of ReferenceMBPolTwoBodyForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/MBPolTwoBodyForce.h"
#include "MBPolReferenceTwoBodyForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "openmm/reference/RealVec.h"

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};

using namespace  OpenMM;
using namespace MBPolPlugin;

const double TOL = 1e-4;

void testTwoBody( double boxDimension ) {

    std::string testName      = "testMBPol2BodyInteraction";

    System system;
    int numberOfParticles          = 6;
    MBPolTwoBodyForce* mbpolTwoBodyForce = new MBPolTwoBodyForce();
    double cutoff = 10;
    mbpolTwoBodyForce->setCutoff( cutoff );

    unsigned int particlesPerMolecule = 3;

    if( boxDimension > 0.0 ){
        Vec3 a( boxDimension, 0.0, 0.0 );
        Vec3 b( 0.0, boxDimension, 0.0 );
        Vec3 c( 0.0, 0.0, boxDimension );
        system.setDefaultPeriodicBoxVectors( a, b, c );
        mbpolTwoBodyForce->setNonbondedMethod(MBPolTwoBodyForce::CutoffPeriodic);
    } else {
        mbpolTwoBodyForce->setNonbondedMethod(MBPolTwoBodyForce::CutoffNonPeriodic);
    }

    std::vector<int> particleIndices(particlesPerMolecule);
    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        particleIndices[0] = jj;
        particleIndices[1] = jj+1;
        particleIndices[2] = jj+2;
        mbpolTwoBodyForce->addParticle( particleIndices);
    }


    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);
    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy;

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );

    positions[3]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[4]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[5]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    expectedForces[0]     = Vec3( -4.85337479, -4.47836379 ,-20.08989563);
    expectedForces[1]     = Vec3( -0.31239868,  0.52518586 , -1.88893830);
    expectedForces[2]     = Vec3(  0.00886712,  0.73323536 , -1.81715325);
    expectedForces[3]     = Vec3( -0.65181727, -0.72947395 ,  5.88973293);
    expectedForces[4]     = Vec3(  4.82340981,  3.20090213 , 16.49522051);
    expectedForces[5]     = Vec3(  0.98531382,  0.74851439 ,  1.41103374);

    expectedEnergy        = 6.14207815;

    system.addForce(mbpolTwoBodyForce);
    std::string platformName;
    #define AngstromToNm 0.1    
    #define CalToJoule   4.184

    platformName = "Reference";
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        forces[ii] /= CalToJoule*10;
        expectedForces[ii] *= -1; // gradient -> forces
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

void testImageMolecules( bool runTestWithAtomImaging) {

    double boxDimension = 10.;
    RealVec box( boxDimension, boxDimension, boxDimension );

    unsigned int numberOfParticles = 6;
    std::vector<RealVec> particlePositions(numberOfParticles);

    particlePositions[0]             = RealVec( 0., 0., 0. );
    particlePositions[1]             = RealVec( 0., 1., 0. );
    particlePositions[2]             = RealVec( 0., 0., 0. );

    particlePositions[3]             = RealVec( 4.5, 0., 0. );
    particlePositions[4]             = RealVec( 5.5, 0., 0. );
    particlePositions[5]             = RealVec( 4.5, 0., 0. );

    std::vector<std::vector<int> > allParticleIndices(6);

    allParticleIndices[0].push_back(0);
    allParticleIndices[0].push_back(1);
    allParticleIndices[0].push_back(2);

    allParticleIndices[3].push_back(3);
    allParticleIndices[3].push_back(4);
    allParticleIndices[3].push_back(5);

    double imagedPositions[numberOfParticles*2];
    imageMolecules(box, 0, 3, particlePositions, allParticleIndices, imagedPositions);

    if (runTestWithAtomImaging)
    {
        // Check that making periodic images of everything with respect to the first oxygen fails
        int siteI = 0;
        int siteJ = 3;

        // Image the hydrogens of the second molecule with respect to the oxygen of the FIRST MOLECULE instead of SECOND
        for(unsigned int a = 1; a < 3; ++a)
            imageParticles(box, imagedPositions, particlePositions[allParticleIndices[siteJ][a]], imagedPositions + 9 + 3*a);
    }

    RealVec tempPosition;
    for (int i=0; i<numberOfParticles; i++) {
           std::cout << "Position atom " << i << ": " << particlePositions[i] << " A" << std::endl;
           for (int xyz=0; xyz<3; xyz++) {
               tempPosition[xyz] = imagedPositions[3*i + xyz];
           }
           std::cout << "Imaged   atom " << i << ": " << tempPosition << " A" << std::endl;
    }
}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceMBPolTwoBodyForce running test..." << std::endl;

        double boxDimension = 0;
        std::cout << "TestReferenceMBPolTwoBodyForce Cluster" << std::endl;
        testTwoBody( boxDimension );

        std::cout << "TestReferenceMBPolTwoBodyForce  Periodic boundary conditions Huge Box" << std::endl;
        boxDimension = 50;
        testTwoBody( boxDimension);

        bool runTestWithAtomImaging = false;
        testImageMolecules(runTestWithAtomImaging);

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }

    std::cout << "Done" << std::endl;
    return 0;
}
