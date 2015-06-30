/* -------------------------------------------------------------------------- *
 *                                   OpenMMMBPol                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,  *
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
 * This tests the Reference implementation of ReferenceMBPolElectrostaticsForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "MBPolReferenceElectrostaticsForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VirtualSite.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include <fenv.h>


#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};


using namespace  OpenMM;
using namespace MBPolPlugin;
using namespace std;
const double TOL = 1e-4;
const double cal2joule = 4.184;

static void testWater3VirtualSite( ) {

    std::string testName      = "testWater3VirtualSiteChloride";
    std::cout << "Test START: " << testName << std::endl;

    int particlesPerMolecule = 4;
    int numberOfParticles     = particlesPerMolecule*3+1;
    double cutoff             = 0.70;

    std::vector<double> outputElectrostaticsMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;


    // beginning of Electrostatics setup
    MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod = MBPolElectrostaticsForce::NoCutoff;

    System system;
    // box dimensions

    // double boxDimension                               = 1.8643;
    // Vec3 a( boxDimension, 0.0, 0.0 );
    // Vec3 b( 0.0, boxDimension, 0.0 );
    // Vec3 c( 0.0, 0.0, boxDimension );
    // system.setDefaultPeriodicBoxVectors( a, b, c );

    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setNonbondedMethod( nonbondedMethod );
    //mbpolElectrostaticsForce->setPolarizationType( polarizationType );
    //mbpolElectrostaticsForce->setCutoffDistance( cutoff );
    //mbpolElectrostaticsForce->setMutualInducedTargetEpsilon( 1.0e-06 );
    //mbpolElectrostaticsForce->setMutualInducedMaxIterations( 500 );
    //mbpolElectrostaticsForce->setAEwald( 5.4459052e+00 );
    //mbpolElectrostaticsForce->setEwaldErrorTolerance( 1.0e-04 );
    if (particlesPerMolecule == 3)
    	mbpolElectrostaticsForce->setIncludeChargeRedistribution(false);

    double virtualSiteWeightO = 0.573293118;
    double virtualSiteWeightH = 0.213353441;
    for( unsigned int jj = 0; jj < numberOfParticles-1; jj += particlesPerMolecule){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );

        if (particlesPerMolecule > 3) {
        system.addParticle( 0. ); // Virtual Site
        system.setVirtualSite(jj+3, new ThreeParticleAverageSite(jj, jj+1, jj+2,
                                                           virtualSiteWeightO, virtualSiteWeightH,virtualSiteWeightH));
        }

    }
    system.addParticle( 35.5 );

    std::vector<double> thole(5);

    thole[TCC] = 0.4;
    thole[TCD] = 0.4;
    thole[TDD] = 0.055;
    thole[TDDOH]  = 0.626;
    thole[TDDHH] = 0.055;

    for( unsigned int jj = 0; jj < numberOfParticles-1; jj += particlesPerMolecule ){
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                            thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                            thole, 0.000294, 0.000294 );
        if (particlesPerMolecule > 3) {
        	mbpolElectrostaticsForce->addElectrostatics(  0., jj, jj+1, jj+2,
                                                    thole,  0.001310,  0.);
        }
    }
    mbpolElectrostaticsForce->addElectrostatics(  -1., -1,-1,-1,
                                                thole,  0.0053602, 0.0053602);

    system.addForce(mbpolElectrostaticsForce);

    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    int i = 0;
    positions[i]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  ); i++;
    positions[i]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  ); i++;
    positions[i]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  ); i++;

    if (particlesPerMolecule==4){
    	positions[i]             = Vec3( -1.43230412, -0.33360265,  1.64727446 ); i++;
    }

    if (numberOfParticles > 5) {

    positions[i]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  ); i++;
    positions[i]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  ); i++;
    positions[i]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  ); i++;
    if (particlesPerMolecule==4){
    	positions[i]             = Vec3( -1.95661974, -0.48654484, -1.18917052 ); i++;
    }
    positions[i]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  ); i++;
    positions[i]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  ); i++;
    positions[i]            = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  ); i++;
    if (particlesPerMolecule==4){
    	positions[i]            = Vec3( -0.73151769,  1.8136042 , -0.13676332 ); i++;
    }
    }
    positions[i]            = Vec3( 0,0,0 ); i++;

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    std::string platformName;
    platformName = "Reference";
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
    double cal2joule = 4.184;

    double expectedEnergy = (-41.1134)*cal2joule;
    std::cout << "Energy: " << energy/cal2joule << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy/cal2joule << " Kcal/mol "<< std::endl;

    std::vector<Vec3> expectedForces(numberOfParticles);
    i = 0;
    expectedForces[i]         = Vec3( -26.1368, -5.40827, 16.2952  ); i++;
    expectedForces[i]         = Vec3( 22.1692, 25.9079, -62.9639  ); i++;
    expectedForces[i]         = Vec3( -1.5177, 3.04609, -3.8373  ); i++;
    if (particlesPerMolecule > 3)
        i++;
    if (numberOfParticles >5) {

    expectedForces[i]         = Vec3( -16.3635, 4.3361, -19.7297  ); i++;
    expectedForces[i]         = Vec3( 36.4086, 8.18161, 28.0008  ); i++;
    expectedForces[i]         = Vec3( 1.65623, 1.10512, -0.705989 ); i++;
    if (particlesPerMolecule > 3)
        i++;
    expectedForces[i]         = Vec3( -9.97784, 32.2664, -8.60591  ); i++;
    expectedForces[i]         = Vec3( 15.3607, -33.4848, -8.22473  ); i++;
    expectedForces[i]        = Vec3( 2.03685, -13.2406, 9.84  ); i++;
    if (particlesPerMolecule > 3)
        i++;
    }
    expectedForces[i]        = Vec3( -23.6357, -22.7096, 49.9315  ); i++;

    // gradient -> forces
    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            forces[i][j] /= cal2joule*10;
            if ((particlesPerMolecule>3) && ((i+1) % 4 == 0)) { // Set virtual site force to 0
                forces[i][j] = 0;
            }
           }
       }
    std::cout  << std::endl << "Forces:" << std::endl;


    const double eps = 1.0e-4;

    double x_orig;

    std::vector<Vec3> finiteDifferenceForces(numberOfParticles);
    for (int i=0; i<numberOfParticles; i++) {
        finiteDifferenceForces.push_back(Vec3( 0.,  0., 0.  ));
    }
//    for (int i=0; i<numberOfParticles; i++) {
//        for (int xyz=0; xyz<3; xyz++) {
//            x_orig = positions[i][xyz];
//
//            positions[i][xyz] = x_orig + eps;
//            context.setPositions(positions);
//            context.applyConstraints(1e-4); // update position of virtual site
//            state                = context.getState(State::Energy);
//            const double Ep  = state.getPotentialEnergy();
//
//            positions[i][xyz] = x_orig + 2*eps;
//            context.setPositions(positions);
//            context.applyConstraints(1e-4); // update position of virtual site
//            state                = context.getState(State::Energy);
//            const double E2p  = state.getPotentialEnergy();
//
//            positions[i][xyz] = x_orig - eps;
//            context.setPositions(positions);
//            context.applyConstraints(1e-4); // update position of virtual site
//            state                = context.getState(State::Energy);
//            const double Em   = state.getPotentialEnergy();
//
//            positions[i][xyz] = x_orig - 2*eps;
//            context.setPositions(positions);
//            context.applyConstraints(1e-4); // update position of virtual site
//            state                = context.getState(State::Energy);
//            const double E2m   = state.getPotentialEnergy();
//
//            finiteDifferenceForces[i][xyz] = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
//            positions[i][xyz] = x_orig;
//        }
//
//    }
//
//    for (int i=0; i<numberOfParticles; i++) {
//           for (int j=0; j<3; j++) {
//            finiteDifferenceForces[i][j] /= -1*cal2joule*10;
//           }
//
//       }

    for (int i=0; i<numberOfParticles; i++) {
        std::cout << "Force atom " << i << ": " << expectedForces[i] << " Kcal/mol/A <openmm-mbpol precomputed finite differences>" << std::endl;
        std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl<< std::endl;
//        std::cout << "Force atom " << i << ": " << finiteDifferenceForces[i] << " Kcal/mol/A <openmm-mbpol finite differences>" << std::endl
    }

    std::cout << "Comparison of energy and forces with tolerance: " << tolerance << std::endl << std::endl;

    ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    }

    std::cout << "Test Successful: " << testName << std::endl << std::endl;


    return;
}

static void testWater3VirtualSitePMEHugeBox( ) {

    std::string testName      = "testWater3VirtualSitePMEHugeBox";
    std::cout << "Test START: " << testName << std::endl;

    unsigned int particlesPerMolecule = 4;

    int numberOfParticles     = particlesPerMolecule*3+1;
    double cutoff             = 10;

    std::vector<double> outputElectrostaticsMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;


    // beginning of Electrostatics setup
    MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod = MBPolElectrostaticsForce::PME;

    System system;

    double boxDimension                               = 50;
    Vec3 a( boxDimension, 0.0, 0.0 );
    Vec3 b( 0.0, boxDimension, 0.0 );
    Vec3 c( 0.0, 0.0, boxDimension );
    system.setDefaultPeriodicBoxVectors( a, b, c );

    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setNonbondedMethod( nonbondedMethod );
    mbpolElectrostaticsForce->setCutoffDistance( cutoff );

    if (particlesPerMolecule == 3) {
    	mbpolElectrostaticsForce->setIncludeChargeRedistribution(false);
    }
    // disable Ewald by setting alpha to very low value
    mbpolElectrostaticsForce->setAEwald( 1e-15 );

    std::vector<int> pmeGridDimension( 3 );
    int inputPmeGridDimension = 20;
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
    mbpolElectrostaticsForce->setPmeGridDimensions( pmeGridDimension );

    double virtualSiteWeightO = 0.573293118;
    double virtualSiteWeightH = 0.213353441;
    for( unsigned int jj = 0; jj < numberOfParticles-1; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        if (particlesPerMolecule == 4) {
        system.addParticle( 0. ); // Virtual Site
        system.setVirtualSite(jj+3, new ThreeParticleAverageSite(jj, jj+1, jj+2,
                                                           virtualSiteWeightO, virtualSiteWeightH,virtualSiteWeightH));
        }
    }
    system.addParticle( 35.5 );
    // system.addParticle( 0. );

    std::vector<double> thole(5);

    thole[TCC] = 0.4;
    thole[TCD] = 0.4;
    thole[TDD] = 0.055;
    thole[TDDOH]  = 0.626;
    thole[TDDHH] = 0.055;

    for( unsigned int jj = 0; jj < numberOfParticles-1; jj += 4 ){
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                            thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                            thole, 0.000294, 0.000294 );
        if (particlesPerMolecule == 4)
        mbpolElectrostaticsForce->addElectrostatics(  0., jj, jj+1, jj+2,
                                                    thole,  0.001310,  0.);

//        if (particlesPerMolecule == 4) {
//        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
//                                            thole, 0.001310, 0. );
//        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
//                                            thole, 0.000294, 0. );
//        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
//                                            thole, 0.000294, 0. );
//        mbpolElectrostaticsForce->addElectrostatics(  0., jj, jj+1, jj+2,
//                                                    thole,  0.001310,  0.);
//    } else {
//        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+2,
//                                            thole, 0.001310, 0. );
//        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+2,
//                                            thole, 0.000294, 0. );
//        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+1,
//                                            thole, 0.000294, 0. );
//    }
    }
//    mbpolElectrostaticsForce->addElectrostatics(  -1, -1, -1, -1,
//                                                thole,  0.0053602, 0.0053602);
    mbpolElectrostaticsForce->addElectrostatics(  -1, -1,-1,-1,
                                                thole,  0.0053602, 0.0053602);

    system.addForce(mbpolElectrostaticsForce);

    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );
    if (particlesPerMolecule == 4)
    positions[3]             = Vec3( -1.43230412, -0.33360265,  1.64727446 );

    if (numberOfParticles > 5) {
    positions[4]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[5]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[6]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );
    positions[7]             = Vec3( -1.95661974, -0.48654484, -1.18917052 );
    }

    if (numberOfParticles > 9) {
    positions[8]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
    positions[9]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
    positions[10]            = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );
    positions[11]            = Vec3( -0.73151769,  1.8136042 , -0.13676332 );
    }
    positions[numberOfParticles-1]            = Vec3( 2,0,0 );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    std::string platformName;
    platformName = "Reference";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    context.applyConstraints(1e-4); // update position of virtual site

    double tolerance          = 1.0e-02;

//    // test energy and forces
//
    State state                = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces   = state.getForces();
    double energy              = state.getPotentialEnergy();
    double cal2joule = 4.184;

    double expectedEnergy = (-13.4724)*cal2joule;
    std::cout << "Energy: " << energy/cal2joule << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy/cal2joule << " Kcal/mol "<< std::endl;

    std::cout  << std::endl << "AEwald:" << mbpolElectrostaticsForce->getAEwald() << std::endl;
    mbpolElectrostaticsForce->getPmeGridDimensions(pmeGridDimension);
    std::cout  << std::endl << "PmeGridDimensions:" << pmeGridDimension[0] << std::endl;


    std::vector<Vec3> expectedForces(numberOfParticles);

    // forces from finite differences
    unsigned int i = 0;
    expectedForces[i]         = Vec3(-13.5221, -2.12167, -3.5848 ); i++;
    expectedForces[i]         = Vec3(18.7817, 3.59001, -11.592 ); i++;
    expectedForces[i]         = Vec3( -0.079764, 3.07032, -4.39644  ); i++;
    if (particlesPerMolecule > 3)
        i++;
    if (numberOfParticles >5) {

    expectedForces[i]         = Vec3( -9.9527, 2.38737, -7.39074  ); i++;
    expectedForces[i]         = Vec3( 9.19241, 5.01003, 18.8356  ); i++;
    expectedForces[i]         = Vec3(3.27407, 1.228, -0.103192  ); i++;
    if (particlesPerMolecule > 3)
        i++;
    expectedForces[i]         = Vec3( -9.97872, 14.3853, 1.29869 ); i++;
    expectedForces[i]         = Vec3( 2.50169, -9.47539, -1.8169 ); i++;
    expectedForces[i]        = Vec3( 2.8797, -11.3631, 1.09643  ); i++;
    if (particlesPerMolecule > 3)
        i++;
    }
    expectedForces[i]        = Vec3(-3.09631, -6.71087, 7.65339 ); i++;

    // gradient -> forces
    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            forces[i][j] /= cal2joule*10;
            if (particlesPerMolecule == 4) {
            if ((i+1) % 4 == 0) { // Set virtual site force to 0
                forces[i][j] = 0;
            }}
           }
       }

    // ALREADY MINUS
//    for (int i=0; i<numberOfParticles; i++) {
//           for (int j=0; j<3; j++) {
//               expectedForces[i][j] *= -1;
//           }
//       }
    std::cout  << std::endl << "Forces:" << std::endl;


    const double eps = 1.0e-4;

    double x_orig;

    std::vector<Vec3> finiteDifferenceForces(numberOfParticles);
    for (int i=0; i<numberOfParticles; i++) {
        finiteDifferenceForces.push_back(Vec3( 0.,  0., 0.  ));
    }
//    for (int i=0; i<numberOfParticles; i=i+1) {
//        for (int xyz=0; xyz<3; xyz++) {
//            x_orig = positions[i][xyz];
//
//            positions[i][xyz] = x_orig + eps;
//            context.setPositions(positions);
//            context.applyConstraints(1e-4); // update position of virtual site
//            state                = context.getState(State::Energy);
//            const double Ep  = state.getPotentialEnergy();
//
//            positions[i][xyz] = x_orig + 2*eps;
//            context.setPositions(positions);
//            context.applyConstraints(1e-4); // update position of virtual site
//            state                = context.getState(State::Energy);
//            const double E2p  = state.getPotentialEnergy();
//
//            positions[i][xyz] = x_orig - eps;
//            context.setPositions(positions);
//            context.applyConstraints(1e-4); // update position of virtual site
//            state                = context.getState(State::Energy);
//            const double Em   = state.getPotentialEnergy();
//
//            positions[i][xyz] = x_orig - 2*eps;
//            context.setPositions(positions);
//            context.applyConstraints(1e-4); // update position of virtual site
//            state                = context.getState(State::Energy);
//            const double E2m   = state.getPotentialEnergy();
//
//            finiteDifferenceForces[i][xyz] = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
//            positions[i][xyz] = x_orig;
//        }
//
//        std::cout << "Force atom " << i << std::endl;
//
////        std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl;
////        std::cout << "Force atom " << i << ": " << finiteDifferenceForces[i] << " Kcal/mol/A <openmm-mbpol finite differences>" << std::endl << std::endl;
//    }
//
//    for (int i=0; i<numberOfParticles; i++) {
//           for (int j=0; j<3; j++) {
//            finiteDifferenceForces[i][j] /= -1*cal2joule*10;
//           }
//
//       }


    for (int i=0; i<numberOfParticles; i++) {
        std::cout << "Force atom " << i << ": " << expectedForces[i] << " Kcal/mol/A <openmm-mbpol precomputed finite differences>" << std::endl;
        std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl;
//        std::cout << "Force atom " << i << ": " << finiteDifferenceForces[i] << " Kcal/mol/A <openmm-mbpol finite differences>" << std::endl;
        std::cout << std::endl;
    }

//    std::cout << "Comparison of energy and forces with tolerance: " << tolerance << std::endl << std::endl;
//
    ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    }

    std::cout << "Test Successful: " << testName << std::endl << std::endl;


    return;
}


int main( int numberOfArguments, char* argv[] ) {

	feenableexcept(FE_INVALID | FE_OVERFLOW);

    try {
        std::cout << "TestReferenceMBPolElectrostaticsForce running test..." << std::endl;

        testWater3VirtualSite( );
        testWater3VirtualSitePMEHugeBox();

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
