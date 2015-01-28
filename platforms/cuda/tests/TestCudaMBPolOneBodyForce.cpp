/* -------------------------------------------------------------------------- *
 *                                   OpenMMMBPol                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace  OpenMM;
using namespace MBPolPlugin;

extern "C" void registerMBPolCudaKernelFactories();

const double TOL = 1e-4;
#define PI_M               3.141592653589
#define RADIAN            57.29577951308
const double DegreesToRadians = PI_M/180.0;

/* ---------------------------------------------------------------------------------------

   Compute cross product of two 3-vectors and place in 3rd vector

   vectorZ = vectorX x vectorY

   @param vectorX             x-vector
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
     
static void crossProductVector3( double* vectorX, double* vectorY, double* vectorZ ){

    vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
    vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
    vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];

    return;
}

static double dotVector3( double* vectorX, double* vectorY ){
    return vectorX[0]*vectorY[0] + vectorX[1]*vectorY[1] + vectorX[2]*vectorY[2];
}


static void computeMBPolStretchBendForce(int bondIndex,  std::vector<Vec3>& positions, MBPolStretchBendForce& mbpolStretchBendForce,
                                          std::vector<Vec3>& forces, double* energy, FILE* log ) {

    int particle1, particle2, particle3;
    double abBondLength, cbBondLength, angleStretchBend, kStretchBend;

    mbpolStretchBendForce.getStretchBendParameters(bondIndex, particle1, particle2, particle3, abBondLength, cbBondLength, angleStretchBend, kStretchBend);
    angleStretchBend *= RADIAN;
#ifdef MBPOL_DEBUG
    if( log ){
        (void) fprintf( log, "computeMBPolStretchBendForce: bond %d [%d %d %d] ab=%10.3e cb=%10.3e angle=%10.3e k=%10.3e\n", 
                             bondIndex, particle1, particle2, particle3, abBondLength, cbBondLength, angleStretchBend, kStretchBend );
        (void) fflush( log );
    }
#endif

    enum { A, B, C, LastAtomIndex };
    enum { AB, CB, CBxAB, ABxP, CBxP, LastDeltaIndex };
 
    // ---------------------------------------------------------------------------------------
 
    // get deltaR between various combinations of the 3 atoms
    // and various intermediate terms
 
    double deltaR[LastDeltaIndex][3];
    double rAB2 = 0.0;
    double rCB2 = 0.0;
    for( int ii = 0; ii < 3; ii++ ){
         deltaR[AB][ii]  = positions[particle1][ii] - positions[particle2][ii];
         rAB2           += deltaR[AB][ii]*deltaR[AB][ii];

         deltaR[CB][ii]  = positions[particle3][ii] - positions[particle2][ii];
         rCB2           += deltaR[CB][ii]*deltaR[CB][ii];
    }
    double rAB   = sqrt( rAB2 );
    double rCB   = sqrt( rCB2 );

    crossProductVector3( deltaR[CB], deltaR[AB], deltaR[CBxAB] );
    double  rP   = dotVector3( deltaR[CBxAB], deltaR[CBxAB] );
            rP   = sqrt( rP );
 
    if( rP <= 0.0 ){
       return;
    }
    double dot    = dotVector3( deltaR[CB], deltaR[AB] );
    double cosine = dot/(rAB*rCB);
 
    double angle;
    if( cosine >= 1.0 ){
       angle = 0.0;
    } else if( cosine <= -1.0 ){
       angle = PI_M;
    } else {
       angle = RADIAN*acos(cosine);
    }
 
    double termA = -RADIAN/(rAB2*rP);
    double termC =  RADIAN/(rCB2*rP);
 
    // P = CBxAB
 
    crossProductVector3( deltaR[AB], deltaR[CBxAB], deltaR[ABxP] );
    crossProductVector3( deltaR[CB], deltaR[CBxAB], deltaR[CBxP] );
    for( int ii = 0; ii < 3; ii++ ){
       deltaR[ABxP][ii] *= termA;
       deltaR[CBxP][ii] *= termC;
    }
 
    double dr    = rAB - abBondLength + rCB - cbBondLength;
 
    termA        = 1.0/rAB;
    termC        = 1.0/rCB;
 
    double term  = kStretchBend;
 
    // ---------------------------------------------------------------------------------------
 
    // forces
 
    // calculate forces for atoms a, b, c
    // the force for b is then -( a + c)
 
    double subForce[LastAtomIndex][3];
    double dt = angle - angleStretchBend;
    for( int jj = 0; jj < 3; jj++ ){
        subForce[A][jj] = term*(dt*termA*deltaR[AB][jj] + dr*deltaR[ABxP][jj] );
        subForce[C][jj] = term*(dt*termC*deltaR[CB][jj] + dr*deltaR[CBxP][jj] );
        subForce[B][jj] = -( subForce[A][jj] + subForce[C][jj] );
    }
 
    // ---------------------------------------------------------------------------------------
 
    // accumulate forces and energy
 
    forces[particle1][0]       -= subForce[0][0];
    forces[particle1][1]       -= subForce[0][1];
    forces[particle1][2]       -= subForce[0][2];

    forces[particle2][0]       -= subForce[1][0];
    forces[particle2][1]       -= subForce[1][1];
    forces[particle2][2]       -= subForce[1][2];

    forces[particle3][0]       -= subForce[2][0];
    forces[particle3][1]       -= subForce[2][1];
    forces[particle3][2]       -= subForce[2][2];

    *energy                    += term*dt*dr;

#ifdef MBPOL_DEBUG
    if( log ){
        (void) fprintf( log, "computeMBPolStretchBendForce: angle=%10.3e dt=%10.3e dr=%10.3e\n", angle, dt, dr ); 
        (void) fflush( log );
    }
#endif

    return;
}
 
static void computeMBPolStretchBendForces( Context& context, MBPolStretchBendForce& mbpolStretchBendForce,
                                          std::vector<Vec3>& expectedForces, double* expectedEnergy, FILE* log ) {

    // get positions and zero forces

    State state                 = context.getState(State::Positions);
    std::vector<Vec3> positions = state.getPositions();
    expectedForces.resize( positions.size() );
    
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii][0] = expectedForces[ii][1] = expectedForces[ii][2] = 0.0;
    }

    // calculates forces/energy

    *expectedEnergy = 0.0;
    for( int ii = 0; ii < mbpolStretchBendForce.getNumStretchBends(); ii++ ){
        computeMBPolStretchBendForce(ii, positions, mbpolStretchBendForce, expectedForces, expectedEnergy, log );
    }

#ifdef MBPOL_DEBUG
    if( log ){
        (void) fprintf( log, "computeMBPolStretchBendForces: expected energy=%14.7e\n", *expectedEnergy );
        for( unsigned int ii = 0; ii < positions.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]\n", ii, expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2] );
        }
        (void) fflush( log );
    }
#endif
    return;

}

void compareWithExpectedForceAndEnergy( Context& context, MBPolStretchBendForce& mbpolStretchBendForce,
                                        double tolerance, const std::string& idString, FILE* log) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeMBPolStretchBendForces( context, mbpolStretchBendForce, expectedForces, &expectedEnergy, log );
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();

#ifdef MBPOL_DEBUG
    if( log ){
        (void) fprintf( log, "computeMBPolStretchBendForces: expected energy=%14.7e %14.7e\n", expectedEnergy, state.getPotentialEnergy() );
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%14.7e %14.7e %14.7e]   [%14.7e %14.7e %14.7e]\n", ii,
                            expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2], forces[ii][0], forces[ii][1], forces[ii][2] );
        }
        (void) fflush( log );
    }
#endif

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
    }
    ASSERT_EQUAL_TOL( expectedEnergy, state.getPotentialEnergy(), tolerance );
}

void testOneStretchBend( FILE* log ) {

    System system;
    int numberOfParticles = 3;
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    MBPolStretchBendForce* mbpolStretchBendForce = new MBPolStretchBendForce();

    double abLength         = 0.144800000E+01;
    double cbLength         = 0.101500000E+01;
    double angleStretchBend = 0.108500000E+03*DegreesToRadians;
    //double kStretchBend     = 0.750491578E-01;
    double kStretchBend     = 1.0;

    mbpolStretchBendForce->addStretchBend(0, 1, 2, abLength, cbLength, angleStretchBend, kStretchBend );

    system.addForce(mbpolStretchBendForce);
    Context context(system, integrator, Platform::getPlatformByName( "CUDA"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0] = Vec3( 0.262660000E+02,  0.254130000E+02,  0.284200000E+01 );
    positions[1] = Vec3( 0.273400000E+02,  0.244300000E+02,  0.261400000E+01 );
    positions[2] = Vec3( 0.269573220E+02,  0.236108860E+02,  0.216376800E+01 );

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy( context, *mbpolStretchBendForce, TOL, "testOneStretchBend", log );
    
    // Try changing the stretch-bend parameters and make sure it's still correct.
    
    mbpolStretchBendForce->setStretchBendParameters(0, 0, 1, 2, 1.1*abLength, 1.2*cbLength, 1.3*angleStretchBend, 1.4*kStretchBend);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareWithExpectedForceAndEnergy( context, *mbpolStretchBendForce, TOL, "testOneStretchBend", log );
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    mbpolStretchBendForce->updateParametersInContext(context);
    compareWithExpectedForceAndEnergy( context, *mbpolStretchBendForce, TOL, "testOneStretchBend", log );
}

int main(int argc, char* argv[]) {
    try {
        std::cout << "TestCudaMBPolStretchBendForce running test..." << std::endl;
        registerMBPolCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("CudaPrecision", std::string(argv[1]));
        FILE* log = NULL;
        testOneStretchBend( log );
    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}



