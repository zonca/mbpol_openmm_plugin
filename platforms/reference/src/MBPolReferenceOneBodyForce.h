
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

#ifndef __MBPolReferenceOneBodyForce_H__
#define __MBPolReferenceOneBodyForce_H__

#include "openmm/reference/RealVec.h"
#include <vector>

using OpenMM::RealVec;
// ---------------------------------------------------------------------------------------

class MBPolReferenceOneBodyForce {

public:
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    MBPolReferenceOneBodyForce( ){};
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
          --------------------------------------------------------------------------------------- */
 
    ~MBPolReferenceOneBodyForce( ){};

     /**---------------------------------------------------------------------------------------
     
        Calculate MBPol stretch bend ixns (force and energy)
     
        @param numBonds                number of angles
        @param posData                 particle positions
        @param particle1               particle 1 indices
        @param particle2               particle 2 indices
        @param particle3               particle 3 indices
        @param lengthABParameters      ideal AB bond length 
        @param lengthCBParameters      ideal CB bond length 
        @param angle                   ideal angle 
        @param kQuadratic              force constant
        @param forces                  output force vector
     
        @return total energy

     
        --------------------------------------------------------------------------------------- */

    RealOpenMM calculateForceAndEnergy( int numOneBodys, const std::vector<RealVec>& particlePositions, const std::vector<std::vector<int> >& allParticleIndices,
                                                                           std::vector<RealVec>& forces) const;


private:

    /**---------------------------------------------------------------------------------------
    
       Calculate MBPol stretch bend angle ixn (force and energy)
    
    
       @return energy
    
       --------------------------------------------------------------------------------------- */

    double calculateOneBodyIxn(const OpenMM::RealVec& positionO, const OpenMM::RealVec& positionH1, const OpenMM::RealVec& positionH2,
            OpenMM::RealVec& forceO, OpenMM::RealVec& forceH1, OpenMM::RealVec& forceH2) const;
 
};

// ---------------------------------------------------------------------------------------

#endif // _MBPolReferenceOneBodyForce___
