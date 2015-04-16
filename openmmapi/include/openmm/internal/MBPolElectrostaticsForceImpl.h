#ifndef OPENMM_MBPOL_MULTIPOLE_FORCE_IMPL_H_
#define OPENMM_MBPOL_MULTIPOLE_FORCE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                                OpenMMMBPol                                *
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

#include "openmm/internal/ForceImpl.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/Kernel.h"
#include "openmm/Vec3.h"
#include <utility>
#include <string>

namespace MBPolPlugin {

/**
 * This is the internal implementation of MBPolElectrostaticsForce.
 */

class OPENMM_EXPORT_MBPOL MBPolElectrostaticsForceImpl : public ForceImpl {
public:
    MBPolElectrostaticsForceImpl(const MBPolElectrostaticsForce& owner);
    ~MBPolElectrostaticsForceImpl();
    void initialize(ContextImpl& context);
    const MBPolElectrostaticsForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>(); // This force field doesn't define any parameters.
    }
    std::vector<std::string> getKernelNames();

    /**
     * Get the CovalentMap for an atom
     * 
     * @param force                MBPolElectrostaticsForce force reference
     * @param index                the index of the atom for which to set parameters
     * @param minCovalentIndex     minimum covalent index
     * @param maxCovalentIndex     maximum covalent index
     */
    static void getCovalentRange( const MBPolElectrostaticsForce& force, int index,
                                  const std::vector< MBPolElectrostaticsForce::CovalentType>& lists,
                                  int* minCovalentIndex, int* maxCovalentIndex );

    /**
     * Get the covalent degree for the  CovalentEnd lists
     * 
     * @param force                MBPolElectrostaticsForce force reference
     * @param covalentDegree      covalent degrees for the CovalentEnd lists
     */
    static void getCovalentDegree( const MBPolElectrostaticsForce& force, std::vector<int>& covalentDegree );

    void getElectrostaticPotential( ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                    std::vector< double >& outputElectrostaticPotential );

    void updateParametersInContext(ContextImpl& context);
 

private:
    const MBPolElectrostaticsForce& owner;
    Kernel kernel;

    static int CovalentDegrees[MBPolElectrostaticsForce::CovalentEnd];
    static bool initializedCovalentDegrees;
    static const int* getCovalentDegrees( void );
};

} // namespace MBPolPlugin

#endif /*OPENMM_MBPOL_MULTIPOLE_FORCE_IMPL_H_*/
