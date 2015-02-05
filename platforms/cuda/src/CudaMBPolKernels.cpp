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

#include "CudaMBPolKernels.h"
#include "CudaMBPolKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"

using namespace MBPolPlugin;
using namespace OpenMM;
using namespace std;

class CudaMBPolOneBodyForceInfo : public CudaForceInfo {
public:
    CudaMBPolOneBodyForceInfo(const MBPolOneBodyForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumOneBodys();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        force.getOneBodyParameters(index, particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        vector<int> particleIndices;
        force.getOneBodyParameters(group1, particleIndices);
        force.getOneBodyParameters(group2, particleIndices);
        // there are no per-molecule parameters, so groups are always identical
        return true;
    }
private:
    const MBPolOneBodyForce& force;
};

CudaCalcMBPolOneBodyForceKernel::~CudaCalcMBPolOneBodyForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcMBPolOneBodyForceKernel::initialize(const System& system, const MBPolOneBodyForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumOneBodys()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumOneBodys()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms;
    for (int i = 0; i < numBonds; i++) {
        vector<int> particleIndices;
        force.getOneBodyParameters(startIndex+i, particleIndices);
        atoms.push_back(particleIndices);
    }
    cu.getBondedUtilities().addInteraction(atoms, CudaMBPolKernelSources::onebodyForce, force.getForceGroup());
    cu.addForce(new CudaMBPolOneBodyForceInfo(force));
}

double CudaCalcMBPolOneBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcMBPolOneBodyForceKernel::copyParametersToContext(ContextImpl& context, const MBPolOneBodyForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumOneBodys()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumOneBodys()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}
