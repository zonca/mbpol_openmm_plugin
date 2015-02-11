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
#include "openmm/cuda/CudaNonbondedUtilities.h"
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

///////////////////////////////////////////// MBPolTwoBodyForce ////////////////////////////////////

class CudaMBPolTwoBodyForceInfo : public CudaForceInfo {
public:
    CudaMBPolTwoBodyForceInfo(const MBPolTwoBodyForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumMolecules();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        force.getParticleParameters(index, particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        return true;
    }
private:
    const MBPolTwoBodyForce& force;
};

CudaCalcMBPolTwoBodyForceKernel::~CudaCalcMBPolTwoBodyForceKernel() {
    cu.setAsCurrent();
}

void CudaCalcMBPolTwoBodyForceKernel::initialize(const System& system, const MBPolTwoBodyForce& force) {
    cu.setAsCurrent();

    // device array
    particleIndices = CudaArray::create<float4>(cu, cu.getPaddedNumAtoms(), "particleIndices");

    // suffix Vec is used for host arrays
    // FIXME forced to convert to float, otherwise type error in real_shfl
    // how to use ints?
    vector<float4> particleIndicesVec(cu.getPaddedNumAtoms());
    for (int i=0; i <  force.getNumMolecules(); i++) {
        std::vector<int> singleParticleIndices;
        force.getParticleParameters(i, singleParticleIndices );
        particleIndicesVec[i] = make_float4((float) singleParticleIndices[0], (float) singleParticleIndices[1], (float) singleParticleIndices[2], (float) singleParticleIndices[3]);
    }

    particleIndices->upload(particleIndicesVec);

    // a parameter is defined per mulecule
    // particleIndices as a parameter fails with an error on read_shfl
    cu.getNonbondedUtilities().addParameter(CudaNonbondedUtilities::ParameterInfo("particleIndices", "float", 4, sizeof(float4), particleIndices->getDevicePointer()));
    map<string, string> replacements;
    // replacements["PARAMS"] = cu.getNonbondedUtilities().addArgument(particleIndices->getDevicePointer(), "int4");
    
    // using an argument instead
   // posq is already on the device, format is float4 (x, y, z, charge) 
   // so, we can just pass as parameters the indices of the particles as we do
   // in the reference platform
   // after we copy params to the device
   //  replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float2");
// 
   // then we also add another argument with the pointer to posq
   // the cu.getPosq().getDevicePointer()
   // so we can then access the position of all particles on the device
   //
   
    // replacements["POSQ"] = cu.getBondedUtilities().addArgument( cu.getPosq().getDevicePointer(), "float3");

    bool useCutoff = (force.getNonbondedMethod() != MBPolTwoBodyForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() == MBPolTwoBodyForce::CutoffPeriodic);
    vector< vector<int> > exclusions;
    cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, false, force.getCutoff(), exclusions, cu.replaceStrings(CudaMBPolKernelSources::twobodyForce, replacements), force.getForceGroup());
    cu.addForce(new CudaMBPolTwoBodyForceInfo(force));

}

double CudaCalcMBPolTwoBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcMBPolTwoBodyForceKernel::copyParametersToContext(ContextImpl& context, const MBPolTwoBodyForce& force) {
    cu.setAsCurrent();
    throw OpenMMException(" CudaCalcMBPolTwoBodyForceKernel::copyParametersToContext not implemented");
    
    cu.invalidateMolecules();
}
