#ifndef CUDA_MBPOL_KERNELS_H_
#define CUDA_MBPOL_KERNELS_H_

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

#include "openmm/mbpolKernels.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/cuda/CudaArray.h"

namespace MBPolPlugin {

/**
 * This kernel is invoked by MBPolForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolOneBodyForceKernel : public CalcMBPolOneBodyForceKernel {
public:
    CudaCalcMBPolOneBodyForceKernel(std::string name, const OpenMM::Platform& platform, OpenMM::CudaContext& cu, const OpenMM::System& system) :
            CalcMBPolOneBodyForceKernel(name, platform), hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }
    ~CudaCalcMBPolOneBodyForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolOneBodyForce this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const MBPolOneBodyForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolOneBodyForce to copy the parameters from
     */
    void copyParametersToContext(OpenMM::ContextImpl& context, const MBPolOneBodyForce& force);
private:
    int numBonds;
    bool hasInitializedKernel;
    OpenMM::CudaContext& cu;
    const OpenMM::System& system;
    OpenMM::CudaArray* params;
};

class CudaCalcMBPolTwoBodyForceKernel : public CalcMBPolTwoBodyForceKernel {
public:
    CudaCalcMBPolTwoBodyForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);

    ~CudaCalcMBPolTwoBodyForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolTwoBodyForce this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const MBPolTwoBodyForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolTwoBodyForce to copy the parameters from
     */
    void copyParametersToContext(OpenMM::ContextImpl& context, const MBPolTwoBodyForce& force);
private:
    int numMolecules;
    CudaArray* particleIndices;
    bool hasInitializedKernel;
    OpenMM::CudaContext& cu;
    const OpenMM::System& system;
    OpenMM::CudaArray* params;
    CUfunction computeTwoBodyForceKernel;

};

class CudaCalcMBPolThreeBodyForceKernel : public CalcMBPolThreeBodyForceKernel {
public:
    CudaCalcMBPolThreeBodyForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);

    ~CudaCalcMBPolThreeBodyForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolThreeBodyForce this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const MBPolThreeBodyForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolThreeBodyForce to copy the parameters from
     */
    void copyParametersToContext(OpenMM::ContextImpl& context, const MBPolThreeBodyForce& force);
private:
    int numMolecules;
    CudaArray* particleIndices;
    bool hasInitializedKernel;
    OpenMM::CudaContext& cu;
    const OpenMM::System& system;
    OpenMM::CudaArray* params;
    CUfunction computeThreeBodyForceKernel;



//    ///////////// things added for neighbor list //////////////////////////
	CudaArray* blockCenter;
	CudaArray* blockBoundingBox;
	CudaArray* neighborPairs;
	CudaArray* numNeighborPairs;
	CudaArray* neighborStartIndex;
	CudaArray* numNeighborsForAtom;
	CudaArray* neighbors;
	int maxNeighborPairs;
//    CUfunction copyPairsKernel, startIndicesKernel, neighborsKernel;
//    std::vector<void*> forceArgs, blockBoundsArgs, neighborsArgs, startIndicesArgs, copyPairsArgs;
//    OpenMM::CudaArray* neighbors, *neighborStartIndex;
//    ///////////////////////////////////////////////////////////////////////
};

} // namespace MBPolPlugin

#endif /*CUDA_MBPOL_KERNELS_H_*/
