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
#include "CudaKernelSources.h"
#include "openmm/cuda/CudaParameterSet.h"

#include <cuda_runtime.h>

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

CudaCalcMBPolTwoBodyForceKernel::CudaCalcMBPolTwoBodyForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) :
        CalcMBPolTwoBodyForceKernel(name, platform), cu(cu), system(system) {
}


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
    // cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, false, force.getCutoff(), exclusions, cu.replaceStrings(CudaMBPolKernelSources::twobodyForce, replacements), force.getForceGroup());
    // cu.addForce(new CudaMBPolTwoBodyForceInfo(force));
    
    // create an explicit CUDA kernel, this is necessary because we need access to
    // position and forces of all atoms in each molecule
    //
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
    defines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
    defines["THREAD_BLOCK_SIZE"] = cu.intToString(cu.getNonbondedUtilities().getNumForceThreadBlocks());
    //
    // tiles with exclusions setup
    
    int numContexts = cu.getPlatformData().contexts.size();
    // nb.initialize(system);
    // int numExclusionTiles = nb.getExclusionTiles().getSize();
    int numExclusionTiles = 1;

    defines["NUM_TILES_WITH_EXCLUSIONS"] = cu.intToString(numExclusionTiles);
    int startExclusionIndex = cu.getContextIndex()*numExclusionTiles/numContexts;
    int endExclusionIndex = (cu.getContextIndex()+1)*numExclusionTiles/numContexts;
        defines["FIRST_EXCLUSION_TILE"] = cu.intToString(startExclusionIndex);
    defines["LAST_EXCLUSION_TILE"] = cu.intToString(endExclusionIndex);
    // end of tiles with exclusions setup
    //
    if (useCutoff)
        defines["USE_CUTOFF"] = "1";
    double cutoff = force.getCutoff();
    defines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);

    if (usePeriodic)
        defines["USE_PERIODIC"] = "1";

    CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaMBPolKernelSources::multibodyLibrary + CudaMBPolKernelSources::twobodyForcePolynomial + CudaMBPolKernelSources::twobodyForce, defines);
    computeTwoBodyForceKernel = cu.getKernel(module, "computeTwoBodyForce");

    // Add an interaction to the default nonbonded kernel.  This doesn't actually do any calculations.  It's
    // just so that CudaNonbondedUtilities will build the exclusion flags and maintain the neighbor list.

    cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, false, force.getCutoff(), exclusions, "", force.getForceGroup());
    // cu.getNonbondedUtilities().setUsePadding(false);
    cu.addForce(new CudaMBPolTwoBodyForceInfo(force));

}

double CudaCalcMBPolTwoBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();

    int startTileIndex = nb.getStartTileIndex();
    int numTileIndices = nb.getNumTiles();
    unsigned int maxTiles;
    if (nb.getUseCutoff()) {
        maxTiles = nb.getInteractingTiles().getSize();
    }

    void* args[] = {&cu.getForce().getDevicePointer(),
        &cu.getEnergyBuffer().getDevicePointer(),
        &cu.getPosq().getDevicePointer(),
        &nb.getExclusionTiles().getDevicePointer(),
        &startTileIndex,
        &numTileIndices,
        &cu.getNonbondedUtilities().getInteractingTiles().getDevicePointer(),
        &cu.getNonbondedUtilities().getInteractionCount().getDevicePointer(),
        cu.getPeriodicBoxSizePointer(),
        cu.getInvPeriodicBoxSizePointer(),
        &maxTiles,
       // &cu.getNonbondedUtilities().getBlock().getDevicePointer(),
        &cu.getNonbondedUtilities().getInteractingAtoms().getDevicePointer()
    };
    cu.executeKernel(computeTwoBodyForceKernel, args, cu.getPaddedNumAtoms());
    return 0.0;
}

void CudaCalcMBPolTwoBodyForceKernel::copyParametersToContext(ContextImpl& context, const MBPolTwoBodyForce& force) {
    cu.setAsCurrent();
    throw OpenMMException(" CudaCalcMBPolTwoBodyForceKernel::copyParametersToContext not implemented");
    
    cu.invalidateMolecules();
}

///////////////////////////////////////////// MBPolThreeBodyForce ////////////////////////////////////

class CudaMBPolThreeBodyForceInfo : public CudaForceInfo {
public:
    CudaMBPolThreeBodyForceInfo(const MBPolThreeBodyForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumMolecules();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        force.getParticleParameters(index, particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        return false;
    }
    bool areParticlesIdentical(int particle1, int particle2) {
    	return false;
    }
private:
    const MBPolThreeBodyForce& force;
};

CudaCalcMBPolThreeBodyForceKernel::CudaCalcMBPolThreeBodyForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) :
        CalcMBPolThreeBodyForceKernel(name, platform), cu(cu), system(system) {
}


CudaCalcMBPolThreeBodyForceKernel::~CudaCalcMBPolThreeBodyForceKernel() {
    cu.setAsCurrent();
}

void CudaCalcMBPolThreeBodyForceKernel::initialize(const System& system, const MBPolThreeBodyForce& force) {
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

	    bool useCutoff = (force.getNonbondedMethod() != MBPolThreeBodyForce::NoCutoff);
	    bool usePeriodic = (force.getNonbondedMethod() == MBPolThreeBodyForce::CutoffPeriodic);
	    vector< vector<int> > exclusions;
	    // cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, false, force.getCutoff(), exclusions, cu.replaceStrings(CudaMBPolKernelSources::threebodyForce, replacements), force.getForceGroup());
	    // cu.addForce(new CudaMBPolThreeBodyForceInfo(force));

	    // create an explicit CUDA kernel, this is necessary because we need access to
	    // position and forces of all atoms in each molecule
	    //


	    // Build data structures for the neighbor list.
	    int numParticles = force.getNumParticles();
//	    std::cout << "NumParticles = " << numParticles << std::endl;
		if (useCutoff) {
			int numAtomBlocks = cu.getNumAtomBlocks();
			int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
			blockCenter = new CudaArray(cu, numAtomBlocks, 4*elementSize, "blockCenter");
			blockBoundingBox = new CudaArray(cu, numAtomBlocks, 4*elementSize, "blockBoundingBox");
			numNeighborPairs = CudaArray::create<int>(cu, 1, "customManyParticleNumNeighborPairs");
			neighborStartIndex = CudaArray::create<int>(cu, numParticles+1, "customManyParticleNeighborStartIndex");
			numNeighborsForAtom = CudaArray::create<int>(cu, numParticles, "customManyParticleNumNeighborsForAtom");
			//CHECK_RESULT(cuEventCreate(&event, CU_EVENT_DISABLE_TIMING), "Error creating event for CustomManyParticleForce");

			// Select a size for the array that holds the neighbor list.  We have to make a fairly
			// arbitrary guess, but if this turns out to be too small we'll increase it later.

			maxNeighborPairs = 150*numParticles/3;
			neighborPairs = CudaArray::create<int2>(cu, maxNeighborPairs, "customManyParticleNeighborPairs");
			neighbors = CudaArray::create<int>(cu, maxNeighborPairs, "customManyParticleNeighbors");
		}



	    map<string, string> defines;
	    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
	    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
	    defines["NUM_BLOCKS"] = "1";
	    defines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
	    defines["THREAD_BLOCK_SIZE"] = "1";
	    //
	    // tiles with exclusions setup

	    int numContexts = cu.getPlatformData().contexts.size();
	    // nb.initialize(system);
	    // int numExclusionTiles = nb.getExclusionTiles().getSize();
	    int numExclusionTiles = 1;

	    defines["NUM_TILES_WITH_EXCLUSIONS"] = cu.intToString(numExclusionTiles);
	    int startExclusionIndex = cu.getContextIndex()*numExclusionTiles/numContexts;
	    int endExclusionIndex = (cu.getContextIndex()+1)*numExclusionTiles/numContexts;
	        defines["FIRST_EXCLUSION_TILE"] = cu.intToString(startExclusionIndex);
	    defines["LAST_EXCLUSION_TILE"] = cu.intToString(endExclusionIndex);
	    // end of tiles with exclusions setup
	    //
	    if (useCutoff)
	        defines["USE_CUTOFF"] = "1";
	    double cutoff = force.getCutoff();
	    defines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);
//	    cout << "CUTOFF_SQUARED = " << cu.doubleToString(cutoff*cutoff) << std::endl;
//	    std::cout << "usePeriodic = " << usePeriodic << std::endl;
	    if (usePeriodic) //usePeriodic not being set
	        defines["USE_PERIODIC"] = "1";

	    findNeighborsWorkgroupSize = 128; //TODO: error when running findNeighbors

	    forceWorkgroupSize = 128;
	    defines["FIND_NEIGHBORS_WORKGROUP_SIZE"] = cu.intToString(findNeighborsWorkgroupSize);

	    CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaMBPolKernelSources::multibodyLibrary + CudaMBPolKernelSources::threebodyForcePolynomial + CudaMBPolKernelSources::threebodyForce, defines);
	    computeThreeBodyForceKernel = cu.getKernel(module, "computeThreeBodyForce");
	    blockBoundsKernel = cu.getKernel(module, "findBlockBounds");
		neighborsKernel = cu.getKernel(module, "findNeighbors");
		startIndicesKernel = cu.getKernel(module, "computeNeighborStartIndices");
		copyPairsKernel = cu.getKernel(module, "copyPairsToNeighborList");
	    // Add an interaction to the default nonbonded kernel.  This doesn't actually do any calculations.  It's
	    // just so that CudaNonbondedUtilities will build the exclusion flags and maintain the neighbor list.

	    cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, false, force.getCutoff(), exclusions, "", force.getForceGroup());
//	    cu.getNonbondedUtilities().setUsePadding(true);
	    cu.addForce(new CudaMBPolThreeBodyForceInfo(force));
}


#define CHECK_RESULT(result, prefix) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

double CudaCalcMBPolThreeBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//	//print positins
//    //void getPositions(std::vector<Vec3>& positions);
//	std::vector<Vec3> positions;
//	context.getPositions(positions);
//	for(std::vector<Vec3>::iterator it = positions.begin(); it != positions.end(); ++it) {
//		std::cout<< (*it)[0] << std::endl;
//	}
	///////////////////////////From Two Body/////////////////////////////////////

    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();

//    int startTileIndex = nb.getStartTileIndex();
//    int numTileIndices = nb.getNumTiles();
//    unsigned int maxTiles;
//    if (nb.getUseCutoff()) {
//        maxTiles = nb.getInteractingTiles().getSize();
//    }
//
//    void* args[] = {&cu.getForce().getDevicePointer(),
//        &cu.getEnergyBuffer().getDevicePointer(),
//        &cu.getPosq().getDevicePointer(),
//        &nb.getExclusionTiles().getDevicePointer(),
//        &startTileIndex,
//        &numTileIndices,
//        &cu.getNonbondedUtilities().getInteractingTiles().getDevicePointer(),
//        &cu.getNonbondedUtilities().getInteractionCount().getDevicePointer(),
//        cu.getPeriodicBoxSizePointer(),
//        cu.getInvPeriodicBoxSizePointer(),
//        &maxTiles,
//       // &cu.getNonbondedUtilities().getBlock().getDevicePointer(),
//        &cu.getNonbondedUtilities().getInteractingAtoms().getDevicePointer()
//    };
    ////////////////////////////////////////////////////////////////////////////////////////


    if (!hasInitializedKernel) {
        hasInitializedKernel = true;

        forceArgs.push_back(&cu.getForce().getDevicePointer());
		forceArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
		forceArgs.push_back(&cu.getPosq().getDevicePointer());
		forceArgs.push_back(cu.getPeriodicBoxSizePointer());
		forceArgs.push_back(cu.getInvPeriodicBoxSizePointer());
		forceArgs.push_back(cu.getPeriodicBoxVecXPointer());
		forceArgs.push_back(cu.getPeriodicBoxVecYPointer());
		forceArgs.push_back(cu.getPeriodicBoxVecZPointer());
		if (nb.getUseCutoff()) {
			forceArgs.push_back(&neighbors->getDevicePointer());
			forceArgs.push_back(&neighborStartIndex->getDevicePointer());
		}


		if (nb.getUseCutoff()) {
			// Set arguments for the block bounds kernel.

			blockBoundsArgs.push_back(cu.getPeriodicBoxSizePointer());
			blockBoundsArgs.push_back(cu.getInvPeriodicBoxSizePointer());
			blockBoundsArgs.push_back(cu.getPeriodicBoxVecXPointer());
			blockBoundsArgs.push_back(cu.getPeriodicBoxVecYPointer());
			blockBoundsArgs.push_back(cu.getPeriodicBoxVecZPointer());
			blockBoundsArgs.push_back(&cu.getPosq().getDevicePointer());
			blockBoundsArgs.push_back(&blockCenter->getDevicePointer());
			blockBoundsArgs.push_back(&blockBoundingBox->getDevicePointer());
			blockBoundsArgs.push_back(&numNeighborPairs->getDevicePointer());

			// Set arguments for the neighbor list kernel.

			neighborsArgs.push_back(cu.getPeriodicBoxSizePointer());
			neighborsArgs.push_back(cu.getInvPeriodicBoxSizePointer());
			neighborsArgs.push_back(cu.getPeriodicBoxVecXPointer());
			neighborsArgs.push_back(cu.getPeriodicBoxVecYPointer());
			neighborsArgs.push_back(cu.getPeriodicBoxVecZPointer());
			neighborsArgs.push_back(&cu.getPosq().getDevicePointer());
			neighborsArgs.push_back(&blockCenter->getDevicePointer());
			neighborsArgs.push_back(&blockBoundingBox->getDevicePointer());
			neighborsArgs.push_back(&neighborPairs->getDevicePointer());
			neighborsArgs.push_back(&numNeighborPairs->getDevicePointer());
			neighborsArgs.push_back(&numNeighborsForAtom->getDevicePointer());
			neighborsArgs.push_back(&maxNeighborPairs);
//			if (exclusions != NULL) {
//				neighborsArgs.push_back(&exclusions->getDevicePointer());
//				neighborsArgs.push_back(&exclusionStartIndex->getDevicePointer());
//			}

			// Set arguments for the kernel to find neighbor list start indices.

			startIndicesArgs.push_back(&numNeighborsForAtom->getDevicePointer());
			startIndicesArgs.push_back(&neighborStartIndex->getDevicePointer());
			startIndicesArgs.push_back(&numNeighborPairs->getDevicePointer());
			startIndicesArgs.push_back(&maxNeighborPairs);

			// Set arguments for the kernel to assemble the final neighbor list.

			copyPairsArgs.push_back(&neighborPairs->getDevicePointer());
			copyPairsArgs.push_back(&neighbors->getDevicePointer());
			copyPairsArgs.push_back(&numNeighborPairs->getDevicePointer());
			copyPairsArgs.push_back(&maxNeighborPairs);
			copyPairsArgs.push_back(&numNeighborsForAtom->getDevicePointer());
			copyPairsArgs.push_back(&neighborStartIndex->getDevicePointer());
	   }
    }
    CUevent event;

    while (true) {
		int* numPairs = (int*) cu.getPinnedBuffer();
		if (nb.getUseCutoff()) {
			cu.executeKernel(blockBoundsKernel, &blockBoundsArgs[0], cu.getNumAtomBlocks());
			//CHECK_RESULT(cuEventRecord(event, 0), "Error recording event for CustomManyParticleForce");

			cu.executeKernel(neighborsKernel, &neighborsArgs[0], cu.getNumAtoms(), findNeighborsWorkgroupSize);
			// We need to make sure there was enough memory for the neighbor list.  Download the
			// information asynchronously so kernels can be running at the same time.
			// download originally was set false,
		    cudaDeviceSynchronize();

			numNeighborPairs->download(numPairs, true);

			cu.executeKernel(startIndicesKernel, &startIndicesArgs[0], 256, 256, 256*sizeof(int));

			cu.executeKernel(copyPairsKernel, &copyPairsArgs[0], maxNeighborPairs);
		}
		int maxThreads = min(cu.getNumAtoms()*forceWorkgroupSize, cu.getEnergyBuffer().getSize());

		//TODO: implement the computation of the three body kernels
		cu.executeKernel(computeThreeBodyForceKernel, &forceArgs[0], maxThreads, forceWorkgroupSize);
		if (nb.getUseCutoff()) {
			 //Make sure there was enough memory for the neighbor list.

			//CHECK_RESULT(cuEventSynchronize(event), "Error synchronizing on event for CustomManyParticleForce");
			if (*numPairs > maxNeighborPairs) {
				// Resize the arrays and run the calculation again.

				delete neighborPairs;
				neighborPairs = NULL;
				delete neighbors;
				neighbors = NULL;
				maxNeighborPairs = (int) (1.1*(*numPairs));
				neighborPairs = CudaArray::create<int2>(cu, maxNeighborPairs, "customManyParticleNeighborPairs");
				neighbors = CudaArray::create<int>(cu, maxNeighborPairs, "customManyParticleNeighbors");
				forceArgs[5] = &neighbors->getDevicePointer();
				neighborsArgs[5] = &neighborPairs->getDevicePointer();
				copyPairsArgs[0] = &neighborPairs->getDevicePointer();
				copyPairsArgs[1] = &neighbors->getDevicePointer();
				continue;
			}
		}
		std::cout << "Number of Neighbor Pairs: " << *numPairs << std::endl;
//		std::cout<< "cutoff = " << nb.getCutoffDistance() << std::endl;

		break;
	}
	return 0.0;

}

void CudaCalcMBPolThreeBodyForceKernel::copyParametersToContext(ContextImpl& context, const MBPolThreeBodyForce& force) {
    cu.setAsCurrent();
    throw OpenMMException(" CudaCalcMBPolThreeBodyForceKernel::copyParametersToContext not implemented");

    cu.invalidateMolecules();
}
