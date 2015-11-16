typedef struct {
    real x, y, z;
    real fx, fy, fz;
} AtomData;

#define Oa  0
#define Ha1 1
#define Ha2 2
#define Ob  3
#define Hb1 4
#define Hb2 5
#define Oc  6
#define Hc1 7
#define Hc2 8

#define r3i 0.000000000000000e+00 // A
#define r3f 4.500000000000000e+00 // A
#define kHH_intra -2.254303257872797e+00 // A^(-1)
#define kOH_intra -2.564771901404151e-01 // A^(-1)
#define kHH 4.247920544074333e-01 // A^(-1)
#define kOH 8.128985941165371e-01 // A^(-1)
#define kOO 3.749457984616480e-02 // A^(-1)
#define dHH_intra 1.690594510379166e+00 // A
#define dOH_intra -2.999999868517452e+00 // A
#define dHH 3.499031358429095e+00 // A
#define dOH 4.854042963514281e+00 // A
#define dOO 4.816312044947604e-08 // A

extern "C" __device__ void computeVar(real k, real r0, real3 * a1, real3 * a2, real * var)
{
    real3 dx = (*a1-*a2)*NM_TO_A;
    real dsq = dx.x*dx.x + dx.y*dx.y + dx.z*dx.z;
    real d = SQRT(dsq);

    *var = EXP(-k*(d - r0));
}

extern "C" __device__ void computeGVar(real * g, real k, real r0, real3 * a1, real3 * a2, real3 * g1, real3 * g2)
{
    real3 dx = (*a1-*a2)*NM_TO_A;

    real dsq = dot(dx,dx);
    real d = SQRT(dsq);

    real gg = - k*(*g)*EXP(-k*(d - r0))/d;

    real cal2joule = 4.184;

    gg *= cal2joule * 10.*-1;
    
    g1->x += gg*dx.x;
    g2->x -= gg*dx.x;
    
    g1->y += gg*dx.y;
    g2->y -= gg*dx.y;
    
    g1->z += gg*dx.z;
    g2->z -= gg*dx.z;
}

extern "C" __device__ void evaluateSwitchFunc(real r, real * s, real * g)
{
    if (r > r3f) {
        *g = 0.0;
        *s = 0.0;
    } else if (r > r3i) {
        real t1 = M_PI/(r3f - r3i);
        real x = (r - r3i)*t1;
        *g = - SIN(x)*t1/2.0;
        *s = (1.0 + COS(x))/2.0;
    } else {
        *g = 0.0;
        *s = 1.0;
    }
}
extern "C" __device__ real computeInteraction(
		const unsigned int atom1,
        const unsigned int atom2,
        const unsigned int atom3,
        const real4* __restrict__ posq,
        const real4* periodicBoxSize,
        real3 * forces) {
        			real tempEnergy = 0.0f;
        			// 3 water
      				real3 positions[9];
   			     	// first water
    			    for (int i = 0; i < 3; i++) {
     			       	positions[Oa + i] = make_real3( posq[atom1+i].x * NM_TO_A,
                       	                   			   	posq[atom1+i].y * NM_TO_A,
                             		               	   	posq[atom1+i].z * NM_TO_A);
    			      	positions[Ob + i] = make_real3( posq[atom2+i].x * NM_TO_A,
                       			                       	posq[atom2+i].y * NM_TO_A,
                                                       	posq[atom2+i].z * NM_TO_A);
     			   	   	positions[Oc + i] = make_real3( posq[atom3+i].x * NM_TO_A,
                       			                       	posq[atom3+i].y * NM_TO_A,
                                			           	posq[atom3+i].z * NM_TO_A);
      				}
#if USE_PERIODIC
         			imageMolecules(periodicBoxSize, positions);
#endif

		real3 rab, rac, rbc;
		real drab(0), drac(0), drbc(0);

		rab = (positions[Oa] - positions[Ob])*NM_TO_A;
		drab += dot(rab, rab);

		rac = (positions[Oa] - positions[Oc])*NM_TO_A;
		drac += dot(rac, rac);

		rbc = (positions[Ob] - positions[Oc])*NM_TO_A;
		drbc += dot(rbc, rbc);

		drab = SQRT(drab);
		drac = SQRT(drac);
		drbc = SQRT(drbc);
		
		real cal2joule = 4.184;	


		if ((drab < 2) or (drac < 2) or (drbc < 2))
             tempEnergy = 0.;
        else {
        	real x[36];
        	int i = 0;
        	computeVar(kHH_intra, dHH_intra, positions +Ha1, positions +Ha2, x+i); ++i;
        	computeVar(kHH_intra, dHH_intra, positions +Hb1, positions +Hb2, x+i); ++i;
        	computeVar(kHH_intra, dHH_intra, positions +Hc1, positions +Hc2, x+i); ++i;
        	computeVar(kOH_intra, dOH_intra, positions +Oa,	 positions +Ha1, x+i); ++i;
        	computeVar(kOH_intra, dOH_intra, positions +Oa,	 positions +Ha2, x+i); ++i;
        	computeVar(kOH_intra, dOH_intra, positions +Ob,	 positions +Hb1, x+i); ++i; //5
        	computeVar(kOH_intra, dOH_intra, positions +Ob,	 positions +Hb2, x+i); ++i;
        	computeVar(kOH_intra, dOH_intra, positions +Oc,	 positions +Hc1, x+i); ++i;
        	computeVar(kOH_intra, dOH_intra, positions +Oc,	 positions +Hc2, x+i); ++i;

        	computeVar(kHH, dHH, positions +Ha1, positions +Hb1, x+i); ++i;
        	computeVar(kHH, dHH, positions +Ha1, positions +Hb2, x+i); ++i; //10
        	computeVar(kHH, dHH, positions +Ha1, positions +Hc1, x+i); ++i;
        	computeVar(kHH, dHH, positions +Ha1, positions +Hc2, x+i); ++i;
        	computeVar(kHH, dHH, positions +Ha2, positions +Hb1, x+i); ++i;
        	computeVar(kHH, dHH, positions +Ha2, positions +Hb2, x+i); ++i;
        	computeVar(kHH, dHH, positions +Ha2, positions +Hc1, x+i); ++i; //15
        	computeVar(kHH, dHH, positions +Ha2, positions +Hc2, x+i); ++i;
        	computeVar(kHH, dHH, positions +Hb1, positions +Hc1, x+i); ++i;
        	computeVar(kHH, dHH, positions +Hb1, positions +Hc2, x+i); ++i;
        	computeVar(kHH, dHH, positions +Hb2, positions +Hc1, x+i); ++i;
        	computeVar(kHH, dHH, positions +Hb2, positions +Hc2, x+i); ++i; //20
        	computeVar(kOH, dOH, positions +Oa, positions +Hb1, x+i); ++i;
        	computeVar(kOH, dOH, positions +Oa, positions +Hb2, x+i); ++i;
        	computeVar(kOH, dOH, positions +Oa, positions +Hc1, x+i); ++i;
        	computeVar(kOH, dOH, positions +Oa, positions +Hc2, x+i); ++i;
        	computeVar(kOH, dOH, positions +Ob, positions +Ha1, x+i); ++i; //25
        	computeVar(kOH, dOH, positions +Ob, positions +Ha2, x+i); ++i;
        	computeVar(kOH, dOH, positions +Ob, positions +Hc1, x+i); ++i;
        	computeVar(kOH, dOH, positions +Ob, positions +Hc2, x+i); ++i;
        	computeVar(kOH, dOH, positions +Oc, positions +Ha1, x+i); ++i;
        	computeVar(kOH, dOH, positions +Oc, positions +Ha2, x+i); ++i; //30
        	computeVar(kOH, dOH, positions +Oc, positions +Hb1, x+i); ++i;
        	computeVar(kOH, dOH, positions +Oc, positions +Hb2, x+i); ++i;
        	computeVar(kOO, dOO, positions +Oa, positions +Ob, x+i); ++i;
        	computeVar(kOO, dOO, positions +Oa, positions +Oc, x+i); ++i;
        	computeVar(kOO, dOO, positions +Ob, positions +Oc, x+i); ++i; //35

        	real g[36];
            tempEnergy = poly_3b_v2x_eval(x, g);

			real gab, gac, gbc;
			real sab, sac, sbc;
			evaluateSwitchFunc(drab, &gab, &sab);
			evaluateSwitchFunc(drac, &gac, &sac);
			evaluateSwitchFunc(drbc, &gbc, &sbc);

			real s = sab*sac + sab*sbc + sac*sbc;

			for (int n = 0; n < 36; ++n)
				g[n] *= s;

			//extern "C" __device__ void computeGVar(real g, real k, real r0, real3 * a1, real3 * a2, real3 * g1, real3 * g2)
        	i = 0;
        	computeGVar(g+i, kHH_intra, dHH_intra, positions +Ha1, positions +Ha2, forces +Ha1, forces +Ha2); ++i; //0
        	computeGVar(g+i, kHH_intra, dHH_intra, positions +Hb1, positions +Hb2, forces +Hb1, forces +Hb2); ++i;
        	computeGVar(g+i, kHH_intra, dHH_intra, positions +Hc1, positions +Hc2, forces +Hc1, forces +Hc2); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions +Oa, positions +Ha1, forces +Oa, forces +Ha1); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions +Oa, positions +Ha2, forces +Oa, forces +Ha1); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions +Ob, positions +Hb1, forces +Ob, forces +Hb1); ++i; //5
        	computeGVar(g+i, kOH_intra, dOH_intra, positions +Ob, positions +Hb2, forces +Ob, forces +Hb2); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions +Oc, positions +Hc1, forces +Oc, forces +Hc1); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions +Oc, positions +Hc2, forces +Oc, forces +Hc2); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Ha1, positions +Hb1, forces +Ha1, forces +Hb1); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Ha1, positions +Hb2, forces +Ha1, forces +Hb2); ++i; //10
        	computeGVar(g+i, kHH, dHH, positions +Ha1, positions +Hc1, forces +Ha1, forces +Hc1); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Ha1, positions +Hc2, forces +Ha1, forces +Hc2); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Ha2, positions +Hb1, forces +Ha2, forces +Hb1); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Ha2, positions +Hb2, forces +Ha2, forces +Hb2); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Ha2, positions +Hc1, forces +Ha2, forces +Hc1); ++i; //15
        	computeGVar(g+i, kHH, dHH, positions +Ha2, positions +Hc2, forces +Ha2, forces +Hc2); ++i;
			computeGVar(g+i, kHH, dHH, positions +Hb1, positions +Hc1, forces +Hb1, forces +Hc1); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Hb1, positions +Hc2, forces +Hb1, forces +Hc2); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Hb2, positions +Hc1, forces +Hb2, forces +Hc1); ++i;
        	computeGVar(g+i, kHH, dHH, positions +Hb2, positions +Hc2, forces +Hb2, forces +Hc2); ++i; //20
        	computeGVar(g+i, kOH, dOH, positions +Oa, positions +Hb1, forces +Oa, forces +Hb1); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Oa, positions +Hb2, forces +Oa, forces +Hb2); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Oa, positions +Hc1, forces +Oa, forces +Hc1); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Oa, positions +Hc2, forces +Oa, forces +Hc2); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Ob, positions +Ha1, forces +Ob, forces +Ha1); ++i; //25
        	computeGVar(g+i, kOH, dOH, positions +Ob, positions +Ha2, forces +Ob, forces +Ha2); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Ob, positions +Hc1, forces +Ob, forces +Hc1); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Ob, positions +Hc2, forces +Ob, forces +Hc2); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Oc, positions +Ha1, forces +Oc, forces +Ha1); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Oc, positions +Ha2, forces +Oc, forces +Ha2); ++i; //30
        	computeGVar(g+i, kOH, dOH, positions +Oc, positions +Hb1, forces +Oc, forces +Hb1); ++i;
        	computeGVar(g+i, kOH, dOH, positions +Oc, positions +Hb2, forces +Oc, forces +Hb2); ++i;
        	computeGVar(g+i, kOO, dOO, positions +Oa, positions +Ob, forces +Oa, forces +Ob); ++i;
        	computeGVar(g+i, kOO, dOO, positions +Oa, positions +Oc, forces +Oa, forces +Oc); ++i;
        	computeGVar(g+i, kOO, dOO, positions +Ob, positions +Oc, forces +Ob, forces +Oc); ++i; //35

        	gab *= (sac + sbc)*tempEnergy/drab;
            gac *= (sab + sbc)*tempEnergy/drac;
            gbc *= (sab + sac)*tempEnergy/drbc;

			tempEnergy *= s;
			
			forces[Oa].x += (gab*rab.x + gac*rac.x) * cal2joule * -NM_TO_A;
			forces[Ob].x += (gbc*rbc.x - gab*rab.x) * cal2joule * -NM_TO_A;
			forces[Oc].x -= (gac*rac.x + gbc*rbc.x) * cal2joule * -NM_TO_A;

            forces[Oa].y += (gab*rab.y + gac*rac.y) * cal2joule * -NM_TO_A;            
            forces[Ob].y += (gbc*rbc.y - gab*rab.y) * cal2joule * -NM_TO_A;
            forces[Oc].y -= (gac*rac.y + gbc*rbc.y) * cal2joule * -NM_TO_A;

            forces[Oa].z += (gab*rab.z + gac*rac.z) * cal2joule * -NM_TO_A;
            forces[Ob].z += (gbc*rbc.z - gab*rab.z) * cal2joule * -NM_TO_A;
            forces[Oc].z -= (gac*rac.z + gbc*rbc.z) * cal2joule * -NM_TO_A;

// Is it okay to calculate the force in the shared variable like in cuda 2body
// or should we calculate in a seperate variable and add to the actual at one
// time at the end like in refrence code ??
        }
        real energy = tempEnergy * cal2joule;

        return energy;

}


extern "C" __global__ void computeThreeBodyForce(

        // const unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ tempForceBuffers) {
        unsigned long long* __restrict__ forceBuffers,
        real* __restrict__ energyBuffer,
        const real4* __restrict__ posq,
        const ushort2* __restrict__ exclusionTiles,
        unsigned int startTileIndex,
        unsigned int numTileIndices
#ifdef USE_CUTOFF
        , const int* __restrict__ tiles,
        const unsigned int* __restrict__ interactionCount,
        real4 periodicBoxSize,
        real4 invPeriodicBoxSize,
        unsigned int maxTiles,
        // const real4* __restrict__ blockSize,
        const unsigned int* __restrict__ interactingAtoms
#endif
        ) {
    real energy = 0.0f;

    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;

}




///////////// things added for neighbor list //////////////////////////

/**
 * Convert a real4 to a real3 by removing its last element.
 */
inline __device__ real3 trim(real4 v) {
    return make_real3(v.x, v.y, v.z);
}
/**
 * Compute the difference between two vectors, taking periodic boundary conditions into account
 * and setting the fourth component to the squared magnitude.
 */
inline __device__ real4 delta(real3 vec1, real3 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = make_real4(vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0.0f);
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_DELTA(result)
#endif
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}
/**
 * Find a list of neighbors for each atom.
 */
extern "C" __global__ void findNeighbors(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, const real4* __restrict__ blockCenter, const real4* __restrict__ blockBoundingBox, int2* __restrict__ neighborPairs,
        int* __restrict__ numNeighborPairs, int* __restrict__ numNeighborsForAtom, int maxNeighborPairs
#ifdef USE_EXCLUSIONS
        , int* __restrict__ exclusions, int* __restrict__ exclusionStartIndex
#endif
        ) {
    __shared__ real3 positionCache[/*FIND_NEIGHBORS_WORKGROUP_SIZE = 128*/ 128];
    int indexInWarp = threadIdx.x%32;
    for (int atom1 = blockIdx.x*blockDim.x+threadIdx.x; atom1 < PADDED_NUM_ATOMS; atom1 += blockDim.x*gridDim.x) {
        // Load data for this atom.  Note that all threads in a warp are processing atoms from the same block.

        real3 pos1 = trim(posq[atom1]);
        int block1 = atom1/TILE_SIZE;
        real4 blockCenter1 = blockCenter[block1];
        real4 blockSize1 = blockBoundingBox[block1];
        int totalNeighborsForAtom1 = 0;

        // Loop over atom blocks to search for neighbors.  The threads in a warp compare block1 against 32
        // other blocks in parallel.

#ifdef USE_CENTRAL_PARTICLE
        int startBlock = 0;
#else
        int startBlock = block1;
#endif
        for (int block2Base = startBlock; block2Base < NUM_BLOCKS; block2Base += 32) {
            int block2 = block2Base+indexInWarp;
            bool includeBlock2 = (block2 < NUM_BLOCKS);
            if (includeBlock2) {
                real4 blockCenter2 = blockCenter[block2];
                real4 blockSize2 = blockBoundingBox[block2];
                real4 blockDelta = blockCenter1-blockCenter2;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                blockDelta.x = max(0.0f, fabs(blockDelta.x)-blockSize1.x-blockSize2.x);
                blockDelta.y = max(0.0f, fabs(blockDelta.y)-blockSize1.y-blockSize2.y);
                blockDelta.z = max(0.0f, fabs(blockDelta.z)-blockSize1.z-blockSize2.z);
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < CUTOFF_SQUARED);
            }

            // Loop over any blocks we identified as potentially containing neighbors.

            int includeBlockFlags = __ballot(includeBlock2);
            while (includeBlockFlags != 0) {
                int i = __ffs(includeBlockFlags)-1;
                includeBlockFlags &= includeBlockFlags-1;
                int block2 = block2Base+i;

                // Loop over atoms in this block.

                int start = block2*TILE_SIZE;
                int included[TILE_SIZE];
                int numIncluded = 0;
                positionCache[threadIdx.x] = trim(posq[start+indexInWarp]);
                if (atom1 < NUM_ATOMS) {
                    for (int j = 0; j < 32; j++) {
                        int atom2 = start+j;
                        real3 pos2 = positionCache[threadIdx.x-indexInWarp+j];

                        // Decide whether to include this atom pair in the neighbor list.

                        real4 atomDelta = delta(pos1, pos2, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
#ifdef USE_CENTRAL_PARTICLE
                        bool includeAtom = (atom2 != atom1 && atom2 < NUM_ATOMS && atomDelta.w < CUTOFF_SQUARED);
#else
                        bool includeAtom = (atom2 > atom1 && atom2 < NUM_ATOMS && atomDelta.w < CUTOFF_SQUARED);
#endif
#ifdef USE_EXCLUSIONS
                        if (includeAtom)
                            includeAtom &= !isInteractionExcluded(atom1, atom2, exclusions, exclusionStartIndex);
#endif
                        if (includeAtom)
                            included[numIncluded++] = atom2;
                    }
                }

                // If we found any neighbors, store them to the neighbor list.

                if (numIncluded > 0) {
                    int baseIndex = atomicAdd(numNeighborPairs, numIncluded);
                    if (baseIndex+numIncluded <= maxNeighborPairs)
                        for (int j = 0; j < numIncluded; j++)
                            neighborPairs[baseIndex+j] = make_int2(atom1, included[j]);
                    totalNeighborsForAtom1 += numIncluded;
                }
            }
        }
        numNeighborsForAtom[atom1] = totalNeighborsForAtom1;
    }
}

/**
 * Sum the neighbor counts to compute the start position of each atom.  This kernel
 * is executed as a single work group.
 */
extern "C" __global__ void computeNeighborStartIndices(int* __restrict__ numNeighborsForAtom, int* __restrict__ neighborStartIndex,
            int* __restrict__ numNeighborPairs, int maxNeighborPairs) {
    extern __shared__ unsigned int posBuffer[];
    if (*numNeighborPairs > maxNeighborPairs) {
        // There wasn't enough memory for the neighbor list, so we'll need to rebuild it.  Set the neighbor start
        // indices to indicate no neighbors for any atom.

        for (int i = threadIdx.x; i <= NUM_ATOMS; i += blockDim.x)
            neighborStartIndex[i] = 0;
        return;
    }
    unsigned int globalOffset = 0;
    for (unsigned int startAtom = 0; startAtom < NUM_ATOMS; startAtom += blockDim.x) {
        // Load the neighbor counts into local memory.

        unsigned int globalIndex = startAtom+threadIdx.x;
        posBuffer[threadIdx.x] = (globalIndex < NUM_ATOMS ? numNeighborsForAtom[globalIndex] : 0);
        __syncthreads();

        // Perform a parallel prefix sum.

        for (unsigned int step = 1; step < blockDim.x; step *= 2) {
            unsigned int add = (threadIdx.x >= step ? posBuffer[threadIdx.x-step] : 0);
            __syncthreads();
            posBuffer[threadIdx.x] += add;
            __syncthreads();
        }

        // Write the results back to global memory.

        if (globalIndex < NUM_ATOMS) {
            neighborStartIndex[globalIndex+1] = posBuffer[threadIdx.x]+globalOffset;
            numNeighborsForAtom[globalIndex] = 0; // Clear this so the next kernel can use it as a counter
        }
        globalOffset += posBuffer[blockDim.x-1];
    }
    if (threadIdx.x == 0)
        neighborStartIndex[0] = 0;
}

/**
 * Assemble the final neighbor list.
 */
extern "C" __global__ void copyPairsToNeighborList(const int2* __restrict__ neighborPairs, int* __restrict__ neighbors, int* __restrict__ numNeighborPairs,
            int maxNeighborPairs, int* __restrict__ numNeighborsForAtom, const int* __restrict__ neighborStartIndex) {
    int actualPairs = *numNeighborPairs;
    if (actualPairs > maxNeighborPairs)
        return; // There wasn't enough memory for the neighbor list, so we'll need to rebuild it.
    for (unsigned int index = blockDim.x*blockIdx.x+threadIdx.x; index < actualPairs; index += blockDim.x*gridDim.x) {
        int2 pair = neighborPairs[index];
        int startIndex = neighborStartIndex[pair.x];
        int offset = atomicAdd(numNeighborsForAtom+pair.x, 1);
        neighbors[startIndex+offset] = pair.y;
    }
}

/**
 * Find a bounding box for the atoms in each block.
 */
extern "C" __global__ void findBlockBounds(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, real4* __restrict__ blockCenter, real4* __restrict__ blockBoundingBox, int* __restrict__ numNeighborPairs) {
	int index = blockIdx.x*blockDim.x+threadIdx.x;
    int base = index*TILE_SIZE;
    while (base < NUM_ATOMS) {
        real4 pos = posq[base];
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_POS(pos)
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, NUM_ATOMS);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
            minPos = make_real4(min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
            maxPos = make_real4(max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
        }
        real4 blockSize = 0.5f*(maxPos-minPos);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = 0.5f*(maxPos+minPos);
        index += blockDim.x*gridDim.x;
        base = index*TILE_SIZE;
    }
    if (blockIdx.x == 0 && threadIdx.x == 0)
        *numNeighborPairs = 0;
}

    
       
       
       
       
       
       
       
       
       
       
