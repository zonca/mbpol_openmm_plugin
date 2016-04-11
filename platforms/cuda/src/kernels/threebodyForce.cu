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

extern "C" __device__ void imageOxygens(const real4 * box, real3 * allPositions)
{
    // Take first oxygen as central atom

    // Now image the oxygen of the second &molecule
    imageParticles(box, &allPositions[0], &allPositions[1]);
}

/**
 * Find a list of neighbors for each atom.
 */
extern "C" __global__ void findNeighbors(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, const real4* __restrict__ blockCenter, const real4* __restrict__ blockBoundingBox, int2* __restrict__ neighborPairs,
        int* __restrict__ numNeighborPairs, int* __restrict__ numNeighborsForAtom, int maxNeighborPairs
//#ifdef USE_EXCLUSIONS
//        , int* __restrict__ exclusions, int* __restrict__ exclusionStartIndex
//#endif
        ) {
//	printf("threadIdx.x = %d\n", threadIdx.x);
//	//if (threadIdx.x == 3) {
//		printf("FIND_NEIGHBORS_WORKGROUP_SIZE = %d\n", FIND_NEIGHBORS_WORKGROUP_SIZE);
//		printf("PADDED_NUM_ATOMS = %d\n", PADDED_NUM_ATOMS); //320 for 100 mol
//		printf("NUM_BLOCKS = %d\n", NUM_BLOCKS);
//		printf("TILE_SIZE = %d\n", TILE_SIZE);
//		printf("CUTOFF_SQUARED = %d\n", CUTOFF_SQUARED);
//		printf("NUM_ATOMS = %d\n",NUM_ATOMS);
////		for (int i = 0; i<NUM_BLOCKS; i++) {
////			printf("blockCenter[%d] = {%lf, %lf, %lf, %lf}\n", i, blockCenter[i].x, blockCenter[i].y, blockCenter[i].z, blockCenter[i].w);
////			printf("blockBoundingBox[%d] = {%lf, %lf, %lf, %lf}\n", i, blockBoundingBox[i].x, blockBoundingBox[i].y, blockBoundingBox[i].z, blockBoundingBox[i].w);
////		//}
//		printf("blockDim.x*gridDim.x = %d\n", blockDim.x*gridDim.x);
//		printf("blockIdx.x*blockDim.x+threadIdx.x = %d\n",blockIdx.x*blockDim.x+threadIdx.x);
//	//}
    __shared__ real3 positionCache[FIND_NEIGHBORS_WORKGROUP_SIZE];
    int indexInWarp = threadIdx.x%32;
    for (int atom1 = blockIdx.x*blockDim.x+threadIdx.x; atom1 < PADDED_NUM_ATOMS; atom1 += blockDim.x*gridDim.x) {
        // Load data for this atom.  Note that all threads in a warp are processing atoms from the same block.
    	
    	//FIXME: temporary fix to num neigbors for atom is not paddded so the size is 9
    	// but needs to be executed for paddedNumAtoms which is 32
    	if (atom1 >= NUM_ATOMS)
        	continue;
    	//FIXME: temp fix to building neighborlist to include only oxygen particles
    	if (atom1%3 != 0)
    		continue;
        
        real3 pos1 = trim(posq[atom1]);
        //printf("find neighbors pos1 = <%10lf, %10lf, %10lf>\n",  pos1.x, pos1.y, pos1.z);
        
        int block1 = atom1/TILE_SIZE;
        real4 blockCenter1 = blockCenter[block1];
        real4 blockSize1 = blockBoundingBox[block1];
        int totalNeighborsForAtom1 = 0;
        
        // Loop over atom blocks to search for neighbors.  The threads in a warp compare block1 against 32
        // other blocks in parallel.

        int startBlock = block1;

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
                    	//FIXME: temp fix to building neighborlist to include only oxygen particles
                    	if (atom2%3 != 0)
                    		continue;
                        real3 pos2 = positionCache[threadIdx.x-indexInWarp+j];
//                        printf("threadIdx.x-indexInWarp+j = %d\n", threadIdx.x-indexInWarp+j);
                        // Decide whether to include this atom pair in the neighbor list.
                                      
                        real4 atomDelta = delta(pos1, pos2, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);

                        bool includeAtom = (atom2 > atom1 && atom2 < NUM_ATOMS && atomDelta.w < CUTOFF_SQUARED);

                        if (includeAtom) {
                            included[numIncluded++] = atom2;
//                            printf("Tid %d, pair found: %d, %d\n", threadIdx.x, atom1, atom2);
                        }
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
    //printf("completed call of findNeighbors: %d\n", threadIdx.x);
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
    //printf("actualPairs = %d\n", actualPairs);
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
	//printf("Thread %d\n", index);
    int base = index*TILE_SIZE;
    while (base < NUM_ATOMS) {
        real4 pos = posq[base];
//        printf("Thread %d blockbounds pos = <%10lf, %10lf, %10lf>\n", index,  pos.x, pos.y, pos.z);
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

////////////////////////////////////////////////////////////////


extern "C" __device__ void computeVar(real k, real r0, real3 * a1, real3 * a2, real * var)
{
    real3 dx = (*a1-*a2)*NM_TO_A;
    real dsq = dx.x*dx.x + dx.y*dx.y + dx.z*dx.z;
    real d = SQRT(dsq);

    *var = EXP(-k*(d - r0));
}

extern "C" __device__ void computeGVar(real * g, 
									   real k, 
									   real r0, 
									   real3 * a1, real3 * a2, 
									   real3 * g1, real3 * g2)
{
//	printf ("parameter check of computeGVar: g = %10lf, k = %10lf, r0 = %10lf, a1 = < %10lf, %10lf, %10lf>, a2 = < %10lf, %10lf, %10lf>, g1 = < %10lf, %10lf, %10lf> g2 = < %10lf, %10lf, %10lf>\n",
//			*g, k, r0, a1->x, a1->y, a1->z, a2->x, a2->y, a2->z, g1->x, g1->y, g1->z, g2->x, g2->y, g2->z);
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

extern "C" __device__ void evaluateSwitchFunc(real r, real * g, real * s)
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

extern "C" __device__ void imageMolecules(const real4 * box, real3 * allPositions)
{

    // Take first oxygen as central atom

    // image its two hydrogens with respect of the first oxygen

    imageParticles(box, &allPositions[Oa], &allPositions[Ha1]);
    imageParticles(box, &allPositions[Oa], &allPositions[Ha2]);

    // Now image the oxygen of the second molecule

    imageParticles(box, &allPositions[Oa], &allPositions[Ob]);

    // Image the hydrogen of the second molecule with respect to the oxygen of the second molecule
    imageParticles(box, &allPositions[Ob], &allPositions[Hb1]);
    imageParticles(box, &allPositions[Ob], &allPositions[Hb2]);
    
    // Now image the oxygen of the third molecule

    imageParticles(box, &allPositions[Oa], &allPositions[Oc]);

    // Image the hydrogen of the third mo&lecule with respect to the oxygen of the second third
    imageParticles(box, &allPositions[Oc], &allPositions[Hc1]);
    imageParticles(box, &allPositions[Oc], &allPositions[Hc2]);

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
		//ordering to match ref
		positions[Oa + i] = make_real3( posq[atom3+i].x,
										posq[atom3+i].y,
										posq[atom3+i].z);
		positions[Ob + i] = make_real3( posq[atom2+i].x,
										posq[atom2+i].y,
										posq[atom2+i].z);
		positions[Oc + i] = make_real3( posq[atom1+i].x,
										posq[atom1+i].y,
										posq[atom1+i].z);
	}
	
//	int temp_int = atom1;
//	atom1 = atom2;
//	atom2 = temp_int;
	
//	real3 temp = positions[0];
//	positions[0] = positions[3];
//	positions[3] = temp;
//	
//	temp = positions[1];
//	positions[1] = positions[4];
//	positions[4] = temp;
//	
//	temp = positions[2];
//	positions[2] = positions[5];
//	positions[5] = temp;
	
	
	
	
	
//	for (int i = 0; i<9; i++) {
//		printf("before imaging positions[%d] = <%10lf, %10lf, %10lf>\n", i, 
//				(positions[i].x > 49 ? positions[i].x-50 : positions[i].x),
//				(positions[i].y > 49 ? positions[i].y-50 : positions[i].y),
//				(positions[i].z > 49 ? positions[i].z-50 : positions[i].z));
//	}
//	for (int i = 0; i<9; i++) {
//		printf("before imaging positions[%d] = <%10lf, %10lf, %10lf>\n", i, positions[i].x, positions[i].y, positions[i].z);
//	}
#ifdef USE_PERIODIC
         			imageMolecules(periodicBoxSize, positions);
#endif
         			
//	for (int i = 0; i<9; i++) {
//		printf(" after imaging positions[%d] = <%10lf, %10lf, %10lf>\n", i, positions[i].x, positions[i].y, positions[i].z);
//	}
//	for (int i = 0; i<9; i++) {
//		printf(" after imaging positions[%d] = <%10lf, %10lf, %10lf>\n", i, 
//				(positions[i].x > 49 ? positions[i].x-50 : positions[i].x),
//				(positions[i].y > 49 ? positions[i].y-50 : positions[i].y),
//				(positions[i].z > 49 ? positions[i].z-50 : positions[i].z));
//	}

		real3 rab, rac, rbc;
		real drab(0), drac(0), drbc(0);

		rab = (positions[Oa] - positions[Ob])* NM_TO_A;
		drab += dot(rab, rab);

		rac = (positions[Oa] - positions[Oc])* NM_TO_A;
		drac += dot(rac, rac);

		rbc = (positions[Ob] - positions[Oc])* NM_TO_A;
		drbc += dot(rbc, rbc);
//		printf("rab = <%lf, %lf, %lf>\n",rab.x, rab.y, rab.z);
//		printf("rac = <%lf, %lf, %lf>\n",rac.x, rac.y, rac.z);
//		printf("rbc = <%lf, %lf, %lf>\n",rbc.x, rbc.y, rbc.z);
		/* rab ra rbc good
		 ref
		 rab = <1.204804, 2.388359, 1.161075>
		 rac = <0.957227, 2.209016, -1.593952>
		 rbc = <-0.247577, -0.179343, -2.755027>
		 cuda
		rab = <1.204804, 2.388359, 1.161075>
		rac = <0.957227, 2.209016, -1.593952>
		rbc = <-0.247577, -0.179343, -2.755027>
		 */

		drab = SQRT(drab);
		drac = SQRT(drac);
		drbc = SQRT(drbc);
		real cal2joule = 4.184;	
		
		//printf("drab = %lf, drac = %lf, drbc = %lf\n", drab, drac, drbc);
		// drab drac drbc good
		//ref drab = 2.91615, drac = 2.88734, drbc = 2.77194
		//cuda drab = 2.771936, drac = 2.887337, drbc = 2.916146


		if ((drab < 2) or (drac < 2) or (drbc < 2)) {
             tempEnergy = 0.;
		}
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
            
//            for (int i = 0; i<36; i++)
//            	printf("before computeGVar g[%d] = %lf\n", i, g[i]);
            
			real gab, gac, gbc;
			real sab, sac, sbc;

			evaluateSwitchFunc(drab, &gab, &sab);
			evaluateSwitchFunc(drac, &gac, &sac);
			evaluateSwitchFunc(drbc, &gbc, &sbc);
			
			real s = sab*sac + sab*sbc + sac*sbc;
			//printf("s = %lf\n", s);
			// s is good
			for (int n = 0; n < 36; ++n)
				g[n] *= s;
			
			// zero all forces
			for (int n = 0; n < 9; ++n){
				forces[n].x = 0;
				forces[n].y = 0;
				forces[n].z = 0;
			}
				
			//extern "C" __device__ void computeGVar(real g, real k, real r0, real3 * a1, real3 * a2, real3 * g1, real3 * g2)
			for (int n = 0; n < 9; ++n)
				printf("b4 forces[%d] = <%lf, %lf, %lf>\n",n, forces[n].x ,forces[n].y, forces[n].z);

			i = 0;
        	computeGVar(g+i, kHH_intra, dHH_intra, positions+ Ha1, positions+ Ha2, forces+ Ha1, forces+ Ha2); ++i; //0
        	computeGVar(g+i, kHH_intra, dHH_intra, positions+ Hb1, positions+ Hb2, forces+ Hb1, forces+ Hb2); ++i;
        	computeGVar(g+i, kHH_intra, dHH_intra, positions+ Hc1, positions+ Hc2, forces+ Hc1, forces+ Hc2); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions+  Oa, positions+ Ha1, forces+  Oa, forces+ Ha1); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions+  Oa, positions+ Ha2, forces+  Oa, forces+ Ha2); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions+  Ob, positions+ Hb1, forces+  Ob, forces+ Hb1); ++i; //5
        	computeGVar(g+i, kOH_intra, dOH_intra, positions+  Ob, positions+ Hb2, forces+  Ob, forces+ Hb2); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions+  Oc, positions+ Hc1, forces+  Oc, forces+ Hc1); ++i;
        	computeGVar(g+i, kOH_intra, dOH_intra, positions+  Oc, positions+ Hc2, forces+  Oc, forces+ Hc2); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Ha1, positions+ Hb1, forces+ Ha1, forces+ Hb1); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Ha1, positions+ Hb2, forces+ Ha1, forces+ Hb2); ++i; //10
        	computeGVar(g+i,       kHH,       dHH, positions+ Ha1, positions+ Hc1, forces+ Ha1, forces+ Hc1); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Ha1, positions+ Hc2, forces+ Ha1, forces+ Hc2); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Ha2, positions+ Hb1, forces+ Ha2, forces+ Hb1); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Ha2, positions+ Hb2, forces+ Ha2, forces+ Hb2); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Ha2, positions+ Hc1, forces+ Ha2, forces+ Hc1); ++i; //15
        	computeGVar(g+i,       kHH,       dHH, positions+ Ha2, positions+ Hc2, forces+ Ha2, forces+ Hc2); ++i;
			computeGVar(g+i,       kHH,       dHH, positions+ Hb1, positions+ Hc1, forces+ Hb1, forces+ Hc1); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Hb1, positions+ Hc2, forces+ Hb1, forces+ Hc2); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Hb2, positions+ Hc1, forces+ Hb2, forces+ Hc1); ++i;
        	computeGVar(g+i,       kHH,       dHH, positions+ Hb2, positions+ Hc2, forces+ Hb2, forces+ Hc2); ++i; //20
        	computeGVar(g+i,       kOH,       dOH, positions+  Oa, positions+ Hb1, forces+  Oa, forces+ Hb1); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Oa, positions+ Hb2, forces+  Oa, forces+ Hb2); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Oa, positions+ Hc1, forces+  Oa, forces+ Hc1); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Oa, positions+ Hc2, forces+  Oa, forces+ Hc2); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Ob, positions+ Ha1, forces+  Ob, forces+ Ha1); ++i; //25
        	computeGVar(g+i,       kOH,       dOH, positions+  Ob, positions+ Ha2, forces+  Ob, forces+ Ha2); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Ob, positions+ Hc1, forces+  Ob, forces+ Hc1); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Ob, positions+ Hc2, forces+  Ob, forces+ Hc2); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Oc, positions+ Ha1, forces+  Oc, forces+ Ha1); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Oc, positions+ Ha2, forces+  Oc, forces+ Ha2); ++i; //30
        	computeGVar(g+i,       kOH,       dOH, positions+  Oc, positions+ Hb1, forces+  Oc, forces+ Hb1); ++i;
        	computeGVar(g+i,       kOH,       dOH, positions+  Oc, positions+ Hb2, forces+  Oc, forces+ Hb2); ++i;
        	computeGVar(g+i,       kOO,       dOO, positions+  Oa, positions+  Ob, forces+  Oa, forces+  Ob); ++i;
        	computeGVar(g+i,       kOO,       dOO, positions+  Oa, positions+  Oc, forces+  Oa, forces+  Oc); ++i;
        	computeGVar(g+i,       kOO,       dOO, positions+  Ob, positions+  Oc, forces+  Ob, forces+  Oc); ++i; //35
			//degbuging gradients
//            for (int i = 0; i<36; i++)
//          	  printf("after g_var g[%d] = %lf\n", i, g[i]);
			for (int n = 0; n < 9; ++n)
				printf("forces[%d] = <%lf, %lf, %lf>\n",n, forces[n].x ,forces[n].y, forces[n].z);
        	gab *= (sac + sbc)*tempEnergy/drab;
            gac *= (sab + sbc)*tempEnergy/drac;
            gbc *= (sab + sac)*tempEnergy/drbc;
            
//            printf("gab = %lf\n", gab);
//            printf("gac = %lf\n", gac);
//            printf("gbc = %lf\n", gbc);
            
            /* ref
             gab gac gbc good
            gab = -0.0390668
			gac = -0.0392534
			gbc = -0.0397021
            
             cuda
            gab = -0.039702
			gac = -0.039253
			gbc = -0.039067
            
             */

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

//            printf("forces[Oa] = <%lf, %lf, %lf>\n",forces[Oa].x ,forces[Oa].y, forces[Oa].z);
//            printf("forces[Ob] = <%lf, %lf, %lf>\n",forces[Ob].x ,forces[Ob].y, forces[Ob].z);
//            printf("forces[Oc] = <%lf, %lf, %lf>\n",forces[Oc].x ,forces[Oc].y, forces[Oc].z);
            
            /* ref
			allForces[Oa] = <-23.7328, -19.2615, 16.736>
			allForces[Ob] = <0.457798, 15.1341, 2.14656>
			allForces[Oc] = <-12.5181, 14.6274, 6.79418>

             
             cuda
			forces[Oa] = <189.473160, -111.636002, 77.668442>
			forces[Ob] = <13.362639, 5.023933, -245.217773>
			forces[Oc] = <-96.027763, -86.165901, -109.821854>
             */

        }
        real energy = tempEnergy * cal2joule;
		
		//printf("TempEnergy (before return) = %lf\n", energy);
		
        return energy;

}

extern "C" __global__ void computeThreeBodyForce(
        unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
#ifdef USE_CUTOFF
        , const int* __restrict__ neighbors, const int* __restrict__ neighborStartIndex
#endif
        /*PARAMETER_ARGUMENTS*/) {
    real energy = 0;
    real3 forces[9];


    // Loop over particles to be the first one in the set.

    for (int p1 = blockIdx.x; p1 < NUM_ATOMS; p1 += gridDim.x) {
        const int a1 = 0;
#ifdef USE_CUTOFF
        int firstNeighbor = neighborStartIndex[p1];
        int numNeighbors = neighborStartIndex[p1+1]-firstNeighbor;
#else
        int numNeighbors = NUM_ATOMS-p1-1;
#endif
        int numCombinations = numNeighbors*numNeighbors;
        for (int index = threadIdx.x; index < numCombinations; index += blockDim.x) {
#ifdef USE_CUTOFF
        	//FIND_ATOMS_FOR_COMBINATION_INDEX;
			int tempIndex = index;
			int a2 = 1+tempIndex%numNeighbors;
			tempIndex /= numNeighbors;
			int a3 = 1+tempIndex%numNeighbors;
			a2 = (a3%2 == 0 ? a2 : numNeighbors-a2+1);
			int p2 = neighbors[firstNeighbor-1+a2];
			int p3 = neighbors[firstNeighbor-1+a3];
#else
			//FIND_ATOMS_FOR_COMBINATION_INDEX;
        	int tempIndex = index;
        	int a2 = 1+tempIndex%numNeighbors;
        	tempIndex /= numNeighbors;
        	int a3 = 1+tempIndex%numNeighbors;
        	a2 = (a3%2 == 0 ? a2 : numNeighbors-a2+1);
        	int p2 = p1+a2;
        	int p3 = p1+a3;
#endif
			//bool includeInteraction = IS_VALID_COMBINATION;
            bool includeInteraction = (a3>a2);
#ifdef USE_CUTOFF
            if (includeInteraction) {
                //VERIFY_CUTOFF;
            	real3 pos2 = trim(posq[p2]);
            	real3 pos3 = trim(posq[p3]);
            	includeInteraction &= (delta(pos2, pos3, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ).w < CUTOFF_SQUARED);
            }
#endif
           if (includeInteraction) {
//             PERMUTE_ATOMS;
        	   int atom1 = p1;
        	   int atom2 = p2;
        	   int atom3 = p3;

//             LOAD_PARTICLE_DATA;
//             COMPUTE_INTERACTION;
        	   real computed_energy = computeInteraction(atom1, atom2, atom3, posq, &periodicBoxSize, forces);
//        	   printf("computed energy = %lf for atoms { %d, %d, %d } in thread: %d\n", computed_energy, atom1, atom2, atom3, threadIdx.x);
        	   energy += computed_energy;
        	   int oxygens[] = {atom3, atom2, atom1}; // ordered to match ref
        	   for (int j = 0, k = 0; j<3; j++){// j used to select oxygen index, k for index in force array
				   for (int i=0, atom = oxygens[j]; i<3; i++, k++) {// i used to index of each particle associated with the oxygen
					   atomicAdd(&forceBuffers[atom + i], static_cast<unsigned long long>((long long) (forces[k].x*0x100000000)));
					   atomicAdd(&forceBuffers[atom + i+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (forces[k].y*0x100000000)));
					   atomicAdd(&forceBuffers[atom + i+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (forces[k].z*0x100000000)));
				   }
        	   }
           }
        }
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}

    
       
       
       
       
       
       
       
       
       
       
