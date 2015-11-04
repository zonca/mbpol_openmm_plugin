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
#define kHH  4.247920544074333e-01 // A^(-1)
#define kOH  8.128985941165371e-01 // A^(-1)
#define kOO  3.749457984616480e-02 // A^(-1)
#define dHH_intra  1.690594510379166e+00 // A
#define dOH_intra -2.999999868517452e+00 // A
#define dHH  3.499031358429095e+00 // A
#define dOH  4.854042963514281e+00 // A
#define dOO  4.816312044947604e-08 // A

extern "C" __device__ void computeVar(real k, real r0, real3 * a1, real3 * a2, real * var)
{
    real3 dx = {(a1[0] - a2[0])*NM_TO_A,
                (a1[1] - a2[1])*NM_TO_A,
                (a1[2] - a2[2])*NM_TO_A};

    real dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    real d = SQRT(dsq);

    *var = EXP(-k*(d - r0));
}

extern "C" __device__ void computeGVar(real g, real k, real r0, real3 * a1, real3 * a2, real3 * g1, real3 * g2)
{
    real3 dx = {(a1[0] - a2[0])*NM_TO_A,
                (a1[1] - a2[1])*NM_TO_A,
                (a1[2] - a2[2])*NM_TO_A};

    real dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    real d = SQRT(dsq);

    real gg = - k*g*EXP(-k*(d - r0))/d;

    real cal2joule = 4.184;

    gg *= cal2joule * 10.*-1;

    for (int i = 0; i < 3; ++i) {
        g1[i] += gg*dx[i];
        g2[i] -= gg*dx[i];
    }
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
            tempEnergy = poly_3b_v2x_eval(exp, g);
			
			real gab, gac, gbc;
			real sab, sac, sbc;
			evaluateSwitchFunc(drab, gab, &sab);
			evaluateSwitchFunc(drac, gac, &sac);
			evaluateSwitchFunc(drbc, gbc, &sbc);
			
			real s = sab*sac + sab*sbc + sac*sbc;
				
			for (int j = 0; j < 36; ++j)
				g[n] *= s;
			
			
				
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
			
			real cal2joule = 4.184;
			
			for (int n = 0; n < 3; ++n) {
              	forces[Oa][n] += (gab*rab[n] + gac*rac[n]) * cal2joule * -NM_TO_A;
              	forces[Ob][n] += (gbc*rbc[n] - gab*rab[n]) * cal2joule * -NM_TO_A;
              	forces[Oc][n] -= (gac*rac[n] + gbc*rbc[n]) * cal2joule * -NM_TO_A;
          	 }
          	 
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
        const ushort3* __restrict__ exclusionTiles,
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
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int globtx = (blockIdx.x*blockDim.x+threadIdx.x); // global index
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex
    real energy = 0.0f;
    // used shared memory if the device cannot shuffle
    // localData contains positions and forces for all atoms, also H
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort3 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        
        const unsigned int z = tileIndices.z;
        
        real3 forces[10];
        for (int i=0; i<3; i++) {
           forces[i] = make_real3(0);
        }
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        const bool hasExclusions = true;
        if (x == y) {
            // This tile is on the diagonal.
            // Diagonal tile means that a block is interacting with itself.
            // For simplicity here OpenMM just repeats twice the calculation
            // of the interaction, i.e. atom0 with atom3 and atom3 with atom0
            // TODO improve this by only computing interactions once,
            // this would require to write results to localData and to
            // use the same scanning technique used in other tiles.
            localData[threadIdx.x].x = posq1.x;
            localData[threadIdx.x].y = posq1.y;
            localData[threadIdx.x].z = posq1.z;

            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                // second atom is always changing so need to zero out
                for (int i=3; i<10; i++) {
                   forces[i] = make_real3(0);
                }
                int atom2 = tbx+j;
                int atom3 = 0; //TODO: determine the actual value of atom3
                
                real3 posq2;
                posq2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                atom2 = y*TILE_SIZE+j;
                real dEdR = 0.0f;
                real tempEnergy = 0.0f;
                if ((atom1 % 3 == 0) && (atom2 % 3 == 0) && (NUM_ATOMS > atom2) && (atom1 < NUM_ATOMS) && (atom1 != atom2)) {
                    // this computes both atom0-atom3 and atom3-atom0
                    // COMPUTE_INTERACTION exclusions diagonal tile
                    energy += computeInteraction(atom1, atom2, atom3, posq, &periodicBoxSize, forces)/2.;
                }
            }
        }
        else {
            // This is an off-diagonal tile.
            unsigned int j = y*TILE_SIZE + tgx;
            real4 shflPosq = posq[j];
            localData[threadIdx.x].x = shflPosq.x;
            localData[threadIdx.x].y = shflPosq.y;
            localData[threadIdx.x].z = shflPosq.z;
            localData[threadIdx.x].fx = 0.0f;
            localData[threadIdx.x].fy = 0.0f;
            localData[threadIdx.x].fz = 0.0f;
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+tj;
                int atom3 = 0; //TODO: determine the actual value of atom3
                
                real3 posq2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                atom2 = y*TILE_SIZE+tj;
                real dEdR = 0.0f;
                real tempEnergy = 0.0f;
                if ((atom1 % 3 == 0) && (atom2 % 3 == 0) && (atom1 > atom2) && (atom1 < NUM_ATOMS)) {
                    // COMPUTE_INTERACTION exclusions off diagonal tile
                    // this computes only atom3-atom0
                    energy += computeInteraction(atom1, atom2, atom3, posq, &periodicBoxSize, forces);
                }
                for (int i=0; i<3; i++) {
                    localData[tbx+tj+i].fx += forces[Ob + i].x;
                    localData[tbx+tj+i].fy += forces[Ob + i].y;
                    localData[tbx+tj+i].fz += forces[Ob + i].z;
                }
                // cycles the indices
                // 0 1 2 3 4 5 6 7 -> 1 2 3 4 5 6 7 0
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            const unsigned int offset = y*TILE_SIZE + tgx;
            // write results for off diagonal tiles
            if (offset < PADDED_NUM_ATOMS) {
                atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) ((CAL2JOULE * -10 * localData[threadIdx.x].fx)*0x100000000)));
                atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) ((CAL2JOULE * -10 * localData[threadIdx.x].fy)*0x100000000)));
                atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) ((CAL2JOULE * -10 * localData[threadIdx.x].fz)*0x100000000)));
            }
        }
        // Write results for on and off diagonal tiles
        const unsigned int offset = x*TILE_SIZE + tgx;
        for (int i=0; i<3; i++) {
            forces[i] *= CAL2JOULE * -10;
        }

        // Write results.
        for (int i=0; i<3; i++) {
            atomicAdd(&forceBuffers[offset + i], static_cast<unsigned long long>((long long) (forces[i].x*0x100000000)));
            atomicAdd(&forceBuffers[offset + i+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (forces[i].y*0x100000000)));
            atomicAdd(&forceBuffers[offset + i+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (forces[i].z*0x100000000)));
        }
    }
    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).
#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    int pos = (int) (numTiles > maxTiles ? startTileIndex+warp*(long long)numTileIndices/totalWarps : warp*(long long)numTiles/totalWarps);
    int end = (int) (numTiles > maxTiles ? startTileIndex+(warp+1)*(long long)numTileIndices/totalWarps : (warp+1)*(long long)numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    // atomIndices can probably be shuffled as well
    // but it probably wouldn't make things any faster
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;

    while (pos < end) {
        real3 forces[10];
        // set only forces for current water to 0
        // forces is a variable local to the thread,
        // forces [0:3] contains the local water forces,
        // those forces are accumulated for each interaction with
        // other molecules
        // forces[4:6] contains the second water molecule that is
        // different for every interaction so we need to set it to
        // zero in the inner loop.
        // then those forces are added to the localData array
        // that is in shared memory and accumulates the forces as we
        // go through the grid of interactions.
        // need to make sure that localData has complete water molecules

        for (int i=0; i<3; i++) {
           forces[i] = make_real3(0);
        }
        bool includeTile = true;

        // Extract the coordinates of this tile.
        int x, y;
#ifdef USE_CUTOFF
        if (numTiles <= maxTiles) {
            x = tiles[pos];
        }
        else
#endif
        {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }

            // Skip over tiles that have exclusions, since they were already processed.

            while (skipTiles[tbx+TILE_SIZE-1] < pos) {
                if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                    ushort3 tile = exclusionTiles[skipBase+tgx];
                    skipTiles[threadIdx.x] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
                }
                else
                    skipTiles[threadIdx.x] = end;
                skipBase += TILE_SIZE;
                currentSkipIndex = tbx;
            }
            while (skipTiles[currentSkipIndex] < pos)
                currentSkipIndex++;
            includeTile = (skipTiles[currentSkipIndex] != pos);
        }
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;
            // Load atom data for this tile.
            real4 posq1 = posq[atom1];
            // LOAD_ATOM1_PARAMETERS
            //const unsigned int localAtomIndex = threadIdx.x;
#ifdef USE_CUTOFF
            unsigned int j = (numTiles <= maxTiles ? interactingAtoms[pos*TILE_SIZE+tgx] : y*TILE_SIZE + tgx);
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
            if (j < PADDED_NUM_ATOMS) {
                // Load position of atom j from from global memory
                localData[threadIdx.x].x = posq[j].x;
                localData[threadIdx.x].y = posq[j].y;
                localData[threadIdx.x].z = posq[j].z;
                localData[threadIdx.x].fx = 0.0f;
                localData[threadIdx.x].fy = 0.0f;
                localData[threadIdx.x].fz = 0.0f;
                // LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            else {
                localData[threadIdx.x].x = 0;
                localData[threadIdx.x].y = 0;
                localData[threadIdx.x].z = 0;
            }
            // We need to apply periodic boundary conditions separately for each interaction.
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                // second atom is always changing so need to zero out
                for (int i=3; i<10; i++) {
                   forces[i] = make_real3(0);
                }
                unsigned int atom2 = tbx+tj;
                int atom3 = 0; //TODO: determine the actual value of atom3
                
                real3 posq2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                // LOAD_ATOM2_PARAMETERS
                atom2 = atomIndices[tbx+tj];
                real dEdR = 0.0f;
                real tempEnergy = 0.0f;
                //if ((atom1 < NUM_ATOMS) && (atom2 < NUM_ATOMS))
                //tempEnergy = 1.;
                // FIXME temporary implementation of filtering O atoms.
                // here we are using the Neighbor list implementation available in NonBondedUtilities.
                // Therefore we have a Neighbor list of all atoms,
                // then in this loop we filter out only the Oxygens.
                // Better implementation would be to write our own implemenation of a O-only Neighbor
                // list based either on NonBondedUtilities or on CustomManyParticleForce
                if ((atom1 % 3 == 0) && (atom2 % 3 == 0) && (atom1 > atom2) && (atom1 < NUM_ATOMS)) {
                    // COMPUTE_INTERACTION no exclusions
                    // this computes only atom3-atom0
                    energy += computeInteraction(atom1, atom2, atom3, posq, &periodicBoxSize, forces);

                    // write forces of second molecule to shared memory

                    for (int i=0; i<3; i++) {
                        localData[tbx+tj+i].fx += forces[Ob + i].x;
                        localData[tbx+tj+i].fy += forces[Ob + i].y;
                        localData[tbx+tj+i].fz += forces[Ob + i].z;
                    }

                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }


            for (int i=0; i<3; i++) {
                forces[i] *= CAL2JOULE * -10;
            }

            // Write results.
            for (int i=0; i<3; i++) {
                atomicAdd(&forceBuffers[atom1 + i], static_cast<unsigned long long>((long long) (forces[i].x*0x100000000)));
                atomicAdd(&forceBuffers[atom1 + i+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (forces[i].y*0x100000000)));
                atomicAdd(&forceBuffers[atom1 + i+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (forces[i].z*0x100000000)));
            }
#ifdef USE_CUTOFF
            unsigned int atom2 = atomIndices[threadIdx.x];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
            if (atom2 < PADDED_NUM_ATOMS) {
                atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long) ((CAL2JOULE * -10 * localData[threadIdx.x].fx)*0x100000000)));
                atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) ((CAL2JOULE * -10 * localData[threadIdx.x].fy)*0x100000000)));
                atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) ((CAL2JOULE * -10 * localData[threadIdx.x].fz)*0x100000000)));
            }
        }
        pos++;
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;

}       
       
       
       
       
       
       
       
       
       
       