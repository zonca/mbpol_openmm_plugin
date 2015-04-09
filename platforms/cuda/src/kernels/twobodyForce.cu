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
#define Xa1 6
#define Xa2 7
#define Xb1 8
#define Xb2 9

#define k_HH_intra -6.480884773303821e-01 // A^(-1)
#define k_OH_intra  1.674518993682975e+00 // A^(-1)
#define k_HH_coul 1.148231864355956e+00 // A^(-1)
#define k_OH_coul 1.205989761123099e+00 // A^(-1)
#define k_OO_coul 1.395357065790959e+00 // A^(-1)
#define k_XH_main 7.347036852042255e-01 // A^(-1)
#define k_XO_main 7.998249864422826e-01 // A^(-1)
#define k_XX_main 7.960663960630585e-01 // A^(-1)
#define in_plane_gamma  -9.721486914088159e-02
#define out_of_plane_gamma  9.859272078406150e-02
#define r2i 4.500000000000000e+00 // A
#define r2f 6.500000000000000e+00 // A
#define d_intra 1.0
#define d_inter 4.0

extern "C" __device__ void computeExtraPoint(real3 * O, real3 * H1, real3 * H2, real3 * X1, real3 * X2) {
    // TODO save oh1 and oh2 to be used later?
    real3 oh1 = *H1 - *O;
    real3 oh2 = *H2 - *O;

    real3 v = cross(oh1, oh2);
    real3 in_plane = (*O) + (oh1 + oh2) * 0.5 * in_plane_gamma;
    real3 out_of_plane = v * out_of_plane_gamma;

    *X1 = in_plane + out_of_plane;
    *X2 = in_plane - out_of_plane;
}

extern "C" __device__ void computeExp(real r0, real k, real3 * O1, real3 * O2, real * exp1, real3 * g) {
    *g = *O1 - *O2;

    real r = SQRT(dot(*g, *g));
    *exp1 = EXP(k*(r0 - r));
    *g *= -k * (*exp1) / r;
}

extern "C" __device__ void computeCoul(real r0, real k, real3 * O1, real3 * O2, real * val, real3 * g) {
    *g = *O1 - *O2;

    real r = SQRT(dot(*g, *g));
    real exp1 = EXP(k * (r0 - r));
    real rinv = 1.0/r;
    *val = exp1*rinv;
    *g *=  - (k + rinv) * (*val) * rinv;
}

extern "C" __device__ void computeGrads(real * g, real3 * gOO, real3 * force1, real3 * force2, real sw) {

    real3 d = *g * (*gOO);
    *force1 += sw * d;
    *force2 -= sw * d;
}

extern "C" __device__ void distributeXpointGrad(real3 * O, real3 * H1, real3 * H2, real3 * forceX1, real3 * forceX2, real3 * forceO, real3 * forceH1, real3 * forceH2, real sw) {

    // TODO save oh1 and oh2 to be used later?
    real3 oh1 = *H1 - *O;
    real3 oh2 = *H2 - *O;

    real3 gm = *forceX1-*forceX2;

    real3 t1 = cross(oh2, gm);

    real3 t2 = cross(oh1, gm);

    real3 gsum = *forceX1 + *forceX2;
    real3 in_plane = gsum*0.5*in_plane_gamma;

    real3 gh1 = in_plane + t1*out_of_plane_gamma;
    real3 gh2 = in_plane - t2*out_of_plane_gamma;

    *forceO +=  sw * (gsum - (gh1 + gh2)); // O
    *forceH1 += sw * gh1; // H1
    *forceH2 += sw * gh2; // H2

}

extern "C" __device__ void evaluateSwitchFunc(real r, real * sw, real * gsw) {

    if (r > r2f) {
        *gsw = 0.0;
        *sw  = 0.0;
    } else if (r > r2i) {
        real t1 = M_PI/(r2f - r2i);
        real x = (r - r2i)*t1;
        *gsw = - SIN(x)*t1/2.0;
        *sw  = (1.0 + COS(x))/2.0;
    } else {
        *gsw = 0.0;
        *sw = 1.0;
    }
}

extern "C" __device__ void imageMolecules(const real4 * box, real3 * allPositions)
{

    // Take first oxygen as central atom

    // image its two hydrogens with respect of the first oxygen

    imageParticles(box, &allPositions[Oa], &allPositions[Ha1]);
    imageParticles(box, &allPositions[Oa], &allPositions[Ha2]);

    // Now image the oxygen of the second &molecule

    imageParticles(box, &allPositions[Oa], &allPositions[Ob]);

    // Image the hydrogen of the second mo&lecule with respect to the oxygen of the second molecule
    imageParticles(box, &allPositions[Ob], &allPositions[Hb1]);
    imageParticles(box, &allPositions[Ob], &allPositions[Hb2]);

}

extern "C" __device__ real computeInteraction(
        const unsigned int atom1,
        const unsigned int atom2,
        const real4* __restrict__ posq,
        const real4* periodicBoxSize,
        real3 * forces) {
                    real tempEnergy = 0.0f;
                    // 2 water molecules and extra positions
                    real3 positions[10];
                    // first water
                    for (int i = 0; i < 3; i++) {
                        positions[Oa + i] = make_real3( posq[atom1+i].x * NM_TO_A,
                                                        posq[atom1+i].y * NM_TO_A,
                                                        posq[atom1+i].z * NM_TO_A);
                        positions[Ob + i] = make_real3( posq[atom2+i].x * NM_TO_A,
                                                        posq[atom2+i].y * NM_TO_A,
                                                        posq[atom2+i].z * NM_TO_A);
                    }

#if USE_PERIODIC
                    imageMolecules(periodicBoxSize, positions);
#endif


                    real3 delta = make_real3(positions[Ob].x-positions[Oa].x, positions[Ob].y-positions[Oa].y, positions[Ob].z-positions[Oa].z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real invR = RSQRT(r2);
                    real rOO = r2*invR;
                    real sw = 1.;
                    real gsw = 1.;

                    if ((rOO > r2f) || (rOO < 2.)) {
                        tempEnergy = 0.;
                    } else {

                        evaluateSwitchFunc(rOO, &sw, &gsw);

                        computeExtraPoint(positions + Oa, positions + Ha1, positions + Ha2,
                               positions + Xa1, positions + Xa2);
                        computeExtraPoint(positions + Ob, positions + Hb1, positions + Hb2,
                                positions + Xb1, positions + Xb2);

                        real exp[31];
                        real3 gOO[31];
                        int i = 0;
                        computeExp(d_intra, k_HH_intra, positions +Ha1, positions +Ha2, exp+i, gOO+i); i++;
                        computeExp(d_intra, k_HH_intra, positions +Hb1, positions +Hb2, exp+i, gOO+i); i++;
                        computeExp(d_intra, k_OH_intra, positions +Oa,  positions +Ha1, exp+i, gOO+i); i++;
                        computeExp(d_intra, k_OH_intra, positions +Oa,  positions +Ha2, exp+i, gOO+i); i++;
                        computeExp(d_intra, k_OH_intra, positions +Ob,  positions +Hb1, exp+i, gOO+i); i++;
                        computeExp(d_intra, k_OH_intra, positions +Ob,  positions +Hb2, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_HH_coul, positions +Ha1, positions +Hb1, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_HH_coul, positions +Ha1, positions +Hb2, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_HH_coul, positions +Ha2, positions +Hb1, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_HH_coul, positions +Ha2, positions +Hb2, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_OH_coul, positions +Oa,  positions +Hb1, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_OH_coul, positions +Oa,  positions +Hb2, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_OH_coul, positions +Ob,  positions +Ha1, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_OH_coul, positions +Ob,  positions +Ha2, exp+i, gOO+i); i++;
                        computeCoul(d_inter, k_OO_coul, positions +Oa,  positions +Ob , exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XH_main,  positions +Xa1, positions +Hb1, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XH_main,  positions +Xa1, positions +Hb2, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XH_main,  positions +Xa2, positions +Hb1, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XH_main,  positions +Xa2, positions +Hb2, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XH_main,  positions +Xb1, positions +Ha1, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XH_main,  positions +Xb1, positions +Ha2, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XH_main,  positions +Xb2, positions +Ha1, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XH_main,  positions +Xb2, positions +Ha2, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XO_main,  positions +Oa , positions +Xb1, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XO_main,  positions +Oa , positions +Xb2, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XO_main,  positions +Ob , positions +Xa1, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XO_main,  positions +Ob , positions +Xa2, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XX_main,  positions +Xa1, positions +Xb1, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XX_main,  positions +Xa1, positions +Xb2, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XX_main,  positions +Xa2, positions +Xb1, exp+i, gOO+i); i++;
                        computeExp(d_inter, k_XX_main,  positions +Xa2, positions +Xb2, exp+i, gOO+i); i++;

                        real g[31];
                        tempEnergy = poly_2b_v6x_eval(exp, g);

                        computeGrads(g+0,  gOO+0,  forces + Ha1, forces + Ha2, sw);
                        computeGrads(g+1,  gOO+1,  forces + Hb1, forces + Hb2, sw);
                        computeGrads(g+2,  gOO+2,  forces + Oa , forces + Ha1, sw);
                        computeGrads(g+3,  gOO+3,  forces + Oa , forces + Ha2, sw);
                        computeGrads(g+4,  gOO+4,  forces + Ob , forces + Hb1, sw);
                        computeGrads(g+5,  gOO+5,  forces + Ob , forces + Hb2, sw);
                        computeGrads(g+6,  gOO+6,  forces + Ha1, forces + Hb1, sw);
                        computeGrads(g+7,  gOO+7,  forces + Ha1, forces + Hb2, sw);
                        computeGrads(g+8,  gOO+8,  forces + Ha2, forces + Hb1, sw);
                        computeGrads(g+9,  gOO+9,  forces + Ha2, forces + Hb2, sw);
                        computeGrads(g+10, gOO+10, forces + Oa , forces + Hb1, sw);
                        computeGrads(g+11, gOO+11, forces + Oa , forces + Hb2, sw);
                        computeGrads(g+12, gOO+12, forces + Ob , forces + Ha1, sw);
                        computeGrads(g+13, gOO+13, forces + Ob , forces + Ha2, sw);
                        computeGrads(g+14, gOO+14, forces + Oa , forces + Ob , sw);
                        computeGrads(g+15, gOO+15, forces + Xa1, forces + Hb1, sw);
                        computeGrads(g+16, gOO+16, forces + Xa1, forces + Hb2, sw);
                        computeGrads(g+17, gOO+17, forces + Xa2, forces + Hb1, sw);
                        computeGrads(g+18, gOO+18, forces + Xa2, forces + Hb2, sw);
                        computeGrads(g+19, gOO+19, forces + Xb1, forces + Ha1, sw);
                        computeGrads(g+20, gOO+20, forces + Xb1, forces + Ha2, sw);
                        computeGrads(g+21, gOO+21, forces + Xb2, forces + Ha1, sw);
                        computeGrads(g+22, gOO+22, forces + Xb2, forces + Ha2, sw);
                        computeGrads(g+23, gOO+23, forces + Oa , forces + Xb1, sw);
                        computeGrads(g+24, gOO+24, forces + Oa , forces + Xb2, sw);
                        computeGrads(g+25, gOO+25, forces + Ob , forces + Xa1, sw);
                        computeGrads(g+26, gOO+26, forces + Ob , forces + Xa2, sw);
                        computeGrads(g+27, gOO+27, forces + Xa1, forces + Xb1, sw);
                        computeGrads(g+28, gOO+28, forces + Xa1, forces + Xb2, sw);
                        computeGrads(g+29, gOO+29, forces + Xa2, forces + Xb1, sw);
                        computeGrads(g+30, gOO+30, forces + Xa2, forces + Xb2, sw);


                    distributeXpointGrad(positions + Oa, positions + Ha1, positions + Ha2,
                            forces + Xa1, forces + Xa2,
                            forces + Oa, forces + Ha1, forces + Ha2, sw);

                    distributeXpointGrad(positions + Ob, positions + Hb1, positions + Hb2,
                            forces + Xb1, forces + Xb2,
                            forces + Ob, forces + Hb1, forces + Hb2, sw);

                    }

                    // gradient of the switch
                    gsw *= tempEnergy/rOO;
                    real3 d = gsw * delta;
                    forces[Oa] += d;
                    forces[Ob] -= d;

                    return sw * tempEnergy * CAL2JOULE;
}

/**
 * Compute nonbonded interactions. The kernel is separated into two parts,
 * tiles with exclusions and tiles without exclusions. It relies heavily on
 * implicit warp-level synchronization. A tile is defined by two atom blocks
 * each of warpsize. Each warp computes a range of tiles.
 *
 * Tiles with exclusions compute the entire set of interactions across
 * atom blocks, equal to warpsize*warpsize. In order to avoid access conflicts
 * the forces are computed and accumulated diagonally in the manner shown below
 * where, suppose
 *
 * [a-h] comprise atom block 1, [i-p] comprise atom block 2
 *
 * 1 denotes the first set of calculations within the warp
 * 2 denotes the second set of calculations within the warp
 * ... etc.
 *
 *        threads
 *     0 1 2 3 4 5 6 7
 *         atom1
 * L    a b c d e f g h
 * o  i 1 2 3 4 5 6 7 8
 * c  j 8 1 2 3 4 5 6 7
 * a  k 7 8 1 2 3 4 5 6
 * l  l 6 7 8 1 2 3 4 5
 * D  m 5 6 7 8 1 2 3 4
 * a  n 4 5 6 7 8 1 2 3
 * t  o 3 4 5 6 7 8 1 2
 * a  p 2 3 4 5 6 7 8 1
 *
 * AZ: so thread 0 is responsible for atom a, and accumulates the force in the `force` local
 * variable, only at the end this is written to the force buffer.
 * then at each step it has exclusive access to one of the index in localData in the order
 * detailed above.
 * For MBPol, the difference is that we work with a full molecule at a time, so now the local\
 * force variable is instead an array of 10 vectors, 3 for each molecule and 2 extras.
 * The first 3 components act exactly like `force`, the rest is used temporarily and then
 * copied to localData at each step.
 *
 *
 * Tiles without exclusions read off directly from the neighbourlist interactingAtoms
 * and follows the same force accumulation method. If more there are more interactingTiles
 * than the size of the neighbourlist initially allocated, the neighbourlist is rebuilt
 * and the full tileset is computed. This should happen on the first step, and very rarely
 * afterwards.
 *
 * On CUDA devices that support the shuffle intrinsic, on diagonal exclusion tiles use
 * __shfl to broadcast. For all other types of tiles __shfl is used to pass around the
 * forces, positions, and parameters when computing the forces.
 *
 * [out]forceBuffers    - forces on each atom to eventually be accumulated
 * [out]energyBuffer    - energyBuffer to eventually be accumulated
 * [in]posq             - x,y,z,charge
 * [in]exclusions       - 1024-bit flags denoting atom-atom exclusions for each tile
 * [in]exclusionTiles   - x,y denotes the indices of tiles that have an exclusion
 * [in]startTileIndex   - index into first tile to be processed
 * [in]numTileIndices   - number of tiles this context is responsible for processing
 * [in]int tiles        - the atom block for each tile
 * [in]interactionCount - total number of tiles that have an interaction
 * [in]maxTiles         - stores the size of the neighbourlist in case it needs
 *                      - to be expanded
 * [in]periodicBoxSize  - size of the Periodic Box, last dimension (w) not used
 * [in]invPeriodicBox   - inverse of the periodicBoxSize, pre-computed for speed
 * [in]blockCenter      - the center of each block in euclidean coordinates
 * [in]blockSize        - size of the each block, radiating from the center
 *                      - x is half the distance of total length
 *                      - y is half the distance of total width
 *                      - z is half the distance of total height
 *                      - w is not used
 * [in]interactingAtoms - a list of interactions within a given tile
 *
 */

extern "C" __global__ void computeTwoBodyForce(

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
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
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
                real3 posq2;
                posq2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                atom2 = y*TILE_SIZE+j;
                real dEdR = 0.0f;
                real tempEnergy = 0.0f;
                if ((atom1 % 3 == 0) && (atom2 % 3 == 0) && (NUM_ATOMS > atom2) && (atom1 < NUM_ATOMS) && (atom1 != atom2)) {
                    // this computes both atom0-atom3 and atom3-atom0
                    // COMPUTE_INTERACTION exclusions diagonal tile
                    energy += computeInteraction(atom1, atom2, posq, &periodicBoxSize, forces)/2.;
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
                real3 posq2 = make_real3(localData[atom2].x, localData[atom2].y, localData[atom2].z);
                atom2 = y*TILE_SIZE+tj;
                real dEdR = 0.0f;
                real tempEnergy = 0.0f;
                if ((atom1 % 3 == 0) && (atom2 % 3 == 0) && (atom1 > atom2) && (atom1 < NUM_ATOMS)) {
                    // COMPUTE_INTERACTION exclusions off diagonal tile
                    // this computes only atom3-atom0
                    energy += computeInteraction(atom1, atom2, posq, &periodicBoxSize, forces);
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
                    ushort2 tile = exclusionTiles[skipBase+tgx];
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
                    energy += computeInteraction(atom1, atom2, posq, &periodicBoxSize, forces);

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

