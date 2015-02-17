extern "C" __global__ void computeTwoBodyForce(const unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ tempForceBuffers) {
    for (unsigned int atom1 = blockIdx.x*blockDim.x+threadIdx.x; atom1 < 3; atom1 += blockDim.x*gridDim.x) {
        //int atom2 = bondReductionAtoms[atom1];
        //long long fx1 = forceBuffers[atom1];
        //long long fy1 = forceBuffers[atom1+PADDED_NUM_ATOMS];
        //long long fz1 = forceBuffers[atom1+PADDED_NUM_ATOMS*2];
        //if (atom1 != atom2) {
        //    double factor = (double) bondReductionFactors[atom1];
        //    long long fx2 = (long long) ((1-factor)*fx1);
        //    long long fy2 = (long long) ((1-factor)*fy1);
        //    long long fz2 = (long long) ((1-factor)*fz1);
        //    atomicAdd(&tempForceBuffers[atom2], static_cast<unsigned long long>(fx2));
        //    atomicAdd(&tempForceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>(fy2));
        //    atomicAdd(&tempForceBuffers[atom2+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>(fz2));
        //    fx1 = (long long) (factor*fx1);
        //    fy1 = (long long) (factor*fy1);
        //    fz1 = (long long) (factor*fz1);
        //}
        atomicAdd(&tempForceBuffers[atom1], static_cast<unsigned long long>(1000.));
        atomicAdd(&tempForceBuffers[atom1+1], static_cast<unsigned long long>(2.));
    }
}
