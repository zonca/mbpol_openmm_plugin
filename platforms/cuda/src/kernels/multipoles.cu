extern "C" __global__ void recordInducedDipoles(const long long* __restrict__ fieldBuffers, const long long* __restrict__ fieldPolarBuffers,
		real* __restrict__ inducedDipole, real* __restrict__ inducedDipolePolar, const float* __restrict__ polarizability) {
	for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += gridDim.x*blockDim.x) {
		real scale = polarizability[atom]/(real) 0x100000000;
		inducedDipole[3*atom] = scale*fieldBuffers[atom];
		inducedDipole[3*atom+1] = scale*fieldBuffers[atom+PADDED_NUM_ATOMS];
		inducedDipole[3*atom+2] = scale*fieldBuffers[atom+PADDED_NUM_ATOMS*2];
		inducedDipolePolar[3*atom] = scale*fieldPolarBuffers[atom];
		inducedDipolePolar[3*atom+1] = scale*fieldPolarBuffers[atom+PADDED_NUM_ATOMS];
		inducedDipolePolar[3*atom+2] = scale*fieldPolarBuffers[atom+PADDED_NUM_ATOMS*2];
	}
}

/**
 * Normalize a vector and return what its magnitude was.
 */
inline __device__ real normVector(real3& v) {
	real n = SQRT(dot(v, v));
	v *= (n > 0 ? RECIP(n) : 0);
	return n;
}

/**
 * Compute the electrostatic potential at each of a set of points.
 */
extern "C" __global__ void computePotentialAtPoints(const real4* __restrict__ posq,
		const real* __restrict__ inducedDipole, const real4* __restrict__ points,
		real* __restrict__ potential, int numPoints, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
	extern __shared__ real4 localPosq[];
	real3* localInducedDipole = (real3*)&localPosq[blockDim.x];
	for (int basePoint = blockIdx.x*blockDim.x; basePoint < numPoints; basePoint += gridDim.x*blockDim.x) {
		int point = basePoint+threadIdx.x;
		real4 pointPos = points[point];
		real p = 0;
		for (int baseAtom = 0; baseAtom < NUM_ATOMS; baseAtom += blockDim.x) {
			int atom = baseAtom+threadIdx.x;

			// Load data into shared memory.

			if (atom < NUM_ATOMS) {
				localPosq[threadIdx.x] = posq[atom];
				localInducedDipole[threadIdx.x] = make_real3(inducedDipole[3*atom], inducedDipole[3*atom+1], inducedDipole[3*atom+2]);
			}
			__syncthreads();

			// Loop over atoms and compute the potential at this point.

			if (point < numPoints) {
				int end = min(blockDim.x, NUM_ATOMS-baseAtom);
				for (int i = 0; i < end; i++) {
					real3 delta = trimTo3(localPosq[i]-pointPos);
#ifdef USE_PERIODIC
					APPLY_PERIODIC_TO_DELTA(delta)
#endif
					real r2 = dot(delta, delta);
					real rInv = RSQRT(r2);
					p += localPosq[i].w*rInv;
					real rr2 = rInv*rInv;
					real rr3 = rInv*rr2;
					real scu = dot(localInducedDipole[i], delta);
					p -= (scu)*rr3;
				}
			}
			__syncthreads();
		}
		potential[point] = p;
	}
}
