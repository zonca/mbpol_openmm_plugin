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
#define Ob  6
#define Hb1 7
#define Hb2 8

real r3i =  0.000000000000000e+00; // A
real r3f =  4.500000000000000e+00; // A

real kHH_intra = -2.254303257872797e+00; // A^(-1)
real kOH_intra = -2.564771901404151e-01; // A^(-1)

real kHH =  4.247920544074333e-01; // A^(-1)
real kOH =  8.128985941165371e-01; // A^(-1)
real kOO =  3.749457984616480e-02; // A^(-1)

real dHH_intra =  1.690594510379166e+00; // A
real dOH_intra = -2.999999868517452e+00; // A

real dHH =  3.499031358429095e+00; // A
real dOH =  4.854042963514281e+00; // A
real dOO =  4.816312044947604e-08; // A

extern "C" __device__ void computeVar(real k, real r0, real3 * a1, real3 * a2, real * var)
{
    real3 dx = {(a1[0] - a2[0])*nm_to_A,
                (a1[1] - a2[1])*nm_to_A,
                (a1[2] - a2[2])*nm_to_A};

    real dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    real d = SQRT(dsq);

    *var = EXP(-k*(d - r0));
}

extern "C" __device__ void computeGVar(real g, real k, real r0, real3 * a1, real3 * a2, real3 * g1, real3 * g2)
{
    real3 dx = {(a1[0] - a2[0])*nm_to_A,
                (a1[1] - a2[1])*nm_to_A,
                (a1[2] - a2[2])*nm_to_A};

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

extern "C" __device__ void evaluateSwitchFunc(real r, real g, real * s)
{
    if (r > r3f) {
        g = 0.0;
        *s = 0.0;
    } else if (r > r3i) {
        real t1 = M_PI/(r3f - r3i);
        real x = (r - r3i)*t1;
        g = - SIN(x)*t1/2.0;
        *s = (1.0 + COS(x))/2.0;
    } else {
        g = 0.0;
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
		
		rab = (positions[Oa] - positions[Ob])*nm_to_A;
		drab += dot(rab, rab);
			
		rac = (positions[Oa] - positions[Oc])*nm_to_A;
		drac += dot(rac, rac);
		
		rab = (positions[Ob] - positions[Oc])*nm_to_A;
		drbc += dot(drbc, drdc);
		
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
            tempEnergy = poly-3b-v2x_eval(exp, g);
			
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
              	forces[Oa][n] += (gab*rab[n] + gac*rac[n]) * cal2joule * -nm_to_A;
              	forces[Ob][n] += (gbc*rbc[n] - gab*rab[n]) * cal2joule * -nm_to_A;
              	forces[Oc][n] -= (gac*rac[n] + gbc*rbc[n]) * cal2joule * -nm_to_A;
          	 }
          	 
// Is it okay to calculate the force in the shared variable like in cuda 2body
// or should we calculate in a seperate variable and add to the actual at one 
// time at the end like in refrence code ??

          	 real energy = tempEnergy * cal2joule;
          	 
          	 return energy;
          	 
        }
       
       
       
       
       
       
       
       
       
       
       