////////////////////////////////////////////////////////////////////////////////

__device__ const real EPS = 2.2204460492503131E-16; //std::numeric_limits<real>::epsilon();
__device__ const real FPMIN = 2.2250738585072014e-308/EPS; //std::numeric_limits<real>::min()/EPS;

const int ngau = 18;

__device__ const real y[18] = {0.0021695375159141994,
    0.011413521097787704,0.027972308950302116,0.051727015600492421,
    0.082502225484340941, 0.12007019910960293,0.16415283300752470,
    0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
    0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
    0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
    0.87126389619061517, 0.95698180152629142};

__device__ const real w[18] = {0.0055657196642445571,
    0.012915947284065419,0.020181515297735382,0.027298621498568734,
    0.034213810770299537,0.040875750923643261,0.047235083490265582,
    0.053244713977759692,0.058860144245324798,0.064039797355015485,
    0.068745323835736408,0.072941885005653087,0.076598410645870640,
    0.079687828912071670,0.082187266704339706,0.084078218979661945,
    0.085346685739338721,0.085983275670394821};

////////////////////////////////////////////////////////////////////////////////

__device__ void gammln_device(const real& xx, real& retval){
    const real cof[14] = {57.1562356658629235,-59.5979603554754912,
        14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
        .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
        -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
        .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};

    assert(xx > 0.0);

    real x = xx;
    real y = xx;

    real tmp = x + 5.24218750000000000;
    tmp = (x + 0.5)*std::log(tmp) - tmp;

    real ser = 0.999999999999997092;
    for (int j = 0; j < 14; ++j)
        ser += cof[j]/++y;

    retval = tmp + std::log(2.5066282746310005*ser/x);
}


///////////////////////////////////////////////////////////////////////////////

__device__ real gammpapprox(const real& a, const real& x, int psig)
{
    real gln;
    gammln_device(a, gln);

    const real a1 = a-1.0;
    const real lna1 = log(a1);
    const real sqrta1 = sqrt(a1);

    real xu, t, sum, ans;

    if (x > a1)
        xu = fmax(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
    else
        xu = fmax(0.0, fmin(a1 - 7.5*sqrta1, x - 5.0*sqrta1));

    sum = 0.0;
    for (int j = 0; j < ngau; ++j) {
        t = x + (xu - x)*y[j];
        sum += w[j]*std::exp(-(t - a1) + a1*(std::log(t) - lna1));
    }

    ans = sum*(xu - x)*std::exp(a1*(lna1 - 1.0) - gln);

    return (psig ? (ans > 0.0 ? 1.0 - ans : -ans)
            : (ans >= 0.0 ? ans : 1.0 + ans));
}

////////////////////////////////////////////////////////////////////////////////

__device__ real gser(const real& a, const real& x) {
    real gln;
    gammln_device(a, gln);


    real ap = a;
    real sum = 1.0/a;
    real del = sum;

    for (int i=0; i<100; i++) {
        ++ap;
        del *= x/ap;
        sum += del;
        if (fabs(del) < fabs(sum)*EPS) {
            break;
        }
    }
    return sum*exp(- x + a*log(x) - gln);
}
////////////////////////////////////////////////////////////////////////////////

__device__ real gcf(const real& a, const real& x) {
    real gln;
    gammln_device(a, gln);

    real b = x + 1.0 - a;
    real c = 1.0/FPMIN;
    real d = 1.0/b;
    real h = d;

    for (int i = 1;i<100; ++i) {
        const real an = -i*(i - a);

        b += 2.0;
        d = an*d + b;
        if (fabs(d) < FPMIN)
            d = FPMIN;

        c = b + an/c;

        if (fabs(c) < FPMIN)
            c = FPMIN;

        d = 1.0/d;
        const real del = d*c;
        h *= del;

        if (fabs(del - 1.0) <= EPS)
            break;
    }

    return exp( - x + a*std::log(x) - gln)*h;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
__device__ real gammq(const real& a, const real &x) {
     const int ASWITCH = 100;

     if (!(x >= 0.0 && a > 0.0)) {
         printf("gammq: x = %lf, a = %lf\n", x, a);
     }

     assert(x >= 0.0 && a > 0.0);

     if (x == 0.0)
      return 1.0;
     else if (int(a) >= ASWITCH)
         return gammpapprox(a, x, 0);
     else if (x < a + 1.0)
         return 1.0 - gser(a,x);
     else
         return gcf(a,x);
}
__device__ void computeOneInteractionF1(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real& energy, real3& outputForce) {

    // FIXME thole copy in unique location
    const enum TholeIndices { TCC, TCD, TDD, TDDOH, TDDHH };
    const float thole[5] =  { 0.4, 0.4, 0.4,   0.4,   0.4 };
	// thole[TDD] = 0.055;
	// thole[TDDOH] = 0.626;
	// thole[TDDHH] = 0.055;

	// deltas
    real3 delta = make_real3(atom2.posq.x - atom1.posq.x, atom2.posq.y - atom1.posq.y, atom2.posq.z - atom1.posq.z);

    real r2 = dot(delta, delta);
    real r = SQRT(r2);
    real rr1 = RECIP(r);
    real rr2 = rr1*rr1;
    real rr3 = rr1*rr2;
    real rr5 = 3*rr3*rr2;
    real rr7 = 5*rr5*rr2;

    real sci3 = dot(atom1.inducedDipole, delta);
    real sci4 = dot(atom2.inducedDipole, delta);

    real scip2 = dot(atom1.inducedDipole,atom2.inducedDipolePolar) +
                 dot(atom2.inducedDipole,atom1.inducedDipolePolar);

    real scip3 = dot(atom1.inducedDipolePolar, delta);

    real scip4 = dot(atom2.inducedDipolePolar, delta);

    real gli1 = atom2.posq.w*sci3 - atom1.posq.w*sci4;

    real glip1 = atom2.posq.w*scip3 - atom1.posq.w*scip4;

    real gl0 = atom1.posq.w*atom2.posq.w;

    // if isSameWater set gl0 to zero
    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    gl0  *= !isSameWater;
    gli1 *= !isSameWater;
    glip1 *= !isSameWater;

    // real scale1CC = getAndScaleInverseRs( particleI, particleK, r, true, 1, TCC);
    // real scale3CD = getAndScaleInverseRs( particleI, particleK, r, true, 3, TCD);

    float damp      = POW(atom1.damp*atom2.damp, 1.0f/6.0f); // AA in MBPol

    real do_scaling = (damp != 0.0) & ( damp > -50.0 ); // damp or not

    real ratio       = POW(r/damp, 4.); // rA4 in MBPol

    real dampForExpCC = -1 * thole[TCC] * ratio;
    // EXP(ttm::gammln(3.0/4.0)) = 1.2254167024651776
    real scale3CC = ( 1.0 - do_scaling*EXP(dampForExpCC) ); // needed for force
    real scale1CC = scale3CC + do_scaling*(POW(thole[TCC], 1.0f/4.0f)*(r/damp)*1.2254167024651776*gammq(3.0/4.0, -dampForExpCC));

    real dampForExpCD = -1 * thole[TCD] * ratio;
    real scale3CD = ( 1.0 - do_scaling*EXP(dampForExpCD) );

    real em = rr1 * gl0 * scale1CC;
    real ei = 0.5f * gli1 * rr3 * scale3CD;
    energy = em+ei;

    if ((atom1.moleculeIndex == 0) & (atom2.moleculeIndex == 1) & (atom1.atomType==0) & (atom2.atomType==1))
    {
        printf("%d\n", isSameWater);
        printf("%.10f\n", scale1CC);
    }
    energy *= 138.9354558456;

    // RealOpenMM scale3CC = getAndScaleInverseRs( particleI, particleK, r, true, 3, TCC);
    // RealOpenMM scale5CD = getAndScaleInverseRs( particleI, particleK, r, true, 5, TCD);
    // RealOpenMM scale5DD = getAndScaleInverseRs( particleI, particleK, r, true, 5, TDD);
    // RealOpenMM scale7DD = getAndScaleInverseRs( particleI, particleK, r, true, 7, TDD);

    real scale5CD = scale3CD - do_scaling * (4./3.) * thole[TCD] * EXP(dampForExpCD) * ratio;

    int tdd = TDD;
    if (isSameWater) {
        if ((atom1.atomType == 0) | (atom2.atomType == 0)) { // one is oxygen
            tdd = TDDOH;
        } else { // both hydrogens
            tdd = TDDHH;
        }
    }
    real dampForExpDD = thole[tdd];
    real scale5DD =  1.0 - do_scaling * EXP(dampForExpDD) *  (1. + (4./3.) * thole[tdd] * ratio);
    real scale7DD = scale5DD - do_scaling * ((4./15.) * thole[tdd] * (4. * thole[tdd] * ratio - 1.) * EXP(dampForExpDD) / POW(damp, 4.0f) * POW(r, 4));

    real gf1 = rr3*gl0*scale3CC;

    real3 ftm2 = gf1*delta;

    // intermediate variables for the induced components

    real gfi1 = 0.5 * rr5 * gli1 * scale5CD +  // charge - induced dipole
                0.5 * rr5 * glip1* scale5CD  +  // charge - induced dipole
                0.5 * rr5 * scip2* scale5DD  +  //induced dipole - induced dipole
              - 0.5 * rr7 * (sci3*scip4+scip3*sci4) *scale7DD; // induced dipole - induced dipole

    real3 ftm2i = gfi1*delta;

    ftm2i += ( atom1.inducedDipolePolar *  sci4 + // iPdipole_i * idipole_k
               atom1.inducedDipole      * scip4 +
               atom2.inducedDipolePolar *  sci3 + // iPdipole_k * idipole_i
               atom2.inducedDipole * scip3  ) * 0.5 * rr5 * scale5DD;

    // Same water atoms have no induced-dipole/charge interaction
    if (not( isSameWater )) {
    ftm2i += (
                   ( atom1.inducedDipole +
             atom1.inducedDipolePolar )*-atom2.posq.w +
                   ( atom2.inducedDipole +
             atom2.inducedDipolePolar )* atom1.posq.w
                 ) * 0.5 * rr3 * scale3CD;
    }

#ifdef DIRECT_POLARIZATION
    // real gfd = 0.5*(3*rr2*scip2*scale3i - 5*rr2*(scip3*sci4+sci3*scip4)*scale5i);
    // real temp5 = 0.5*scale5i;
    // real fdir_0 = gfd*xr + temp5*(sci4*atom1.inducedDipolePolar.x + scip4*atom1.inducedDipole.x + sci3*atom2.inducedDipolePolar.x + scip3*atom2.inducedDipole.x);
    // real fdir_1 = gfd*yr + temp5*(sci4*atom1.inducedDipolePolar.y + scip4*atom1.inducedDipole.y + sci3*atom2.inducedDipolePolar.y + scip3*atom2.inducedDipole.y);
    // real fdir_2 = gfd*zr + temp5*(sci4*atom1.inducedDipolePolar.z + scip4*atom1.inducedDipole.z + sci3*atom2.inducedDipolePolar.z + scip3*atom2.inducedDipole.z);
    // ftm2i_0 -= fdir_0;
    // ftm2i_1 -= fdir_1;
    // ftm2i_2 -= fdir_2;
#else
#endif

    outputForce = -(ftm2+ftm2i);
}
