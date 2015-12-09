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

    for (;;) {
        ++ap;
        del *= x/ap;
        sum += del;
        if (fabs(del) < fabs(sum)*EPS) {
            return sum*exp(- x + a*log(x) - gln);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

__device__ real gcf(const real& a, const real& x) {
    real gln;
    gammln_device(a, gln);

    real b = x + 1.0 - a;
    real c = 1.0/FPMIN;
    real d = 1.0/b;
    real h = d;

    for (int i = 1;; ++i) {
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
    const float uScale = 1;
    real ddsc3_0 = 0;
    real ddsc3_1 = 0;
    real ddsc3_2 = 0;

    real ddsc5_0 = 0;
    real ddsc5_1 = 0;
    real ddsc5_2 = 0;

    real ddsc7_0 = 0;
    real ddsc7_1 = 0;
    real ddsc7_2 = 0;
    
	// deltas 
    real xr = atom2.posq.x - atom1.posq.x;
    real yr = atom2.posq.y - atom1.posq.y;
    real zr = atom2.posq.z - atom1.posq.z;
    
    real r2 = xr*xr + yr*yr + zr*zr;
    real r = SQRT(r2);
    real rr1 = RECIP(r);
    real rr2 = rr1*rr1;
    real rr3 = rr1*rr2;
    real rr5 = 3*rr3*rr2;
    real rr7 = 5*rr5*rr2;

    real scale3 = 1;
    real scale5 = 1;
    real scale7 = 1;
    
	// get and scale move to diffrent kernel func
    real pdamp = atom1.damp*atom2.damp;
    if (pdamp != 0) {
   
        real ratio = r/pdamp;
        //FIXME: 
        float thole = 0.4;
        float pGamma = thole;

        real damp = ratio*ratio*ratio*pGamma;
        real dampExp = EXP(-damp);
        real damp1 = damp + 1;
        real damp2 = damp*damp;

        scale3 = 1 - dampExp;
        scale5 = 1 - damp1*dampExp;
        scale7 = 1 - (damp1 + 0.6f*damp2)*dampExp;

        real factor = 3*damp*dampExp*rr2;
        real factor7 = -0.2f + 0.6f*damp;
        
        ddsc3_0 = factor*xr;
        ddsc5_0 = ddsc3_0*damp;
        ddsc7_0 = ddsc5_0*factor7;

        ddsc3_1 = factor*yr;
        ddsc5_1 = ddsc3_1*damp;
        ddsc7_1 = ddsc5_1*factor7;

        ddsc3_2 = factor*zr;
        ddsc5_2 = ddsc3_2*damp;
        ddsc7_2 = ddsc5_2*factor7;


    }
      
    real scale3i = rr3*scale3*uScale;
    real scale5i = rr5*scale5*uScale;

    real dsc3 = rr3*scale3*dScale;
    real psc3 = rr3*scale3*pScale;

    real dsc5 = rr5*scale5*dScale;
    real psc5 = rr5*scale5*pScale;

    real dsc7 = rr7*scale7*dScale;
    real psc7 = rr7*scale7*pScale;

    real sc2 = atom1.dipole.x*atom2.dipole.x + atom1.dipole.y*atom2.dipole.y + atom1.dipole.z*atom2.dipole.z;

    real sc4 = atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr;

    real sc3 = atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr;
    
    real sci3 = atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr;
    real sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    
    real scip1 = atom1.inducedDipolePolar.x*atom2.dipole.x + atom1.inducedDipolePolar.y*atom2.dipole.y + atom1.inducedDipolePolar.z*atom2.dipole.z +
                 atom2.inducedDipolePolar.x*atom1.dipole.x + atom2.inducedDipolePolar.y*atom1.dipole.y + atom2.inducedDipolePolar.z*atom1.dipole.z;

    real scip2 = atom1.inducedDipole.x*atom2.inducedDipolePolar.x + atom1.inducedDipole.y*atom2.inducedDipolePolar.y + atom1.inducedDipole.z*atom2.inducedDipolePolar.z +
                 atom2.inducedDipole.x*atom1.inducedDipolePolar.x + atom2.inducedDipole.y*atom1.inducedDipolePolar.y + atom2.inducedDipole.z*atom1.inducedDipolePolar.z;
    
    real scip3 = ((atom1.inducedDipolePolar.x)*(xr) + (atom1.inducedDipolePolar.y)*(yr) + (atom1.inducedDipolePolar.z)*(zr));

    real scip4 = ((atom2.inducedDipolePolar.x)*(xr) + (atom2.inducedDipolePolar.y)*(yr) + (atom2.inducedDipolePolar.z)*(zr));

    real gli1 = atom2.posq.w*sci3 - atom1.posq.w*sci4;
    
    real glip1 = atom2.posq.w*scip3 - atom1.posq.w*scip4;
    real glip6 = scip1;
    real gli2 = -sc3*sci4 - sci3*sc4;
    
    real glip2 = -sc3*scip4 - scip3*sc4;
    real factor3 = rr3*((gli1)*pScale + (glip1  + glip6)*dScale);
    real factor5 = rr5*(gli2*pScale + glip2*dScale);
    
    real ftm2i_0 = -0.5f*(factor3*ddsc3_0 + factor5*ddsc5_0);
    real ftm2i_1 = -0.5f*(factor3*ddsc3_1 + factor5*ddsc5_1);
    real ftm2i_2 = -0.5f*(factor3*ddsc3_2 + factor5*ddsc5_2);
      
    real gl0 = atom1.posq.w*atom2.posq.w;
    real gl1 = atom2.posq.w*sc3 - atom1.posq.w*sc4;
    real gl2 = -sc3*sc4;
    real gl6 = sc2;
    
    // if isSameWater set gl0 to zero
    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    gl0  *= !isSameWater;
    gli1 *= !isSameWater;

    // real scale1CC = getAndScaleInverseRs( particleI, particleK, r, true, 1, TCC);
    // real scale3CD = getAndScaleInverseRs( particleI, particleK, r, true, 3, TCD);

    real damp      = pow(atom1.damp*atom2.damp, 1.0f/6.0f); // AA in MBPol

    real do_scaling = (damp != 0.0) & ( damp > -50.0 ); // damp or not

    real ratio       = pow(r/damp, 4); // rA4 in MBPol
    real pgamma = thole[TCC];
    real dampForExp = -1 * pgamma * ratio;

    real scale3CD = ( 1.0 - do_scaling*EXP(dampForExp) );

    // EXP(ttm::gammln(3.0/4.0)) = 1.2254167024651776
    real scale1CC = scale3CD + do_scaling*(pow(pgamma, 1.0/4.0)*(r/damp)*1.2254167024651776*gammq(3.0/4.0, -dampForExp));
    // end getAndScaleInverseRs

    real em = rr1 * gl0 * scale1CC;
    real ei = 0.5f * gli1 * rr3 * scale3CD;
    energy = em+ei;

    energy *= 138.9354558456;

    real gf1 = rr3*gl0 + rr5*(gl1+gl6) + rr7*gl2;
    real gf2 = -atom2.posq.w*rr3 + sc4*rr5;
    real gf3 =  atom1.posq.w*rr3 + sc3*rr5;
    
    real ftm2_0 = mScale*(gf1*xr + gf2*atom1.dipole.x + gf3*atom2.dipole.x);
    real ftm2_1 = mScale*(gf1*yr + gf2*atom1.dipole.y + gf3*atom2.dipole.y);
    real ftm2_2 = mScale*(gf1*zr + gf2*atom1.dipole.z + gf3*atom2.dipole.z);

    real gfi1 = rr2*(1.5f*((gli1)*psc3 + (glip1+glip6)*dsc3 + scip2*scale3i) + 2.5f*(gli2*psc5 + glip2*dsc5 - (sci3*scip4+scip3*sci4)*scale5i));
    ftm2i_0 += gfi1*xr;
    ftm2i_1 += gfi1*yr;
    ftm2i_2 += gfi1*zr;

    ftm2i_0 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.x*psc3 + atom1.inducedDipolePolar.x*dsc3) +
               sc4*(atom1.inducedDipole.x*psc5 + atom1.inducedDipolePolar.x*dsc5)) +
      
               0.5f*(atom1.posq.w*(atom2.inducedDipole.x*psc3+atom2.inducedDipolePolar.x*dsc3) +
               sc3*(atom2.inducedDipole.x*psc5 +atom2.inducedDipolePolar.x*dsc5)) +

               scale5i*(sci4*atom1.inducedDipolePolar.x+scip4*atom1.inducedDipole.x +
                        sci3*atom2.inducedDipolePolar.x+scip3*atom2.inducedDipole.x)*0.5f +
      
               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.x +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.x;
      
    ftm2i_1 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.y*psc3 + atom1.inducedDipolePolar.y*dsc3) +
               sc4*(atom1.inducedDipole.y*psc5 + atom1.inducedDipolePolar.y*dsc5)) +

               (atom1.posq.w*(atom2.inducedDipole.y*psc3+atom2.inducedDipolePolar.y*dsc3) +
                    sc3*(atom2.inducedDipole.y*psc5+atom2.inducedDipolePolar.y*dsc5))*0.5f +
                    scale5i*(sci4*atom1.inducedDipolePolar.y+scip4*atom1.inducedDipole.y + sci3*atom2.inducedDipolePolar.y+scip3*atom2.inducedDipole.y)*0.5f +

               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.y +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.y;
      
    ftm2i_2 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.z*psc3 + atom1.inducedDipolePolar.z*dsc3) +
               sc4*(atom1.inducedDipole.z*psc5 + atom1.inducedDipolePolar.z*dsc5)) +

               (atom1.posq.w*(atom2.inducedDipole.z*psc3+atom2.inducedDipolePolar.z*dsc3) +
                    sc3*(atom2.inducedDipole.z*psc5+atom2.inducedDipolePolar.z*dsc5))*0.5f +
                    scale5i*(sci4*atom1.inducedDipolePolar.z+scip4*atom1.inducedDipole.z +
                    sci3*atom2.inducedDipolePolar.z+scip3*atom2.inducedDipole.z)*0.5f +

               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.z +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.z;

#ifdef DIRECT_POLARIZATION
    real gfd = 0.5*(3*rr2*scip2*scale3i - 5*rr2*(scip3*sci4+sci3*scip4)*scale5i);
    real temp5 = 0.5*scale5i;
    real fdir_0 = gfd*xr + temp5*(sci4*atom1.inducedDipolePolar.x + scip4*atom1.inducedDipole.x + sci3*atom2.inducedDipolePolar.x + scip3*atom2.inducedDipole.x);
    real fdir_1 = gfd*yr + temp5*(sci4*atom1.inducedDipolePolar.y + scip4*atom1.inducedDipole.y + sci3*atom2.inducedDipolePolar.y + scip3*atom2.inducedDipole.y);
    real fdir_2 = gfd*zr + temp5*(sci4*atom1.inducedDipolePolar.z + scip4*atom1.inducedDipole.z + sci3*atom2.inducedDipolePolar.z + scip3*atom2.inducedDipole.z);
    ftm2i_0 -= fdir_0;
    ftm2i_1 -= fdir_1;
    ftm2i_2 -= fdir_2;
#else
    real scaleF = 0.5f*uScale;
    real inducedFactor3 = scip2*rr3*scaleF;
    real inducedFactor5 = (sci3*scip4+scip3*sci4)*rr5*scaleF;
    real findmp_0 = inducedFactor3*ddsc3_0 - inducedFactor5*ddsc5_0;
    real findmp_1 = inducedFactor3*ddsc3_1 - inducedFactor5*ddsc5_1;
    real findmp_2 = inducedFactor3*ddsc3_2 - inducedFactor5*ddsc5_2;
    ftm2i_0 -= findmp_0;
    ftm2i_1 -= findmp_1;
    ftm2i_2 -= findmp_2;
#endif

    outputForce.x = -(ftm2_0+ftm2i_0);
    outputForce.y = -(ftm2_1+ftm2i_1);
    outputForce.z = -(ftm2_2+ftm2i_2);
}
