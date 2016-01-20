__device__ void
computeOneInteractionF1(
        AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, real bn5, float forceFactor,
        float dScale, float pScale, float mScale,
        real3& force, real& energy) {
    // FIXME thole copy in unique location
    const enum TholeIndices { TCC, TCD, TDD, TDDOH, TDDHH };
    const float thole[5] =  { 0.4, 0.4, 0.4,   0.4,   0.4 };
	// thole[TDD] = 0.055;
	// thole[TDDOH] = 0.626;
	// thole[TDDHH] = 0.055;
    real xr = delta.x;
    real yr = delta.y;
    real zr = delta.z;
    real rr1 = delta.w;
    real r = SQRT(xr*xr + yr*yr + zr*zr);

    // set the permanent multipole and induced dipole values;

    real ci = atom1.q;

    real ck = atom2.q;

    real bn1 = bn.x;
    real bn2 = bn.y;
    real bn3 = bn.z;
    real bn4 = bn.w;

    real rr3 = rr1*rr1*rr1;
    real3 ftm2 = make_real3(0);

    // calculate the scalar products for permanent components

    real gl0 = ci*ck;

    // RealOpenMM scale1CC =getAndScaleInverseRs(particleI,particleJ,r,true,1,TCC);
    // RealOpenMM scale3CD =getAndScaleInverseRs(particleI,particleJ,r,true,3,TCD);

    float damp      = POW(atom1.damp*atom2.damp, 1.0f/6.0f); // AA in MBPol

    bool do_scaling = (damp != 0.0) & ( damp > -50.0 );

    real ratio       = POW(r/damp, 4.); // rA4 in MBPol

    real dampForExpCC = -1 * thole[TCC] * ratio;
    // EXP(ttm::gammln(3.0/4.0)) = 1.2254167024651776
    real scale3CC = 1.0;
    if (do_scaling)
        scale3CC -= EXP(dampForExpCC); // needed for force
    real scale1CC = scale3CC;
    if (do_scaling)
        scale1CC += POW(thole[TCC], 1.0f/4.0f)*(r/damp)*1.2254167024651776*gammq(3.0/4.0, -dampForExpCC);

    real dampForExpCD = -1 * thole[TCD] * ratio;
    real scale3CD = 1.0;
    if (do_scaling)
        scale3CD -= do_scaling*EXP(dampForExpCD);

    // in PME same water interactions are not excluded,
    // but the scale factors are set to 0.

    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    scale1CC *= !isSameWater;
    scale3CD *= !isSameWater;
    scale3CC *= !isSameWater;

    energy += -forceFactor*(rr1*gl0*(1 - scale1CC));
    if ((atom1.moleculeIndex ==0) & (atom1.atomType == 0) & (abs(atom2.pos.x-50+0.0621) < 0.001))
    //if ((atom1.moleculeIndex ==0) & (atom1.atomType == 0) & (atom2.moleculeIndex ==0) & (atom2.atomType == 1))
    {
        printf("atom1.pos.x: %.8g\n", atom1.pos.x-50);
        printf("atom2.pos.x: %.8g\n", atom2.pos.x-50);
        printf("gl0: %.8g\n", gl0);
        printf("scale1CC: %.8g\n", scale1CC);
        printf("erl: %.8g\n", forceFactor*(rr1*gl0*(1 - scale1CC))/4.184*ENERGY_SCALE_FACTOR);
    }

    real3 delta3 = trimTo3(delta);
    real gf1 = bn1*gl0;
    real offset = 1.;
    gf1 -= (1 - scale3CC) * rr3 * gl0;
    ftm2 += gf1*delta3;

    force = ftm2;
}


__device__ void
computeOneInteractionF2(
        AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, float forceFactor,
        float dScale, float pScale, float mScale,
        real3& force, real& energy) {
    // FIXME thole copy in unique location
    const enum TholeIndices { TCC, TCD, TDD, TDDOH, TDDHH };
    const float thole[5] =  { 0.4, 0.4, 0.4,   0.4,   0.4 };
	// thole[TDD] = 0.055;
	// thole[TDDOH] = 0.626;
	// thole[TDDHH] = 0.055;
    const float uScale = 1;
    real xr = delta.x;
    real yr = delta.y;
    real zr = delta.z;
    real3 delta3 = trimTo3(delta);
    real rr1 = delta.w;
    real psc3 = 1.;
    real dsc3 = 1.;
    real usc5 = 1.;
    real usc3 = 1.;

    // set the permanent multipole and induced dipole values;

    real ci = atom1.q;

    real bn1 = bn.x;
    real bn2 = bn.y;
    real bn3 = bn.z;
    real bn4 = bn.w;

    real rr5 = rr1*rr1;
    rr5 = 3*rr1*rr5*rr5;
    real rr7 = 5.0 * rr5 * rr1 * rr1;

    real3 ftm2 = make_real3(0);

    real rr3 = rr1*rr1*rr1;

    // RealOpenMM scale3CD =getAndScaleInverseRs(particleI,particleJ,r,true,3,TCD);

    float damp      = POW(atom1.damp*atom2.damp, 1.0f/6.0f); // AA in MBPol

    bool do_scaling = (damp != 0.0) & ( damp > -50.0 );

    real r = SQRT(xr*xr + yr*yr + zr*zr);
    real ratio       = POW(r/damp, 4.); // rA4 in MBPol

    real dampForExpCD = -1 * thole[TCD] * ratio;
    real scale3CD = 1.0;
    if (do_scaling)
        scale3CD -= EXP(dampForExpCD);
    real scale5CD = scale3CD;
    if (do_scaling)
        scale5CD -= (4./3.) * thole[TCD] * EXP(dampForExpCD) * ratio;

    // in PME same water interactions are not excluded,
    // but the scale factors are set to 0.

    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    scale3CD *= !isSameWater;
    scale5CD *= !isSameWater;

    real sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    energy += forceFactor*0.5f*sci4*(rr3 * (1 - scale3CD) - bn1)*ci;

    real scip4 = atom2.inducedDipolePolar.x*xr + atom2.inducedDipolePolar.y*yr + atom2.inducedDipolePolar.z*zr;
//#ifndef DIRECT_POLARIZATION
//    prefactor1 = 0.5f*(bn2 );
//    ftm21 += prefactor1*((sci4*atom1.inducedDipolePolar.x + scip4*atom1.inducedDipole.x));
//    ftm22 += prefactor1*((sci4*atom1.inducedDipolePolar.y + scip4*atom1.inducedDipole.y));
//    ftm23 += prefactor1*((sci4*atom1.inducedDipolePolar.z + scip4*atom1.inducedDipole.z));
//#endif

    real ck = atom2.q;
    real sci3 = dot(atom1.inducedDipole, delta3);
    real scip3 = dot(atom1.inducedDipolePolar, delta3);
    energy += forceFactor*0.5f*sci3*(ck*(bn1-rr3 * (1 - scale3CD)));

    real gli1 = ck*sci3 - ci*sci4;
    real glip1 = ck * scip3 -ci*scip4;
    // get the induced force with screening


    real scip2 = dot(atom1.inducedDipole, atom2.inducedDipolePolar) +
                 dot(atom2.inducedDipole, atom1.inducedDipolePolar);

    real gfi1 = bn2*(gli1+glip1+scip2) - bn3*(scip3*sci4+sci3*scip4);
    //gfi1 -= (rr1*rr1)*(3*(gli1*psc3 + glip1*dsc3) + 5*(gli2*psc5 + glip2*dsc5));
    gfi1 *= 0.5f;
    ftm2 += gfi1 * delta3;

    real gfi2 = (-ck*bn1 );
    real gfi3 = ci*bn1;

    ftm2 += 0.5f * gfi2 * (atom1.inducedDipole + atom1.inducedDipolePolar);
    ftm2 += 0.5f * gfi3 * (atom2.inducedDipole + atom2.inducedDipolePolar);

    ftm2 += 0.5f*bn2*(sci3*atom2.inducedDipolePolar + scip3*atom2.inducedDipole);
    ftm2 += 0.5f*bn2*(sci4*atom1.inducedDipolePolar + scip4*atom1.inducedDipole);

    {
        real expdamp = EXP(damp);
        real temp3 = -1.5f*damp*expdamp*rr1*rr1;
        real temp5 = -damp;
        real temp7 = -0.2f - 0.6f*damp;

        real3 ddsc3 = temp3*delta3;

        real rr3 = rr1*rr1*rr1;
        temp3 = (gli1*pScale + glip1*dScale);
        ftm2 -= rr3*temp3*ddsc3;
    }

    // get the induced force without screening

    int t = TDD;
    if ((isSameWater) & ((atom1.atomType == 0) | (atom2.atomType == 0)))
        t = TDDOH;

    if ((isSameWater) & ((atom1.atomType == 1) & (atom2.atomType == 1)))
        t = TDDHH;

    real dampForExpDD = -1 * thole[t] * ratio;
    real scale5DD = 1.0;
    if (do_scaling)
        scale5DD -= EXP(dampForExpDD) + (4./3.) * thole[t] * EXP(dampForExpDD) * ratio;
    real scale7DD = scale5DD;
    if (do_scaling)
        scale7DD -= (4./15.) * thole[t] * (4. * thole[t] * ratio - 1.) * EXP(dampForExpDD) / POW(damp, 4) * POW(r, 4);

    real gfri1 = 0.5f*(rr5 * ( gli1  * (1 - scale5CD)   // charge - inddip
                        + glip1 * (1 - scale5CD)   // charge - inddip
                  + scip2 * (1 - scale5DD) ) // inddip - inddip
                      - rr7 * (sci3*scip4+scip3*sci4)
                              * (1 - scale7DD)   // inddip - inddip
               );

    ftm2 -= gfri1 * delta3;

    ftm2 -= 0.5f*rr5*(1 - scale5DD)*((sci3*atom2.inducedDipolePolar + scip3*atom2.inducedDipole) +
                      (sci4*atom1.inducedDipolePolar + scip4*atom1.inducedDipole));

    ftm2 -= 0.5f * rr3*(1 - scale3CD) *(-(atom1.inducedDipole + atom1.inducedDipolePolar) * ck +
                                         (atom2.inducedDipole + atom2.inducedDipolePolar) * ci);

    if ((atom1.moleculeIndex ==0) & (atom1.atomType == 0) & (abs(atom2.pos.x-50+0.19) < 0.001))
    {
        printf("scale3CD: %.8g\n", scale3CD);
        printf("first part: %.8g\n", forceFactor*0.5f*sci4*(rr3 * (1 - scale3CD) - bn1)*ci);
        printf("second part: %.8g\n", forceFactor*0.5f*sci3*(ck*(bn1-rr3 * (1 - scale3CD))));
    }

    force += ftm2;
}
