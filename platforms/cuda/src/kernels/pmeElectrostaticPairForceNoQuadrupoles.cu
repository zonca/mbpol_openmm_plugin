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

    real do_scaling = (damp != 0.0) & ( damp > -50.0 ); // damp or not

    real ratio       = POW(r/damp, 4.); // rA4 in MBPol

    real dampForExpCC = -1 * thole[TCC] * ratio;
    // EXP(ttm::gammln(3.0/4.0)) = 1.2254167024651776
    real scale3CC = ( 1.0 - do_scaling*EXP(dampForExpCC) ); // needed for force
    real scale1CC = scale3CC + do_scaling*(POW(thole[TCC], 1.0f/4.0f)*(r/damp)*1.2254167024651776*gammq(3.0/4.0, -dampForExpCC));

    real dampForExpCD = -1 * thole[TCD] * ratio;
    real scale3CD = ( 1.0 - do_scaling*EXP(dampForExpCD) );

    // in PME same water interactions are not excluded,
    // but the scale factors are set to 0.

    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    scale1CC *= !isSameWater;
    scale3CD *= !isSameWater;

    energy += -forceFactor*(rr1*gl0*(1 - scale1CC));

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

    real3 ftm2 = make_real3(0);

    real rr3 = rr1*rr1*rr1;

    real gfi3 = ci*bn1;

    real prefactor1;
    prefactor1 = 0.5f;//*(ci*psc3 + sc3*psc5 - gfi3);
    ftm2 -= prefactor1*atom2.inducedDipole;

    prefactor1 = 0.5f;//*(ci*dsc3 + sc3*dsc5 - gfi3);
    ftm2 -= prefactor1*atom2.inducedDipolePolar;

    // RealOpenMM scale3CD =getAndScaleInverseRs(particleI,particleJ,r,true,3,TCD);

    float damp      = POW(atom1.damp*atom2.damp, 1.0f/6.0f); // AA in MBPol

    real do_scaling = (damp != 0.0) & ( damp > -50.0 ); // damp or not

    real r = SQRT(xr*xr + yr*yr + zr*zr);
    real ratio       = POW(r/damp, 4.); // rA4 in MBPol

    real dampForExpCD = -1 * thole[TCD] * ratio;
    real scale3CD = ( 1.0 - do_scaling*EXP(dampForExpCD) );
    real scale5CD = scale3CD - do_scaling * (4./3.) * thole[TCD] * EXP(dampForExpCD) * ratio;

    // in PME same water interactions are not excluded,
    // but the scale factors are set to 0.

    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    scale3CD *= !isSameWater;

    real sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    energy += forceFactor*0.5f*sci4*(rr3 * (1 - scale3CD) - bn1)*ci;

    real scip4 = atom2.inducedDipolePolar.x*xr + atom2.inducedDipolePolar.y*yr + atom2.inducedDipolePolar.z*zr;
//#ifndef DIRECT_POLARIZATION
//    prefactor1 = 0.5f*(bn2 );
//    ftm21 += prefactor1*((sci4*atom1.inducedDipolePolar.x + scip4*atom1.inducedDipole.x));
//    ftm22 += prefactor1*((sci4*atom1.inducedDipolePolar.y + scip4*atom1.inducedDipole.y));
//    ftm23 += prefactor1*((sci4*atom1.inducedDipolePolar.z + scip4*atom1.inducedDipole.z));
//#endif

    real gli1 = -ci*sci4;
    real glip1 = -ci*scip4;

    real gfi1 = rr5 * (gli1+glip1) * (1 - scale5CD); // charge - inddip
    gfi1 -= (rr1*rr1)*(3*(gli1*psc3 + glip1*dsc3));
    gfi1 *= 0.5f;
    ftm2 += gfi1 * delta3;

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

//K

    real ck = atom2.q;
    real gfi2 = (-ck*bn1 );

    prefactor1 = 0.5f*(ck*psc3 + gfi2);
    ftm2 += prefactor1*atom1.inducedDipole;

    prefactor1 = 0.5f*(ck*dsc3 + gfi2);
    ftm2 += prefactor1*atom1.inducedDipolePolar;

    real sci3 = atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr;
    energy += forceFactor*0.5f*sci3*(ck*(bn1-rr3 * (1 - scale3CD)));
    real scip3 = atom1.inducedDipolePolar.x*xr + atom1.inducedDipolePolar.y*yr + atom1.inducedDipolePolar.z*zr;

#ifndef DIRECT_POLARIZATION
    prefactor1 = 0.5f*(bn2 );

    ftm2 += prefactor1*(sci3*atom2.inducedDipolePolar + scip3*atom2.inducedDipole);
    
    real sci34;
    sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    scip4 = atom2.inducedDipolePolar.x*xr + atom2.inducedDipolePolar.y*yr + atom2.inducedDipolePolar.z*zr;
    sci34 = (sci3*scip4+scip3*sci4);

    gfi1 = sci34*(usc5*(5*rr1*rr1) -bn3);
//#else
//    gfi1 = 0;
//#endif
    

    real scip2 = atom1.inducedDipole.x*atom2.inducedDipolePolar.x +
                                  atom1.inducedDipole.y*atom2.inducedDipolePolar.y +
                                  atom1.inducedDipole.z*atom2.inducedDipolePolar.z +
                                  atom2.inducedDipole.x*atom1.inducedDipolePolar.x +
                                  atom2.inducedDipole.y*atom1.inducedDipolePolar.y +
                                  atom2.inducedDipole.z*atom1.inducedDipolePolar.z;

           gli1 = ck*sci3;
          glip1 = ck*scip3;


    gfi1 += (bn2*(gli1+glip1));
    gfi1 -= (rr1*rr1)*(3*(gli1*psc3 + glip1*dsc3));
//#ifndef DIRECT_POLARIZATION
    gfi1 += scip2*(bn2 - (3*rr1*rr1)*usc3);
//#endif
    
    gfi1 *= 0.5f;

    ftm2 += gfi1*delta3;

    {
        real expdamp = EXP(damp);
        real temp3 = -1.5f*damp*expdamp*rr1*rr1;
        real temp5 = -damp;
        real temp7 = -0.2f - 0.6f*damp;

        real3 ddsc3 = temp3*delta3;

        real rr3 = rr1*rr1*rr1;

        temp3 = gli1*pScale + glip1*dScale;

        ftm2 -= rr3*temp3*ddsc3;

#ifndef DIRECT_POLARIZATION
        temp3 =  uScale*scip2;
        ftm2 -= rr3*(temp3*ddsc3);
#endif
    }

    force += ftm2;
}
