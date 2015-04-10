const real nm_to_A = 10.;
const real kcal_permol_Aminus6_to_kJ_permol_nmminus6 = 4.184e6;
const real cal2joule = 4.184;

#ifdef USE_CUTOFF
if (!isExcluded && atom1 != atom2 && r2 < CUTOFF_SQUARED) {
#else
if (!isExcluded && atom1 != atom2) {
#endif

    real c6 = 3.;
    real d6 = 10.;

    real d6r = r * d6; // r in nm and d6 is 1/nm so units compensate
    // tang_toennies function
    int nn = 6;

    real sum = 1.0 + d6r/nn;

    while (--nn != 0)
        sum = 1.0 + sum*d6r/nn;

    real tt6 = 1.0 - sum*EXP(-d6r);

    if (abs(tt6) < 1.0e-8) {

        real term(1);
        for (nn = 6; nn != 0; --nn)
            term *= d6r/nn;

        sum = 0.0;
        for (nn = 6 + 1; nn < 1000; ++nn) {
            term *= d6r/nn;
            sum += term;

            if (abs(term/sum) < 1.0e-8)
                break;
        }

        tt6 = sum*std::exp(-d6r);
    }
    // end of tang_toennies function

    real e6 = c6*tt6*POW(1/(r * nm_to_A), 6)/kcal_permol_Aminus6_to_kJ_permol_nmminus6;

    real if6 = 1.0/720; // factorial(6)

    dEdR += (6*e6/(r * r * nm_to_A) - 
      c6/kcal_permol_Aminus6_to_kJ_permol_nmminus6*POW(d6/nm_to_A, 7)*if6*EXP(-d6r)/r) * cal2joule;
    // factored out 1 nm_to_A factor

    tempEnergy += -e6 * cal2joule;
}
