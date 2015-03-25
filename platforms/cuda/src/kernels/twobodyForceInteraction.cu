
if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
    real sigma = 1.3;
    real epsilon = 2.;
    real x = sigma/r;
    real eps = SQRT(2);
    dEdR += 2 * epsilon*eps*(12*POW(x, 12.0)-6*POW(x, 6.0)) * invR * invR;
    tempEnergy += 4.0*eps*(POW(x, 12.0)-POW(x, 6.0));
}
