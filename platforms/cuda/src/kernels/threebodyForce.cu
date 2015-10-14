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