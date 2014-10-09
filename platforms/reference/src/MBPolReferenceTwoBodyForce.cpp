
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#include "MBPolReferenceForce.h"
#include "MBPolReferenceTwoBodyForce.h"
#include <algorithm>
#include <cctype>
#include "mbpol_2body_constants.h"
#include "poly-2b-v6x.h"

using std::vector;
using OpenMM::RealVec;

MBPolReferenceTwoBodyForce::MBPolReferenceTwoBodyForce( ) : _nonbondedMethod(NoCutoff), _cutoff(1.0e+10) {

    _periodicBoxDimensions = RealVec( 0.0, 0.0, 0.0 );
}

MBPolReferenceTwoBodyForce::NonbondedMethod MBPolReferenceTwoBodyForce::getNonbondedMethod( void ) const {
    return _nonbondedMethod;
}

void MBPolReferenceTwoBodyForce::setNonbondedMethod( MBPolReferenceTwoBodyForce::NonbondedMethod nonbondedMethod ){
    _nonbondedMethod = nonbondedMethod;
}

void MBPolReferenceTwoBodyForce::setCutoff( double cutoff ){
    _cutoff  = cutoff;
}

double MBPolReferenceTwoBodyForce::getCutoff( void ) const {
    return _cutoff;
}

void MBPolReferenceTwoBodyForce::setPeriodicBox( const RealVec& box ){
    _periodicBoxDimensions = box;
}

RealVec MBPolReferenceTwoBodyForce::getPeriodicBox( void ) const {
    return _periodicBoxDimensions;
}

void imageParticles(const RealVec& box, const RealVec & referenceParticle, RealVec& particleToImage)
{
    // Periodic boundary conditions imaging of particleToImage with respect to referenceParticle

    double distance, factor;
    for (unsigned int i=0; i < 3; i++) {
        distance = referenceParticle[i] - particleToImage[i];
        factor = std::floor(distance/box[i] + 0.5);
        particleToImage[i] += box[i]*factor;
    }
}

void imageMolecules(const RealVec& box, std::vector<RealVec> allPositions)
{

    // Take first oxygen as central atom

    // image its two hydrogens with respect of the first oxygen

    imageParticles(box, allPositions[Oa], allPositions[Ha1]);
    imageParticles(box, allPositions[Oa], allPositions[Ha2]);

    // Now image the oxygen of the second molecule

    imageParticles(box, allPositions[Oa], allPositions[Ob]);

    // Image the hydrogen of the second molecule with respect to the oxygen of the second molecule
    imageParticles(box, allPositions[Ob], allPositions[Hb1]);
    imageParticles(box, allPositions[Ob], allPositions[Hb2]);

    if (allPositions.size() >= Hc1) // Three molecules
    {
        // Now image the oxygen of the third molecule
        imageParticles(box, allPositions[Oa], allPositions[Oc]);

        // Image the hydrogen of the third molecule with respect to the oxygen of the third molecule
        imageParticles(box, allPositions[Oc], allPositions[Hc1]);
        imageParticles(box, allPositions[Oc], allPositions[Hc2]);
    }

}
RealOpenMM MBPolReferenceTwoBodyForce::calculatePairIxn( int siteI, int siteJ,
                                                      const std::vector<RealVec>& particlePositions,
                                                      const std::vector<std::vector<int> >& allParticleIndices,
                                                      vector<RealVec>& forces ) const {

        // siteI and siteJ are indices in a oxygen-only array, in order to get the position of an oxygen, we need:
        // allParticleIndices[siteI][0]
        // first hydrogen: allParticleIndices[siteI][1]
        // second hydrogen: allParticleIndices[siteI][2]
        // same for the second water molecule
        const double nm_to_A = 10.;

        // offsets

        const int Oa  = 0;
        const int Ha1 = 3;
        const int Ha2 = 6;

        const int Ob  = 9;
        const int Hb1 = 12;
        const int Hb2 = 15;

        const int Xa1 = 18;
        const int Xa2 = 21;

        const int Xb1 = 24;
        const int Xb2 = 27;

        double xcrd[30]; // coordinates including extra-points

        if( _nonbondedMethod == CutoffPeriodic ){
            imageMolecules(_periodicBoxDimensions, siteI, siteJ, particlePositions, allParticleIndices, xcrd);
        } else {
            for (unsigned int i=0; i < 3; i++) {
                // first water molecule
                xcrd[Oa + i] =  particlePositions[allParticleIndices[siteI][0]][i];
                xcrd[Ha1 + i] = particlePositions[allParticleIndices[siteI][1]][i];
                xcrd[Ha2 + i] = particlePositions[allParticleIndices[siteI][2]][i];
                // second water molecule
                xcrd[Ob + i] =  particlePositions[allParticleIndices[siteJ][0]][i];
                xcrd[Hb1 + i] = particlePositions[allParticleIndices[siteJ][1]][i];
                xcrd[Hb2 + i] = particlePositions[allParticleIndices[siteJ][2]][i];
            }
        }


        for (unsigned int i=0; i < 6*3; i++) {
            xcrd[i] *= nm_to_A;
        }

        const double dOO[3] = {xcrd[Oa] - xcrd[Ob],
                               xcrd[Oa + 1] - xcrd[Ob + 1],
                               xcrd[Oa + 2] - xcrd[Ob + 2],};

        const double rOOsq = dOO[0]*dOO[0] + dOO[1]*dOO[1] + dOO[2]*dOO[2];
        const double rOO = std::sqrt(rOOsq);

        if (rOO > r2f)
            return 0.0;

        if (rOO < 2.)
            return 0.0;


        // the extra-points

        monomer ma, mb;

        ma.setup(xcrd + Oa,
                 in_plane_gamma, out_of_plane_gamma,
                 xcrd + Xa1, xcrd + Xa2);

        mb.setup(xcrd + Ob,
                 in_plane_gamma, out_of_plane_gamma,
                 xcrd + Xb1, xcrd + Xb2);

        // variables

        const double d0_intra = 1.0;
        const double d0_inter = 4.0;

        double v[31]; // stored separately (gets passed to poly::eval)

        variable ctxt[31];

        v[0] = ctxt[0].v_exp(d0_intra, k_HH_intra, xcrd, Ha1, Ha2);
        v[1] = ctxt[1].v_exp(d0_intra, k_HH_intra, xcrd, Hb1, Hb2);

        v[2] = ctxt[2].v_exp(d0_intra, k_OH_intra, xcrd, Oa, Ha1);
        v[3] = ctxt[3].v_exp(d0_intra, k_OH_intra, xcrd, Oa, Ha2);
        v[4] = ctxt[4].v_exp(d0_intra, k_OH_intra, xcrd, Ob, Hb1);
        v[5] = ctxt[5].v_exp(d0_intra, k_OH_intra, xcrd, Ob, Hb2);

        v[6] = ctxt[6].v_coul(d0_inter, k_HH_coul, xcrd, Ha1, Hb1);
        v[7] = ctxt[7].v_coul(d0_inter, k_HH_coul, xcrd, Ha1, Hb2);
        v[8] = ctxt[8].v_coul(d0_inter, k_HH_coul, xcrd, Ha2, Hb1);
        v[9] = ctxt[9].v_coul(d0_inter, k_HH_coul, xcrd, Ha2, Hb2);

        v[10] = ctxt[10].v_coul(d0_inter, k_OH_coul, xcrd, Oa, Hb1);
        v[11] = ctxt[11].v_coul(d0_inter, k_OH_coul, xcrd, Oa, Hb2);
        v[12] = ctxt[12].v_coul(d0_inter, k_OH_coul, xcrd, Ob, Ha1);
        v[13] = ctxt[13].v_coul(d0_inter, k_OH_coul, xcrd, Ob, Ha2);

        v[14] = ctxt[14].v_coul(d0_inter, k_OO_coul, xcrd, Oa, Ob);

        v[15] = ctxt[15].v_exp(d0_inter, k_XH_main, xcrd, Xa1, Hb1);
        v[16] = ctxt[16].v_exp(d0_inter, k_XH_main, xcrd, Xa1, Hb2);
        v[17] = ctxt[17].v_exp(d0_inter, k_XH_main, xcrd, Xa2, Hb1);
        v[18] = ctxt[18].v_exp(d0_inter, k_XH_main, xcrd, Xa2, Hb2);
        v[19] = ctxt[19].v_exp(d0_inter, k_XH_main, xcrd, Xb1, Ha1);
        v[20] = ctxt[20].v_exp(d0_inter, k_XH_main, xcrd, Xb1, Ha2);
        v[21] = ctxt[21].v_exp(d0_inter, k_XH_main, xcrd, Xb2, Ha1);
        v[22] = ctxt[22].v_exp(d0_inter, k_XH_main, xcrd, Xb2, Ha2);

        v[23] = ctxt[23].v_exp(d0_inter, k_XO_main, xcrd, Oa, Xb1);
        v[24] = ctxt[24].v_exp(d0_inter, k_XO_main, xcrd, Oa, Xb2);
        v[25] = ctxt[25].v_exp(d0_inter, k_XO_main, xcrd, Ob, Xa1);
        v[26] = ctxt[26].v_exp(d0_inter, k_XO_main, xcrd, Ob, Xa2);

        v[27] = ctxt[27].v_exp(d0_inter, k_XX_main, xcrd, Xa1, Xb1);
        v[28] = ctxt[28].v_exp(d0_inter, k_XX_main, xcrd, Xa1, Xb2);
        v[29] = ctxt[29].v_exp(d0_inter, k_XX_main, xcrd, Xa2, Xb1);
        v[30] = ctxt[30].v_exp(d0_inter, k_XX_main, xcrd, Xa2, Xb2);

        double g[31];
        const double E_poly = poly_2b_v6x_eval(thefit, v, g);

        double xgrd[30];
        std::fill(xgrd, xgrd + 30, 0.0);

        ctxt[0].grads(g[0], xgrd, Ha1, Ha2);
        ctxt[1].grads(g[1], xgrd, Hb1, Hb2);

        ctxt[2].grads(g[2], xgrd, Oa, Ha1);
        ctxt[3].grads(g[3], xgrd, Oa, Ha2);
        ctxt[4].grads(g[4], xgrd, Ob, Hb1);
        ctxt[5].grads(g[5], xgrd, Ob, Hb2);

        ctxt[6].grads(g[6], xgrd, Ha1, Hb1);
        ctxt[7].grads(g[7], xgrd, Ha1, Hb2);
        ctxt[8].grads(g[8], xgrd, Ha2, Hb1);
        ctxt[9].grads(g[9], xgrd, Ha2, Hb2);

        ctxt[10].grads(g[10], xgrd, Oa, Hb1);
        ctxt[11].grads(g[11], xgrd, Oa, Hb2);
        ctxt[12].grads(g[12], xgrd, Ob, Ha1);
        ctxt[13].grads(g[13], xgrd, Ob, Ha2);

        ctxt[14].grads(g[14], xgrd, Oa, Ob);

        ctxt[15].grads(g[15], xgrd, Xa1, Hb1);
        ctxt[16].grads(g[16], xgrd, Xa1, Hb2);
        ctxt[17].grads(g[17], xgrd, Xa2, Hb1);
        ctxt[18].grads(g[18], xgrd, Xa2, Hb2);
        ctxt[19].grads(g[19], xgrd, Xb1, Ha1);
        ctxt[20].grads(g[20], xgrd, Xb1, Ha2);
        ctxt[21].grads(g[21], xgrd, Xb2, Ha1);
        ctxt[22].grads(g[22], xgrd, Xb2, Ha2);

        ctxt[23].grads(g[23], xgrd, Oa, Xb1);
        ctxt[24].grads(g[24], xgrd, Oa, Xb2);
        ctxt[25].grads(g[25], xgrd, Ob, Xa1);
        ctxt[26].grads(g[26], xgrd, Ob, Xa2);

        ctxt[27].grads(g[27], xgrd, Xa1, Xb1);
        ctxt[28].grads(g[28], xgrd, Xa1, Xb2);
        ctxt[29].grads(g[29], xgrd, Xa2, Xb1);
        ctxt[30].grads(g[30], xgrd, Xa2, Xb2);

        // distribute gradients w.r.t. the X-points

        ma.grads(xgrd + Xa1, xgrd + Xa2,
                 in_plane_gamma, out_of_plane_gamma,
                 xgrd + Oa);

        mb.grads(xgrd + Xb1, xgrd + Xb2,
                 in_plane_gamma, out_of_plane_gamma,
                 xgrd + Ob);

        // the switch

        double gsw;
        const double sw = f_switch(rOO, gsw);

        double cal2joule = 4.184;

        for (int i = 0; i < 3; ++i) {
            // first water molecule
            forces[allParticleIndices[siteI][0]][i] += sw*xgrd[Oa + i]  * cal2joule * -10.;
            forces[allParticleIndices[siteI][1]][i] += sw*xgrd[Ha1 + i] * cal2joule * -10.;
            forces[allParticleIndices[siteI][2]][i] += sw*xgrd[Ha2 + i] * cal2joule * -10.;
            // second water molecule
            forces[allParticleIndices[siteJ][0]][i] += sw*xgrd[Ob + i]  * cal2joule * -10.;
            forces[allParticleIndices[siteJ][1]][i] += sw*xgrd[Hb1 + i] * cal2joule * -10.;
            forces[allParticleIndices[siteJ][2]][i] += sw*xgrd[Hb2 + i] * cal2joule * -10.;
        }

        // gradient of the switch
        gsw *= E_poly/rOO;
        for (int i = 0; i < 3; ++i) {
            const double d = gsw*dOO[i];
            forces[allParticleIndices[siteI][0]][i] += d * cal2joule * -10.;
            forces[allParticleIndices[siteJ][0]][i] -= d * cal2joule * -10.;
        }

    RealOpenMM energy=sw*E_poly * cal2joule;

    return energy;

}

RealOpenMM MBPolReferenceTwoBodyForce::calculateForceAndEnergy( int numParticles,
                                                             const vector<RealVec>& particlePositions,
                                                             const std::vector<std::vector<int> >& allParticleIndices,
                                                             const NeighborList& neighborList,
                                                             vector<RealVec>& forces ) const {

    // loop over neighbor list
    //    (1) calculate pair TwoBody ixn
    //    (2) accumulate forces: if particle is a site where interaction position != particle position,
    //        then call addReducedForce() to apportion force to particle and its covalent partner
    //        based on reduction factor

    RealOpenMM energy = 0.;
    for( unsigned int ii = 0; ii < neighborList.size(); ii++ ){

        OpenMM::AtomPair pair       = neighborList[ii];
        int siteI                   = pair.first;
        int siteJ                   = pair.second;

        energy                     += calculatePairIxn( siteI, siteJ,
                particlePositions, allParticleIndices, forces );

    }

    return energy;
}
