/* -------------------------------------------------------------------------- *
 *                               OpenMMMBPol                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2013 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "MBPolCudaKernels.h"
#include "CudaMBPolKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/MBPolElectrostaticsForceImpl.h"
#include "openmm/internal/MBPolTwoBodyForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "CudaBondedUtilities.h"
#include "CudaForceInfo.h"
#include "CudaKernelSources.h"
#include "CudaNonbondedUtilities.h"

#include <algorithm>
#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace  OpenMM;
using namespace MBPolPlugin;
using namespace std;

#define CHECK_RESULT(result) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<errorMessage<<": "<<cu.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

/* -------------------------------------------------------------------------- *
 *                            MBPolBondForce                                 *
 * -------------------------------------------------------------------------- */

class CudaCalcMBPolOneBodyForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const MBPolBondForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumOneBodys();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        force.getOneBodyParameters(index, particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        // They can be identical only if they are the same
        // FIXME we could optionally check all distances.
        return (group1 == group2)
    }
private:
    const MBPolOneBodyForce& force;
};

CudaCalcMBPolBondForceKernel::CudaCalcMBPolBondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : 
                CalcMBPolBondForceKernel(name, platform), cu(cu), system(system), params(NULL) {
}

CudaCalcMBPolBondForceKernel::~CudaCalcMBPolBondForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcMBPolBondForceKernel::initialize(const System& system, const MBPolBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = CudaArray::create<float2>(cu, numBonds, "bondParams");
    vector<float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = make_float2((float) length, (float) k);
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = CudaMBPolKernelSources::mbpolBondForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float2");
    replacements["CUBIC_K"] = cu.doubleToString(force.getMBPolGlobalBondCubic());
    replacements["QUARTIC_K"] = cu.doubleToString(force.getMBPolGlobalBondQuartic());
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::bondForce, replacements), force.getForceGroup());
    cu.addForce(new ForceInfo(force));
}

double CudaCalcMBPolBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcMBPolBondForceKernel::copyParametersToContext(ContextImpl& context, const MBPolBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int atom1, atom2;
        double length, k;
        force.getBondParameters(startIndex+i, atom1, atom2, length, k);
        paramVector[i] = make_float2((float) length, (float) k);
    }
    params->upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

/* -------------------------------------------------------------------------- *
 *                             MBPolMultipole                                *
 * -------------------------------------------------------------------------- */

class CudaCalcMBPolMultipoleForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const MBPolMultipoleForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, thole1, thole2, damping1, damping2, polarity1, polarity2;
        int axis1, axis2, multipole11, multipole12, multipole21, multipole22, multipole31, multipole32;
        vector<double> dipole1, dipole2, quadrupole1, quadrupole2;
        force.getMultipoleParameters(particle1, charge1, dipole1, quadrupole1, axis1, multipole11, multipole21, multipole31, thole1, damping1, polarity1);
        force.getMultipoleParameters(particle2, charge2, dipole2, quadrupole2, axis2, multipole12, multipole22, multipole32, thole2, damping2, polarity2);
        if (charge1 != charge2 || thole1 != thole2 || damping1 != damping2 || polarity1 != polarity2 || axis1 != axis2){
            return false;
        }
        for (int i = 0; i < (int) dipole1.size(); ++i){
            if (dipole1[i] != dipole2[i]){
                return false;
            }
        }
        for (int i = 0; i < (int) quadrupole1.size(); ++i){
            if (quadrupole1[i] != quadrupole2[i]){
                return false;
            }
        }
        return true;
    }
    int getNumParticleGroups() {
        return 7*force.getNumMultipoles();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle = index/7;
        int type = index-7*particle;
        force.getCovalentMap(particle, MBPolMultipoleForce::CovalentType(type), particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        return ((group1%7) == (group2%7));
    }
private:
    const MBPolMultipoleForce& force;
};

CudaCalcMBPolMultipoleForceKernel::CudaCalcMBPolMultipoleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : 
        CalcMBPolMultipoleForceKernel(name, platform), cu(cu), system(system), hasInitializedScaleFactors(false), hasInitializedFFT(false), multipolesAreValid(false),
        multipoleParticles(NULL), molecularDipoles(NULL), molecularQuadrupoles(NULL), labFrameDipoles(NULL), labFrameQuadrupoles(NULL),
        field(NULL), fieldPolar(NULL), inducedField(NULL), inducedFieldPolar(NULL), torque(NULL), dampingAndThole(NULL),
        inducedDipole(NULL), inducedDipolePolar(NULL), inducedDipoleErrors(NULL), polarizability(NULL), covalentFlags(NULL), polarizationGroupFlags(NULL),
        pmeGrid(NULL), pmeBsplineModuliX(NULL), pmeBsplineModuliY(NULL), pmeBsplineModuliZ(NULL), pmeIgrid(NULL), pmePhi(NULL),
        pmePhid(NULL), pmePhip(NULL), pmePhidp(NULL), pmeAtomGridIndex(NULL), lastPositions(NULL), sort(NULL), gkKernel(NULL) {
}

CudaCalcMBPolMultipoleForceKernel::~CudaCalcMBPolMultipoleForceKernel() {
    cu.setAsCurrent();
    if (multipoleParticles != NULL)
        delete multipoleParticles;
    if (molecularDipoles != NULL)
        delete molecularDipoles;
    if (molecularQuadrupoles != NULL)
        delete molecularQuadrupoles;
    if (labFrameDipoles != NULL)
        delete labFrameDipoles;
    if (labFrameQuadrupoles != NULL)
        delete labFrameQuadrupoles;
    if (field != NULL)
        delete field;
    if (fieldPolar != NULL)
        delete fieldPolar;
    if (inducedField != NULL)
        delete inducedField;
    if (inducedFieldPolar != NULL)
        delete inducedFieldPolar;
    if (torque != NULL)
        delete torque;
    if (dampingAndThole != NULL)
        delete dampingAndThole;
    if (inducedDipole != NULL)
        delete inducedDipole;
    if (inducedDipolePolar != NULL)
        delete inducedDipolePolar;
    if (inducedDipoleErrors != NULL)
        delete inducedDipoleErrors;
    if (polarizability != NULL)
        delete polarizability;
    if (covalentFlags != NULL)
        delete covalentFlags;
    if (polarizationGroupFlags != NULL)
        delete polarizationGroupFlags;
    if (pmeGrid != NULL)
        delete pmeGrid;
    if (pmeBsplineModuliX != NULL)
        delete pmeBsplineModuliX;
    if (pmeBsplineModuliY != NULL)
        delete pmeBsplineModuliY;
    if (pmeBsplineModuliZ != NULL)
        delete pmeBsplineModuliZ;
    if (pmeIgrid != NULL)
        delete pmeIgrid;
    if (pmePhi != NULL)
        delete pmePhi;
    if (pmePhid != NULL)
        delete pmePhid;
    if (pmePhip != NULL)
        delete pmePhip;
    if (pmePhidp != NULL)
        delete pmePhidp;
    if (pmeAtomGridIndex != NULL)
        delete pmeAtomGridIndex;
    if (lastPositions != NULL)
        delete lastPositions;
    if (sort != NULL)
        delete sort;
    if (hasInitializedFFT)
        cufftDestroy(fft);
}

/**
 * Select a size for an FFT that is a multiple of 2, 3, 5, and 7.
 */
static int findFFTDimension(int minimum) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        int unfactored = minimum;
        for (int factor = 2; factor < 8; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1)
            return minimum;
        minimum++;
    }
}

void CudaCalcMBPolMultipoleForceKernel::initialize(const System& system, const MBPolMultipoleForce& force) {
    cu.setAsCurrent();

    // Initialize multipole parameters.

    numMultipoles = force.getNumMultipoles();
    CudaArray& posq = cu.getPosq();
    float4* posqf = (float4*) cu.getPinnedBuffer();
    double4* posqd = (double4*) cu.getPinnedBuffer();
    vector<float2> dampingAndTholeVec;
    vector<float> polarizabilityVec;
    vector<float> molecularDipolesVec;
    vector<float> molecularQuadrupolesVec;
    vector<int4> multipoleParticlesVec;
    for (int i = 0; i < numMultipoles; i++) {
        double charge, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getMultipoleParameters(i, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        if (cu.getUseDoublePrecision())
            posqd[i] = make_double4(0, 0, 0, charge);
        else
            posqf[i] = make_float4(0, 0, 0, (float) charge);
        dampingAndTholeVec.push_back(make_float2((float) damping, (float) thole));
        polarizabilityVec.push_back((float) polarity);
        multipoleParticlesVec.push_back(make_int4(atomX, atomY, atomZ, axisType));
        for (int j = 0; j < 3; j++)
            molecularDipolesVec.push_back((float) dipole[j]);
        molecularQuadrupolesVec.push_back((float) quadrupole[0]);
        molecularQuadrupolesVec.push_back((float) quadrupole[1]);
        molecularQuadrupolesVec.push_back((float) quadrupole[2]);
        molecularQuadrupolesVec.push_back((float) quadrupole[4]);
        molecularQuadrupolesVec.push_back((float) quadrupole[5]);
    }
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    for (int i = numMultipoles; i < paddedNumAtoms; i++) {
        dampingAndTholeVec.push_back(make_float2(0, 0));
        polarizabilityVec.push_back(0);
        multipoleParticlesVec.push_back(make_int4(0, 0, 0, 0));
        for (int j = 0; j < 3; j++)
            molecularDipolesVec.push_back(0);
        for (int j = 0; j < 5; j++)
            molecularQuadrupolesVec.push_back(0);
    }
    dampingAndThole = CudaArray::create<float2>(cu, paddedNumAtoms, "dampingAndThole");
    polarizability = CudaArray::create<float>(cu, paddedNumAtoms, "polarizability");
    multipoleParticles = CudaArray::create<int4>(cu, paddedNumAtoms, "multipoleParticles");
    molecularDipoles = CudaArray::create<float>(cu, 3*paddedNumAtoms, "molecularDipoles");
    molecularQuadrupoles = CudaArray::create<float>(cu, 5*paddedNumAtoms, "molecularQuadrupoles");
    lastPositions = new CudaArray(cu, cu.getPosq().getSize(), cu.getPosq().getElementSize(), "lastPositions");
    dampingAndThole->upload(dampingAndTholeVec);
    polarizability->upload(polarizabilityVec);
    multipoleParticles->upload(multipoleParticlesVec);
    molecularDipoles->upload(molecularDipolesVec);
    molecularQuadrupoles->upload(molecularQuadrupolesVec);
    posq.upload(cu.getPinnedBuffer());
    
    // Create workspace arrays.
    
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    labFrameDipoles = new CudaArray(cu, 3*paddedNumAtoms, elementSize, "labFrameDipoles");
    labFrameQuadrupoles = new CudaArray(cu, 9*paddedNumAtoms, elementSize, "labFrameQuadrupoles");
    field = new CudaArray(cu, 3*paddedNumAtoms, sizeof(long long), "field");
    fieldPolar = new CudaArray(cu, 3*paddedNumAtoms, sizeof(long long), "fieldPolar");
    torque = new CudaArray(cu, 3*paddedNumAtoms, sizeof(long long), "torque");
    inducedDipole = new CudaArray(cu, 3*paddedNumAtoms, elementSize, "inducedDipole");
    inducedDipolePolar = new CudaArray(cu, 3*paddedNumAtoms, elementSize, "inducedDipolePolar");
    inducedDipoleErrors = new CudaArray(cu, cu.getNumThreadBlocks(), sizeof(float2), "inducedDipoleErrors");
    cu.addAutoclearBuffer(*field);
    cu.addAutoclearBuffer(*fieldPolar);
    cu.addAutoclearBuffer(*torque);
    
    // Record which atoms should be flagged as exclusions based on covalent groups, and determine
    // the values for the covalent group flags.
    
    vector<vector<int> > exclusions(numMultipoles);
    for (int i = 0; i < numMultipoles; i++) {
        vector<int> atoms;
        set<int> allAtoms;
        allAtoms.insert(i);
        force.getCovalentMap(i, MBPolMultipoleForce::Covalent12, atoms);
        allAtoms.insert(atoms.begin(), atoms.end());
        force.getCovalentMap(i, MBPolMultipoleForce::Covalent13, atoms);
        allAtoms.insert(atoms.begin(), atoms.end());
        for (set<int>::const_iterator iter = allAtoms.begin(); iter != allAtoms.end(); ++iter)
            covalentFlagValues.push_back(make_int3(i, *iter, 0));
        force.getCovalentMap(i, MBPolMultipoleForce::Covalent14, atoms);
        allAtoms.insert(atoms.begin(), atoms.end());
        for (int j = 0; j < (int) atoms.size(); j++)
            covalentFlagValues.push_back(make_int3(i, atoms[j], 1));
        force.getCovalentMap(i, MBPolMultipoleForce::Covalent15, atoms);
        for (int j = 0; j < (int) atoms.size(); j++)
            covalentFlagValues.push_back(make_int3(i, atoms[j], 2));
        allAtoms.insert(atoms.begin(), atoms.end());
        force.getCovalentMap(i, MBPolMultipoleForce::PolarizationCovalent11, atoms);
        allAtoms.insert(atoms.begin(), atoms.end());
        exclusions[i].insert(exclusions[i].end(), allAtoms.begin(), allAtoms.end());

        // Workaround for bug in TINKER: if an atom is listed in both the PolarizationCovalent11
        // and PolarizationCovalent12 maps, the latter takes precedence.

        vector<int> atoms12;
        force.getCovalentMap(i, MBPolMultipoleForce::PolarizationCovalent12, atoms12);
        for (int j = 0; j < (int) atoms.size(); j++)
            if (find(atoms12.begin(), atoms12.end(), atoms[j]) == atoms12.end())
                polarizationFlagValues.push_back(make_int2(i, atoms[j]));
    }
    set<pair<int, int> > tilesWithExclusions;
    for (int atom1 = 0; atom1 < (int) exclusions.size(); ++atom1) {
        int x = atom1/CudaContext::TileSize;
        for (int j = 0; j < (int) exclusions[atom1].size(); ++j) {
            int atom2 = exclusions[atom1][j];
            int y = atom2/CudaContext::TileSize;
            tilesWithExclusions.insert(make_pair(max(x, y), min(x, y)));
        }
    }
    
    // Record other options.
    
    if (force.getPolarizationType() == MBPolMultipoleForce::Mutual) {
        maxInducedIterations = force.getMutualInducedMaxIterations();
        inducedEpsilon = force.getMutualInducedTargetEpsilon();
        inducedField = new CudaArray(cu, 3*paddedNumAtoms, sizeof(long long), "inducedField");
        inducedFieldPolar = new CudaArray(cu, 3*paddedNumAtoms, sizeof(long long), "inducedFieldPolar");
    }
    else
        maxInducedIterations = 0;
    bool usePME = (force.getNonbondedMethod() == MBPolMultipoleForce::PME);
    
    // See whether there's an MBPolGeneralizedKirkwoodForce in the System.

    const MBPolGeneralizedKirkwoodForce* gk = NULL;
    for (int i = 0; i < system.getNumForces() && gk == NULL; i++)
        gk = dynamic_cast<const MBPolGeneralizedKirkwoodForce*>(&system.getForce(i));
    double innerDielectric = (gk == NULL ? 1.0 : gk->getSoluteDielectric());
    
    // Create the kernels.

    bool useShuffle = (cu.getComputeCapability() >= 3.0 && !cu.getUseDoublePrecision());
    double fixedThreadMemory = 19*elementSize+2*sizeof(float)+3*sizeof(int)/(double) cu.TileSize;
    double inducedThreadMemory = 15*elementSize+2*sizeof(float);
    double electrostaticsThreadMemory = 0;
    if (!useShuffle)
        fixedThreadMemory += 3*elementSize;
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(numMultipoles);
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
    defines["ENERGY_SCALE_FACTOR"] = cu.doubleToString(138.9354558456/innerDielectric);
    if (force.getPolarizationType() == MBPolMultipoleForce::Direct)
        defines["DIRECT_POLARIZATION"] = "";
    if (useShuffle)
        defines["USE_SHUFFLE"] = "";
    defines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
    int numExclusionTiles = tilesWithExclusions.size();
    defines["NUM_TILES_WITH_EXCLUSIONS"] = cu.intToString(numExclusionTiles);
    int numContexts = cu.getPlatformData().contexts.size();
    int startExclusionIndex = cu.getContextIndex()*numExclusionTiles/numContexts;
    int endExclusionIndex = (cu.getContextIndex()+1)*numExclusionTiles/numContexts;
    defines["FIRST_EXCLUSION_TILE"] = cu.intToString(startExclusionIndex);
    defines["LAST_EXCLUSION_TILE"] = cu.intToString(endExclusionIndex);
    double alpha = force.getAEwald();
    int gridSizeX, gridSizeY, gridSizeZ;
    if (usePME) {
        vector<int> pmeGridDimension;
        force.getPmeGridDimensions(pmeGridDimension);
        if (pmeGridDimension[0] == 0 || alpha == 0.0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            NonbondedForceImpl::calcPMEParameters(system, nb, alpha, gridSizeX, gridSizeY, gridSizeZ);
            gridSizeX = findFFTDimension(gridSizeX);
            gridSizeY = findFFTDimension(gridSizeY);
            gridSizeZ = findFFTDimension(gridSizeZ);
        } else {
            gridSizeX = pmeGridDimension[0];
            gridSizeY = pmeGridDimension[1];
            gridSizeZ = pmeGridDimension[2];
        }
        defines["EWALD_ALPHA"] = cu.doubleToString(alpha);
        defines["SQRT_PI"] = cu.doubleToString(sqrt(M_PI));
        defines["USE_EWALD"] = "";
        defines["USE_CUTOFF"] = "";
        defines["USE_PERIODIC"] = "";
        defines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    if (gk != NULL) {
        defines["USE_GK"] = "";
        defines["GK_C"] = cu.doubleToString(2.455);
        double solventDielectric = gk->getSolventDielectric();
        defines["GK_FC"] = cu.doubleToString(1*(1-solventDielectric)/(0+1*solventDielectric));
        defines["GK_FD"] = cu.doubleToString(2*(1-solventDielectric)/(1+2*solventDielectric));
        defines["GK_FQ"] = cu.doubleToString(3*(1-solventDielectric)/(2+3*solventDielectric));
        fixedThreadMemory += 4*elementSize;
        inducedThreadMemory += 13*elementSize;
    }
    int maxThreads = cu.getNonbondedUtilities().getForceThreadBlockSize();
    fixedFieldThreads = min(maxThreads, cu.computeThreadBlockSize(fixedThreadMemory));
    inducedFieldThreads = min(maxThreads, cu.computeThreadBlockSize(inducedThreadMemory));
    CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaMBPolKernelSources::multipoles, defines);
    computeMomentsKernel = cu.getKernel(module, "computeLabFrameMoments");
    recordInducedDipolesKernel = cu.getKernel(module, "recordInducedDipoles");
    mapTorqueKernel = cu.getKernel(module, "mapTorqueToForce");
    computePotentialKernel = cu.getKernel(module, "computePotentialAtPoints");
    defines["THREAD_BLOCK_SIZE"] = cu.intToString(fixedFieldThreads);
    module = cu.createModule(CudaKernelSources::vectorOps+CudaMBPolKernelSources::multipoleFixedField, defines);
    computeFixedFieldKernel = cu.getKernel(module, "computeFixedField");
    if (maxInducedIterations > 0) {
        defines["THREAD_BLOCK_SIZE"] = cu.intToString(inducedFieldThreads);
        module = cu.createModule(CudaKernelSources::vectorOps+CudaMBPolKernelSources::multipoleInducedField, defines);
        computeInducedFieldKernel = cu.getKernel(module, "computeInducedField");
        updateInducedFieldKernel = cu.getKernel(module, "updateInducedFieldBySOR");
    }
    stringstream electrostaticsSource;
    if (usePME) {
        electrostaticsSource << CudaKernelSources::vectorOps;
        electrostaticsSource << CudaMBPolKernelSources::pmeMultipoleElectrostatics;
        electrostaticsSource << CudaMBPolKernelSources::pmeElectrostaticPairForce;
        electrostaticsSource << "#define APPLY_SCALE\n";
        electrostaticsSource << CudaMBPolKernelSources::pmeElectrostaticPairForce;
        electrostaticsThreadMemory = 24*elementSize+3*sizeof(float)+3*sizeof(int)/(double) cu.TileSize;
        if (!useShuffle)
            electrostaticsThreadMemory += 3*elementSize;
    }
    else {
        electrostaticsSource << CudaKernelSources::vectorOps;
        electrostaticsSource << CudaMBPolKernelSources::multipoleElectrostatics;
        electrostaticsSource << "#define F1\n";
        electrostaticsSource << CudaMBPolKernelSources::electrostaticPairForce;
        electrostaticsSource << "#undef F1\n";
        electrostaticsSource << "#define T1\n";
        electrostaticsSource << CudaMBPolKernelSources::electrostaticPairForce;
        electrostaticsSource << "#undef T1\n";
        electrostaticsSource << "#define T3\n";
        electrostaticsSource << CudaMBPolKernelSources::electrostaticPairForce;
        electrostaticsThreadMemory = 21*elementSize+2*sizeof(float)+3*sizeof(int)/(double) cu.TileSize;
        if (!useShuffle)
            electrostaticsThreadMemory += 3*elementSize;
        if (gk != NULL)
            electrostaticsThreadMemory += 4*elementSize;
    }
    electrostaticsThreads = min(maxThreads, cu.computeThreadBlockSize(electrostaticsThreadMemory));
    defines["THREAD_BLOCK_SIZE"] = cu.intToString(electrostaticsThreads);
    module = cu.createModule(electrostaticsSource.str(), defines);
    electrostaticsKernel = cu.getKernel(module, "computeElectrostatics");

    // Set up PME.
    
    if (usePME) {
        // Create the PME kernels.

        map<string, string> pmeDefines;
        pmeDefines["EWALD_ALPHA"] = cu.doubleToString(alpha);
        pmeDefines["PME_ORDER"] = cu.intToString(PmeOrder);
        pmeDefines["NUM_ATOMS"] = cu.intToString(numMultipoles);
        pmeDefines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        pmeDefines["EPSILON_FACTOR"] = cu.doubleToString(138.9354558456);
        pmeDefines["GRID_SIZE_X"] = cu.intToString(gridSizeX);
        pmeDefines["GRID_SIZE_Y"] = cu.intToString(gridSizeY);
        pmeDefines["GRID_SIZE_Z"] = cu.intToString(gridSizeZ);
        pmeDefines["M_PI"] = cu.doubleToString(M_PI);
        pmeDefines["SQRT_PI"] = cu.doubleToString(sqrt(M_PI));
        if (force.getPolarizationType() == MBPolMultipoleForce::Direct)
            pmeDefines["DIRECT_POLARIZATION"] = "";
        CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaMBPolKernelSources::multipolePme, pmeDefines);
        pmeGridIndexKernel = cu.getKernel(module, "findAtomGridIndex");
        pmeSpreadFixedMultipolesKernel = cu.getKernel(module, "gridSpreadFixedMultipoles");
        pmeSpreadInducedDipolesKernel = cu.getKernel(module, "gridSpreadInducedDipoles");
        pmeFinishSpreadChargeKernel = cu.getKernel(module, "finishSpreadCharge");
        pmeConvolutionKernel = cu.getKernel(module, "reciprocalConvolution");
        pmeFixedPotentialKernel = cu.getKernel(module, "computeFixedPotentialFromGrid");
        pmeInducedPotentialKernel = cu.getKernel(module, "computeInducedPotentialFromGrid");
        pmeFixedForceKernel = cu.getKernel(module, "computeFixedMultipoleForceAndEnergy");
        pmeInducedForceKernel = cu.getKernel(module, "computeInducedDipoleForceAndEnergy");
        pmeRecordInducedFieldDipolesKernel = cu.getKernel(module, "recordInducedFieldDipoles");
        cuFuncSetCacheConfig(pmeSpreadFixedMultipolesKernel, CU_FUNC_CACHE_PREFER_L1);
        cuFuncSetCacheConfig(pmeSpreadInducedDipolesKernel, CU_FUNC_CACHE_PREFER_L1);
        cuFuncSetCacheConfig(pmeFixedPotentialKernel, CU_FUNC_CACHE_PREFER_L1);
        cuFuncSetCacheConfig(pmeInducedPotentialKernel, CU_FUNC_CACHE_PREFER_L1);

        // Create required data structures.

        int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        pmeGrid = new CudaArray(cu, gridSizeX*gridSizeY*gridSizeZ, 2*elementSize, "pmeGrid");
        cu.addAutoclearBuffer(*pmeGrid);
        pmeBsplineModuliX = new CudaArray(cu, gridSizeX, elementSize, "pmeBsplineModuliX");
        pmeBsplineModuliY = new CudaArray(cu, gridSizeY, elementSize, "pmeBsplineModuliY");
        pmeBsplineModuliZ = new CudaArray(cu, gridSizeZ, elementSize, "pmeBsplineModuliZ");
        pmeIgrid = CudaArray::create<int4>(cu, numMultipoles, "pmeIgrid");
        pmePhi = new CudaArray(cu, 20*numMultipoles, elementSize, "pmePhi");
        pmePhid = new CudaArray(cu, 10*numMultipoles, elementSize, "pmePhid");
        pmePhip = new CudaArray(cu, 10*numMultipoles, elementSize, "pmePhip");
        pmePhidp = new CudaArray(cu, 20*numMultipoles, elementSize, "pmePhidp");
        pmeAtomRange = CudaArray::create<int>(cu, gridSizeX*gridSizeY*gridSizeZ+1, "pmeAtomRange");
        pmeAtomGridIndex = CudaArray::create<int2>(cu, numMultipoles, "pmeAtomGridIndex");
        sort = new CudaSort(cu, new SortTrait(), cu.getNumAtoms());
        cufftResult result = cufftPlan3d(&fft, gridSizeX, gridSizeY, gridSizeZ, cu.getUseDoublePrecision() ? CUFFT_Z2Z : CUFFT_C2C);
        if (result != CUFFT_SUCCESS)
            throw OpenMMException("Error initializing FFT: "+cu.intToString(result));
        hasInitializedFFT = true;

        // Initialize the b-spline moduli.

        double data[PmeOrder];
        double x = 0.0;
        data[0] = 1.0 - x;
        data[1] = x;
        for (int i = 2; i < PmeOrder; i++) {
            double denom = 1.0/i;
            data[i] = x*data[i-1]*denom;
            for (int j = 1; j < i; j++)
                data[i-j] = ((x+j)*data[i-j-1] + ((i-j+1)-x)*data[i-j])*denom;
            data[0] = (1.0-x)*data[0]*denom;
        }
        int maxSize = max(max(gridSizeX, gridSizeY), gridSizeZ);
        vector<double> bsplines_data(maxSize+1, 0.0);
        for (int i = 2; i <= PmeOrder+1; i++)
            bsplines_data[i] = data[i-2];
        for (int dim = 0; dim < 3; dim++) {
            int ndata = (dim == 0 ? gridSizeX : dim == 1 ? gridSizeY : gridSizeZ);
            vector<double> moduli(ndata);

            // get the modulus of the discrete Fourier transform

            double factor = 2.0*M_PI/ndata;
            for (int i = 0; i < ndata; i++) {
                double sc = 0.0;
                double ss = 0.0;
                for (int j = 1; j <= ndata; j++) {
                    double arg = factor*i*(j-1);
                    sc += bsplines_data[j]*cos(arg);
                    ss += bsplines_data[j]*sin(arg);
                }
                moduli[i] = sc*sc+ss*ss;
            }

            // Fix for exponential Euler spline interpolation failure.

            double eps = 1.0e-7;
            if (moduli[0] < eps)
                moduli[0] = 0.9*moduli[1];
            for (int i = 1; i < ndata-1; i++)
                if (moduli[i] < eps)
                    moduli[i] = 0.9*(moduli[i-1]+moduli[i+1]);
            if (moduli[ndata-1] < eps)
                moduli[ndata-1] = 0.9*moduli[ndata-2];

            // Compute and apply the optimal zeta coefficient.

            int jcut = 50;
            for (int i = 1; i <= ndata; i++) {
                int k = i - 1;
                if (i > ndata/2)
                    k = k - ndata;
                double zeta;
                if (k == 0)
                    zeta = 1.0;
                else {
                    double sum1 = 1.0;
                    double sum2 = 1.0;
                    factor = M_PI*k/ndata;
                    for (int j = 1; j <= jcut; j++) {
                        double arg = factor/(factor+M_PI*j);
                        sum1 += pow(arg, PmeOrder);
                        sum2 += pow(arg, 2*PmeOrder);
                    }
                    for (int j = 1; j <= jcut; j++) {
                        double arg = factor/(factor-M_PI*j);
                        sum1 += pow(arg, PmeOrder);
                        sum2 += pow(arg, 2*PmeOrder);
                    }
                    zeta = sum2/sum1;
                }
                moduli[i-1] = moduli[i-1]*zeta*zeta;
            }
            if (cu.getUseDoublePrecision()) {
                if (dim == 0)
                    pmeBsplineModuliX->upload(moduli);
                else if (dim == 1)
                    pmeBsplineModuliY->upload(moduli);
                else
                    pmeBsplineModuliZ->upload(moduli);
            }
            else {
                vector<float> modulif(ndata);
                for (int i = 0; i < ndata; i++)
                    modulif[i] = (float) moduli[i];
                if (dim == 0)
                    pmeBsplineModuliX->upload(modulif);
                else if (dim == 1)
                    pmeBsplineModuliY->upload(modulif);
                else
                    pmeBsplineModuliZ->upload(modulif);
            }
        }
    }

    // Add an interaction to the default nonbonded kernel.  This doesn't actually do any calculations.  It's
    // just so that CudaNonbondedUtilities will build the exclusion flags and maintain the neighbor list.
    
    cu.getNonbondedUtilities().addInteraction(usePME, usePME, true, force.getCutoffDistance(), exclusions, "", force.getForceGroup());
    cu.getNonbondedUtilities().setUsePadding(false);
    cu.addForce(new ForceInfo(force));
}

void CudaCalcMBPolMultipoleForceKernel::initializeScaleFactors() {
    hasInitializedScaleFactors = true;
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    
    // Figure out the covalent flag values to use for each atom pair.

    vector<ushort2> exclusionTiles;
    nb.getExclusionTiles().download(exclusionTiles);
    map<pair<int, int>, int> exclusionTileMap;
    for (int i = 0; i < (int) exclusionTiles.size(); i++) {
        ushort2 tile = exclusionTiles[i];
        exclusionTileMap[make_pair(tile.x, tile.y)] = i;
    }
    covalentFlags = CudaArray::create<uint2>(cu, nb.getExclusions().getSize(), "covalentFlags");
    vector<uint2> covalentFlagsVec(nb.getExclusions().getSize(), make_uint2(0, 0));
    for (int i = 0; i < (int) covalentFlagValues.size(); i++) {
        int atom1 = covalentFlagValues[i].x;
        int atom2 = covalentFlagValues[i].y;
        int value = covalentFlagValues[i].z;
        int x = atom1/CudaContext::TileSize;
        int offset1 = atom1-x*CudaContext::TileSize;
        int y = atom2/CudaContext::TileSize;
        int offset2 = atom2-y*CudaContext::TileSize;
        int f1 = (value == 0 || value == 1 ? 1 : 0);
        int f2 = (value == 0 || value == 2 ? 1 : 0);
        if (x == y) {
            int index = exclusionTileMap[make_pair(x, y)]*CudaContext::TileSize;
            covalentFlagsVec[index+offset1].x |= f1<<offset2;
            covalentFlagsVec[index+offset1].y |= f2<<offset2;
            covalentFlagsVec[index+offset2].x |= f1<<offset1;
            covalentFlagsVec[index+offset2].y |= f2<<offset1;
        }
        else if (x > y) {
            int index = exclusionTileMap[make_pair(x, y)]*CudaContext::TileSize;
            covalentFlagsVec[index+offset1].x |= f1<<offset2;
            covalentFlagsVec[index+offset1].y |= f2<<offset2;
        }
        else {
            int index = exclusionTileMap[make_pair(y, x)]*CudaContext::TileSize;
            covalentFlagsVec[index+offset2].x |= f1<<offset1;
            covalentFlagsVec[index+offset2].y |= f2<<offset1;
        }
    }
    covalentFlags->upload(covalentFlagsVec);
    
    // Do the same for the polarization flags.
    
    polarizationGroupFlags = CudaArray::create<unsigned int>(cu, nb.getExclusions().getSize(), "polarizationGroupFlags");
    vector<unsigned int> polarizationGroupFlagsVec(nb.getExclusions().getSize(), 0);
    for (int i = 0; i < (int) polarizationFlagValues.size(); i++) {
        int atom1 = polarizationFlagValues[i].x;
        int atom2 = polarizationFlagValues[i].y;
        int x = atom1/CudaContext::TileSize;
        int offset1 = atom1-x*CudaContext::TileSize;
        int y = atom2/CudaContext::TileSize;
        int offset2 = atom2-y*CudaContext::TileSize;
        if (x == y) {
            int index = exclusionTileMap[make_pair(x, y)]*CudaContext::TileSize;
            polarizationGroupFlagsVec[index+offset1] |= 1<<offset2;
            polarizationGroupFlagsVec[index+offset2] |= 1<<offset1;
        }
        else if (x > y) {
            int index = exclusionTileMap[make_pair(x, y)]*CudaContext::TileSize;
            polarizationGroupFlagsVec[index+offset1] |= 1<<offset2;
        }
        else {
            int index = exclusionTileMap[make_pair(y, x)]*CudaContext::TileSize;
            polarizationGroupFlagsVec[index+offset2] |= 1<<offset1;
        }
    }
    polarizationGroupFlags->upload(polarizationGroupFlagsVec);
}

double CudaCalcMBPolMultipoleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (!hasInitializedScaleFactors) {
        initializeScaleFactors();
        for (int i = 0; i < (int) context.getForceImpls().size() && gkKernel == NULL; i++) {
            MBPolGeneralizedKirkwoodForceImpl* gkImpl = dynamic_cast<MBPolGeneralizedKirkwoodForceImpl*>(context.getForceImpls()[i]);
            if (gkImpl != NULL)
                gkKernel = dynamic_cast<CudaCalcMBPolGeneralizedKirkwoodForceKernel*>(&gkImpl->getKernel().getImpl());
        }
    }
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    
    // Compute the lab frame moments.

    void* computeMomentsArgs[] = {&cu.getPosq().getDevicePointer(), &multipoleParticles->getDevicePointer(),
        &molecularDipoles->getDevicePointer(), &molecularQuadrupoles->getDevicePointer(),
        &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer()};
    cu.executeKernel(computeMomentsKernel, computeMomentsArgs, cu.getNumAtoms());
    int startTileIndex = nb.getStartTileIndex();
    int numTileIndices = nb.getNumTiles();
    int numForceThreadBlocks = nb.getNumForceThreadBlocks();
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    void* npt = NULL;
    if (pmeGrid == NULL) {
        // Compute induced dipoles.
        
        if (gkKernel == NULL) {
            void* computeFixedFieldArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(), &cu.getPosq().getDevicePointer(),
                &covalentFlags->getDevicePointer(), &polarizationGroupFlags->getDevicePointer(), &nb.getExclusionTiles().getDevicePointer(), &startTileIndex, &numTileIndices,
                &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer(), &dampingAndThole->getDevicePointer()};
            cu.executeKernel(computeFixedFieldKernel, computeFixedFieldArgs, numForceThreadBlocks*fixedFieldThreads, fixedFieldThreads);
            void* recordInducedDipolesArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(),
                &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(), &polarizability->getDevicePointer()};
            cu.executeKernel(recordInducedDipolesKernel, recordInducedDipolesArgs, cu.getNumAtoms());
        }
        else {
            gkKernel->computeBornRadii();
            void* computeFixedFieldArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(), &cu.getPosq().getDevicePointer(),
                &covalentFlags->getDevicePointer(), &polarizationGroupFlags->getDevicePointer(), &nb.getExclusionTiles().getDevicePointer(), &startTileIndex, &numTileIndices,
                &gkKernel->getBornRadii()->getDevicePointer(), &gkKernel->getField()->getDevicePointer(),
                &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer(), &dampingAndThole->getDevicePointer()};
            cu.executeKernel(computeFixedFieldKernel, computeFixedFieldArgs, numForceThreadBlocks*fixedFieldThreads, fixedFieldThreads);
            void* recordInducedDipolesArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(),
                &gkKernel->getField()->getDevicePointer(), &gkKernel->getInducedDipoles()->getDevicePointer(),
                &gkKernel->getInducedDipolesPolar()->getDevicePointer(), &inducedDipole->getDevicePointer(),
                &inducedDipolePolar->getDevicePointer(), &polarizability->getDevicePointer()};
            cu.executeKernel(recordInducedDipolesKernel, recordInducedDipolesArgs, cu.getNumAtoms());
        }
        
        // Iterate until the dipoles converge.
        
        vector<float2> errors;
        for (int i = 0; i < maxInducedIterations; i++) {
            cu.clearBuffer(*inducedField);
            cu.clearBuffer(*inducedFieldPolar);
            if (gkKernel == NULL) {
                void* computeInducedFieldArgs[] = {&inducedField->getDevicePointer(), &inducedFieldPolar->getDevicePointer(), &cu.getPosq().getDevicePointer(),
                    &nb.getExclusionTiles().getDevicePointer(), &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(), &startTileIndex, &numTileIndices,
                    &dampingAndThole->getDevicePointer()};
                cu.executeKernel(computeInducedFieldKernel, computeInducedFieldArgs, numForceThreadBlocks*inducedFieldThreads, inducedFieldThreads);
            }
            else {
                cu.clearBuffer(*gkKernel->getInducedField());
                cu.clearBuffer(*gkKernel->getInducedFieldPolar());
                void* computeInducedFieldArgs[] = {&inducedField->getDevicePointer(), &inducedFieldPolar->getDevicePointer(), &cu.getPosq().getDevicePointer(),
                    &nb.getExclusionTiles().getDevicePointer(), &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(), &startTileIndex, &numTileIndices,
                    &gkKernel->getInducedField()->getDevicePointer(), &gkKernel->getInducedFieldPolar()->getDevicePointer(),
                    &gkKernel->getInducedDipoles()->getDevicePointer(), &gkKernel->getInducedDipolesPolar()->getDevicePointer(),
                    &gkKernel->getBornRadii()->getDevicePointer(), &dampingAndThole->getDevicePointer()};
                cu.executeKernel(computeInducedFieldKernel, computeInducedFieldArgs, numForceThreadBlocks*inducedFieldThreads, inducedFieldThreads);
                void* updateInducedGkFieldArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(),
                    &gkKernel->getField()->getDevicePointer(), &gkKernel->getInducedField()->getDevicePointer(),
                    &gkKernel->getInducedFieldPolar()->getDevicePointer(), &gkKernel->getInducedDipoles()->getDevicePointer(),
                    &gkKernel->getInducedDipolesPolar()->getDevicePointer(), &polarizability->getDevicePointer(), &inducedDipoleErrors->getDevicePointer()};
                cu.executeKernel(updateInducedFieldKernel, updateInducedGkFieldArgs, cu.getNumThreadBlocks()*cu.ThreadBlockSize, cu.ThreadBlockSize, cu.ThreadBlockSize*elementSize*2);
            }
            void* updateInducedFieldArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(), &npt, &inducedField->getDevicePointer(),
                &inducedFieldPolar->getDevicePointer(), &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(),
                &polarizability->getDevicePointer(), &inducedDipoleErrors->getDevicePointer()};
            cu.executeKernel(updateInducedFieldKernel, updateInducedFieldArgs, cu.getNumThreadBlocks()*cu.ThreadBlockSize, cu.ThreadBlockSize, cu.ThreadBlockSize*elementSize*2);
            inducedDipoleErrors->download(errors);
            double total1 = 0.0, total2 = 0.0;
            for (int j = 0; j < (int) errors.size(); j++) {
                total1 += errors[j].x;
                total2 += errors[j].y;
            }
            if (48.033324*sqrt(max(total1, total2)/cu.getNumAtoms()) < inducedEpsilon)
                break;
        }
        
        // Compute electrostatic force.
        
        void* electrostaticsArgs[] = {&cu.getForce().getDevicePointer(), &torque->getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(),
            &cu.getPosq().getDevicePointer(), &covalentFlags->getDevicePointer(), &polarizationGroupFlags->getDevicePointer(),
            &nb.getExclusionTiles().getDevicePointer(), &startTileIndex, &numTileIndices,
            &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer(), &inducedDipole->getDevicePointer(),
            &inducedDipolePolar->getDevicePointer(), &dampingAndThole->getDevicePointer()};
        cu.executeKernel(electrostaticsKernel, electrostaticsArgs, numForceThreadBlocks*electrostaticsThreads, electrostaticsThreads);
        if (gkKernel != NULL)
            gkKernel->finishComputation(*torque, *labFrameDipoles, *labFrameQuadrupoles, *inducedDipole, *inducedDipolePolar, *dampingAndThole, *covalentFlags, *polarizationGroupFlags);
    }
    else {
        // Reciprocal space calculation.
        
        unsigned int maxTiles = nb.getInteractingTiles().getSize();
        void* gridIndexArgs[] = {&cu.getPosq().getDevicePointer(), &pmeAtomGridIndex->getDevicePointer(),
            cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer()};
        cu.executeKernel(pmeGridIndexKernel, gridIndexArgs, cu.getNumAtoms(), cu.ThreadBlockSize, cu.ThreadBlockSize*PmeOrder*PmeOrder*elementSize);
        sort->sort(*pmeAtomGridIndex);
        void* pmeSpreadFixedMultipolesArgs[] = {&cu.getPosq().getDevicePointer(), &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer(),
            &pmeGrid->getDevicePointer(), &pmeAtomGridIndex->getDevicePointer(),  cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer()};
        cu.executeKernel(pmeSpreadFixedMultipolesKernel, pmeSpreadFixedMultipolesArgs, cu.getNumAtoms());
        void* finishSpreadArgs[] = {&pmeGrid->getDevicePointer()};
        if (cu.getUseDoublePrecision())
            cu.executeKernel(pmeFinishSpreadChargeKernel, finishSpreadArgs, pmeGrid->getSize());
        if (cu.getUseDoublePrecision())
            cufftExecZ2Z(fft, (double2*) pmeGrid->getDevicePointer(), (double2*) pmeGrid->getDevicePointer(), CUFFT_FORWARD);
        else
            cufftExecC2C(fft, (float2*) pmeGrid->getDevicePointer(), (float2*) pmeGrid->getDevicePointer(), CUFFT_FORWARD);
        void* pmeConvolutionArgs[] = {&pmeGrid->getDevicePointer(), &pmeBsplineModuliX->getDevicePointer(), &pmeBsplineModuliY->getDevicePointer(),
            &pmeBsplineModuliZ->getDevicePointer(), cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer()};
        cu.executeKernel(pmeConvolutionKernel, pmeConvolutionArgs, cu.getNumAtoms());
        if (cu.getUseDoublePrecision())
            cufftExecZ2Z(fft, (double2*) pmeGrid->getDevicePointer(), (double2*) pmeGrid->getDevicePointer(), CUFFT_INVERSE);
        else
            cufftExecC2C(fft, (float2*) pmeGrid->getDevicePointer(), (float2*) pmeGrid->getDevicePointer(), CUFFT_INVERSE);
        void* pmeFixedPotentialArgs[] = {&pmeGrid->getDevicePointer(), &pmePhi->getDevicePointer(), &field->getDevicePointer(),
            &fieldPolar ->getDevicePointer(), &cu.getPosq().getDevicePointer(), &labFrameDipoles->getDevicePointer(),
            cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer(), &pmeAtomGridIndex->getDevicePointer()};
        cu.executeKernel(pmeFixedPotentialKernel, pmeFixedPotentialArgs, cu.getNumAtoms());
        void* pmeFixedForceArgs[] = {&cu.getPosq().getDevicePointer(), &cu.getForce().getDevicePointer(), &torque->getDevicePointer(),
            &cu.getEnergyBuffer().getDevicePointer(), &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer(),
            &pmePhi->getDevicePointer(), cu.getInvPeriodicBoxSizePointer()};
        cu.executeKernel(pmeFixedForceKernel, pmeFixedForceArgs, cu.getNumAtoms());
        
        // Direct space calculation.
        
        void* computeFixedFieldArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(), &cu.getPosq().getDevicePointer(),
            &covalentFlags->getDevicePointer(), &polarizationGroupFlags->getDevicePointer(), &nb.getExclusionTiles().getDevicePointer(), &startTileIndex, &numTileIndices,
            &nb.getInteractingTiles().getDevicePointer(), &nb.getInteractionCount().getDevicePointer(), cu.getPeriodicBoxSizePointer(),
            cu.getInvPeriodicBoxSizePointer(), &maxTiles, &nb.getBlockCenters().getDevicePointer(), &nb.getInteractingAtoms().getDevicePointer(),
            &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer(), &dampingAndThole->getDevicePointer()};
        cu.executeKernel(computeFixedFieldKernel, computeFixedFieldArgs, numForceThreadBlocks*fixedFieldThreads, fixedFieldThreads);
        void* recordInducedDipolesArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(),
            &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(), &polarizability->getDevicePointer()};
        cu.executeKernel(recordInducedDipolesKernel, recordInducedDipolesArgs, cu.getNumAtoms());

        // Reciprocal space calculation for the induced dipoles.

        cu.clearBuffer(*pmeGrid);
        void* pmeSpreadInducedDipolesArgs[] = {&cu.getPosq().getDevicePointer(), &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(),
            &pmeGrid->getDevicePointer(), &pmeAtomGridIndex->getDevicePointer(), cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer()};
        cu.executeKernel(pmeSpreadInducedDipolesKernel, pmeSpreadInducedDipolesArgs, cu.getNumAtoms());
        if (cu.getUseDoublePrecision())
            cu.executeKernel(pmeFinishSpreadChargeKernel, finishSpreadArgs, pmeGrid->getSize());
        if (cu.getUseDoublePrecision())
            cufftExecZ2Z(fft, (double2*) pmeGrid->getDevicePointer(), (double2*) pmeGrid->getDevicePointer(), CUFFT_FORWARD);
        else
            cufftExecC2C(fft, (float2*) pmeGrid->getDevicePointer(), (float2*) pmeGrid->getDevicePointer(), CUFFT_FORWARD);
        cu.executeKernel(pmeConvolutionKernel, pmeConvolutionArgs, cu.getNumAtoms());
        if (cu.getUseDoublePrecision())
            cufftExecZ2Z(fft, (double2*) pmeGrid->getDevicePointer(), (double2*) pmeGrid->getDevicePointer(), CUFFT_INVERSE);
        else
            cufftExecC2C(fft, (float2*) pmeGrid->getDevicePointer(), (float2*) pmeGrid->getDevicePointer(), CUFFT_INVERSE);
        void* pmeInducedPotentialArgs[] = {&pmeGrid->getDevicePointer(), &pmePhid->getDevicePointer(), &pmePhip->getDevicePointer(),
            &pmePhidp->getDevicePointer(), &cu.getPosq().getDevicePointer(), cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer(),
            &pmeAtomGridIndex->getDevicePointer()};
        cu.executeKernel(pmeInducedPotentialKernel, pmeInducedPotentialArgs, cu.getNumAtoms());
        
        // Iterate until the dipoles converge.
        
        vector<float2> errors;
        for (int i = 0; i < maxInducedIterations; i++) {
            cu.clearBuffer(*inducedField);
            cu.clearBuffer(*inducedFieldPolar);
            void* computeInducedFieldArgs[] = {&inducedField->getDevicePointer(), &inducedFieldPolar->getDevicePointer(), &cu.getPosq().getDevicePointer(),
                &nb.getExclusionTiles().getDevicePointer(), &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(), &startTileIndex, &numTileIndices,
                &nb.getInteractingTiles().getDevicePointer(), &nb.getInteractionCount().getDevicePointer(), cu.getPeriodicBoxSizePointer(),
                cu.getInvPeriodicBoxSizePointer(), &maxTiles, &nb.getBlockCenters().getDevicePointer(), &nb.getInteractingAtoms().getDevicePointer(),
                &dampingAndThole->getDevicePointer()};
            cu.executeKernel(computeInducedFieldKernel, computeInducedFieldArgs, numForceThreadBlocks*inducedFieldThreads, inducedFieldThreads);
            cu.clearBuffer(*pmeGrid);
            cu.executeKernel(pmeSpreadInducedDipolesKernel, pmeSpreadInducedDipolesArgs, cu.getNumAtoms());
            if (cu.getUseDoublePrecision())
                cu.executeKernel(pmeFinishSpreadChargeKernel, finishSpreadArgs, pmeGrid->getSize());
            if (cu.getUseDoublePrecision())
                cufftExecZ2Z(fft, (double2*) pmeGrid->getDevicePointer(), (double2*) pmeGrid->getDevicePointer(), CUFFT_FORWARD);
            else
                cufftExecC2C(fft, (float2*) pmeGrid->getDevicePointer(), (float2*) pmeGrid->getDevicePointer(), CUFFT_FORWARD);
            cu.executeKernel(pmeConvolutionKernel, pmeConvolutionArgs, cu.getNumAtoms());
            if (cu.getUseDoublePrecision())
                cufftExecZ2Z(fft, (double2*) pmeGrid->getDevicePointer(), (double2*) pmeGrid->getDevicePointer(), CUFFT_INVERSE);
            else
                cufftExecC2C(fft, (float2*) pmeGrid->getDevicePointer(), (float2*) pmeGrid->getDevicePointer(), CUFFT_INVERSE);
            cu.executeKernel(pmeInducedPotentialKernel, pmeInducedPotentialArgs, cu.getNumAtoms());
            void* pmeRecordInducedFieldDipolesArgs[] = {&pmePhid->getDevicePointer(), &pmePhip->getDevicePointer(),
                &inducedField->getDevicePointer(), &inducedFieldPolar->getDevicePointer(), cu.getInvPeriodicBoxSizePointer()};
            cu.executeKernel(pmeRecordInducedFieldDipolesKernel, pmeRecordInducedFieldDipolesArgs, cu.getNumAtoms());
            void* updateInducedFieldArgs[] = {&field->getDevicePointer(), &fieldPolar->getDevicePointer(), &npt, &inducedField->getDevicePointer(),
                &inducedFieldPolar->getDevicePointer(), &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(),
                &polarizability->getDevicePointer(), &inducedDipoleErrors->getDevicePointer()};
            cu.executeKernel(updateInducedFieldKernel, updateInducedFieldArgs, cu.getNumThreadBlocks()*cu.ThreadBlockSize, cu.ThreadBlockSize, cu.ThreadBlockSize*elementSize*2);
            inducedDipoleErrors->download(errors);
            double total1 = 0.0, total2 = 0.0;
            for (int j = 0; j < (int) errors.size(); j++) {
                total1 += errors[j].x;
                total2 += errors[j].y;
            }
            if (48.033324*sqrt(max(total1, total2)/cu.getNumAtoms()) < inducedEpsilon)
                break;
        }
        
        // Compute electrostatic force.
        
        void* electrostaticsArgs[] = {&cu.getForce().getDevicePointer(), &torque->getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(),
            &cu.getPosq().getDevicePointer(), &covalentFlags->getDevicePointer(), &polarizationGroupFlags->getDevicePointer(),
            &nb.getExclusionTiles().getDevicePointer(), &startTileIndex, &numTileIndices,
            &nb.getInteractingTiles().getDevicePointer(), &nb.getInteractionCount().getDevicePointer(),
            cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer(), &maxTiles, &nb.getBlockCenters().getDevicePointer(), &nb.getInteractingAtoms().getDevicePointer(),
            &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer(), &inducedDipole->getDevicePointer(),
            &inducedDipolePolar->getDevicePointer(), &dampingAndThole->getDevicePointer()};
        cu.executeKernel(electrostaticsKernel, electrostaticsArgs, numForceThreadBlocks*electrostaticsThreads, electrostaticsThreads);
        void* pmeInducedForceArgs[] = {&cu.getPosq().getDevicePointer(), &cu.getForce().getDevicePointer(), &torque->getDevicePointer(),
            &cu.getEnergyBuffer().getDevicePointer(), &labFrameDipoles->getDevicePointer(), &labFrameQuadrupoles->getDevicePointer(),
            &inducedDipole->getDevicePointer(), &inducedDipolePolar->getDevicePointer(), &pmePhi->getDevicePointer(), &pmePhid->getDevicePointer(),
            &pmePhip->getDevicePointer(), &pmePhidp->getDevicePointer(), cu.getInvPeriodicBoxSizePointer()};
        cu.executeKernel(pmeInducedForceKernel, pmeInducedForceArgs, cu.getNumAtoms());
    }

    // Map torques to force.

    void* mapTorqueArgs[] = {&cu.getForce().getDevicePointer(), &torque->getDevicePointer(),
        &cu.getPosq().getDevicePointer(), &multipoleParticles->getDevicePointer()};
    cu.executeKernel(mapTorqueKernel, mapTorqueArgs, cu.getNumAtoms());
    
    // Record the current atom positions so we can tell later if they have changed.
    
    cu.getPosq().copyTo(*lastPositions);
    multipolesAreValid = true;
    return 0.0;
}

void CudaCalcMBPolMultipoleForceKernel::ensureMultipolesValid(ContextImpl& context) {
    if (multipolesAreValid) {
        int numParticles = cu.getNumAtoms();
        if (cu.getUseDoublePrecision()) {
            vector<double4> pos1, pos2;
            cu.getPosq().download(pos1);
            lastPositions->download(pos2);
            for (int i = 0; i < numParticles; i++)
                if (pos1[i].x != pos2[i].x || pos1[i].y != pos2[i].y || pos1[i].z != pos2[i].z) {
                    multipolesAreValid = false;
                    break;
                }
        }
        else {
            vector<float4> pos1, pos2;
            cu.getPosq().download(pos1);
            lastPositions->download(pos2);
            for (int i = 0; i < numParticles; i++)
                if (pos1[i].x != pos2[i].x || pos1[i].y != pos2[i].y || pos1[i].z != pos2[i].z) {
                    multipolesAreValid = false;
                    break;
                }
        }
    }
    if (!multipolesAreValid)
        context.calcForcesAndEnergy(false, false, -1);
}

void CudaCalcMBPolMultipoleForceKernel::getElectrostaticPotential(ContextImpl& context, const vector<Vec3>& inputGrid, vector<double>& outputElectrostaticPotential) {
    ensureMultipolesValid(context);
    int numPoints = inputGrid.size();
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    CudaArray points(cu, numPoints, 4*elementSize, "points");
    CudaArray potential(cu, numPoints, elementSize, "potential");
    
    // Copy the grid points to the GPU.
    
    if (cu.getUseDoublePrecision()) {
        vector<double4> p(numPoints);
        for (int i = 0; i < numPoints; i++)
            p[i] = make_double4(inputGrid[i][0], inputGrid[i][1], inputGrid[i][2], 0);
        points.upload(p);
    }
    else {
        vector<float4> p(numPoints);
        for (int i = 0; i < numPoints; i++)
            p[i] = make_float4((float) inputGrid[i][0], (float) inputGrid[i][1], (float) inputGrid[i][2], 0);
        points.upload(p);
    }
    
    // Compute the potential.
    
    void* computePotentialArgs[] = {&cu.getPosq().getDevicePointer(), &labFrameDipoles->getDevicePointer(),
        &labFrameQuadrupoles->getDevicePointer(), &inducedDipole->getDevicePointer(), &points.getDevicePointer(),
        &potential.getDevicePointer(), &numPoints, cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer()};
    int blockSize = 128;
    cu.executeKernel(computePotentialKernel, computePotentialArgs, numPoints, blockSize, blockSize*15*elementSize);
    outputElectrostaticPotential.resize(numPoints);
    if (cu.getUseDoublePrecision())
        potential.download(outputElectrostaticPotential);
    else {
        vector<float> p(numPoints);
        potential.download(p);
        for (int i = 0; i < numPoints; i++)
            outputElectrostaticPotential[i] = p[i];
    }
}

template <class T, class T4, class M4>
void CudaCalcMBPolMultipoleForceKernel::computeSystemMultipoleMoments(ContextImpl& context, vector<double>& outputMultipoleMoments) {
    // Compute the local coordinates relative to the center of mass.
    int numAtoms = cu.getNumAtoms();
    vector<T4> posq;
    vector<M4> velm;
    cu.getPosq().download(posq);
    cu.getVelm().download(velm);
    double totalMass = 0.0;
    Vec3 centerOfMass(0, 0, 0);
    for (int i = 0; i < numAtoms; i++) {
        double mass = (velm[i].w > 0 ? 1.0/velm[i].w : 0.0);
        totalMass += mass;
        centerOfMass[0] += mass*posq[i].x;
        centerOfMass[1] += mass*posq[i].y;
        centerOfMass[2] += mass*posq[i].z;
    }
    if (totalMass > 0.0) {
        centerOfMass[0] /= totalMass;
        centerOfMass[1] /= totalMass;
        centerOfMass[2] /= totalMass;
    }
    vector<double4> posqLocal(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        posqLocal[i].x = posq[i].x - centerOfMass[0];
        posqLocal[i].y = posq[i].y - centerOfMass[1];
        posqLocal[i].z = posq[i].z - centerOfMass[2];
        posqLocal[i].w = posq[i].w;
    }

    // Compute the multipole moments.
    
    double totalCharge = 0.0;
    double xdpl = 0.0;
    double ydpl = 0.0;
    double zdpl = 0.0;
    double xxqdp = 0.0;
    double xyqdp = 0.0;
    double xzqdp = 0.0;
    double yxqdp = 0.0;
    double yyqdp = 0.0;
    double yzqdp = 0.0;
    double zxqdp = 0.0;
    double zyqdp = 0.0;
    double zzqdp = 0.0;
    vector<T> labDipoleVec, inducedDipoleVec, quadrupoleVec;
    labFrameDipoles->download(labDipoleVec);
    inducedDipole->download(inducedDipoleVec);
    labFrameQuadrupoles->download(quadrupoleVec);
    for (int i = 0; i < numAtoms; i++) {
        totalCharge += posqLocal[i].w;
        double netDipoleX = (labDipoleVec[3*i]    + inducedDipoleVec[3*i]);
        double netDipoleY = (labDipoleVec[3*i+1]  + inducedDipoleVec[3*i+1]);
        double netDipoleZ = (labDipoleVec[3*i+2]  + inducedDipoleVec[3*i+2]);
        xdpl += posqLocal[i].x*posqLocal[i].w + netDipoleX;
        ydpl += posqLocal[i].y*posqLocal[i].w + netDipoleY;
        zdpl += posqLocal[i].z*posqLocal[i].w + netDipoleZ;
        xxqdp += posqLocal[i].x*posqLocal[i].x*posqLocal[i].w + 2*posqLocal[i].x*netDipoleX;
        xyqdp += posqLocal[i].x*posqLocal[i].y*posqLocal[i].w + posqLocal[i].x*netDipoleY + posqLocal[i].y*netDipoleX;
        xzqdp += posqLocal[i].x*posqLocal[i].z*posqLocal[i].w + posqLocal[i].x*netDipoleZ + posqLocal[i].z*netDipoleX;
        yxqdp += posqLocal[i].y*posqLocal[i].x*posqLocal[i].w + posqLocal[i].y*netDipoleX + posqLocal[i].x*netDipoleY;
        yyqdp += posqLocal[i].y*posqLocal[i].y*posqLocal[i].w + 2*posqLocal[i].y*netDipoleY;
        yzqdp += posqLocal[i].y*posqLocal[i].z*posqLocal[i].w + posqLocal[i].y*netDipoleZ + posqLocal[i].z*netDipoleY;
        zxqdp += posqLocal[i].z*posqLocal[i].x*posqLocal[i].w + posqLocal[i].z*netDipoleX + posqLocal[i].x*netDipoleZ;
        zyqdp += posqLocal[i].z*posqLocal[i].y*posqLocal[i].w + posqLocal[i].z*netDipoleY + posqLocal[i].y*netDipoleZ;
        zzqdp += posqLocal[i].z*posqLocal[i].z*posqLocal[i].w + 2*posqLocal[i].z*netDipoleZ;
    }

    // Convert the quadrupole from traced to traceless form.
 
    double qave = (xxqdp + yyqdp + zzqdp)/3;
    xxqdp = 1.5*(xxqdp-qave);
    xyqdp = 1.5*xyqdp;
    xzqdp = 1.5*xzqdp;
    yxqdp = 1.5*yxqdp;
    yyqdp = 1.5*(yyqdp-qave);
    yzqdp = 1.5*yzqdp;
    zxqdp = 1.5*zxqdp;
    zyqdp = 1.5*zyqdp;
    zzqdp = 1.5*(zzqdp-qave);

    // Add the traceless atomic quadrupoles to the total quadrupole moment.

    for (int i = 0; i < numAtoms; i++) {
        xxqdp = xxqdp + 3*quadrupoleVec[5*i];
        xyqdp = xyqdp + 3*quadrupoleVec[5*i+1];
        xzqdp = xzqdp + 3*quadrupoleVec[5*i+2];
        yxqdp = yxqdp + 3*quadrupoleVec[5*i+1];
        yyqdp = yyqdp + 3*quadrupoleVec[5*i+3];
        yzqdp = yzqdp + 3*quadrupoleVec[5*i+4];
        zxqdp = zxqdp + 3*quadrupoleVec[5*i+2];
        zyqdp = zyqdp + 3*quadrupoleVec[5*i+4];
        zzqdp = zzqdp + -3*(quadrupoleVec[5*i]+quadrupoleVec[5*i+3]);
    }
 
    double debye = 4.80321;
    outputMultipoleMoments.resize(13);
    outputMultipoleMoments[0] = totalCharge;
    outputMultipoleMoments[1] = 10.0*xdpl*debye;
    outputMultipoleMoments[2] = 10.0*ydpl*debye;
    outputMultipoleMoments[3] = 10.0*zdpl*debye;
    outputMultipoleMoments[4] = 100.0*xxqdp*debye;
    outputMultipoleMoments[5] = 100.0*xyqdp*debye;
    outputMultipoleMoments[6] = 100.0*xzqdp*debye;
    outputMultipoleMoments[7] = 100.0*yxqdp*debye;
    outputMultipoleMoments[8] = 100.0*yyqdp*debye;
    outputMultipoleMoments[9] = 100.0*yzqdp*debye;
    outputMultipoleMoments[10] = 100.0*zxqdp*debye;
    outputMultipoleMoments[11] = 100.0*zyqdp*debye;
    outputMultipoleMoments[12] = 100.0*zzqdp*debye;
}

void CudaCalcMBPolMultipoleForceKernel::getSystemMultipoleMoments(ContextImpl& context, vector<double>& outputMultipoleMoments) {
    ensureMultipolesValid(context);
    if (cu.getUseDoublePrecision())
        computeSystemMultipoleMoments<double, double4, double4>(context, outputMultipoleMoments);
    else if (cu.getUseMixedPrecision())
        computeSystemMultipoleMoments<float, float4, double4>(context, outputMultipoleMoments);
    else
        computeSystemMultipoleMoments<float, float4, float4>(context, outputMultipoleMoments);
}

void CudaCalcMBPolMultipoleForceKernel::copyParametersToContext(ContextImpl& context, const MBPolMultipoleForce& force) {
    // Make sure the new parameters are acceptable.
    
    cu.setAsCurrent();
    if (force.getNumMultipoles() != cu.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of multipoles has changed");
    
    // Record the per-multipole parameters.
    
    cu.getPosq().download(cu.getPinnedBuffer());
    float4* posqf = (float4*) cu.getPinnedBuffer();
    double4* posqd = (double4*) cu.getPinnedBuffer();
    vector<float2> dampingAndTholeVec;
    vector<float> polarizabilityVec;
    vector<float> molecularDipolesVec;
    vector<float> molecularQuadrupolesVec;
    vector<int4> multipoleParticlesVec;
    for (int i = 0; i < force.getNumMultipoles(); i++) {
        double charge, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getMultipoleParameters(i, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        if (cu.getUseDoublePrecision())
            posqd[i].w = charge;
        else
            posqf[i].w = (float) charge;
        dampingAndTholeVec.push_back(make_float2((float) damping, (float) thole));
        polarizabilityVec.push_back((float) polarity);
        multipoleParticlesVec.push_back(make_int4(atomX, atomY, atomZ, axisType));
        for (int j = 0; j < 3; j++)
            molecularDipolesVec.push_back((float) dipole[j]);
        molecularQuadrupolesVec.push_back((float) quadrupole[0]);
        molecularQuadrupolesVec.push_back((float) quadrupole[1]);
        molecularQuadrupolesVec.push_back((float) quadrupole[2]);
        molecularQuadrupolesVec.push_back((float) quadrupole[4]);
        molecularQuadrupolesVec.push_back((float) quadrupole[5]);
    }
    for (int i = force.getNumMultipoles(); i < cu.getPaddedNumAtoms(); i++) {
        dampingAndTholeVec.push_back(make_float2(0, 0));
        polarizabilityVec.push_back(0);
        multipoleParticlesVec.push_back(make_int4(0, 0, 0, 0));
        for (int j = 0; j < 3; j++)
            molecularDipolesVec.push_back(0);
        for (int j = 0; j < 5; j++)
            molecularQuadrupolesVec.push_back(0);
    }
    dampingAndThole->upload(dampingAndTholeVec);
    polarizability->upload(polarizabilityVec);
    multipoleParticles->upload(multipoleParticlesVec);
    molecularDipoles->upload(molecularDipolesVec);
    molecularQuadrupoles->upload(molecularQuadrupolesVec);
    cu.getPosq().upload(cu.getPinnedBuffer());
    cu.invalidateMolecules();
    multipolesAreValid = false;
}
