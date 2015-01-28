#ifndef MBPOL_OPENMM_CUDAKERNELS_H_
#define MBPOL_OPENMM_CUDAKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMMBPol                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2013 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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

#include "openmm/mbpolKernels.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "CudaArray.h"
#include "CudaContext.h"
#include "CudaSort.h"
#include <cufft.h>

namespace MBPolPlugin {

/**
 * This kernel is invoked by MBPolOneBodyForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolOneBodyForceKernel : public CalcMBPolOneBodyForceKernel {
public:
    CudaCalcMBPolOneBodyForceKernel(std::string name, 
                                          const Platform& platform,
                                          CudaContext& cu,
                                          const System& system);
    ~CudaCalcMBPolOneBodyForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolOneBodyForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolOneBodyForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolOneBodyForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolOneBodyForce& force);
private:
    class ForceInfo;
    int numOneBodys;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MBPolAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolAngleForceKernel : public CalcMBPolAngleForceKernel {
public:
    CudaCalcMBPolAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolAngleForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolAngleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolAngleForce& force);
private:
    class ForceInfo;
    int numAngles;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MBPolInPlaneAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolInPlaneAngleForceKernel : public CalcMBPolInPlaneAngleForceKernel {
public:
    CudaCalcMBPolInPlaneAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolInPlaneAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolInPlaneAngleForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolInPlaneAngleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolInPlaneAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolInPlaneAngleForce& force);
private:
    class ForceInfo;
    int numAngles;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MBPolPiTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolPiTorsionForceKernel : public CalcMBPolPiTorsionForceKernel {
public:
    CudaCalcMBPolPiTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolPiTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolPiTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolPiTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolPiTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolPiTorsionForce& force);
private:
    class ForceInfo;
    int numPiTorsions;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MBPolStretchBendForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolStretchBendForceKernel : public CalcMBPolStretchBendForceKernel {
public:
    CudaCalcMBPolStretchBendForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolStretchBendForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolStretchBendForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolStretchBendForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolStretchBendForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolStretchBendForce& force);
private:
    class ForceInfo;
    int numStretchBends;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MBPolOutOfPlaneBendForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolOutOfPlaneBendForceKernel : public CalcMBPolOutOfPlaneBendForceKernel {
public:
    CudaCalcMBPolOutOfPlaneBendForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolOutOfPlaneBendForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolOutOfPlaneBendForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolOutOfPlaneBendForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolOutOfPlaneBendForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolOutOfPlaneBendForce& force);
private:
    class ForceInfo;
    int numOutOfPlaneBends;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by MBPolTorsionTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolTorsionTorsionForceKernel : public CalcMBPolTorsionTorsionForceKernel {
public:
    CudaCalcMBPolTorsionTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolTorsionTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolTorsionTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolTorsionTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
private:
    class ForceInfo;
    int numTorsionTorsions;
    int numTorsionTorsionGrids;
    CudaContext& cu;
    const System& system;
    CudaArray* gridValues;
    CudaArray* gridParams;
    CudaArray* torsionParams;
};

/**
 * This kernel is invoked by MBPolMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolMultipoleForceKernel : public CalcMBPolMultipoleForceKernel {
public:
    CudaCalcMBPolMultipoleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolMultipoleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolMultipoleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Execute the kernel to calculate the electrostatic potential
     *
     * @param context        the context in which to execute this kernel
     * @param inputGrid      input grid coordinates
     * @param outputElectrostaticPotential output potential 
     */
    void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                   std::vector< double >& outputElectrostaticPotential );

   /** 
     * Get the system multipole moments
     *
     * @param context      context
     * @param outputMultipoleMoments (charge,
     *                                dipole_x, dipole_y, dipole_z,
     *                                quadrupole_xx, quadrupole_xy, quadrupole_xz,
     *                                quadrupole_yx, quadrupole_yy, quadrupole_yz,
     *                                quadrupole_zx, quadrupole_zy, quadrupole_zz )
     */
    void getSystemMultipoleMoments(ContextImpl& context, std::vector<double>& outputMultipoleMoments);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolMultipoleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolMultipoleForce& force);
private:
    class ForceInfo;
    class SortTrait : public CudaSort::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "INT_MIN";}
        const char* getMaxKey() const {return "INT_MAX";}
        const char* getMaxValue() const {return "make_int2(INT_MAX, INT_MAX)";}
        const char* getSortKey() const {return "value.y";}
    };
    void initializeScaleFactors();
    void ensureMultipolesValid(ContextImpl& context);
    template <class T, class T4, class M4> void computeSystemMultipoleMoments(ContextImpl& context, std::vector<double>& outputMultipoleMoments);
    int numMultipoles, maxInducedIterations;
    int fixedFieldThreads, inducedFieldThreads, electrostaticsThreads;
    double inducedEpsilon;
    bool hasInitializedScaleFactors, hasInitializedFFT, multipolesAreValid;
    CudaContext& cu;
    const System& system;
    std::vector<int3> covalentFlagValues;
    std::vector<int2> polarizationFlagValues;
    CudaArray* multipoleParticles;
    CudaArray* molecularDipoles;
    CudaArray* molecularQuadrupoles;
    CudaArray* labFrameDipoles;
    CudaArray* labFrameQuadrupoles;
    CudaArray* field;
    CudaArray* fieldPolar;
    CudaArray* inducedField;
    CudaArray* inducedFieldPolar;
    CudaArray* torque;
    CudaArray* dampingAndThole;
    CudaArray* inducedDipole;
    CudaArray* inducedDipolePolar;
    CudaArray* inducedDipoleErrors;
    CudaArray* polarizability;
    CudaArray* covalentFlags;
    CudaArray* polarizationGroupFlags;
    CudaArray* pmeGrid;
    CudaArray* pmeBsplineModuliX;
    CudaArray* pmeBsplineModuliY;
    CudaArray* pmeBsplineModuliZ;
    CudaArray* pmeIgrid;
    CudaArray* pmePhi;
    CudaArray* pmePhid;
    CudaArray* pmePhip;
    CudaArray* pmePhidp;
    CudaArray* pmeAtomRange;
    CudaArray* pmeAtomGridIndex;
    CudaArray* lastPositions;
    CudaSort* sort;
    cufftHandle fft;
    CUfunction computeMomentsKernel, recordInducedDipolesKernel, computeFixedFieldKernel, computeInducedFieldKernel, updateInducedFieldKernel, electrostaticsKernel, mapTorqueKernel;
    CUfunction pmeGridIndexKernel, pmeSpreadFixedMultipolesKernel, pmeSpreadInducedDipolesKernel, pmeFinishSpreadChargeKernel, pmeConvolutionKernel;
    CUfunction pmeFixedPotentialKernel, pmeInducedPotentialKernel, pmeFixedForceKernel, pmeInducedForceKernel, pmeRecordInducedFieldDipolesKernel, computePotentialKernel;
    CudaCalcMBPolGeneralizedKirkwoodForceKernel* gkKernel;
    static const int PmeOrder = 5;
};

/**
 * This kernel is invoked by MBPolMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolGeneralizedKirkwoodForceKernel : public CalcMBPolGeneralizedKirkwoodForceKernel {
public:
    CudaCalcMBPolGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolGeneralizedKirkwoodForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolGeneralizedKirkwoodForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Perform the computation of Born radii.
     */
    void computeBornRadii();
    /**
     * Perform the final parts of the force/energy computation.
     */
    void finishComputation(CudaArray& torque, CudaArray& labFrameDipoles, CudaArray& labFrameQuadrupoles, CudaArray& inducedDipole, CudaArray& inducedDipolePolar, CudaArray& dampingAndThole, CudaArray& covalentFlags, CudaArray& polarizationGroupFlags);
    CudaArray* getBornRadii() {
        return bornRadii;
    }
    CudaArray* getField() {
        return field;
    }
    CudaArray* getInducedField() {
        return inducedField;
    }
    CudaArray* getInducedFieldPolar() {
        return inducedFieldPolar;
    }
    CudaArray* getInducedDipoles() {
        return inducedDipoleS;
    }
    CudaArray* getInducedDipolesPolar() {
        return inducedDipolePolarS;
    }
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolGeneralizedKirkwoodForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolGeneralizedKirkwoodForce& force);
private:
    class ForceInfo;
    CudaContext& cu;
    const System& system;
    bool includeSurfaceArea, hasInitializedKernels;
    int computeBornSumThreads, gkForceThreads, chainRuleThreads, ediffThreads;
    std::map<std::string, std::string> defines;
    CudaArray* params;
    CudaArray* bornSum;
    CudaArray* bornRadii;
    CudaArray* bornForce;
    CudaArray* field;
    CudaArray* inducedField;
    CudaArray* inducedFieldPolar;
    CudaArray* inducedDipoleS;
    CudaArray* inducedDipolePolarS;
    CUfunction computeBornSumKernel, reduceBornSumKernel, surfaceAreaKernel, gkForceKernel, chainRuleKernel, ediffKernel;
};

/**
 * This kernel is invoked to calculate the vdw forces acting on the system and the energy of the system.
 */
class CudaCalcMBPolVdwForceKernel : public CalcMBPolVdwForceKernel {
public:
    CudaCalcMBPolVdwForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system);
    ~CudaCalcMBPolVdwForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolMultipoleForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolVdwForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolVdwForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolVdwForce& force);
private:
    class ForceInfo;
    CudaContext& cu;
    const System& system;
    bool hasInitializedNonBonded;
    double dispersionCoefficient;
    CudaArray* sigmaEpsilon;
    CudaArray* BondReductionAtoms;
    CudaArray* BondReductionFactors;
    CudaArray* tempPosq;
    CudaArray* tempForces;
    CudaNonBondedUtilities* nonBonded;
    CUfunction prepareKernel, spreadKernel;
};

} // namespace MBPolPlugin

#endif /*MBPOL_OPENMM_CUDAKERNELS_H*/
