/* -------------------------------------------------------------------------- *
 *                              OpenMMMBPol                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
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

#include "MBPolCpuKernelFactory.h"
#include "MBPolCpuKernels.h"
#include "openmm/cpu/CpuPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace  OpenMM;
using namespace MBPolPlugin;

// need to add registerKernelFactories from https://github.com/peastman/openmmexampleplugin/blob/master/platforms/reference/src/ReferenceExampleKernelFactory.cpp

//extern "C" OPENMM_EXPORT void registerKernelFactories() {
//    for (int i = 0; i < Platform::getNumPlatforms(); i++) {
//        Platform& platform = Platform::getPlatform(i);
//        if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
//            ReferenceExampleKernelFactory* factory = new ReferenceExampleKernelFactory();
//            platform.registerKernelFactory(CalcExampleForceKernel::Name(), factory);
//        }
//    }
//}

extern "C" void initMBPolCpuKernels() {
    for( int ii = 0; ii < Platform::getNumPlatforms(); ii++ ){
        Platform& platform = Platform::getPlatform(ii);
        if( platform.getName() == "Cpu" ){

             MBPolCpuKernelFactory* factory = new MBPolCpuKernelFactory();

             //platform.registerKernelFactory(CalcMBPolOneBodyForceKernel::Name(),           factory);
             //platform.registerKernelFactory(CalcMBPolTwoBodyForceKernel::Name(),                   factory);
             //platform.registerKernelFactory(CalcMBPolThreeBodyForceKernel::Name(),                   factory);
             //platform.registerKernelFactory(CalcMBPolElectrostaticsForceKernel::Name(),             factory);
        }
    }
}

KernelImpl* MBPolCpuKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CpuPlatform::PlatformData& data = *static_cast<CpuPlatform::PlatformData*>(context.getPlatformData());

    // create MBPolCpuData object if contextToMBPolDataMap does not contain
    // key equal to current context
    //if (name == CalcMBPolOneBodyForceKernel::Name())
    //    return new ReferenceCalcMBPolOneBodyForceKernel(name, platform, context.getSystem());

    //if (name == CalcMBPolTwoBodyForceKernel::Name())
    //    return new ReferenceCalcMBPolTwoBodyForceKernel(name, platform, context.getSystem());

    //if (name == CalcMBPolThreeBodyForceKernel::Name())
    //        return new ReferenceCalcMBPolThreeBodyForceKernel(name, platform, context.getSystem());

    //if (name == CalcMBPolElectrostaticsForceKernel::Name())
    //    return new ReferenceCalcMBPolElectrostaticsForceKernel(name, platform, context.getSystem());

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
