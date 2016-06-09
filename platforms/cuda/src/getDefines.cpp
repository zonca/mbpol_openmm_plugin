#include <sstream>
#include <iostream>
using namespace std;
int main() {
    bool hasTypeFilters = false;
    int particlesPerSet = 3;
    bool centralParticleMode = false;
    int nonbondedMethod = 1, NoCutoff = 1;
    // Create other replacements that depend on the number of particles per set.

       stringstream numCombinations, atomsForCombination, isValidCombination, permute, loadData, verifyCutoff, verifyExclusions;
       if (hasTypeFilters) {
           permute<<"int particleSet[] = {";
           for (int i = 0; i < particlesPerSet; i++) {
               permute<<"p"<<(i+1);
               if (i < particlesPerSet-1)
                   permute<<", ";
           }
           permute<<"};\n";
       }
       for (int i = 0; i < particlesPerSet; i++) {
           if (hasTypeFilters)
               permute<<"int atom"<<(i+1)<<" = particleSet[particleOrder["<<particlesPerSet<<"*order+"<<i<<"]];\n";
           else
               permute<<"int atom"<<(i+1)<<" = p"<<(i+1)<<";\n";
           loadData<<"real3 pos"<<(i+1)<<" = trim(posq[atom"<<(i+1)<<"]);\n";
//           for (int j = 0; j < (int) params->getBuffers().size(); j++)
//               loadData<<params->getBuffers()[j].getType()<<" params"<<(j+1)<<(i+1)<<" = global_params"<<(j+1)<<"[atom"<<(i+1)<<"];\n";
       }
       if (centralParticleMode) {
           for (int i = 1; i < particlesPerSet; i++) {
               if (i > 1)
                   isValidCombination<<" && p"<<(i+1)<<">p"<<i<<" && ";
               isValidCombination<<"p"<<(i+1)<<"!=p1";
           }
       }
       else {
           for (int i = 2; i < particlesPerSet; i++) {
               if (i > 2)
                   isValidCombination<<" && ";
               isValidCombination<<"a"<<(i+1)<<">a"<<i;
           }
       }
       atomsForCombination<<"int tempIndex = index;\n";
       for (int i = 1; i < particlesPerSet; i++) {
           if (i > 1)
               numCombinations<<"*";
           numCombinations<<"numNeighbors";
           if (centralParticleMode)
               atomsForCombination<<"int a"<<(i+1)<<" = tempIndex%numNeighbors;\n";
           else
               atomsForCombination<<"int a"<<(i+1)<<" = 1+tempIndex%numNeighbors;\n";
           if (i < particlesPerSet-1)
               atomsForCombination<<"tempIndex /= numNeighbors;\n";
       }
       if (particlesPerSet > 2) {
           if (centralParticleMode)
               atomsForCombination<<"a2 = (a3%2 == 0 ? a2 : numNeighbors-a2-1);\n";
           else
               atomsForCombination<<"a2 = (a3%2 == 0 ? a2 : numNeighbors-a2+1);\n";
       }
       for (int i = 1; i < particlesPerSet; i++) {
           if (nonbondedMethod == NoCutoff) {
               if (centralParticleMode)
                   atomsForCombination<<"int p"<<(i+1)<<" = a"<<(i+1)<<";\n";
               else
                   atomsForCombination<<"int p"<<(i+1)<<" = p1+a"<<(i+1)<<";\n";
           }
           else {
               if (centralParticleMode)
                   atomsForCombination<<"int p"<<(i+1)<<" = neighbors[firstNeighbor+a"<<(i+1)<<"];\n";
               else
                   atomsForCombination<<"int p"<<(i+1)<<" = neighbors[firstNeighbor-1+a"<<(i+1)<<"];\n";
           }
       }
       if (nonbondedMethod != NoCutoff) {
           for (int i = 1; i < particlesPerSet; i++)
               verifyCutoff<<"real3 pos"<<(i+1)<<" = trim(posq[p"<<(i+1)<<"]);\n";
           if (!centralParticleMode) {
               for (int i = 1; i < particlesPerSet; i++) {
                   for (int j = i+1; j < particlesPerSet; j++)
                       verifyCutoff<<"includeInteraction &= (delta(pos"<<(i+1)<<", pos"<<(j+1)<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ).w < CUTOFF_SQUARED);\n";
               }
           }
       }
//       if (force.getNumExclusions() > 0) {
//           int startCheckFrom = (nonbondedMethod == NoCutoff ? 0 : 1);
//           for (int i = startCheckFrom; i < particlesPerSet; i++)
//               for (int j = i+1; j < particlesPerSet; j++)
//                   verifyExclusions<<"includeInteraction &= !isInteractionExcluded(p"<<(i+1)<<", p"<<(j+1)<<", exclusions, exclusionStartIndex);\n";
//       }
//       string computeTypeIndex = "particleTypes[p"+cu.intToString(particlesPerSet)+"]";
//       for (int i = particlesPerSet-2; i >= 0; i--)
//           computeTypeIndex = "particleTypes[p"+cu.intToString(i+1)+"]+"+cu.intToString(numTypes)+"*("+computeTypeIndex+")";

    std::cout << atomsForCombination.str() << std::endl;
}
