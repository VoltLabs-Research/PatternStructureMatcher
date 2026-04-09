#include <volt/analysis/pattern_matching_shared.h>

#include <volt/analysis/cna_classifier.h>

namespace Volt {

static bool assignGenericCnaPermutationRecursive(
    const CompiledPatternLocalMatcher& matcher,
    const std::vector<PatternCnaSignature>& runtimeSignatures,
    const NeighborBondArray& runtimeNeighborArray,
    std::vector<int>& canonicalToRuntime,
    std::vector<unsigned char>& usedRuntimeSlots,
    int canonicalSlot
){
    if(canonicalSlot == matcher.coordinationNumber){
        return true;
    }

    const PatternCnaSignature& expectedSignature = matcher.cnaSignatures[static_cast<std::size_t>(canonicalSlot)];
    for(int runtimeSlot = 0; runtimeSlot < matcher.coordinationNumber; ++runtimeSlot){
        if(usedRuntimeSlots[static_cast<std::size_t>(runtimeSlot)]){
            continue;
        }
        if(!(runtimeSignatures[static_cast<std::size_t>(runtimeSlot)] == expectedSignature)){
            continue;
        }

        bool bondsMatch = true;
        for(int previousCanonicalSlot = 0; previousCanonicalSlot < canonicalSlot; ++previousCanonicalSlot){
            const int previousRuntimeSlot = canonicalToRuntime[static_cast<std::size_t>(previousCanonicalSlot)];
            if(previousRuntimeSlot < 0){
                continue;
            }

            const bool runtimeBond = runtimeNeighborArray.neighborBond(runtimeSlot, previousRuntimeSlot);
            const bool expectedBond =
                (matcher.neighborBondRows[static_cast<std::size_t>(canonicalSlot)] & (1u << previousCanonicalSlot)) != 0;
            if(runtimeBond != expectedBond){
                bondsMatch = false;
                break;
            }
        }
        if(!bondsMatch){
            continue;
        }

        canonicalToRuntime[static_cast<std::size_t>(canonicalSlot)] = runtimeSlot;
        usedRuntimeSlots[static_cast<std::size_t>(runtimeSlot)] = 1;
        if(assignGenericCnaPermutationRecursive(
            matcher,
            runtimeSignatures,
            runtimeNeighborArray,
            canonicalToRuntime,
            usedRuntimeSlots,
            canonicalSlot + 1
        )){
            return true;
        }
        usedRuntimeSlots[static_cast<std::size_t>(runtimeSlot)] = 0;
        canonicalToRuntime[static_cast<std::size_t>(canonicalSlot)] = -1;
    }

    return false;
}

static bool atomTypesAvailable(const std::vector<int>* atomTypes, std::size_t atomCount){
    return atomTypes != nullptr && atomTypes->size() >= atomCount;
}

PatternCnaSignature computePatternCnaSignature(
    const NeighborBondArray& neighborArray,
    int neighborIndex,
    int coordinationNumber
){
    unsigned int commonNeighbors = 0;
    const int numCommonNeighbors = CommonNeighborAnalysis::findCommonNeighbors(
        neighborArray,
        neighborIndex,
        commonNeighbors,
        coordinationNumber
    );

    CommonNeighborAnalysis::CNAPairBond neighborBonds[MAX_NEIGHBORS * MAX_NEIGHBORS];
    const int numBonds = CommonNeighborAnalysis::findNeighborBonds(
        neighborArray,
        commonNeighbors,
        coordinationNumber,
        neighborBonds
    );

    PatternCnaSignature signature;
    signature.numCommonNeighbors = numCommonNeighbors;
    signature.numBonds = numBonds;
    signature.maxChainLength = CommonNeighborAnalysis::calcMaxChainLength(neighborBonds, numBonds);
    return signature;
}

bool matchSpeciesWithSymmetry(
    const CompiledPatternLocalMatcher& matcher,
    const std::vector<int>& canonicalNeighborAtomIndices,
    const std::vector<int>* atomTypes,
    int centerAtomIndex,
    bool allowEmptySymmetryFallback
){
    if(!matcher.requiresSpecies){
        return true;
    }
    if(!atomTypesAvailable(atomTypes, static_cast<std::size_t>(centerAtomIndex) + 1)){
        return false;
    }
    if((*atomTypes)[static_cast<std::size_t>(centerAtomIndex)] != matcher.centerSpecies){
        return false;
    }

    auto symmetryMatches = [&](const PatternSymmetryPermutation& symmetry){
        for(int canonicalSlot = 0; canonicalSlot < matcher.coordinationNumber; ++canonicalSlot){
            const int neighborAtomIndex = canonicalNeighborAtomIndices[static_cast<std::size_t>(canonicalSlot)];
            if(neighborAtomIndex < 0 || !atomTypesAvailable(atomTypes, static_cast<std::size_t>(neighborAtomIndex) + 1)){
                return false;
            }

            const int mappedSlot = symmetry.permutation[static_cast<std::size_t>(canonicalSlot)];
            if(mappedSlot < 0 || mappedSlot >= static_cast<int>(matcher.neighborSpeciesByCanonicalSlot.size())){
                return false;
            }
            if((*atomTypes)[static_cast<std::size_t>(neighborAtomIndex)] !=
               matcher.neighborSpeciesByCanonicalSlot[static_cast<std::size_t>(mappedSlot)]){
                return false;
            }
        }
        return true;
    };

    for(const auto& symmetry : matcher.symmetries){
        if(symmetryMatches(symmetry)){
            return true;
        }
    }

    return allowEmptySymmetryFallback && matcher.symmetries.empty() && matcher.neighborSpeciesByCanonicalSlot.empty();
}

bool assignGenericCnaPermutation(
    const CompiledPatternLocalMatcher& matcher,
    const std::vector<PatternCnaSignature>& runtimeSignatures,
    const NeighborBondArray& runtimeNeighborArray,
    std::vector<int>& canonicalToRuntime
){
    canonicalToRuntime.assign(static_cast<std::size_t>(matcher.coordinationNumber), -1);
    std::vector<unsigned char> usedRuntimeSlots(static_cast<std::size_t>(matcher.coordinationNumber), 0);
    return assignGenericCnaPermutationRecursive(
        matcher,
        runtimeSignatures,
        runtimeNeighborArray,
        canonicalToRuntime,
        usedRuntimeSlots,
        0
    );
}

bool equivalentLocalMatchers(
    const CompiledPatternLocalMatcher& lhs,
    const CompiledPatternLocalMatcher& rhs
){
    return lhs.kind == rhs.kind &&
        lhs.coordinationNumber == rhs.coordinationNumber &&
        lhs.centerSpecies == rhs.centerSpecies &&
        lhs.requiresSpecies == rhs.requiresSpecies &&
        lhs.neighborSpeciesByCanonicalSlot == rhs.neighborSpeciesByCanonicalSlot &&
        lhs.sortedCnaSignatures == rhs.sortedCnaSignatures;
}

}
