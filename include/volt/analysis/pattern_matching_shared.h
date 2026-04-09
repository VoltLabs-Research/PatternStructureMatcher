#pragma once

#include <volt/analysis/pattern_catalog.h>
#include <volt/structures/neighbor_bond_array.h>

#include <vector>

namespace Volt {

PatternCnaSignature computePatternCnaSignature(
    const NeighborBondArray& neighborArray,
    int neighborIndex,
    int coordinationNumber
);

bool matchSpeciesWithSymmetry(
    const CompiledPatternLocalMatcher& matcher,
    const std::vector<int>& canonicalNeighborAtomIndices,
    const std::vector<int>* atomTypes,
    int centerAtomIndex,
    bool allowEmptySymmetryFallback
);

bool assignGenericCnaPermutation(
    const CompiledPatternLocalMatcher& matcher,
    const std::vector<PatternCnaSignature>& runtimeSignatures,
    const NeighborBondArray& runtimeNeighborArray,
    std::vector<int>& canonicalToRuntime
);

bool equivalentLocalMatchers(
    const CompiledPatternLocalMatcher& lhs,
    const CompiledPatternLocalMatcher& rhs
);

}
