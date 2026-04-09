#include <volt/analysis/pattern_classifier.h>

#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/analysis/crystal_topology_library.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/analysis/pattern_dxa_topology_provider.h>
#include <volt/analysis/pattern_matching_shared.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Volt {

constexpr double kPatternCanonicalVectorTolerance = 1e-4;

bool matchGenericLocalMatcherImpl(
    const CompiledPattern& pattern,
    const CompiledPatternLocalMatcher& matcher,
    int localMatcherIndex,
    const NearestNeighborFinder& neighFinder,
    int atomIndex,
    const std::vector<int>* atomTypes,
    bool allowCollapsedNeighborShell,
    bool* outRejectedByNeighborGap,
    PatternAtomMatch& outMatch
){
    if(outRejectedByNeighborGap){
        *outRejectedByNeighborGap = false;
    }

    NearestNeighborFinder::Query<MAX_NEIGHBORS> query(neighFinder);
    query.findNeighbors(neighFinder.particlePos(atomIndex));
    if(static_cast<int>(query.results().size()) < matcher.coordinationNumber){
        return false;
    }

    double localScale = 0.0;
    for(int neighborSlot = 0; neighborSlot < matcher.coordinationNumber; ++neighborSlot){
        localScale += std::sqrt(query.results()[neighborSlot].distanceSq) *
            matcher.scaleFactors[static_cast<std::size_t>(neighborSlot)];
    }
    if(localScale <= EPSILON){
        return false;
    }

    if(static_cast<int>(query.results().size()) > matcher.coordinationNumber){
        const double localGap =
            std::sqrt(query.results()[matcher.coordinationNumber].distanceSq) -
            std::sqrt(query.results()[matcher.coordinationNumber - 1].distanceSq);
        if(localGap <= matcher.cutoffGap * localScale){
            if(outRejectedByNeighborGap){
                *outRejectedByNeighborGap = true;
            }
            if(!allowCollapsedNeighborShell){
                return false;
            }
        }
    }

    NeighborBondArray runtimeNeighborArray;
    const double cutoffSquared = matcher.localCutoff * matcher.localCutoff * localScale * localScale;
    for(int ni1 = 0; ni1 < matcher.coordinationNumber; ++ni1){
        runtimeNeighborArray.setNeighborBond(ni1, ni1, false);
        for(int ni2 = ni1 + 1; ni2 < matcher.coordinationNumber; ++ni2){
            const bool bonded = (
                query.results()[ni1].delta - query.results()[ni2].delta
            ).squaredLength() < cutoffSquared;
            runtimeNeighborArray.setNeighborBond(ni1, ni2, bonded);
        }
    }

    std::vector<PatternCnaSignature> runtimeSignatures;
    runtimeSignatures.reserve(static_cast<std::size_t>(matcher.coordinationNumber));
    for(int neighborSlot = 0; neighborSlot < matcher.coordinationNumber; ++neighborSlot){
        runtimeSignatures.push_back(
            computePatternCnaSignature(
                runtimeNeighborArray,
                neighborSlot,
                matcher.coordinationNumber
            )
        );
    }

    std::vector<PatternCnaSignature> sortedRuntimeSignatures = runtimeSignatures;
    std::sort(sortedRuntimeSignatures.begin(), sortedRuntimeSignatures.end());
    if(sortedRuntimeSignatures != matcher.sortedCnaSignatures){
        return false;
    }

    std::vector<int> canonicalToRuntime;
    if(!assignGenericCnaPermutation(
        matcher,
        runtimeSignatures,
        runtimeNeighborArray,
        canonicalToRuntime
    )){
        return false;
    }

    std::vector<int> canonicalNeighborAtomIndices(static_cast<std::size_t>(matcher.coordinationNumber), -1);
    for(int canonicalSlot = 0; canonicalSlot < matcher.coordinationNumber; ++canonicalSlot){
        canonicalNeighborAtomIndices[static_cast<std::size_t>(canonicalSlot)] =
            query.results()[canonicalToRuntime[static_cast<std::size_t>(canonicalSlot)]].index;
    }
    if(!matchSpeciesWithSymmetry(
        matcher,
        canonicalNeighborAtomIndices,
        atomTypes,
        atomIndex,
        true
    )){
        return false;
    }

    outMatch = PatternAtomMatch{};
    outMatch.patternId = pattern.id;
    outMatch.structureType = pattern.structureType;
    outMatch.localMatcherIndex = localMatcherIndex;
    outMatch.localCutoff = matcher.localCutoff * localScale;
    outMatch.coordinationNumber = matcher.coordinationNumber;
    outMatch.allowedSymmetryMask = AnalysisSymmetryUtils::fullSymmetryMask(static_cast<int>(matcher.symmetries.size()));
    outMatch.symmetryPermutation = -1;
    for(int canonicalSlot = 0; canonicalSlot < matcher.coordinationNumber; ++canonicalSlot){
        outMatch.idealNeighborVectors[static_cast<std::size_t>(canonicalSlot)] =
            matcher.canonicalNeighborVectors[static_cast<std::size_t>(canonicalSlot)];
        outMatch.orderedNeighborIndices[static_cast<std::size_t>(canonicalSlot)] =
            canonicalNeighborAtomIndices[static_cast<std::size_t>(canonicalSlot)];
    }
    return true;
}

bool matchGenericLocalMatcher(
    const CompiledPattern& pattern,
    const CompiledPatternLocalMatcher& matcher,
    int localMatcherIndex,
    const NearestNeighborFinder& neighFinder,
    int atomIndex,
    const std::vector<int>* atomTypes,
    PatternAtomMatch& outMatch
){
    bool rejectedByNeighborGap = false;
    if(matchGenericLocalMatcherImpl(
        pattern,
        matcher,
        localMatcherIndex,
        neighFinder,
        atomIndex,
        atomTypes,
        false,
        &rejectedByNeighborGap,
        outMatch
    )){
        return true;
    }
    if(!rejectedByNeighborGap){
        return false;
    }
    return matchGenericLocalMatcherImpl(
        pattern,
        matcher,
        localMatcherIndex,
        neighFinder,
        atomIndex,
        atomTypes,
        true,
        nullptr,
        outMatch
    );
}

void canonicalizeSharedTopologyMatch(const CompiledPattern& pattern, PatternAtomMatch& match){
    if(match.structureType == StructureType::OTHER){
        return;
    }

    const SharedCrystalTopology* topology = sharedCrystalTopology(pattern.name);
    if(!topology || topology->coordinationNumber <= 0 || match.coordinationNumber != topology->coordinationNumber){
        return;
    }

    std::array<Vector3, MAX_NEIGHBORS> canonicalVectors;
    std::array<int, MAX_NEIGHBORS> canonicalNeighbors;
    canonicalVectors.fill(Vector3::Zero());
    canonicalNeighbors.fill(-1);
    std::array<bool, MAX_NEIGHBORS> consumed{};
    consumed.fill(false);

    for(int canonicalSlot = 0; canonicalSlot < topology->coordinationNumber; ++canonicalSlot){
        const Vector3& expectedVector = topology->neighborVectors[static_cast<std::size_t>(canonicalSlot)];
        int matchedSlot = -1;
        for(int slot = 0; slot < match.coordinationNumber; ++slot){
            if(consumed[static_cast<std::size_t>(slot)]){
                continue;
            }
            if((match.idealNeighborVectors[static_cast<std::size_t>(slot)] - expectedVector)
                .isZero(kPatternCanonicalVectorTolerance)){
                matchedSlot = slot;
                break;
            }
        }
        if(matchedSlot < 0){
            return;
        }
        consumed[static_cast<std::size_t>(matchedSlot)] = true;
        canonicalVectors[static_cast<std::size_t>(canonicalSlot)] = expectedVector;
        canonicalNeighbors[static_cast<std::size_t>(canonicalSlot)] =
            match.orderedNeighborIndices[static_cast<std::size_t>(matchedSlot)];
    }

    match.idealNeighborVectors = canonicalVectors;
    match.orderedNeighborIndices = canonicalNeighbors;
    match.allowedSymmetryMask = AnalysisSymmetryUtils::fullSymmetryMask(
        static_cast<int>(topology->symmetries.size())
    );
    const int identityLikeSymmetry = findClosestSharedCrystalSymmetryPermutation(
        *topology,
        Matrix3::Identity()
    );
    match.symmetryPermutation = identityLikeSymmetry >= 0 ? identityLikeSymmetry : 0;
}

bool matchPatternAtom(
    const CompiledPattern& pattern,
    const NearestNeighborFinder& neighFinder,
    int atomIndex,
    const std::vector<int>* atomTypes,
    PatternAtomMatch& outMatch
){
    for(std::size_t matcherIndex = 0; matcherIndex < pattern.localMatchers.size(); ++matcherIndex){
        if(matchGenericLocalMatcher(
            pattern,
            pattern.localMatchers[matcherIndex],
            static_cast<int>(matcherIndex),
            neighFinder,
            atomIndex,
            atomTypes,
            outMatch
        )){
            canonicalizeSharedTopologyMatch(pattern, outMatch);
            return true;
        }
    }
    return false;
}

int patternSelectionPriority(const CompiledPattern& pattern){
    if(pattern.requiresAtomTypes){
        return 1;
    }
    if(!pattern.localMatchers.empty()){
        return 2;
    }
    return 3;
}

PatternClassifier::PatternClassifier(std::shared_ptr<const PatternCatalog> catalog)
    : _catalog(std::move(catalog)),
      _atomDxaStates(std::make_shared<std::vector<PatternDxaAtomState>>()){}

void PatternClassifier::setAtomTypes(const std::vector<int>* atomTypes){
    _atomTypes = atomTypes;
}

void PatternClassifier::setSelectedPatternIds(std::vector<int> selectedPatternIds){
    _configuredSelectedPatternIds = std::move(selectedPatternIds);
}

void PatternClassifier::classify(StructureAnalysis& analysis){
    if(!_catalog){
        throw std::runtime_error("Pattern classifier catalog is not configured");
    }

    _selectedPatternIds = _configuredSelectedPatternIds.empty()
        ? _catalog->defaultSelection()
        : _configuredSelectedPatternIds;
    for(int patternId : _selectedPatternIds){
        if(patternId < 0 || patternId >= static_cast<int>(_catalog->patterns().size())){
            throw std::runtime_error("Configured pattern id out of range.");
        }
    }
    std::stable_sort(_selectedPatternIds.begin(), _selectedPatternIds.end(), [&](int lhs, int rhs){
        return patternSelectionPriority(_catalog->patternById(lhs)) < patternSelectionPriority(_catalog->patternById(rhs));
    });
    if(_selectedPatternIds.empty()){
        throw std::runtime_error("No supported lattice patterns available for pattern matching.");
    }

    StructureContext& context = analysis.context();
    analysis.setCrystalInfoProvider(std::make_shared<PatternDxaTopologyProvider>(_catalog, _selectedPatternIds));

    std::fill(
        context.structureTypes->dataInt(),
        context.structureTypes->dataInt() + context.structureTypes->size(),
        static_cast<int>(StructureType::OTHER)
    );
    std::fill(
        context.neighborCounts->dataInt(),
        context.neighborCounts->dataInt() + context.neighborCounts->size(),
        0
    );
    if(context.atomAllowedSymmetryMasks){
        std::fill(
            context.atomAllowedSymmetryMasks->dataInt64(),
            context.atomAllowedSymmetryMasks->dataInt64() + context.atomAllowedSymmetryMasks->size(),
            0
        );
    }
    if(context.atomSymmetryPermutations){
        std::fill(
            context.atomSymmetryPermutations->dataInt(),
            context.atomSymmetryPermutations->dataInt() + context.atomSymmetryPermutations->size(),
            -1
        );
    }

    _atomPatternIds.assign(context.atomCount(), -1);
    _atomDxaStates = std::make_shared<std::vector<PatternDxaAtomState>>(context.atomCount());

    NearestNeighborFinder neighborFinder(MAX_NEIGHBORS);
    if(!neighborFinder.prepare(context.positions, context.simCell, context.particleSelection)){
        throw std::runtime_error("Error preparing nearest-neighbor finder for pattern matching.");
    }

    const std::size_t atomCount = context.atomCount();
    std::vector<int> localCounts(atomCount, 0);
    std::vector<std::array<int, MAX_NEIGHBORS>> orderedNeighborIndices(atomCount);
    for(auto& neighbors : orderedNeighborIndices){
        neighbors.fill(-1);
    }
    auto storeMatch = [&](std::size_t atomIndex, const PatternAtomMatch& match, double& localMaxDistance){
        context.structureTypes->setInt(atomIndex, match.structureType);
        _atomPatternIds[atomIndex] = match.patternId;
        auto& dxaState = (*_atomDxaStates)[atomIndex];
        dxaState.patternId = match.patternId;
        dxaState.structureType = match.structureType;
        dxaState.localMatcherIndex = match.localMatcherIndex;
        dxaState.coordinationNumber = match.coordinationNumber;
        dxaState.allowedSymmetryMask = match.allowedSymmetryMask;
        dxaState.symmetryPermutation = match.symmetryPermutation;
        dxaState.neighborAtomIndices = match.orderedNeighborIndices;
        dxaState.idealNeighborVectors = match.idealNeighborVectors;
        localCounts[atomIndex] = match.coordinationNumber;
        orderedNeighborIndices[atomIndex] = match.orderedNeighborIndices;
        localMaxDistance = std::max(localMaxDistance, match.localCutoff);
    };

    context.maximumNeighborDistance = tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, atomCount),
        0.0,
        [&](const tbb::blocked_range<std::size_t>& range, double localMaxDistance) -> double {
            for(std::size_t atomIndex = range.begin(); atomIndex != range.end(); ++atomIndex){
                for(int patternId : _selectedPatternIds){
                    const CompiledPattern& pattern = _catalog->patternById(patternId);
                    PatternAtomMatch match;
                    if(!matchPatternAtom(pattern, neighborFinder, static_cast<int>(atomIndex), _atomTypes, match)){
                        continue;
                    }

                    storeMatch(atomIndex, match, localMaxDistance);
                    break;
                }
            }
            return localMaxDistance;
        },
        [](double lhs, double rhs) -> double {
            return std::max(lhs, rhs);
        }
    );
    auto* offsets = context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        offsets[atomIndex + 1] = offsets[atomIndex] + localCounts[atomIndex];
    }

    const std::size_t totalNeighbors = static_cast<std::size_t>(offsets[atomCount]);
    context.neighborIndices = std::make_shared<ParticleProperty>(totalNeighbors, DataType::Int, 1, 0, false);

    int* neighborIndices = context.neighborIndices->dataInt();
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, atomCount), [&](const tbb::blocked_range<std::size_t>& range){
        for(std::size_t atomIndex = range.begin(); atomIndex != range.end(); ++atomIndex){
            const int count = localCounts[atomIndex];
            context.atomAllowedSymmetryMasks->setInt64(
                atomIndex,
                count > 0 ? static_cast<std::int64_t>((*_atomDxaStates)[atomIndex].allowedSymmetryMask) : 0
            );
            context.atomSymmetryPermutations->setInt(
                atomIndex,
                count > 0 ? (*_atomDxaStates)[atomIndex].symmetryPermutation : -1
            );
            if(count <= 0){
                continue;
            }

            const int start = offsets[atomIndex];
            for(int neighborSlot = 0; neighborSlot < count; ++neighborSlot){
                neighborIndices[start + neighborSlot] = orderedNeighborIndices[atomIndex][static_cast<std::size_t>(neighborSlot)];
            }
            context.neighborCounts->setInt(atomIndex, count);
        }
    });
}

void PatternClassifier::configureDxaClustering(StructureAnalysis& analysis) const{
    analysis.setCrystalInfoProvider(std::make_shared<PatternDxaTopologyProvider>(_catalog, _selectedPatternIds));
}

}
