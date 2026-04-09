#include <volt/analysis/pattern_classifier.h>

#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/analysis/pattern_dxa_topology_provider.h>
#include <volt/analysis/pattern_matching_shared.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace Volt {

namespace{

bool orthonormalizeMatrix(const Matrix3& input, Matrix3& output){
    Vector3 c0 = input.column(0);
    Vector3 c1 = input.column(1);
    Vector3 c2 = input.column(2);

    const double l0 = c0.length();
    if(l0 <= EPSILON){
        return false;
    }
    c0 /= l0;

    c1 -= c0 * c0.dot(c1);
    const double l1 = c1.length();
    if(l1 <= EPSILON){
        return false;
    }
    c1 /= l1;

    c2 -= c0 * c0.dot(c2);
    c2 -= c1 * c1.dot(c2);
    const double l2 = c2.length();
    if(l2 <= EPSILON){
        return false;
    }
    c2 /= l2;

    output = Matrix3(c0, c1, c2);
    if(output.determinant() < 0.0){
        output.column(2) = -output.column(2);
    }
    return true;
}

bool lexicographicallyGreaterMatrix(const Matrix3& lhs, const Matrix3& rhs){
    constexpr double epsilon = 1e-8;
    for(int row = 0; row < 3; ++row){
        for(int column = 0; column < 3; ++column){
            const double delta = lhs(row, column) - rhs(row, column);
            if(std::abs(delta) <= epsilon){
                continue;
            }
            return delta > 0.0;
        }
    }
    return false;
}

void permuteCanonicalToRuntimeBySymmetry(
    const std::vector<int>& canonicalToRuntime,
    const PatternSymmetryPermutation& symmetry,
    int coordinationNumber,
    std::vector<int>& outCanonicalToRuntime
){
    outCanonicalToRuntime.assign(static_cast<std::size_t>(coordinationNumber), -1);

    std::array<int, MAX_NEIGHBORS> inversePermutation{};
    inversePermutation.fill(-1);
    for(int slot = 0; slot < coordinationNumber; ++slot){
        const int mappedSlot = symmetry.permutation[static_cast<std::size_t>(slot)];
        if(mappedSlot >= 0 && mappedSlot < coordinationNumber){
            inversePermutation[static_cast<std::size_t>(mappedSlot)] = slot;
        }
    }

    for(int canonicalSlot = 0; canonicalSlot < coordinationNumber; ++canonicalSlot){
        const int inverseSlot = inversePermutation[static_cast<std::size_t>(canonicalSlot)];
        if(inverseSlot < 0 || inverseSlot >= coordinationNumber){
            continue;
        }
        outCanonicalToRuntime[static_cast<std::size_t>(canonicalSlot)] =
            canonicalToRuntime[static_cast<std::size_t>(inverseSlot)];
    }
}

double orientationAssignmentError(
    const CompiledPatternLocalMatcher& matcher,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& query,
    const std::vector<int>& canonicalToRuntime,
    const Matrix3& orientation
){
    double error = 0.0;
    for(int canonicalSlot = 0; canonicalSlot < matcher.coordinationNumber; ++canonicalSlot){
        const int runtimeSlot = canonicalToRuntime[static_cast<std::size_t>(canonicalSlot)];
        if(runtimeSlot < 0 || runtimeSlot >= static_cast<int>(query.results().size())){
            return std::numeric_limits<double>::infinity();
        }

        Vector3 expected = orientation * matcher.canonicalNeighborVectors[static_cast<std::size_t>(canonicalSlot)];
        Vector3 actual = query.results()[runtimeSlot].delta;
        const double expectedLength = expected.length();
        const double actualLength = actual.length();
        if(expectedLength <= EPSILON || actualLength <= EPSILON){
            return std::numeric_limits<double>::infinity();
        }

        expected /= expectedLength;
        actual /= actualLength;
        const double dot = std::clamp(expected.dot(actual), -1.0, 1.0);
        error += 1.0 - dot;
    }
    return error;
}

bool computeMatchOrientation(
    const CompiledPatternLocalMatcher& matcher,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& query,
    const std::vector<int>& canonicalToRuntime,
    Matrix3& outOrientation
){
    if(static_cast<int>(canonicalToRuntime.size()) < matcher.coordinationNumber || matcher.coordinationNumber < 3){
        return false;
    }

    Matrix3 orientationV = Matrix3::Zero();
    Matrix3 orientationW = Matrix3::Zero();

    for(int canonicalSlot = 0; canonicalSlot < matcher.coordinationNumber; ++canonicalSlot){
        const int runtimeSlot = canonicalToRuntime[static_cast<std::size_t>(canonicalSlot)];
        if(runtimeSlot < 0 || runtimeSlot >= static_cast<int>(query.results().size())){
            return false;
        }

        const Vector3& idealVector = matcher.canonicalNeighborVectors[static_cast<std::size_t>(canonicalSlot)];
        const Vector3& spatialVector = query.results()[runtimeSlot].delta;
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 3; ++j){
                orientationV(i, j) += idealVector[j] * idealVector[i];
                orientationW(i, j) += idealVector[j] * spatialVector[i];
            }
        }
    }

    Matrix3 orientationVInverse;
    if(!orientationV.inverse(orientationVInverse)){
        return false;
    }

    const Matrix3 rawOrientation = Matrix3(orientationW * orientationVInverse);
    return orthonormalizeMatrix(rawOrientation, outOrientation);
}

double matrixTrace(const Matrix3& matrix){
    return matrix(0, 0) + matrix(1, 1) + matrix(2, 2);
}

bool reduceOrientationToFundamentalZone(
    const CompiledPatternLocalMatcher& matcher,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& query,
    std::vector<int>& inOutCanonicalToRuntime,
    Matrix3& inOutOrientation
){
    if(matcher.symmetries.empty()){
        return false;
    }

    std::vector<int> bestCanonicalToRuntime = inOutCanonicalToRuntime;
    Matrix3 bestOrientation = inOutOrientation;
    double bestTrace = matrixTrace(inOutOrientation);
    double bestError = orientationAssignmentError(
        matcher,
        query,
        inOutCanonicalToRuntime,
        inOutOrientation
    );
    bool changed = false;

    for(const PatternSymmetryPermutation& symmetry : matcher.symmetries){
        std::vector<int> candidateCanonicalToRuntime;
        permuteCanonicalToRuntimeBySymmetry(
            inOutCanonicalToRuntime,
            symmetry,
            matcher.coordinationNumber,
            candidateCanonicalToRuntime
        );

        Matrix3 candidateOrientation;
        if(!computeMatchOrientation(
            matcher,
            query,
            candidateCanonicalToRuntime,
            candidateOrientation
        )){
            continue;
        }

        const double candidateTrace = matrixTrace(candidateOrientation);
        const double candidateError = orientationAssignmentError(
            matcher,
            query,
            candidateCanonicalToRuntime,
            candidateOrientation
        );

        const bool preferCandidate =
            candidateTrace > bestTrace + 1e-8 ||
            (std::abs(candidateTrace - bestTrace) <= 1e-8 && candidateError < bestError - 1e-8) ||
            (std::abs(candidateTrace - bestTrace) <= 1e-8 &&
             std::abs(candidateError - bestError) <= 1e-8 &&
             lexicographicallyGreaterMatrix(candidateOrientation, bestOrientation));

        if(!preferCandidate){
            continue;
        }

        bestCanonicalToRuntime = std::move(candidateCanonicalToRuntime);
        bestOrientation = candidateOrientation;
        bestTrace = candidateTrace;
        bestError = candidateError;
        changed = true;
    }

    if(changed){
        inOutCanonicalToRuntime = std::move(bestCanonicalToRuntime);
        inOutOrientation = bestOrientation;
    }
    return changed;
}

bool assignBestOrientationFallbackPermutation(
    const CompiledPatternLocalMatcher& matcher,
    const std::vector<PatternCnaSignature>& runtimeSignatures,
    const NeighborBondArray& runtimeNeighborArray,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& query,
    std::vector<int>& outCanonicalToRuntime,
    Matrix3& outOrientation
){
    std::vector<int> currentCanonicalToRuntime(static_cast<std::size_t>(matcher.coordinationNumber), -1);
    std::vector<int> bestCanonicalToRuntime = currentCanonicalToRuntime;
    std::vector<unsigned char> usedRuntimeSlots(static_cast<std::size_t>(matcher.coordinationNumber), 0);
    Matrix3 bestOrientation = Matrix3::Identity();
    double bestError = std::numeric_limits<double>::infinity();
    bool foundBest = false;

    std::function<void(int)> visit = [&](int canonicalSlot){
        if(canonicalSlot == matcher.coordinationNumber){
            Matrix3 candidateOrientation;
            if(!computeMatchOrientation(
                matcher,
                query,
                currentCanonicalToRuntime,
                candidateOrientation
            )){
                return;
            }

            const double candidateError = orientationAssignmentError(
                matcher,
                query,
                currentCanonicalToRuntime,
                candidateOrientation
            );
            const bool preferCandidate =
                !foundBest ||
                candidateError < bestError - 1e-8 ||
                (std::abs(candidateError - bestError) <= 1e-8 &&
                 lexicographicallyGreaterMatrix(candidateOrientation, bestOrientation));
            if(!preferCandidate){
                return;
            }

            foundBest = true;
            bestError = candidateError;
            bestOrientation = candidateOrientation;
            bestCanonicalToRuntime = currentCanonicalToRuntime;
            return;
        }

        const PatternCnaSignature& expectedSignature =
            matcher.cnaSignatures[static_cast<std::size_t>(canonicalSlot)];

        for(int runtimeSlot = 0; runtimeSlot < matcher.coordinationNumber; ++runtimeSlot){
            if(usedRuntimeSlots[static_cast<std::size_t>(runtimeSlot)]){
                continue;
            }
            if(!(runtimeSignatures[static_cast<std::size_t>(runtimeSlot)] == expectedSignature)){
                continue;
            }

            bool bondsMatch = true;
            for(int previousCanonicalSlot = 0; previousCanonicalSlot < canonicalSlot; ++previousCanonicalSlot){
                const int previousRuntimeSlot =
                    currentCanonicalToRuntime[static_cast<std::size_t>(previousCanonicalSlot)];
                if(previousRuntimeSlot < 0){
                    continue;
                }

                const bool runtimeBond = runtimeNeighborArray.neighborBond(runtimeSlot, previousRuntimeSlot);
                const bool expectedBond =
                    (matcher.neighborBondRows[static_cast<std::size_t>(canonicalSlot)] &
                     (1u << previousCanonicalSlot)) != 0;
                if(runtimeBond != expectedBond){
                    bondsMatch = false;
                    break;
                }
            }
            if(!bondsMatch){
                continue;
            }

            currentCanonicalToRuntime[static_cast<std::size_t>(canonicalSlot)] = runtimeSlot;
            usedRuntimeSlots[static_cast<std::size_t>(runtimeSlot)] = 1;
            visit(canonicalSlot + 1);
            usedRuntimeSlots[static_cast<std::size_t>(runtimeSlot)] = 0;
            currentCanonicalToRuntime[static_cast<std::size_t>(canonicalSlot)] = -1;
        }
    };

    visit(0);
    if(!foundBest){
        return false;
    }

    outCanonicalToRuntime = std::move(bestCanonicalToRuntime);
    outOrientation = bestOrientation;
    return true;
}

}

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

    if(matcher.referenceNeighborCount <= 0 ||
       matcher.referenceNeighborOffset < 0 ||
       matcher.referenceNeighborOffset + matcher.referenceNeighborCount > matcher.coordinationNumber ||
       matcher.cutoffMultiplier <= EPSILON){
        return false;
    }

    double localScale = 0.0;
    const int referenceBegin = matcher.referenceNeighborOffset;
    const int referenceEnd = referenceBegin + matcher.referenceNeighborCount;
    for(int neighborSlot = referenceBegin; neighborSlot < referenceEnd; ++neighborSlot){
        localScale += std::sqrt(query.results()[neighborSlot].distanceSq);
    }
    localScale /= static_cast<double>(matcher.referenceNeighborCount);
    if(localScale <= EPSILON){
        return false;
    }

    const double localCutoff = matcher.cutoffMultiplier * localScale;
    const double cutoffSquared = localCutoff * localCutoff;

    if(matcher.extraNeighborRejectIndex >= 0 &&
       static_cast<int>(query.results().size()) > matcher.extraNeighborRejectIndex &&
       query.results()[matcher.extraNeighborRejectIndex].distanceSq <= cutoffSquared){
        if(outRejectedByNeighborGap){
            *outRejectedByNeighborGap = true;
        }
        if(!allowCollapsedNeighborShell){
            return false;
        }
    }

    NeighborBondArray runtimeNeighborArray;
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

    constexpr int kExhaustiveOrientationFallbackCoordinationLimit = 8;

    std::vector<int> canonicalToRuntime;
    Matrix3 canonicalOrientation = Matrix3::Identity();
    bool canonicalOrientationValid = false;
    int canonicalSymmetryPermutation = -1;
    const bool useExhaustiveOrientationFallback =
        matcher.requiresOrientationFallback &&
        matcher.coordinationNumber <= kExhaustiveOrientationFallbackCoordinationLimit;

    if(useExhaustiveOrientationFallback){
        if(!assignBestOrientationFallbackPermutation(
            matcher,
            runtimeSignatures,
            runtimeNeighborArray,
            query,
            canonicalToRuntime,
            canonicalOrientation
        )){
            return false;
        }
        canonicalOrientationValid = true;
    }else{
        if(!assignGenericCnaPermutation(
            matcher,
            runtimeSignatures,
            runtimeNeighborArray,
            canonicalToRuntime
        )){
            return false;
        }
    }

    auto buildCanonicalNeighborAtomIndices = [&](const std::vector<int>& slotMapping){
        std::vector<int> indices(static_cast<std::size_t>(matcher.coordinationNumber), -1);
        for(int canonicalSlot = 0; canonicalSlot < matcher.coordinationNumber; ++canonicalSlot){
            indices[static_cast<std::size_t>(canonicalSlot)] =
                query.results()[slotMapping[static_cast<std::size_t>(canonicalSlot)]].index;
        }
        return indices;
    };

    std::vector<int> canonicalNeighborAtomIndices = buildCanonicalNeighborAtomIndices(canonicalToRuntime);
    if(!matchSpeciesWithSymmetry(
        matcher,
        canonicalNeighborAtomIndices,
        atomTypes,
        atomIndex,
        true
    )){
        return false;
    }

    if(!useExhaustiveOrientationFallback && matcher.requiresOrientationFallback && !matcher.symmetries.empty()){
        std::vector<int> bestCanonicalToRuntime = canonicalToRuntime;
        Matrix3 bestOrientation = Matrix3::Identity();
        double bestError = std::numeric_limits<double>::infinity();
        bool foundBest = false;

        for(const PatternSymmetryPermutation& symmetry : matcher.symmetries){
            std::vector<int> candidateCanonicalToRuntime;
            permuteCanonicalToRuntimeBySymmetry(
                canonicalToRuntime,
                symmetry,
                matcher.coordinationNumber,
                candidateCanonicalToRuntime
            );

            Matrix3 candidateOrientation;
            if(!computeMatchOrientation(
                matcher,
                query,
                candidateCanonicalToRuntime,
                candidateOrientation
            )){
                continue;
            }

            const double candidateError = orientationAssignmentError(
                matcher,
                query,
                candidateCanonicalToRuntime,
                candidateOrientation
            );
            const bool preferCandidate =
                !foundBest ||
                candidateError < bestError - 1e-8 ||
                (std::abs(candidateError - bestError) <= 1e-8 &&
                 lexicographicallyGreaterMatrix(candidateOrientation, bestOrientation));

            if(!preferCandidate){
                continue;
            }

            foundBest = true;
            bestError = candidateError;
            bestOrientation = candidateOrientation;
            bestCanonicalToRuntime = std::move(candidateCanonicalToRuntime);
        }

        if(foundBest){
            canonicalToRuntime = std::move(bestCanonicalToRuntime);
            canonicalNeighborAtomIndices = buildCanonicalNeighborAtomIndices(canonicalToRuntime);
            canonicalOrientation = bestOrientation;
            canonicalOrientationValid = true;
        }
    }

    if(matcher.requiresOrientationFallback && canonicalOrientationValid){
        reduceOrientationToFundamentalZone(
            matcher,
            query,
            canonicalToRuntime,
            canonicalOrientation
        );
        canonicalNeighborAtomIndices = buildCanonicalNeighborAtomIndices(canonicalToRuntime);
    }

    outMatch = PatternAtomMatch{};
    outMatch.patternId = pattern.id;
    outMatch.structureType = pattern.structureType;
    outMatch.localMatcherIndex = localMatcherIndex;
    outMatch.localCutoff = localCutoff;
    outMatch.coordinationNumber = matcher.coordinationNumber;
    outMatch.allowedSymmetryMask = AnalysisSymmetryUtils::fullSymmetryMask(static_cast<int>(matcher.symmetries.size()));
    outMatch.symmetryPermutation = canonicalSymmetryPermutation;
    outMatch.orientationValid = canonicalOrientationValid;
    outMatch.orientation = canonicalOrientation;
    if(!outMatch.orientationValid){
        outMatch.symmetryPermutation = -1;
        outMatch.orientationValid = computeMatchOrientation(
            matcher,
            query,
            canonicalToRuntime,
            outMatch.orientation
        );
    }
    if(outMatch.orientationValid && !matcher.symmetries.empty() && !matcher.requiresOrientationFallback){
        if(outMatch.symmetryPermutation < 0){
            outMatch.symmetryPermutation = AnalysisSymmetryUtils::findClosestSymmetryPermutation(
                matcher.symmetries,
                outMatch.orientation
            );
        }
    }
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
        dxaState.orientation = match.orientation;
        dxaState.orientationValid = match.orientationValid;
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
