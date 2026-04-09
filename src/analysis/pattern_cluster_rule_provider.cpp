#include <volt/analysis/pattern_cluster_rule_provider.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/structures/cluster.h>

#include <algorithm>

namespace Volt {

constexpr double kPatternClusterOrientationTolerance = 0.2;
constexpr double kPatternNeighborVectorTolerance = 1e-4;

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

bool findMatchingIdealOffset(
    const PatternDxaAtomState& sourceState,
    int sourceSlot,
    const PatternDxaAtomState& targetState,
    int excludedTargetSlot,
    int& outTargetSlot
){
    const Vector3 expectedOffset = sourceState.idealNeighborVectors[static_cast<std::size_t>(sourceSlot)] -
        sourceState.idealNeighborVectors[static_cast<std::size_t>(excludedTargetSlot)];

    for(int targetSlot = 0; targetSlot < targetState.coordinationNumber; ++targetSlot){
        if(targetSlot == excludedTargetSlot){
            continue;
        }
        if((expectedOffset - targetState.idealNeighborVectors[static_cast<std::size_t>(targetSlot)]).isZero(kPatternNeighborVectorTolerance)){
            outTargetSlot = targetSlot;
            return true;
        }
    }

    return false;
}

bool computeIdealTransition(
    const PatternDxaAtomState& currentState,
    int neighborSlot,
    const PatternDxaAtomState& neighborState,
    int reverseSlot,
    Matrix3& outTransition
){
    std::array<std::array<int, 2>, 2> commonNeighborPairs{};
    int commonNeighborCount = 0;

    for(int currentSlot = 0; currentSlot < currentState.coordinationNumber; ++currentSlot){
        if(currentSlot == neighborSlot){
            continue;
        }

        const Vector3 expectedFromNeighbor =
            currentState.idealNeighborVectors[static_cast<std::size_t>(currentSlot)] -
            currentState.idealNeighborVectors[static_cast<std::size_t>(neighborSlot)];

        for(int neighborSlotCandidate = 0; neighborSlotCandidate < neighborState.coordinationNumber; ++neighborSlotCandidate){
            if(neighborSlotCandidate == reverseSlot){
                continue;
            }
            if(!(expectedFromNeighbor - neighborState.idealNeighborVectors[static_cast<std::size_t>(neighborSlotCandidate)]).isZero(kPatternNeighborVectorTolerance)){
                continue;
            }

            if(commonNeighborCount < 2){
                commonNeighborPairs[static_cast<std::size_t>(commonNeighborCount)] = { currentSlot, neighborSlotCandidate };
            }
            ++commonNeighborCount;
            break;
        }
    }

    if(commonNeighborCount < 2){
        return false;
    }

    Matrix3 tm1 = Matrix3::Zero();
    Matrix3 tm2 = Matrix3::Zero();
    for(int axis = 0; axis < 2; ++axis){
        const int currentSlot = commonNeighborPairs[static_cast<std::size_t>(axis)][0];
        const int neighborSlotCandidate = commonNeighborPairs[static_cast<std::size_t>(axis)][1];
        tm1.column(axis) =
            currentState.idealNeighborVectors[static_cast<std::size_t>(currentSlot)] -
            currentState.idealNeighborVectors[static_cast<std::size_t>(neighborSlot)];
        tm2.column(axis) = neighborState.idealNeighborVectors[static_cast<std::size_t>(neighborSlotCandidate)];
    }
    tm1.column(2) = -currentState.idealNeighborVectors[static_cast<std::size_t>(neighborSlot)];
    tm2.column(2) = neighborState.idealNeighborVectors[static_cast<std::size_t>(reverseSlot)];

    Matrix3 tm1Inverse;
    if(!tm1.inverse(tm1Inverse)){
        return false;
    }
    outTransition = Matrix3(tm2 * tm1Inverse);
    return outTransition.isOrthogonalMatrix(kPatternClusterOrientationTolerance);
}

std::vector<std::array<int, 3>> enumerateIndependentBases(const std::vector<Vector3>& vectors){
    std::vector<std::array<int, 3>> bases;
    for(int i = 0; i < static_cast<int>(vectors.size()); ++i){
        for(int j = i + 1; j < static_cast<int>(vectors.size()); ++j){
            for(int k = j + 1; k < static_cast<int>(vectors.size()); ++k){
                Matrix3 basis(vectors[static_cast<std::size_t>(i)], vectors[static_cast<std::size_t>(j)], vectors[static_cast<std::size_t>(k)]);
                if(std::abs(basis.determinant()) <= EPSILON){
                    continue;
                }
                bases.push_back({i, j, k});
            }
        }
    }
    return bases;
}

PatternClusterRuleProvider::PatternClusterRuleProvider(
    std::shared_ptr<const PatternCatalog> catalog,
    std::shared_ptr<const std::vector<PatternDxaAtomState>> atomStates
) :
    _catalog(std::move(catalog)),
    _atomStates(std::move(atomStates))
{
    if(!_catalog){
        return;
    }

    for(const auto& lhsPattern : _catalog->patterns()){
        for(std::size_t lhsMatcherIndex = 0; lhsMatcherIndex < lhsPattern.localMatchers.size(); ++lhsMatcherIndex){
            for(const auto& rhsPattern : _catalog->patterns()){
                for(std::size_t rhsMatcherIndex = 0; rhsMatcherIndex < rhsPattern.localMatchers.size(); ++rhsMatcherIndex){
                    const MatcherPairKey key = {
                        {lhsPattern.id, static_cast<int>(lhsMatcherIndex)},
                        {rhsPattern.id, static_cast<int>(rhsMatcherIndex)}
                    };
                    _transitionCandidates.emplace(
                        key,
                        buildCandidateTransitions(
                            lhsPattern.localMatchers[lhsMatcherIndex],
                            rhsPattern.localMatchers[rhsMatcherIndex]
                        )
                    );
                }
            }
        }
    }
}

std::vector<Matrix3> PatternClusterRuleProvider::buildCandidateTransitions(
    const CompiledPatternLocalMatcher& lhs,
    const CompiledPatternLocalMatcher& rhs
){
    std::vector<Matrix3> candidates;
    const auto lhsBases = enumerateIndependentBases(lhs.canonicalNeighborVectors);
    const auto rhsBases = enumerateIndependentBases(rhs.canonicalNeighborVectors);

    for(const auto& lhsBasisIndices : lhsBases){
        Matrix3 lhsBasis(
            lhs.canonicalNeighborVectors[static_cast<std::size_t>(lhsBasisIndices[0])],
            lhs.canonicalNeighborVectors[static_cast<std::size_t>(lhsBasisIndices[1])],
            lhs.canonicalNeighborVectors[static_cast<std::size_t>(lhsBasisIndices[2])]
        );
        Matrix3 lhsBasisInverse;
        if(!lhsBasis.inverse(lhsBasisInverse)){
            continue;
        }

        for(const auto& rhsBasisIndices : rhsBases){
            Matrix3 rhsBasis(
                rhs.canonicalNeighborVectors[static_cast<std::size_t>(rhsBasisIndices[0])],
                rhs.canonicalNeighborVectors[static_cast<std::size_t>(rhsBasisIndices[1])],
                rhs.canonicalNeighborVectors[static_cast<std::size_t>(rhsBasisIndices[2])]
            );
            Matrix3 transition = Matrix3(rhsBasis * lhsBasisInverse);
            if(!transition.isOrthogonalMatrix(kPatternClusterOrientationTolerance)){
                continue;
            }

            Matrix3 canonicalTransition;
            if(orthonormalizeMatrix(transition, canonicalTransition)){
                transition = canonicalTransition;
            }

            const auto duplicate = std::find_if(candidates.begin(), candidates.end(), [&](const Matrix3& existing){
                return existing.equals(transition, CA_TRANSITION_MATRIX_EPSILON);
            });
            if(duplicate == candidates.end()){
                candidates.push_back(transition);
            }
        }
    }

    if(candidates.empty()){
        candidates.push_back(Matrix3::Identity());
    }
    return candidates;
}

Matrix3 PatternClusterRuleProvider::quantizeTransition(
    const Matrix3& rawTransition,
    const std::vector<Matrix3>& candidates
){
    Matrix3 best = Matrix3::Identity();
    double bestError = DOUBLE_MAX;
    for(const Matrix3& candidate : candidates){
        double error = 0.0;
        for(int row = 0; row < 3; ++row){
            for(int col = 0; col < 3; ++col){
                const double diff = candidate(row, col) - rawTransition(row, col);
                error += diff * diff;
            }
        }
        if(error < bestError){
            bestError = error;
            best = candidate;
        }
    }
    return best;
}

bool orthogonalizeTransition(const Matrix3& input, Matrix3& output){
    return orthonormalizeMatrix(input, output);
}

void PatternClusterRuleProvider::initializeClusterSeed(
    const StructureAnalysis&,
    const AnalysisContext&,
    Cluster& cluster,
    int,
    int
) const{
    cluster.symmetryTransformation = 0;
}

bool PatternClusterRuleProvider::computeAtomOrientation(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    const PatternDxaAtomState& state,
    int atomIndex,
    Matrix3& outOrientation
) const{
    if(!state.isMatched() || state.coordinationNumber < 3){
        return false;
    }

    Matrix3 orientationV = Matrix3::Zero();
    Matrix3 orientationW = Matrix3::Zero();
    int vectorCount = 0;

    const int symmetryIndex = std::max(0, context.atomSymmetryPermutations->getInt(static_cast<std::size_t>(atomIndex)));

    for(int slot = 0; slot < state.coordinationNumber; ++slot){
        const int neighborAtomIndex = state.neighborAtomIndices[static_cast<std::size_t>(slot)];
        if(neighborAtomIndex < 0){
            continue;
        }

        const Vector3 idealVector = transformedIdealVector(analysis, state, symmetryIndex, slot);
        const Vector3 spatialVector = context.simCell.wrapVector(
            context.positions->getPoint3(static_cast<std::size_t>(neighborAtomIndex)) -
            context.positions->getPoint3(static_cast<std::size_t>(atomIndex))
        );

        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 3; ++j){
                orientationV(i, j) += idealVector[j] * idealVector[i];
                orientationW(i, j) += idealVector[j] * spatialVector[i];
            }
        }
        ++vectorCount;
    }

    if(vectorCount < 3){
        return false;
    }

    Matrix3 orientationVInverse;
    if(!orientationV.inverse(orientationVInverse)){
        return false;
    }

    const Matrix3 rawOrientation = Matrix3(orientationW * orientationVInverse);
    return orthonormalizeMatrix(rawOrientation, outOrientation);
}

const CompiledPatternLocalMatcher* PatternClusterRuleProvider::matcherFor(const PatternDxaAtomState& state) const{
    if(!_catalog || state.patternId < 0){
        return nullptr;
    }
    const auto& pattern = _catalog->patternById(state.patternId);
    if(state.localMatcherIndex < 0 || state.localMatcherIndex >= static_cast<int>(pattern.localMatchers.size())){
        return nullptr;
    }
    return &pattern.localMatchers[static_cast<std::size_t>(state.localMatcherIndex)];
}

Vector3 PatternClusterRuleProvider::transformedIdealVector(
    const StructureAnalysis& analysis,
    const PatternDxaAtomState& state,
    int symmetryIndex,
    int slot
) const{
    if(slot < 0 || slot >= state.coordinationNumber){
        return Vector3::Zero();
    }

    if(!isSyntheticPatternStructureType(state.structureType)){
        const int mappedSlot = analysis.symmetryPermutationEntry(state.structureType, symmetryIndex, slot);
        if(mappedSlot < 0 || mappedSlot >= state.coordinationNumber){
            return Vector3::Zero();
        }
        return analysis.latticeVector(state.structureType, mappedSlot);
    }

    const auto* matcher = matcherFor(state);
    if(!matcher || symmetryIndex < 0 || symmetryIndex >= static_cast<int>(matcher->symmetries.size())){
        return state.idealNeighborVectors[static_cast<std::size_t>(slot)];
    }

    const int mappedSlot = matcher->symmetries[static_cast<std::size_t>(symmetryIndex)].permutation[static_cast<std::size_t>(slot)];
    if(mappedSlot < 0 || mappedSlot >= state.coordinationNumber){
        return Vector3::Zero();
    }
    return state.idealNeighborVectors[static_cast<std::size_t>(mappedSlot)];
}

int PatternClusterRuleProvider::findReverseNeighborSlot(const PatternDxaAtomState& state, int atomIndex) const{
    for(int slot = 0; slot < state.coordinationNumber; ++slot){
        if(state.neighborAtomIndices[static_cast<std::size_t>(slot)] == atomIndex){
            return slot;
        }
    }
    return -1;
}

bool PatternClusterRuleProvider::validateLocalOverlap(
    const StructureAnalysis& analysis,
    const PatternDxaAtomState& currentState,
    int currentSymmetry,
    int neighborSlot,
    const PatternDxaAtomState& neighborState,
    int neighborSymmetry,
    int reverseSlot
) const{
    if(!(transformedIdealVector(analysis, currentState, currentSymmetry, neighborSlot) +
         transformedIdealVector(analysis, neighborState, neighborSymmetry, reverseSlot)).isZero(kPatternNeighborVectorTolerance)){
        return false;
    }

    int commonNeighborCount = 0;
    for(int currentSlot = 0; currentSlot < currentState.coordinationNumber; ++currentSlot){
        if(currentSlot == neighborSlot){
            continue;
        }

        int matchingNeighborSlot = -1;
        const Vector3 expectedFromNeighbor =
            transformedIdealVector(analysis, currentState, currentSymmetry, currentSlot) -
            transformedIdealVector(analysis, currentState, currentSymmetry, neighborSlot);

        for(int neighborCandidateSlot = 0; neighborCandidateSlot < neighborState.coordinationNumber; ++neighborCandidateSlot){
            if(neighborCandidateSlot == reverseSlot){
                continue;
            }
            if((expectedFromNeighbor - transformedIdealVector(analysis, neighborState, neighborSymmetry, neighborCandidateSlot)).isZero(kPatternNeighborVectorTolerance)){
                matchingNeighborSlot = neighborCandidateSlot;
                break;
            }
        }

        if(matchingNeighborSlot < 0){
            continue;
        }
        ++commonNeighborCount;
    }

    return commonNeighborCount >= 2;
}

bool PatternClusterRuleProvider::finalizeClusterOrientation(
    const StructureAnalysis&,
    const AnalysisContext& context,
    Cluster& cluster,
    int,
    int
) const{
    (void)context;
    (void)cluster;
    return false;
}

ClusterRuleDecision PatternClusterRuleProvider::tryAssignNeighbor(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    const Cluster&,
    int currentAtomIndex,
    int neighborAtomIndex,
    int neighborIndex,
    int structureType,
    int& outNeighborSymmetry
) const{
    const auto& currentState = (*_atomStates)[static_cast<std::size_t>(currentAtomIndex)];
    const auto& neighborState = (*_atomStates)[static_cast<std::size_t>(neighborAtomIndex)];
    if(!currentState.isMatched() || !neighborState.isMatched()){
        return ClusterRuleDecision::Unhandled;
    }
    if(currentState.structureType != structureType || neighborState.structureType != structureType){
        return ClusterRuleDecision::Rejected;
    }
    if(neighborIndex < 0 || neighborIndex >= currentState.coordinationNumber){
        return ClusterRuleDecision::Rejected;
    }
    if(currentState.neighborAtomIndices[static_cast<std::size_t>(neighborIndex)] != neighborAtomIndex){
        return ClusterRuleDecision::Rejected;
    }

    const int reverseSlot = findReverseNeighborSlot(neighborState, currentAtomIndex);
    if(reverseSlot < 0){
        return ClusterRuleDecision::Rejected;
    }

    const int currentSymmetry = std::max(0, context.atomSymmetryPermutations->getInt(static_cast<std::size_t>(currentAtomIndex)));
    const auto* neighborMatcher = matcherFor(neighborState);
    const std::uint64_t allowedMask = static_cast<std::uint64_t>(
        context.atomAllowedSymmetryMasks->getInt64(static_cast<std::size_t>(neighborAtomIndex))
    );
    const int neighborSymmetryCount = isSyntheticPatternStructureType(neighborState.structureType)
        ? (neighborMatcher ? static_cast<int>(neighborMatcher->symmetries.size()) : 1)
        : std::max(1, analysis.symmetryPermutationCount(structureType));

    for(int neighborSymmetry = 0; neighborSymmetry < std::max(1, neighborSymmetryCount); ++neighborSymmetry){
        if(neighborSymmetry < 63 && allowedMask != 0 && (allowedMask & (std::uint64_t{1} << neighborSymmetry)) == 0){
            continue;
        }
        if(!validateLocalOverlap(analysis, currentState, currentSymmetry, neighborIndex, neighborState, neighborSymmetry, reverseSlot)){
            continue;
        }
        outNeighborSymmetry = neighborSymmetry;
        return ClusterRuleDecision::Accepted;
    }

    return ClusterRuleDecision::Rejected;
}

ClusterRuleDecision PatternClusterRuleProvider::tryCalculateTransition(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    int atomIndex,
    int neighborAtomIndex,
    int neighborIndex,
    Matrix3& outTransition
) const{
    const auto& atomState = (*_atomStates)[static_cast<std::size_t>(atomIndex)];
    const auto& neighborState = (*_atomStates)[static_cast<std::size_t>(neighborAtomIndex)];
    if(!atomState.isMatched() || !neighborState.isMatched()){
        return ClusterRuleDecision::Unhandled;
    }
    if(
        !isSyntheticPatternStructureType(atomState.structureType) &&
        !isSyntheticPatternStructureType(neighborState.structureType)
    ){
        return ClusterRuleDecision::Unhandled;
    }

    const int reverseSlot = findReverseNeighborSlot(neighborState, atomIndex);
    if(reverseSlot >= 0 && computeIdealTransition(atomState, neighborIndex, neighborState, reverseSlot, outTransition)){
        const MatcherPairKey key = {
            {atomState.patternId, atomState.localMatcherIndex},
            {neighborState.patternId, neighborState.localMatcherIndex}
        };
        const auto it = _transitionCandidates.find(key);
        if(it != _transitionCandidates.end() && !it->second.empty()){
            outTransition = quantizeTransition(outTransition, it->second);
        }else{
            Matrix3 orthogonalized;
            if(orthogonalizeTransition(outTransition, orthogonalized)){
                outTransition = orthogonalized;
            }
        }
        if(outTransition.isOrthogonalMatrix(CA_TRANSITION_MATRIX_EPSILON)){
            return ClusterRuleDecision::Accepted;
        }
    }

    Matrix3 atomOrientation;
    Matrix3 neighborOrientation;
    if(computeAtomOrientation(analysis, context, atomState, atomIndex, atomOrientation) &&
       computeAtomOrientation(analysis, context, neighborState, neighborAtomIndex, neighborOrientation)){
        const Matrix3 transition = Matrix3(neighborOrientation.transposed() * atomOrientation);
        if(transition.isOrthogonalMatrix(kPatternClusterOrientationTolerance)){
            const MatcherPairKey key = {
                {atomState.patternId, atomState.localMatcherIndex},
                {neighborState.patternId, neighborState.localMatcherIndex}
            };
            const auto it = _transitionCandidates.find(key);
            Matrix3 quantized = it == _transitionCandidates.end()
                ? transition
                : quantizeTransition(transition, it->second);
            outTransition = quantized;
            if(outTransition.isOrthogonalMatrix(CA_TRANSITION_MATRIX_EPSILON)){
                return ClusterRuleDecision::Accepted;
            }
        }
    }

    outTransition = Matrix3::Identity();
    return ClusterRuleDecision::Accepted;
}

}
