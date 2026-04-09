#pragma once

#include <volt/analysis/cluster_rule_provider.h>
#include <volt/analysis/pattern_catalog.h>

#include <map>
#include <memory>
#include <vector>

namespace Volt {

class PatternClusterRuleProvider final : public ClusterRuleProvider {
public:
    PatternClusterRuleProvider(
        std::shared_ptr<const PatternCatalog> catalog,
        std::shared_ptr<const std::vector<PatternDxaAtomState>> atomStates
    );

    void initializeClusterSeed(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        Cluster& cluster,
        int seedAtomIndex,
        int structureType
    ) const override;

    bool finalizeClusterOrientation(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        Cluster& cluster,
        int seedAtomIndex,
        int structureType
    ) const override;

    ClusterRuleDecision tryAssignNeighbor(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        const Cluster& cluster,
        int currentAtomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        int structureType,
        int& outNeighborSymmetry
    ) const override;

    ClusterRuleDecision tryCalculateTransition(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        int atomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        Matrix3& outTransition
    ) const override;

private:
    using MatcherKey = std::pair<int, int>;
    using MatcherPairKey = std::pair<MatcherKey, MatcherKey>;

    static std::vector<Matrix3> buildCandidateTransitions(
        const CompiledPatternLocalMatcher& lhs,
        const CompiledPatternLocalMatcher& rhs
    );
    static Matrix3 quantizeTransition(
        const Matrix3& rawTransition,
        const std::vector<Matrix3>& candidates
    );

    bool computeAtomOrientation(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        const PatternDxaAtomState& state,
        int atomIndex,
        Matrix3& outOrientation
    ) const;
    const CompiledPatternLocalMatcher* matcherFor(const PatternDxaAtomState& state) const;
    Vector3 transformedIdealVector(
        const StructureAnalysis& analysis,
        const PatternDxaAtomState& state,
        int symmetryIndex,
        int slot
    ) const;

    int findReverseNeighborSlot(const PatternDxaAtomState& state, int atomIndex) const;
    bool validateLocalOverlap(
        const StructureAnalysis& analysis,
        const PatternDxaAtomState& currentState,
        int currentSymmetry,
        int neighborSlot,
        const PatternDxaAtomState& neighborState,
        int neighborSymmetry,
        int reverseSlot
    ) const;

    std::shared_ptr<const PatternCatalog> _catalog;
    std::shared_ptr<const std::vector<PatternDxaAtomState>> _atomStates;
    std::map<MatcherPairKey, std::vector<Matrix3>> _transitionCandidates;
};

}
