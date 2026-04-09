#pragma once

#include <volt/analysis/pattern_catalog.h>
#include <volt/analysis/structure_analysis.h>

#include <memory>
#include <unordered_map>

namespace Volt {

class PatternDxaTopologyProvider final : public StructureAnalysisCrystalInfo {
public:
    struct TopologyData {
        int coordinationNumber = 0;
        std::array<std::array<int, 2>, MAX_NEIGHBORS> commonNeighbors{};
        std::vector<PatternSymmetryPermutation> symmetries;
        std::vector<Vector3> neighborVectors;
        std::string name;
    };

    explicit PatternDxaTopologyProvider(
        std::shared_ptr<const PatternCatalog> catalog,
        std::vector<int> selectedPatternIds = {}
    );

    int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const override;
    int coordinationNumber(int structureType) const override;
    int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const override;
    int symmetryPermutationCount(int structureType) const override;
    int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const override;
    const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const override;
    int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const override;
    const Vector3& latticeVector(int structureType, int latticeVectorIndex) const override;
    std::string_view topologyName(int structureType) const override;

private:
    const TopologyData& dataFor(int structureType) const;
    static TopologyData buildTopologyData(const CompiledPattern& pattern);

    std::shared_ptr<const PatternCatalog> _catalog;
    std::unordered_map<int, TopologyData> _topologies;
};

}
