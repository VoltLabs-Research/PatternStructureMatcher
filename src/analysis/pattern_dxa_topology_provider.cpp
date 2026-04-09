#include <volt/analysis/pattern_dxa_topology_provider.h>

#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/structures/neighbor_bond_array.h>

#include <stdexcept>

namespace Volt {

namespace{

const Vector3& zeroVector(){
    static const Vector3 vector = Vector3::Zero();
    return vector;
}

const Matrix3& identityMatrix(){
    static const Matrix3 matrix = Matrix3::Identity();
    return matrix;
}

}

PatternDxaTopologyProvider::PatternDxaTopologyProvider(
    std::shared_ptr<const PatternCatalog> catalog,
    std::vector<int> selectedPatternIds
)
    : _catalog(std::move(catalog))
{
    if(!_catalog){
        return;
    }

    const auto visitPattern = [&](const CompiledPattern& pattern){
        const TopologyData topology = buildTopologyData(pattern);
        const auto [it, inserted] = _topologies.emplace(pattern.structureType, topology);
        if(!inserted && it->second.coordinationNumber <= 0 && topology.coordinationNumber > 0){
            it->second = topology;
        }
    };

    if(selectedPatternIds.empty()){
        for(const CompiledPattern& pattern : _catalog->patterns()){
            visitPattern(pattern);
        }
        return;
    }

    for(const int patternId : selectedPatternIds){
        if(patternId < 0 || patternId >= static_cast<int>(_catalog->patterns().size())){
            continue;
        }
        visitPattern(_catalog->patternById(patternId));
    }
}

PatternDxaTopologyProvider::TopologyData PatternDxaTopologyProvider::buildTopologyData(const CompiledPattern& pattern){
    TopologyData data;
    data.commonNeighbors.fill({-1, -1});
    data.name = pattern.name;

    const CompiledPatternLocalMatcher* referenceMatcher = nullptr;
    for(const auto& matcher : pattern.localMatchers){
        if(!referenceMatcher || matcher.coordinationNumber > referenceMatcher->coordinationNumber){
            referenceMatcher = &matcher;
        }
    }

    if(!referenceMatcher){
        return data;
    }

    data.coordinationNumber = referenceMatcher->coordinationNumber;
    data.neighborVectors = referenceMatcher->canonicalNeighborVectors;
    data.symmetries = referenceMatcher->symmetries;
    if(data.symmetries.empty()){
        PatternSymmetryPermutation identity;
        for(int slot = 0; slot < data.coordinationNumber && slot < static_cast<int>(identity.permutation.size()); ++slot){
            identity.permutation[static_cast<std::size_t>(slot)] = slot;
        }
        identity.inverseProduct = {0};
        data.symmetries.push_back(std::move(identity));
    }

    NeighborBondArray neighborArray;
    for(int row = 0; row < static_cast<int>(referenceMatcher->neighborBondRows.size()); ++row){
        neighborArray.neighborArray[row] = referenceMatcher->neighborBondRows[static_cast<std::size_t>(row)];
    }

    for(int neighborIndex = 0; neighborIndex < data.coordinationNumber; ++neighborIndex){
        Matrix3 basis;
        basis.column(0) = data.neighborVectors[static_cast<std::size_t>(neighborIndex)];
        bool found = false;

        for(int i1 = 0; i1 < data.coordinationNumber && !found; ++i1){
            if(!neighborArray.neighborBond(neighborIndex, i1)){
                continue;
            }
            basis.column(1) = data.neighborVectors[static_cast<std::size_t>(i1)];

            for(int i2 = i1 + 1; i2 < data.coordinationNumber; ++i2){
                if(!neighborArray.neighborBond(neighborIndex, i2)){
                    continue;
                }
                basis.column(2) = data.neighborVectors[static_cast<std::size_t>(i2)];

                if(std::abs(basis.determinant()) <= EPSILON){
                    continue;
                }

                data.commonNeighbors[static_cast<std::size_t>(neighborIndex)][0] = i1;
                data.commonNeighbors[static_cast<std::size_t>(neighborIndex)][1] = i2;
                found = true;
                break;
            }
        }
    }

    return data;
}

const PatternDxaTopologyProvider::TopologyData& PatternDxaTopologyProvider::dataFor(int structureType) const{
    const auto it = _topologies.find(structureType);
    if(it == _topologies.end()){
        throw std::out_of_range("Unknown pattern structure type");
    }
    return it->second;
}

int PatternDxaTopologyProvider::findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const{
    const auto& data = dataFor(structureType);
    return AnalysisSymmetryUtils::findClosestSymmetryPermutation(data.symmetries, rotation);
}

int PatternDxaTopologyProvider::coordinationNumber(int structureType) const{
    return dataFor(structureType).coordinationNumber;
}

int PatternDxaTopologyProvider::commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const{
    const auto& data = dataFor(structureType);
    if(neighborIndex < 0 || neighborIndex >= data.coordinationNumber || commonNeighborSlot < 0 || commonNeighborSlot > 1){
        return -1;
    }
    return data.commonNeighbors[static_cast<std::size_t>(neighborIndex)][static_cast<std::size_t>(commonNeighborSlot)];
}

int PatternDxaTopologyProvider::symmetryPermutationCount(int structureType) const{
    return static_cast<int>(dataFor(structureType).symmetries.size());
}

int PatternDxaTopologyProvider::symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size()) || neighborIndex < 0 || neighborIndex >= data.coordinationNumber){
        return neighborIndex;
    }
    return data.symmetries[static_cast<std::size_t>(symmetryIndex)].permutation[static_cast<std::size_t>(neighborIndex)];
}

const Matrix3& PatternDxaTopologyProvider::symmetryTransformation(int structureType, int symmetryIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size())){
        return identityMatrix();
    }
    return data.symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
}

int PatternDxaTopologyProvider::symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size())){
        return 0;
    }

    const auto& inverseProduct = data.symmetries[static_cast<std::size_t>(symmetryIndex)].inverseProduct;
    if(transformationIndex < 0 || transformationIndex >= static_cast<int>(inverseProduct.size())){
        return 0;
    }
    return inverseProduct[static_cast<std::size_t>(transformationIndex)];
}

const Vector3& PatternDxaTopologyProvider::latticeVector(int structureType, int latticeVectorIndex) const{
    const auto& data = dataFor(structureType);
    if(latticeVectorIndex < 0 || latticeVectorIndex >= static_cast<int>(data.neighborVectors.size())){
        return zeroVector();
    }
    return data.neighborVectors[static_cast<std::size_t>(latticeVectorIndex)];
}

std::string_view PatternDxaTopologyProvider::topologyName(int structureType) const{
    return dataFor(structureType).name;
}

}
