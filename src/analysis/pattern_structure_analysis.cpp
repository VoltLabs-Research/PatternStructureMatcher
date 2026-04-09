#include <volt/analysis/pattern_structure_analysis.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>

namespace Volt {

void computeMaximumNeighborDistanceFromPatterns(StructureAnalysis& analysis){
    StructureContext& context = analysis.context();
    const std::size_t atomCount = context.atomCount();
    if(atomCount == 0 || !context.neighborCounts || !context.neighborOffsets || !context.neighborIndices){
        context.maximumNeighborDistance = 0.0;
        return;
    }

    const auto* positions = context.positions->constDataPoint3();
    const auto& inverseMatrix = context.simCell.inverseMatrix();
    const auto& directMatrix = context.simCell.matrix();
    const int* counts = context.neighborCounts->constDataInt();
    const int* offsets = context.neighborOffsets->constDataInt();
    const int* indices = context.neighborIndices->constDataInt();

    context.maximumNeighborDistance = tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, atomCount),
        0.0,
        [&](const tbb::blocked_range<std::size_t>& range, double maxDistance) -> double {
            for(std::size_t atomIndex = range.begin(); atomIndex != range.end(); ++atomIndex){
                const int neighborCount = counts[atomIndex];
                const int start = offsets[atomIndex];
                for(int neighborSlot = 0; neighborSlot < neighborCount; ++neighborSlot){
                    const int neighborIndex = indices[start + neighborSlot];
                    Vector3 delta = positions[neighborIndex] - positions[atomIndex];

                    double fractional[3] = {0.0, 0.0, 0.0};
                    for(int dimension = 0; dimension < 3; ++dimension){
                        fractional[dimension] = inverseMatrix.prodrow(delta, dimension);
                        fractional[dimension] -= std::round(fractional[dimension]);
                    }

                    Vector3 minimumImage =
                        directMatrix.column(0) * fractional[0] +
                        directMatrix.column(1) * fractional[1] +
                        directMatrix.column(2) * fractional[2];

                    maxDistance = std::max(maxDistance, minimumImage.length());
                }
            }
            return maxDistance;
        },
        [](double lhs, double rhs) -> double {
            return std::max(lhs, rhs);
        }
    );
}

}
