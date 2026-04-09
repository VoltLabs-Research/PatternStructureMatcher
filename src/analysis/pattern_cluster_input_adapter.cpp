#include <volt/analysis/pattern_cluster_input_adapter.h>
#include <volt/analysis/cluster_input_preparation.h>

namespace Volt {

void PatternClusterInputAdapter::prepare(StructureAnalysis& analysis, AnalysisContext& context){
    ClusterInputAdapterUtils::prepareSymmetryAwareClusterInputs(
        analysis,
        context,
        false,
        [&](std::size_t atomIndex, int structureType) {
            if(structureType == LATTICE_OTHER){
                return false;
            }
            if(analysis.numberOfNeighbors(static_cast<int>(atomIndex)) == 0){
                return false;
            }
            return context.atomAllowedSymmetryMasks->getInt64(atomIndex) == 0;
        }
    );
}

}
