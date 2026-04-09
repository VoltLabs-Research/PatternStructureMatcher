#include <volt/analysis/pattern_service.h>

#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/analysis/pattern_catalog.h>
#include <volt/analysis/pattern_classifier.h>
#include <volt/analysis/pattern_cluster_input_adapter.h>
#include <volt/analysis/pattern_dxa_topology_provider.h>
#include <volt/analysis/orientation_cluster_rule_provider.h>
#include <volt/analysis/pattern_structure_analysis.h>
#include <volt/analysis/reconstructed_state_canonicalizer.h>
#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/json_utils.h>
#include <volt/structures/crystal_topology_registry.h>

#include <algorithm>
#include <cctype>
#include <map>
#include <unordered_set>

namespace Volt {

std::string trimAscii(std::string_view value){
    std::size_t first = 0;
    while(first < value.size() && std::isspace(static_cast<unsigned char>(value[first]))){
        ++first;
    }
    std::size_t last = value.size();
    while(last > first && std::isspace(static_cast<unsigned char>(value[last - 1]))){
        --last;
    }
    return std::string(value.substr(first, last - first));
}

std::vector<std::string> splitPatternCsv(std::string_view csv){
    std::vector<std::string> values;
    std::size_t start = 0;
    while(start <= csv.size()){
        const std::size_t comma = csv.find(',', start);
        const std::size_t end = comma == std::string_view::npos ? csv.size() : comma;
        std::string token = trimAscii(csv.substr(start, end - start));
        if(!token.empty()){
            values.push_back(std::move(token));
        }
        if(comma == std::string_view::npos){
            break;
        }
        start = comma + 1;
    }
    return values;
}

bool resolveSelectedPatternIds(
    const PatternCatalog& catalog,
    const std::string& csv,
    std::vector<int>& outPatternIds,
    std::string& outError
){
    outPatternIds.clear();
    const auto tokens = splitPatternCsv(csv);
    if(tokens.empty()){
        outError = "--patterns was provided but no valid pattern names were found.";
        return false;
    }

    std::unordered_set<int> seenIds;
    for(const std::string& rawName : tokens){
        const std::string normalized = normalizePatternName(rawName);
        if(normalized == "all"){
            outPatternIds = catalog.defaultSelection();
            return !outPatternIds.empty();
        }

        const CompiledPattern* pattern = catalog.findPatternByName(normalized);
        if(!pattern){
            outError = "Unknown pattern in --patterns: '" + rawName + "'";
            return false;
        }
        if(!pattern->supportedForLocalMatching){
            continue;
        }
        if(seenIds.insert(pattern->id).second){
            outPatternIds.push_back(pattern->id);
        }
    }

    if(outPatternIds.empty()){
        outError = "No supported lattice patterns remained after applying --patterns filter.";
        return false;
    }
    return true;
}

int structureTypeForPattern(const CompiledPattern& pattern){
    return pattern.structureType;
}

json buildMainListing(
    const AnalysisContext& context,
    const PatternCatalog& catalog,
    const std::vector<int>& atomPatternIds
){
    struct ListingEntry {
        std::string name;
        int label = static_cast<int>(StructureType::OTHER);
        int count = 0;
    };

    std::map<int, ListingEntry> counts;
    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int patternId = atomIndex < atomPatternIds.size() ? atomPatternIds[atomIndex] : -1;
        if(patternId >= 0){
            const auto& pattern = catalog.patternById(patternId);
            const int structureLabel = context.structureTypes->getInt(atomIndex);
            auto& entry = counts[structureLabel];
            entry.label = structureLabel;
            entry.name = pattern.name;
            ++entry.count;
            continue;
        }

        auto& entry = counts[static_cast<int>(StructureType::OTHER)];
        entry.label = static_cast<int>(StructureType::OTHER);
        entry.name = "OTHER";
        ++entry.count;
    }

    json listing = json::array();
    for(const auto& [_, entry] : counts){
        listing.push_back({
            {"structure_type", entry.label},
            {"structure_name", entry.name},
            {"count", entry.count},
        });
    }

    std::sort(listing.begin(), listing.end(), [](const json& lhs, const json& rhs){
        return lhs.value("structure_name", "") < rhs.value("structure_name", "");
    });
    return listing;
}

json buildPatternListing(const PatternCatalog& catalog, const std::vector<int>& atomPatternIds){
    std::vector<int> counts(catalog.patterns().size(), 0);
    for(int patternId : atomPatternIds){
        if(patternId < 0 || patternId >= static_cast<int>(counts.size())){
            continue;
        }
        counts[static_cast<std::size_t>(patternId)]++;
    }

    json listing = json::array();
    for(const CompiledPattern& pattern : catalog.patterns()){
        const int count = counts[static_cast<std::size_t>(pattern.id)];
        if(count == 0){
            continue;
        }
        listing.push_back({
            {"pattern_id", pattern.id},
            {"pattern_name", pattern.name},
            {"count", count},
        });
    }

    std::sort(listing.begin(), listing.end(), [](const json& lhs, const json& rhs){
        return lhs.value("pattern_name", "") < rhs.value("pattern_name", "");
    });
    return listing;
}

json buildPerAtomProperties(
    const LammpsParser::Frame& frame,
    const AnalysisContext& context,
    const PatternCatalog& catalog,
    const std::vector<int>& atomPatternIds
){
    json perAtom = json::array();
    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int structureType = context.structureTypes->getInt(atomIndex);
        const int patternId = atomPatternIds[atomIndex];
        const CompiledPattern* pattern = nullptr;
        if(patternId >= 0 && patternId < static_cast<int>(catalog.patterns().size())){
            pattern = &catalog.patternById(patternId);
        }

        json atom;
        atom["id"] = atomIndex < frame.ids.size() ? frame.ids[atomIndex] : static_cast<int>(atomIndex);
        atom["structure_type"] = structureType;
        atom["structure_name"] = pattern ? pattern->name : structureTypeNameForExport(structureType);
        atom["pattern_id"] = patternId;
        atom["pattern_name"] = pattern ? pattern->name : "";
        atom["cluster_id"] = context.atomClusters ? context.atomClusters->getInt(atomIndex) : 0;

        if(atomIndex < frame.positions.size()){
            const auto& position = frame.positions[atomIndex];
            atom["pos"] = {position.x(), position.y(), position.z()};
        }else{
            atom["pos"] = {0.0, 0.0, 0.0};
        }

        perAtom.push_back(std::move(atom));
    }
    return perAtom;
}

json buildSelectedPatterns(const PatternCatalog& catalog, const std::vector<int>& selectedPatternIds){
    json selected = json::array();
    for(int patternId : selectedPatternIds){
        const CompiledPattern& pattern = catalog.patternById(patternId);
        selected.push_back({
            {"pattern_id", pattern.id},
            {"structure_type", structureTypeForPattern(pattern)},
            {"pattern_name", pattern.name},
        });
    }
    return selected;
}

void applyPatternNeighborVectorOverrides(
    StructureAnalysis& analysis,
    const AnalysisContext& context,
    const std::shared_ptr<const std::vector<PatternDxaAtomState>>& atomStates
){
    if(!atomStates || atomStates->size() != context.atomCount()){
        return;
    }

    std::vector<Vector3> neighborVectorOverrides(
        context.atomCount() * static_cast<std::size_t>(MAX_NEIGHBORS),
        Vector3::Zero()
    );

    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const PatternDxaAtomState& state = (*atomStates)[atomIndex];
        const int count = std::min(state.coordinationNumber, static_cast<int>(MAX_NEIGHBORS));
        for(int neighborSlot = 0; neighborSlot < count; ++neighborSlot){
            neighborVectorOverrides[
                atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) +
                static_cast<std::size_t>(neighborSlot)
            ] = state.idealNeighborVectors[static_cast<std::size_t>(neighborSlot)];
        }
    }

    analysis.setNeighborLatticeVectorOverrides(
        std::move(neighborVectorOverrides),
        static_cast<std::size_t>(MAX_NEIGHBORS)
    );
}

std::shared_ptr<std::vector<OrientationClusterAtomState>> buildOrientationClusterStates(
    const std::shared_ptr<const std::vector<PatternDxaAtomState>>& atomStates
){
    if(!atomStates){
        return nullptr;
    }

    auto states = std::make_shared<std::vector<OrientationClusterAtomState>>(atomStates->size());
    for(std::size_t atomIndex = 0; atomIndex < atomStates->size(); ++atomIndex){
        const PatternDxaAtomState& source = (*atomStates)[atomIndex];
        auto& target = (*states)[atomIndex];
        target.orientation = source.orientation;
        target.valid = source.orientationValid;
        target.preferredSymmetry = source.symmetryPermutation;
    }
    return states;
}

bool orthonormalizeOrientationMatrix(const Matrix3& input, Matrix3& output){
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

std::shared_ptr<const std::vector<PatternDxaAtomState>> smoothOrientationFallbackStates(
    const PatternCatalog& catalog,
    const std::shared_ptr<const std::vector<PatternDxaAtomState>>& atomStates
){
    if(!atomStates){
        return atomStates;
    }

    constexpr int kSmoothingIterations = 2;
    auto current = std::make_shared<std::vector<PatternDxaAtomState>>(*atomStates);
    auto next = std::make_shared<std::vector<PatternDxaAtomState>>(*atomStates);

    for(int iteration = 0; iteration < kSmoothingIterations; ++iteration){
        *next = *current;

        for(std::size_t atomIndex = 0; atomIndex < current->size(); ++atomIndex){
            const PatternDxaAtomState& state = (*current)[atomIndex];
            if(!state.isMatched() || !state.orientationValid || state.patternId < 0){
                continue;
            }

            const CompiledPattern& pattern = catalog.patternById(state.patternId);
            if(state.localMatcherIndex < 0 ||
               state.localMatcherIndex >= static_cast<int>(pattern.localMatchers.size())){
                continue;
            }

            const CompiledPatternLocalMatcher& matcher =
                pattern.localMatchers[static_cast<std::size_t>(state.localMatcherIndex)];
            if(!matcher.requiresOrientationFallback || matcher.symmetries.empty()){
                continue;
            }

            Matrix3 accumulated = state.orientation;
            int contributionCount = 1;

            const int neighborCount = std::min(state.coordinationNumber, static_cast<int>(MAX_NEIGHBORS));
            for(int neighborSlot = 0; neighborSlot < neighborCount; ++neighborSlot){
                const int neighborAtomIndex = state.neighborAtomIndices[static_cast<std::size_t>(neighborSlot)];
                if(neighborAtomIndex < 0 ||
                   neighborAtomIndex >= static_cast<int>(current->size()) ||
                   neighborAtomIndex == static_cast<int>(atomIndex)){
                    continue;
                }

                const PatternDxaAtomState& neighborState = (*current)[static_cast<std::size_t>(neighborAtomIndex)];
                if(!neighborState.isMatched() ||
                   !neighborState.orientationValid ||
                   neighborState.structureType != state.structureType){
                    continue;
                }

                const Matrix3 relative = Matrix3(neighborState.orientation.transposed() * state.orientation);
                const int symmetryIndex = AnalysisSymmetryUtils::findClosestSymmetryPermutation(
                    matcher.symmetries,
                    relative
                );
                const Matrix3 alignedNeighborOrientation = Matrix3(
                    neighborState.orientation *
                    matcher.symmetries[static_cast<std::size_t>(symmetryIndex)].transformation
                );

                accumulated = Matrix3(accumulated + alignedNeighborOrientation);
                ++contributionCount;
            }

            if(contributionCount <= 1){
                continue;
            }

            Matrix3 smoothedOrientation;
            if(!orthonormalizeOrientationMatrix(accumulated, smoothedOrientation)){
                continue;
            }

            (*next)[atomIndex].orientation = smoothedOrientation;
            (*next)[atomIndex].orientationValid = true;
        }

        std::swap(current, next);
    }

    return current;
}

void recomputeClusterOrientations(StructureAnalysis& analysis, AnalysisContext& context){
    for(Cluster* cluster : analysis.clusterGraph().clusters()){
        if(!cluster || cluster->id == 0 || cluster->structure == LATTICE_OTHER){
            continue;
        }

        Matrix3 orientationV = Matrix3::Zero();
        Matrix3 orientationW = Matrix3::Zero();
        int vectorCount = 0;

        for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
            if(context.atomClusters->getInt(atomIndex) != cluster->id){
                continue;
            }

            const int symmetryPermutation = context.atomSymmetryPermutations->getInt(atomIndex);
            if(symmetryPermutation < 0){
                continue;
            }

            const int neighborCount = context.neighborCounts->getInt(atomIndex);
            for(int neighborSlot = 0; neighborSlot < neighborCount; ++neighborSlot){
                const int neighborAtomIndex = analysis.getNeighbor(static_cast<int>(atomIndex), neighborSlot);
                if(neighborAtomIndex < 0){
                    continue;
                }

                const int latticeVectorIndex = analysis.symmetryPermutationEntry(
                    cluster->structure,
                    symmetryPermutation,
                    neighborSlot
                );
                const Vector3& idealVector = analysis.latticeVector(cluster->structure, latticeVectorIndex);
                const Vector3 spatialVector = context.simCell.wrapVector(
                    context.positions->getPoint3(static_cast<std::size_t>(neighborAtomIndex)) -
                    context.positions->getPoint3(atomIndex)
                );

                for(int i = 0; i < 3; ++i){
                    for(int j = 0; j < 3; ++j){
                        orientationV(i, j) += idealVector[j] * idealVector[i];
                        orientationW(i, j) += idealVector[j] * spatialVector[i];
                    }
                }
                ++vectorCount;
            }
        }

        if(vectorCount < 3){
            continue;
        }

        Matrix3 orientationVInverse;
        if(!orientationV.inverse(orientationVInverse)){
            continue;
        }

        cluster->orientation = Matrix3(orientationW * orientationVInverse);
    }
}

PatternStructureMatchingService::PatternStructureMatchingService() = default;

void PatternStructureMatchingService::setLatticeDirectory(std::string latticeDirectory){
    _latticeDirectory = std::move(latticeDirectory);
}

void PatternStructureMatchingService::setReferenceLatticeDirectory(std::string referenceLatticeDirectory){
    _referenceLatticeDirectory = std::move(referenceLatticeDirectory);
}

void PatternStructureMatchingService::setSelectedPatterns(std::string selectedPatternsCsv){
    _selectedPatternsCsv = std::move(selectedPatternsCsv);
}

void PatternStructureMatchingService::setDissolveSmallClusters(bool dissolveSmallClusters){
    _dissolveSmallClusters = dissolveSmallClusters;
}

json PatternStructureMatchingService::compute(
    const LammpsParser::Frame& frame,
    const std::string& outputBase,
    const std::string& inputDumpPath
){
    if(_latticeDirectory.empty()){
        return AnalysisResult::failure("Lattice directory is required. Use --lattice-dir <path>.");
    }
    if(_referenceLatticeDirectory.empty()){
        return AnalysisResult::failure(
            "Reference lattice directory is required. Use --reference-lattice-dir <path>."
        );
    }

    std::string frameError;
    auto session = AnalysisPipelineUtils::prepareAnalysisSession(
        frame,
        LATTICE_OTHER,
        &frameError
    );
    if(!session){
        return AnalysisResult::failure(frameError);
    }

    AnalysisContext& context = session->context;
    const std::string annotatedDumpPath = outputBase.empty()
        ? inputDumpPath + ".annotated.dump"
        : outputBase + "_annotated.dump";

    try{
        setCrystalTopologySearchRoot(_referenceLatticeDirectory);
        const auto catalog = PatternCatalog::loadFromDirectory(_latticeDirectory);
        std::vector<int> selectedPatternIds;
        if(_selectedPatternsCsv.empty()){
            selectedPatternIds = catalog->defaultSelection();
        }else{
            if(!resolveSelectedPatternIds(*catalog, _selectedPatternsCsv, selectedPatternIds, frameError)){
                return AnalysisResult::failure(frameError);
            }
        }
        if(selectedPatternIds.empty()){
            return AnalysisResult::failure("No supported lattice patterns available for pattern matching.");
        }

        StructureAnalysis analysis(context);
        std::vector<int> atomPatternIds;
        PatternClassifier classifier(catalog);
        classifier.setAtomTypes(frame.types.empty() ? nullptr : &frame.types);
        classifier.setSelectedPatternIds(selectedPatternIds);
        classifier.classify(analysis);
        selectedPatternIds = classifier.selectedPatternIds();
        atomPatternIds = classifier.atomPatternIds();
        const auto atomDxaStates = classifier.atomDxaStates();

        if(!selectedPatternIds.empty()){
            const CompiledPattern& referencePattern = catalog->patternById(selectedPatternIds.front());
            context.inputCrystalType = static_cast<LatticeStructureType>(referencePattern.structureType);
        }else{
            context.inputCrystalType = LATTICE_OTHER;
        }
        classifier.configureDxaClustering(analysis);
        computeMaximumNeighborDistanceFromPatterns(analysis);

        PatternClusterInputAdapter clusterInputAdapter;
        clusterInputAdapter.prepare(analysis, context);
        const auto smoothedAtomDxaStates = smoothOrientationFallbackStates(*catalog, atomDxaStates);
        analysis.setClusterRuleProvider(
            std::make_shared<OrientationClusterRuleProvider>(buildOrientationClusterStates(smoothedAtomDxaStates))
        );
        ClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);
        normalizeReconstructedClusterGraphForExport(analysis, context);
        recomputeClusterOrientations(analysis, context);
        ReconstructedStateCanonicalizer::canonicalizeNeighborShellsToExportConvention(analysis, context);

        json result = AnalysisResult::success();
        result["lattice_directory"] = catalog->latticeDirectory().string();
        result["reference_lattice_directory"] = _referenceLatticeDirectory;
        result["selected_patterns"] = buildSelectedPatterns(*catalog, selectedPatternIds);
        result["cluster_mode"] = "cluster-builder";
        result["main_listing"] = buildMainListing(context, *catalog, atomPatternIds);
        result["pattern_listing"] = buildPatternListing(*catalog, atomPatternIds);
        result["per-atom-properties"] = buildPerAtomProperties(frame, context, *catalog, atomPatternIds);
        if(!AnalysisPipelineUtils::appendClusterOutputs(
            frame,
            outputBase,
            annotatedDumpPath,
            context,
            analysis,
            result,
            &frameError
        )){
            return AnalysisResult::failure(frameError);
        }

        if(!outputBase.empty()){
            const std::string summaryPath = outputBase + "_pattern_analysis.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(result, summaryPath, false)){
                return AnalysisResult::failure("Failed to write " + summaryPath);
            }
            result["pattern_analysis"] = summaryPath;
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("Pattern matching analysis failed: ") + error.what());
    }
}

}
