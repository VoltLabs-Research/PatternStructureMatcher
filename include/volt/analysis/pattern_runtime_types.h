#pragma once

#include <volt/core/volt.h>
#include <volt/math/point3.h>
#include <volt/structures/crystal_structure_types.h>

#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace Volt {

struct PatternTemplateSource {
    std::string name;
    std::filesystem::path definitionPath;
    int coordinationNumber = 0;
};

struct TemplateStructureData {
    Matrix3 cell = Matrix3::Zero();
    std::vector<Point3> positions;
    std::vector<int> species;
};

struct PatternSymmetryPermutation {
    Matrix3 transformation = Matrix3::Identity();
    std::array<int, MAX_NEIGHBORS> permutation{};
    std::vector<int> inverseProduct;
};

struct PatternCnaSignature {
    int numCommonNeighbors = 0;
    int numBonds = 0;
    int maxChainLength = 0;

    auto tie() const{
        return std::tie(numCommonNeighbors, numBonds, maxChainLength);
    }

    bool operator==(const PatternCnaSignature& other) const{
        return tie() == other.tie();
    }

    bool operator<(const PatternCnaSignature& other) const{
        return tie() < other.tie();
    }
};

enum class PatternLocalMatcherKind {
    GenericCna,
};

struct CompiledPatternLocalMatcher {
    PatternLocalMatcherKind kind = PatternLocalMatcherKind::GenericCna;
    int coordinationNumber = 0;
    double localCutoff = 0.0;
    int referenceNeighborOffset = 0;
    int referenceNeighborCount = 0;
    double cutoffMultiplier = 0.0;
    int extraNeighborRejectIndex = -1;
    int centerSpecies = 0;
    bool requiresSpecies = false;
    std::vector<Vector3> canonicalNeighborVectors;
    std::vector<int> neighborSpeciesByCanonicalSlot;
    std::array<unsigned int, 32> neighborBondRows{};
    std::vector<PatternCnaSignature> cnaSignatures;
    std::vector<PatternCnaSignature> sortedCnaSignatures;
    std::vector<PatternSymmetryPermutation> symmetries;
    bool requiresOrientationFallback = false;
};

struct CompiledPattern {
    int id = -1;
    std::string name;
    int coordinationNumber = 0;
    bool supportedForLocalMatching = false;
    bool requiresAtomTypes = false;
    int structureType = static_cast<int>(StructureType::OTHER);
    std::vector<CompiledPatternLocalMatcher> localMatchers;
};

struct PatternAtomMatch {
    PatternAtomMatch(){
        orderedNeighborIndices.fill(-1);
        idealNeighborVectors.fill(Vector3::Zero());
    }

    int patternId = -1;
    int structureType = static_cast<int>(StructureType::OTHER);
    int localMatcherIndex = -1;
    double localCutoff = 0.0;
    int coordinationNumber = 0;
    std::uint64_t allowedSymmetryMask = 0;
    int symmetryPermutation = -1;
    Matrix3 orientation = Matrix3::Identity();
    bool orientationValid = false;
    std::array<int, MAX_NEIGHBORS> orderedNeighborIndices;
    std::array<Vector3, MAX_NEIGHBORS> idealNeighborVectors;
};

struct PatternDxaAtomState {
    PatternDxaAtomState(){
        neighborAtomIndices.fill(-1);
        idealNeighborVectors.fill(Vector3::Zero());
    }

    int patternId = -1;
    int structureType = static_cast<int>(StructureType::OTHER);
    int localMatcherIndex = -1;
    int coordinationNumber = 0;
    std::uint64_t allowedSymmetryMask = 0;
    int symmetryPermutation = -1;
    Matrix3 orientation = Matrix3::Identity();
    bool orientationValid = false;
    std::array<int, MAX_NEIGHBORS> neighborAtomIndices;
    std::array<Vector3, MAX_NEIGHBORS> idealNeighborVectors;

    bool isMatched() const{
        return patternId >= 0 && coordinationNumber > 0;
    }
};

inline std::string normalizePatternName(std::string_view text){
    std::string normalized;
    normalized.reserve(text.size());
    for(char character : text){
        normalized.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(character))));
    }
    return normalized;
}

inline const char* structureTypeNameForExport(int structureType){
    if(structureType == static_cast<int>(StructureType::OTHER)){
        return "OTHER";
    }
    if(structureType >= 1000){
        return "PATTERN";
    }
    return "STRUCTURE";
}

inline int syntheticStructureTypeForPattern(int patternId){
    return 1000 + patternId;
}

inline bool isSyntheticPatternStructureType(int structureType){
    return structureType >= 1000;
}

inline int patternIdFromSyntheticStructureType(int structureType){
    return isSyntheticPatternStructureType(structureType) ? (structureType - 1000) : -1;
}

}
