#include <volt/analysis/pattern_catalog.h>

#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/analysis/pattern_matching_shared.h>

#include <algorithm>

namespace Volt {

namespace{

constexpr double kPatternShellDistanceTolerance = 1e-6;

struct PatternNeighborCandidate {
    Vector3 vector = Vector3::Zero();
    int species = 0;
    double distance = 0.0;
};

struct NeighborShellRange {
    int begin = 0;
    int count = 0;
    double averageDistance = 0.0;
};

bool approximatelySameDistance(double lhs, double rhs){
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= kPatternShellDistanceTolerance * scale;
}

std::vector<NeighborShellRange> buildNeighborShellRanges(
    const std::vector<PatternNeighborCandidate>& neighbors,
    int coordinationNumber
){
    std::vector<NeighborShellRange> shells;
    if(coordinationNumber <= 0){
        return shells;
    }

    int shellBegin = 0;
    double shellDistanceSum = neighbors.front().distance;
    for(int neighborIndex = 1; neighborIndex < coordinationNumber; ++neighborIndex){
        if(approximatelySameDistance(
            neighbors[static_cast<std::size_t>(neighborIndex)].distance,
            neighbors[static_cast<std::size_t>(neighborIndex - 1)].distance
        )){
            shellDistanceSum += neighbors[static_cast<std::size_t>(neighborIndex)].distance;
            continue;
        }

        const int shellCount = neighborIndex - shellBegin;
        shells.push_back({
            shellBegin,
            shellCount,
            shellDistanceSum / static_cast<double>(shellCount)
        });
        shellBegin = neighborIndex;
        shellDistanceSum = neighbors[static_cast<std::size_t>(neighborIndex)].distance;
    }

    const int finalShellCount = coordinationNumber - shellBegin;
    shells.push_back({
        shellBegin,
        finalShellCount,
        shellDistanceSum / static_cast<double>(finalShellCount)
    });
    return shells;
}

NeighborShellRange selectReferenceShell(const std::vector<NeighborShellRange>& shells){
    NeighborShellRange best;
    bool initialized = false;
    for(const NeighborShellRange& shell : shells){
        if(!initialized ||
           shell.count > best.count ||
           (shell.count == best.count && shell.averageDistance < best.averageDistance)){
            best = shell;
            initialized = true;
        }
    }
    return best;
}

Vector3 pointToVector(const Point3& point){
    return Vector3(point.x(), point.y(), point.z());
}

std::vector<PatternNeighborCandidate> gatherPeriodicNeighbors(
    const TemplateStructureData& templateData,
    int centerAtomIndex,
    int imageRadius
){
    std::vector<PatternNeighborCandidate> neighbors;
    neighbors.reserve(
        static_cast<std::size_t>((imageRadius * 2 + 1) * (imageRadius * 2 + 1) * (imageRadius * 2 + 1)) *
        templateData.positions.size()
    );

    const Vector3 center = pointToVector(templateData.positions[static_cast<std::size_t>(centerAtomIndex)]);
    for(int ix = -imageRadius; ix <= imageRadius; ++ix){
        for(int iy = -imageRadius; iy <= imageRadius; ++iy){
            for(int iz = -imageRadius; iz <= imageRadius; ++iz){
                const Vector3 translation =
                    templateData.cell.column(0) * static_cast<double>(ix) +
                    templateData.cell.column(1) * static_cast<double>(iy) +
                    templateData.cell.column(2) * static_cast<double>(iz);

                for(std::size_t atomIndex = 0; atomIndex < templateData.positions.size(); ++atomIndex){
                    if(ix == 0 && iy == 0 && iz == 0 && atomIndex == static_cast<std::size_t>(centerAtomIndex)){
                        continue;
                    }

                    PatternNeighborCandidate candidate;
                    candidate.vector = pointToVector(templateData.positions[atomIndex]) + translation - center;
                    candidate.distance = candidate.vector.length();
                    candidate.species = templateData.species[atomIndex];
                    if(candidate.distance <= EPSILON){
                        continue;
                    }
                    neighbors.push_back(candidate);
                }
            }
        }
    }

    std::sort(neighbors.begin(), neighbors.end(), [](const auto& lhs, const auto& rhs){
        if(std::abs(lhs.distance - rhs.distance) > 1e-8){
            return lhs.distance < rhs.distance;
        }
        if(std::abs(lhs.vector.x() - rhs.vector.x()) > 1e-8){
            return lhs.vector.x() < rhs.vector.x();
        }
        if(std::abs(lhs.vector.y() - rhs.vector.y()) > 1e-8){
            return lhs.vector.y() < rhs.vector.y();
        }
        return lhs.vector.z() < rhs.vector.z();
    });

    return neighbors;
}

int appendOrResolveLocalMatcherIndex(
    std::vector<CompiledPatternLocalMatcher>& localMatchers,
    CompiledPatternLocalMatcher matcher
){
    const auto it = std::find_if(localMatchers.begin(), localMatchers.end(), [&](const auto& existing){
        return equivalentLocalMatchers(existing, matcher);
    });
    if(it != localMatchers.end()){
        return static_cast<int>(std::distance(localMatchers.begin(), it));
    }
    localMatchers.push_back(std::move(matcher));
    return static_cast<int>(localMatchers.size()) - 1;
}

int findIdentitySymmetryIndex(const std::vector<PatternSymmetryPermutation>& symmetries){
    if(symmetries.empty()){
        return -1;
    }
    return AnalysisSymmetryUtils::findClosestSymmetryPermutation(symmetries, Matrix3::Identity());
}

const Vector3& transformedCanonicalNeighborVector(
    const CompiledPatternLocalMatcher& matcher,
    int symmetryIndex,
    int slot
){
    if(slot < 0 || slot >= matcher.coordinationNumber){
        static const Vector3 zero = Vector3::Zero();
        return zero;
    }
    const auto& symmetry = matcher.symmetries[static_cast<std::size_t>(symmetryIndex)];
    const int mappedSlot = symmetry.permutation[static_cast<std::size_t>(slot)];
    if(mappedSlot < 0 || mappedSlot >= matcher.coordinationNumber){
        static const Vector3 zero = Vector3::Zero();
        return zero;
    }
    return matcher.canonicalNeighborVectors[static_cast<std::size_t>(mappedSlot)];
}

int countMatchingLocalOverlap(
    const CompiledPatternLocalMatcher& matcher,
    int currentSymmetry,
    int neighborSlot,
    int neighborSymmetry,
    int reverseSlot
){
    int commonNeighborCount = 0;
    for(int currentSlot = 0; currentSlot < matcher.coordinationNumber; ++currentSlot){
        if(currentSlot == neighborSlot){
            continue;
        }

        const Vector3 expectedFromNeighbor =
            transformedCanonicalNeighborVector(matcher, currentSymmetry, currentSlot) -
            transformedCanonicalNeighborVector(matcher, currentSymmetry, neighborSlot);

        for(int neighborCandidateSlot = 0; neighborCandidateSlot < matcher.coordinationNumber; ++neighborCandidateSlot){
            if(neighborCandidateSlot == reverseSlot){
                continue;
            }
            if(!(expectedFromNeighbor - transformedCanonicalNeighborVector(
                matcher,
                neighborSymmetry,
                neighborCandidateSlot
            )).isZero(kPatternShellDistanceTolerance)){
                continue;
            }

            ++commonNeighborCount;
            break;
        }
    }

    return commonNeighborCount;
}

bool requiresOrientationFallback(const CompiledPatternLocalMatcher& matcher){
    if(matcher.symmetries.empty() || matcher.coordinationNumber <= 0){
        return false;
    }

    const int identitySymmetry = findIdentitySymmetryIndex(matcher.symmetries);
    if(identitySymmetry < 0){
        return false;
    }

    for(int neighborIndex = 0; neighborIndex < matcher.coordinationNumber; ++neighborIndex){
        const Vector3 currentBond = transformedCanonicalNeighborVector(
            matcher,
            identitySymmetry,
            neighborIndex
        );
        bool hasUsableOverlap = false;

        for(int neighborSymmetry = 0; neighborSymmetry < static_cast<int>(matcher.symmetries.size()) && !hasUsableOverlap; ++neighborSymmetry){
            for(int reverseSlot = 0; reverseSlot < matcher.coordinationNumber; ++reverseSlot){
                const Vector3 neighborBond = transformedCanonicalNeighborVector(
                    matcher,
                    neighborSymmetry,
                    reverseSlot
                );
                if(!(currentBond + neighborBond).isZero(kPatternShellDistanceTolerance)){
                    continue;
                }
                if(countMatchingLocalOverlap(
                    matcher,
                    identitySymmetry,
                    neighborIndex,
                    neighborSymmetry,
                    reverseSlot
                ) >= 2){
                    hasUsableOverlap = true;
                    break;
                }
            }
        }

        if(!hasUsableOverlap){
            return true;
        }
    }

    return false;
}

void initializeIdentitySymmetry(
    const std::vector<Vector3>& canonicalNeighborVectors,
    std::vector<PatternSymmetryPermutation>& symmetries
){
    PatternSymmetryPermutation identity;
    for(std::size_t i = 0; i < canonicalNeighborVectors.size() && i < identity.permutation.size(); ++i){
        identity.permutation[i] = static_cast<int>(i);
    }
    identity.inverseProduct = {0};
    symmetries.clear();
    symmetries.push_back(std::move(identity));
}

void retainProperRotations(std::vector<PatternSymmetryPermutation>& symmetries){
    symmetries.erase(
        std::remove_if(symmetries.begin(), symmetries.end(), [](const PatternSymmetryPermutation& symmetry){
            return symmetry.transformation.determinant() <= 0.0;
        }),
        symmetries.end()
    );
}

void buildGenericSymmetries(
    const std::vector<Vector3>& canonicalNeighborVectors,
    std::vector<PatternSymmetryPermutation>& symmetries
){
    symmetries.clear();
    try{
        AnalysisSymmetryUtils::generateSymmetryPermutations(
            canonicalNeighborVectors,
            static_cast<int>(canonicalNeighborVectors.size()),
            canonicalNeighborVectors,
            symmetries
        );
        retainProperRotations(symmetries);
        AnalysisSymmetryUtils::calculateSymmetryProducts(symmetries);
        if(symmetries.empty()){
            initializeIdentitySymmetry(canonicalNeighborVectors, symmetries);
        }
    }catch(...){
        initializeIdentitySymmetry(canonicalNeighborVectors, symmetries);
    }
}

bool compileGenericLocalMatcherForAtom(
    const PatternTemplateSource& source,
    const TemplateStructureData& templateData,
    int atomIndex,
    CompiledPatternLocalMatcher& outMatcher
){
    const auto neighbors = gatherPeriodicNeighbors(templateData, atomIndex, 1);
    if(static_cast<int>(neighbors.size()) < source.coordinationNumber){
        return false;
    }

    const int coordinationNumber = source.coordinationNumber;
    if(coordinationNumber <= 0 || coordinationNumber > MAX_NEIGHBORS){
        return false;
    }

    const double lastInnerRadius = neighbors[static_cast<std::size_t>(coordinationNumber - 1)].distance;
    double firstOuterRadius = lastInnerRadius * 1.1;
    if(coordinationNumber < static_cast<int>(neighbors.size())){
        firstOuterRadius = neighbors[static_cast<std::size_t>(coordinationNumber)].distance;
    }

    if(firstOuterRadius <= lastInnerRadius + 1e-8){
        firstOuterRadius = lastInnerRadius * 1.1;
    }

    CompiledPatternLocalMatcher matcher;
    matcher.kind = PatternLocalMatcherKind::GenericCna;
    matcher.coordinationNumber = coordinationNumber;
    matcher.localCutoff = 0.5 * (lastInnerRadius + firstOuterRadius);
    matcher.centerSpecies = templateData.species[static_cast<std::size_t>(atomIndex)];
    matcher.requiresSpecies = std::any_of(
        templateData.species.begin(),
        templateData.species.end(),
        [&](int species){ return species != matcher.centerSpecies; }
    );
    matcher.canonicalNeighborVectors.reserve(static_cast<std::size_t>(coordinationNumber));
    matcher.neighborSpeciesByCanonicalSlot.reserve(static_cast<std::size_t>(coordinationNumber));
    matcher.cnaSignatures.reserve(static_cast<std::size_t>(coordinationNumber));

    for(int neighborIndex = 0; neighborIndex < coordinationNumber; ++neighborIndex){
        const auto& candidate = neighbors[static_cast<std::size_t>(neighborIndex)];
        matcher.canonicalNeighborVectors.push_back(candidate.vector);
        matcher.neighborSpeciesByCanonicalSlot.push_back(candidate.species);
    }

    const std::vector<NeighborShellRange> shells = buildNeighborShellRanges(neighbors, coordinationNumber);
    const NeighborShellRange referenceShell = selectReferenceShell(shells);
    if(referenceShell.count <= 0 || referenceShell.averageDistance <= EPSILON){
        return false;
    }
    matcher.referenceNeighborOffset = referenceShell.begin;
    matcher.referenceNeighborCount = referenceShell.count;
    matcher.cutoffMultiplier = matcher.localCutoff / referenceShell.averageDistance;
    matcher.extraNeighborRejectIndex = coordinationNumber < static_cast<int>(neighbors.size())
        ? coordinationNumber
        : -1;

    NeighborBondArray neighborArray;
    const double cutoffSquared = matcher.localCutoff * matcher.localCutoff;
    for(int ni1 = 0; ni1 < coordinationNumber; ++ni1){
        neighborArray.setNeighborBond(ni1, ni1, false);
        for(int ni2 = ni1 + 1; ni2 < coordinationNumber; ++ni2){
            const bool bonded = (
                matcher.canonicalNeighborVectors[static_cast<std::size_t>(ni1)] -
                matcher.canonicalNeighborVectors[static_cast<std::size_t>(ni2)]
            ).squaredLength() < cutoffSquared;
            neighborArray.setNeighborBond(ni1, ni2, bonded);
        }
    }

    for(int row = 0; row < 32; ++row){
        matcher.neighborBondRows[static_cast<std::size_t>(row)] = neighborArray.neighborArray[row];
    }

    for(int neighborIndex = 0; neighborIndex < coordinationNumber; ++neighborIndex){
        matcher.cnaSignatures.push_back(
            computePatternCnaSignature(
                neighborArray,
                neighborIndex,
                coordinationNumber
            )
        );
    }
    matcher.sortedCnaSignatures = matcher.cnaSignatures;
    std::sort(matcher.sortedCnaSignatures.begin(), matcher.sortedCnaSignatures.end());
    buildGenericSymmetries(matcher.canonicalNeighborVectors, matcher.symmetries);
    matcher.requiresOrientationFallback = requiresOrientationFallback(matcher);

    outMatcher = std::move(matcher);
    return true;
}

void compileGenericLocalMatchers(
    CompiledPattern& pattern,
    const PatternTemplateSource& source,
    const TemplateStructureData& templateData
){
    pattern.localMatchers.clear();
    for(int atomIndex = 0; atomIndex < static_cast<int>(templateData.positions.size()); ++atomIndex){
        CompiledPatternLocalMatcher matcher;
        if(!compileGenericLocalMatcherForAtom(source, templateData, atomIndex, matcher)){
            continue;
        }
        appendOrResolveLocalMatcherIndex(pattern.localMatchers, std::move(matcher));
    }

    pattern.supportedForLocalMatching = !pattern.localMatchers.empty();
    pattern.requiresAtomTypes = std::any_of(
        pattern.localMatchers.begin(),
        pattern.localMatchers.end(),
        [](const auto& matcher){ return matcher.requiresSpecies; }
    );
}

}

CompiledPattern compilePattern(
    const PatternTemplateSource& source,
    const TemplateStructureData* templateData
){
    CompiledPattern pattern;
    pattern.name = source.name;
    pattern.coordinationNumber = source.coordinationNumber;
    pattern.structureType = static_cast<int>(StructureType::OTHER);
    pattern.supportedForLocalMatching = false;
    pattern.requiresAtomTypes = false;

    if(templateData){
        compileGenericLocalMatchers(pattern, source, *templateData);
    }

    return pattern;
}

}
