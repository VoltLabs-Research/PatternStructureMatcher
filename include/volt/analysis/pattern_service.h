#pragma once

#include <volt/core/lammps_parser.h>

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace Volt {

using json = nlohmann::json;

class PatternStructureMatchingService {
public:
    PatternStructureMatchingService();

    void setLatticeDirectory(std::string latticeDirectory);
    void setReferenceLatticeDirectory(std::string referenceLatticeDirectory);
    void setSelectedPatterns(std::string selectedPatternsCsv);
    void setDissolveSmallClusters(bool dissolveSmallClusters);

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase,
        const std::string& inputDumpPath
    );

private:
    std::string _latticeDirectory;
    std::string _referenceLatticeDirectory;
    std::string _selectedPatternsCsv;
    bool _dissolveSmallClusters = false;
};

}
