#include <volt/analysis/pattern_catalog.h>
#include <volt/analysis/pattern_service.h>
#include <volt/cli/common.h>

using namespace Volt;
using namespace Volt::CLI;

namespace Volt {

void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Pattern Structure Matching");
    std::cerr
        << "  --lattice-dir <path>         Directory containing PatternStructureMatching lattice YAMLs.\n"
        << "  --reference-lattice-dir <path> Directory containing OpenDXA reference lattice YAMLs.\n"
        << "  --patterns <csv>             Optional lattice filter (e.g. fcc,bcc). Default: all.\n"
        << "  --dissolveSmallClusters      Mark small clusters as OTHER after building clusters.\n";
    printHelpOption();
}

}

int main(int argc, char* argv[]){
    std::string filename;
    std::string outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    if(const int startupStatus = handleHelpOrMissingInput(argc, argv, opts, filename, Volt::showUsage);
       startupStatus >= 0){
        return startupStatus;
    }

    initLogging("volt-pattern-structure-matching");

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)){
        return 1;
    }

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    const std::string latticeDirectory = getString(opts, "--lattice-dir", "");
    if(latticeDirectory.empty()){
        spdlog::error("No lattice directory resolved. Use --lattice-dir <path>.");
        return 1;
    }
    const std::string referenceLatticeDirectory = getString(opts, "--reference-lattice-dir", "");
    if(referenceLatticeDirectory.empty()){
        spdlog::error("No reference lattice directory resolved. Use --reference-lattice-dir <path>.");
        return 1;
    }

    const std::string selectedPatterns = getString(opts, "--patterns", "");

    PatternStructureMatchingService analyzer;
    analyzer.setLatticeDirectory(latticeDirectory);
    analyzer.setReferenceLatticeDirectory(referenceLatticeDirectory);
    analyzer.setSelectedPatterns(selectedPatterns);
    analyzer.setDissolveSmallClusters(hasOption(opts, "--dissolveSmallClusters"));

    spdlog::info("Using lattice directory: {}", latticeDirectory);
    if(!referenceLatticeDirectory.empty()){
        spdlog::info("Using reference lattice directory: {}", referenceLatticeDirectory);
    }
    spdlog::info("Starting pattern structure matching...");

    json result = analyzer.compute(frame, outputBase, filename);
    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("Pattern structure matching completed.");
    return 0;
}
