#include <volt/plugin/plugin_entry.h>
#include <volt/analysis/pattern_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = PatternStructureMatchingService;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--lattice_dir", "PatternStructureMatching lattice YAMLs directory", "", &S::setLatticeDirectory),
    opt("--reference_lattice_dir", "OpenDXA reference lattice YAMLs directory", "", &S::setReferenceLatticeDirectory),
    opt("--patterns", "Optional lattice filter (e.g. fcc,bcc)", "", &S::setSelectedPatterns),
    opt("--dissolve_small_clusters", "Mark small clusters as OTHER", false, &S::setDissolveSmallClusters),
};

VOLT_SERVICE_PLUGIN("volt-pattern-structure-matching", "Pattern Structure Matching", S, bindings)
