#pragma once

#include <volt/analysis/pattern_catalog.h>
#include <volt/analysis/structure_analysis.h>

#include <memory>
#include <string>
#include <vector>

namespace Volt {

class PatternClassifier {
public:
    explicit PatternClassifier(std::shared_ptr<const PatternCatalog> catalog);

    void setAtomTypes(const std::vector<int>* atomTypes);
    void setSelectedPatternIds(std::vector<int> selectedPatternIds);
    void classify(StructureAnalysis& analysis);
    void configureDxaClustering(StructureAnalysis& analysis) const;

    const std::shared_ptr<const PatternCatalog>& catalog() const{
        return _catalog;
    }

    const std::vector<int>& selectedPatternIds() const{
        return _selectedPatternIds;
    }

    const std::vector<int>& atomPatternIds() const{
        return _atomPatternIds;
    }

    std::shared_ptr<const std::vector<PatternDxaAtomState>> atomDxaStates() const{
        return _atomDxaStates;
    }

private:
    std::shared_ptr<const PatternCatalog> _catalog;
    std::vector<int> _selectedPatternIds;
    std::vector<int> _configuredSelectedPatternIds;
    std::vector<int> _atomPatternIds;
    std::shared_ptr<std::vector<PatternDxaAtomState>> _atomDxaStates;
    const std::vector<int>* _atomTypes = nullptr;
};

}
