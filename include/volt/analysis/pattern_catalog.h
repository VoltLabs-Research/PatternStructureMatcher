#pragma once

#include <volt/analysis/pattern_runtime_types.h>

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace Volt {

bool loadPatternTemplateSources(
    const std::filesystem::path& latticeDirectory,
    std::vector<PatternTemplateSource>& outSources,
    std::string* errorMessage
);

bool readTemplateStructureFile(
    const std::filesystem::path& filePath,
    TemplateStructureData& outData,
    std::string* errorMessage
);

CompiledPattern compilePattern(
    const PatternTemplateSource& source,
    const TemplateStructureData* templateData
);

class PatternCatalog {
public:
    static std::shared_ptr<const PatternCatalog> loadFromDirectory(const std::string& latticeDirectory);

    const std::filesystem::path& latticeDirectory() const{
        return _latticeDirectory;
    }

    const std::vector<CompiledPattern>& patterns() const{
        return _patterns;
    }

    const CompiledPattern& patternById(int patternId) const;
    const CompiledPattern* findPatternByName(std::string_view patternName) const;

    std::vector<int> defaultSelection() const;

private:
    PatternCatalog(std::filesystem::path latticeDirectory, std::vector<CompiledPattern> patterns);

    std::filesystem::path _latticeDirectory;
    std::vector<CompiledPattern> _patterns;
    std::unordered_map<std::string, int> _nameToId;
};

}
