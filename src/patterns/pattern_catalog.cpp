#include <volt/analysis/pattern_catalog.h>

#include <spdlog/spdlog.h>

#include <stdexcept>

namespace Volt {

PatternCatalog::PatternCatalog(
    std::filesystem::path latticeDirectory,
    std::vector<CompiledPattern> patterns
) :
    _latticeDirectory(std::move(latticeDirectory)),
    _patterns(std::move(patterns))
{
    for(const CompiledPattern& pattern : _patterns){
        _nameToId.emplace(normalizePatternName(pattern.name), pattern.id);
    }
}

std::shared_ptr<const PatternCatalog> PatternCatalog::loadFromDirectory(const std::string& latticeDirectory){
    std::vector<PatternTemplateSource> sources;
    std::string errorMessage;
    if(!loadPatternTemplateSources(latticeDirectory, sources, &errorMessage)){
        throw std::runtime_error(errorMessage);
    }

    std::vector<CompiledPattern> compiledPatterns;
    compiledPatterns.reserve(sources.size());

    for(const PatternTemplateSource& source : sources){
        TemplateStructureData templateData;
        const TemplateStructureData* templatePtr = nullptr;

        if(!std::filesystem::exists(source.definitionPath)){
            spdlog::warn(
                "Lattice '{}' is missing its definition file '{}'",
                source.name,
                source.definitionPath.string()
            );
        }else if(readTemplateStructureFile(source.definitionPath, templateData, &errorMessage)){
            templatePtr = &templateData;
        }else{
            spdlog::warn("Skipping template payload for lattice '{}': {}", source.name, errorMessage);
        }

        CompiledPattern compiled = compilePattern(source, templatePtr);
        compiled.id = static_cast<int>(compiledPatterns.size());
        compiled.structureType = syntheticStructureTypeForPattern(compiled.id);
        compiledPatterns.push_back(std::move(compiled));
    }

    return std::shared_ptr<const PatternCatalog>(
        new PatternCatalog(std::filesystem::path(latticeDirectory), std::move(compiledPatterns))
    );
}

const CompiledPattern& PatternCatalog::patternById(int patternId) const{
    if(patternId < 0 || patternId >= static_cast<int>(_patterns.size())){
        throw std::out_of_range("Pattern id out of range");
    }
    return _patterns[static_cast<std::size_t>(patternId)];
}

const CompiledPattern* PatternCatalog::findPatternByName(std::string_view patternName) const{
    const auto it = _nameToId.find(normalizePatternName(patternName));
    if(it == _nameToId.end()){
        return nullptr;
    }
    return &_patterns[static_cast<std::size_t>(it->second)];
}

std::vector<int> PatternCatalog::defaultSelection() const{
    std::vector<int> selection;
    for(const CompiledPattern& pattern : _patterns){
        if(!pattern.supportedForLocalMatching){
            continue;
        }
        selection.push_back(pattern.id);
    }
    return selection;
}

}
