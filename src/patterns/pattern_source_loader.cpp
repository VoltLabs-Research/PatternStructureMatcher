#include <volt/analysis/pattern_catalog.h>

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <system_error>
#include <stdexcept>

namespace Volt {

namespace PatternSourceLoaderDetail{

void setError(std::string* errorMessage, const std::string& message){
    if(errorMessage){
        *errorMessage = message;
    }
}

PatternTemplateSource parsePatternSourceFile(const std::filesystem::path& filePath){
    YAML::Node document;
    try{
        document = YAML::LoadFile(filePath.string());
    }catch(const std::exception& error){
        throw std::runtime_error(
            "Unable to parse lattice YAML '" + filePath.string() + "': " + error.what()
        );
    }

    if(!document || !document.IsMap()){
        throw std::runtime_error("Lattice YAML must contain a mapping root.");
    }

    PatternTemplateSource source;
    source.definitionPath = filePath;
    source.name = document["name"].as<std::string>();
    source.coordinationNumber = document["coordination_number"].as<int>();

    if(source.name.empty()){
        throw std::runtime_error("Lattice name cannot be empty.");
    }
    if(source.coordinationNumber <= 0 || source.coordinationNumber > MAX_NEIGHBORS){
        throw std::runtime_error("coordination_number must be between 1 and MAX_NEIGHBORS.");
    }

    return source;
}

}

using namespace PatternSourceLoaderDetail;

bool loadPatternTemplateSources(
    const std::filesystem::path& latticeDirectory,
    std::vector<PatternTemplateSource>& outSources,
    std::string* errorMessage
){
    std::error_code error;
    if(latticeDirectory.empty()){
        setError(errorMessage, "Lattice directory cannot be empty.");
        return false;
    }
    if(!std::filesystem::exists(latticeDirectory, error) || !std::filesystem::is_directory(latticeDirectory, error)){
        setError(
            errorMessage,
            "Lattice directory '" + latticeDirectory.string() + "' does not exist or is not a directory."
        );
        return false;
    }

    std::vector<std::filesystem::path> latticeFiles;
    for(const auto& entry : std::filesystem::recursive_directory_iterator(latticeDirectory, error)){
        if(error){
            break;
        }
        if(!entry.is_regular_file()){
            continue;
        }

        const std::string extension = normalizePatternName(entry.path().extension().string());
        if(extension != ".yml" && extension != ".yaml"){
            continue;
        }
        latticeFiles.push_back(entry.path().lexically_normal());
    }
    std::sort(latticeFiles.begin(), latticeFiles.end());

    if(error){
        setError(
            errorMessage,
            "Failed to enumerate lattice YAML files under '" + latticeDirectory.string() + "'."
        );
        return false;
    }
    if(latticeFiles.empty()){
        setError(
            errorMessage,
            "No lattice YAML files were found under '" + latticeDirectory.string() + "'."
        );
        return false;
    }

    outSources.clear();
    outSources.reserve(latticeFiles.size());

    try{
        for(const auto& latticePath : latticeFiles){
            outSources.push_back(parsePatternSourceFile(latticePath));
        }
    }catch(const std::exception& error){
        setError(errorMessage, error.what());
        outSources.clear();
        return false;
    }

    return true;
}

}
