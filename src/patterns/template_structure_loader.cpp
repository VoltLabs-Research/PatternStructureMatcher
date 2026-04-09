#include <volt/analysis/pattern_catalog.h>

#include <yaml-cpp/yaml.h>

#include <stdexcept>

namespace Volt {

namespace TemplateStructureLoaderDetail{

Vector3 parseVector3(const YAML::Node& value){
    if(!value || !value.IsSequence() || value.size() != 3){
        throw std::runtime_error("Expected a 3-component vector.");
    }
    return Vector3(
        value[0].as<double>(),
        value[1].as<double>(),
        value[2].as<double>()
    );
}

void setError(std::string* errorMessage, const std::string& message){
    if(errorMessage){
        *errorMessage = message;
    }
}

}

using namespace TemplateStructureLoaderDetail;

bool readTemplateStructureFile(
    const std::filesystem::path& filePath,
    TemplateStructureData& outData,
    std::string* errorMessage
){
    YAML::Node document;
    try{
        document = YAML::LoadFile(filePath.string());
    }catch(const std::exception& error){
        setError(
            errorMessage,
            "Unable to parse template lattice YAML '" + filePath.string() + "': " + error.what()
        );
        return false;
    }

    try{
        if(!document || !document.IsMap()){
            throw std::runtime_error("Template lattice YAML must contain a mapping root.");
        }

        const double scale = document["scale"] ? document["scale"].as<double>() : 1.0;
        const YAML::Node cellNode = document["cell"];
        if(!cellNode || !cellNode.IsSequence() || cellNode.size() != 3){
            throw std::runtime_error("cell must contain exactly three vectors.");
        }

        Matrix3 cell = Matrix3::Zero();
        cell.column(0) = parseVector3(cellNode[0]) * scale;
        cell.column(1) = parseVector3(cellNode[1]) * scale;
        cell.column(2) = parseVector3(cellNode[2]) * scale;

        const std::string coordinateMode = normalizePatternName(
            document["coordinate_mode"] ? document["coordinate_mode"].as<std::string>() : "fractional"
        );
        const bool cartesianCoordinates =
            coordinateMode == "cartesian" || coordinateMode == "cart" || coordinateMode == "c";
        if(!cartesianCoordinates && coordinateMode != "fractional" && coordinateMode != "reduced" && coordinateMode != "direct"){
            throw std::runtime_error("coordinate_mode must be either 'fractional' or 'cartesian'.");
        }

        const YAML::Node basisNode = document["basis"];
        if(!basisNode || !basisNode.IsSequence() || basisNode.size() == 0){
            throw std::runtime_error("basis must contain at least one atom.");
        }

        outData = TemplateStructureData{};
        outData.cell = cell;
        outData.positions.reserve(basisNode.size());
        outData.species.reserve(basisNode.size());

        for(const auto& atomNode : basisNode){
            if(!atomNode || !atomNode.IsMap()){
                throw std::runtime_error("basis entries must be mappings.");
            }

            const int species = atomNode["species"] ? atomNode["species"].as<int>() : 1;
            Vector3 coordinates = parseVector3(atomNode["position"]);
            if(cartesianCoordinates){
                coordinates *= scale;
            }else{
                coordinates =
                    cell.column(0) * coordinates.x() +
                    cell.column(1) * coordinates.y() +
                    cell.column(2) * coordinates.z();
            }

            outData.positions.emplace_back(coordinates.x(), coordinates.y(), coordinates.z());
            outData.species.push_back(species);
        }
    }catch(const std::exception& error){
        setError(errorMessage, std::string(error.what()) + " In '" + filePath.string() + "'.");
        return false;
    }

    return true;
}

}
