cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(pattern-structure-matching VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)

find_package(coretoolkit REQUIRED)
find_package(common-neighbor-analysis REQUIRED)
find_package(structure-identification REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(yaml-cpp REQUIRED)

set(PATTERN_STRUCTURE_MATCHING_BOOST_TARGET "")
if(TARGET Boost::headers)
    set(PATTERN_STRUCTURE_MATCHING_BOOST_TARGET Boost::headers)
elseif(TARGET Boost::boost)
    set(PATTERN_STRUCTURE_MATCHING_BOOST_TARGET Boost::boost)
elseif(TARGET boost::headers)
    set(PATTERN_STRUCTURE_MATCHING_BOOST_TARGET boost::headers)
elseif(TARGET boost::boost)
    set(PATTERN_STRUCTURE_MATCHING_BOOST_TARGET boost::boost)
else()
    message(FATAL_ERROR "PatternStructureMatching requires a Boost headers target")
endif()

set(LIB_SOURCES
    ${CMAKE_SOURCE_DIR}/src/analysis/pattern_service.cpp
    ${CMAKE_SOURCE_DIR}/src/analysis/pattern_structure_analysis.cpp
    ${CMAKE_SOURCE_DIR}/src/analysis/pattern_cluster_input_adapter.cpp
    ${CMAKE_SOURCE_DIR}/src/analysis/pattern_classifier.cpp
    ${CMAKE_SOURCE_DIR}/src/analysis/pattern_matching_shared.cpp
    ${CMAKE_SOURCE_DIR}/src/analysis/pattern_dxa_topology_provider.cpp
    ${CMAKE_SOURCE_DIR}/src/patterns/pattern_catalog.cpp
    ${CMAKE_SOURCE_DIR}/src/patterns/pattern_source_loader.cpp
    ${CMAKE_SOURCE_DIR}/src/patterns/template_structure_loader.cpp
    ${CMAKE_SOURCE_DIR}/src/patterns/pattern_compiler.cpp
)

add_library(${PROJECT_NAME}_lib STATIC ${LIB_SOURCES})
set_target_properties(${PROJECT_NAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

get_target_property(
    PATTERN_STRUCTURE_MATCHING_CORETOOLKIT_INCLUDE_DIRS
    coretoolkit::coretoolkit
    INTERFACE_INCLUDE_DIRECTORIES
)
get_target_property(
    PATTERN_STRUCTURE_MATCHING_COMMON_NEIGHBOR_ANALYSIS_INCLUDE_DIRS
    common-neighbor-analysis::common-neighbor-analysis
    INTERFACE_INCLUDE_DIRECTORIES
)
get_target_property(
    PATTERN_STRUCTURE_MATCHING_STRUCTURE_IDENTIFICATION_INCLUDE_DIRS
    structure-identification::structure-identification
    INTERFACE_INCLUDE_DIRECTORIES
)

target_include_directories(${PROJECT_NAME}_lib BEFORE PUBLIC
    ${PATTERN_STRUCTURE_MATCHING_CORETOOLKIT_INCLUDE_DIRS}
    ${PATTERN_STRUCTURE_MATCHING_COMMON_NEIGHBOR_ANALYSIS_INCLUDE_DIRS}
    ${PATTERN_STRUCTURE_MATCHING_STRUCTURE_IDENTIFICATION_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
)
target_precompile_headers(${PROJECT_NAME}_lib PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC
    ${PATTERN_STRUCTURE_MATCHING_BOOST_TARGET}
    TBB::tbb
    common-neighbor-analysis::common-neighbor-analysis
    coretoolkit::coretoolkit
    structure-identification::structure-identification
    spdlog::spdlog
    nlohmann_json::nlohmann_json
    yaml-cpp::yaml-cpp
)

add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME} PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME} PRIVATE
    TBB::tbb
    ${PROJECT_NAME}_lib
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(TARGETS ${PROJECT_NAME}_lib DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/lattices/ DESTINATION lattices)
