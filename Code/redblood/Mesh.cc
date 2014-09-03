#include <fstream>
#include <assert.h>

#include "redblood/Mesh.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "Exception.h"

namespace hemelb { namespace redblood {

boost::shared_ptr<redblood::MeshData> read_mesh(std::string const &_filename) {
    log::Logger::Log<log::Debug, log::Singleton>(
            "Reading red blood cell from %s", _filename.c_str());

    // Open file if it exists
    std::string line;
    std::ifstream file;
    if (!util::file_exists(_filename.c_str()))
        throw Exception() << "Red-blood-cell mesh file '"
            << _filename.c_str() << "' does not exist";
    file.open(_filename.c_str());

    // Drop header
    for(int i(0); i < 4; ++i) std::getline(file, line);

    // Number of vertices
    unsigned int num_vertices;
    file >> num_vertices;


    // Create Mesh data
    boost::shared_ptr<MeshData> result(new MeshData);
    result->vertices.resize(num_vertices);

    // Then read in first and subsequent lines
    unsigned int offset;
    file >> offset >> result->vertices[0].x
        >> result->vertices[0].y
        >> result->vertices[0].z;
    log::Logger::Log<log::Trace, log::Singleton>(
        "Vertex 0 at %d, %d, %d", result->vertices[0].x,
        result->vertices[0].y, result->vertices[0].z
    );
    for(unsigned int i(1), index(0); i < num_vertices; ++i) {
        file >> index >> result->vertices[i].x
            >> result->vertices[i].y
            >> result->vertices[i].z;
        log::Logger::Log<log::Trace, log::Singleton>(
            "Vertex %i at %d, %d, %d", i, result->vertices[i].x,
            result->vertices[i].y, result->vertices[i].z
        );
    }

    // Drop mid-file headers
    for(int i(0); i < 3; ++i) std::getline(file, line);

    // Read facet indices
    unsigned int num_facets, indices[3];
    file >> num_facets;
    result->facets.resize(num_facets);
    for(unsigned int i(0), dummy(0); i < num_facets; ++i) {
        file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
            >> indices[0] >> indices[1] >> indices[2];
        result->facets[i].x = indices[0] - offset;
        result->facets[i].y = indices[1] - offset;
        result->facets[i].z = indices[2] - offset;
        log::Logger::Log<log::Trace, log::Singleton>(
            "Facet %i with %i, %i, %i", i, result->facets[i].x,
            result->facets[i].y, result->facets[i].z
        );
    }

    log::Logger::Log<log::Debug, log::Singleton>(
            "Read %i vertices and %i triangular facets", num_vertices,
            num_facets);

    return result;
}

}} // hemelb::rbc
