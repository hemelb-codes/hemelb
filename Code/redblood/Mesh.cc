#include <fstream>
#include <assert.h>
#include <numeric>

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
        for(unsigned j(0); j < 3; ++j)
            result->facets[i].insert(indices[j] - offset);
        log::Logger::Log<log::Trace, log::Singleton>(
            "Facet %i with %i, %i, %i", i, indices[0] - offset,
            indices[1] - offset, indices[2] - offset
        );
    }

    log::Logger::Log<log::Debug, log::Singleton>(
            "Read %i vertices and %i triangular facets", num_vertices,
            num_facets);

    return result;
}

t_Real3D Mesh::barycenter() const {
    typedef MeshData::t_Vertices::value_type t_Vertex;
    return std::accumulate(
        mesh_->vertices.begin(), mesh_->vertices.end(),
        t_Vertex(0, 0, 0)
    ) / t_Vertex::value_type(mesh_->vertices.size());
}

namespace {
    bool edge_sharing(std::set<unsigned int> const &_a, 
        std::set<unsigned int> const &_b ) {
        unsigned int result(0);
        std::set<unsigned int> :: const_iterator i_neigh = _a.begin();
        for(; i_neigh != _a.end(); ++i_neigh)
            if(_b.count(*i_neigh) == 1) ++result;
        return result == 2;
    }
    // Adds value as first non-negative number, if value not in array yet
    void insert(boost::array<unsigned, 3> &_container,
        unsigned _value, unsigned _max) {
        for(unsigned i(0); i < _container.size(); ++i)
            if(_container[i] == _max) { _container[i] = _value; return; }
            else if(_container[i] == _value) return;
    }
}

MeshTopology::MeshTopology(MeshData const &_mesh) {
    vertex_to_facets.resize(_mesh.vertices.size());
    facet_neighbors.resize(_mesh.facets.size());

    // Loop over facets to create map from vertices to facets
    MeshData::t_Facets::const_iterator i_facet = _mesh.facets.begin();
    MeshData::t_Facets::const_iterator const i_facet_end = _mesh.facets.end();
    for(unsigned int i(0); i_facet != i_facet_end; ++i_facet, ++i) {
        t_Indices::const_iterator i_vertex = i_facet->begin();
        for(; i_vertex != i_facet->end(); ++i_vertex)
            vertex_to_facets.at(*i_vertex).insert(i);
    }

    // Now creates map of neighboring facets
    unsigned int const Nmax = _mesh.facets.size();
    boost::array<unsigned int, 3> const neg = {{ Nmax, Nmax, Nmax }};
    for(unsigned int i(0); i < facet_neighbors.size(); ++i)
        facet_neighbors[i] = neg;
    i_facet = _mesh.facets.begin();
    for(unsigned int i(0); i_facet != i_facet_end; ++i_facet, ++i) {
        t_Indices::const_iterator i_vertex = i_facet->begin();
        for(; i_vertex != i_facet->end(); ++i_vertex) {
            // check facets that this node is attached to
            std::set<unsigned int> const &facets
                = vertex_to_facets.at(*i_vertex);
            std::set<unsigned int> :: const_iterator i_neigh = facets.begin();
            for(; i_neigh != facets.end(); ++i_neigh) {
                if(edge_sharing(*i_facet, _mesh.facets.at(*i_neigh)))
                  insert(facet_neighbors.at(i), *i_neigh, Nmax);
            }
        }
    }
#   ifndef NDEBUG
    // Checks there are no uninitialized values
    for(unsigned int i(0); i < facet_neighbors.size(); ++i)
        for(unsigned int j(0); j < 3; ++j)
            assert(facet_neighbors[i][j] < Nmax);
#   endif
}

}} // hemelb::rbc
