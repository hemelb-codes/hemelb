//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <fstream>
#include <cassert>
#include <numeric>
#include <map>

#include "redblood/Mesh.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "Exception.h"
#include "constants.h"

namespace hemelb
{
  namespace redblood
  {
    std::shared_ptr<MeshData> read_mesh(std::string const &filename)
    {
      log::Logger::Log<log::Debug, log::Singleton>("Reading red blood cell from %s",
                                                   filename.c_str());

      // Open file if it exists
      std::ifstream file;

      if (!util::file_exists(filename.c_str()))
        throw Exception() << "Red-blood-cell mesh file '" << filename.c_str() << "' does not exist";

      file.open(filename.c_str());
      return read_mesh(file);
    }

    std::shared_ptr<MeshData> read_mesh(std::istream &stream)
    {
      log::Logger::Log<log::Debug, log::Singleton>("Reading red blood cell from stream");

      std::string line;

      // Drop header
      for (int i(0); i < 4; ++i)
      {
        std::getline(stream, line);
      }

      // Number of vertices
      unsigned int num_vertices;
      stream >> num_vertices;

      // Create Mesh data
      std::shared_ptr<MeshData> result(new MeshData);
      result->vertices.resize(num_vertices);

      // Then read in first and subsequent lines
      MeshData::Facet::value_type offset;
      stream >> offset >> result->vertices[0].x >> result->vertices[0].y >> result->vertices[0].z;
      log::Logger::Log<log::Trace, log::Singleton>("Vertex 0 at %d, %d, %d",
                                                   result->vertices[0].x,
                                                   result->vertices[0].y,
                                                   result->vertices[0].z);

      for (unsigned int i(1), index(0); i < num_vertices; ++i)
      {
        stream >> index >> result->vertices[i].x >> result->vertices[i].y >> result->vertices[i].z;
        log::Logger::Log<log::Trace, log::Singleton>("Vertex %i at %d, %d, %d",
                                                     i,
                                                     result->vertices[i].x,
                                                     result->vertices[i].y,
                                                     result->vertices[i].z);
      }

      // Drop mid-file headers
      for (int i(0); i < 3; ++i)
      {
        std::getline(stream, line);
      }

      // Read facet indices
      unsigned int num_facets;
      MeshData::Facet indices;
      stream >> num_facets;
      result->facets.resize(num_facets);

      for (unsigned int i(0), dummy(0); i < num_facets; ++i)
      {
        stream >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> indices[0] >> indices[1]
            >> indices[2];
        result->facets[i][0] = indices[0] - offset;
        result->facets[i][1] = indices[1] - offset;
        result->facets[i][2] = indices[2] - offset;
        log::Logger::Log<log::Trace, log::Singleton>("Facet %i with %i, %i, %i",
                                                     i,
                                                     indices[0] - offset,
                                                     indices[1] - offset,
                                                     indices[2] - offset);
      }

      log::Logger::Log<log::Debug, log::Singleton>("Read %i vertices and %i triangular facets",
                                                   num_vertices,
                                                   num_facets);

      return result;
    }

    void write_mesh(std::string const &filename, MeshData const &data)
    {
      log::Logger::Log<log::Debug, log::Singleton>("Writing red blood cell from %s",
                                                   filename.c_str());
      std::ofstream file(filename.c_str());
      write_mesh(file, data);
    }

    void write_vtkmesh(std::string const &filename, MeshData const &data)
    {
      log::Logger::Log<log::Debug, log::Singleton>("Writing red blood cell from %s",
                                                   filename.c_str());
      std::ofstream file(filename.c_str());
      write_vtkmesh(file, data);
    }

    void write_mesh(std::ostream &stream, MeshData const &data)
    {
      // Write Header
      stream << "$MeshFormat\n2 0 8\n$EndMeshFormat\n" << "$Nodes\n" << data.vertices.size()
          << "\n";

      typedef MeshData::Vertices::const_iterator VertexIterator;
      VertexIterator i_vertex = data.vertices.begin();
      VertexIterator const i_vertex_end = data.vertices.end();

      for (unsigned i(1); i_vertex != i_vertex_end; ++i_vertex, ++i)
        stream << i << " " << (*i_vertex)[0] << " " << (*i_vertex)[1] << " " << (*i_vertex)[2]
            << "\n";

      stream << "$EndNode\n" << "$Elements\n" << data.facets.size() << "\n";

      typedef MeshData::Facets::const_iterator FacetIterator;
      FacetIterator i_facet = data.facets.begin();
      FacetIterator const i_facet_end = data.facets.end();

      for (unsigned i(1); i_facet != i_facet_end; ++i_facet, ++i)
        stream << i << " 1 2 3 4 5 " << (*i_facet)[0] + 1 << " " << (*i_facet)[1] + 1 << " "
            << (*i_facet)[2] + 1 << "\n";

      stream << "$EndElement\n";
    }

    void write_vtkmesh(std::ostream &stream, MeshData const &data)
    {
      // Write Header
      stream << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          << "  <PolyData>\n" << "    <Piece NumberOfPoints=\"" << data.vertices.size()
          << "\" NumberOfPolys=\"" << data.facets.size() << "\">\n" << "      <Points>\n"
          << "        <DataArray NumberOfComponents=\"3\" type=\"Float32\">\n";

      typedef MeshData::Vertices::const_iterator VertexIterator;
      VertexIterator i_vertex = data.vertices.begin();
      VertexIterator const i_vertex_end = data.vertices.end();

      for (unsigned i(1); i_vertex != i_vertex_end; ++i_vertex, ++i)
        stream << (*i_vertex)[0] << " " << (*i_vertex)[1] << " " << (*i_vertex)[2] << " ";

      stream << "\n        </DataArray>\n" << "      </Points>\n" << "      <Polys>\n"
          << "        <DataArray type=\"Int32\" Name=\"connectivity\">\n";

      typedef MeshData::Facets::const_iterator FacetIterator;
      FacetIterator i_facet = data.facets.begin();
      FacetIterator const i_facet_end = data.facets.end();

      for (; i_facet != i_facet_end; ++i_facet)
        stream << (*i_facet)[0] << " " << (*i_facet)[1] << " " << (*i_facet)[2] << " ";

      stream << "\n        </DataArray>\n"
          << "        <DataArray type=\"Int32\" Name=\"offsets\">\n";

      for (unsigned i(0); i < data.facets.size(); ++i)
      {
        stream << (i + 1) * 3 << " ";
      }

      stream << "\n        </DataArray>\n" << "      </Polys>\n";
      stream << "    </Piece>\n  </PolyData>\n</VTKFile>\n";
    }

    LatticePosition barycenter(MeshData::Vertices const &vertices)
    {
      typedef MeshData::Vertices::value_type Vertex;
      return std::accumulate(vertices.begin(), vertices.end(), Vertex(0, 0, 0))
          / Vertex::value_type(vertices.size());
    }
    LatticePosition barycenter(MeshData const &mesh)
    {
      return barycenter(mesh.vertices);
    }
    PhysicalVolume volume(MeshData::Vertices const &vertices, MeshData::Facets const &facets)
    {
      MeshData::Facets::const_iterator i_facet = facets.begin();
      MeshData::Facets::const_iterator const i_facet_end(facets.end());
      PhysicalVolume result(0);

      for (; i_facet != i_facet_end; ++i_facet)
      {
        LatticePosition const &v0(vertices[ (*i_facet)[0]]);
        LatticePosition const &v1(vertices[ (*i_facet)[1]]);
        LatticePosition const &v2(vertices[ (*i_facet)[2]]);
        result += v0.Cross(v1).Dot(v2);
      }

      // Minus sign comes from outward facing facet orientation
      return -result / PhysicalVolume(6);
    }
    PhysicalVolume volume(MeshData const &mesh)
    {
      return volume(mesh.vertices, mesh.facets);
    }

    PhysicalVolume area(MeshData::Vertices const &vertices, MeshData::Facets const &facets)
    {
      MeshData::Facets::const_iterator i_facet = facets.begin();
      MeshData::Facets::const_iterator const i_facet_end(facets.end());
      PhysicalVolume result(0);

      for (; i_facet != i_facet_end; ++i_facet)
      {
        LatticePosition const &v0(vertices[ (*i_facet)[0]]);
        LatticePosition const &v1(vertices[ (*i_facet)[1]]);
        LatticePosition const &v2(vertices[ (*i_facet)[2]]);
        result += (v0 - v1).Cross(v2 - v1).GetMagnitude();
      }

      return result * 0.5;
    }
    PhysicalArea area(MeshData const &mesh)
    {
      return area(mesh.vertices, mesh.facets);
    }

    namespace
    {
      bool contains(MeshData::Facet const &a, MeshData::Facet::value_type v)
      {
        return a[0] == v or a[1] == v or a[2] == v;
      }
      bool edge_sharing(MeshData::Facet const &a, MeshData::Facet const &b)
      {
        return (contains(b, a[0]) ?
          1 :
          0) + (contains(b, a[1]) ?
          1 :
          0) + (contains(b, a[2]) ?
          1 :
          0) >= 2;
      }
      // Adds value as first non-negative number, if value not in array yet
      void insert(MeshData::Facet &container, MeshData::Facet::value_type value,
                  MeshData::Facet::value_type max)
      {
        for (size_t i(0); i < container.size(); ++i)
          if (container[i] >= max)
          {
            container[i] = value;
            return;
          }
          else if (container[i] == value)
          {
            return;
          }
      }
    }

    MeshTopology::MeshTopology(MeshData const &mesh)
    {
      vertexToFacets.resize(mesh.vertices.size());
      facetNeighbors.resize(mesh.facets.size());

      // Loop over facets to create map from vertices to facets
      MeshData::Facets::const_iterator i_facet = mesh.facets.begin();
      MeshData::Facets::const_iterator const i_facet_end = mesh.facets.end();

      for (unsigned int i(0); i_facet != i_facet_end; ++i_facet, ++i)
      {
        vertexToFacets.at( (*i_facet)[0]).insert(i);
        vertexToFacets.at( (*i_facet)[1]).insert(i);
        vertexToFacets.at( (*i_facet)[2]).insert(i);
      }

      // Now creates map of neighboring facets
      size_t const Nmax = mesh.facets.size();
      FacetNeighbors::value_type intelcompiler;
      // intel compiler without parameter list initialization
      intelcompiler[0] = Nmax; intelcompiler[1] = Nmax; intelcompiler[2] = Nmax; 
      std::fill(facetNeighbors.begin(), facetNeighbors.end(), intelcompiler);
      i_facet = mesh.facets.begin();

      for (unsigned int i(0); i_facet != i_facet_end; ++i_facet, ++i)
      {
        for (size_t node(0); node != i_facet->size(); ++node)
        {
          // check facets that this node is attached to
          MeshTopology::VertexToFacets::const_reference neighboringFacets =
              vertexToFacets.at( (*i_facet)[node]);
          MeshTopology::VertexToFacets::value_type::const_iterator i_neigh =
              neighboringFacets.begin();

          for (; i_neigh != neighboringFacets.end(); ++i_neigh)
          {
            if (i == *i_neigh)
            {
              continue;
            }

            if (edge_sharing(*i_facet, mesh.facets.at(*i_neigh)))
            {
              insert(facetNeighbors.at(i), *i_neigh, Nmax);
            }
          }
        }
      }

#ifndef NDEBUG

      // Checks there are no uninitialized values
      for (unsigned int i(0); i < facetNeighbors.size(); ++i)
        for (unsigned int j(0); j < 3; ++j)
        {
          assert(facetNeighbors[i][j] < Nmax);
        }

#endif
    }

    void Mesh::operator*=(Dimensionless const &scale)
    {
      auto const barycenter = GetBarycenter();

      for (auto &vertex : mesh->vertices)
      {
        vertex = (vertex - barycenter) * scale + barycenter;
      }
    }

    void Mesh::operator*=(util::Matrix3D const &scale)
    {
      auto const barycenter = GetBarycenter();

      for (auto &vertex : mesh->vertices)
      {
        scale.timesVector(vertex - barycenter, vertex);
        vertex += barycenter;
      }
    }

    void Mesh::operator+=(LatticePosition const &offset)
    {
      for (auto &vertex : mesh->vertices)
      {
        vertex += offset;
      }
    }

    void Mesh::operator+=(std::vector<LatticePosition> const &displacements)
    {
      assert(displacements.size() == mesh->vertices.size());
      auto i_disp = displacements.begin();

      for (auto &vertex : mesh->vertices)
      {
        vertex += * (i_disp++);
      }
    }

    namespace
    {
      std::shared_ptr<MeshData> initial_tetrahedron()
      {
        std::shared_ptr<MeshData> data(new MeshData);

        // facets at something degrees from one another
        data->vertices.push_back(LatticePosition(0, 0, 0));
        data->vertices.push_back(LatticePosition(1, 0, 1));
        data->vertices.push_back(LatticePosition(1, 1, 0));
        data->vertices.push_back(LatticePosition(0, 1, 1));

        redblood::MeshData::Facet indices;
        indices[0] = 0;
        indices[1] = 1;
        indices[2] = 2;
        data->facets.push_back(indices);
        indices[0] = 0;
        indices[1] = 2;
        indices[2] = 3;
        data->facets.push_back(indices);
        indices[0] = 0;
        indices[1] = 3;
        indices[2] = 1;
        data->facets.push_back(indices);
        indices[0] = 1;
        indices[1] = 3;
        indices[2] = 2;
        data->facets.push_back(indices);
        return data;
      }

      size_t vertex(std::shared_ptr<MeshData> &data,
                    std::map<std::pair<size_t, size_t>, size_t> &vertices, size_t const &i0,
                    size_t const &i1)
      {
        std::pair<size_t, size_t> const indices(i0, i1);
        std::map<std::pair<size_t, size_t>, size_t>::const_iterator i_found =
            vertices.find(indices);

        if (i_found == vertices.end())
        {
          i_found = vertices.find(std::pair<size_t, size_t>(i1, i0));
        }

        if (i_found == vertices.end())
        {
          data->vertices.push_back( (data->vertices[i0] + data->vertices[i1]) * 0.5);
          vertices[indices] = data->vertices.size() - 1;
          return data->vertices.size() - 1;
        }

        return i_found->second;
      }

      void refine(std::shared_ptr<MeshData> &data)
      {
        MeshData::Facets const facets(data->facets);
        data->facets.clear();
        data->facets.resize(facets.size() * 4);
        data->vertices.reserve(data->facets.size() * 3 + data->vertices.size());

        // Container with midpoint indices, so midpoints are only added once
        typedef std::pair<size_t, size_t> Pair;
        std::map<Pair, size_t> new_vertices;

        MeshData::Facets::const_iterator i_orig_facet(facets.begin());
        MeshData::Facets::const_iterator const i_orig_facet_end(facets.end());
        MeshData::Facets::iterator i_facet = data->facets.begin();

        for (; i_orig_facet != i_orig_facet_end; ++i_orig_facet)
        {
          MeshData::Facet::value_type const i0 = (*i_orig_facet)[0];
          MeshData::Facet::value_type const i1 = (*i_orig_facet)[1];
          MeshData::Facet::value_type const i2 = (*i_orig_facet)[2];

          // Adds new vertices halfway through edges
          MeshData::Facet::value_type mid0(vertex(data, new_vertices, i0, i1)),
              mid1(vertex(data, new_vertices, i1, i2)), mid2(vertex(data, new_vertices, i2, i0));

          // Adds all four new faces
          (*i_facet)[0] = i0;
          (*i_facet)[1] = mid0;
          (*i_facet)[2] = mid2;
          ++i_facet;

          (*i_facet)[0] = mid0;
          (*i_facet)[1] = i1;
          (*i_facet)[2] = mid1;
          ++i_facet;

          (*i_facet)[0] = mid1;
          (*i_facet)[1] = i2;
          (*i_facet)[2] = mid2;
          ++i_facet;

          (*i_facet)[0] = mid0;
          (*i_facet)[1] = mid1;
          (*i_facet)[2] = mid2;
          ++i_facet;
        }
      }
    }

    Mesh refine(Mesh data, unsigned int depth)
    {
      if (depth == 0)
      {
        return data.clone();
      }

      std::shared_ptr<MeshData> newData(new MeshData(*data.GetData()));

      for (unsigned int i(0); i < depth; ++i)
      {
        refine(newData);
      }

      return Mesh(newData);
    }

    Mesh tetrahedron(unsigned int depth)
    {
      std::shared_ptr<MeshData> result(initial_tetrahedron());

      for (unsigned int i(0); i < depth; ++i)
      {
        refine(result);
      }

      return Mesh(result);
    }

    Mesh pancakeSamosa(unsigned int depth)
    {
      std::shared_ptr<redblood::MeshData> mesh(new redblood::MeshData);

      // facets at something degrees from one another
      mesh->vertices.push_back(LatticePosition(0, 0, 0));
      mesh->vertices.push_back(LatticePosition(1, 0, 1));
      mesh->vertices.push_back(LatticePosition(1, 1, 0));

      redblood::MeshData::Facet indices;
      indices[0] = 0;
      indices[1] = 1;
      indices[2] = 2;
      mesh->facets.push_back(indices);
      indices[0] = 2;
      indices[1] = 1;
      indices[2] = 0;
      mesh->facets.push_back(indices);

      // Create topology by hand cos we generally don't allow for this kind of
      // ambiguous self-referencing shape.
      std::shared_ptr<redblood::MeshTopology> topo(new redblood::MeshTopology);
      MeshTopology::VertexToFacets::value_type v2f;
      v2f.insert(0);
      v2f.insert(1);
      topo->vertexToFacets.resize(3, v2f);

      MeshTopology::FacetNeighbors::value_type neighbors[2] = { { { 0, 0, 0 } }, { { 1, 1, 1 } } };
      topo->facetNeighbors.push_back(neighbors[1]);
      topo->facetNeighbors.push_back(neighbors[0]);

      return refine(Mesh(mesh, topo), depth);
    }

    Mesh icoSphere(unsigned int depth)
    {
      std::shared_ptr<redblood::MeshData> mesh(new redblood::MeshData);
      // First creates the simplest icosphere
      // Dumbly copied from
      // http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
      PhysicalDistance const t = 0.5 * (1. + std::sqrt(5.0));
      mesh->vertices.clear();
      mesh->vertices.emplace_back(-1, t, 0);
      mesh->vertices.emplace_back( 1, t, 0);
      mesh->vertices.emplace_back(-1, -t, 0);
      mesh->vertices.emplace_back( 1, -t, 0);
      mesh->vertices.emplace_back( 0, -1, t);
      mesh->vertices.emplace_back( 0, 1, t);
      mesh->vertices.emplace_back( 0, -1, -t);
      mesh->vertices.emplace_back( 0, 1, -t);
      mesh->vertices.emplace_back( t, 0, -1);
      mesh->vertices.emplace_back(t, 0, 1);
      mesh->vertices.emplace_back(-t, 0, -1);
      mesh->vertices.emplace_back(-t, 0, 1);
      auto intel_compiler = [](size_t i, size_t j, size_t k) -> MeshData::Facet
      {
        MeshData::Facet result;
        result[0] = i; result[1] = j; result[2] = k;
        return result;
      };
      mesh->facets.clear();
      mesh->facets.push_back(intel_compiler(0, 11, 5));
      mesh->facets.push_back(intel_compiler(0, 5, 1));
      mesh->facets.push_back(intel_compiler(0, 1, 7));
      mesh->facets.push_back(intel_compiler(0, 7, 10));
      mesh->facets.push_back(intel_compiler(0, 10, 11));
      mesh->facets.push_back(intel_compiler(1, 5, 9));
      mesh->facets.push_back(intel_compiler(5, 11, 4));
      mesh->facets.push_back(intel_compiler(11, 10, 2));
      mesh->facets.push_back(intel_compiler(10, 7, 6));
      mesh->facets.push_back(intel_compiler(7, 1, 8));
      mesh->facets.push_back(intel_compiler(3, 9, 4));
      mesh->facets.push_back(intel_compiler(3, 4, 2));
      mesh->facets.push_back(intel_compiler(3, 2, 6));
      mesh->facets.push_back(intel_compiler(3, 6, 8));
      mesh->facets.push_back(intel_compiler(3, 8, 9));
      mesh->facets.push_back(intel_compiler(4, 9, 5));
      mesh->facets.push_back(intel_compiler(2, 4, 11));
      mesh->facets.push_back(intel_compiler(6, 2, 10));
      mesh->facets.push_back(intel_compiler(8, 6, 7));
      mesh->facets.push_back(intel_compiler(9, 8, 1));
      // Then refines it
      for (unsigned int i(0); i < depth; ++i)
      {
        refine(mesh);
      }
      // Finally, makes sure vertices are on unit-sphere
      for (auto &vertex : mesh->vertices)
      {
        vertex.Normalise();
      }
      return mesh;
    }
  }
} // hemelb::redblood
