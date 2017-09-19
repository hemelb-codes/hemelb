#ifndef HEMELBSETUPTOOL_TESTRESOURCES_MESHES_HPP
#define HEMELBSETUPTOOL_TESTRESOURCES_MESHES_HPP

#include <vector>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#include <cppunit/extensions/HelperMacros.h>

#include "../Vector.h"
#include "../Iolet.h"

#include "helpers.h"

struct MeshData {
  std::vector<Vector> points;
  std::vector<Index> triangles;
  std::vector<Vector> normals;
  std::vector<int> labels;
  std::vector<Iolet> iolets;
};

class SimpleMeshFactory {
public:
  static std::shared_ptr<MeshData> MkTrivial() {
    auto ans = std::make_shared<MeshData>();
    ans->points = {{1.2, 1.2, 1.2},
		   {1.2, 1.2, 2.2},
		   {1.2, 2.2, 1.2},
		   {1.2, 2.2, 2.2}};
    ans->triangles = {{0,1,2},
		      {2,1,3}};
    ans->normals = {{-1, 0, 0},
		    {-1, 0, 0}};
    ans->labels = {-1,
		   -1};
    return ans;
  }

  static std::shared_ptr<MeshData> MkSphere() {
    if (!sphere_mesh) {
      auto vtk = ReadSphere();
      sphere_mesh = VtkToMesh(vtk);
    }
    return sphere_mesh;
  }
  
  static std::shared_ptr<MeshData> MkDuct() {
    auto ans = std::make_shared<MeshData>();
    ans->points = {
      {0.5, 0.5, 2.5},
      {0.5, 3.5, 2.5},
      {3.5, 3.5, 2.5},
      {3.5, 0.5, 2.5},
      
      {0.5, 0.5, 14.5},
      {0.5, 3.5, 14.5},
      {3.5, 3.5, 14.5},
      {3.5, 0.5, 14.5}
    };
    
    ans->triangles = {
      // -x
      {0,4,1},
      {4,5,1},
      // +y
      {1,5,2},
      {5,6,2},
      // +x
      {2,6,3},
      {6,7,3},
      // -y
      {3,7,0},
      {7,4,0},
      // inlet
      {0,1,3},
      {1,2,3},
      // outlet
      {5,4,7},
      {6,5,7}
    };
    
    ans->normals = {
      {-1,0,0},
      {-1,0,0},
      {0,+1,0},
      {0,+1,0},
      {+1,0,0},
      {+1,0,0},
      {0,-1,0},
      {0,-1,0},
      {0,0,-1},
      {0,0,-1},
      {0,0,+1},
      {0,0,+1}
    };
    
    ans->labels = {-1,-1, -1,-1, -1,-1, -1,-1, 0,0, 1,1};
    ans->iolets = {Iolet{0, true}, Iolet{0, false}};
    return ans;
  };

private:
  typedef vtkSmartPointer<vtkXMLPolyDataReader> PDR_ptr;
  typedef vtkSmartPointer<vtkPolyData> PD_ptr;

  
  static PD_ptr ReadSphere() {
    auto reader = PDR_ptr::New();
    const std::string&& sphere = GetResource("sphere.vtp");
    reader->SetFileName(sphere.c_str());
    reader->Update();
    CPPUNIT_ASSERT(reader->GetOutput() != nullptr);
    return reader->GetOutput();
  }

  static std::shared_ptr<MeshData> sphere_mesh;
  
  static std::shared_ptr<MeshData> VtkToMesh(PD_ptr pd) {
    auto ans = std::make_shared<MeshData>();
    
    // Points
    auto nPoints = pd->GetNumberOfPoints();
    ans->points.resize(nPoints);
    for (auto i = 0; i < nPoints; ++i) {
      auto pt = pd->GetPoint(i);
      for (auto d = 0; d<3; ++d)
	ans->points[i][d] = pt[d];
    }

    
    auto nTri = pd->GetNumberOfCells();
    ans->triangles.resize(nTri);
    auto raw_tris = pd->GetPolys()->GetData();
    
    ans->normals.resize(nTri);
    auto normals = pd->GetCellData()->GetNormals();
    
    ans->labels.resize(nTri);
    auto labels = pd->GetCellData()->GetScalars();

    auto cell_idx = 0;
    for (auto i = 0; i < nTri; ++i) {
      CPPUNIT_ASSERT(raw_tris->GetTuple1(cell_idx) == 3);
      cell_idx++;
      
      auto n = normals->GetTuple(i);
      
      for (auto d = 0; d < 3; ++d) {
	ans->triangles[i][d] = raw_tris->GetTuple1(cell_idx + d);
	ans->normals[i][d] = n[d];
      }
      ans->labels[i] = labels->GetTuple1(i);
      
      cell_idx += 3;
    }

    return ans;
  }

};
std::shared_ptr<MeshData> SimpleMeshFactory::sphere_mesh;
#endif

