#ifndef HEMELBSETUPTOOL_CYLINDERGENERATOR_H
#define HEMELBSETUPTOOL_CYLINDERGENERATOR_H

#include <string>
#include <vector>

// VTK bits we need
class vtkPolyData;
class vtkPoints;
class vtkIdList;
class vtkIntArray;

#include "Index.h"
#include "GetSet.h"
#include "Iolet.h"
#include "GenerationError.h"

class GeometryWriter;
class Site;
class BlockWriter;

class CylinderGenerator {
public:
	CylinderGenerator();
	~CylinderGenerator();
	void Execute() throw (GenerationError);

	inline double GetVoxelSizeMetres(void) {
		return this->VoxelSizeMetres;
	}
	inline void SetVoxelSizeMetres(double val) {
		this->VoxelSizeMetres = val;
	}

	inline double GetVoxelSizeWorking(void) {
		return 1;
	}

	inline std::string GetOutputGeometryFile(void) {
		return this->OutputGeometryFile;
	}
	inline void SetOutputGeometryFile(std::string val) {
		this->OutputGeometryFile = val;
	}

	inline std::vector<Iolet*>& GetIolets() {
		return this->Iolets;
	}
	inline void SetIolets(std::vector<Iolet*> iv) {
		this->Iolets = std::vector<Iolet*>(iv);
	}

	inline void SetCylinderCentre(Vector v) {
		this->Cylinder->Centre = v;
	}
	inline void SetCylinderAxis(Vector n) {
		this->Cylinder->Axis = n;
	}
	inline void SetCylinderRadius(double r) {
		this->Cylinder->Radius = r;
	}
	inline void SetCylinderLength(double l) {
		this->Cylinder->Length = l;
	}

private:
	void ClassifySite(Site& site);
	void WriteSolidSite(BlockWriter& blockWriter, Site& site);
	void WriteFluidSite(BlockWriter& blockWriter, Site& site);
	int ComputeIntersections(Site& from, Site& to);
	// Members set from outside to initialise
	double VoxelSizeMetres;
	std::string OutputGeometryFile;
	std::vector<Iolet*> Iolets;

	struct CylinderData {
		// Cylinder parameters
		Vector Centre;
		Vector Axis;
		double Radius;
		double Length;
	};
	CylinderData* Cylinder;
	// Members used internally
	vtkPoints* hitPoints;
	vtkIdList* hitCellIds;
	vtkIntArray* IoletIdArray;
};

#endif // HEMELBSETUPTOOL_CYLINDERGENERATOR_H
