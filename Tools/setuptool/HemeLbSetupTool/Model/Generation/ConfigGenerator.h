#ifndef HEMELBSETUPTOOL_CONFIGGENERATOR_H
#define HEMELBSETUPTOOL_CONFIGGENERATOR_H

#include <string>
#include <vector>

// VTK bits we need
//#include "vtkSmartPointer.h"
class vtkPolyDataAlgorithm;
class vtkOBBTree;
//typedef vtkSmartPointer<vtkPolyDataAlgorithm> vtkPolyDataAlgorithmPtr;
//typedef vtkSmartPointer<vtkOBBTree> vtkOBBTreePtr;
typedef vtkPolyDataAlgorithm* vtkPolyDataAlgorithmPtr;
typedef vtkOBBTree* vtkOBBTreePtr;

#include "GetSet.h"

class ConfigWriter;
class Site;

class ConfigGenerator {
public:
	ConfigGenerator();
	~ConfigGenerator();
	void Execute();
	void ClassifySite(Site& site);

	inline double GetVoxelSize(void) {
		return this->VoxelSize;
	}
	;
	inline void SetVoxelSize(double val) {
		this->VoxelSize = val;
	}
	;

	inline std::string GetOutputConfigFile(void) {
		return this->OutputConfigFile;
	}
	inline void SetOutputConfigFile(std::string val) {
		this->OutputConfigFile = val;
	}

	inline std::string GetOutputXmlFile(void) {
		return this->OutputXmlFile;
	}
	inline void SetOutputXmlFile(std::string val) {
		this->OutputXmlFile = val;
	}

	inline int GetStressType(void) {
		return this->StressType;
	}
	inline void SetStressType(int val) {
		this->StressType = val;
	}

	inline std::vector<int>& GetIoletIds() {
		return this->IoletIds;
	}
	inline void SetIoletIds(std::vector<int> ids) {
		this->IoletIds = std::vector<int>(ids);
	}

	inline void GetSeedPoint(double out[3]) {
		for (unsigned int i = 0; i < 3; ++i)
			out[i] = this->SeedPoint[i];
		return;
	}
	inline void SetSeedPoint(double out[3]) {
		for (unsigned int i = 0; i < 3; ++i)
			this->SeedPoint[i] = out[i];
	}
	inline void SetSeedPoint(double x, double y, double z) {
		this->SeedPoint[0] = x;
		this->SeedPoint[1] = y;
		this->SeedPoint[2] = z;
	}

	inline vtkPolyDataAlgorithmPtr GetSurfaceSource(void) {
		return this->SurfaceSource;
	}
	inline void SetSurfaceSource(vtkPolyDataAlgorithmPtr val) {
		this->SurfaceSource = val;
	}

	inline vtkPolyDataAlgorithmPtr GetClippedSurface(void) {
		return this->ClippedSurface;
	}
	inline void SetClippedSurface(vtkPolyDataAlgorithmPtr val) {
		this->ClippedSurface = val;
	}

	inline vtkOBBTreePtr GetLocator(void) {
		return this->Locator;
	}
	inline void SetLocator(vtkOBBTreePtr val) {
		this->Locator = val;
	}

protected:
	double VoxelSize;
	std::string OutputConfigFile;
	std::string OutputXmlFile;
	int StressType;

	std::vector<int> IoletIds;
	double SeedPoint[3];
	vtkPolyDataAlgorithmPtr SurfaceSource;
	vtkPolyDataAlgorithmPtr ClippedSurface;
	vtkOBBTreePtr Locator;
};

#endif // HEMELBSETUPTOOL_CONFIGGENERATOR_H
