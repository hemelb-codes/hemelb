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

#include "GetSet.h"
#include "Iolet.h"

class ConfigWriter;
class Site;

class ConfigGenerator {
public:
//	ConfigGenerator();
//	~ConfigGenerator();
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

//	inline std::string GetOutputXmlFile(void) {
//		return this->OutputXmlFile;
//	}
//	inline void SetOutputXmlFile(std::string val) {
//		this->OutputXmlFile = val;
//	}

	inline int GetStressType(void) {
		return this->StressType;
	}
	inline void SetStressType(int val) {
		this->StressType = val;
	}

	inline std::vector<Iolet>& GetIolets() {
		return this->Iolets;
	}
	inline void SetIolets(std::vector<Iolet> iv) {
		this->Iolets = std::vector<Iolet>(iv);
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

	inline vtkPolyDataAlgorithm* GetSurfaceSource(void) {
		return this->SurfaceSource;
	}
	inline void SetSurfaceSource(vtkPolyDataAlgorithm* val) {
		this->SurfaceSource = val;
	}

	inline vtkPolyDataAlgorithm* GetClippedSurface(void) {
		return this->ClippedSurface;
	}
	inline void SetClippedSurface(vtkPolyDataAlgorithm* val) {
		this->ClippedSurface = val;
	}

	inline vtkOBBTree* GetLocator(void) {
		return this->Locator;
	}
	inline void SetLocator(vtkOBBTree* val) {
		this->Locator = val;
	}

protected:
	double VoxelSize;
	std::string OutputConfigFile;
//	std::string OutputXmlFile;
	int StressType;

	std::vector<Iolet> Iolets;
	double SeedPoint[3];
	vtkPolyDataAlgorithm* SurfaceSource;
	vtkPolyDataAlgorithm* ClippedSurface;
	vtkOBBTree* Locator;
};

#endif // HEMELBSETUPTOOL_CONFIGGENERATOR_H
