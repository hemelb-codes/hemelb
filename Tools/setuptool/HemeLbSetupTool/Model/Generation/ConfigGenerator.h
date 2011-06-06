#include <string>
#include <vector>

#include "vtkSmartPointer.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkOBBTree.h"

typedef vtkSmartPointer<vtkPolyDataAlgorithm> vtkPolyDataAlgorithmPtr;
typedef vtkSmartPointer<vtkOBBTree> vtkOBBTreePtr;

/*    _Args = {'StlFile': None,
             'StlFileUnitId': 1,
             'Iolets': ObservableList(),
             'VoxelSize': 0.,
             'SeedPoint': Vector(),
             'OutputConfigFile': None,
             'OutputXmlFile': None,
             'StressType': 1}
                     self.SurfaceSource = profile.SurfaceSource

        self.ClippedSurface = self.ConstructClipPipeline()
        self.Locator = vtkOBBTree()
        self.Locator.SetTolerance(0.)

*/
#define GETTER(name, type) inline type Get##name (void) {return this->name;}
#define SETTER(name, type) inline void Set##name(type val) {this->name = val;}


class ConfigGenerator {
public:
	GETTER(VoxelSize, double);
	SETTER(VoxelSize, double);

	GETTER(OutputConfigFile, std::string)
	SETTER(OutputConfigFile, std::string)

	GETTER(OutputXmlFile, std::string)
	SETTER(OutputXmlFile, std::string)

	GETTER(StressType, int)
	SETTER(StressType, int)

	inline const std::vector<int>& GetIoletIds() {
		return this->IoletIds;
	}
	inline std::vector<int>& GetIoletIds() {
			return this->IoletIds;
	}
	inline void SetIoletIds(std::vector<int> ids) {
		this->IoletIds = std::vector<int>(ids);
	}

	inline void GetSeedPoint(double out[3]) {
		for (unsigned int i=0; i<3; ++i)
			out[i] = this->SeedPoint[i];
		return;
	}
	inline void SetSeedPoint(double out[3]) {
		for (unsigned int i=0; i<3; ++i)
			this->SeedPoint[i] = out[i];
	}
	inline void SetSeedPoint(double x, double y, double z) {
		this->SeedPoint[0] = x;
		this->SeedPoint[1] = y;
		this->SeedPoint[2] = z;
	}

	GETTER(SurfaceSource, vtkPolyDataAlgorithmPtr)
	SETTER(SurfaceSource, vtkPolyDataAlgorithmPtr)

	GETTER(ClippedSurface, vtkPolyDataAlgorithmPtr)
	SETTER(ClippedSurface, vtkPolyDataAlgorithmPtr)

	GETTER(Locator, vtkOBBTreePtr)
	SETTER(Locator, vtkOBBTreePtr)

protected:
	double VoxelSize;
	std::string OutputConfigFile;
	std::string OutputXmlFile;
	int StressType;

	std::vector<int> IoletIds;
	double[3] SeedPoint;
	vtkPolyDataAlgorithmPtr SurfaceSource;
	vtkPolyDataAlgorithmPtr ClippedSurface;
	vtkOBBTreePtr Locator;
};
