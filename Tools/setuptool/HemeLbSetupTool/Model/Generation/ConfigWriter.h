#include <string>
#include <cstdio>

class ConfigWriter {
public:
	ConfigWriter(std::string OutputConfigFile, int StressType, int BlockSize,
			int BlockCounts[3], double VoxelSize, double Origin[3]);

protected:
	std::string OutputConfigFile;
	int StressType;
	int BlockSize;
	int BlockCounts[3];
	double VoxelSize;
	double Origin[3];

	std::FILE* file;
	int headerStart;
	int bodyStart;
};
