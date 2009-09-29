#ifndef MESH
#ifndef DCM_INPUT
#define DCM_INPUT

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGDCMImageIO.h"

#include "itkOrientedImage.h" 
#include "itkGDCMSeriesFileNames.h" 
#include "itkImageSeriesReader.h" 

#include <list>
#include <fstream>

class DICOM_read { 

	public:
		DICOM_read(char* dir_path);
		void DICOM_volume_size(int &x, int &y, int &z);
		void DICOM_copy_data(signed short *a, int length);
		// Method to get the data into a format for SegTool

	private:

		typedef signed short PixelType;

		// Statically define the dimension of the data as 3 here
		typedef itk::Image< PixelType, 3 > ImageType;

		typedef itk::ImageSeriesReader< ImageType > ReaderType;

		typedef itk::GDCMImageIO ImageIOType;

		typedef itk::GDCMSeriesFileNames NamesGeneratorType; 
		NamesGeneratorType::Pointer nameGenerator;

		typedef std::vector< std::string >    SeriesIdContainer; 

		typedef std::vector< std::string >   FileNamesContainer; 

		FileNamesContainer fileNames; 

	  	std::string seriesIdentifier;
		ImageType::Pointer image;

		int x, y, z;

		typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;

};

#endif // DCM_INPUT
#endif // MESH
