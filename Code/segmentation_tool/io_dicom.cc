#ifndef MESH
#include "io_dicom.h"

DICOM_read :: DICOM_read(char* dir_path) {

	ReaderType::Pointer reader = ReaderType::New();

	ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
	
	
	reader->SetImageIO( gdcmImageIO );
  
	nameGenerator = NamesGeneratorType::New(); 

	nameGenerator->SetDirectory( dir_path ); 

	try {
		std::cout << std::endl << "The directory: " << std::endl; 
		std::cout << std::endl << dir_path << std::endl << std::endl; 
		std::cout << "Contains the following DICOM Series: "; 
		std::cout << std::endl << std::endl; 

		const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

		SeriesIdContainer::const_iterator seriesItr = seriesUID.begin(); 
		SeriesIdContainer::const_iterator seriesEnd = seriesUID.end(); 

		while ( seriesItr != seriesEnd ) {
			std::cout << seriesItr->c_str() << std::endl; 
			seriesItr++; 
		}
		seriesIdentifier = seriesUID.begin()->c_str();
		
		std::cout << std::endl << std::endl; 
		std::cout << "Now reading series: " << std::endl << std::endl; 
		std::cout << seriesIdentifier << std::endl; 
		
		fileNames = nameGenerator->GetFileNames( seriesIdentifier ); 
		
		std::cout << fileNames.size() << " files in series, names: " << std::endl;
		
		for (int i = 0; i < fileNames.size(); i++)
			std::cout << fileNames[i] << std::endl;
		std::cout << std::endl;
		
		reader->SetFileNames( fileNames ); 
		
		try {
		        reader->Update(); 
		} catch (itk::ExceptionObject &ex) {
			std::cout << ex << std::endl; 
			// return EXIT_FAILURE; 
		}
		image = reader->GetOutput();
		
	} catch (itk::ExceptionObject &ex) {
		std::cout << ex << std::endl; 
		// return EXIT_FAILURE; 
	}
	// return EXIT_SUCCESS; 
}


void DICOM_read :: DICOM_volume_size(int &x, int &y, int &z) {

	x = image->GetLargestPossibleRegion().GetSize()[0];
	y = image->GetLargestPossibleRegion().GetSize()[1];
	z = image->GetLargestPossibleRegion().GetSize()[2];
}


void DICOM_read :: DICOM_copy_data(signed short *a, int length) {

	ConstIteratorType constIterator( image, image->GetRequestedRegion() );
	
	constIterator.GoToBegin();
	
	for (int i = 0; i < length; i++) {
		a[i] = constIterator.Get();
		++constIterator;
	}
}
#endif // MESH
