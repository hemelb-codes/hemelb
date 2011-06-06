#include "ConfigWriter.h"

ConfigWriter::ConfigWriter(std::string OutputConfigFile, int StressType,
		int BlockSize, int BlockCounts[3], double VoxelSize, double Origin[3]) {
	this->OutputConfigFile = OutputConfigFile;
	this->StressType = StressType;
	this->BlockSize = BlockSize;
	this->VoxelSize = VoxelSize;

	for (unsigned int i = 0; i < 3; ++i) {
		this->BlockCounts[i] = BlockCounts[i];
		this->Origin[i] = Origin[i];
	}
	// Open for read/write, creating if non-existent, truncating otherwise.
	this->file = std::open(this->OutputConfigFile, O_RDWR | O_CREAT | O_TRUNC);
}
/*
 # Truncate, and open for read & write in binary mode
 self.file = file(OutputConfigFile, 'w+b')

 encoder = xdrlib.Packer()
 # Write the preamble, starting with the stress type
 encoder.pack_int(StressType)

 # Blocks in each dimension
 for count in BlockCounts:
 encoder.pack_uint(int(count))
 continue
 # Sites along 1 dimension of a block
 encoder.pack_uint(BlockSize)

 # Voxel Size, in metres
 encoder.pack_double(VoxelSize)
 # Position of site index (0,0,0) in block index (0,0,0), in
 # metres in the STL file's coordinate system
 for ori in Origin:
 encoder.pack_double(ori)
 continue

 # Write this to the file
 self.file.write(encoder.get_buffer())

 # Reset, we're going to write a dummy header now
 encoder.reset()
 # For each block
 for ignored in xrange(np.prod(BlockCounts)):
 encoder.pack_uint(0) # n fluid sites
 encoder.pack_uint(0) # n bytes
 continue
 # Note the start of the header
 self.headerStart = self.file.tell()
 # Write it
 self.file.write(encoder.get_buffer())
 # Note the start of the body
 self.bodyStart = self.file.tell()

 self.HeaderEncoder = xdrlib.Packer()
 return

 def Close(self):
 """Rewrites the header map and closes the file.
 """
 head = self.HeaderEncoder.get_buffer()
 assert len(head) == (self.bodyStart - self.headerStart)

 self.file.seek(self.headerStart)
 self.file.write(head)
 self.file.close()
 return

 @contextmanager
 def BlockStarted(self):
 """A context manager (see PEP 343) for use in a 'with'
 statement.

 This context manager will return an enhanced xdrlib.Packer
 instance which should be used to write the data for a single
 MacroBlock within the body of the 'with' statement. Every time
 a fluid site is written, the IncrementFluidSitesCount method
 of the Packer should be called. This is used at the end of the
 'with' block to write a record in the header.

 """
 # The encode we will yiled
 encoder = xdrlib.Packer()

 # This function will be set to the IncrementFluidSitesCount
 # attribute of the Packer
 def incrementor():
 """Increase the count of fluid sites that have been
 written by the active BlockStarted context.
 """
 incrementor.count += 1
 return
 incrementor.count = 0

 encoder.IncrementFluidSitesCount = incrementor
 # Give the altered XDR encoder back to our caller
 yield encoder

 # Write our record into the header buffer
 self.HeaderEncoder.pack_uint(incrementor.count)

 if incrementor.count == 0:
 # We mustn't write anything to the file, so note that it
 # takes zero bytes
 self.HeaderEncoder.pack_uint(0)
 else:
 # Write the block to the main file

 blockStart = self.file.tell()
 self.file.write(encoder.get_buffer())
 blockEnd = self.file.tell()
 self.HeaderEncoder.pack_uint(blockEnd - blockStart)
 pass

 return

 pass
 *
 */
