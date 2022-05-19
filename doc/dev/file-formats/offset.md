# Offset files

These store where each MPI rank of the simulation wrote its contiguous
chunk of data. This makes resuming a checkpoint on the same number of
ranks much easier and faster.

The file has a header and a body.

## Header
This contains, encoded as uint32, in order:
 - HemeLbMagicNumber
 - OffsetMagicNumber (0x6F666604 == 'xtr\eof')
 - OffsetVersionNumber (currently 1)
 - Number of ranks

# Body

This holds number of ranks + 1 entries, XDR encoded as uint64. Each
gives the offset (in number of bytes) into the extraction file of the
start of that rank's chunk. The final value holds the
just-past-the-end value.

The length of a single timestep of pure data in an extraction file is
given by `data[n_ranks] - data[0]`.
