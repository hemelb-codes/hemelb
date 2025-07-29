// LocalPropertyHdf5Output.cc
#include "LocalPropertyHdf5Output.h"
#ifdef USE_HDF5
#include <vector>
#include <cassert>

namespace hemelb::io {

#pragma pack(push,1)
struct PointRec
{
  uint32_t x, y, z;
  double   pressure;
};
#pragma pack(pop)

LocalPropertyHdf5Output::LocalPropertyHdf5Output(IterableDataSource& src,
                                                 const PropertyOutputFile& spec,
                                                 const net::IOCommunicator& comms)
  : LocalPropertyOutputBase(src, spec, comms),
    comm_(comms.GetCommunicator())
{
  // ---- 定义 compound datatype ----
  memType_ = H5Tcreate(H5T_COMPOUND, sizeof(PointRec));
  H5Tinsert(memType_, "x", HOFFSET(PointRec, x), H5T_NATIVE_UINT32);
  H5Tinsert(memType_, "y", HOFFSET(PointRec, y), H5T_NATIVE_UINT32);
  H5Tinsert(memType_, "z", HOFFSET(PointRec, z), H5T_NATIVE_UINT32);
  H5Tinsert(memType_, "p", HOFFSET(PointRec, pressure), H5T_NATIVE_DOUBLE);
}

void LocalPropertyHdf5Output::StartFile(const std::string& fn)
{
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, comm_, MPI_INFO_NULL);
  file_ = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);

  const hsize_t global = globalSiteCount;      // 这个成员继承自基类
  fileSpace_ = H5Screate_simple(1, &global, nullptr);
  dset_      = H5Dcreate(file_, "points", memType_, fileSpace_,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void LocalPropertyHdf5Output::Write(unsigned long /*ts*/, unsigned long /*total*/)
{
  std::vector<PointRec> buf;
  buf.reserve(localSiteCount);

  while (dataSource.ReadNext()) {
    const auto& pos = dataSource.GetPosition();
    buf.push_back({pos.x, pos.y, pos.z, dataSource.GetPressure()});
  }
  assert(buf.size() == localSiteCount);

  // 选择 hyperslab
  const hsize_t start = localSiteStartRow;
  const hsize_t count = localSiteCount;
  H5Sselect_hyperslab(fileSpace_, H5S_SELECT_SET, &start, nullptr, &count, nullptr);

  hid_t memSpace = H5Screate_simple(1, &count, nullptr);
  hid_t xfer     = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);

  H5Dwrite(dset_, memType_, memSpace, fileSpace_, xfer, buf.data());

  H5Sclose(memSpace);
  H5Pclose(xfer);
}

void LocalPropertyHdf5Output::Close()
{
  if (dset_      >= 0) H5Dclose(dset_);
  if (fileSpace_ >= 0) H5Sclose(fileSpace_);
  if (memType_   >= 0) H5Tclose(memType_);
  if (file_      >= 0) H5Fclose(file_);
}

} // namespace hemelb::io
#endif // USE_HDF5
