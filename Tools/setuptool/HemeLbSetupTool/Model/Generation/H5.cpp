// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "H5.h"

namespace H5
{
  Object::Object() :
    m_Id(H5I_INVALID_HID)
  {

  }
  Object::Object(hid_t id) :
    m_Id(id)
  {

  }
  Object::~Object()
  {
  }

  PList::PList() :
    Object(H5P_DEFAULT)
  {
  }
  PList::PList(hid_t cls) :
    Object()
  {
    H5_CONSTRUCT(m_Id, H5Pcreate, (cls));
  }
  PList::~PList()
  {
    Close();
  }
  void PList::Close()
  {
    H5_CALL(H5Pclose, (m_Id));
    m_Id = H5I_INVALID_HID;
  }

  /// Default options
  PListSharedPtr PList::Default()
  {
    return PListSharedPtr(new PList());
  }

  /// Properties for object creation
  PListSharedPtr PList::ObjectCreate()
  {
    return PListSharedPtr(new PList(H5P_OBJECT_CREATE));
  }

  /// Properties for file creation
  PListSharedPtr PList::FileCreate()
  {
    return PListSharedPtr(new PList(H5P_FILE_CREATE));
  }

  /// Properties for file access
  PListSharedPtr PList::FileAccess()
  {
    return PListSharedPtr(new PList(H5P_FILE_ACCESS));
  }

  /// Properties for dataset creation
  PListSharedPtr PList::DatasetCreate()
  {
    return PListSharedPtr(new PList(H5P_DATASET_CREATE));
  }

  /// Properties for dataset access
  PListSharedPtr PList::DatasetAccess()
  {
    return PListSharedPtr(new PList(H5P_DATASET_ACCESS));
  }

  /// Properties for raw data transfer
  PListSharedPtr PList::DatasetXfer()
  {
    return PListSharedPtr(new PList(H5P_DATASET_XFER));
  }

  /// Properties for file mounting
  PListSharedPtr PList::FileMount()
  {
    return PListSharedPtr(new PList(H5P_FILE_MOUNT));
  }

  /// Properties for group creation
  PListSharedPtr PList::GroupCreate()
  {
    return PListSharedPtr(new PList(H5P_GROUP_CREATE));
  }

  /// Properties for group access
  PListSharedPtr PList::GroupAccess()
  {
    return PListSharedPtr(new PList(H5P_GROUP_ACCESS));
  }

  /// Properties for datatype creation
  PListSharedPtr PList::DatatypeCreate()
  {
    return PListSharedPtr(new PList(H5P_DATATYPE_CREATE));
  }

  /// Properties for datatype access
  PListSharedPtr PList::DatatypeAccess()
  {
    return PListSharedPtr(new PList(H5P_DATATYPE_ACCESS));
  }

  /// Properties for character encoding when encoding strings or object names
  PListSharedPtr PList::StringCreate()
  {
    return PListSharedPtr(new PList(H5P_STRING_CREATE));
  }

  /// Properties for attribute creation
  PListSharedPtr PList::AttributeCreate()
  {
    return PListSharedPtr(new PList(H5P_ATTRIBUTE_CREATE));
  }

  /// Properties governing the object copying process
  PListSharedPtr PList::ObjectCopy()
  {
    return PListSharedPtr(new PList(H5P_OBJECT_COPY));
  }

  /// Properties governing link creation
  PListSharedPtr PList::LinkCreate()
  {
    return PListSharedPtr(new PList(H5P_LINK_CREATE));
  }

  /// Properties governing link traversal when accessing objects
  PListSharedPtr PList::LinkAccess()
  {
    return PListSharedPtr(new PList(H5P_LINK_ACCESS));
  }
  void PList::SetChunk(const std::vector<hsize_t>& dims)
  {
    H5_CALL(H5Pset_chunk, (m_Id, dims.size(), &dims[0]));
  }
  void PList::SetDeflate(const unsigned level)
  {
    H5_CALL(H5Pset_deflate, (m_Id, level));
  }
  
  GroupSharedPtr CanHaveGroupsDataSets::CreateGroup(
						    const std::string& name, PListSharedPtr createPL,
						    PListSharedPtr accessPL)
  {
    hid_t grp;
    H5_CONSTRUCT(grp, H5Gcreate,
		 (m_Id, name.c_str(), H5P_DEFAULT, createPL->GetId(), accessPL->GetId()));
    GroupSharedPtr ans(new Group(grp));
    return ans;
  }

  DataSetSharedPtr CanHaveGroupsDataSets::CreateDataSet(
							const std::string& name, DataTypeSharedPtr type,
							DataSpaceSharedPtr space, PListSharedPtr createPL,
							PListSharedPtr accessPL)
  {
    hid_t ds;
    H5_CONSTRUCT(ds, H5Dcreate,
		 (m_Id, name.c_str(), type->GetId(), space->GetId(), H5P_DEFAULT, createPL->GetId(), accessPL->GetId()));

    DataSetSharedPtr ans(new DataSet(ds));
    return ans;
  }

  // Open an existing group.
  // The accessPL can be omitted to use the defaults
  GroupSharedPtr CanHaveGroupsDataSets::OpenGroup(
						  const std::string& name, PListSharedPtr accessPL) const
  {
    hid_t grp;
    H5_CONSTRUCT(grp, H5Gopen2,
		 (m_Id, name.c_str(), accessPL->GetId()));
    GroupSharedPtr ans(new Group(grp));
    return ans;
  }

  // Open an existing dataset
  // The accessPL can be omitted to use the defaults
  DataSetSharedPtr CanHaveGroupsDataSets::OpenDataSet(
						      const std::string& name, PListSharedPtr accessPL) const
  {
    hid_t ds;
    H5_CONSTRUCT(ds, H5Dopen2,
		 (m_Id, name.c_str(), accessPL->GetId()));
    DataSetSharedPtr ans(new DataSet(ds));
    return ans;
  }

  CanHaveGroupsDataSets::LinkIterator CanHaveGroupsDataSets::begin()
  {
    // Have to use dynamic because of virtual inheritance
    CanHaveGroupsDataSetsSharedPtr thisSh =
      std::dynamic_pointer_cast < CanHaveGroupsDataSets
				    > (shared_from_this());
    return CanHaveGroupsDataSets::LinkIterator(thisSh);
  }

  CanHaveGroupsDataSets::LinkIterator CanHaveGroupsDataSets::end()
  {
    // Have to use dynamic because of virtual inheritance
    CanHaveGroupsDataSetsSharedPtr thisSh =
      std::dynamic_pointer_cast < CanHaveGroupsDataSets
				    > (shared_from_this());
    return CanHaveGroupsDataSets::LinkIterator(thisSh,
					       GetNumElements());
  }

  CanHaveGroupsDataSets::LinkIterator::LinkIterator(
						    CanHaveGroupsDataSetsSharedPtr grp, hsize_t idx) :
    m_grp(grp), m_idx(-1), m_next(idx), m_size(grp->GetNumElements())
  {
    ++*this;
  }

  const std::string& CanHaveGroupsDataSets::LinkIterator::operator*()
  {
    return m_currentName;
  }
  CanHaveGroupsDataSets::LinkIterator& CanHaveGroupsDataSets::LinkIterator::operator++()
  {
    m_idx = m_next;
    if (m_idx < m_size)
      {
	H5_CALL(H5Literate,
		(m_grp->GetId(), H5_INDEX_NAME, H5_ITER_NATIVE, &m_next, LinkIterator::helper, this));
      }
    return *this;
  }
  bool CanHaveGroupsDataSets::LinkIterator::operator==(
						       const CanHaveGroupsDataSets::LinkIterator& other) const
  {
    if (m_grp == other.m_grp)
      {
	if (m_idx == other.m_idx)
	  {
	    return true;
	  }
      }
    return false;
  }
  herr_t CanHaveGroupsDataSets::LinkIterator::helper(hid_t g_id,
						     const char *name, const H5L_info_t *info, void *op_data)
  {
    CanHaveGroupsDataSets::LinkIterator* iter =
      static_cast<CanHaveGroupsDataSets::LinkIterator*>(op_data);
    iter->m_currentName = name;
    return 1;
  }

  AttributeSharedPtr CanHaveAttributes::CreateAttribute(
							const std::string& name, DataTypeSharedPtr type,
							DataSpaceSharedPtr space)
  {
    return Attribute::Create(m_Id, name, type, space);
  }

  AttributeSharedPtr CanHaveAttributes::OpenAttribute(
						      const std::string& name)
  {
    return Attribute::Open(m_Id, name);
  }

  int CanHaveAttributes::GetNumAttr() const
  {
    H5O_info_t info;
    H5_CALL(H5Oget_info, (m_Id, &info));
    return info.num_attrs;
  }

  CanHaveAttributes::AttrIterator CanHaveAttributes::attr_begin()
  {
    // Have to use dynamic because of virtual inheritance
    CanHaveAttributesSharedPtr thisSh = std::dynamic_pointer_cast
      < CanHaveAttributes > (shared_from_this());
    return CanHaveAttributes::AttrIterator(thisSh);

  }
  CanHaveAttributes::AttrIterator CanHaveAttributes::attr_end()
  {
    // Have to use dynamic because of virtual inheritance
    CanHaveAttributesSharedPtr thisSh = std::dynamic_pointer_cast
      < CanHaveAttributes > (shared_from_this());
    return CanHaveAttributes::AttrIterator(thisSh, GetNumAttr());
  }

  CanHaveAttributes::AttrIterator::AttrIterator(
						CanHaveAttributesSharedPtr obj, hsize_t idx) :
    m_obj(obj), m_idx(-1), m_next(idx), m_size(obj->GetNumAttr())
  {
    ++*this;
  }

  const std::string& CanHaveAttributes::AttrIterator::operator*()
  {
    return m_currentName;
  }
  CanHaveAttributes::AttrIterator& CanHaveAttributes::AttrIterator::operator++()
  {
    m_idx = m_next;
    if (m_next < m_size)
      {
	H5_CALL(H5Aiterate2,
		(m_obj->GetId(), H5_INDEX_CRT_ORDER, H5_ITER_INC, &m_next, AttrIterator::helper, this));
      }
    return *this;
  }
  bool CanHaveAttributes::AttrIterator::operator==(
						   const CanHaveAttributes::AttrIterator& other) const
  {
    if (m_obj == other.m_obj)
      {
	if (m_idx == other.m_idx)
	  {
	    return true;
	  }
      }
    return false;
  }

  herr_t CanHaveAttributes::AttrIterator::helper(hid_t g_id,
						 const char *name, const H5A_info_t *info, void *op_data)
  {
    CanHaveAttributes::AttrIterator* iter =
      static_cast<CanHaveAttributes::AttrIterator*>(op_data);
    iter->m_currentName = name;
    return 1;
  }

  DataSpaceSharedPtr DataSpace::Null()
  {
    DataSpaceSharedPtr ans(new DataSpace);
    H5_CONSTRUCT(ans->m_Id, H5Screate, (H5S_NULL));
    return ans;
  }

  DataSpaceSharedPtr DataSpace::Scalar()
  {
    DataSpaceSharedPtr ans(new DataSpace);
    H5_CONSTRUCT(ans->m_Id, H5Screate, (H5S_SCALAR));
    return ans;
  }
  DataSpaceSharedPtr DataSpace::OneD(hsize_t size)
  {
    DataSpaceSharedPtr ans(new DataSpace);
    H5_CONSTRUCT(ans->m_Id, H5Screate_simple, (1, &size, NULL));
    return ans;
  }

  DataSpace::DataSpace() :
    Object()
  {
  }

  DataSpace::DataSpace(hid_t id) :
    Object(id)
  {
  }

  DataSpace::DataSpace(const std::vector<hsize_t>& dims) :
    Object()
  {
    int rank = dims.size();
    H5_CONSTRUCT(m_Id, H5Screate_simple, (rank, &dims[0], NULL));
  }

  DataSpace::DataSpace(const hsize_t size, const hsize_t max) :
    Object()
  {
    const hsize_t* max_p = NULL;
    if (max != (H5S_UNLIMITED - 1))
      max_p = &max;
    H5_CONSTRUCT(m_Id, H5Screate_simple, (1, &size, max_p));
  }

  DataSpace::DataSpace(const std::vector<hsize_t>& dims,
		       const std::vector<hsize_t>& max_dims) :
    Object()
  {
    int rank = dims.size();
    H5_CONSTRUCT(m_Id, H5Screate_simple,
		 (rank, &dims[0], &max_dims[0]));
  }

  DataSpace::~DataSpace()
  {
    Close();
  }

  void DataSpace::Close()
  {
    H5_CALL(H5Sclose, (m_Id));
    m_Id = H5I_INVALID_HID;
  }
  void DataSpace::SelectRange(const hsize_t start, const hsize_t count)
  {
    H5_CALL(H5Sselect_hyperslab,
	    (m_Id, H5S_SELECT_SET, &start, NULL, &count, NULL));
  }

  DataType::DataType(hid_t id) :
    Object(id)
  {
  }
  DataTypeSharedPtr DataType::String(size_t len)
  {
    DataTypeSharedPtr s1 = PredefinedDataType::CS1();
    DataTypeSharedPtr ans = s1->Copy();
    if (len == 0)
      len = H5T_VARIABLE;
    H5_CALL(H5Tset_size, (ans->GetId(), len));
    return ans;
  }

  void DataType::Close()
  {
    H5_CALL(H5Tclose, (m_Id));
    m_Id = H5I_INVALID_HID;
  }

  DataTypeSharedPtr DataType::Copy() const
  {
    hid_t ans_id = H5I_INVALID_HID;
    H5_CONSTRUCT(ans_id, H5Tcopy, (m_Id));
    DataTypeSharedPtr ans(new DataType(ans_id));
    return ans;
  }

  DataTypeSharedPtr PredefinedDataType::CS1()
  {
    return DataTypeSharedPtr(new PredefinedDataType(H5T_C_S1));
  }

  PredefinedDataType::PredefinedDataType(hid_t id) :
    DataType(id)
  {
  }

  void PredefinedDataType::Close()
  {
    // No-op
    m_Id = H5I_INVALID_HID;
  }

  template<>
  const hid_t DataTypeTraits<char>::NativeType = H5T_NATIVE_CHAR;

  template<>
  const hid_t DataTypeTraits<int>::NativeType = H5T_NATIVE_INT;

  template<>
  const hid_t DataTypeTraits<unsigned int>::NativeType =
    H5T_NATIVE_UINT;

  template<>
  const hid_t DataTypeTraits<unsigned long long>::NativeType =
    H5T_NATIVE_ULLONG;

  template<>
  const hid_t DataTypeTraits<float>::NativeType = H5T_NATIVE_FLOAT;
  template<>
  const hid_t DataTypeTraits<double>::NativeType = H5T_NATIVE_DOUBLE;

  AttributeSharedPtr Attribute::Create(hid_t parent,
				       const std::string& name, DataTypeSharedPtr type,
				       DataSpaceSharedPtr space)
  {
    hid_t id;
    H5_CONSTRUCT(id, H5Acreate,
		 (parent, name.c_str(), type->GetId(), space->GetId(), H5P_DEFAULT, H5P_DEFAULT));
    return AttributeSharedPtr(new Attribute(id));
  }

  AttributeSharedPtr Attribute::Open(hid_t parent,
				     const std::string& name)
  {
    hid_t id;
    H5_CONSTRUCT(id, H5Aopen, (parent, name.c_str(), H5P_DEFAULT));
    return AttributeSharedPtr(new Attribute(id));
  }

  Attribute::~Attribute()
  {
    Close();
  }
  void Attribute::Close()
  {
    H5_CALL(H5Aclose, (m_Id));
    m_Id = H5I_INVALID_HID;
  }

  DataSpaceSharedPtr Attribute::GetSpace() const
  {
    return DataSpaceSharedPtr(new DataSpace(H5Aget_space(m_Id)));
  }

  File::File(hid_t id) :
    Object(id)
  {
  }
  FileSharedPtr File::Create(const std::string& filename,
			     unsigned mode, PListSharedPtr createPL,
			     PListSharedPtr accessPL)
  {
    hid_t id;
    H5_CONSTRUCT(id, H5Fcreate,
		 (filename.c_str(), mode, createPL->GetId(), accessPL->GetId()));
    return FileSharedPtr(new File(id));
  }
  FileSharedPtr File::Open(const std::string& filename, unsigned mode,
			   PListSharedPtr accessPL)
  {
    hid_t id;
    H5_CONSTRUCT(id, H5Fopen,
		 (filename.c_str(), mode, accessPL->GetId()));
    return FileSharedPtr(new File(id));
  }

  File::~File()
  {
    Close();
  }

  void File::Close()
  {
    H5_CALL(H5Fclose, (m_Id));
    m_Id = H5I_INVALID_HID;
  }
  hsize_t File::GetNumElements()
  {
    GroupSharedPtr root = OpenGroup("/");
    return root->GetNumElements();
  }
  Group::Group(hid_t id) :
    Object(id)
  {
  }

  Group::~Group()
  {
    Close();
  }

  void Group::Close()
  {
    H5_CALL(H5Gclose, (m_Id));
    m_Id = H5I_INVALID_HID;
  }

  hsize_t Group::GetNumElements()
  {
    H5G_info_t info;
    H5_CALL(H5Gget_info, (m_Id, &info));
    return info.nlinks;
  }

  DataSet::DataSet(hid_t id) :
    Object(id)
  {
  }

  DataSet::~DataSet()
  {
    Close();
  }

  void DataSet::Close()
  {
    H5_CALL(H5Dclose, (m_Id));
    m_Id = H5I_INVALID_HID;
  }

  DataSpaceSharedPtr DataSet::GetSpace() const
  {
    return DataSpaceSharedPtr(new DataSpace(H5Dget_space(m_Id)));
  }

}
