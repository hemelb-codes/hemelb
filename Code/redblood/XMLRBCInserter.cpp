#include "XMLRBCInserter.h"
#include "io/xml/XmlAbstractionLayer.h"

namespace hemelb
{
  namespace redblood
  {

XMLRBCInserter::XMLRBCInserter(const std::string & xml_path,
                               const std::string & mesh_path) {
  this->shape = read_mesh(mesh_path);
  read_position_and_scale_from_xml(xml_path, this->position, this->scale);
}

XMLRBCInserter::XMLRBCInserter(const std::string & xml_path,
                               std::istream & mesh_stream) {
  this->shape = read_mesh(mesh_stream);
  read_position_and_scale_from_xml(xml_path, this->position, this->scale);
}

XMLRBCInserter::XMLRBCInserter(const std::string & xml_path,
                               const MeshData & shape) :
  shape(new MeshData(shape)) {
  read_position_and_scale_from_xml(xml_path, this->position, this->scale);
}

XMLRBCInserter::XMLRBCInserter(const MeshData::Vertices::value_type & position,
                               Dimensionless scale, const MeshData & shape) :
    position(position), scale(scale), shape(new MeshData(shape))
{
}

void XMLRBCInserter::InsertCell(std::function<void(CellContainer::value_type)> insertFn) {
  Cell cell(this->shape->vertices, this->shape, this->scale);
  cell += this->position;
  insertFn(cell);
}

void XMLRBCInserter::SetShape(const MeshData & shape)
{
  this->shape.reset(new MeshData(shape));
}
std::shared_ptr<MeshData> XMLRBCInserter::GetShape()
{
  return this->shape;
}
std::shared_ptr<const MeshData> XMLRBCInserter::GetShape() const
{
  return this->shape;
}

void XMLRBCInserter::SetPosition(const MeshData::Vertices::value_type & position)
{
  this->position = position;
}
MeshData::Vertices::value_type & XMLRBCInserter::GetPosition()
{
  return this->position;
}
const MeshData::Vertices::value_type & XMLRBCInserter::GetPosition() const
{
  return this->position;
}

void XMLRBCInserter::SetScale(Dimensionless scale)
{
  this->scale = scale;
}
Dimensionless XMLRBCInserter::GetScale() const
{
  return this->scale;
}

void XMLRBCInserter::read_position_and_scale_from_xml(const std::string & xml_path,
                                                      MeshData::Vertices::value_type & position,
                                                      Dimensionless & scale) {
  using hemelb::io::xml::Document;
  using hemelb::io::xml::Element;

  Document document(xml_path);
  Element root = document.GetRoot();
  Element rbc = root.GetChildOrThrow("redbloodcells");
  Element ins = rbc.GetChildOrThrow("insertion");
  ins.GetAttributeOrThrow("position", position);
  ins.GetAttributeOrThrow("scale", scale);
}

  }
}
