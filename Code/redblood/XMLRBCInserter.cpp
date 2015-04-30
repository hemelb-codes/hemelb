#include "XMLRBCInserter.h"
#include "io/xml/XmlAbstractionLayer.h"

namespace hemelb
{
  namespace redblood
  {

XMLRBCInserter::XMLRBCInserter(std::function<bool()> condition,
                               const std::string & mesh_path,
                               const lb::iolets::InOutLet * inlet,
                               Dimensionless scale) :
    condition(condition), inlet(inlet), scale(scale)
{
  this->shape = read_mesh(mesh_path);
}

XMLRBCInserter::XMLRBCInserter(std::function<bool()> condition,
                               std::istream & mesh_stream,
                               const lb::iolets::InOutLet * inlet,
                               Dimensionless scale) :
    condition(condition), inlet(inlet), scale(scale)
{
  this->shape = read_mesh(mesh_stream);
}

XMLRBCInserter::XMLRBCInserter(std::function<bool()> condition, const MeshData & shape,
                               const lb::iolets::InOutLet * inlet,
                               Dimensionless scale) :
    condition(condition), shape(new MeshData(shape)), inlet(inlet), scale(scale)
{
}

void XMLRBCInserter::operator()(std::function<void(CellContainer::value_type)> insertFn) {
  while (condition()) {
    std::shared_ptr<Cell> cell = std::make_shared<Cell>(this->shape->vertices, Mesh(*this->shape), this->scale);
    const std::shared_ptr<FlowExtension> flowExt = this->inlet->GetFlowExtension();
    *cell += this->inlet->GetPosition();
    if (flowExt)
      *cell += flowExt->normal * flowExt->length;
    insertFn(cell);
  }
}

void XMLRBCInserter::SetShape(const MeshData & shape)
{
  this->shape.reset(new MeshData(shape));
}
std::shared_ptr<const MeshData> XMLRBCInserter::GetShape() const
{
  return this->shape;
}

void XMLRBCInserter::SetInLet(const lb::iolets::InOutLet * inlet) {
  this->inlet = inlet;
}
const lb::iolets::InOutLet * XMLRBCInserter::GetInLet() const {
  return this->inlet;
}

void XMLRBCInserter::SetScale(Dimensionless scale)
{
  this->scale = scale;
}
Dimensionless XMLRBCInserter::GetScale() const
{
  return this->scale;
}

void XMLRBCInserter::SetCondition(std::function<bool()> condition)
{
  this->condition = condition;
}
std::function<bool()> XMLRBCInserter::GetCondition() const
{
  return this->condition;
}

  }
}
