#include "io/xml/XmlAbstractionLayer.h"
#include "RBCInserter.h"

namespace hemelb
{
  namespace redblood
  {

RBCInserter::RBCInserter(std::function<bool()> condition,
                         const std::string & mesh_path,
                         const lb::iolets::InOutLet * inlet,
                         Dimensionless scale) :
    condition(condition), inlet(inlet), scale(scale)
{
  this->shape = read_mesh(mesh_path);
}

RBCInserter::RBCInserter(std::function<bool()> condition,
                         std::istream & mesh_stream,
                         const lb::iolets::InOutLet * inlet,
                         Dimensionless scale) :
    condition(condition), inlet(inlet), scale(scale)
{
  this->shape = read_mesh(mesh_stream);
}

RBCInserter::RBCInserter(std::function<bool()> condition,
                         const MeshData & shape,
                         const lb::iolets::InOutLet * inlet,
                         Dimensionless scale) :
    condition(condition), shape(new MeshData(shape)), inlet(inlet), scale(scale)
{
}

void RBCInserter::operator()(std::function<void(CellContainer::value_type)> insertFn) {
  while (condition()) {
    std::shared_ptr<Cell> cell = std::make_shared<Cell>(this->shape->vertices,
                                                        Mesh(*this->shape),
                                                        this->scale);
    const std::shared_ptr<FlowExtension> flowExt = this->inlet->GetFlowExtension();
    *cell += this->inlet->GetPosition();
    if (flowExt)
      *cell += flowExt->normal * flowExt->length;
    insertFn(cell);
  }
}

void RBCInserter::SetShape(const MeshData & shape)
{
  this->shape.reset(new MeshData(shape));
}
std::shared_ptr<const MeshData> RBCInserter::GetShape() const
{
  return this->shape;
}

void RBCInserter::SetInLet(const lb::iolets::InOutLet * inlet) {
  this->inlet = inlet;
}
const lb::iolets::InOutLet * RBCInserter::GetInLet() const {
  return this->inlet;
}

void RBCInserter::SetScale(Dimensionless scale)
{
  this->scale = scale;
}
Dimensionless RBCInserter::GetScale() const
{
  return this->scale;
}

void RBCInserter::SetCondition(std::function<bool()> condition)
{
  this->condition = condition;
}
std::function<bool()> RBCInserter::GetCondition() const
{
  return this->condition;
}

  }
}
