#include "io/xml/XmlAbstractionLayer.h"
#include "RBCInserter.h"

namespace hemelb
{
  namespace redblood
  {

    RBCInserter::RBCInserter(std::function<bool()> condition, const std::string & mesh_path,
                             std::vector<lb::iolets::InOutLet *> inlets, Cell::Moduli moduli,
                             Dimensionless scale) :
        condition(condition), inlets(inlets), moduli(moduli), scale(scale)
    {
      shape = read_mesh(mesh_path);
    }

    RBCInserter::RBCInserter(std::function<bool()> condition, std::istream & mesh_stream,
                             std::vector<lb::iolets::InOutLet *> inlets, Cell::Moduli moduli,
                             Dimensionless scale) :
        condition(condition), inlets(inlets), moduli(moduli), scale(scale)
    {
      shape = read_mesh(mesh_stream);
    }

    RBCInserter::RBCInserter(std::function<bool()> condition, const MeshData & shape,
                             std::vector<lb::iolets::InOutLet *> inlets, Cell::Moduli moduli,
                             Dimensionless scale) :
        condition(condition), shape(new MeshData(shape)), inlets(inlets), moduli(moduli),
            scale(scale)
    {
    }

    void RBCInserter::operator()(std::function<void(CellContainer::value_type)> insertFn) const
    {
      for (auto inlet : inlets)
      {
        if (condition())
        {
          std::shared_ptr < Cell > cell = std::make_shared < Cell
              > (shape->vertices, Mesh(*shape), scale);
          cell->moduli = moduli;
          const std::shared_ptr<FlowExtension> flowExt = inlet->GetFlowExtension();
          *cell += inlet->GetPosition();
          if (flowExt)
            *cell += flowExt->normal * flowExt->length;
          insertFn(cell);
        }
        else
          break;
      }
    }

    void RBCInserter::SetShape(const MeshData & s)
    {
      shape.reset(new MeshData(s));
    }
    std::shared_ptr<const MeshData> RBCInserter::GetShape() const
    {
      return shape;
    }

    void RBCInserter::AddInLet(lb::iolets::InOutLet * inlet)
    {
      inlets.push_back(inlet);
    }
    void RBCInserter::RemoveInLet(lb::iolets::InOutLet * inlet)
    {
      inlets.erase(std::remove(std::begin(inlets), std::end(inlets), inlet), std::end(inlets));
    }

    void RBCInserter::SetScale(Dimensionless s)
    {
      scale = s;
    }
    Dimensionless RBCInserter::GetScale() const
    {
      return scale;
    }

    void RBCInserter::SetModuli(Cell::Moduli & mod)
    {
      moduli = mod;
    }
    const Cell::Moduli & RBCInserter::GetModuli() const
    {
      return moduli;
    }

    void RBCInserter::SetCondition(std::function<bool()> cond)
    {
      condition = cond;
    }
    std::function<bool()> RBCInserter::GetCondition() const
    {
      return condition;
    }

  }
}
