#include <vector>
#include <sstream>
#include <string>
#include <iterator>
#include <exception>
#include "redblood/Mesh.h"

int get_depth(std::vector<std::string> const &_args) {
  std::vector<std::string>::const_iterator i_found
    = std::find(_args.begin(), _args.end(), std::string("--depth"));
  if(i_found == _args.end())
    i_found = std::find(_args.begin(), _args.end(), std::string("-d"));
  if(i_found == _args.end()) return 0;
  if((++i_found) == _args.end()) {
    std::cerr << "Incorrect command-line argument";
    throw std::exception();
  }
  std::istringstream sstr(i_found->c_str());
  int result;
  sstr >> result;
  return result;
}

std::string get_output(std::vector<std::string> const &_args) {
  std::vector<std::string>::const_iterator i_current = _args.begin();
  std::vector<std::string>::const_iterator const i_end = _args.end();
  for(;i_current != i_end; ++i_current)
    if((*i_current)[0] != '-') return *i_current;
    else if((++i_current) == i_end) break;
  return "";
}

bool pop_option(
    std::vector<std::string> &_args, std::string const &_option) {
  std::vector<std::string>::const_iterator i_found
    = std::find(_args.begin(), _args.end(), _option);
  if(i_found == _args.end()) return false;
  _args.erase(i_found);
  return true;
}

// Usage is $ tetrahedron [filename] [-d, --depth positive integer] [--vtk]
int main(int argc, char const *argv[]) {

  std::vector<std::string> args;
  std::copy(argv+1, argv + argc, std::back_inserter(args));

  bool const dovtk = pop_option(args, "--vtk");
  int const depth(get_depth(args));
  std::string const output(get_output(args));

  hemelb::redblood::Mesh result = hemelb::redblood::tetrahedron(depth);
  if(output != "" and not dovtk)
    hemelb::redblood::write_mesh(output, *result.GetData());
  else if(output == "" and not dovtk)
    hemelb::redblood::write_mesh(std::cout, *result.GetData());
  else if(output != "")
    hemelb::redblood::write_vtkmesh(output, *result.GetData());
  else
    hemelb::redblood::write_vtkmesh(std::cout, *result.GetData());

  return 0;
}
