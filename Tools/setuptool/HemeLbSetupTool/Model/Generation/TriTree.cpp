#include "TriTree.h"

std::ostream& operator<<(std::ostream& os, const IdList& lst) {
  bool first = true;
  os << '[';
  for (auto i: lst) {
    if (!first)
      os << ',';
    os << i;
    first = false;
  }
  os << ']';
  return os;
}
