#ifndef HEMELBSETUPTOOL_IOLET_H
#define HEMELBSETUPTOOL_IOLET_H

#include <vector>

struct Iolet {
  std::vector<double> Centre;
  std::vector<double> Normal;
  double Radius;
  int Id;
};

//typedef std::vector<Iolet> IoletVector;

#endif // HEMELBSETUPTOOL_IOLET_H
