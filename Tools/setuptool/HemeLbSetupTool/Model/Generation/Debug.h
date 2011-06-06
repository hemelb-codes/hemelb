#ifndef HEMELBSETUPTOOL_DEBUG_H
#define HEMELBSETUPTOOL_DEBUG_H

#include <iostream>

class DummyStream {
};

template<typename T>
DummyStream& operator<<(DummyStream &ds, const T& val) {
	return ds;
}
inline DummyStream& operator<<(DummyStream& ds, std::ostream& (*func) ( std::ostream& os )) {
	return ds;
}

//#define Log() std::cout
DummyStream& Log();

#endif // HEMELBSETUPTOOL_DEBUG_H
//DummyStream() << "Domain size " << this->BlockCounts << std::endl;
