#include "Debug.h"

namespace {
	DummyStream ds;
}

DummyStream& Log() {
	return ds;
}
