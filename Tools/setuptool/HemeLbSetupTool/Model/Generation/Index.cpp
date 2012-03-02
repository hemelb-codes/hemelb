#include "Index.h"

namespace {
	void DefaultHandlerFunction(int direction) {
		throw IndexError();
	}
}

namespace hemelb {
	namespace util {
		Vector3DBase::HandlerFunction* Vector3DBase::handler = DefaultHandlerFunction;
	}
}
