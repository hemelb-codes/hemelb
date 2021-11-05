# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import string
import struct
import xdrlib
from typing import List

# Here we generate a header that contains a class template test_data
# that has static member functions to return some values to be packed
# with XDR and the encoded data. We provide explicit specialisations
# for types we actually use in HemeLB.

# XDR works on 4 byte words
# Every type is coded in whole words

# For each integer type generate the values with all zero bits and
# then with each bit flipped. Same for floating types. For string, we
# try as few arbitrary things.

def indent(spaces: int, text: str) -> str:
    return text.replace("\n", "\n" + (" " * spaces))

spec_template = string.Template("""template<>
const std::vector<$typename>& test_data<$typename>::unpacked() {
  static const std::vector<$typename> unpacked = {
    $values_code
  };
  return unpacked;
}

template<>
const std::vector<char>& test_data<$typename>::packed() {
  static const std::vector<char> packed = {$buffer_data};
  return packed;
};""")

def mk_ints(nbits: int) -> List[int]:
    values = [0]
    for i in range(nbits):
        values.append(1 << i)
    return values

def mk_specialisation_int(signed: bool, nbits: int) -> str:
    p = xdrlib.Packer()
    if nbits == 32:
        if signed:
            typename = "int32_t"
            pack = p.pack_int
            suffix = ""
        else:
            typename = "uint32_t"
            pack = p.pack_uint
            suffix = "U"
    elif nbits == 64:
        if signed:
            typename = "int64_t"
            pack = p.pack_hyper
            suffix = "L"
        else:
            typename = "uint64_t"
            pack = p.pack_uhyper
            suffix = "UL"
    else:
        raise ValueError("only deal with 32/64 bit ints")
    
    values = mk_ints(nbits)
    values_code_lines = [f"0x{v:x}" for v in values]
    if signed:
        values_code_lines[-1] = "static_cast<{}>({})".format(typename, values_code_lines[-1])

    values_code = ",\n    ".join(values_code_lines)
    if signed:
        values[-1] *= -1
    for v in values:
        pack(v)

    buffer_data = ", ".join(f"'\\x{x:02x}'" for x in p.get_buffer())
    return spec_template.substitute(locals())

def mk_specialisation_float(nbits: int) -> str:
    p = xdrlib.Packer()
    if nbits == 32:
        typename = "float"
        pack = p.pack_float
        suffix = "U"
        flt_type_str = ">f"
        int_type_str = ">I"
        func_name = "binflt"
    elif nbits == 64:
        typename = "double"
        pack = p.pack_double
        suffix = "UL"
        flt_type_str = ">d"
        int_type_str = ">Q"
        func_name = "bindbl"
    else:
        raise ValueError("only deal with 32/64 bit floats")
    intvals =  mk_ints(nbits)
    fltvals = [struct.unpack(flt_type_str, struct.pack(int_type_str, v))[0] for v in intvals]
    for v in fltvals:
        pack(v)
    values_code = ",\n    ".join(f"{func_name}(0x{v:x})" for v in intvals)
    buffer_data = ", ".join(f"'\\x{x:02x}'" for x in p.get_buffer())
    return spec_template.substitute(locals())

def mk_specialisation_string() -> str:
    typename = "std::string";
    p = xdrlib.Packer()
    values = [
        "",
        "a",
        "test",
        "This is a slightly longer string whose length is a multiple of four.",
        "This is a slightly longer string whose length isn't a multiple of four.",
    ]
    for v in values:
        p.pack_string(v.encode())
    values_code = ",\n    ".join(f'"{v}"' for v in values)
    buffer_data = ", ".join(f"'\\x{x:02x}'" for x in p.get_buffer())
    return spec_template.substitute(locals())

def mk_header():
    header = """// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_IO_XDR_TEST_DATA_H
#define HEMELB_TESTS_IO_XDR_TEST_DATA_H
#include <vector>

namespace hemelb
{
  namespace tests
  {
    template <typename T>
    struct test_data {
      static const std::vector<T>& unpacked();
      static const std::vector<char>& packed();
    };
  }
}
#endif
"""

    impl_template = string.Template("""
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/io/xdr_test_data.h"

#include <cstdint>
#include <string>

namespace hemelb
{
  namespace tests
  {

    namespace {
      float binflt(uint32_t x) {
        return *reinterpret_cast<float*>(&x);
      }
      double bindbl(uint64_t x) {
        return *reinterpret_cast<double*>(&x);
      }
    }
    $specialisations

  }
}
""")
    specs = [
        mk_specialisation_int(True, 32),
        mk_specialisation_int(True, 64),
        mk_specialisation_int(False, 32),
        mk_specialisation_int(False, 64),
        
        mk_specialisation_float(32),
        mk_specialisation_float(64),

        mk_specialisation_string(),
    ]
    spec_text = "\n".join(indent(4, s) for s in specs)
    impl = impl_template.substitute(specialisations=spec_text)

    with open("xdr_test_data.h", "w") as f:
        f.write(header)
    
    with open("xdr_test_data.cc", "w") as f:
        f.write(impl)

if __name__ == "__main__":
    mk_header()
