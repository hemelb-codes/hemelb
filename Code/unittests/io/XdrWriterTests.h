// -*- mode: C++ -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_IO_XDRWRITERTESTS_H
#define HEMELB_UNITTESTS_IO_XDRWRITERTESTS_H

#include <cppunit/TestFixture.h>
#include <type_traits>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/writers/xdr/XdrWriter.h"

namespace hemelb
{
  namespace unittests
  {
    namespace io
    {
      
      class XdrWriterTests : public CppUnit::TestFixture
      {
	CPPUNIT_TEST_SUITE(XdrWriterTests);
	CPPUNIT_TEST(TestInt16);
	CPPUNIT_TEST(TestUInt16);
	CPPUNIT_TEST(TestInt32);
	CPPUNIT_TEST(TestUInt32);
	CPPUNIT_TEST(TestInt64);
	CPPUNIT_TEST(TestUInt64);

	CPPUNIT_TEST(TestFloat);
	CPPUNIT_TEST(TestDouble);

	CPPUNIT_TEST(TestString);
	CPPUNIT_TEST_SUITE_END();

	using buf_t = std::vector<char>;
	class xdrmem {
	  XDR mXdr;
	public:
	  xdrmem(char* buf, size_t len) {
	    xdrmem_create(&mXdr, buf, len, XDR_ENCODE);
	  }
	  ~xdrmem() {
	    xdr_destroy(&mXdr);
	  }
	  operator XDR*() {
	    return &mXdr;
	  }
	};

	template<typename T, typename XDRFUNC, typename GENERATOR>
	void TestBasic(XDRFUNC xdr_func, GENERATOR g) {
	  // Figure out sizes and counts
	  constexpr auto sz = sizeof(T);
	  static_assert(sz > 0 && sz <= 8,
			"Only works on 8--64 bit types");
	  constexpr auto type_bits = sz*8;
	  // XDR works on 32 bit words
	  constexpr auto coded_words = (sz - 1)/4 + 1;
	  constexpr auto coded_bytes = coded_words * 4;
	  // Need a value for each bit being on + zero
	  constexpr auto n_vals = type_bits + 1;
	  // Buffers for the encoded data
	  constexpr auto buf_size = n_vals * coded_bytes;
	  char xdr_buf[buf_size];
	  char our_buf[buf_size];
	  // Fill with binary ones to trigger failure if we're not
	  // writing enough zeros
	  std::fill(xdr_buf, xdr_buf+buf_size, ~'\0');
	  std::fill(our_buf, xdr_buf+buf_size, ~'\0');
	  // Type for comparing as a single operation
	  using CUINT = typename std::conditional<coded_words == 1, uint32_t, uint64_t>::type;
	  const CUINT* xdr_check_buf = reinterpret_cast<CUINT*>(xdr_buf);
	  const CUINT* our_check_buf = reinterpret_cast<CUINT*>(our_buf);

	  // Make the encoders: our reference and the libC one
	  auto xdr_coder = xdrmem(xdr_buf, buf_size);
	  auto our_coder = hemelb::io::MakeXdrWriter(our_buf, our_buf + buf_size);

	  // Check every bit in turn to ensure endian issues are dealt with
	  for (size_t i_bit = 0; i_bit <= type_bits; ++i_bit) {
	    T val = g(i_bit);
	    xdr_func(xdr_coder, &val);
	    our_coder << val;
	    CPPUNIT_ASSERT_EQUAL(xdr_check_buf[i_bit], our_check_buf[i_bit]);
	  }
	}

	// Flip every bit in the range in turn and check it encodes to
	// the same as libc's xdr_encode. Works only for 1- or 2-word
	// encoded size types.
	template<typename INT, typename XDRFUNC>
	void TestInt(XDRFUNC xdr_func) {
	  static_assert(std::is_integral<INT>::value, "only works on (u)ints");
	  auto gen = [](size_t i_bit) -> INT {
	    return i_bit ? (INT(1) << INT(i_bit-1)) : INT(0);
	  };
	  CPPUNIT_ASSERT_EQUAL(INT(0), gen(0));
	  CPPUNIT_ASSERT_EQUAL(INT(1), gen(1));
	  CPPUNIT_ASSERT_EQUAL(INT(2), gen(2));
	  CPPUNIT_ASSERT_EQUAL(INT(4), gen(3));
	  CPPUNIT_ASSERT_EQUAL(INT(8), gen(4));
	  TestBasic<INT>(xdr_func, gen);
	}

	template<typename FLOAT, typename XDRFUNC>
	void TestFloating(XDRFUNC xdr_func) {
	  static_assert(std::is_floating_point<FLOAT>::value, "Floats only please!");
	  constexpr bool is_double = sizeof(FLOAT) == 8;
	  
	  auto gen = [&](size_t i_bit) -> FLOAT {
	    using INT = typename std::conditional<is_double, uint64_t, uint32_t>::type;
	    INT data = i_bit ? (INT(1) << INT(i_bit-1)) : INT(0);
	    auto ans_p = reinterpret_cast<FLOAT*>(&data);
	    return *ans_p;
	  };

	  CPPUNIT_ASSERT_EQUAL(FLOAT(0.0), gen(0));

	  constexpr auto sign_bit = sizeof(FLOAT)*8 - 1;
	  CPPUNIT_ASSERT_EQUAL(-FLOAT(0.0), gen(sign_bit + 1));

	  CPPUNIT_ASSERT_EQUAL(std::numeric_limits<FLOAT>::denorm_min(), gen(1));

	  constexpr size_t mantissa_bits = is_double ? 52 : 23;
	  CPPUNIT_ASSERT_EQUAL(std::numeric_limits<FLOAT>::min(), gen(mantissa_bits+1));
	  TestBasic<FLOAT>(xdr_func, gen);
	}

      public:
	void TestInt16() {
	  TestInt<int16_t>(xdr_int16_t);
	}
	void TestUInt16() {
	  TestInt<uint16_t>(xdr_u_int16_t);
	}
	void TestInt32() {
	  TestInt<int32_t>(xdr_int);
	}
	void TestUInt32() {
	  TestInt<uint32_t>(xdr_u_int32_t);
	}
	void TestInt64() {
	  TestInt<int64_t>(xdr_int64_t);
	}
	void TestUInt64() {
	  TestInt<uint64_t>(xdr_u_int64_t);
	}

	void TestFloat() {
	  TestFloating<float>(xdr_float);
	}
	void TestDouble() {
	  TestFloating<double>(xdr_double);
	}

	void TestString() {
	  using UPC = std::unique_ptr<char[]>;
	  // XDR strings are serialised as length, data (0-padded to a word)
	  auto coded_length = [](const std::string& s) {
	    return s.size() ? 4U * (2U + (s.size() - 1U) / 4U) : 4U;
	  };
	  auto compare = [&](const std::string& s) {
	    auto n = coded_length(s);
	    auto make_ones = [&](size_t n) -> UPC {
	      // Fill with binary ones to trigger failure if we're not
	      // writing enough zeros
	      auto ans = UPC(new char[n]);
	      std::fill(ans.get(), ans.get() + n, ~'\0');
	      return ans;
	    };

	    auto ref = make_ones(n);
	    auto test = make_ones(n);
	    // Make the encoders: our reference and the libC one
	    auto ref_coder = xdrmem(ref.get(), n);
	    auto test_coder = hemelb::io::MakeXdrWriter(test.get(), test.get() + n);

	    char* tmp = const_cast<char*>(s.c_str());
	    xdr_string(ref_coder, &tmp, s.size());
	    test_coder << s;

	    for (int i = 0; i < n; ++i) {
	      CPPUNIT_ASSERT_EQUAL(ref[i], test[i]);
	    }
	  };
	  compare("");
	  compare("a");
	  compare("test");
	  compare("This is a slightly longer string whose length is a multiple of four.");
	  compare("This is a slightly longer string whose length isn't a multiple of four.");
	}
      private:
      };
      CPPUNIT_TEST_SUITE_REGISTRATION(XdrWriterTests);
    }
  }
}

#endif
