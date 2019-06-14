// -*- mode: C++ -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_IO_XDRREADERTESTS_H
#define HEMELB_UNITTESTS_IO_XDRREADERTESTS_H

#include <cppunit/TestFixture.h>
#include <type_traits>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/writers/xdr/XdrReader.h"

namespace hemelb
{
  namespace unittests
  {
    namespace io
    {      
      class XdrReaderTests : public CppUnit::TestFixture
      {
	CPPUNIT_TEST_SUITE(XdrReaderTests);
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
	  xdrmem(const char* buf, size_t len) {
	    xdrmem_create(&mXdr, const_cast<char*>(buf), len, XDR_DECODE);
	  }
	  ~xdrmem() {
	    xdr_destroy(&mXdr);
	  }
	  operator XDR*() {
	    return &mXdr;
	  }
	};

	template<typename T, unsigned N = 0, typename XDRFUNC, typename GENERATOR>
	void TestBasic(XDRFUNC xdr_func, GENERATOR g) {
	  // Figure out sizes and counts
	  constexpr auto sz = sizeof(T);
	  static_assert(sz > 0U && sz <= 8U,
			"Only works on 8--64 bit types");
	  constexpr auto type_bits = sz*8U;
	  // XDR works on 32 bit words
	  constexpr auto coded_words = (sz - 1U)/4U + 1U;
	  constexpr auto coded_bytes = coded_words * 4U;
	  // Need a value for each bit being on + zero
	  constexpr auto n_vals = N ? N : type_bits + 1U;
	  // Buffers for the encoded data
	  constexpr auto buf_size = n_vals * coded_bytes;
	  char xdr_buf[buf_size];
	  char our_buf[buf_size];
	  // Fill with zeros
	  std::fill(xdr_buf, xdr_buf+buf_size, '\0');
	  std::fill(our_buf, xdr_buf+buf_size, '\0');
	  for (auto i = 0; i < n_vals; ++i) {
	    reinterpret_cast<T*>(xdr_buf)[i] = reinterpret_cast<T*>(our_buf)[i] = g(i);
	  }
	  // Make the decoders: our reference and the libC one
	  auto xdr_coder = xdrmem(xdr_buf, buf_size);
	  auto our_coder = hemelb::io::writers::xdr::XdrMemReader(our_buf, buf_size);

	  // Check they are equal
	  for (auto i = 0; i < n_vals; ++i) {
	    T xdr_val, our_val;
	    xdr_func(xdr_coder, &xdr_val);
	    our_coder.read(our_val);
	    CPPUNIT_ASSERT_EQUAL(xdr_val, our_val);
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

	template<typename FLOAT, typename UINT, typename XDRFUNC>
	void CheckFloatDecode(FLOAT val, UINT coded, XDRFUNC xdr_func) {
	  static_assert(sizeof(FLOAT) == sizeof(UINT), "size mis-match");
	  const char* buf = reinterpret_cast<char*>(&coded);
	  constexpr auto len = sizeof(UINT);
	  auto xdr_coder = xdrmem(buf, len);
	  FLOAT xdr_val;
	  xdr_func(xdr_coder, &xdr_val);
	  CPPUNIT_ASSERT_EQUAL(val, xdr_val);

	  auto our_coder = hemelb::io::writers::xdr::XdrMemReader(buf, len);
	  auto our_val = our_coder.read<FLOAT>();
	  CPPUNIT_ASSERT_EQUAL(val, our_val);
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
	  // +/-0
	  CheckFloatDecode(0.0f, 0x00000000U, xdr_float);
	  CheckFloatDecode(-0.0f, 0x00000080U, xdr_float);
	  // +/-1
	  CheckFloatDecode(1.0f, 0x0000803fU, xdr_float);
	  CheckFloatDecode(-1.0f, 0x000080bfU, xdr_float);
	  // +/- inf
	  CheckFloatDecode(std::numeric_limits<float>::infinity(), 0x0000807fU, xdr_float);
	  CheckFloatDecode(-std::numeric_limits<float>::infinity(), 0x000080ffU, xdr_float);
	  // 1.00...01 (binary)
	  CheckFloatDecode(1.0000001f, 0x0100803fU, xdr_float);	  
	}
	void TestDouble() {
	  // +/-0
	  CheckFloatDecode(0.0, 0x0000000000000000UL, xdr_double);
	  CheckFloatDecode(-0.0, 0x0000000000000080UL, xdr_double);
	  // +/-1
	  CheckFloatDecode(1.0, 0x000000000000f03fUL, xdr_double);
	  CheckFloatDecode(-1.0, 0x000000000000f0bfUL, xdr_double);
	  // +/- inf
	  CheckFloatDecode(std::numeric_limits<double>::infinity(), 0x000000000000f07fUL, xdr_double);
	  CheckFloatDecode(-std::numeric_limits<double>::infinity(), 0x000000000000f0ffUL, xdr_double);
	  // 1.00...01 (binary)
	  CheckFloatDecode(1.0000000000000002, 0x010000000000f03fUL, xdr_double);
	}

	void TestString() {
	  auto compare = [&](const char* buf, unsigned len) {
	    // Make the decoders: our reference and the libC one
	    auto xdr_coder = xdrmem(buf, len);
	    auto our_coder = hemelb::io::writers::xdr::XdrMemReader(buf, len);

	    auto xdr_val = [&]() {
	      // Remember to delete the memory allocated by XDR
	      char* xdr_raw = nullptr;
	      xdr_wrapstring(xdr_coder, &xdr_raw);
	      auto ans = std::string{xdr_raw};
	      std::free(xdr_raw);
	      return ans;
	    }();

	    auto our_val = our_coder.read<std::string>();

	    CPPUNIT_ASSERT_EQUAL(xdr_val, our_val);
	  };
	  // Null string
	  compare("\x00\x00\x00\x00",
		  4 + 0);
	  compare("\x00\x00\x00\x01"
		  "a\x00\x00\x00",
		  4 + 4);
	  compare("\x00\x00\x00\x04"
		  "test",
		  4 + 4);
	  compare("\x00\x00\x00\x44"
		  "This is a slightly longer string whose length is a multiple of four.",
		  4 + 68);
	  compare("\x00\x00\x00\x47"
		  "This is a slightly longer string whose length isn't a multiple of four.",
		  4 + 72);
	}
      private:
      };
      CPPUNIT_TEST_SUITE_REGISTRATION(XdrReaderTests);
    }
  }
}

#endif
