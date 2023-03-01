// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "util/Vector3D.h"


namespace hemelb::tests
{
    // Aliases to save typing
    using VD = util::Vector3D<double>;
    using VF = util::Vector3D<float>;
    using VI = util::Vector3D<int>;
    using VU = util::Vector3D<unsigned>;
    using VI8 = util::Vector3D<std::int8_t>;
    using VI16 = util::Vector3D<std::int16_t>;
    using VI32 = util::Vector3D<std::int32_t>;
    using VI64 = util::Vector3D<std::int64_t>;
    using VU8 = util::Vector3D<std::uint8_t>;
    using VU16 = util::Vector3D<std::uint16_t>;
    using VU32 = util::Vector3D<std::uint32_t>;
    using VU64 = util::Vector3D<std::uint64_t>;

    TEST_CASE("Vector3D construction", "[util]") {
        SECTION("default") {
            VD a;
            a[0] = 0;
        }

        SECTION("three value") {
            auto a = VI{1, 2, 3};
            REQUIRE(*a.begin() == 1);
            REQUIRE(a.y() == 2);
            REQUIRE(a[2] == 3);
        }

        SECTION("one value") {
            auto a = VI8{'a'};
            REQUIRE(a[0] == 97);
            REQUIRE(a[1] == 97);
            REQUIRE(a[2] == 97);
        }

        SECTION("zero") {
            auto a = VI32::Zero();
            REQUIRE(a[0] == 0);
            REQUIRE(a[1] == 0);
            REQUIRE(a[2] == 0);
        }

        SECTION("one") {
            auto a = VU32::Ones();
            REQUIRE(a[0] == 1U);
            REQUIRE(a[1] == 1U);
            REQUIRE(a[2] == 1U);
        }

        SECTION("min/max") {
            auto big = VI8::Largest();
            REQUIRE(big[0] == 127);
            auto small = VI8::Lowest();
            REQUIRE(small[0] == -128);
        }

        SECTION("copy") {
            auto orig = VU16::Ones();
            {
                VU16 copy(orig);
                REQUIRE(*copy.begin() == 1U);
                REQUIRE(copy.y() == 1U);
                REQUIRE(copy[2] == 1U);
            }
            {
                VI64 copy(orig);
                REQUIRE(*copy.begin() == 1U);
                REQUIRE(copy.y() == 1U);
                REQUIRE(copy[2] == 1U);
            }
        }

        SECTION("move") {
            auto src = VI16{-1};
            VI16 dest(std::move(src));
            REQUIRE(dest[0] == -1);
        }

        SECTION("as") {
            auto vec = VU8::Largest().as<int>();
            STATIC_REQUIRE(std::is_same_v<decltype(vec), VI>);
            REQUIRE(vec[1] == 255);
        }
    }

    TEST_CASE("Vector3D in place arithmetic", "[util]") {
        {
            // language loves turning things into ints
            VI sum(0);
            sum += VI::Ones();
            sum += VI8{2};
            REQUIRE(sum[0] == 3);
        }
        {
            VU8 sum(0U);
            sum += VU::Ones();
            sum += VU8{2};
            REQUIRE(sum[0] == 3);
        }

        {
            auto v = VI::Ones();
            v *= 2;
            REQUIRE(v[2] == 2);
        }

        {
            auto v = VU64{2, 4, 8};
            v /= 2;
            REQUIRE(v == VU64{1, 2, 4});
        }
        {
            auto v = VU{4, 5, 6};
            v %= 5;
            REQUIRE(v == VU{4, 0, 1});
        }
    }

#define TYPECHECK(TYPE, expr) STATIC_REQUIRE(std::same_as<TYPE, decltype(expr)>)

    TEST_CASE("Vector3D arithmetic", "[util]") {
        auto a = VI::Ones() + VI::Zero();
        TYPECHECK(VI, a);
        TYPECHECK(VI, VI::Ones() + VI16::Ones());
        TYPECHECK(VI, VI16::Ones() + VI::Ones());

        SECTION("addition") {
            REQUIRE(a == VI::Ones());

            a += VI::Ones();
            REQUIRE(a == VI{2, 2, 2});

            auto b = a + VD::Ones();
            TYPECHECK(VD, b);
            REQUIRE(b == VD{3, 3, 3});

            auto c = VD::Ones() + b;
            TYPECHECK(VD, c);
            REQUIRE(c == VD{4, 4, 4});
        }

        SECTION("subtraction") {
            a -= VI::Ones();
            REQUIRE(a == VI::Zero());

            auto b = VI32::Largest() - VI32{0x7ffffffe};
            TYPECHECK(VI32, b);
            REQUIRE(b == VI32::Ones());
        }

        SECTION("multiplication by scalar") {
            a *= 0;
            REQUIRE(a == VI::Zero());

            auto b = VU::Ones() * 2;
            TYPECHECK(VU, b);
            REQUIRE(b == VU{2, 2, 2});

            auto c = 3 * VU::Ones();
            TYPECHECK(VU, c);
            REQUIRE(c == VU{3, 3, 3});
        }

        SECTION("dot product") {
            auto b = VU32{1, 2, 3};

            auto b2 = Dot(b, b);
            TYPECHECK(std::uint32_t, b2);
            REQUIRE(b2 == (1 + 4 + 9));

            auto c = -VI64::Ones();
            TYPECHECK(VI64, c);
            REQUIRE(c == VI64{-1, -1, -1});

            auto bc = Dot(b, c);
            TYPECHECK(std::int64_t, bc);
            auto cb = Dot(c, b);
            TYPECHECK(std::int64_t, cb);
            REQUIRE(bc == cb);
            REQUIRE(bc == -6);
        }

        SECTION("cross product") {
            auto x = VI{1, 0, 0};
            auto y = VI{0, 1, 0};
            auto z = VI{0, 0, 1};

            auto x_y = Cross(x, y);
            TYPECHECK(VI, x_y);
            REQUIRE(x_y == z);

            REQUIRE(Cross(y, x) == VI{0,0, -1});
        }
    }

    TEST_CASE("Vector3D CastsInVector3DProduct", "[util]") {
	const double dblMax = std::numeric_limits<double>::max();
	const unsigned uintMax = std::numeric_limits<unsigned>::max();

	SECTION("(int, float) -> float") {
	  VI foo(-1, 0, 1);
	  float bar = 0.8;
	  VF baz = foo * bar;

	  REQUIRE(-0.8f == baz[0]);
	  REQUIRE(0.0f == baz[1]);
	  REQUIRE(0.8f == baz[2]);
	}

	SECTION("(float, double) -> double"){
	  VF foo(0.0f, 1.0f, -1.0f);
	  double bar = dblMax;
	  VD baz = foo * bar;

	  REQUIRE(0.0 == baz[0]);
	  REQUIRE(dblMax == baz[1]);
	  REQUIRE(-dblMax == baz[2]);
	}

	SECTION("(int, unsigned) -> unsigned") {
	  VI foo(0, 2, 2);
	  unsigned bar = uintMax / 2u;
	  VU baz = foo * bar;

	 REQUIRE(0u == baz[0]);
	 REQUIRE(uintMax == baz[1] + uintMax % 2);
	 REQUIRE(uintMax == baz[2] + uintMax % 2);
	}
    }

    TEST_CASE("Vector3D works at compile time", "[util]") {
      constexpr auto z = VI::Zero();
      STATIC_REQUIRE(z.x() == 0);
      STATIC_REQUIRE(z.y() == 0);
      STATIC_REQUIRE(z.z() == 0);

      constexpr auto one = VI::Ones();
      STATIC_REQUIRE(one.x() == 1);
      STATIC_REQUIRE(one.y() == 1);
      STATIC_REQUIRE(one.z() == 1);

    }
}
