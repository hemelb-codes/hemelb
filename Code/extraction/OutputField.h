// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_OUTPUTFIELD_H
#define HEMELB_EXTRACTION_OUTPUTFIELD_H

#include <variant>

#include "Exception.h"
#include "io/formats/extraction.h"

namespace hemelb::extraction
{
  namespace detail {
    // Helpers for allowing one to pass a number of lambdas to visit a
    // std:variant (see visit functions in namespaces below).
    template<typename... Ts> struct make_overload: Ts... { using Ts::operator()...; };
    template<typename... Ts> make_overload(Ts...) -> make_overload<Ts...>;
  }

  // Namespace holding tag types and variant for the data source.
  namespace source {
    // Tag types
    struct Pressure {};
    struct Velocity {};
    struct ShearStress {};
    struct VonMisesStress {};
    struct ShearRate {};
    struct StressTensor {};
    struct Traction {};
    struct TangentialProjectionTraction {};
    struct Distributions {};
    struct MpiRank {};

    // Variant of tags. This is a finite sum type so can get compiler
    // to exhaustively check for use (much neater than an enum!).
    using Type = std::variant<
      Pressure,
      Velocity,
      ShearStress,
      VonMisesStress,
      ShearRate,
      StressTensor,
      Traction,
      TangentialProjectionTraction,
      Distributions,
      MpiRank
    >;

    // Helper to visit the variant with a bunch of lamdbas. E.g.:
    //
    // source::visit(src_var,
    //   [](Pressure unused) {/* pressure specific stuff */},
    //   [](Velocity) { /* velocity specifi stuff */},
    //   ... callables for other possible types ...
    // );
    template <typename... Callables>
    auto visit(Type const& t, Callables&&... calls) {
      return std::visit(detail::make_overload{std::forward<Callables>(calls)...}, t);
    }
  }

  // Namespace holding variant and helpers for the type to be saved to
  // the file.
  namespace code {
    using io::formats::extraction::TypeCode;
    using Type = std::variant<
      float,
      double,
      std::int32_t,
      std::uint32_t,
      std::int64_t,
      std::uint64_t
    >;

    // Visit helper func (see source::visit above)
    template <typename... Callables>
    auto visit(Type const& t, Callables&&... calls) {
      return std::visit(detail::make_overload{std::forward<Callables>(calls)...}, t);
    }

    // Convert the enum to the variant with the appropriate type
    // active.
    inline Type enum_to_type(TypeCode tc) {
      switch(tc) {
      case TypeCode::FLOAT:
	return 0.0f;
      case TypeCode::DOUBLE:
	return 0.0;
      case TypeCode::INT32:
	return std::int32_t{0};
      case TypeCode::UINT32:
	return std::uint32_t{0};
      case TypeCode::INT64:
	return std::int64_t{0};
      case TypeCode::UINT64:
	return std::uint64_t{0};
      default:
	throw Exception() << "Invalid type";
      }
    }

    // Convert the variant to an enum
    inline TypeCode type_to_enum(Type const& tvar) {
      return visit(tvar,
	[](float) { return TypeCode::FLOAT; },
	[](double) { return TypeCode::DOUBLE; },
	[](std::int32_t) { return TypeCode::INT32; },
	[](std::uint32_t) { return TypeCode::UINT32; },
	[](std::int64_t) { return TypeCode::INT64; },
	[](std::uint64_t) { return TypeCode::UINT64; }
      );
    }

    // Compute sizeof the active member of the variant.
    inline std::size_t type_to_size(Type const& tvar) {
      return std::visit(
	[](auto&& t) { return sizeof(decltype(t)); },
	tvar
      );
    }

  }

  struct OutputField
  {
    std::string name;
    source::Type src;
    code::Type typecode;
    std::uint32_t noffsets;
    // Data sources are double so do subtraction at full precision
    // before converting.
    std::vector<double> offset;
  };
}

#endif /* HEMELB_EXTRACTION_OUTPUTFIELD_H */
