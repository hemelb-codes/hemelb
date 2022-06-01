// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <string>

#include <catch2/catch.hpp>

#include "util/clone_ptr.h"

namespace hemelb::tests {

  // Set up a simple polymorphic hierarchy with `get_name` and `clone`
  // members. Need to track construction and destruction.
  namespace {
    std::vector<std::pair<std::string, std::string>> event_log;
  }

  class Base {
  public:
    Base() {
      event_log.emplace_back("Base", "default construct");
    }
    Base(Base const& b) {
      event_log.emplace_back("Base", "copy construct");
    }

    virtual ~Base() noexcept {
      event_log.emplace_back("Base", "destruct");
    }
    virtual std::string const& get_name() const = 0;
    virtual Base* clone() const = 0;
  };

  // This type all instances have the same name
  class ConstantName : public Base {
    static std::string const NAME;
  public:
    ConstantName() {
      event_log.emplace_back("ConstantName", "default construct");
    }      
    ConstantName(ConstantName const&) {
      event_log.emplace_back("ConstantName", "copy construct");
    }

    ~ConstantName() override {
      event_log.emplace_back("ConstantName", "destruct");
    }
    std::string const& get_name() const override {
      return NAME;
    }
    Base* clone() const override {
      return new ConstantName{};
    }
  };
  std::string const ConstantName::NAME = "Constant";

  // These instances have their own name
  class InstanceName : public Base {
    std::string m_name;
  public:
    InstanceName(std::string n) : m_name{n} {
      event_log.emplace_back("InstanceName", "construct");
    }
    InstanceName(InstanceName const& o) : m_name{o.m_name} {
      event_log.emplace_back("InstanceName", "copy construct");
    }
    ~InstanceName() override {
      event_log.emplace_back("InstanceName", "destruct");
    }
    std::string const& get_name() const override {
      return m_name;
    }
    Base* clone() const override {
      return new InstanceName{*this};
    }
  };

  TEST_CASE("clone_ptr_tests") {
    using namespace util;

    SECTION("Default constructible") {
      clone_ptr<Base> ptr;
      REQUIRE(ptr.get() == nullptr);
    }
    SECTION("Constructible from instances") {
      event_log.clear();
      {
	auto c = clone_ptr<Base>{new ConstantName{}};
	REQUIRE(c->get_name() == "Constant");
	REQUIRE(event_log[0].first == "Base");
	REQUIRE(event_log[1].first == "ConstantName");
	auto two = clone_ptr<Base>{new InstanceName{"two"}};
	REQUIRE(event_log[2].first == "Base");
	REQUIRE(event_log[3].first == "InstanceName");
	REQUIRE(two->get_name() == "two");
      }
      REQUIRE(event_log[4].first == "InstanceName");
      REQUIRE(event_log[5].first == "Base");
      REQUIRE(event_log[6].first == "ConstantName");
      REQUIRE(event_log[7].first == "Base");
      REQUIRE(event_log.size() == 8);
    }

    SECTION("Cloning on constant type works") {
      event_log.clear();
      {
	auto two = clone_ptr<Base>{new InstanceName{"two"}};
	REQUIRE(event_log[0].first == "Base");
	REQUIRE(event_log[1].first == "InstanceName");
	auto two_ptr = two->clone();
	REQUIRE(event_log[2].first == "Base");
	REQUIRE(event_log[3].first == "InstanceName");
	static_assert(std::is_same_v<decltype(two_ptr), Base*>);
	REQUIRE(two_ptr->get_name() == "two");
	delete two_ptr;
	REQUIRE(event_log[4].first == "InstanceName");
	REQUIRE(event_log[5].first == "Base");
      }
      REQUIRE(event_log[6].first == "InstanceName");
      REQUIRE(event_log[7].first == "Base");
      REQUIRE(event_log.size() == 8);
    }

    SECTION("Move construct works, including up hierarchy") {
      event_log.clear();
      {
	auto three = clone_ptr<InstanceName>{new InstanceName{"three"}};
	REQUIRE(event_log[0].first == "Base");
	REQUIRE(event_log[1].first == "InstanceName");
	auto also_three = clone_ptr<InstanceName>{std::move(three)};
	INFO("Must null the moved-from object");
	REQUIRE(!three);
	REQUIRE(event_log.size() == 2);
	clone_ptr<Base> base3{std::move(also_three)};
	REQUIRE(!also_three);
	REQUIRE(event_log.size() == 2);

	auto must_fail = clone_dynamic_cast<ConstantName>(std::move(base3));
	REQUIRE(!must_fail);
	auto derived3 = clone_dynamic_cast<InstanceName>(std::move(base3));
	static_assert(std::is_same_v<decltype(derived3), clone_ptr<InstanceName>>);
      }
      REQUIRE(event_log[2].first == "InstanceName");
      REQUIRE(event_log[3].first == "Base");
      REQUIRE(event_log.size() == 4);
    }

    SECTION("Move assign works, including up hierarchy") {
      event_log.clear();
      {
	auto four = clone_ptr<InstanceName>{new InstanceName{"four"}};
	REQUIRE(event_log[0].first == "Base");
	REQUIRE(event_log[1].first == "InstanceName");
	clone_ptr<InstanceName> also_four;
	also_four = std::move(four);
	INFO("Must null the moved-from object");
	REQUIRE(!four);
	REQUIRE(event_log.size() == 2);

	clone_ptr<Base> base4;
	base4 = std::move(also_four);
	REQUIRE(event_log.size() == 2);

	clone_ptr<InstanceName> derived4;
	derived4 = clone_dynamic_cast<InstanceName>(std::move(base4));
	REQUIRE(event_log.size() == 2);

	// This must destruct
	derived4 = nullptr;
	REQUIRE(event_log[2].first == "InstanceName");
	REQUIRE(event_log[3].first == "Base");
	REQUIRE(event_log.size() == 4);
	REQUIRE(!derived4);
      }
    }

    SECTION("Copy construct clones") {
      event_log.clear();
      {
	auto const five = clone_ptr<InstanceName>{new InstanceName{"five"}};
	REQUIRE(event_log[0].first == "Base");
	REQUIRE(event_log[1].first == "InstanceName");
	REQUIRE(event_log.size() == 2);
	auto const cp5 = clone_ptr<InstanceName>{five};
	REQUIRE(event_log[2].first == "Base");
	REQUIRE(event_log[3].first == "InstanceName");
	REQUIRE(event_log.size() == 4);

	clone_ptr<Base> const b5{cp5};
	REQUIRE(event_log[4].first == "Base");
	REQUIRE(event_log[5].first == "InstanceName");
	REQUIRE(event_log.size() == 6);
	REQUIRE(b5->get_name() == "five");

	clone_ptr<InstanceName> const d5{clone_dynamic_cast<InstanceName>(b5)};
	REQUIRE(event_log.size() == 8);
      }
      REQUIRE(event_log.size() == 16);
    }

    SECTION("Copy assign clones") {
      event_log.clear();
      {
	auto six_1 = clone_ptr<ConstantName>{new ConstantName};
	REQUIRE(event_log[0].first == "Base");
	REQUIRE(event_log[1].first == "ConstantName");
	REQUIRE(event_log.size() == 2);

	clone_ptr<ConstantName> six_2;
	six_2 = six_1;
	REQUIRE(event_log.size() == 4);

	clone_ptr<Base> six_3;
	six_3 = six_2;
	REQUIRE(event_log.size() == 6);
      }
      REQUIRE(event_log.size() == 12);
    }

    SECTION("make_clone_ptr works") {
      event_log.clear();
      auto c = make_clone_ptr<ConstantName>();
      static_assert(std::is_same_v<decltype(c), clone_ptr<ConstantName>>);

      clone_ptr<Base> b = make_clone_ptr<InstanceName>("wibble");
      REQUIRE(event_log.size() == 4);
    }

    SECTION("clone member works") {
      event_log.clear();
      clone_ptr<Base> b = make_clone_ptr<InstanceName>("wibble");
      {
	auto copy = b.clone();
	static_assert(std::is_same_v<decltype(copy), clone_ptr<Base>>);
	REQUIRE(event_log.size() == 4);
	auto derived_cp = clone_dynamic_cast<InstanceName>(std::move(copy));
	REQUIRE(!copy);
	static_assert(std::is_same_v<decltype(derived_cp), clone_ptr<InstanceName>>);
	REQUIRE(event_log.size() == 4);
      }
      
    }
    
  }
}
