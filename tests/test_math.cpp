#include <catch2/catch_test_macros.hpp>
#include "sim_math.h"

TEST_CASE("Basic math operations", "[math]") {
    REQUIRE(sim_math::add(2, 3) == 5);
    REQUIRE(sim_math::sub(5, 3) == 2);
    REQUIRE(sim_math::mul(4, 3) == 12);
}

TEST_CASE("Division", "[math]") {
    REQUIRE(sim_math::div(10, 2) == 5);
    REQUIRE_THROWS(sim_math::div(10, 0));
}

