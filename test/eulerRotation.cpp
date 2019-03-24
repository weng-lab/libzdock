#include "Eigen/Dense"
#include "TransformUtil.hpp"
#include "Test.hpp"

TEST_CASE("TransformUtil basic checks", "[eulerRotation]") {
  Eigen::Matrix<double, 4, 4> x, y;
  double epsilon = 1e-30;
  SECTION("eulerRotation 1/7") {
    // expected result
    x << 0.9716366168288380, 0.114650906882496, -0.2068271123146840, 0,
        0.2259211030107070, -0.708457854264245, 0.6686158268734950, 0,
        -0.0698708812870948, -0.696378229328232, -0.7142656520272000, 0, 0, 0,
        0, 1;
    y = zdock::TransformUtil::eulerRotation({0.3, 10.2, 0.1}).matrix();
    REQUIRE((x - y).squaredNorm() < epsilon);
  }
  SECTION("eulerRotation 2/7") {
    // expected result
    x << 0.4709650619398320, -0.412444527669252, 0.7797957566104720, 0,
        -0.0970453222046322, 0.854380927745508, 0.5105050790569320, 0,
        -0.8767976481892560, -0.316105586632715, 0.3623577544766740, 0, 0, 0, 0,
        1;
    y = zdock::TransformUtil::eulerRotation({8.2, -1.2, -21}, true).matrix();
    REQUIRE((x - y).squaredNorm() < epsilon);
  }
  SECTION("eulerRotation 3/7") {
    // expected result
    x << 0.9716366168288380, 0.225921103010707, -0.0698708812870948, 0,
        0.1146509068824960, -0.708457854264245, -0.6963782293282320, 0,
        -0.2068271123146840, 0.668615826873495, -0.7142656520272000, 0, 0, 0, 0,
        1;
    y = zdock::TransformUtil::eulerRotation({0.3, 10.2, 0.1}, true).matrix();
    REQUIRE((x - y).squaredNorm() < epsilon);
  }
  SECTION("eulerRotation 4/7") {
    // expected result
    x << 1, 0, 0, 0, 0, 0.540302305868140, -0.8414709848078970, 0, 0,
        0.841470984807897, 0.5403023058681400, 0, 0, 0, 0, 1;
    y = zdock::TransformUtil::eulerRotation({0, 1, 0}).matrix();
    REQUIRE((x - y).squaredNorm() < epsilon);
  }
  SECTION("eulerRotation 5/7") {
    // expected result
    x << -1, -1.22464679914735e-16, 0, 0, 1.22464679914735e-16, -1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 1;
    y = zdock::TransformUtil::eulerRotation({zdock::TransformUtil::PI, 0, 0})
            .matrix();
    REQUIRE((x - y).squaredNorm() < epsilon);
  }
  SECTION("eulerRotation 6/7") {
    // expected result
    x << 1, 0, 0, 0, 0, -1, -1.22464679914735e-16, 0, 0, 1.22464679914735e-16,
        -1, 0, 0, 0, 0, 1;
    y = zdock::TransformUtil::eulerRotation({0, zdock::TransformUtil::PI, 0})
            .matrix();
    REQUIRE((x - y).squaredNorm() < epsilon);
  }
  SECTION("eulerRotation 7/7") {
    // expected result
    x << -1, -1.22464679914735e-16, 0, 0, 1.22464679914735e-16, -1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 1;
    y = zdock::TransformUtil::eulerRotation({0, 0, zdock::TransformUtil::PI})
            .matrix();
    REQUIRE((x - y).squaredNorm() < epsilon);
  }
}
