#include "Eigen/Dense"
#include "Test.hpp"
#include "TransformMultimer.hpp"
#include "ZDOCK.hpp"

TEST_CASE("M-ZDOCK transformations are rigid", "[multimer]") {
  // Invariant: every generated subunit preserves all pairwise atom distances.
  const zdock::ZDOCK docking(test::getpath("ZDOCK/mzdock.out"));
  const zdock::TransformMultimer transform(docking);
  zdock::PDB::Matrix input(3, 3);
  input << 0.0, 1.0, 0.0,
           0.0, 0.0, 2.0,
           0.0, 0.0, 0.0;

  const auto first = transform.txMultimer(input, docking.predictions()[0], 0);
  const auto second = transform.txMultimer(input, docking.predictions()[0], 1);

  zdock::PDB::Matrix expected(3, 3);
  expected << -19.855310217231672, -20.543392706773322,
      -20.711382348244648, 17.358073888770218, 17.338265216204988,
      18.994614060996650, 0.666941280683854, 1.392303331026556,
      -0.100441786103513;

  // Invariant: the documented M-ZDOCK rotations and translations produce this
  // independently composed reference pose, including their order and direction.
  REQUIRE(first.isApprox(expected, 1e-12));

  for (int left = 0; left < input.cols(); ++left) {
    for (int right = 0; right < input.cols(); ++right) {
      REQUIRE((first.col(left) - first.col(right)).norm() ==
              Catch::Approx((input.col(left) - input.col(right)).norm()));
      REQUIRE((second.col(left) - second.col(right)).norm() ==
              Catch::Approx((input.col(left) - input.col(right)).norm()));
    }
  }
  REQUIRE_FALSE(first.isApprox(second));
}

TEST_CASE("TransformMultimer accepts an M-ZDOCK filename", "[multimer]") {
  // Invariant: filename and parsed-object constructors produce identical poses.
  const std::string filename = test::getpath("ZDOCK/mzdock.out");
  const zdock::ZDOCK docking(filename);
  const zdock::TransformMultimer fromFile(filename);
  const zdock::TransformMultimer fromObject(docking);
  zdock::PDB::Matrix input(3, 1);
  input << 1.0, 2.0, 3.0;

  REQUIRE(fromFile.txMultimer(input, docking.predictions()[1], 2).isApprox(
      fromObject.txMultimer(input, docking.predictions()[1], 2)));
}

TEST_CASE("Multimer chain identifiers cover the complete PDB alphabet",
          "[multimer]") {
  // Invariant: each of the 52 supported components receives a unique chain ID.
  std::string chains;
  for (const auto &chain : zdock::TransformMultimer::CHAINS) {
    chains += chain;
  }

  REQUIRE(chains == "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
}

TEST_CASE("Multimer transformations reject incompatible requests",
          "[multimer]") {
  // Invariant: wrong formats and components outside symmetry fail explicitly.
  const zdock::ZDOCK regular(test::getpath("ZDOCK/6GWC.zd.out"));
  const zdock::ZDOCK multimer(test::getpath("ZDOCK/mzdock.out"));
  const zdock::TransformMultimer wrongFormat(regular);
  const zdock::TransformMultimer transform(multimer);
  zdock::PDB::Matrix input(3, 1);
  input << 1.0, 2.0, 3.0;

  REQUIRE_THROWS_AS(
      wrongFormat.txMultimer(input, regular.predictions()[0], 0),
      zdock::ZDOCKUnsupported);
  REQUIRE_THROWS_AS(
      transform.txMultimer(input, multimer.predictions()[0], -1),
      zdock::ZDOCKUnsupported);
  REQUIRE_THROWS_AS(
      transform.txMultimer(input, multimer.predictions()[0], 24),
      zdock::ZDOCKUnsupported);
}
