// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

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

TEST_CASE("C2 multimer subunits are related by a twofold axis", "[multimer]") {
  // Invariant: for symmetry two, the rigid operator mapping subunit A onto
  // subunit B is a 180 degree rotation: symmetric, self-inverse, trace -1.
  const zdock::ZDOCK docking(test::getpath("ZDOCK/mzdock_c2.out"));
  const zdock::TransformMultimer transform(docking);
  zdock::PDB::Matrix input(3, 4);
  input << 0.0, 4.0, 0.0, 1.0,
           0.0, 0.0, 3.0, 2.0,
           0.0, 0.0, 0.0, 5.0;

  const auto a = transform.txMultimer(input, docking.predictions()[0], 0);
  const auto b = transform.txMultimer(input, docking.predictions()[0], 1);

  // Solve b - b0 = R (a - a0) exactly from three independent edge vectors.
  Eigen::Matrix3d edgesA;
  Eigen::Matrix3d edgesB;
  for (int i = 0; i < 3; ++i) {
    edgesA.col(i) = a.col(i + 1) - a.col(0);
    edgesB.col(i) = b.col(i + 1) - b.col(0);
  }
  const Eigen::Matrix3d rotation = edgesB * edgesA.inverse();

  REQUIRE(rotation.isApprox(rotation.transpose(), 1e-9));
  REQUIRE((rotation * rotation).isApprox(Eigen::Matrix3d::Identity(), 1e-9));
  REQUIRE(rotation.trace() == Catch::Approx(-1.0).margin(1e-9));
  REQUIRE(rotation.determinant() == Catch::Approx(1.0).margin(1e-9));
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
