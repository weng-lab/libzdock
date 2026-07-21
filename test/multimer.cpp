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
