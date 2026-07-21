#include "Exception.hpp"
#include "PDB.hpp"
#include "Test.hpp"

#include <algorithm>
#include <string>

TEST_CASE("PDB atoms can be located by stable identity", "[pdb]") {
  // Invariant: serial and full-coordinate lookups identify the same atom.
  const zdock::PDB pdb(test::getpath("2OOB/ligand.pdb"));
  const auto &bySerial = pdb[352];
  zdock::RecordCoord coordinate;
  coordinate.serialNum = 352;
  coordinate.atomName = "N";
  coordinate.resName = "MET";
  coordinate.chain = 'b';
  coordinate.resNum = 1;

  REQUIRE(pdb[coordinate] == bySerial);
  REQUIRE(bySerial->atom.serialNum == 352);
}

TEST_CASE("PDB atom lookup rejects absent or mismatched identities", "[pdb]") {
  // Invariant: lookup never returns an atom whose requested identity is wrong.
  const zdock::PDB pdb(test::getpath("2OOB/ligand.pdb"));
  zdock::RecordCoord coordinate;
  coordinate.serialNum = 352;
  coordinate.atomName = "CA";
  coordinate.resName = "MET";
  coordinate.chain = 'b';
  coordinate.resNum = 1;

  REQUIRE_THROWS_AS(pdb[-1], zdock::AtomNotFoundException);
  REQUIRE_THROWS_AS(pdb[coordinate], zdock::AtomNotFoundException);
}

TEST_CASE("PDB filters affect atoms but preserve input records", "[pdb]") {
  // Invariant: filtering selects coordinates without discarding parsed records.
  const std::string filename = test::getpath("2OOB/ligand.pdb");
  const zdock::PDB unfiltered(filename);
  const zdock::PDB alphaCarbons(filename, [](const libpdb::PDB &record) {
    return zdock::Utils::trim_copy(record.atom.name) == "CA";
  });

  REQUIRE(alphaCarbons.atoms().size() < unfiltered.atoms().size());
  REQUIRE(alphaCarbons.records().size() == unfiltered.records().size());
  REQUIRE(alphaCarbons.matrix().cols() ==
          static_cast<long>(alphaCarbons.atoms().size()));
  for (const auto &atom : alphaCarbons.atoms()) {
    REQUIRE(zdock::Utils::trim_copy(atom->atom.name) == "CA");
  }
}

TEST_CASE("PDB value operations preserve observable coordinates", "[pdb]") {
  // Invariant: copies preserve values without sharing mutable record ownership.
  zdock::PDB original(test::getpath("2OOB/ligand.pdb"));
  zdock::PDB copied(original);
  zdock::PDB assigned;
  assigned = original;

  REQUIRE(copied.matrix().isApprox(original.matrix()));
  REQUIRE(copied.atoms().size() == original.atoms().size());
  REQUIRE(assigned.matrix().isApprox(original.matrix()));
  REQUIRE(assigned.records().size() == original.records().size());

  const double originalX = original.atoms()[0]->atom.xyz[0];
  copied.atoms()[0]->atom.xyz[0] += 10.0;
  assigned.atoms()[0]->atom.xyz[0] += 20.0;
  REQUIRE(original.atoms()[0]->atom.xyz[0] == originalX);
  REQUIRE(copied.atoms()[0] != original.atoms()[0]);
  REQUIRE(assigned.atoms()[0] != original.atoms()[0]);
}

TEST_CASE("PDB copies preserve a filter's selected record graph", "[pdb]") {
  // Invariant: copies preserve prior selections without running the filter again.
  auto firstAtomOnly = [seen = size_t{0}](const libpdb::PDB &) mutable {
    return seen++ == 0;
  };
  const zdock::PDB original(test::getpath("2OOB/ligand.pdb"), firstAtomOnly);
  const zdock::PDB copied(original);
  zdock::PDB assigned;
  assigned = original;

  REQUIRE(original.atoms().size() == 1);
  REQUIRE(copied.atoms().size() == original.atoms().size());
  REQUIRE(assigned.atoms().size() == original.atoms().size());
  REQUIRE(copied.records().size() == original.records().size());
  REQUIRE(assigned.records().size() == original.records().size());
  REQUIRE(std::find(copied.records().begin(), copied.records().end(),
                    copied.atoms()[0]) != copied.records().end());
  REQUIRE(std::find(assigned.records().begin(), assigned.records().end(),
                    assigned.atoms()[0]) != assigned.records().end());
}

TEST_CASE("PDB copies preserve model record identity", "[pdb]") {
  // Invariant: a copied model and its parent refer to the same cloned atom.
  const zdock::PDB original(test::getpath("PDB/models.pdb"));
  const zdock::PDB copied(original);

  REQUIRE(copied.nmodels() == original.nmodels());
  REQUIRE(copied.matrix().isApprox(original.matrix()));
  REQUIRE(std::find(copied.records().begin(), copied.records().end(),
                    copied.models()[0]->atoms()[0]) != copied.records().end());
  REQUIRE(copied.models()[0]->atoms()[0] != original.models()[0]->atoms()[0]);
}

TEST_CASE("PDB reports files that cannot be opened", "[pdb]") {
  // Invariant: opening a missing PDB fails with the public PDB exception type.
  REQUIRE_THROWS_AS(zdock::PDB("/path/that/does/not/exist.pdb"),
                    zdock::PDBOpenException);
}

TEST_CASE("PDB models expose and update the first model", "[pdb]") {
  // Invariant: model-aware matrix access and updates consistently target model 1.
  zdock::PDB pdb(test::getpath("PDB/models.pdb"));

  REQUIRE(pdb.nmodels() == 2);
  REQUIRE(pdb.models()[0]->modelNum() == 1);
  REQUIRE(pdb.models()[1]->modelNum() == 2);
  REQUIRE(pdb.matrix().cols() == 2);
  REQUIRE(pdb.matrix()(0, 0) == Catch::Approx(1.0));
  REQUIRE(pdb.centroid().isApprox(zdock::PDB::Coord(2.5, 3.5, 4.5)));

  zdock::PDB::Matrix updated = pdb.matrix();
  updated.array() += 5.0;
  REQUIRE(pdb.setMatrix(updated).isApprox(updated));
  REQUIRE(pdb.matrix().isApprox(updated));
  REQUIRE(pdb.models()[1]->matrix()(0, 0) == Catch::Approx(7.0));
}

TEST_CASE("PDB rejects invalid matrix and model operations", "[pdb]") {
  // Invariant: public mutations validate dimensions and model indices at runtime.
  zdock::PDB pdb(test::getpath("2OOB/ligand.pdb"));
  zdock::PDB::Matrix wrongColumns(3, pdb.matrix().cols() - 1);

  REQUIRE_THROWS_AS(pdb.setMatrix(wrongColumns), zdock::Exception);
  REQUIRE_THROWS_AS(pdb.append(pdb.atoms()[0], 1), zdock::Exception);
}

TEST_CASE("Empty PDB structures have no centroid", "[pdb]") {
  // Invariant: centroid calculation rejects an empty coordinate collection.
  const zdock::PDB empty;
  REQUIRE_THROWS_AS(empty.centroid(), zdock::Exception);
}
