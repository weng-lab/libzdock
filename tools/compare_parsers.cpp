// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

// Legacy-vs-modern parser comparison harness. Walks a file or a directory
// tree of PDB files and pushes every line through both parsers, comparing
// the record type and the rendered text byte-for-byte. The PDBRUN version
// state is mirrored between the two parsers line by line: the legacy
// parser keeps it in process globals (flipped by USER PDBRUN lines and by
// rendering USER PDBRUN records), the modern parser carries it in
// ParseContext/WriteContext, and the harness synchronizes the two after
// every line.
//
// Build (or use the compare-parsers make target):
//   g++ -O2 -std=c++17 -Icontrib/eigen -Iinclude -Isrc/libpdb++
//       -Isrc/record tools/compare_parsers.cpp lib/libzdock.a
//       -o bin/compare_parsers
//
// Usage:
//   bin/compare_parsers <file-or-directory> [--max-files N] [--max-lines N]
//                        [--all-extensions] [--quiet]
//
// Exit status: 0 = full parity, 1 = mismatches found, 2 = usage or I/O
// error. Files are visited in sorted order for reproducibility; only
// *.pdb / *.PDB files are visited unless --all-extensions is given.

#include "PDB.hpp"
#include "Record.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

struct Options {
  std::string root;
  std::size_t maxFiles = 0; // 0 = unlimited
  std::size_t maxLines = 0; // 0 = unlimited
  bool allExtensions = false;
  bool quiet = false;
};

struct Counters {
  std::size_t files = 0;
  std::size_t lines = 0;
  std::size_t typeMismatches = 0;
  std::size_t textMismatches = 0;
  std::size_t rejections = 0; // new parser returned nullopt
};

bool pdbExtension(const std::filesystem::path &path) {
  std::string ext = path.extension().string();
  for (char &c : ext)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return ext == ".pdb";
}

std::string shorten(const std::string &text) {
  if (text.size() <= 72)
    return text;
  return text.substr(0, 72) + "...";
}

void report(const std::string &where, std::size_t lineNo,
            const std::string &what, const std::string &detail) {
  fprintf(stderr, "%s:%zu: %s: %s\n", where.c_str(), lineNo, what.c_str(),
          detail.c_str());
}

// Compares one line through both parsers and updates the mirrored PDBRUN
// version state. Returns true on parity.
bool compareLine(const std::string &line, const std::string &where,
                 std::size_t lineNo, int &inVersion, int &outVersion,
                 Counters &counters, bool quiet) {
  libpdb::PDB::PdbrunInputVersion(inVersion);
  libpdb::PDB::PdbrunOutputVersion(outVersion);

  const libpdb::PDB legacy(line.c_str());
  const std::string legacyText = legacy.chars();
  const int inAfterParse = libpdb::PDB::PdbrunInputVersion();
  const int outAfterRender = libpdb::PDB::PdbrunOutputVersion();

  zdock::pdb::ParseContext parseContext;
  parseContext.pdbrunVersion = inVersion;
  const std::optional<zdock::pdb::Record> parsed =
      zdock::pdb::Record::parse(line, parseContext);

  bool parity = true;
  if (!parsed.has_value()) {
    ++counters.rejections;
    parity = false;
    if (!quiet)
      report(where, lineNo, "new parser rejected",
             "legacy type " + std::to_string((int)legacy.type()) + " | " +
                 shorten(line));
  } else {
    // mirror the legacy render side effect: rendering a USER PDBRUN
    // record moves the legacy output version to the record's version
    int outForRender = outAfterRender;
    if (parsed->type() == zdock::pdb::RecordType::USER_PDBRUN)
      outForRender = parsed->as<zdock::pdb::UserPdbrunRecord>().version;
    zdock::pdb::WriteContext writeContext;
    writeContext.pdbrunVersion = outForRender;
    const std::string newText = parsed->chars(writeContext);

    if (static_cast<int>(parsed->type()) != static_cast<int>(legacy.type())) {
      ++counters.typeMismatches;
      parity = false;
      if (!quiet)
        report(where, lineNo, "type mismatch",
               "legacy " + std::to_string((int)legacy.type()) + " new " +
                   std::to_string((int)parsed->type()) + " | " +
                   shorten(line));
    }
    if (newText != legacyText) {
      ++counters.textMismatches;
      parity = false;
      if (!quiet)
        report(where, lineNo, "text mismatch",
               "legacy '" + shorten(legacyText) + "' new '" +
                   shorten(newText) + "'");
    }
  }

  // synchronize the mirrored version state with the legacy globals
  inVersion = inAfterParse;
  outVersion = outAfterRender;
  return parity;
}

} // namespace

int main(int argc, char **argv) {
  Options options;
  if (argc < 2) {
    fprintf(stderr,
            "usage: %s <file-or-directory> [--max-files N] [--max-lines N] "
            "[--all-extensions] [--quiet]\n",
            argv[0]);
    return 2;
  }
  options.root = argv[1];
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--max-files" && i + 1 < argc)
      options.maxFiles = static_cast<std::size_t>(std::atoi(argv[++i]));
    else if (arg == "--max-lines" && i + 1 < argc)
      options.maxLines = static_cast<std::size_t>(std::atoi(argv[++i]));
    else if (arg == "--all-extensions")
      options.allExtensions = true;
    else if (arg == "--quiet")
      options.quiet = true;
    else {
      fprintf(stderr, "unknown argument: %s\n", arg.c_str());
      return 2;
    }
  }

  std::vector<std::filesystem::path> files;
  try {
    const std::filesystem::path root(options.root);
    if (std::filesystem::is_regular_file(root)) {
      files.push_back(root);
    } else {
      for (const auto &entry :
           std::filesystem::recursive_directory_iterator(root))
        if (entry.is_regular_file() &&
            (options.allExtensions || pdbExtension(entry.path())))
          files.push_back(entry.path());
      std::sort(files.begin(), files.end());
    }
  } catch (const std::filesystem::filesystem_error &error) {
    fprintf(stderr, "filesystem error: %s\n", error.what());
    return 2;
  }
  if (options.maxFiles > 0 && files.size() > options.maxFiles)
    files.resize(options.maxFiles);
  if (files.empty()) {
    fprintf(stderr, "no input files under %s\n", options.root.c_str());
    return 2;
  }

  Counters counters;
  int inVersion = 6;
  int outVersion = 6;
  for (const auto &file : files) {
    std::ifstream in(file);
    if (!in) {
      fprintf(stderr, "cannot open %s\n", file.string().c_str());
      continue;
    }
    ++counters.files;
    const std::string where = file.string();
    std::string line;
    std::size_t lineNo = 0;
    while (std::getline(in, line)) {
      ++lineNo;
      ++counters.lines;
      compareLine(line, where, lineNo, inVersion, outVersion, counters,
                  options.quiet);
      if (options.maxLines > 0 && counters.lines >= options.maxLines)
        break;
    }
    if (options.maxLines > 0 && counters.lines >= options.maxLines)
      break;
  }

  printf("files %zu | lines %zu | type mismatches %zu | text mismatches %zu "
         "| rejections %zu\n",
         counters.files, counters.lines, counters.typeMismatches,
         counters.textMismatches, counters.rejections);
  return (counters.typeMismatches == 0 && counters.textMismatches == 0 &&
          counters.rejections == 0)
             ? 0
             : 1;
}
