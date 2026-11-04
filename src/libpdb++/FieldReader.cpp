// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "FieldReader.hpp"

#include <ctype.h>
#include <stdlib.h>

namespace libpdb {
namespace {
constexpr int kMaxFieldSize = 80; // legacy scratch size for numeric fields
}

FieldReader::FieldReader(const char *line) noexcept : cursor_(line), nmatch_(0) {}

bool FieldReader::atEnd() const noexcept {
  return *cursor_ == '\0' || *cursor_ == '\n';
}

int FieldReader::nmatch() const noexcept { return nmatch_; }

bool FieldReader::matchLiteral(char expected) noexcept {
  if (*cursor_ == expected) {
    ++cursor_;
    return true;
  }
  return atEnd(); // a literal past end-of-line is tolerated
}

bool FieldReader::beginField() {
  if (atEnd()) {
    return false;
  }
  ++nmatch_;
  return true;
}

bool FieldReader::readInt(int width, int *value) {
  if (!beginField()) {
    *value = 0;
    return true;
  }
  char tmp[kMaxFieldSize];
  char *s = tmp;
  for (int i = 0; i < width && !atEnd() && s - tmp < kMaxFieldSize - 1; ++i) {
    *s++ = *cursor_++;
  }
  *s = '\0';
  while (s > tmp && isspace(static_cast<unsigned char>(*(s - 1)))) {
    *--s = '\0';
  }
  char *t = nullptr;
  *value = static_cast<int>(strtol(tmp, &t, 10));
  return t == s;
}

bool FieldReader::readReal(int width, double *value) {
  if (!beginField()) {
    *value = 0.0;
    return true;
  }
  char tmp[kMaxFieldSize];
  char *s = tmp;
  for (int i = 0; i < width && !atEnd() && s - tmp < kMaxFieldSize - 1; ++i) {
    *s++ = *cursor_++;
  }
  *s = '\0';
  while (s > tmp && isspace(static_cast<unsigned char>(*(s - 1)))) {
    *--s = '\0';
  }
  char *t = nullptr;
  *value = strtod(tmp, &t);
  return t == s;
}

void FieldReader::readString(int width, char *value) {
  if (!beginField()) {
    *value = '\0';
    return;
  }
  char *s = value;
  for (int i = 0; i < width && !atEnd(); ++i) {
    *s++ = *cursor_++;
  }
  *s = '\0';
  while (s > value && isspace(static_cast<unsigned char>(*--s))) {
    *s = '\0';
  }
}

void FieldReader::readChars(int width, char *value) {
  for (int i = 0; i < width; ++i) {
    value[i] = ' ';
  }
  if (!beginField()) {
    return;
  }
  for (int i = 0; i < width && !atEnd(); ++i) {
    value[i] = *cursor_++;
  }
}

void FieldReader::skipChars(int width) {
  if (!beginField()) {
    return;
  }
  for (int i = 0; i < width && !atEnd(); ++i) {
    ++cursor_;
  }
}

} // namespace libpdb
