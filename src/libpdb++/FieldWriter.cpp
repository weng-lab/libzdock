// Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
// SPDX-License-Identifier: BSD-2-Clause

#include "FieldWriter.hpp"

#include <string.h>

#define OVERFLOW_CHAR '*'

namespace libpdb {

FieldWriter::FieldWriter(char *out) noexcept : out_(out), cursor_(out) {}

void FieldWriter::writeChar(char c) { *cursor_++ = c; }

void FieldWriter::writePaddedChar(char c, int width, char fill,
                                  bool leftJustify) {
  if (leftJustify)
    *cursor_++ = c;
  while (--width > 0)
    *cursor_++ = fill;
  if (!leftJustify)
    *cursor_++ = c;
}

void FieldWriter::overflow(int width) {
  while (width-- > 0)
    *cursor_++ = OVERFLOW_CHAR;
}

void FieldWriter::writeInt(int value, int width, int radix, char fill,
                           char hex, bool leftJustify, char zero) {
  char *s;
  int n;
  bool negative;

  if (value < 0)
    negative = true, value = -value, width--;
  else
    negative = false;
  s = scratch_;
  if (value)
    do {
      n = value % radix;
      *s++ = n < 10 ? '0' + n : hex + n - 10;
      value /= radix;
    } while (value);
  else
    *s++ = zero;
  n = s - scratch_;
  if (width != -1 && n > width) {
    overflow(width + negative);
    return;
  }

  if (negative && fill == '0')
    *cursor_++ = '-';
  if (!leftJustify)
    while (width-- > n)
      *cursor_++ = fill;
  if (negative && fill == ' ')
    *cursor_++ = '-';
  while (--s >= scratch_)
    *cursor_++ = *s;
  if (leftJustify)
    while (width-- > n)
      *cursor_++ = fill;
}

void FieldWriter::writeUnsigned(unsigned value, int width, char fill,
                                bool leftJustify) {
  char *s;
  int n;

  s = scratch_;
  while (value) {
    *s++ = value % 10 + '0';
    value /= 10;
  }
  n = s - scratch_;
  if (n == 0)
    *s++ = '0', n = 1;
  if (width != -1 && n > width) {
    overflow(width);
    return;
  }

  if (!leftJustify)
    while (width-- > n)
      *cursor_++ = fill;
  while (--s >= scratch_)
    *cursor_++ = *s;
  if (leftJustify)
    while (width-- > n)
      *cursor_++ = fill;
}

void FieldWriter::writeString(const char *s, int width, int maxLen, char fill,
                              bool leftJustify) {
  int len;

  len = strlen(s);
  if (maxLen >= 0 && len > maxLen)
    len = maxLen;
  if (width != -1 && len > width) {
    overflow(width);
    return;
  }

  if (!leftJustify)
    while (width-- > len)
      *cursor_++ = fill;
  else
    width -= len;
  while (len--)
    *cursor_++ = *s++;
  if (leftJustify)
    while (width-- > 0)
      *cursor_++ = fill;
}

void FieldWriter::writeFloat(double value, int width, int places, char fill,
                             bool leftJustify) {
  int i, intval;
  char *place, *to, *from;
  bool negative;
  char *const p = cursor_;

  negative = value < 0.0 ? true : false;

  if (negative)
    value = -value;

  for (i = 0; i < places; i++)
    value *= 10.0;

  intval = (int)(value + 0.5);

  if (width == -1)
    width = places + 4; /* TODO: fix */
  else if (places + (places == 0 ? 1 : 2) > width) {
    overflow(width);
    return;
  }

  for (place = p + width - 1; place >= p + width - places; place--) {
    *place = '0' + intval % 10;
    intval /= 10;
  }

  if (places > 0)
    *place-- = '.';

  if (intval == 0)
    *place-- = '0';

  for (; place >= p; place--) {
    if (intval == 0)
      break;
    else {
      *place = '0' + intval % 10;
      intval /= 10;
    }
  }

  if (intval != 0) {
    overflow(width);
    return;
  }

  if (place < p && negative) {
    overflow(width);
    return;
  }

  if (leftJustify) {
    for (from = place + 1, to = (negative ? p + 1 : p); from < p + width;
         from++, to++)
      *to = *from;
    for (; to < p + width; to++)
      *to = fill;
    if (negative)
      *p = '-';
  } else {
    for (to = place; to >= p; to--)
      *to = fill;
    if (negative) {
      if (fill == ' ')
        *place = '-';
      else
        *p = '-';
    }
  }

  cursor_ = p + width;
}

void FieldWriter::writeExp(double value, int width, int places, char fill,
                           bool leftJustify) {
  int n;
  char *s;
  bool negative;
  double fraction;

  if (value < 0)
    negative = true, value = -value, width--;
  else
    negative = false;

  n = 0;
  while (value > 10)
    n++, value /= 10;
  if (value)
    while (value < 1)
      n--, value *= 10;

  s = scratch_;
  if (n < 0) {
    n = -n;
    *s++ = n % 10 + '0';
    *s++ = n / 10 + '0';
    *s++ = '-';
  } else {
    *s++ = n % 10 + '0';
    *s++ = n / 10 + '0';
    *s++ = '+';
  }
  *s = 'e';

  s = scratch_ + places + 4; /* 4 == strlen("e+00") */
  fraction = value - (int)value;
  for (n = 0; n < places; n++) {
    fraction *= 10.0;
    *--s = '0' + (int)fraction;
    fraction -= (int)fraction;
  }

  s = scratch_ + places + 4;
  if (places)
    *s++ = '.';
  n = (int)value;
  if (n)
    *s++ = n % 10 + '0';
  else
    *s++ = '0';
  n = s - scratch_;
  if (width != -1 && n > width) {
    overflow(width + negative);
    return;
  }

  if (negative && fill == '0')
    *cursor_++ = '-';
  if (!leftJustify)
    while (width-- > n)
      *cursor_++ = fill;
  if (negative && fill == ' ')
    *cursor_++ = '-';
  while (--s >= scratch_)
    *cursor_++ = *s;
  if (leftJustify)
    while (width-- > n)
      *cursor_++ = fill;
}

void FieldWriter::terminate() { *cursor_ = '\0'; }

int FieldWriter::count() const noexcept { return cursor_ - out_; }

} // namespace libpdb
