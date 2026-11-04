//
//	Copyright (c) 1989,1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
//	$Id: pdb_sprntf.cc,v 1.4 94/09/06 15:03:07 gregc Exp $
//

#include "pdb++.h"
#include "FieldWriter.hpp"
#include <ctype.h>
#include <stdarg.h>

namespace libpdb {

int PDB::sprintf(char *outbuf, const char *fmt, ...) {
  va_list argv;
  const char *f;
  int field1, field2;
  char fill_char;
  double fnum;
  bool left_justify;

  va_start(argv, fmt);
  f = fmt;
  FieldWriter writer(outbuf);
  while (*f) {
    if (*f == '%') {
      f++;
      if (*f == '-')
        left_justify = true, f++;
      else
        left_justify = false;

      if (*f == '0')
        fill_char = '0', f++;
      else
        fill_char = ' ';

      if (isdigit(static_cast<unsigned char>(*f))) {
        field1 = *f++ - '0';
        while (isdigit(static_cast<unsigned char>(*f)))
          field1 = field1 * 10 + *f++ - '0';
      } else
        field1 = -1;

      if (*f == '.') {
        f++;
        field2 = 0;
        while (isdigit(static_cast<unsigned char>(*f)))
          field2 = field2 * 10 + *f++ - '0';
      } else
        field2 = -1;

      if (*f == 'l' || *f == 'h')
        f++;

      while (isspace(static_cast<unsigned char>(*f)))
        f++;
      switch (*f) {
      case 'c': {
        char c = (char)va_arg(argv, int);
        if (c == '\0')
          c = ' ';
        writer.writePaddedChar(c, field1, fill_char, left_justify);
        break;
      }
      case 'd':
      case 'D':
        writer.writeInt(va_arg(argv, int), field1, 10, fill_char, 'a',
                        left_justify, (*f == 'D') ? ' ' : '0');
        break;
      case 'e':
        fnum = va_arg(argv, double);
        if (field2 < 0)
          field2 = 6;
        writer.writeExp(fnum, field1, field2, fill_char, left_justify);
        break;
      case 'f':
        fnum = va_arg(argv, double);
        if (field2 < 0)
          field2 = 6;
        writer.writeFloat(fnum, field1, field2, fill_char, left_justify);
        break;
      case 'o':
        writer.writeInt(va_arg(argv, int), field1, 8, fill_char, 'a',
                        left_justify, '0');
        break;
      case 's':
        writer.writeString(va_arg(argv, char *), field1, field2, fill_char,
                           left_justify);
        break;
      case 'u':
        writer.writeUnsigned(va_arg(argv, unsigned), field1, fill_char,
                             left_justify);
        break;
      case 'x':
        writer.writeInt(va_arg(argv, int), field1, 16, fill_char, 'a',
                        left_justify, '0');
        break;
      case 'X':
        writer.writeInt(va_arg(argv, int), field1, 16, fill_char, 'A',
                        left_justify, '0');
        break;
      default:
        writer.writePaddedChar(*f, field1, fill_char, left_justify);
        break;
      }
      f++;
    } else if (*f == '\\') { /* Special character */
      switch (*++f) {
      case 'n':
        writer.writeChar('\n');
        break;
      case 'r':
        writer.writeChar('\r');
        break;
      case 'b':
        writer.writeChar('\b');
        break;
      case 't':
        writer.writeChar('\t');
        break;
      case 'f':
        writer.writeChar('\f');
        break;
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7': {
        int inum = *f++ - '0';
        if (*f >= '0' && *f <= '7') {
          inum = inum * 8 + *f++ - '0';
          if (*f >= '0' && *f <= '7')
            inum = inum * 8 + *f++ - '0';
        }
        f--;
        writer.writeChar((char)inum);
        break;
      }
      default:
        writer.writeChar(*f);
      }
      f++;
    } else /* Normal character */
      writer.writeChar(*f++);
  }
  writer.terminate();
  va_end(argv);
  return writer.count();
}

} // namespace libpdb
