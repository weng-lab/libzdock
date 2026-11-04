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
//	$Id: pdb_sscanf.cc,v 1.8 94/09/06 15:03:23 gregc Exp $
//

#include "pdb++.h"
#include "FieldReader.hpp"
#include <ctype.h>
#include <stdarg.h>

//
//	pdb_sscanf performs similarly to sscanf, execept that fields are of
//	fixed length and a complete line is always consumed.  The field
//	width defaults to one.  If the line is shorter than expected then
//	the default is returned.
//
//		d	get an integer.  Default:  0.
//		f	get a floating point number (C double).  Default:  0.0.
//		(space) ignore characters within field
//		s	get a C string, leading and trailing spaces are
//			stripped; the field width is used as a limit on
//			the string length, the null character is appended
//			to the end of the string.  Default:  empty string.
//		c	get a character(s); no stripping of spaces, nor is
//			a null character appended.  Default:  space(s).
//
//	The field mechanics live in FieldReader; this wrapper only walks
//	the format string.
//

namespace libpdb {

int PDB::sscanf(const char *buffer, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  FieldReader reader(buffer);
  for (; *fmt != '\0'; fmt++) {
    if (*fmt != '%') {
      if (!reader.matchLiteral(*fmt)) {
        va_end(ap);
        return -1;
      }
      continue;
    }

    // calculate field_width
    int field_width = 0;
    for (++fmt; isdigit(static_cast<unsigned char>(*fmt)); fmt++)
      field_width = field_width * 10 + *fmt - '0';
    if (field_width == 0)
      field_width = 1; // default

    bool matched = true;
    switch (*fmt) {

    case 'd': // integer
      matched = reader.readInt(field_width, va_arg(ap, int *));
      break;

    case 'f': // floating point
      matched = reader.readReal(field_width, va_arg(ap, double *));
      break;

    case 's': // string
      reader.readString(field_width, va_arg(ap, char *));
      break;

    case 'c': // character(s)
      reader.readChars(field_width, va_arg(ap, char *));
      break;

    case ' ': // space (ignore)
      reader.skipChars(field_width);
      break;

    default:
      va_end(ap);
      return -1;
    }
    if (!matched) {
      va_end(ap);
      return -1;
    }
  }
  va_end(ap);
  return reader.nmatch();
}

} // namespace libpdb
