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
//	$Id: pdb_chars.cc,v 1.8 95/02/28 14:09:52 gregc Exp $
//
//	subroutine for writing PDB format files
//

#include "pdb++.h"
#include "PdbFormat.hpp"
#include <ctype.h>

namespace libpdb {

const char *PDB::chars(void) const {
  static char buf[BufLen];

  // convert C structure to pdb record

  int count =
      PdbFormat::formatRecord(*this, buf, rType, pdbrunOutputVersion);

  // Bounded system formatting reports the required size, so clamp before
  // scanning the buffer. The legacy formatter is addressed separately.
  if (count < 0)
    count = 0;
  else if (count >= BufLen)
    count = BufLen - 1;

  // find last non-blank in buf, and shorten it
  while (count > 1 && isspace(static_cast<unsigned char>(buf[count - 1])))
    count -= 1;
  buf[count] = '\0';
  return buf;
}

} // namespace libpdb
