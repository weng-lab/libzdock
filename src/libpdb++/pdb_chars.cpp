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
  // The legacy fixed-width formatter is unbounded, and its widest layout,
  // the 82-column ZDATOM record, needs two more cells than BufLen. Every
  // record is clamped to the 80-character record length except ZDATOM,
  // which keeps all 82 columns of its layout.
  static char buf[BufLen + 8];

  // convert C structure to pdb record

  int count =
      PdbFormat::formatRecord(*this, buf, rType, pdbrunOutputVersion);
  const bool zdatom = (rType == ATOM || rType == HETATM ||
                       rType == SIGATM) && atom.iszdatom;

  // Bounded system formatting reports the required size, so clamp before
  // scanning the buffer. Non-ZDATOM counts (including the would-be lengths
  // of the v5 USER libc paths) clamp to the 80-character record length as
  // they have since the bounds were introduced; a ZDATOM render writes at
  // most its 82-column layout, well inside the buffer.
  if (count < 0)
    count = 0;
  else if (!zdatom) {
    if (count >= BufLen)
      count = BufLen - 1;
  } else if (count >= static_cast<int>(sizeof(buf)))
    count = static_cast<int>(sizeof(buf)) - 1;

  // find last non-blank in buf, and shorten it
  while (count > 1 && isspace(static_cast<unsigned char>(buf[count - 1])))
    count -= 1;
  buf[count] = '\0';
  return buf;
}

} // namespace libpdb
