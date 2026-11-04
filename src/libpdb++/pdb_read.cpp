//
//	Copyright (c) 1989, 1992 The Regents of the University of California.
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
//	$Id: pdb_read.cc,v 1.8 95/02/28 14:09:28 gregc Exp $
//
//	subroutine for reading PDB format files
//

#include "pdb++.h"
#include "PdbFormat.hpp"
#include <cstring>

namespace libpdb {

PDB::PDB(const char *buf) {
  // convert pdb record to C structure

  std::memset(reinterpret_cast<unsigned char *>(this), 0, sizeof *this);
  rType = getType(buf);
  rType = PdbFormat::parse(*this, buf, rType, pdbrunInputVersion);
}

} // namespace libpdb
