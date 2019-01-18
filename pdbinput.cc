//
//	Copyright (c) 1992 The Regents of the University of California.
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
//	$Id: pdbinput.cc,v 1.1 94/02/22 17:09:27 gregc Exp $
//

#include "pdb++.h"

std::istream &
operator>>(std::istream &s, PDB &p)
{
	char	buf[4 * PDB::BufLen];

	s.getline(buf, 4 * PDB::BufLen);
	if (s)
		p = PDB(buf);
	return s;
}
