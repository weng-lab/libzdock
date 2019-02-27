# PDB and ZDOCK

SYNOPSIS
--------

```C++
#include "PDB.hpp"
#include "ZDOCK.hpp"
#include "Utils.hpp"

....

// read pdb file (CA only)
PDB pdb("filename.pdb", [](const auto &r) { return Utils::trim_copy(r.atom.name) == "CA"; });

// read zdock file
ZDOCK z("zdock.out");

// transform class for ZDOCK (i.e. not M-ZDOCK)
TransformLigand txl(z);

// grab a prediction and transform the PDB atom coordinatea
const Prediction pred = z.predictions()[0];
pdb.setMatrix(txl.txLigand(pdb.matrix(), pred));

// print updated PDB contents.
for (const auto& x : pdb.records()) {
 std::cout << *x << '\n';
}

....
```


# libpdb++
libpdb++ from http://www.cgl.ucsf.edu/Overview/software.html#pdbio

PDB++
=====

[NAME](#NAME)\
 [SYNOPSIS](#SYNOPSIS)\
 [DESCRIPTION](#DESCRIPTION)\
 [MEMBER CONSTANTS](#MEMBER%20CONSTANTS)\
 [MEMBER TYPES](#MEMBER%20TYPES)\
 [MEMBER FUNCTIONS](#MEMBER%20FUNCTIONS)\
 [I/O FUNCTIONS](#I/O%20FUNCTIONS)\
 [SEE ALSO](#SEE%20ALSO)\
 [NOTES](#NOTES)\
 [BUGS](#BUGS)\
 [COPYRIGHT](#COPYRIGHT)


* * * * *

NAME
----

pdb++ − A C++ class for manipulating Brookhaven Protein DataBank records

SYNOPSIS
--------

```C++
#include <pdb++.h>

....

PDB record;

while (cin >> record) {
  switch (record.type()) {
    case PDB::ATOM:
      cout << record.atom.xyz[0] << ' ' << record.atom.xyz[0]
            << ' ' << record.atom.xyz[0] << endl;
      ....
      break;
  }
}
....
```

DESCRIPTION
-----------

The routines listed above are available in the pdb++ library,
*-L/usr/local/lib/midas -lpdb++*.

The **PDB** class provides methods for parsing Brookhaven Protein
DataBank (PDB) records (lines from a PDB file) into structures and
expanding those structures back into PDB records. Rather than provide
access functions for each possible field, the structure containing the
parsed record is publicly available. The field names are listed in the
header file.

The **PDB** class has two enhancements to the Brookhaven Protein
DataBank specification: four character residue names, and the PDBRUN set
of scene annotation records. Four character residue names work because
everywhere in the specification a three character residue name appears,
there is a blank afterwards.

MEMBER CONSTANTS
----------------

**BufLen**

The maximum length of a generated PDB record string (including the null
byte).

**PDBRUNVersion**

The default version of the PDBRUN scene annotation records.

There are also constants for each known PDB record type (*e.g.*,
**ATOM**, **HETATM**, *etc.*), and the constant **UNKNOWN** for an
unknown PDB record.

The following constants are for each PDBRUN scene annotation record
type: **USER\_PDBRUN**, **USER\_EYEPOS**, **USER\_ATPOS**,
**USER\_WINDOW**, **USER\_FOCUS**, **USER\_VIEWPORT**,
**USER\_BGCOLOR**, **USER\_ANGLE**, **USER\_DISTANCE**, **USER\_FILE**,
**USER\_MARKNAME**, **USER\_MARK**, **USER\_CNAME**, **USER\_COLOR**,
**USER\_RADIUS**, **USER\_OBJECT**, **USER\_ENDOBJ**, **USER\_CHAIN**,
**USER\_GFX\_BEGIN**, **USER\_GFX\_END**, **USER\_GFX\_COLOR**,
**USER\_GFX\_NORMAL**, **USER\_GFX\_VERTEX**, **USER\_GFX\_FONT**,
**USER\_GFX\_TEXTPOS**, and **USER\_GFX\_LABEL**.

The following constants are for the various graphics primitives
supported in scenes: **GFX\_UNKNOWN**, **GFX\_POINTS**,
**GFX\_MARKERS**, **GFX\_LINES**, **GFX\_LINE\_STRIP**,
**GFX\_LINE\_LOOP**, **GFX\_TRIANGLES**, **GFX\_TRIANGLE\_STRIP**,
**GFX\_TRIANGLE\_FAN**, **GFX\_QUADS**, **GFX\_QUAD\_STRIP**, and
**GFX\_POLYGON**.

MEMBER TYPES
------------

typedef char **Date**[10]

A text field containing a date, typically *day*−*month*−*year*, where
*day* is numeric, *month* is a three-letter abbreviation, and *year* is
the last two digits of the year.

typedef char **AName**[5]

A PDB atom name, *e.g.*, NO2\*.

typedef char **RName**[5]

Residue name, *e.g.*, ALA.

typedef char **PName**[5]

PDB name, *e.g.*, 9lyz.

typedef char **Id**[4]

Generic short id field.

typedef double **Real**

Size of floating point numbers read and written.

struct **Residue**

A **Residue** consists of a residue name (**name**), a chain identifier
(**chainId**), a sequence number (**seqNum**), and an insertion code
(**insertCode**).

MEMBER FUNCTIONS
----------------

**pdb**() **\
 pdb**(RecordType t) **\
 pdb**(const char \*buf)

Constructors. The first two above construct a zeroed instance of the
given record type (default **UNKNOWN**). The last constructor above
fills in all of the fields of the instance from the given PDB record
text.

RecordType **type**() const

Return the type of PDB instance.

void **type**(RecordType t)

Change the PDB record type of the instance and reinitialize all the
fields to default values (zero in all cases except for an **ATOM**’s
occupancy which defaults to 1.0).

const char \***chars**() const;

Return a string containing the PDB record in textual form.

static int **PdbrunInputVersion**()

Return the current PDBRUN scene annotation version used to parse text
records.

static int **PdbrunOutputVersion**()

Return the current PDBRUN scene annotation version used to create text
records.

static void **PdbrunInputVersion**(int v)

Set the current PDBRUN scene annotation version used to parse text
records.

static void **PdbrunOutputVersion**(int v)

Set the current PDBRUN scene annotation version used to create text
records.

static recordType **getType**(const char \*buf)

Return the PDB record type for the given line of text.

static GfxType **getGfxType**(const char \*buf)\
 static const char \***gfxChars**(GfxType gt)\
 static int **sscanf**(const char \*, const char \*, ...)

A version of **sscanf**(3) whose format’s behave like FORTRAN formats
where field widths are sacrosanct. If the input line is short, then the
rest of the fields are initialized to default values. Any literal
characters in the format must match the input. The format characters
are: space, ignore input character; **c**, character (array), default to
a space; **d**, integer, default zero; **f**, double, default zero;
**s**, get a C string, trailing spaces are stripped and it is null
terminated, default an empty string. **sscanf** returns the number of
input fields converted (may be less than expected if the input line is
short) or −1 if an error is found.

static int **sprintf**(char \*, const char \*, ...)

A version of **sprintf**(3) whose format’s behave like FORTRAN formats
where field widths are sacrosanct. Literal characters are copied as is.
If the text or number to be printed is larger than the given field
width, then the field is filled in with asterisks. The format characters
are: **d**, integer; **D**, integer where zero is written as spaces;
**s**, right-justified string (a negative field width left-justifies);
**c**, character (array), zero characters are converted to spaces;
**f**, floating point, normal **printf** precisions are used.

I/O FUNCTIONS
-------------

ostream &**operator\<\<**(ostream &s, const PDB &p)

Output the current PDB record on the given output stream.

istream &**operator\>\>**(istream &s, PDB &p)

Read a line from the given input stream and convert to a PDB record.

SEE ALSO
--------

‘‘Protein Data Bank Atomic Coordinate and Bibliographic Entry Format
Description,’’ Febuary 1992, Brookhaven National Laboratory, the January
1993 Protein Data Bank Quarterly Newsletter, ‘‘Annotating PDB Files with
Scene Information,’’ Gregory S. Couch, *et. al.*, (submitted for
publication).

NOTES
-----

The subtype field of USERxx structure tells what the *xx* part was. The
rest of the line, up to the card sequence portion, is the text field.

Due to the way Brookhaven encodes their files, atom names often have
leading blanks and sometimes have embedded blanks. Residue names
occasionally have leading blanks too. To be entirely consistent with the
PDB format, the programmer should put those blanks in before using the
**chars** member function.

BUGS
----

Routines are needed to convert to and from PDB typesetting conventions
in **COMPND**, **SOURCE**, **AUTHOR**, and **JRNL** records. Also, some
new record types have appeared since 1994, that currently show up as
**UNKNOWN**, such as **TITLE**.

COPYRIGHT
---------

Copyright © 1994 The Regents of the University of California. All rights
reserved.

Redistribution and use in source and binary forms are permitted provided
that the above copyright notice and this paragraph are duplicated in all
such forms and that any documentation, advertising materials, and other
materials related to such distribution and use acknowledge that the
software was developed by the University of California, San Francisco.
The name of the University may not be used to endorse or promote
products derived from this software without specific prior written
permission. THIS SOFTWARE IS PROVIDED ‘‘AS IS’’ AND WITHOUT ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.

* * * * *
