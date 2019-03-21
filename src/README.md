# Programs

A number of small utilities are built by default, to facilitate basic operations on ZDOCK, M-ZDOCK and PDB files. Each tool is described breifly below.


## centroids
Generates a PDB file with HETATM records indicating the center of mass for the top-_N_ predictions in a ZDOCK output file.

**Usage**
```
usage: centroids [options] <zdock output>

  -n <integer>    number of centroids to generate (top-n) (defaults to 1; the top prediction)
  -l <filename>   ligand PDB filename; defaults to receptor in ZDOCK output
  -c <char>       chain id to use for output (defaults to 'Z')
```

The output looks like this:

```
HETATM    1 N    HOH Z   1      36.995  43.928  73.278  1.00  0.00           N
HETATM    2 N    HOH Z   2      30.233  44.599  71.103  1.00  0.00           N
HETATM    3 N    HOH Z   3      63.936  22.304  57.441  1.00  0.00           N
HETATM    4 N    HOH Z   4      63.456  20.526  57.747  1.00  0.00           N
HETATM    5 N    HOH Z   5      45.213  21.873  26.833  1.00  0.00           N
```


## constraints
Performs constraint based filtering on ZDOCK and M-ZDOCK output, i.e. based on
a minimum or maximum distance between atoms in the predicted complexes. Distance
constraints are specified in a constraints file, of which the format is described below.

  - The format is line-based and each line represents one constraint.
  - Each line is whitespace separated and contains the following columns:
    + Atom serial number in the _first_ structure
    + Atom name (and altLoc)
    + Residue name
    + chain id (exactly _one_ character)
    + residue sequence number
    + Atom serial number in the _second_ structure
    + Atom name (and altLoc)
    + Residue name
    + chain id (exactly _one_ character)
    + residue sequence number
    + the last column is optional and contains "MIN" or "MAX" indicating the type
      of the constraint; minimum or maximum distance. If not specified the default
      constraint type is "MAX".

For ZDOCK constraint based filtering, the two sets of columns refer to atoms in the first and
the second structure respectively, whereas for M-ZDOCK where only only structure is operated on,
boths sets of columns refer to the same single structure.

Example constraints file:
```
13 OE2 GLU A 5    101 OD1 ASP B 12 7.3  MIN
13 OE2 GLU A 5    201 O   GLU B 12 7.5
13 OE2 GLU A 5    234 CA  PRO B 12 20.0 MAX
```

**Usage**
```
usage: constraints [options] <zdock output> <constraints file>

  -r <filename>   receptor PDB filename; defaults to receptor in (M-)ZDOCK output
  -l <filename>   ligand PDB filename; defaults to ligand in ZDOCK output
```

## createlig
Creates a complex or ligand for a ZDOCK prediction. Transformations are documented
[here](https://github.com/weng-lab/libpdb/blob/master/src/zdock/TransformLigand.hpp#L74).

**Usage**
```
usage: createlig [options] <zdock output>

  -n <integer>    index of prediction in ZDOCK file (defaults to 1; the top prediction)
  -c              create complex; by default only ligand is created
  -r <filename>   receptor PDB filename; defaults to receptor in ZDOCK output
  -l <filename>   ligand PDB filename; defaults to ligand in ZDOCK output
  -a              return all records (by default only ATOM and HETATM are returned)
```

## createmultimer
Creates a multimer complex (or a single component) from a M-ZDOCK output file. Relevant
transformations are documented
[here](https://github.com/weng-lab/libpdb/blob/master/src/zdock/TransformMultimer.hpp#L87).

**Usage**
```
usage: createmultimer [options] <zdock output>

  -n <integer>    index of prediction in M-ZDOCK file (defaults to 1; the top prediction)
  -r <filename>   receptor PDB filename; defaults to receptor in M-ZDOCK output
  -m <mer>        component of multimer to output (all if not specified)
  -a              return all records (by default only ATOM and HETATM are returned)
```

## pruning
Performs pruning on ZDOCK and M-ZDOCK output (using the greedy algorthm published here (TODO: add link)

**Usage**
```
usage: pruning [options] <zdock output>

  -c <double>     cutoff RMSD (defaults to 16.00)
  -C              return all prediction, but with score replaced by
                  cluster number.
  -l <filename>   structure PDB filename; defaults to ligand in ZDOCK
```

## zdsplit
Splits a ZDOCK file into equally sized chunks, each retaining the original header (akin to UNIX' _split_)

**Usage**
```
usage: zdsplit [options] <zdock output>

  -n <integer>    chunk size (defaults to input size)
  -p <string>     output filename prefix (defaults to "zdsplit.")
```

## zdunsplit
Reconstitutes a ZDOCK file from multiple _chunks_ (akin to UNIX' _cat_).

**Usage**
```
usage: zdunsplit <zdock output> [file] [...]

```
