ZDOCK
=====

Parser for ZDOCK and M-ZDOCK output. Reads output from all available versions in all variants. This module is part of the [Wenglab](https://zlab.umassmed.edu) [libzdock](https://github.com/weng-lab/libzdock.git) project.

To install this Python package, run the following while in this directory:


system wide:

```
pip install .
```

or, if you'd like to install this for the current user only:

```
pip install --user .
```

This code could be used as follows:

```python
import ZDOCK

# load a ZDOCK output file
myfile = ZDOCK("/path/to/zdock.out")

# figure out whether or not receptor and ligand were swapped
print(myfile.isswitched)

# print receptor
print(myfile.receptor)

# get all predictions
print(myfile)
```
