<!--
Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
SPDX-License-Identifier: BSD-2-Clause
-->

# ZDOCK Python package

Parser for ZDOCK and M-ZDOCK output. Reads output from all available versions in all variants. This module is part of the [Wenglab](https://zlab.umassmed.edu) [libzdock](https://github.com/weng-lab/libzdock.git) project.

## Development

Install the package and its development tools into uv's managed environment:

```console
uv sync --dev
```

Run the tests and quality checks from this directory:

```console
uv run python -m unittest discover -s tests -v
uv run isort --check-only src tests
uv run black --check src tests
uv run mypy src tests
uv run pylint src tests
```

Build source and wheel distributions:

```console
uv build
```

## Usage

```python
from zdock import ZDOCK

docking = ZDOCK("/path/to/zdock.out")

print(docking.isswitched)
print(docking.receptor)
print(docking.predictions)
```

The root GNU Make workflow builds the native command-line tools before running the Python integration tests:

```console
make python-test
```
