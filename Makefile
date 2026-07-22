# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

# BSD make reads this file. GNU make reads GNUmakefile first.

all test clean doc cpp-test python-check python-test check-deps:
	@echo "This project requires GNU make; run 'gmake $@' on BSD." >&2
	@exit 1
