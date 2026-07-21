# BSD make reads this file. GNU make reads GNUmakefile first.

all test clean doc cpp-test python-test check-deps:
	@echo "This project requires GNU make; run 'gmake $@' on BSD." >&2
	@exit 1
