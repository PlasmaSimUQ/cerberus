#!/usr/bin/env python3
"""Path-stable shim: the tool now lives in the eos_tools package.

Kept at this path so the test-suite `run` scripts' invocations
(`python3 ../../python_analysis/eos_table_prep.py synthetic|fpeos|qa ...`)
keep working unchanged. See eos_tools/cli.py for the commands and
eos_tools/formats/eostab.py for the frozen .eostab v1 spec.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from eos_tools.cli import main  # noqa: E402

if __name__ == "__main__":
    main()
