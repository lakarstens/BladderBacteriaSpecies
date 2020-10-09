#!/usr/bin/env python3

import subprocess

checkit=subprocess.run(["python", "run_hw.py"], shell=True, capture_output=True, cwd="../resource_files")

print("\n\targs={}\n\treturn code={}\n\tstdout={}\n\tstderr={}".format(checkit.args, checkit.returncode, checkit.stdout, checkit.stderr))