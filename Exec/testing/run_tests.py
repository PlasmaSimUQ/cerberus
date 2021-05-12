#!/usr/bin/env python3

import os
from subprocess import Popen, PIPE

tests = os.listdir(".")

failed = []

for test in tests:

    if os.path.isfile(test):
        continue

    local_files = os.listdir(test)
    if "check.py" not in local_files:
        continue

    print("running: ",test, flush=True)

    p = Popen("cd %s; sh run"%test, shell=True)
    output, err = p.communicate()
    rc = p.returncode

    if rc > 0:
        print("failed:\n",err,"\n",output)
        failed.append({"test":test, "output":output, "error":err, "rc":rc})

    print("done", flush=True)



print("Completed with {} failed tests:".format(len(failed)))
for fail in failed:
    print("  ",fail["test"])

