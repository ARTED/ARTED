#!/usr/bin/env python
import sys
import os
from subprocess import call
from os.path import abspath, dirname, isfile, join

if os.environ["USER"] == "jenkins":
    build_id = os.environ["BUILD_ID"]
    build_number = os.environ["BUILD_NUMBER"]
    workspace = os.environ["WORKSPACE"]
else:
    build_id, build_number = 0, 0
    workspace = join(dirname(abspath(__file__)), "../../")

for target in "sc", "ms":
    tempdir = join(workspace, "build_temp_%s" % target)
    os.mkdir(tempdir)
    ret = call("../configure.py --target=%s --arch=intel-knl && make" % target, shell=True, cwd=tempdir)
    out = join(workspace, "bin/ARTED_%s.mic" % target)
    if (ret != 0) or (not isfile(out)):
        print("Failure")
        sys.exit(-1)

print("Success id=%s number=%s" % (build_id, build_number))
sys.exit(0)
