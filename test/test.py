#!/usr/bin/env python3

import shutil
import os
import sys
import subprocess
import unittest

py = "./tblastx-synteny.py"

class shell:
    def __init__(self, cmd, strict=True, timeout=None):
        print(cmd)
        cmd = f"set -e; set -u; set -o pipefail\n{cmd}"
        p = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            executable="/bin/bash",
        )
        try:
            stdout, stderr = p.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            p.kill()
            sys.stderr.write(f"Error: Timeout after {timeout} seconds\n")
            stdout, stderr = p.communicate()
        self.returncode = p.returncode
        self.stdout = stdout.decode()
        self.stderr = stderr.decode()
        self.cmd = cmd
        if strict and self.returncode != 0:
            raise subprocess.SubprocessError(
                f"\nSTDOUT:\n{self.stdout}\nSTDERR:\n{self.stderr}\nEXIT CODE: {self.returncode}"
            )


class Test(unittest.TestCase):
    def setUp(self):
        sys.stderr.write("\n" + self.id().split(".")[-1] + "\n")  # Print test name
        if os.path.exists("test_out"):
            shutil.rmtree("test_out")
        os.makedirs("test_out/input")

    def tearDown(self):
        if os.path.exists("test_out"):
            shutil.rmtree("test_out")

    def testCanPrintHelp(self):
        cmd = f"{py} --help"
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue(len(p.stdout) > 4)

    def canCaptureUserError(self):
        cmd = f"{py} --query foo --subject test_out/input/PbANKA_MIT_v3.fa"
        print(cmd)
        shell(cmd)


    def testTblastxSingleChunk(self):
        shell('cp test/data/PbANKA_MIT_v3.fa test/data/Pf3D7_MIT_v3.fa test_out/input/')
        cmd = f"{py} -q test_out/input/PbANKA_MIT_v3.fa -s test_out/input/Pf3D7_MIT_v3.fa -d test_out/"
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)



if __name__ == "__main__":
    unittest.main()
