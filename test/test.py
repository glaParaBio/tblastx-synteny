#!/usr/bin/env python3

import shutil
import os
import sys
import subprocess
import unittest
import pandas
import glob

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

    # def tearDown(self):
    #    if os.path.exists("test_out"):
    #        shutil.rmtree("test_out")

    def testCanPrintHelp(self):
        cmd = f"{py} --help"
        p = shell(cmd)
        self.assertTrue(len(p.stdout) > 4)

    def testCaptureUserError(self):
        cmd = f"{py} --query foo --subject test_out/input/PbANKA_MIT_v3.fa -d test_out"
        p = shell(cmd, strict=False)
        self.assertTrue(p.returncode != 0)
        self.assertTrue("Missing input files" in p.stderr)

        cmd = f"{py} --subject test_out/input/PbANKA_MIT_v3.fa"
        p = shell(cmd, strict=False)
        self.assertTrue(p.returncode != 0)
        self.assertTrue("the following arguments are required" in p.stderr)

    def testTblastxSingleVsMultipleChunks(self):
        shell("cp test/data/PbANKA_MIT_v3.fa test/data/Pf3D7_MIT_v3.fa test_out/input/")
        cmd = f"{py} -q test_out/input/PbANKA_MIT_v3.fa -s test_out/input/Pf3D7_MIT_v3.fa -d test_out/single -S 20000"
        p = shell(cmd)
        out = glob.glob("test_out/single/tblastx/PbANKA_MIT_v3:*.out")
        self.assertEqual(len(out), 1)

        cmd = f"{py} -q test_out/input/PbANKA_MIT_v3.fa -s test_out/input/Pf3D7_MIT_v3.fa -d test_out/multiple -S 2000"
        p = shell(cmd)
        out = glob.glob("test_out/multiple/tblastx/PbANKA_MIT_v3:*.out")
        self.assertEqual(len(out), 3)

        # Test splitting the output gives nearly the same results (can't be the
        # same because of breakpoints)
        single = pandas.read_csv(
            "test_out/single/PbANKA_MIT_v3_vs_Pf3D7_MIT_v3.paf", sep="\t"
        )
        multiple = pandas.read_csv(
            "test_out/multiple/PbANKA_MIT_v3_vs_Pf3D7_MIT_v3.paf", sep="\t"
        )

        # Exclude evalue from comparison
        single = single.drop("blast_evalue", axis=1)
        multiple = multiple.drop("blast_evalue", axis=1)

        # 10 top and bottom rows are the same
        self.assertTrue(single.iloc[0:10,].equals(multiple.iloc[0:10,]))
        self.assertTrue(
            single.tail(10)
            .reset_index(drop=True)
            .equals(multiple.tail(10).reset_index(drop=True))
        )

    def testPerfectMatchSingleVsMultipleChunks(self):
        shell(
            "cp test/data/Pf3D7_MIT_v3.fa test/data/PF3D7_MIT02100.fa test_out/input/"
        )
        cmd = f"{py} -s test_out/input/PF3D7_MIT02100.fa -q test_out/input/Pf3D7_MIT_v3.fa -d test_out/single -S 20000 -b '-seg no -evalue 0.01' "
        p = shell(cmd)
        out = glob.glob("test_out/single/tblastx/Pf3D7_MIT_v3:*.out")
        self.assertEqual(len(out), 1)

        cmd = f"{py} -s test_out/input/PF3D7_MIT02100.fa -q test_out/input/Pf3D7_MIT_v3.fa -d test_out/multiple -S 1900 -b '-seg no -evalue 0.01'"
        p = shell(cmd)
        out = glob.glob("test_out/multiple/tblastx/Pf3D7_MIT_v3:*.out")
        self.assertEqual(len(out), 4)

        # These may be too sensitive:
        p = shell(
            "diff test_out/single/Pf3D7_MIT_v3_vs_PF3D7_MIT02100.paf test_out/multiple/Pf3D7_MIT_v3_vs_PF3D7_MIT02100.paf"
        )
        self.assertEqual(p.stdout, "")
        p = shell(
            "diff test/data/expected_Pf3D7_MIT_v3_vs_PF3D7_MIT02100.paf test_out/multiple/Pf3D7_MIT_v3_vs_PF3D7_MIT02100.paf"
        )
        self.assertEqual(p.stdout, "")

    def testRerunTrigger(self):
        shell("cp test/data/PbANKA_MIT_v3.fa test/data/Pf3D7_MIT_v3.fa test_out/input/")
        cmd = f"{py} -q test_out/input/PbANKA_MIT_v3.fa -s test_out/input/Pf3D7_MIT_v3.fa -d test_out -S 2000"
        shell(cmd)
        cmd = f"{py} -q test_out/input/PbANKA_MIT_v3.fa -s test_out/input/Pf3D7_MIT_v3.fa -d test_out -S 1000"
        p = shell(cmd)
        self.assertTrue("Params have changed since last execution" in p.stderr)
        cmd = f"{py} -q test_out/input/PbANKA_MIT_v3.fa -s test_out/input/Pf3D7_MIT_v3.fa -d test_out -S 1000"
        p = shell(cmd)
        self.assertTrue("Nothing to be done" in p.stderr)
        cmd = f"{py} -q test_out/input/PbANKA_MIT_v3.fa -s test_out/input/Pf3D7_MIT_v3.fa -d test_out -S 1000 --blast-args '-evalue 0.123'"
        p = shell(cmd)
        self.assertTrue("Params have changed since last execution" in p.stderr)

    def testBlastn(self):
        shell("cp test/data/PbANKA_MIT_v3.fa test/data/Pf3D7_MIT_v3.fa test_out/input/")
        cmd = f"{py} -q test_out/input/PbANKA_MIT_v3.fa -s test_out/input/Pf3D7_MIT_v3.fa -d test_out -S 1000 --task blastn"
        shell(cmd)
        paf = pandas.read_csv('test_out/PbANKA_MIT_v3_vs_Pf3D7_MIT_v3.paf', sep='\t')
        top = paf.loc[0]
        self.assertEqual(top['qstart'], 0)
        self.assertEqual(top['qend'], 1000)
        self.assertEqual(top['sstart'], 2)
        self.assertEqual(top['send'], 1009)

##qaccver      | qaccver_length | qstart | qend | strand | saccver      | saccver_length | sstart | send | nident | length | mapq | blast_pident | blast_nident | blast_evalue | blast_bitscore
#PbANKA_MIT_v3 |           5957 |      0 | 1000 | +      | Pf3D7_MIT_v3 |           5967 |      2 | 1009 | *      |   1014 |  255 | pi:f:88.659  | ni:i:899     | ev:f:0.0     | bs:f:1262
#PbANKA_MIT_v3 |           5957 |   1000 | 1998 | +      | Pf3D7_MIT_v3 |           5967 |   1009 | 2001 | *      |   1000 |  255 | pi:f:88.5    | ni:i:885     | ev:f:0.0     | bs:f:1277
#PbANKA_MIT_v3 |           5957 |   2000 | 3000 | +      | Pf3D7_MIT_v3 |           5967 |   2003 | 2999 | *      |   1000 |  255 | pi:f:86.3    | ni:i:863     | ev:f:0.0     | bs:f:1181
#PbANKA_MIT_v3 |           5957 |   3001 | 4000 | +      | Pf3D7_MIT_v3 |           5967 |   3000 | 3992 | *      |    999 |  255 | pi:f:83.684  | ni:i:836     | ev:f:0.0     | bs:f:1069
#PbANKA_MIT_v3 |           5957 |   4001 | 4999 | +      | Pf3D7_MIT_v3 |           5967 |   3993 | 4995 | *      |   1006 |  255 | pi:f:87.674  | ni:i:882     | ev:f:0.0     | bs:f:1249

if __name__ == "__main__":
    unittest.main()
