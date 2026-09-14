"""End-to-end tests against the compiled C executable; Python standard library only."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BINARY = Path(sys.argv[1]).resolve()
sys.argv = [sys.argv[0]]


class SimulatorCLI(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def run_sim(self, *args, success=True):
        result = subprocess.run(
            [str(BINARY), *map(str, args)], cwd=self.root,
            text=True, capture_output=True, timeout=30,
        )
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def rows(self, filename="warp_telemetry.csv"):
        with (self.root / filename).open(newline="") as stream:
            return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]

    def summary(self):
        with (self.root / "warp_summary.json").open() as stream:
            return json.load(stream, parse_constant=lambda text: self.fail(f"non-finite JSON: {text}"))

    def test_default_run(self):
        self.run_sim()
        rows = self.rows()
        self.assertEqual(len(rows), 1601)
        self.assertEqual(rows[0]["time_s"], 0)
        self.assertEqual(rows[0]["x_m"], 0)
        self.assertEqual(rows[0]["throttle"], 0)
        self.assertEqual(rows[-1]["time_s"], 1.6)
        self.assertEqual(rows[-1]["velocity_m_s"], 0)
        self.assertAlmostEqual(rows[-1]["x_m"], 749481145, delta=1e-5)
        self.assertTrue(all(a["time_s"] < b["time_s"] for a, b in zip(rows, rows[1:])))
        self.assertTrue(all(math.isfinite(v) for row in rows for v in row.values()))
        summary = self.summary()
        self.assertEqual(summary["samples"], len(rows))
        self.assertEqual(summary["final_position_m"][0], rows[-1]["x_m"])
        self.assertTrue(summary["quadrature"]["converged"])
        self.assertLess(summary["full_throttle_energy_j"], 0)
        self.assertEqual(summary["config"]["direction"], [1, 0, 0])

    def test_timestep_independent_trajectory(self):
        endpoints = []
        for index, dt in enumerate(("0.001", "0.07", "0.3", "4")):
            self.run_sim("--dt", dt, "--output", f"run{index}.csv",
                         "--summary", f"run{index}.json")
            rows = self.rows(f"run{index}.csv")
            self.assertEqual(rows[0]["time_s"], 0)
            self.assertEqual(rows[-1]["time_s"], 1.6)
            self.assertTrue(all(a["time_s"] < b["time_s"] for a, b in zip(rows, rows[1:])))
            endpoints.append(rows[-1]["x_m"])
        self.assertEqual(len(set(endpoints)), 1)

    def test_direction_and_reverse(self):
        self.run_sim("--warp", "-2", "--direction", "0", "3", "4", "--dt", "0.2")
        final = self.rows()[-1]
        self.assertEqual(final["x_m"], 0)
        self.assertAlmostEqual(final["y_m"], -0.6 * 749481145, delta=1e-5)
        self.assertAlmostEqual(final["z_m"], -0.8 * 749481145, delta=1e-5)
        actual = self.summary()["config"]["direction"]
        for a, b in zip(actual, [0, 0.6, 0.8]):
            self.assertAlmostEqual(a, b, places=14)

    def test_zero_duration(self):
        self.run_sim("--duration", "0")
        self.assertEqual(len(self.rows()), 1)
        self.assertEqual(self.rows()[0]["time_s"], 0)

    def test_nondivisible_endpoint(self):
        self.run_sim("--duration", "0.3", "--dt", "0.07")
        times = [r["time_s"] for r in self.rows()]
        self.assertEqual(len(times), 6)
        self.assertEqual(times[-1], 0.3)
        self.assertTrue(all(a < b for a, b in zip(times, times[1:])))

    def test_zero_warp(self):
        self.run_sim("--warp", "0", "--dt", "0.5")
        self.assertTrue(all(r["x_m"] == 0 and r["integrated_energy_j"] == 0 for r in self.rows()))
        self.assertEqual(self.summary()["full_throttle_energy_j"], 0)

    def test_field_slice(self):
        self.run_sim("--dt", "0.2", "--slice", "field.csv", "--slice-time", "0.75",
                     "--slice-points", "5", "--slice-extent", "100")
        rows = self.rows("field.csv")
        self.assertEqual(len(rows), 25)
        center = next(r for r in rows if r["longitudinal_m"] == r["transverse_m"] == 0)
        self.assertEqual(center["shape"], 1)
        self.assertEqual(center["energy_density_j_m3"], 0)
        axis = [r for r in rows if r["transverse_m"] == 0]
        self.assertTrue(all(r["energy_density_j_m3"] == 0 for r in axis))
        wall = next(r for r in rows if r["longitudinal_m"] == 0 and r["transverse_m"] == 100)
        self.assertLess(wall["energy_density_j_m3"], 0)
        front = next(r for r in rows if r["longitudinal_m"] == 100 and r["transverse_m"] == 0)
        back = next(r for r in rows if r["longitudinal_m"] == -100 and r["transverse_m"] == 0)
        self.assertLess(front["expansion_per_s"], 0)
        self.assertGreater(back["expansion_per_s"], 0)
        self.assertTrue(all(math.isfinite(v) for row in rows for v in row.values()))

    def test_invalid_inputs_create_no_outputs(self):
        cases = [
            ["--dt", "0"], ["--dt", "-1"], ["--dt", "nan"],
            ["--dt", "1e-30"], ["--radius", "-1"], ["--radius", "inf"],
            ["--sigma", "0"], ["--sigma", "1e300"], ["--warp", "1001"],
            ["--warp", "2junk"], ["--direction", "0", "0", "0"],
            ["--direction", "1", "0"], ["--direction", "nan", "1", "0"],
            ["--ramp-up", "0"], ["--ramp-down", "-1"], ["--cruise", "-1"],
            ["--duration", "-1"], ["--energy-intervals", "33"],
            ["--energy-intervals", "-32"], ["--energy-intervals", "1024.0"],
            ["--energy-intervals", "999999999999999999999999"],
            ["--slice-time", "0.5"], ["--slice", "field.csv", "--slice-time", "2"],
            ["--slice", "field.csv", "--slice-points", "2"],
            ["--slice", "field.csv", "--slice-extent", "0"],
            ["--output", ""], ["--output", "-"],
            ["--output", "same.csv", "--summary", "same.csv"],
            ["--unknown", "1"], ["--radius"],
        ]
        for args in cases:
            with self.subTest(args=args):
                self.run_sim(*args, success=False)
                self.assertEqual(list(self.root.iterdir()), [])

    def test_existing_file_preserved(self):
        output = self.root / "warp_telemetry.csv"
        output.write_text("preserve me")
        self.run_sim(success=False)
        self.assertEqual(output.read_text(), "preserve me")
        self.assertFalse((self.root / "warp_summary.json").exists())

    def test_later_open_failure_cleans_new_outputs(self):
        summary = self.root / "warp_summary.json"
        summary.write_text("existing summary")
        self.run_sim(success=False)
        self.assertEqual(summary.read_text(), "existing summary")
        self.assertFalse((self.root / "warp_telemetry.csv").exists())
        self.run_sim("--summary", "missing/summary.json", success=False)
        self.assertFalse((self.root / "warp_telemetry.csv").exists())

    def test_path_alias_cannot_truncate_output(self):
        self.run_sim("--output", "same.csv", "--summary", "./same.csv", success=False)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_help_and_version(self):
        self.assertIn("--slice", self.run_sim("--help").stdout)
        self.assertEqual(self.run_sim("--version").stdout.strip(), "0.2.0")
        self.assertEqual(list(self.root.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
