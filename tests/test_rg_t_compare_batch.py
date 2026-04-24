from __future__ import annotations

import importlib.util
import sys
import unittest
from array import array
from pathlib import Path
from tempfile import TemporaryDirectory


REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "Relaxation-and-Stretching" / "analysis"
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))


def load_module(module_name: str, relative_path: str):
    module_path = REPO_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


rg_t = load_module(
    "rg_T_analyze_force_clamp_aligned_box",
    "Relaxation-and-Stretching/analysis/rg_T_analyze_force_clamp_aligned_box.py",
)
compare_batch = load_module(
    "compare_batch_rg_T",
    "Relaxation-and-Stretching/analysis/compare_batch_rg_T.py",
)


class CompareBatchTests(unittest.TestCase):
    def test_stretch_dataset_schema_check_rejects_old_cache_layout(self) -> None:
        with TemporaryDirectory() as temp_dir_str:
            dataset_path = Path(temp_dir_str) / "stretch_dataset_example.dat"
            dataset_path.write_text(
                "# source_tag timestep time L dL strain contour_fraction Fext Fchain Ftotal Rg Rg2 count total_mass coord_source dump_file\n",
                encoding="utf-8",
            )

            self.assertFalse(rg_t.stretch_dataset_schema_matches(dataset_path))

    def test_compute_mass_weighted_rg_components_reports_total_and_x_projection(self) -> None:
        metrics = rg_t.compute_mass_weighted_rg_components(
            positions=[
                (0.0, 0.0, 0.0),
                (2.0, 0.0, 0.0),
                (4.0, 0.0, 0.0),
            ],
            masses=[1.0, 1.0, 1.0],
        )

        self.assertEqual(metrics["Lx"], 4.0)
        self.assertAlmostEqual(metrics["Rg_x"], (8.0 / 3.0) ** 0.5, places=6)
        self.assertAlmostEqual(metrics["Rg"], metrics["Rg_x"], places=6)

    def test_build_force_summary_rows_uses_last_window_mean_and_std(self) -> None:
        payloads = [
            {
                "source_tag": "restart860000_Tstar_0.50_F0.70_tail_v0_no_thermal_xuyu_zu",
                "table": {
                    "time": array("d", [0.0, 1.0, 2.0, 3.0]),
                    "L": array("d", [10.0, 12.0, 14.0, 16.0]),
                    "contour_fraction": array("d", [0.10, 0.20, 0.30, 0.40]),
                },
                "rg_kinematics": {
                    "time": array("d", [0.0, 1.0, 2.0, 3.0]),
                    "rg": array("d", [3.0, 4.0, 5.0, 6.0]),
                    "rg_x": array("d", [1.0, 2.0, 3.0, 4.0]),
                },
            },
            {
                "source_tag": "restart860000_Tstar_0.50_F1.20_tail_v0_no_thermal_xuyu_zu",
                "table": {
                    "time": array("d", [0.0, 1.0, 2.0, 3.0]),
                    "L": array("d", [20.0, 22.0, 24.0, 26.0]),
                    "contour_fraction": array("d", [0.50, 0.60, 0.70, 0.80]),
                },
                "rg_kinematics": {
                    "time": array("d", [0.0, 1.0, 2.0, 3.0]),
                    "rg": array("d", [7.0, 8.0, 9.0, 10.0]),
                    "rg_x": array("d", [4.0, 5.0, 6.0, 7.0]),
                },
            },
        ]

        summary_rows = compare_batch.build_force_summary_rows(payloads, window_points=3)

        self.assertEqual([row["force"] for row in summary_rows], [0.7, 1.2])
        first = summary_rows[0]
        self.assertEqual(first["mean_L"], 14.0)
        self.assertAlmostEqual(first["std_L"], (8.0 / 3.0) ** 0.5, places=6)
        self.assertEqual(first["mean_contour_fraction"], 0.30)
        self.assertEqual(first["mean_Rg"], 5.0)
        self.assertEqual(first["mean_Rg_x"], 3.0)


if __name__ == "__main__":
    unittest.main()
