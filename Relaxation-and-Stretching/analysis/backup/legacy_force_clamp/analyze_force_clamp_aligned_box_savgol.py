#!/usr/bin/env python3

from analyze_force_clamp_aligned_box import main


if __name__ == "__main__":
    main(default_trend_method="savgol", default_output_dir_name="analysis_savgol", default_sg_window=21)
