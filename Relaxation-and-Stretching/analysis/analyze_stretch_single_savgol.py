#!/usr/bin/env python3

from analyze_stretch_single import main


if __name__ == "__main__":
    main(default_trend_method="savgol", default_output_dir_name="analysis_stretch_savgol", default_sg_window=21)
