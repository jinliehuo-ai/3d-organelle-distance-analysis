# Reproducibility checklist

- Record Fiji, Bio-Formats, 3D ImageJ Suite, Python, NumPy, SciPy, Matplotlib, and OpenPyXL versions.
- Record voxel calibration and channel mapping for each acquisition.
- Keep raw images outside the repository and preserve their checksums in a private lab record.
- Version every change to segmentation, thresholds, inclusion rules, QC rules, or statistical comparisons.
- Back up results before a rerun.
- Report whether distances are raw, signed, floored, or otherwise transformed for plotting.
- Run the figure QA tools after every layout-affecting change.
