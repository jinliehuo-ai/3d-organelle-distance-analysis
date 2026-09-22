# 3D centrosome–nucleus distance analysis

Reproducible Fiji/Python templates for measuring calibrated 3D distances from punctate marker objects to segmented nuclear surfaces in multichannel z-stacks.

This repository contains code and documentation only. It intentionally does not contain experimental images, result tables, microscope metadata, private paths, or derived figures.

## Included

- Fiji/Groovy templates for 3D object segmentation, single-cell crops, and surface-distance measurement.
- Python template for descriptive statistics, two-sided Mann–Whitney tests with Holm correction, and publication figures; group names and comparisons come from a JSON config, not from the code.
- Input/output schema, reproducibility notes, and a public-release safety check.

## Excluded by design

Raw microscopy files (`.nd2`, `.tif`, `.tiff`), result tables (`.xls`, `.xlsx`, `.csv`), metadata, QC images, backups, and local machine paths are excluded. Use a local ignored data directory.

## Minimal workflow

1. Install Fiji with Bio-Formats and 3D ImageJ Suite.
2. Put local input images under an ignored data directory.
3. Edit the clearly marked input/output block in the Fiji template.
4. Run Fiji headless or from the Script Editor.
5. Write a JSON config with your group names, colours and comparisons (see `docs/data_schema.md`).
6. Run the Python figure template against a local summary workbook or tab-delimited summary.

Any change to thresholds, inclusion, outlier handling, or statistical comparisons must be documented as a new version.

## Privacy check

Before publishing code, run:

```bash
python scripts/qa/check_public_release.py .
```

That alone checks only file types, file sizes and private machine paths. Identifiers
in source text - unpublished line or clone names, internal batch numbers, collaborator
names - need a denylist, which must live **outside** this repository, because a
denylist committed here would publish the very terms it hides:

```bash
python scripts/qa/check_public_release.py . --deny-file ../private-denylist.txt
# or export PUBLIC_RELEASE_DENY_FILE=../private-denylist.txt
```

With that variable set, `pytest tests/` gates the same check in CI. Note that removing
an identifier in a new commit does not remove it from git history: if one was ever
committed, rewrite history before the repository becomes public.

## License

MIT. See `LICENSE`.
