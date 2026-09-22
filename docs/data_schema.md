# Local data schema

The public repository contains no observations. A local analysis directory may contain:

```text
data/
├── raw/
│   ├── GROUP_A/*.tif
│   ├── GROUP_B/*.tif
│   └── ...
└── combined_data/       optional multichannel source stacks
```

The distance summary is a tab-delimited table with one row per analyzed image. Required columns:

- `file_name`
- `distance_to_nucleus_surface_um`
- `negative_distance`
- `qc_status`

Additional diagnostic columns may be retained. The Python figure template maps the immediate parent directory to the group label and refuses unknown groups or duplicate basenames.

## Figure config

Group names, colours, the comparison family, the output base name and the workbook
sheet are supplied as JSON, so no experiment-specific label is stored in the code:

```json
{
  "name": "my_experiment_figure",
  "sheet": "Summary",
  "groups": ["Group A", "Group B", "Group C", "Group D"],
  "colors": {"Group A": "#0072B2", "Group B": "#D55E00",
             "Group C": "#CC79A7", "Group D": "#009E73"},
  "comparisons": [["Group A", "Group B"], ["Group B", "Group C"],
                  ["Group B", "Group D"], ["Group A", "Group C"],
                  ["Group A", "Group D"]]
}
```

```bash
python scripts/plotting/make_four_group_figure.py \
  --summary results/summary.xlsx --raw-dir data/raw \
  --config config.local.json --output-dir out/
```

Group names must match the `data/raw/` subdirectory names exactly, and every
comparison must name configured groups; the loader rejects anything else. Keep the
config out of version control if the group names are unpublished (`config.local.*`
is already ignored). Run without `--config` to get neutral Group A-D placeholders.
