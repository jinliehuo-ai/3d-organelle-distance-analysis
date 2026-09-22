"""Four-group distance figure: descriptive statistics, Mann-Whitney tests, publication figure.

This is a template. Group names, colours, the comparison family and the output
base name are supplied by a JSON config (--config); nothing about any particular
experiment is hard-coded here. Run with no config to use neutral placeholder
names (Group A-D) so a fresh clone produces a runnable, identifier-free example.

Input is a per-observation summary table - an .xlsx workbook (a sheet named by
--sheet when present, else the first sheet) or tab-delimited text - carrying at
least file_name, distance_to_nucleus_surface_um and qc_status. Observations are
assigned to groups by the parent directory name of the matching input image.

Measurement policy: every observation is retained. Negative raw distances and
ENGULFED measurements are set to zero only when summarising and plotting; the
source table is never modified.

Statistics: two-sided Mann-Whitney U on the configured comparisons, Holm-adjusted
jointly across that family. Welch t is recorded as a supplementary column only and
never determines the figure annotation.

The unit of analysis is whatever one row of the summary represents. If several
rows come from one field of view, dish or animal, they are not independent and a
per-observation test will overstate significance; aggregate to the replicate level
before testing. State the unit in the methods record either way.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator, PercentFormatter
import numpy as np
import openpyxl
from scipy import stats as scipy_stats
from figure_common import apply_style, describe, write_csv, panel_heading, MM

BASE = Path(__file__).resolve().parents[2]
SKILL = Path(os.environ["NATURE_FIGURE_SKILL_SCRIPTS"]) if os.environ.get("NATURE_FIGURE_SKILL_SCRIPTS") else None
# Neutral defaults. Override every one of these with --config; see docs/data_schema.md.
DEFAULT_CONFIG = {
    "name": "four_group_distance_figure",
    "sheet": None,
    "groups": ["Group A", "Group B", "Group C", "Group D"],
    "colors": {"Group A": "#0072B2", "Group B": "#D55E00",
               "Group C": "#CC79A7", "Group D": "#009E73"},
    "comparisons": [["Group A", "Group B"], ["Group B", "Group C"], ["Group B", "Group D"],
                    ["Group A", "Group C"], ["Group A", "Group D"]],
}
GROUPS = tuple(DEFAULT_CONFIG["groups"])
COLORS = dict(DEFAULT_CONFIG["colors"])
NAME = DEFAULT_CONFIG["name"]
PAIRS = tuple(tuple(pair) for pair in DEFAULT_CONFIG["comparisons"])
SHEET = DEFAULT_CONFIG["sheet"]


def load_config(path):
    """Install group names, colours, comparisons and output name from JSON."""
    global GROUPS, COLORS, NAME, PAIRS, SHEET
    config = dict(DEFAULT_CONFIG)
    if path is not None:
        config.update(json.loads(Path(path).read_text()))
    groups = list(config["groups"])
    if len(groups) != 4 or len(set(groups)) != 4:
        raise ValueError("config 'groups' must list four distinct group names")
    colors = dict(config["colors"])
    missing = [g for g in groups if g not in colors]
    if missing:
        raise ValueError(f"config 'colors' has no entry for {missing}")
    pairs = [tuple(pair) for pair in config["comparisons"]]
    for first, second in pairs:
        if first not in groups or second not in groups:
            raise ValueError(f"comparison ({first}, {second}) names an unknown group")
    if len(set(map(frozenset, pairs))) != len(pairs):
        raise ValueError("config 'comparisons' repeats a pair")
    if not str(config["name"]).strip():
        raise ValueError("config 'name' must be a non-empty output base name")
    GROUPS, COLORS, PAIRS = tuple(groups), colors, tuple(pairs)
    NAME, SHEET = str(config["name"]), config.get("sheet")
    return config


def significance_label(p):
    for threshold, label in ((0.0001, "****"), (0.001, "***"), (0.01, "**"), (0.05, "*")):
        if p < threshold:
            return label
    return "ns"


def pairwise_tests(groups):
    """Every configured comparison forms one Holm family; each pair is counted once."""
    results = []
    for first, second in PAIRS:
        x, y = groups[first], groups[second]
        mw = scipy_stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic", use_continuity=True)
        welch = scipy_stats.ttest_ind(x, y, equal_var=False, alternative="two-sided")
        results.append(dict(group_1=first, group_2=second, n_1=len(x), n_2=len(y),
            mw_u=float(mw.statistic), mw_p_raw=float(mw.pvalue),
            mw_p_holm=0.0, welch_t=float(welch.statistic), welch_p_raw=float(welch.pvalue),
            alternative="two-sided", mw_method="asymptotic; tie correction; continuity correction",
            correction=f"Holm; one family of {len(PAIRS)} MW tests", family_size=len(PAIRS)))
    order = np.argsort([r["mw_p_raw"] for r in results], kind="stable")
    adjusted = np.minimum(1.0, np.maximum.accumulate([
        (len(order)-i)*results[index]["mw_p_raw"] for i,index in enumerate(order)]))
    for index,p in zip(order,adjusted):
        results[index]["mw_p_holm"] = float(p)
    return results


def read_summary_rows(summary):
    """Summary rows as plain strings, from an .xlsx workbook or tab-delimited text."""
    if summary.suffix.lower() in (".xlsx", ".xlsm"):
        book = openpyxl.load_workbook(summary, read_only=True, data_only=True)
        sheet = book[SHEET] if SHEET and SHEET in book.sheetnames else book[book.sheetnames[0]]
        grid = [row for row in sheet.iter_rows(values_only=True)]
        book.close()
        header = [str(cell) for cell in grid[0]]
        assert header[0] == "file_name", f"Unexpected summary header: {header[:3]}"
        rows = []
        for row in grid[1:]:
            if row[0] is None:
                continue
            rows.append({name: ("" if value is None else str(value))
                         for name, value in zip(header, row)})
        return rows
    with summary.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_data(summary, raw_dir):
    mapping = {}
    for path in sorted(raw_dir.rglob("*")):
        # Mirrors run configuration 1.0.3: source fields and candidate crops are
        # not formal four-group inputs.
        if "combined_data" in path.parts or "split_data" in path.parts:
            continue
        if path.suffix.lower() not in (".tif", ".tiff"):
            continue
        if path.name in mapping:
            raise ValueError(f"Duplicate input basename: {path.name}")
        if path.parent.name not in GROUPS:
            raise ValueError(f"Unknown input group: {path.parent.name}")
        mapping[path.name] = path.parent.name
    rows = read_summary_rows(summary)
    names = [r["file_name"] for r in rows]
    assert len(set(names)) == len(names) and set(names) == set(mapping), "Input/summary mismatch"
    for r in rows:
        r["group"] = mapping[r["file_name"]]
        r["raw_distance"] = float(r["distance_to_nucleus_surface_um"])
        assert np.isfinite(r["raw_distance"])
        r["distance"] = 0.0 if "ENGULFED" in r["qc_status"] else max(0.0, r["raw_distance"])
        r["hollow"] = any(s in r["qc_status"] for s in ("MULTIPLE_NUCLEI", "IMPLAUSIBLY_FAR"))
    assert all(any(r["group"] == g for r in rows) for g in GROUPS)
    return rows


def spread_points(ax, x, values, separation_pt=1.8):
    """Greedy horizontal packing in final display coordinates; never alter y."""
    pixels = separation_pt*ax.figure.dpi/72
    rendered_y = ax.transData.transform(np.column_stack((np.zeros(len(values)), values)))[:, 1]
    assigned = np.zeros(len(values))
    done = []
    # Candidates have a stable, data-independent ordering around the group center.
    candidates = [0] + [sign*i*pixels for i in range(1, 24) for sign in (-1, 1)]
    for i in np.argsort(values, kind="stable"):
        for candidate in candidates:
            if all((candidate-assigned[j])**2+(rendered_y[i]-rendered_y[j])**2 >= pixels**2 for j in done):
                assigned[i] = candidate
                done.append(i)
                break
        else:
            raise ValueError("Point packing exceeded candidates")
    center = ax.transData.transform((x, 0))[0]
    positions = ax.transData.inverted().transform(np.column_stack((center+assigned, rendered_y)))[:, 0]
    assert np.max(np.abs(positions-x)) < 0.32, "Widen plot before packing more points"
    return positions


def main(summary, raw_dir, output):
    rows = read_data(summary, raw_dir)
    groups = {g: np.array([r["distance"] for r in rows if r["group"] == g]) for g in GROUPS}
    stats = {g: describe(v) for g,v in groups.items()}
    tests = pairwise_tests(groups)
    output.mkdir(parents=True, exist_ok=True)
    qa = output/"figure_qa"
    qa.mkdir(exist_ok=True)
    source = output/"figure_source_data"
    write_csv(source/f"{NAME}_observations.csv", [dict(file_name=r["file_name"], group=r["group"],
        raw_distance_um=r["raw_distance"], plotted_distance_um=r["distance"],
        qc_status=r["qc_status"], hollow_marker=r["hollow"]) for r in rows])
    write_csv(output/f"descriptive_statistics_{NAME}.csv", [dict(group=g, **stats[g]) for g in GROUPS])
    write_csv(output/f"pairwise_statistics_{NAME}.csv", tests)
    apply_style()
    # Explicit export contract: editable vector text and a 600-dpi raster preview.
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
    width_mm = 183
    fig,(a,b) = plt.subplots(1,2,figsize=(width_mm*MM,175*MM))
    fig.subplots_adjust(left=0.091,right=0.983,bottom=0.325,top=0.741,wspace=0.34)
    a.set_yscale("symlog",linthresh=0.3,linscale=0.65)
    low = min(-0.1,min(s["mean_minus_sd"] for s in stats.values())*1.30)
    high = max(max(v.max(),stats[g]["mean_plus_sd"]) for g,v in groups.items())*1.25
    a.set_ylim(low, high)
    ticks = [v for v in (-10,-5,-2,-1,-0.5,0,0.1,0.3,0.5,1,2,3,5,10,20,50,100) if low <= v <= high]
    a.yaxis.set_major_locator(FixedLocator(ticks))
    a.yaxis.set_major_formatter(FuncFormatter(lambda v,_: f"{v:g}".replace("-","−")))
    a.yaxis.set_minor_locator(NullLocator())
    a.axhspan(low,0,facecolor="#F2F2F2",edgecolor="none",zorder=0)
    a.axhline(0,color="#A0A0A0",linewidth=0.6,zorder=1)
    a.set_xlim(-0.35,3.58)
    a.set_xticks(np.arange(4)+0.10,[f"{g}\nn = {stats[g]['n']}\n{stats[g]['mean']:.2f} ± {stats[g]['sd']:.2f}" for g in GROUPS])
    a.tick_params(axis="x",length=0,pad=6,labelsize=6)
    a.set_ylabel("Centrosome–nuclear surface distance (µm)",labelpad=5)
    # Reserve a physical annotation band above both data rectangles.
    for ax,letter,title,subtitle in ((a,"a","Individual cells","All cells; mean ± SD"),
            (b,"b","Distance distribution","0.5-µm bins; last bin ≥3 µm; labels: counts")):
        ax.set_label(letter)
        ax.annotate(letter,xy=(0,1),xycoords="axes fraction",xytext=(-17,3),
            textcoords="offset points",fontsize=9,fontweight="bold",va="bottom")
        ax.annotate(title,xy=(0,1),xycoords="axes fraction",xytext=(0,104),
            textcoords="offset points",fontsize=8,fontweight="bold",va="bottom")
        ax.annotate(subtitle,xy=(0,1),xycoords="axes fraction",xytext=(0,90),
            textcoords="offset points",fontsize=6.5,va="bottom")
    fig.canvas.draw()
    point_positions = []
    for i,g in enumerate(GROUPS):
        v=groups[g]
        hollow=np.array([r["hollow"] for r in rows if r["group"]==g])
        positions=spread_points(a,i,v)
        point_positions += [{"group":g,"x":float(x),"distance_um":float(y)} for x,y in zip(positions,v)]
        for mask,fill in ((~hollow,COLORS[g]),(hollow,"none")):
            a.scatter(positions[mask],v[mask],s=6.25,facecolors=fill,edgecolors=COLORS[g],linewidths=0.55,zorder=3)
        sx=i+0.40
        a.errorbar(sx,stats[g]["mean"],yerr=stats[g]["sd"],fmt="none",ecolor="#252525",elinewidth=0.8,capsize=2,zorder=4)
        a.plot([sx-0.075,sx+0.075],[stats[g]["mean"]]*2,color="#252525",linewidth=1.6,zorder=5)
    bracket_records=[]
    height_pt=a.get_window_extent().height*72/fig.dpi
    # Shorter intervals below longer intervals; five separate, evenly spaced levels.
    bracket_order=sorted(tests,key=lambda t: (GROUPS.index(t["group_2"])-GROUPS.index(t["group_1"]),GROUPS.index(t["group_1"])))
    for level,test in enumerate(bracket_order):
        x1=GROUPS.index(test["group_1"])+0.1
        x2=GROUPS.index(test["group_2"])+0.1
        y=1+(8+16*level)/height_pt
        a.plot([x1,x1,x2,x2],[y-3/height_pt,y,y,y-3/height_pt],
            transform=a.get_xaxis_transform(),clip_on=False,color="#252525",linewidth=0.7)
        label=significance_label(test["mw_p_holm"])
        a.text((x1+x2)/2,y+2/height_pt,label,transform=a.get_xaxis_transform(),
            ha="center",va="bottom",fontsize=7,clip_on=False)
        bracket_records.append(dict(group_1=test["group_1"],group_2=test["group_2"],
            mw_p_holm=test["mw_p_holm"],label=label,level=level))
    labels=["0–0.5","0.5–1","1–1.5","1.5–2","2–2.5","2.5–3","≥3"]
    edges=np.arange(0,3.01,0.5)
    bins=[]
    max_pct=0
    b.axvspan(-0.5,0.5,facecolor="#F2F2F2",edgecolor="none",zorder=0)
    for i,g in enumerate(GROUPS):
        counts=np.array([np.count_nonzero((groups[g]>=lo)&(groups[g]<hi)) for lo,hi in zip(edges[:-1],edges[1:])]+[np.count_nonzero(groups[g]>=3)])
        assert counts.sum()==len(groups[g])
        pct=counts*100/len(groups[g])
        max_pct=max(max_pct,pct.max())
        bars=b.bar(np.arange(7)+(i-1.5)*0.21,pct,width=0.19,color=COLORS[g],linewidth=0,zorder=3)
        for bar,count,p in zip(bars,counts,pct):
            if count:
                b.text(bar.get_x()+bar.get_width()/2,p+0.8,str(count),ha="center",va="bottom",fontsize=5.5)
        bins += [dict(group=g,bin_um=label,count=int(c),percentage=float(p)) for label,c,p in zip(labels,counts,pct)]
    b.set_xlim(-0.65,6.6)
    b.set_ylim(0,max_pct+9)
    b.set_xticks(np.arange(7),labels,rotation=45,ha="right",rotation_mode="anchor")
    b.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    b.set_xlabel("Distance interval (µm)",labelpad=7)
    b.set_ylabel("Cells within group (%)",labelpad=5)
    # Matplotlib fills legend columns first; reorder handles for row-wise display.
    legend_groups = (GROUPS[0], GROUPS[2], GROUPS[1], GROUPS[3])
    b.legend(handles=[Patch(facecolor=COLORS[g],label=g) for g in legend_groups],loc="upper right",ncol=2,
        handlelength=1,columnspacing=0.8,fontsize=6,borderaxespad=0.3)
    fig.text(0.091,0.228,"Two-sided Mann–Whitney U; Holm correction across five comparisons (exact P values in CSV).",fontsize=6.5)
    fig.text(0.091,0.203,"* P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001; ns P ≥ 0.05 (adjusted P).",fontsize=6.5)
    fig.text(0.091,0.178,"Black bars: mean ± sample SD (µm); values below groups. Brackets show the five specified comparisons.",fontsize=6.5)
    fig.legend(handles=[Line2D([],[],color="#606060",marker="o",markerfacecolor="none",linestyle="none",markersize=3,
        label="Hollow dots: multiple nuclei or distance >5 µm; all cells retained")],loc="lower left",bbox_to_anchor=(0.082,0.132),fontsize=6.5)
    fig.text(0.091,0.105,"a  Grey region: SD below zero, not negative cell distances. Symlog axis (linear near 0).",fontsize=6.5,color="#606060")
    fig.text(0.091,0.077,"b  Grey bin: <0.5 µm (estimated axial resolution); percentages use each group's full n.",fontsize=6.5,color="#606060")
    fig.text(0.091,0.049,f"n denotes analyzed cells/images, not confirmed biological replicates. All {len(rows)} cells are shown.",fontsize=6.5,color="#606060")
    fig.text(0.091,0.021,"Negative / engulfed measurements are set to 0 for display and summaries; raw measurements are retained.",fontsize=6.5,color="#606060")
    write_csv(source/f"{NAME}_bins.csv",bins)
    fig.canvas.draw()
    # The panel-alignment gate is an optional external QA tool. A clone without it still
    # renders the figure; it must never be the reason a fresh checkout cannot run.
    alignment_verdict = "skipped: NATURE_FIGURE_SKILL_SCRIPTS not set"
    if SKILL and SKILL.exists():
        sys.path.insert(0, str(SKILL))
        try:
            from audit_panel_alignment import require_matplotlib_panel_alignment
        except ImportError:
            alignment_verdict = f"skipped: audit_panel_alignment not found in {SKILL}"
        else:
            require_matplotlib_panel_alignment(fig, json_out=qa/f"{NAME}.alignment.json",
                tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=True, strict=True)
            alignment_verdict = "checked"
    print(f"Panel alignment: {alignment_verdict}")
    fig.savefig(output/f"{NAME}.pdf")
    fig.savefig(output/f"{NAME}.svg")
    fig.savefig(output/f"{NAME}.png",dpi=600)
    plt.close(fig)
    report=dict(config_name=NAME,panel_alignment=alignment_verdict,groups=stats,tests=tests,brackets=bracket_records,
        source_summary_sha256=hashlib.sha256(summary.read_bytes()).hexdigest(),rows=len(rows),excluded=0,
        axis_limits=[low,high],axis="symlog: linthresh=0.3, linscale=0.65",point_positions=point_positions,
        hollow_counts={g:sum(r["hollow"] for r in rows if r["group"]==g) for g in GROUPS})
    (qa/f"{NAME}.statistics.json").write_text(json.dumps(report,indent=2)+"\n")
    pair_text = ", ".join(f"{a} vs {b}" for a, b in PAIRS[:-1]) + f" and {PAIRS[-1][0]} vs {PAIRS[-1][1]}"
    group_text = "; ".join(f"{g}, n = {stats[g]['n']}" for g in GROUPS)
    caption = f"""# {NAME} | Distance to segmented surface, four groups

**a**, One marker per analyzed observation, with black bars showing the mean and
full sample SD (ddof = 1). {group_text}. Hollow markers flag observations with
multiple parent objects or distances beyond the plausibility limit; these
observations remain included. The symlog axis is linear near zero. Grey shading
below zero accommodates a negative mean-minus-SD limit and does not represent
negative measured distances. Negative raw distances and ENGULFED measurements are
set to zero only when calculating summaries and plotting; the source table is
unchanged. **b**, Within-group percentages in 0.5-unit intervals, with the final
interval open-ended. Counts appear above bars. Intervals are left-closed and
right-open. No observation is excluded.

The {len(PAIRS)} configured comparisons ({pair_text}) use two-sided Mann-Whitney U
tests. Holm correction is applied jointly across these {len(PAIRS)} tests. Brackets
above panel a show every comparison using the Holm-adjusted P: * P < 0.05,
** P < 0.01, *** P < 0.001, **** P < 0.0001, and ns P >= 0.05; exact adjusted P
values are retained in the CSV statistics table. Unadjusted Mann-Whitney and
supplementary unadjusted Welch t-test results are also retained there; Welch
results do not determine the figure annotation.

Each test treats one summary row as one observation. Replicate structure is a
property of the experiment, not of this script: if several rows share a field of
view, dish, animal or session they are not independent, and the adjusted P values
above will overstate significance. Aggregate to the replicate level before testing
when that applies, and state the unit of analysis and the replicate count in the
methods record.
"""
    (output/f"{NAME}_caption.md").write_text(caption)
    if SKILL and SKILL.exists():
        for tool,args in [
            ("validate_figure.py",[str(Path(__file__).resolve()),"--json"]),
            ("audit_pdf_text.py",[str(output/f"{NAME}.pdf"),"--min-pt","5","--json"]),
            ("audit_figure_collisions.py",[str(output/f"{NAME}.pdf"),"--json-out",str(qa/f"{NAME}.collision-audit.json"),
             "--overlay-pdf",str(qa/f"{NAME}.collision-audit.pdf")])]:
            run=subprocess.run([sys.executable,str(SKILL/tool),*args],text=True,capture_output=True)
            (qa/f"{tool.removesuffix('.py')}.log").write_text(run.stdout+run.stderr)
            if run.returncode:
                raise RuntimeError(f"{tool} failed: {run.stdout} {run.stderr}")
    else:
        (qa/"qa_tools_skipped.log").write_text(
            "Set NATURE_FIGURE_SKILL_SCRIPTS to the directory containing the optional figure QA tools.\n"
        )
    print(json.dumps(stats,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=BASE/"results/summary.xlsx",
                        help="Per-observation summary: .xlsx workbook or tab-delimited text")
    parser.add_argument("--raw-dir", type=Path, default=BASE/"data/raw",
                        help="Input image tree; the parent directory name assigns the group")
    parser.add_argument("--config", type=Path, default=None,
                        help="JSON with name, sheet, groups, colors and comparisons; "
                             "omit to use the neutral Group A-D placeholders")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Use staging; back up existing results before publication")
    args = parser.parse_args()
    load_config(args.config)
    main(args.summary, args.raw_dir, args.output_dir)
