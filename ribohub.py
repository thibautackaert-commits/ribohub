"""
Build a Gwips RiboSeq track hub from sorted bigWig files.

Layout expected on disk:
    {data_dir}/{srr[:6]}/{srr[6:8]}/{filename}

Filenames follow:
    {SRR}_pshifted_{forward,reverse}.{bw,bigWig}            (all reads)
    {SRR}_pshifted_unique_{forward,reverse}.{bw,bigWig}     (unique only)
    {SRR}_pshifted_multimapped_{forward,reverse}.{bw,bigWig} (multimapped only)

bigBed region files (optional):
    {SRR}_*.{bb,bigBed}   (no strand/kind suffix; region label from filename)

Track structure in the hub:
    - One main composite containing every (sample, strand, kind) subtrack,
      with subGroups for UCSC filter-matrix browsing.
    - Per-sample strand-overlay aggregates as siblings of the composite
      (using the bare _forward/_reverse files), enabled via --with-aggregates.
    - An optional region SuperTrack containing bigBed tracks, as a SIBLING
      of the signal composite, enabled via --with-regions. 

Output formats:
    directory   : standard multi-file hub (hub.txt + genomes.txt + hg38/trackDb.txt + symlinks)
    single-file : useOneFile hub (everything in one .hub.txt)
"""
# =============================================================================
# READING THIS FILE TOP-TO-BOTTOM
# -----------------------------------------------------------------------------
# The program flows top-to-bottom in stages, marked by banner comments:
#   INPUT PARSING  -> turn what the user typed into clean data
#   METADATA       -> read + tidy the spreadsheet of sample info
#   FILTER         -> optionally derive/narrow the sample set from metadata
#   DISCOVERY      -> find the actual bigWig/bigBed files on disk
#   COMPOSITE      -> assemble the main filterable track group (the core)
#   AGGREGATE      -> assemble the per-sample overlay tracks
#   REGIONS        -> assemble the bigBed region track group
#   DESCRIPTION    -> build the static HTML explainer page
#   OUTPUT WRITERS -> turn the assembled objects into track-hub text
#   CLI            -> the command-line interface that starts everything
#
# A leading underscore on a name (e.g. _has_value) marks an INTERNAL helper:
# used only within this file, not meant to be called from outside.
#
# For a full plain-English walkthrough, see README.md. #TODO
# =============================================================================

# --- stdlib ---
import csv          
import logging      
import re           
import shutil       
import sys          
from dataclasses import dataclass        
from pathlib import Path                 
from typing import Literal, NamedTuple   

# --- third-party ---
import click            # builds the command-line interface (the --samples flags)
import trackhub  # type: ignore  # the library that writes the track-hub text for us


# ----- type aliases --------------------------------------------------------
# A `Literal` says a value can ONLY be one of these exact strings. So a Strand
# is either the text "forward" or "reverse", nothing else. This documents the
# fixed set of choices and lets tools catch typos like "forwrd".
Strand = Literal["forward", "reverse"]
Kind = Literal["unique", "multimapped"]  # None separately means "all reads / bare file", atleast so it can be a FileKey value

class FileKey(NamedTuple):
    """Dict key identifying one bigWig variant for a sample.

    A NamedTuple is a small, UNCHANGEABLE bundle of named values, like a
    labelled box with fixed compartments. Used here as a "coordinate":
    (strand, kind) together identify one specific kind of file.
    """
    strand: Strand
    kind: Kind | None  # None means "all reads" (no kind suffix in filename)


class BigwigEntry(NamedTuple):
    """One classified bigWig file: which strand, which kind, where it lives.

    Like FileKey but also carries the file's location. FileKey is used as a
    lookup key; BigwigEntry is the full result of inspecting one file.
    """
    strand: Strand
    kind: Kind | None
    path: Path


class BigBedEntry(NamedTuple):
    """One classified bigBed region file: its region label and path.

    bigBed files carry per-feature strand internally, so they have no
    strand/kind suffix. The label is derived from the text between
    "{SRR}_" and the file extension (e.g. "SRR123_orfs.bb" → "orfs").
    """
    label: str
    path: Path


@dataclass(frozen=True)
class BuildContext:
    """Run-wide configuration passed to the composite builder.

    One box holding ALL the settings for a single run, so functions take one
    `ctx` argument instead of ten separate ones. `frozen=True` makes it
    read-only: the settings can't be changed once the run starts.
    """
    base_url: str
    colors: dict[str, str]
    label_fields: tuple[str, ...]          # for individual track long labels
    subgroup_fields: tuple[str, ...]       # metadata columns that become subGroups
    metadata_fields: tuple[str, ...]       # metadata columns that become 'metadata' lines
    kinds: frozenset[str]                  # which kinds to include: subset of {all, unique, multi}
    with_aggregates: bool
    with_regions: bool                     # whether to include bigBed region tracks
    auto_scale: bool                       # True = UCSC autoScale; False = fixed viewLimits

# ----- defaults (informative; runtime values come from CLI) ----------------

DEFAULT_COLORS: dict[str, str] = {
    "fwd":       "#003399",   # navy blue
    "rev":       "#CC5500",   # burnt orange
    "fwd_multi": "#336699",   # steel blue
    "rev_multi": "#E68A00",   # amber
}
DEFAULT_REGION_COLOR: str = "#2E8B57"  # sea green, distinct from signal colors

VIEW_LIMITS: str = "-127:127"

# Fixed metadata columns used for track labels, subgroup filter dimensions,
# and per-track metadata lines. Hardcoded to RiboSeqOrg column names.
# Edit here if the upstream schema changes.
# Possible future improvement: make these configurable via CLI
LABEL_FIELDS: tuple[str, ...] = (
    "ScientificName",
    "LIBRARYTYPE",
    "TISSUE",
    "CELL_LINE",
    "CONDITION",
    "INHIBITOR",
)
SUBGROUP_FIELDS: tuple[str, ...] = (
    "TISSUE",
    "CELL_LINE",
    "CONDITION",
    "INHIBITOR",
    "Stress",
    "GENE",
    "FRACTION",
)
# SRR ID column name in RiboSeqOrg metadata CSV.
SRR_COLUMN: str = "Run"

# Values that RiboSeqOrg uses to mean "no data". Case-insensitive.
NULL_VALUES: frozenset[str] = frozenset({
    "", "nana", "na", "n/a", "none", "null", "nan", "-", "unknown", "0.0", "missing",
})

# Tag values valid in UCSC subGroups must be identifiers. This regex matches
# what _sanitize_tag produces.
# In plain terms: this pattern matches "any run of characters that are NOT
# lowercase letters, digits or underscores", specifically the messy characters we
# want to scrub out of a label. The r"..." is a "raw string" (backslashes are
# left as-is). The ^ inside [...] means "NOT these".
_TAG_SAFE_RE = re.compile(r"[^a-z0-9_]+")

# ----- bigBed filename convention (PROVISIONAL) ---------------------------
# This regex extracts the region label from a bigBed filename.
# Pattern: {SRR}_{label}.bb or {SRR}_{label}.bigBed
# The label is the text between the first underscore after the SRR and the
# extension. PROVISIONAL: change this ONE line if the naming convention changes.
_BIGBED_LABEL_RE = re.compile(r"^[A-Z]{3}\d+_(.+)\.(bb|bigBed)$")

log = logging.getLogger("ribohub")


# ----- input parsing -------------------------------------------------------

def parse_samples(value: str) -> set[str]:
    """Parse the --samples CLI value into a set of SRR IDs.

    Accepts four input formats:
        - Single ID:        "SRR9295900"
        - Comma-separated:  "SRR9295900,SRR9295901"
        - Text file:        path to a .txt with one ID per line
        - CSV file:         path to a .csv; reads the first column, skips the header

    Raises click.BadParameter if no valid IDs are found.
    """
    path = Path(value)
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as f:
            if path.suffix.lower() == ".csv":
                # CSV: first column, skip header row.
                ids = {row[0].strip() for row in csv.reader(f)
                       if row and row[0].strip()}
            else:
                # Plain text file: one ID per line.
                ids = {line.strip() for line in f if line.strip()}
    else:
        # Not a file path, treat as comma-separated IDs.
        ids = {s.strip() for s in value.split(",") if s.strip()}

    # Fail fast if nothing was parsed.
    if not ids:
        raise click.BadParameter("No sample IDs found in --samples value.")
    return ids


def hex_to_trackhub_rgb(value: str) -> str:
    """Validate hex color and convert to trackhub RGB format ('R,G,B')."""
    v = value.lstrip("#")
    if len(v) != 6 or not all(c in "0123456789abcdefABCDEF" for c in v):
        raise click.BadParameter(f"Invalid hex color: {value!r}. Expected #RRGGBB.")
    return trackhub.helpers.hex2rgb(value)


def parse_csv_field(value: str) -> tuple[str, ...]:
    """Split a comma-separated CLI string into a de-duplicated tuple.

    Preserves input order (first occurrence wins). Used for --kinds
    where ordering affects UCSC filter-matrix rendering.
    """
    fields: list[str] = []
    seen: set[str] = set()
    for raw in value.split(","):
        col = raw.strip()
        if not col or col in seen:
            continue
        seen.add(col)
        fields.append(col)
    return tuple(fields)


@dataclass
class FilterResult:
    """Result of applying --filter against loaded metadata.

    Returned by apply_metadata_filter() so every caller (CLI, Galaxy wrapper) can format its own error messages from the same data, rather
    than receiving pre-formatted strings it can't adapt.

    Fields
    ------
    matched         SRR IDs that passed ALL filter conditions.
    per_field       How many SRRs matched each individual field condition.
                    e.g. {"CONDITION": 3, "CELL_LINE": 5}
    zero_fields     Field names where no SRR matched (likely a typo or wrong value).
    available_values
                    For every zero_field, the distinct non-null values that DO
                    exist in the metadata, shown to the user as a "did you mean?"
                    hint.  Only populated for zero_fields to keep it cheap.
    unknown_columns Column names in the filter that don't exist in the CSV at all.
    available_columns
                    The full sorted list of columns in the metadata CSV,
                    shown alongside unknown_columns so the user can self-correct.
    """
    matched: set[str]
    per_field: dict[str, int]
    zero_fields: list[str]
    available_values: dict[str, list[str]]
    unknown_columns: list[str]
    available_columns: list[str]


# Columns the CLI exposes as suggested filter dimensions (shown in --help).
# These are the most likely to vary across samples in RiboSeqOrg.
SUGGESTED_FILTER_COLUMNS: tuple[str, ...] = (
    "CONDITION", "INHIBITOR", "TISSUE", "CELL_LINE",
    "ScientificName", "LIBRARYTYPE", "REPLICATE",
)


def parse_filter(value: str) -> dict[str, list[str]]:
    """Parse a --filter string into {column: [value, ...]} conditions.

        Syntax:
        COL=VAL             single value
        COL=VAL1/VAL2       OR within a field
        COL=X,COL2=Y        AND across fields (comma-separated)

    Examples
    --------
    "CONDITION=High"                  -> {"CONDITION": ["High"]}
    "CONDITION=High/Test"             -> {"CONDITION": ["High", "Test"]}
    "CONDITION=High,CELL_LINE=HEK293" -> {"CONDITION": ["High"], "CELL_LINE": ["HEK293"]}

    Raises click.BadParameter on malformed input so the CLI prints a clean
    error without a traceback.
    See SUGGESTED_FILTER_COLUMNS for commonly useful fields.
    """
    filters: dict[str, list[str]] = {}
    for pair in value.split(","):          # each comma-separated condition
        pair = pair.strip()
        if not pair:
            continue
        if "=" not in pair:
            raise click.BadParameter(
                f"--filter: invalid token {pair!r}. "
                f"Expected COL=VAL or COL=VAL1/VAL2.\n"
                f"Example: --filter \"CONDITION=High/Test,CELL_LINE=HEK293\""
            )
        col, vals = pair.split("=", 1)     # split on first "=" only, values may contain "="
        col = col.strip()
        parsed_vals = [v.strip() for v in vals.split("/") if v.strip()]  # "/" = OR within field
        if not col or not parsed_vals:
            raise click.BadParameter(
                f"--filter: empty column or value in {pair!r}."
            )
        filters[col] = parsed_vals
    if not filters:
        raise click.BadParameter("--filter: no valid conditions found.")
    return filters


def apply_metadata_filter(
    metadata: dict[str, dict[str, str]],
    filters: dict[str, list[str]],
    csv_columns: list[str],
) -> FilterResult:
    """Apply parsed filter conditions against loaded metadata.

    Case-insensitive matching with whitespace stripping,
    "high", "High", " High " all match a filter value of "High".

    Logic: AND across fields, OR within a field.

    Parameters
    ----------
    metadata     {srr: {col: val}} from load_metadata().
    filters      {col: [val, ...]} from parse_filter().
    csv_columns  Raw CSV header; used for unknown-column detection.

    Returns a FilterResult for the caller to format.
    """
    col_set = set(csv_columns)

    # Check for unknown columns up front, before touching any rows.
    unknown = [col for col in filters if col not in col_set]

    # Per-field match counts and available values for zero-match fields.
    per_field: dict[str, int] = {}
    zero_fields: list[str] = []
    available_values: dict[str, list[str]] = {}

    for col, wanted in filters.items():
        if col in unknown:
            continue  # can't count matches for a column that doesn't exist
        wanted_lower = {v.lower() for v in wanted}
        count = sum(
            1 for row in metadata.values()
            if row.get(col, "").strip().lower() in wanted_lower
        )
        per_field[col] = count
        if count == 0:
            zero_fields.append(col)
            # Collect distinct non-null values for the "did you mean?" hint.
            seen: set[str] = set()
            vals: list[str] = []
            for row in metadata.values():
                v = row.get(col, "").strip()
                if v and v.lower() not in NULL_VALUES and v not in seen:
                    seen.add(v)
                    vals.append(v)
            available_values[col] = sorted(vals)

    # Final intersection: rows that pass ALL non-unknown conditions.
    matched: set[str] = set()
    valid_filters = {col: vals for col, vals in filters.items() if col not in unknown}
    for srr, row in metadata.items():
        if all(
            row.get(col, "").strip().lower() in {v.lower() for v in vals}
            for col, vals in valid_filters.items()
        ):
            matched.add(srr)

    return FilterResult(
        matched=matched,
        per_field=per_field,
        zero_fields=zero_fields,
        available_values=available_values,
        unknown_columns=unknown,
        available_columns=sorted(csv_columns),
    )


def srr_to_dir(srr_id: str, data_dir: Path) -> Path:
    """Map an SRR ID to its on-disk directory. Mirrors sorting.sh.

    Layout: {data_dir}/{srr[:6]}/{srr[6:8]}/
    Example: SRR9295900 → data_dir/SRR929/59/
    """
    return data_dir / srr_id[:6] / srr_id[6:8]


def classify_bigwig(path: Path) -> BigwigEntry | None:
    """Classify a bigWig file by strand and kind based on its filename.

    Looks for _forward/_reverse (required) and _unique/_multimapped
    (optional, absent means all reads). Returns None if strand
    can't be determined.
    """
    fname = path.name
    if "_forward" in fname:
        strand: Strand = "forward"
    elif "_reverse" in fname:
        strand = "reverse"
    else:
        return None

    if "_unique" in fname:
        kind: Kind | None = "unique"
    elif "_multimapped" in fname:
        kind = "multimapped"
    else:
        kind = None

    return BigwigEntry(strand=strand, kind=kind, path=path)


def classify_bigbed(path: Path) -> BigBedEntry | None:
    """Recognize a bigBed file and extract its region label.

    Returns None if the filename doesn't match the expected convention.
    The label is the text between "{SRR}_" and the extension:
      SRR123456_orfs.bb      → label = "orfs"
    """
    match = _BIGBED_LABEL_RE.match(path.name)  # pattern defined at module level
    if not match:
        return None
    label = match.group(1)                      # capture group 1 = the region label
    return BigBedEntry(label=label, path=path)


# ----- metadata ------------------------------------------------------------

def _has_value(value: str | None) -> bool:
    """Return True if value is non-empty and not a known null sentinel."""
    if value is None:
        return False
    return value.strip().lower() 
    if stripped in NULL_VALUES: 
        return False
    if stripped.startswith("nana"):
        return False
    return True


def _clean_value(col: str, val: str) -> str:
    """Normalize known quirks in RiboSeqOrg metadata values.

    Currently handles:
        REPLICATE: strips trailing ".0" (e.g. "1.0" → "1"), caused by
        Excel/pandas treating integer columns as float.
    """
    v = val.strip()
    if not v:
        return v
    if col == "REPLICATE" and v.endswith(".0") and v[:-2].isdigit():
        return v[:-2]  # "1.0" → "1"
    return v


def load_metadata(
    path: Path,
    srr_column: str,
    required_fields: tuple[str, ...] | None = None,
) -> dict[str, dict[str, str]]:
    """Load a metadata CSV into {srr_id: {column: value}} format.

    Parameters
    ----------
    path              Path to the CSV file.
    srr_column        Column name containing SRR IDs (e.g. "Run").
    required_fields   Columns the pipeline needs. Raises click.BadParameter
                      if any are missing from the CSV header.

    Returns {srr_id: {col: val}} with values cleaned via _clean_value().
    Duplicate SRR rows are logged as warnings; last occurrence wins.
    """
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []

        # --- validate header before reading any rows ---
        if srr_column not in header:
            raise click.BadParameter(
                f"Column {srr_column!r} not found in metadata CSV. "
                f"Available columns: {', '.join(header)}"
            )
        if required_fields:
            unknown = [c for c in required_fields if c not in header]
            if unknown:
                raise click.BadParameter(
                    f"Required field(s) not in CSV: {', '.join(unknown)}. "
                    f"Available: {', '.join(header)}"
                )

        # --- read and clean rows ---
        rows: dict[str, dict[str, str]] = {}
        duplicates: list[str] = []
        for raw_row in reader:
            srr = (raw_row.get(srr_column) or "").strip()
            if not srr:
                continue
            row = {k: _clean_value(k, v or "") for k, v in raw_row.items()}
            if srr in rows:
                duplicates.append(srr)
            rows[srr] = row                # last occurrence wins on duplicates

    if duplicates:
        log.warning(
            f"Metadata CSV has {len(duplicates)} duplicate SRR(s); "
            f"keeping last occurrence. Examples: {', '.join(duplicates[:5])}"
        )
    log.info(f"Loaded metadata for {len(rows)} sample(s).")
    return rows

def long_label_for_sample(
    srr_id: str,
    metadata_row: dict | None,
    label_fields: tuple[str, ...],
) -> str:
    """Build a descriptive long label for a track from metadata fields.

    Joins the SRR ID with non-null values from label_fields using " | ".
    Example: "SRR9295900 | Homo sapiens | RFP | HEK293 | High"

    Falls back to just the SRR ID if no metadata is available.
    """
    if not metadata_row:
        return srr_id
    parts = [srr_id]
    for col in label_fields:
        val = metadata_row.get(col)
        if _has_value(val):
            parts.append(val.strip())
    return " | ".join(parts)


# ----- discovery -----------------------------------------------------------

def discover_samples(
    sample_set: set[str],
    data_dir: Path,
    with_regions: bool = False,
) -> tuple[
    dict[str, dict[FileKey, str]],
    dict[str, list[tuple[str, str]]],
    set[str],
    set[str],
]:
    """Scan the data directory for bigWig (and optionally bigBed) files.

    For each SRR ID, computes the expected directory via srr_to_dir(),
    then globs for matching files and classifies them.

    Parameters
    ----------
    sample_set    SRR IDs to look for.
    data_dir      Root of the sorted file tree.
    with_regions  Also scan for bigBed region files.

    Returns
    -------
    found         {srr: {FileKey: rel_path}}, bigWig signal files.
    regions       {srr: [(label, rel_path), ...]}, bigBed region files.
    dir_missing   SRR IDs whose expected directory doesn't exist.
    dir_empty     SRR IDs whose directory exists but has no recognized files.
    """
    found: dict[str, dict[FileKey, str]] = {}
    regions: dict[str, list[tuple[str, str]]] = {}
    dir_missing: set[str] = set()
    dir_empty: set[str] = set()

    for srr_id in sorted(sample_set):
        sample_dir = srr_to_dir(srr_id, data_dir)

        if not sample_dir.is_dir():
            dir_missing.add(srr_id)
            continue

        # --- bigWig discovery ---
        bws = sorted(
            list(sample_dir.glob(f"{srr_id}_*.bw")) +
            list(sample_dir.glob(f"{srr_id}_*.bigWig")) # both extensions are valid
        )

        classified: dict[FileKey, str] = {}
        for path in bws:
            entry = classify_bigwig(path)
            if entry is None:
                log.debug(f"Skipping unrecognized bigWig: {path}")
                continue
            rel = path.relative_to(data_dir).as_posix() # store as relative URL path
            classified[FileKey(entry.strand, entry.kind)] = rel

        if classified:
            found[srr_id] = classified

        # --- bigBed region discovery ---
        if with_regions:
            bbs = sorted(
                list(sample_dir.glob(f"{srr_id}_*.bb")) +
                list(sample_dir.glob(f"{srr_id}_*.bigBed"))
            )
            sample_regions: list[tuple[str, str]] = []
            for path in bbs:
                bb_entry = classify_bigbed(path)
                if bb_entry is None:
                    log.debug(f"Skipping unrecognized bigBed: {path}")
                    continue
                rel = path.relative_to(data_dir).as_posix()
                sample_regions.append((bb_entry.label, rel))
            if sample_regions:
                regions[srr_id] = sample_regions

        # A sample is "empty" only if it has neither signal nor region files
        if not classified and srr_id not in regions:
            if srr_id not in dir_missing:
                dir_empty.add(srr_id)

    return found, regions, dir_missing, dir_empty


def report_discovery(
    found: dict[str, dict],
    regions: dict[str, list],
    dir_missing: set[str],
    dir_empty: set[str],
    sample_set: set[str],
) -> None:
    """Log a summary of discovery results to stdout.

    Shows counts for found, missing and empty samples so the user
    can verify before the build proceeds.
    """
    click.echo(f"Requested: {len(sample_set)} samples")
    click.echo(f"  Signal found:    {len(found)}")
    if regions:
        total_bb = sum(len(v) for v in regions.values())
        click.echo(f"  Region files:    {total_bb} across {len(regions)} sample(s)")
    if dir_missing:
        click.echo(f"  Dir missing:     {len(dir_missing)} → "
                   f"{', '.join(sorted(dir_missing))}")
    if dir_empty:
        click.echo(f"  Dir empty:       {len(dir_empty)} → "
                   f"{', '.join(sorted(dir_empty))}")
    click.echo("")


# ----- composite construction ---------------------------------------------
#
# One composite holds every track across all samples. UCSC's filter matrix
# (subGroup1..N + dimensions) lets users filter by sample/strand/kind/condition.
#
# Built in two passes:
#   1. collect_subgroup_vocabularies(), find which values actually appear
#      per dimension; drop single-value dimensions (useless as filters).
#   2. build_composite_trackdb(), create the composite with subGroup
#      definitions, then add one subtrack per (sample, strand, kind).


def _sanitize_tag(value: str) -> str:
    """Convert an arbitrary string to a UCSC-safe subGroup tag.

    UCSC subGroup tags must be valid identifiers, no whitespace, no special
    chars, ideally lowercase. "High dose" → "high_dose", "Replicate 1" → "replicate_1".
    """
    lowered = value.strip().lower()
    cleaned = _TAG_SAFE_RE.sub("_", lowered).strip("_")  # replace non-alphanumeric runs with "_"
    return cleaned or "unknown"                           # fallback if nothing survives


# Tag -> display label for structural subGroups (shown in UCSC filter UI).
_STRAND_TAG_LABELS: dict[str, str] = {"fwd": "Forward", "rev": "Reverse"}
_KIND_TAG_LABELS: dict[str, str] = {
    "all":    "All_reads",
    "unique": "Unique",
    "multi":  "Multimapped",
}
# FileKey.kind → trackDb subgroup tag. None = all reads (bare file).
_KIND_TO_TAG: dict[Kind | None, str] = {
    None:           "all",
    "unique":       "unique",
    "multimapped":  "multi",
}


def collect_subgroup_vocabularies(
    found: dict[str, dict[FileKey, str]],
    metadata: dict[str, dict[str, str]],
    subgroup_fields: tuple[str, ...],
    kinds: frozenset[str],
) -> dict[str, dict[str, str]]:
    """Collect which subgroup values are present in the data (pass 1 of 2).

    Scans discovered files and metadata to build a vocabulary per dimension.
    Dimensions with only one distinct value are dropped, they can't filter
    anything.

    Returns {dimension: {tag: label}}, ordered: structural dims first
    (sample, strand, kind), then metadata dims in subgroup_fields order.
    """
    vocab: dict[str, dict[str, str]] = {}

    # Structural: sample
    samples = sorted(found.keys())

    # Structural: strand (single-value rule below may drop it)
    strand_tags_present: set[str] = set()
    for files in found.values():
        for fk in files:
            strand_tags_present.add("fwd" if fk.strand == "forward" else "rev")
    vocab["strand"] = {
        tag: _STRAND_TAG_LABELS[tag]
        for tag in ("fwd", "rev") if tag in strand_tags_present
    }

    # Structural: kind. Only include kinds the user asked for via --kinds
    # AND that actually exist in the data.
    kind_tags_present: set[str] = set()
    for files in found.values():
        for fk in files:
            tag = _KIND_TO_TAG[fk.kind]
            if tag in kinds:
                kind_tags_present.add(tag)
    vocab["kind"] = {
        tag: _KIND_TAG_LABELS[tag]
        for tag in ("all", "unique", "multi") if tag in kind_tags_present
    }

    # Metadata-derived dimensions
    for col in subgroup_fields:
        dim = col.lower()  # UCSC dimensions are lowercase; CSV column name kept for label lookup
        seen_tags: dict[str, str] = {}  # tag → label (raw value)
        for srr_id in samples:
            row = metadata.get(srr_id)
            if not row:
                continue
            val = row.get(col)
            if not _has_value(val):
                continue
            tag = _sanitize_tag(val)
            seen_tags.setdefault(tag, val.strip())  # first occurrence sets the display label
        if seen_tags:
            has_missing = any(
                not _has_value((metadata.get(srr_id) or {}).get(col))
                for srr_id in samples
            )
            if has_missing:
                seen_tags["na"] = "N/A"
            vocab[dim] = seen_tags

    # Drop dimensions with only one value, they can't filter anything.
    return {dim: vals for dim, vals in vocab.items() if len(vals) > 1}


def vocab_to_subgroup_definitions(
    vocab: dict[str, dict[str, str]],
) -> list["trackhub.SubGroupDefinition"]:
    """Convert the vocab dict into trackhub SubGroupDefinition objects.

    Translates RiboHub's {dimension: {tag: label}} format into the objects
    the trackhub library needs to render subGroup1..N, dimensions and
    filterComposite lines. Order is preserved from the vocab dict.
    """
    defs: list[trackhub.SubGroupDefinition] = []
    for dim, values in vocab.items():
        label = dim.replace("_", " ").title().replace(" ", "_")         # "cell_line" -> "Cell_Line"
        mapping = {tag: lbl.replace(" ", "_") for tag, lbl in values.items()}  # UCSC tags can't have spaces
        defs.append(trackhub.SubGroupDefinition(
            name=dim, label=label, mapping=mapping,
        ))
    return defs


def _subgroup_assignment_for_subtrack(
    srr_id: str,
    file_key: FileKey,
    metadata_row: dict | None,
    vocab: dict[str, dict[str, str]],
    subgroup_fields: tuple[str, ...],
) -> dict[str, str]:
    """Compute the {dimension: tag} subGroups assignment for one subtrack.

    Only includes dimensions present in vocab (i.e. that survived the
    single-value filter in pass 1).

    Example return: {"sample": "SRR123", "strand": "fwd", "kind": "unique",
                     "condition": "high_dose"}
    """
    out: dict[str, str] = {}

    # --- structural dimensions ---
    if "strand" in vocab:
        out["strand"] = "fwd" if file_key.strand == "forward" else "rev"

    if "kind" in vocab:
        out["kind"] = _KIND_TO_TAG[file_key.kind]

    # --- metadata-derived dimensions ---
    for col in subgroup_fields:
        dim = col.lower()
        if dim not in vocab:
            continue  # dropped as single-valued or never had values
        val = (metadata_row or {}).get(col)
        if _has_value(val):
            out[dim] = _sanitize_tag(val)
        elif dim in vocab:
            out[dim] = "na"

    return out


def _metadata_lines_for_subtrack(
    metadata_row: dict | None,
    metadata_fields: tuple[str, ...],
) -> dict[str, str]:
    """Build {key: value} pairs for the UCSC `metadata` trackDb line.

    This line shows extra info in the browser UI but doesn't affect
    filtering. Skips fields that are missing or null-valued.
    """
    if not metadata_row:
        return {}
    out: dict[str, str] = {}
    for col in metadata_fields:
        val = metadata_row.get(col)
        if _has_value(val):
            out[col.lower()] = val.strip()
    return out


def _color_for(strand: Strand, kind: Kind | None, ctx: BuildContext) -> str:
    """Pick the color for a subtrack based on strand and kind."""
    if kind == "multimapped":
        return ctx.colors["fwd_multi"] if strand == "forward" else ctx.colors["rev_multi"]
    # 'all' and 'unique' share the base strand color in this palette.
    return ctx.colors["fwd"] if strand == "forward" else ctx.colors["rev"]


def build_composite_subtrack(
    srr_id: str,
    file_key: FileKey,
    rel_path: str,
    ctx: BuildContext,
    metadata_row: dict | None,
    vocab: dict[str, dict[str, str]],
) -> tuple[trackhub.Track, dict[str, str]]:
    """Build one subtrack of the main composite.

    Returns (track, metadata_dict). The track gets subGroups via the
    trackhub API. The metadata dict is returned separately because
    trackhub 1.0 can't serialize UCSC `metadata` lines, the writer
    injects them later.
    """
    kind_tag = _KIND_TO_TAG[file_key.kind]
    strand_tag = "fwd" if file_key.strand == "forward" else "rev"

    name = trackhub.helpers.sanitize(f"{srr_id}_{file_key.strand}_{kind_tag}")
    short_label = f"{srr_id} {strand_tag} {kind_tag}"
    long_label = (
        f"{long_label_for_sample(srr_id, metadata_row, ctx.label_fields)} "
        f"({file_key.strand}"
        f"{', ' + kind_tag if kind_tag != 'all' else ''})"
    )

    kwargs: dict[str, str] = dict(
        name=name,
        url=f"{ctx.base_url}/{rel_path}",
        tracktype="bigWig",
        short_label=short_label,
        long_label=long_label,
        color=_color_for(file_key.strand, file_key.kind, ctx),
    )
    if file_key.strand == "reverse":
        kwargs["negateValues"] = "on"

    track = trackhub.Track(**kwargs)

    # subGroups assigned via trackhub's native API
    track.add_subgroups(_subgroup_assignment_for_subtrack(
        srr_id, file_key, metadata_row, vocab, ctx.subgroup_fields,
    ))

    # metadata can't go through trackhub, returned separately for the writer to inject
    metadata_line = _metadata_lines_for_subtrack(metadata_row, ctx.metadata_fields)
    return track, metadata_line  # the function now hands back a PAIR


def build_composite_trackdb(
    found: dict[str, dict[FileKey, str]],
    metadata: dict[str, dict[str, str]],
    ctx: BuildContext,
    composite_name: str = "ribohub",
) -> tuple[trackhub.CompositeTrack, dict[str, dict[str, str]]]:
    """Build the composite track and all its subtracks (pass 2 of 2).

    Attaches subgroup definitions, dimensions and filterComposite via
    the trackhub API. Returns (composite, metadata_map) where metadata_map
    is {track_name: {key: value}} for the writer to inject — the one
    field trackhub 1.0 can't serialize.
    """
    vocab = collect_subgroup_vocabularies(
        found, metadata, ctx.subgroup_fields, ctx.kinds,
    )
    subgroup_defs = vocab_to_subgroup_definitions(vocab)

    comp_kwargs: dict[str, str] = dict(
        name=trackhub.helpers.sanitize(composite_name),
        tracktype="bigWig",
        short_label="All RiboSeq tracks",
        long_label="All RiboSeq samples (filter by subgroups)",
        visibility="full",
        viewLimits=VIEW_LIMITS,
        autoScale="on" if ctx.auto_scale else "off",
        maxHeightPixels="100:50:8",
    )
    if not ctx.auto_scale:
        comp_kwargs["viewLimits"] = VIEW_LIMITS
    comp = trackhub.CompositeTrack(**comp_kwargs)

    # Attach subgroup definitions and derived header lines (dimensions, filterComposite)
    if subgroup_defs:
        comp.add_subgroups(subgroup_defs)
        comp.add_params(
            dimensions=trackhub.helpers.dimensions_from_subgroups(subgroup_defs),
        )
        filt = trackhub.helpers.filter_composite_from_subgroups(subgroup_defs)
        if filt:
            comp.add_params(filterComposite=filt)
    # sortOrder isn't derived by the library, set explicitly
    sort_dims = [d for d in ("strand", "kind", "condition", "tissue", "cell_line", "inhibitor")
                 if d in vocab]
    if sort_dims:
        comp.add_params(sortOrder=" ".join(f"{d}=+" for d in sort_dims))

    # metadata per track, keyed by name, passed to the writer for injection
    metadata_map: dict[str, dict[str, str]] = {}

    n_subtracks = 0
    for srr_id in sorted(found.keys()):
        # Sort files consistently; `or ""` handles None (all-reads) in sorting
        for file_key, rel_path in sorted(
            found[srr_id].items(),
            key=lambda kv: (kv[0].strand, kv[0].kind or ""),
        ):
            kind_tag = _KIND_TO_TAG[file_key.kind]
            if kind_tag not in ctx.kinds:
                continue  # user didn't ask for this kind, skip it
            track, meta_line = build_composite_subtrack(
                srr_id, file_key, rel_path, ctx,
                metadata.get(srr_id), vocab,
            )
            comp.add_subtrack(track)
            if meta_line:
                metadata_map[track.name] = meta_line  # remember it for the writer
            n_subtracks += 1

    if n_subtracks == 0:
        return comp, metadata_map

    log.info(f"Built composite with {n_subtracks} subtrack(s) across "
             f"{len(found)} sample(s).")
    return comp, metadata_map


# ----- aggregate construction ---------------------------------------------
# Per-sample strand-overlay aggregates, rendered as siblings of the composite.

def build_aggregate(
    srr_id: str,
    files: dict[FileKey, str],
    ctx: BuildContext,
    metadata_row: dict | None = None,
) -> trackhub.AggregateTrack | None:
    """Strand-overlay aggregate using bare (all-reads) files."""
    fwd = files.get(FileKey("forward", None))
    rev = files.get(FileKey("reverse", None))
    if not (fwd and rev):
        return None

    rich = long_label_for_sample(srr_id, metadata_row, ctx.label_fields)
    agg = trackhub.AggregateTrack(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg"),
        tracktype="bigWig",
        short_label=f"{srr_id} overlay",
        long_label=f"{rich} (strand overlay)",
        aggregate="transparentOverlay",
        showSubtrackColorOnUi="on",
        visibility="full",
        viewLimits=VIEW_LIMITS,
        autoScale="group",
        maxHeightPixels="100:50:8",
    )

    agg.add_subtrack(trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg_forward"),
        url=f"{ctx.base_url}/{fwd}",
        tracktype="bigWig",
        short_label="forward",
        long_label=f"{srr_id} forward",
        color=ctx.colors["fwd"],
    ))

    agg.add_subtrack(trackhub.Track(
        name=trackhub.helpers.sanitize(f"{srr_id}_agg_reverse"),
        url=f"{ctx.base_url}/{rev}",
        tracktype="bigWig",
        short_label="reverse",
        long_label=f"{srr_id} reverse",
        color=ctx.colors["rev"],
        negateValues="on",
    ))
    return agg


# ----- region (bigBed) construction ----------------------------------------
#
# Region tracks live in their OWN SuperTrack, which is a SIBLING of the
# signal composite. This keeps signal and region tracks cleanly separated
# in the UCSC browser, users can toggle the whole region group on/off
# independently.


def build_region_track(
    srr_id: str,
    label: str,
    rel_path: str,
    ctx: BuildContext,
    metadata_row: dict | None = None,
) -> trackhub.Track:
    """Build one bigBed region subtrack.

    bigBed tracks are fundamentally different from bigWig signal tracks:
    no viewLimits, no autoScale, no negateValues. They display genomic
    features (ORFs) as discrete blocks.
    """
    name = trackhub.helpers.sanitize(f"{srr_id}_region_{label}")
    rich = long_label_for_sample(srr_id, metadata_row, ctx.label_fields)
    region_color = hex_to_trackhub_rgb(DEFAULT_REGION_COLOR)

    return trackhub.Track(
        name=name,
        url=f"{ctx.base_url}/{rel_path}",
        tracktype="bigBed",
        short_label=f"{srr_id} {label}",
        long_label=f"{rich} ({label})",
        color=region_color,
        visibility="dense",
    )


def build_region_container(
    regions: dict[str, list[tuple[str, str]]],
    ctx: BuildContext,
    metadata: dict[str, dict[str, str]],
    container_name: str = "ribohub_regions",
) -> trackhub.SuperTrack | None:
    """Build a SuperTrack holding all bigBed region subtracks.

    Returns None if no region files were discovered, the caller should
    simply skip emitting the container in that case.
    """
    if not regions:
        return None

def build_aggregate_container(
    aggregates: list[trackhub.AggregateTrack],
    container_name:str = "ribohub_aggregates",
) -> trackhub.SuperTrack | None:
    if not aggregates:
        return None

    super_track = trackhub.SuperTrack(
        name=trackhub.helpers.sanitize(container_name),
        short_label="Strand overlays",
        long_label="Per-sample strand-overlay aggregates",
    )

    for agg in aggregates:
        super_track.add_tracks(agg)
    
    log.info(f"Built aggregate container with {len(aggregates)} track(s).")
    return super_track
    
    super_track = trackhub.SuperTrack(
        name=trackhub.helpers.sanitize(container_name),
        short_label="Region annotations",
        long_label="Per-sample bigBed region annotations (ORFs)",
    )

    n_tracks = 0
    for srr_id in sorted(regions.keys()):
        for label, rel_path in sorted(regions[srr_id]):
            track = build_region_track(
                srr_id, label, rel_path, ctx, metadata.get(srr_id),
            )
            super_track.add_tracks(track)
            n_tracks += 1

    log.info(f"Built region container with {n_tracks} track(s) across "
             f"{len(regions)} sample(s).")
    return super_track


# ----- description HTML ----------------------------------------------------
#
# A static description.html lives in a fixed location next to ribohub.py.
# It is copied into the output directory at build time so it travels with
# the hub. To update the page content, edit that file directly, no code
# change needed.
#
# Convention: <ribohub_dir>/static/description.html

STATIC_DIR: Path = Path(__file__).parent / "static"


def copy_description_html(output_dir: Path) -> Path | None:
    """Copy static/description.html into output_dir if it exists.

    Returns the destination path, or None if the source file is absent
    (in which case the hub is built without a description page).
    """
    src = STATIC_DIR / "description.html"
    if not src.exists():
        log.warning(
            f"No description HTML found at {src}. "
            "Hub will be built without a description page. "
            "Create static/description.html next to ribohub.py to add one."
        )
        return None
    dst = output_dir / "description.html"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)
    return dst


# ----- output writers ------------------------------------------------------
# trackhub renders most fields natively via str(). The one exception is the
# per-subtrack `metadata` line (rejected by trackhub 1.0 as unknown), so
# writers render via str(comp) then inject metadata lines afterward.


def _inject_metadata(rendered: str,
                     metadata_map: dict[str, dict[str, str]],
                     indent: str = "    ") -> str:
    """Insert `metadata ...` lines into a rendered composite.

    The metadata line is placed immediately after each subtrack's `subGroups`
    line (or, if a subtrack has no subGroups, after its `track` line). The
    map is keyed by track name, matched against the `track <name>` line.
    """
    if not metadata_map:
        return rendered  # nothing to inject; hand back the text unchanged

    out: list[str] = []          # the new text, line by line, as we build it
    current: str | None = None   # name of the track we're currently inside
    pending: dict[str, str] | None = None  # metadata waiting to be written out

    # Small helper: turn a metadata dict into the UCSC text line, indented.
    def fmt(md: dict[str, str]) -> str:
        pairs = " ".join(f'{k}="{v}"' for k, v in md.items())
        return f"{indent}metadata {pairs}"

    # Strategy: scan line by line. On a `track NAME` line, check if that
    # track has metadata (store in `pending`). Write it after the track's
    # `subGroups` line, or before the next `track` line if there's no subGroups.
    for line in rendered.splitlines():
        stripped = line.strip()
        if stripped.startswith("track "):
            # Flush any pending metadata from the previous track before starting the new one.
            if pending:
                out.append(fmt(pending))
                pending = None
            current = stripped.split(None, 1)[1]   # the text after "track "
            out.append(line)
            pending = metadata_map.get(current)    # metadata for this track?
            continue
        out.append(line)
        # Found the anchor: write the pending metadata right after subGroups.
        if pending and stripped.startswith("subGroups "):
            out.append(fmt(pending))
            pending = None
    # Handle the very last track (no following `track` line to trigger a flush).
    if pending:
        out.append(fmt(pending))
    return "\n".join(out)


def _render_composite(comp: trackhub.CompositeTrack,
                      metadata_map: dict[str, dict[str, str]]) -> str:
    """Render the composite via trackhub, then inject metadata lines."""
    return _inject_metadata(str(comp), metadata_map)


def _render_aggregate(agg: trackhub.AggregateTrack) -> str:
    """Render an aggregate (no subGroups/metadata extras needed)."""
    return str(agg)


def _render_region_container(container: trackhub.SuperTrack) -> str:
    """Render the region SuperTrack."""
    return str(container)


def write_single_file_hub(
    hub_name: str,
    short_label: str,
    long_label: str,
    email: str,
    genome: str,
    composite: trackhub.CompositeTrack,
    aggregates: list[trackhub.AggregateTrack],
    output_dir: Path,
    metadata_map: dict[str, dict[str, str]] | None = None,
    region_container: trackhub.SuperTrack | None = None,
) -> Path:
    """Write a useOneFile hub with composite + aggregates + regions as siblings."""
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{hub_name}.hub.txt"

    lines: list[str] = [
        f"hub {hub_name}",
        f"shortLabel {short_label}",
        f"longLabel {long_label}",
        "useOneFile on",
        f"email {email}",
        "",
        f"genome {genome}",
        "",
        _render_composite(composite, metadata_map or {}),
        "",
    ]

    agg_container = build_aggregate_container(aggregates)
    if agg_container is not None:
        lines.append(_render_region_container(agg_container))
        lines.append("")

    if region_container is not None:
        lines.append(_render_region_container(region_container))
        lines.append("")

    out_path.write_text("\n".join(lines))
    return out_path


def write_directory_hub(
    hub_name: str,
    short_label: str,
    long_label: str,
    email: str,
    genome: str,
    composite: trackhub.CompositeTrack,
    aggregates: list[trackhub.AggregateTrack],
    output_dir: Path,
    metadata_map: dict[str, dict[str, str]] | None = None,
    region_container: trackhub.SuperTrack | None = None,
) -> Path:
    """Write a multi-file hub: hub.txt + genomes.txt + <genome>/trackDb.txt.

    Uses str(comp) + metadata injection instead of trackhub.upload.stage_hub(),
    because bigWigs are hosted at an external base_url, not staged locally.
    """
    staging = output_dir / hub_name
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    genome_dir = staging / genome
    genome_dir.mkdir()

    (staging / f"{hub_name}.hub.txt").write_text("\n".join([
        f"hub {hub_name}",
        f"shortLabel {short_label}",
        f"longLabel {long_label}",
        f"genomesFile genomes.txt",
        f"email {email}",
        "",
    ]))
    (staging / "genomes.txt").write_text("\n".join([
        f"genome {genome}",
        f"trackDb {genome}/trackDb.txt",
        "",
    ]))
    trackdb_lines: list[str] = [
        _render_composite(composite, metadata_map or {}),
        "",
    ]
    agg_container = build_aggregate_container(aggregates)
    if agg_container is not None:
        lines.append(_render_region_container(agg_container))
        lines.append("")

    if region_container is not None:
        trackdb_lines.append(_render_region_container(region_container))
        trackdb_lines.append("")

    (genome_dir / "trackDb.txt").write_text("\n".join(trackdb_lines))

    return staging / f"{hub_name}.hub.txt"


# ----- CLI -----------------------------------------------------------------

@click.group()
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
@click.pass_context
def cli(ctx: click.Context, verbose: bool) -> None:
    """RiboHub: Build a Gwips RiboSeq track hub from organized bigWig files."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


@cli.command()
# Required inputs
@click.option("--samples", default=None,
              help="Sample selection. Single ID, comma list, .txt, or .csv (first column). "
                   "Optional when --filter is given (filter derives the sample set from metadata).")
@click.option("--data-dir", required=True, envvar="RIBOHUB_DATA_DIR",
              type=click.Path(exists=True, file_okay=False, path_type=Path), #TODO 
              help="Root directory containing sorted bigWig files.")
@click.option("--output-dir", required=True, envvar="RIBOHUB_OUTPUT_DIR", #TODO
              type=click.Path(file_okay=False, path_type=Path),
              help="Directory where the hub will be written.")
@click.option("--base-url", required=True, envvar="RIBOHUB_BASE_URL", #TODO
              help="Public URL where bigWig files are served (must be reachable by Gwips).")
# Metadata
@click.option("--metadata", "metadata_path", default=None,
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="Optional RiboSeqOrg-style metadata CSV. "
                   "Auto-discovered from --data-dir if named metadata.csv or RiboSeqOrg_Metadata.csv.")
@click.option("--filter", "metadata_filter", default=None,
              help="Filter samples from metadata by column values. "
                   "Format: COL=VAL or COL=VAL1|VAL2, comma-separated for AND logic. "
                   "Example: --filter \"CONDITION=High|Test,CELL_LINE=HEK293\". "
                   f"Suggested columns: {', '.join(SUGGESTED_FILTER_COLUMNS)}. "
                   "Requires --metadata or a metadata file in --data-dir.")
# Hub configuration
@click.option("--genome", default="hg38", show_default=True,
              help="Gwips genome assembly.")
@click.option("--hub-name", default="RiboSeqHub", show_default=True,
              help="Hub directory name (no spaces).")
@click.option("--email", default="your@email.com", show_default=True,
              help="Contact email written into hub.txt.")
# Output
@click.option("--output-format", default="directory", show_default=True,
              type=click.Choice(["directory", "single-file"]),
              help="directory: multi-file hub; single-file: useOneFile hub.")
# Build behavior
@click.option("--kinds", default="all,unique,multi", show_default=True,
              help="Which kinds to include in the composite. Subset of: all, unique, multi.")
@click.option("--with-aggregates/--no-aggregates", default=True, show_default=True,
              help="Whether to include per-sample strand-overlay aggregates.")
@click.option("--with-regions/--no-regions", default=True, show_default=True,
              help="Whether to include per-sample bigBed region tracks.")
@click.option("--strict", is_flag=True,
              help="Exit non-zero if any requested sample is missing or partial.")
@click.option("--dry-run", is_flag=True,
              help="Report what would be built without writing the hub.")
@click.option("--auto-scale/--no-auto-scale", default=True, show_default=True,
              help="Let UCSC auto-scale the y-axis to the visible signel. "
                   "When off, uses fixed viewLimits (-127:127).")
# Colors
@click.option("--color-fwd", default=DEFAULT_COLORS["fwd"], show_default=True,
              help="Forward strand color (hex).")
@click.option("--color-rev", default=DEFAULT_COLORS["rev"], show_default=True,
              help="Reverse strand color (hex).")
@click.option("--color-fwd-multi", default=DEFAULT_COLORS["fwd_multi"], show_default=True,
              help="Forward multimapped color (hex).")
@click.option("--color-rev-multi", default=DEFAULT_COLORS["rev_multi"], show_default=True,
              help="Reverse multimapped color (hex).")
def generate(
    samples: str | None,
    data_dir: Path,
    output_dir: Path,
    base_url: str,
    metadata_path: Path | None,
    metadata_filter: str | None,
    genome: str,
    hub_name: str,
    email: str,
    output_format: str,
    kinds: str,
    with_aggregates: bool,
    with_regions: bool,
    strict: bool,
    dry_run: bool,
    auto_scale : bool,
    color_fwd: str,
    color_rev: str,
    color_fwd_multi: str,
    color_rev_multi: str,
) -> None:
    """Generate a Gwips track hub for the given samples."""
    # ---- This function is the conductor: it runs every stage in order. ----
    #  1. tidy + validate the user's options
    #  2. find the sample files on disk      (discover_samples)
    #  3. read the metadata spreadsheet       (load_metadata)
    #  4. early exits: nothing found / dry-run / strict-mode refusals
    #  5. build the composite + aggregates    (build_composite_trackdb, build_aggregate)
    #  6. build the region container          (build_region_container)
    #  7. write the hub to disk               (write_*_hub)
    #  8. report what happened
    base_url = base_url.rstrip("/")  # drop any trailing "/" so we don't get "//"

    colors: dict[str, str] = {
        "fwd":       hex_to_trackhub_rgb(color_fwd),
        "rev":       hex_to_trackhub_rgb(color_rev),
        "fwd_multi": hex_to_trackhub_rgb(color_fwd_multi),
        "rev_multi": hex_to_trackhub_rgb(color_rev_multi),
    }

    parsed_kinds_raw = parse_csv_field(kinds)
    valid_kinds = {"all", "unique", "multi"}
    bad_kinds = [k for k in parsed_kinds_raw if k not in valid_kinds]
    if bad_kinds:
        raise click.BadParameter(
            f"--kinds: unknown values {bad_kinds}. Must be subset of {sorted(valid_kinds)}."
        )
    parsed_kinds = frozenset(parsed_kinds_raw)
    if not parsed_kinds:
        raise click.BadParameter("--kinds cannot be empty.")

    ctx = BuildContext(
        base_url=base_url,
        colors=colors,
        label_fields=LABEL_FIELDS,
        subgroup_fields=SUBGROUP_FIELDS,
        metadata_fields=LABEL_FIELDS,  # same fields drive both label and metadata lines
        kinds=parsed_kinds,
        with_aggregates=with_aggregates,
        with_regions=with_regions,
        auto_scale=auto_scale,
    )

    # ---- stage 0: validate --samples / --filter mutual requirements ----
    # At least one of --samples or --filter must be given.
    if not samples and not metadata_filter:
        click.echo(
            "ERROR: provide --samples, --filter, or both.",
            err=True,
        )
        sys.exit(1)

    # ---- stage 1: metadata autodiscovery ----
    # Auto-discover metadata CSV in --data-dir if not explicitly provided
    if metadata_path is None:
        for candidate in ("metadata.csv", "RiboSeqOrg_Metadata.csv"):
            candidate_path = data_dir / candidate
            if candidate_path.exists():
                log.info(f"Auto-discovered metadata: {candidate_path}")
                metadata_path = candidate_path
                break

    # --filter requires metadata either explicit or auto-discovered.
    if metadata_filter and metadata_path is None:
        click.echo(
            "ERROR: --filter requires metadata.\n"
            "       Pass --metadata explicitly, or place metadata.csv / "
            "RiboSeqOrg_Metadata.csv in --data-dir.",
            err=True,
        )
        sys.exit(1)

    # ---- stage 2: load metadata (needed for --filter; optional otherwise) ----
    metadata: dict[str, dict[str, str]] = {}
    csv_columns: list[str] = []
    if metadata_path:
        required = tuple(sorted(set(LABEL_FIELDS + SUBGROUP_FIELDS)))
        metadata = load_metadata(metadata_path, SRR_COLUMN, required)
        # Re-open just to grab the header for unknown-column checks.
        with metadata_path.open(newline="", encoding="utf-8") as _f:
            csv_columns = list(csv.DictReader(_f).fieldnames or [])

    # ---- stage 3: resolve sample_set ----
    # Three possible inputs, in priority order:
    #   a) --samples only          → use it directly
    #   b) --filter only           → derive from metadata
    #   c) both                    → intersect (with warning if any dropped)
    sample_set: set[str] = set()

    if samples:
        sample_set = parse_samples(samples)
        log.info(f"Resolved {len(sample_set)} requested samples from --samples.")

    if metadata_filter:
        parsed_filter = parse_filter(metadata_filter)
        result = apply_metadata_filter(metadata, parsed_filter, csv_columns)

        # ---- error: unknown columns ----
        if result.unknown_columns:
            click.echo(
                f"ERROR: --filter references unknown column(s): "
                f"{', '.join(result.unknown_columns)}\n"
                f"       Suggested filter columns: "
                f"{', '.join(SUGGESTED_FILTER_COLUMNS)}\n"
                f"       All available columns: "
                f"{', '.join(result.available_columns)}",
                err=True,
            )
            sys.exit(1)

        # ---- error: zero-match fields ----
        if result.zero_fields:
            lines = ["ERROR: --filter matched 0 samples for the following field(s):\n"]
            for col in result.zero_fields:
                avail = result.available_values.get(col, [])
                avail_str = (
                    ", ".join(avail[:20])
                    + (" ..." if len(avail) > 20 else "")
                ) if avail else "(no non-null values in metadata)"
                lines.append(f"  {col}: no match")
                lines.append(f"       Available values: {avail_str}")
            click.echo("\n".join(lines), err=True)
            sys.exit(1)

        # ---- per-field match summary (informational) ----
        for col, count in result.per_field.items():
            log.info(f"  --filter {col}: {count} match(es) in metadata.")

        if not sample_set:
            # --filter only: use matched set directly
            sample_set = result.matched
            log.info(f"--filter resolved {len(sample_set)} sample(s) from metadata.")
        else:
            # both --samples and --filter: intersect
            before = len(sample_set)
            sample_set = sample_set & result.matched
            dropped = before - len(sample_set)
            if dropped:
                click.echo(
                    f"WARNING: --samples and --filter both given; intersecting.\n"
                    f"         --samples: {before} IDs\n"
                    f"         --filter matched: {len(result.matched)} IDs\n"
                    f"         Intersection: {len(sample_set)} IDs "
                    f"({dropped} dropped).",
                )
            else:
                log.info(
                    f"--samples and --filter intersection: "
                    f"{len(sample_set)} IDs (no samples dropped)."
                )

        if not sample_set:
            click.echo(
                "ERROR: sample set is empty after applying --filter.\n"
                "       Check that the SRR IDs in --samples appear in the metadata.",
                err=True,
            )
            sys.exit(1)

    found, regions, dir_missing, dir_empty = discover_samples(
        sample_set, data_dir, with_regions=with_regions,
    )

    # ---- warn if filter matched metadata but nothing landed on disk ----
    if metadata_filter and not found:
        click.echo(
            f"ERROR: --filter matched {len(sample_set)} sample(s) in metadata "
            f"but 0 have bigWig files in --data-dir.\n"
            f"       Matched: {', '.join(sorted(sample_set))}\n"
            f"       Check that --data-dir points to the right location.",
            err=True,
        )
        sys.exit(1)

    report_discovery(found, regions, dir_missing, dir_empty, sample_set)

    # Validate that the data has at least one of the requested kinds.
    # Error out cleanly if not, better than building an empty composite.
    actual_kinds = {_KIND_TO_TAG[fk.kind]
                    for files in found.values() for fk in files}
    requested_present = parsed_kinds & actual_kinds
    if found and not requested_present:
        click.echo(
            f"ERROR: --kinds requested {sorted(parsed_kinds)} but no matching "
            f"files found in any of {len(found)} samples.",
            err=True,
        )
        click.echo(f"       Available kinds in data: {sorted(actual_kinds)}", err=True)
        sys.exit(1)

    if metadata:
        missing_meta = sorted(set(found) - set(metadata))
        if missing_meta:
            log.warning(
                f"{len(missing_meta)} discovered sample(s) have no metadata row. "
                f"Examples: {', '.join(missing_meta[:5])}"
            )

    if not found:
        click.echo("ERROR: No samples found. Nothing to build.", err=True)
        sys.exit(1)

    if dry_run:
        click.echo(
            f"Dry run: would build a composite of {len(found)} sample(s) "
            f"with kinds={sorted(parsed_kinds)}, aggregates={with_aggregates}, "
            f"regions={with_regions} ({len(regions)} with region files), "
            f"format={output_format}."
        )
        if metadata_filter:
            click.echo(f"  Filter: {metadata_filter} → {len(sample_set)} matched in metadata, "
                       f"{len(found)} have files on disk.")
        if strict and (dir_missing or dir_empty):
            click.echo(f"  Note: --strict would reject this run "
                       f"({len(dir_missing)} missing, {len(dir_empty)} empty).",
                       err=True)
        if metadata:
            with_meta = sum(1 for s in found if s in metadata)
            click.echo(f"Metadata matched: {with_meta}/{len(found)} sample(s).")
        return

    if strict and (dir_missing or dir_empty):
        click.echo("ERROR (--strict): refusing to build a partial hub.", err=True)
        sys.exit(1)

    # ----- build composite -----
    composite, metadata_map = build_composite_trackdb(found, metadata, ctx)
    if not composite.subtracks:
        click.echo(
            "ERROR: no subtracks built (data has no files matching --kinds).",
            err=True,
        )
        sys.exit(1)

    # ----- build aggregates -----
    aggregates: list[trackhub.AggregateTrack] = []
    skipped_aggregates: list[str] = []
    if with_aggregates:
        for srr_id in sorted(found.keys()):
            agg = build_aggregate(srr_id, found[srr_id], ctx, metadata.get(srr_id))
            if agg:
                aggregates.append(agg)
            else:
                skipped_aggregates.append(srr_id)
                log.warning(f"{srr_id}: lacks bare files; aggregate skipped.")

        if skipped_aggregates and len(skipped_aggregates) == len(found):
            click.echo(
                "ERROR: --with-aggregates requested but no samples have both "
                "bare _forward and _reverse files.",
                err=True,
            )
            sys.exit(1)

        if strict and skipped_aggregates:
            click.echo(
                f"ERROR (--strict): {len(skipped_aggregates)} sample(s) had "
                f"aggregate skipped. Refusing to ship a partial hub.",
                err=True,
            )
            sys.exit(1)

    # ----- build region container -----
    region_container: trackhub.SuperTrack | None = None
    if with_regions and regions:
        region_container = build_region_container(regions, ctx, metadata)

    # ----- copy description HTML -----
    desc_path = copy_description_html(output_dir)
    if desc_path:
        log.info(f"Copied description HTML: {desc_path}")

    # ----- write output -----
    if output_format == "directory":
        hub_file = write_directory_hub(
            hub_name=hub_name,
            short_label=hub_name,
            long_label=f"{hub_name} Track Hub",
            email=email,
            genome=genome,
            composite=composite,
            aggregates=aggregates,
            output_dir=output_dir,
            metadata_map=metadata_map,
            region_container=region_container,
        )
        public_url = f"{base_url}/{hub_name}/{hub_name}.hub.txt"
    else:
        hub_file = write_single_file_hub(
            hub_name=hub_name,
            short_label=hub_name,
            long_label=f"{hub_name} Track Hub",
            email=email,
            genome=genome,
            composite=composite,
            aggregates=aggregates,
            output_dir=output_dir,
            metadata_map=metadata_map,
            region_container=region_container,
        )
        public_url = f"{base_url}/{hub_name}.hub.txt"

    click.echo(
        f"Built hub with {len(composite.subtracks)} composite subtrack(s), "
        f"{len(aggregates)} aggregate(s)"
        f"{f' and {len(region_container.subtracks)} region track(s)' if region_container else ''}"
        f" across {len(found)} sample(s)."
    )
    click.echo(f"Wrote: {hub_file}")
    click.echo(f"Hub URL: {public_url}")

if __name__ == "__main__":
    cli()
