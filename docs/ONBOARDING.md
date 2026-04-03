# PhyLabeler Onboarding Guide

## What This Project Does

PhyLabeler checks whether clades in a phylogenetic tree are monophyletic with
respect to the NCBI taxonomy database.

At a high level:

1. Load taxonomy data from `names.dmp` and `nodes.dmp`, or from a cached pickle.
2. Parse a Newick tree into an in-memory node structure.
3. Resolve every tip label in the tree to an NCBI taxonomic identifier.
4. For each internal non-root clade, compute the taxonomic MRCA shared by all
   tips in that clade.
5. Test whether that MRCA is exclusive to the clade or also appears in taxa
   outside the clade.
6. Mark the clade as monophyletic, non-monophyletic, or unresolved.
7. Write a labeled Newick tree and a text report.

## Fast Reading Order

Use this order if you want to understand the modern codebase quickly:

1. `main.py`
   Learn how the project decides between CLI and GUI and where the main
   workflows start.
2. `taxonomy_db.py`
   Understand how NCBI taxonomy is loaded, parsed, cached, and queried.
3. `tree_parser.py`
   Understand the internal tree model and how bipartitions are derived.
4. `monophyly.py`
   Understand the analysis algorithm and the result objects.
5. `gui.py`
   See how the analysis engine is exposed through the Tkinter interface.

Only after that, look at the legacy files:

- `LabelPhy.py`
- `TreeCode.py`
- `Renamer.py`
- `NameNavigate.py`

These older modules are useful for project history, but they are not the main
implementation you should study first.

## Core Modules

### `main.py`

This is the project dispatcher.

- Parses command-line arguments.
- Chooses GUI or CLI mode.
- Loads taxonomy.
- Parses one tree or many trees.
- Runs the monophyly checker.
- Writes labeled tree and report outputs.

Important functions:

- `main()`
- `run_cli(args)`
- `_process_single(args, checker)`
- `_process_batch(args, checker)`
- `run_gui()`

## `taxonomy_db.py`

This file builds the project knowledge base from NCBI taxonomy.

### Main responsibilities

- Parse `names.dmp`
- Parse `nodes.dmp`
- Store lookup dictionaries
- Save and load a pickle cache
- Download fresh taxonomy from NCBI
- Provide lineage lookups used by the monophyly engine

### Key data structures

- `name_to_code`
  Maps taxon name to NCBI tax ID.
- `code_to_names`
  Maps tax ID to all known names.
- `code_to_scientific`
  Maps tax ID to the scientific name.
- `code_to_parent`
  Maps tax ID to parent tax ID.
- `code_to_rank`
  Maps tax ID to taxonomic rank.

### Tip-name resolution behavior

The lookup strategy is intentionally simple:

1. Exact match
2. Replace underscores with spaces
3. Case-insensitive scan

That means naming quality in the input tree matters a lot.

## `tree_parser.py`

This file defines the in-memory tree representation and the Newick parser.

### `Node`

Each node stores:

- `label`
- `branch_length`
- `parent`
- `children`
- `is_tip`
- `monophyly_status`
- `taxonomy_code`
- `taxonomy_name`

### Important tree methods

- `get_tip_labels()`
  Collect all descendant tip labels.
- `get_tip_nodes()`
  Collect all descendant tip node objects.
- `get_bipartitions()`
  Return a post-order list of internal non-root clades.
- `to_newick()`
  Rebuild the tree as a Newick string.
- `count_tips()`, `count_internal()`, `depth()`
  Provide structural summaries.

### Why bipartitions matter here

The monophyly checker analyzes each internal clade independently. In practice,
the clade is represented by the set of tip labels below a node.

## `monophyly.py`

This is the project's analysis engine.

### `MonophylyResult`

Each analyzed clade gets a result object containing:

- the node reference
- the tip labels in the clade
- the taxonomic MRCA code
- the MRCA scientific name
- the MRCA rank
- the status
- the candidate names for that MRCA
- the intruder taxa outside the clade

### Step-by-step algorithm

1. Collect all tree tips.
2. Resolve each tip to an NCBI tax ID.
3. Convert each resolved tip to a lineage of tax IDs from tip to root.
4. Get every internal non-root clade from the tree in post-order.
5. For one clade at a time:
   - gather the tip lineages inside the clade
   - find the first tax ID shared by all of them
   - treat that as the clade's taxonomic MRCA
   - inspect the complement of the clade
   - if any outside lineage also contains the MRCA, record those tips as
     intruders
6. Assign a status:
   - `monophyletic`
   - `not_monophyletic`
   - `unresolved`
7. If monophyletic, label the internal node with the cleaned MRCA name.

### Important conceptual point

Monophyly here is implemented as an exclusivity test relative to the other tips
present in the input tree, not as a separate phylogenetic inference procedure.

## `gui.py`

This file wraps the core engine in a Tkinter desktop interface.

### Main UI sections

- Taxonomy loading controls
- Tree file list
- Tree visualization canvas
- Analysis buttons
- Results notebook

### GUI workflow

1. Start the app.
2. Load taxonomy from cache, local files, or NCBI download.
3. Load one tree or many trees.
4. Parse the selected tree.
5. Run the monophyly analysis.
6. View the results:
   - Summary
   - Non-monophyletic clades
   - Monophyletic clades
   - Unresolved tips
7. Export a labeled tree or report.

### Tree drawing model

The GUI uses a rectangular phylogram-like layout:

- y-position comes from tip order
- x-position comes from node depth
- branches are colored by monophyly status
- internal labels are drawn near labeled clades

## CLI Workflow

### Single-tree mode

Typical path:

1. load taxonomy
2. parse one tree
3. check monophyly across all internal clades
4. label monophyletic internal nodes
5. write:
   - `*.labeled.nwk`
   - `*.report.txt`

### Batch mode

Typical path:

1. scan a directory for Newick-like files
2. parse each tree
3. run the same checker on each tree
4. write one labeled output per tree
5. write a `batch_summary.txt`

## Output Files

### Labeled Newick

The labeled tree preserves the original topology and branch lengths, but
monophyletic internal nodes receive cleaned taxonomy labels suitable for Newick.

### Report

The report summarizes:

- total internal nodes analyzed
- monophyletic count
- non-monophyletic count
- unresolved count
- a sorted list of monophyletic clades
- a sorted list of non-monophyletic clades with intruders
- unresolved tips, if any

## Recommended Study Path For PI Meetings

Use this sequence when explaining the project:

1. Problem statement
   This tool checks whether gene-tree clades agree with NCBI taxonomy.
2. Inputs
   Newick tree plus NCBI taxonomy files or cache.
3. Internal model
   Parsed tree plus taxonomy lookup tables.
4. Core method
   MRCA sharing inside the clade, then exclusivity against the complement.
5. Outputs
   Labeled tree and report.
6. Limitations
   Depends heavily on tip-name matching and the taxa present in the tree.

## Known Caveats

- Tip matching is mostly literal.
- Unresolved names are not force-resolved with fuzzy methods.
- Only internal non-root nodes are analyzed as bipartitions.
- The modern UI and CLI share the same core analysis engine.
- The legacy files reflect an older implementation style and should not be
  treated as the primary architecture.

## Practical Questions To Ask While Reading

- How are tip labels normalized before lookup?
- What exactly counts as a bipartition in this implementation?
- How is the MRCA computed from lineage lists?
- Why does a clade fail monophyly even if its internal members share a common
  taxonomy?
- Which parts are UI-only and which parts are analysis-core?
- What assumptions are made about the input tree labels?

## Suggested Next Improvements

If you want to extend the project later, likely high-value directions are:

- stronger tip-name normalization and alias handling
- search/filter tools in the results UI
- richer tree rendering for large trees
- explicit export of unresolved-name diagnostics
- clearer separation of legacy and modern code paths
