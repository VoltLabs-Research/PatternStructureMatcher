# PatternStructureMatching

`PatternStructureMatching` classifies atoms by matching them against dynamic lattice YAML definitions.

Within VOLT, the plugin also lets the user define a lattice once and derives both:

- the PatternStructureMatching pattern YAML
- the OpenDXA reference lattice YAML

That makes it usable as the upstream structure-identification algorithm for the
`OpenDXA` plugin without forcing the user to maintain two parallel lattice files.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- PatternStructureMatching
```

## CLI

Usage:

```bash
pattern-structure-matching <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--lattice_dir <path>` | Yes | Directory containing PatternStructureMatching lattice YAMLs. | |
| `--reference_lattice_dir <path>` | Yes | Directory containing OpenDXA reference lattice YAMLs. | |
| `--patterns <csv>` | No | Optional lattice filter, for example `fcc,bcc`. | all lattices |
| `--dissolve_small_clusters` | No | Mark small clusters as `OTHER` after clustering. | `false` |
| `--help` | No | Print CLI help. | |

## VOLT Plugin Workflow

Inside VOLT, `PatternStructureMatching` is exposed as a structure-identification
plugin that can be selected upstream of `OpenDXA`.

The plugin does not expose presets in the VOLT UI. Instead, the user defines one
or more custom patterns directly in the form and the wrapper materializes both:

- a PatternStructureMatching lattice YAML for local matching
- an OpenDXA reference lattice YAML with the derived `neighbor_vectors`

The generated artifacts are written to:

- `<output_base>_generated_lattices/pattern-structure-matching`
- `<output_base>_generated_lattices/opendxa`
- `<output_base>_pattern_structure_matching_manifest.json`

### Pattern Definition UX

Each custom pattern includes:

- `Pattern Name`
- `Use As Matrix Phase`
- `Coordination Number`
- `Scale`
- `Coordinate Mode`
- `Reference Basis Atom Index`
- `Cell A`, `Cell B`, `Cell C` as vector tuples with `x`, `y`, `z`
- `Basis Atoms` as a nested editable list of `species`, `x`, `y`, `z`

This makes the UI behave like a guided crystal-template editor instead of a
single free-form YAML or DSL text box.

### Matrix Phase Resolution

The plugin also exposes `Matrix Pattern Name Override`.

Resolution order is:

1. `Matrix Pattern Name Override`, if set
2. the single pattern marked as `Use As Matrix Phase`
3. the only defined pattern, when there is exactly one

If more than one pattern is defined and none is unambiguously selected as the
matrix phase, execution fails with a validation error.

### Using It With OpenDXA

When `PatternStructureMatching` is selected as the upstream algorithm inside the
`OpenDXA` plugin:

- `OpenDXA` passes its matrix phase through `reference_topology`
- that value is mapped into the PSM `Matrix Pattern Name Override`
- the PSM wrapper writes matching OpenDXA YAMLs under `<output_base>_generated_lattices/opendxa`
- the `OpenDXA` wrapper consumes that generated directory automatically

If the matrix phase itself is one of the user-defined custom patterns, the name
seen by OpenDXA must match the corresponding `Pattern Name`.

### Conceptual Example

The UI is equivalent to defining data shaped like this:

```text
Pattern 1:
  name = bi_a7
  is_matrix_phase = true
  coordination_number = 6
  coordinate_mode = fractional
  cell_a = (1, 0, 0)
  cell_b = (0.5, 0.866, 0)
  cell_c = (0, 0, 1.5)
  basis_atoms = [
    { species, x, y, z },
    { species, x, y, z },
    ...
  ]
```
