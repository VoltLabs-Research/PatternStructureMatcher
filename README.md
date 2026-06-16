# Pattern Structure Matching

Matches atomic structures against user-defined patterns and exports the reconstructed-state contract consumed by OpenDXA.

## Install

```bash
vpm install @voltlabs/pattern-structure-matching
```

## CLI

```bash
pattern-structure-matching <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--pattern_definitions <list>` | yes | — | Custom lattice patterns (name, cell vectors, basis atoms, coordination number, scale, coordinate mode). |
| `--matrix_pattern_name <name>` | no | `""` | Override which defined pattern is the matrix phase. |
| `--dissolve_small_clusters` | no | `false` | Mark small clusters as `OTHER` after clustering. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_atoms.parquet` | Structure Identification | AtomisticExporter → glb |
| `{output_base}_pattern_analysis.parquet` | Pattern Analysis | — |
| `{output_base}_neighbor_lattice.parquet` | Neighbor Lattice | — |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/pattern-structure-matcher
