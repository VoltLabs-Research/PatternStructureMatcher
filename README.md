# PatternStructureMatching

`PatternStructureMatching` classifies atoms by matching them against dynamic lattice YAML definitions.

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
| `--lattice-dir <path>` | Yes | Directory containing PatternStructureMatching lattice YAMLs. | |
| `--reference-lattice-dir <path>` | Yes | Directory containing OpenDXA reference lattice YAMLs. | |
| `--patterns <csv>` | No | Optional lattice filter, for example `fcc,bcc`. | all lattices |
| `--dissolveSmallClusters` | No | Mark small clusters as `OTHER` after clustering. | `false` |
| `--help` | No | Print CLI help. | |

## Build With CoreToolkit

```bash
cd /path/to/voltlabs-ecosystem/tools/CoreToolkit
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/StructureIdentification
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/CommonNeighborAnalysis
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/PatternStructureMatching
conan create . -nr
```
