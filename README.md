# PatternStructureMatching

`PatternStructureMatching` classifies atoms by matching them against dynamic lattice YAML definitions.

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
| `--lattice-dir <path>` | Yes | Directory containing PatternStructureMatching lattice YAMLs. | |
| `--reference-lattice-dir <path>` | Yes | Directory containing OpenDXA reference lattice YAMLs. | |
| `--patterns <csv>` | No | Optional lattice filter, for example `fcc,bcc`. | all lattices |
| `--dissolveSmallClusters` | No | Mark small clusters as `OTHER` after clustering. | `false` |
| `--help` | No | Print CLI help. | |
