#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
PLUGIN_ROOT = Path(os.environ.get("PLUGIN_PROJECT_DIR", SCRIPT_DIR.parent)).resolve()
PLUGINS_ROOT = PLUGIN_ROOT.parent if PLUGIN_ROOT.parent.name == "plugins" else None
EMBEDDED_LOADER = (PLUGIN_ROOT / "lib" / "ld-linux-x86-64.so.2").resolve()
EMBEDDED_LIBRARY_DIR = (PLUGIN_ROOT / "lib").resolve()
ENV_PSM_BINARY_OVERRIDES = ("VOLT_PATTERN_STRUCTURE_MATCHING_BINARY", "VOLT_PSM_BINARY")

RESERVED_RUNTIME_FLAGS_WITH_VALUE = {
    "--selectedTimesteps",
    "--selected-timesteps",
}

MAX_IMAGE_RADIUS = 4


@dataclass(frozen=True)
class BasisAtom:
    species: int
    position: tuple[float, float, float]


@dataclass(frozen=True)
class PatternSpec:
    name: str
    coordination_number: int
    scale: float
    cell: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    coordinate_mode: str
    basis: tuple[BasisAtom, ...]


@dataclass(frozen=True)
class PatternDefinition:
    spec: PatternSpec
    reference_basis_index: int
    is_matrix_phase: bool


class WrapperError(RuntimeError):
    pass


def log(message: str) -> None:
    sys.stderr.write(f"[pattern-structure-matching-plugin] {message}\n")
    sys.stderr.flush()


def format_number(value: float) -> str:
    if abs(value) < 1e-12:
        value = 0.0
    text = f"{value:.12f}".rstrip("0").rstrip(".")
    if "." not in text:
        text += ".0"
    return text


def resolve_existing_path(label: str, candidates: Iterable[Path | None]) -> Path:
    for candidate in candidates:
        if candidate and candidate.exists():
            return candidate.resolve()
    pretty = "\n".join(f"  - {candidate}" for candidate in candidates if candidate)
    raise WrapperError(f"No pude resolver {label}. Paths probados:\n{pretty}")


def repo_binary_candidates(binary_name: str, repo_names: Iterable[str]) -> list[Path]:
    candidates: list[Path] = []
    repo_names = list(repo_names)

    if PLUGINS_ROOT is not None:
        for repo_name in repo_names:
            candidates.extend([
                PLUGINS_ROOT / repo_name / "build" / "build" / "Release" / binary_name,
                PLUGINS_ROOT / repo_name / "build-local" / "build" / "Release" / binary_name,
                PLUGINS_ROOT / repo_name / "build-manual" / "build" / "Release" / binary_name,
                PLUGINS_ROOT / repo_name / "build" / "build" / "Release" / "build" / "Release" / binary_name,
                PLUGINS_ROOT / repo_name / "build-local" / "build" / "Release" / "build" / "Release" / binary_name,
            ])

    return candidates


def resolve_binary() -> Path:
    candidates: list[Path | None] = []
    for env_var in ENV_PSM_BINARY_OVERRIDES:
        value = os.environ.get(env_var, "").strip()
        if value:
            candidates.append(Path(value))
    candidates.append(PLUGIN_ROOT / "bin" / "pattern-structure-matching")
    candidates.extend(repo_binary_candidates("pattern-structure-matching", ("PatternStructureMatching", "PatternStructureMatcher")))
    which = shutil.which("pattern-structure-matching")
    if which:
        candidates.append(Path(which))
    return resolve_existing_path("binario pattern-structure-matching", candidates)


def resolve_embedded_runtime_command(command: list[str]) -> list[str]:
    if not command:
        return command
    if not EMBEDDED_LOADER.exists() or not EMBEDDED_LIBRARY_DIR.exists():
        return command

    binary_path = Path(command[0]).resolve()
    try:
        binary_path.relative_to(PLUGIN_ROOT)
    except ValueError:
        return command

    if binary_path == EMBEDDED_LOADER:
        return command

    return [
        str(EMBEDDED_LOADER),
        "--library-path",
        str(EMBEDDED_LIBRARY_DIR),
        str(binary_path),
        *command[1:],
    ]


def ensure_executable(path: Path) -> None:
    try:
        mode = path.stat().st_mode
        if not (mode & 0o111):
            path.chmod(mode | 0o755)
    except OSError:
        pass


def run(command: list[str]) -> None:
    command = resolve_embedded_runtime_command(command)
    log(" ".join(command))
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise WrapperError(f"El comando fallo con exit code {completed.returncode}: {command[0]}")


def ensure_parent_directory(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def require_file(path: Path) -> Path:
    if not path.exists():
        raise WrapperError(f"Falta el archivo requerido: {path}")
    return path


def normalize_token(value: str) -> str:
    return re.sub(r"[\s\-]+", "_", value.strip()).lower()


def sanitize_lattice_name(value: str) -> str:
    token = value.strip().replace("/", "_").replace("\\", "_")
    token = re.sub(r"\s+", "_", token)
    if not token:
        raise WrapperError("El nombre del patron no puede ser vacio.")
    return token


def write_text_file(path: Path, content: str) -> None:
    ensure_parent_directory(path)
    path.write_text(content, encoding="utf-8")


def require_non_empty_string(value: Any, label: str) -> str:
    if not isinstance(value, str):
        raise WrapperError(f"{label} debe ser un string.")
    text = value.strip()
    if not text:
        raise WrapperError(f"{label} no puede ser vacio.")
    return text


def coerce_float(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise WrapperError(f"{label} debe ser numerico.")

    if isinstance(value, (int, float)):
        result = float(value)
    elif isinstance(value, str) and value.strip():
        try:
            result = float(value.strip())
        except ValueError as error:
            raise WrapperError(f"{label} debe ser numerico.") from error
    else:
        raise WrapperError(f"{label} debe ser numerico.")

    if not math.isfinite(result):
        raise WrapperError(f"{label} debe ser finito.")
    return result


def coerce_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise WrapperError(f"{label} debe ser entero.")

    if isinstance(value, int):
        result = value
    elif isinstance(value, float) and value.is_integer():
        result = int(value)
    elif isinstance(value, str) and value.strip():
        try:
            parsed = float(value.strip())
        except ValueError as error:
            raise WrapperError(f"{label} debe ser entero.") from error
        if not parsed.is_integer():
            raise WrapperError(f"{label} debe ser entero.")
        result = int(parsed)
    else:
        raise WrapperError(f"{label} debe ser entero.")

    return result


def coerce_bool(value: Any, label: str) -> bool:
    if isinstance(value, bool):
        return value

    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "y"}:
            return True
        if normalized in {"false", "0", "no", "n", ""}:
            return False

    raise WrapperError(f"{label} debe ser booleano.")


def parse_json_argument(raw_value: str, label: str) -> Any:
    if not raw_value.strip():
        raise WrapperError(f"{label} no puede estar vacio.")

    try:
        return json.loads(raw_value)
    except json.JSONDecodeError as error:
        raise WrapperError(f"{label} no es JSON valido: {error.msg}.") from error


def parse_basis_atoms(raw_value: Any, label: str) -> tuple[BasisAtom, ...]:
    if not isinstance(raw_value, list) or not raw_value:
        raise WrapperError(f"{label} debe contener al menos un atomo.")

    atoms: list[BasisAtom] = []
    for index, raw_atom in enumerate(raw_value):
        atom_label = f"{label}[{index}]"
        if not isinstance(raw_atom, dict):
            raise WrapperError(f"{atom_label} debe ser un objeto.")

        species = coerce_int(raw_atom.get("species", 1), f"{atom_label}.species")
        if species <= 0:
            raise WrapperError(f"{atom_label}.species debe ser mayor que cero.")

        atoms.append(BasisAtom(
            species=species,
            position=(
                coerce_float(raw_atom.get("x"), f"{atom_label}.x"),
                coerce_float(raw_atom.get("y"), f"{atom_label}.y"),
                coerce_float(raw_atom.get("z"), f"{atom_label}.z"),
            ),
        ))

    return tuple(atoms)


def parse_vector3(raw_value: Any, label: str, default: tuple[float, float, float]) -> tuple[float, float, float]:
    if raw_value is None:
        raw_value = {}

    if isinstance(raw_value, dict):
        return (
            coerce_float(raw_value.get("x", default[0]), f"{label}.x"),
            coerce_float(raw_value.get("y", default[1]), f"{label}.y"),
            coerce_float(raw_value.get("z", default[2]), f"{label}.z"),
        )

    if isinstance(raw_value, (list, tuple)) and len(raw_value) == 3:
        return (
            coerce_float(raw_value[0], f"{label}[0]"),
            coerce_float(raw_value[1], f"{label}[1]"),
            coerce_float(raw_value[2], f"{label}[2]"),
        )

    raise WrapperError(f"{label} debe ser un vector de tres componentes.")


def parse_cell_vector(
    raw_pattern: dict[str, Any],
    pattern_label: str,
    key: str,
    legacy_prefix: str,
    default: tuple[float, float, float],
) -> tuple[float, float, float]:
    tuple_value = raw_pattern.get(key)
    if tuple_value is not None:
        return parse_vector3(tuple_value, f"{pattern_label}.{key}", default)

    legacy_value = {
        "x": raw_pattern.get(f"{legacy_prefix}x", default[0]),
        "y": raw_pattern.get(f"{legacy_prefix}y", default[1]),
        "z": raw_pattern.get(f"{legacy_prefix}z", default[2]),
    }
    return parse_vector3(legacy_value, f"{pattern_label}.{key}", default)


def parse_pattern_definitions(raw_value: str) -> list[PatternDefinition]:
    payload = parse_json_argument(raw_value, "pattern_definitions")
    if not isinstance(payload, list) or not payload:
        raise WrapperError("pattern_definitions debe contener al menos un patron.")

    definitions: list[PatternDefinition] = []
    seen_names: set[str] = set()

    for index, raw_pattern in enumerate(payload):
        pattern_label = f"pattern_definitions[{index}]"
        if not isinstance(raw_pattern, dict):
            raise WrapperError(f"{pattern_label} debe ser un objeto.")

        name = sanitize_lattice_name(require_non_empty_string(raw_pattern.get("name"), f"{pattern_label}.name"))
        normalized_name = normalize_token(name)
        if normalized_name in seen_names:
            raise WrapperError(f"El patron '{name}' esta repetido.")
        seen_names.add(normalized_name)

        coordination_number = coerce_int(
            raw_pattern.get("coordination_number"),
            f"{pattern_label}.coordination_number",
        )
        if coordination_number <= 0:
            raise WrapperError(f"{pattern_label}.coordination_number debe ser mayor que cero.")

        coordinate_mode = normalize_token(str(raw_pattern.get("coordinate_mode", "fractional")))
        if coordinate_mode not in {"fractional", "cartesian"}:
            raise WrapperError(f"{pattern_label}.coordinate_mode debe ser fractional o cartesian.")

        reference_basis_index = coerce_int(
            raw_pattern.get("reference_basis_index", 0),
            f"{pattern_label}.reference_basis_index",
        )
        if reference_basis_index < 0:
            raise WrapperError(f"{pattern_label}.reference_basis_index no puede ser negativo.")

        spec = PatternSpec(
            name=name,
            coordination_number=coordination_number,
            scale=coerce_float(raw_pattern.get("scale", 1.0), f"{pattern_label}.scale"),
            cell=(
                parse_cell_vector(raw_pattern, pattern_label, "cell_a", "cell_a", (1.0, 0.0, 0.0)),
                parse_cell_vector(raw_pattern, pattern_label, "cell_b", "cell_b", (0.0, 1.0, 0.0)),
                parse_cell_vector(raw_pattern, pattern_label, "cell_c", "cell_c", (0.0, 0.0, 1.0)),
            ),
            coordinate_mode=coordinate_mode,
            basis=parse_basis_atoms(raw_pattern.get("basis_atoms"), f"{pattern_label}.basis_atoms"),
        )

        definitions.append(PatternDefinition(
            spec=spec,
            reference_basis_index=reference_basis_index,
            is_matrix_phase=coerce_bool(raw_pattern.get("is_matrix_phase", False), f"{pattern_label}.is_matrix_phase"),
        ))

    return definitions


def resolve_matrix_pattern_name(
    matrix_pattern_name: str,
    pattern_definitions: list[PatternDefinition],
) -> str:
    definitions_by_name = {
        normalize_token(pattern_definition.spec.name): pattern_definition.spec.name
        for pattern_definition in pattern_definitions
    }

    override = matrix_pattern_name.strip()
    if override:
        normalized_override = normalize_token(sanitize_lattice_name(override))
        resolved = definitions_by_name.get(normalized_override)
        if not resolved:
            raise WrapperError(
                "matrix_pattern_name debe coincidir con uno de los patrones definidos por el usuario."
            )
        return resolved

    matrix_phase_patterns = [
        pattern_definition.spec.name
        for pattern_definition in pattern_definitions
        if pattern_definition.is_matrix_phase
    ]
    if len(matrix_phase_patterns) == 1:
        return matrix_phase_patterns[0]
    if len(matrix_phase_patterns) > 1:
        raise WrapperError(
            "Solo un patron puede estar marcado como matrix phase. "
            "Usa una sola marca o define matrix_pattern_name."
        )
    if len(pattern_definitions) == 1:
        return pattern_definitions[0].spec.name

    raise WrapperError(
        "Define matrix_pattern_name o marca exactamente un patron como matrix phase."
    )


def build_pattern_yaml(spec: PatternSpec) -> str:
    lines = [
        f"name: {spec.name}",
        f"coordination_number: {spec.coordination_number}",
        "",
    ]
    if abs(spec.scale - 1.0) > 1e-12:
        lines.append(f"scale: {format_number(spec.scale)}")
        lines.append("")

    lines.append("cell:")
    for vector in spec.cell:
        lines.append(f"  - [{format_number(vector[0])}, {format_number(vector[1])}, {format_number(vector[2])}]")

    lines.extend([
        "",
        f"coordinate_mode: {spec.coordinate_mode}",
        "",
        "basis:",
    ])
    for atom in spec.basis:
        lines.append(f"  - species: {atom.species}")
        lines.append(
            "    position: "
            f"[{format_number(atom.position[0])}, {format_number(atom.position[1])}, {format_number(atom.position[2])}]"
        )
    lines.append("")
    return "\n".join(lines)


def build_opendxa_yaml(name: str, coordination_number: int, vectors: list[tuple[float, float, float]]) -> str:
    lines = [
        f"name: {name}",
        f"coordination_number: {coordination_number}",
        "",
        "neighbor_vectors:",
    ]
    for vector in vectors:
        lines.append(f"  - [{format_number(vector[0])}, {format_number(vector[1])}, {format_number(vector[2])}]")
    lines.append("")
    return "\n".join(lines)


def apply_scale(vector: tuple[float, float, float], scale: float) -> tuple[float, float, float]:
    return vector[0] * scale, vector[1] * scale, vector[2] * scale


def cell_basis_to_cartesian(
    cell: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    position: tuple[float, float, float],
) -> tuple[float, float, float]:
    ax, ay, az = cell[0]
    bx, by, bz = cell[1]
    cx, cy, cz = cell[2]
    px, py, pz = position
    return (
        px * ax + py * bx + pz * cx,
        px * ay + py * by + pz * cy,
        px * az + py * bz + pz * cz,
    )


def vector_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


def vector_length_sq(vector: tuple[float, float, float]) -> float:
    return vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2


def derive_neighbor_vectors(spec: PatternSpec, reference_basis_index: int) -> list[tuple[float, float, float]]:
    if reference_basis_index < 0 or reference_basis_index >= len(spec.basis):
        raise WrapperError(
            f"reference_basis_index={reference_basis_index} esta fuera de rango para '{spec.name}'. "
            f"La basis tiene {len(spec.basis)} atomos."
        )

    if spec.coordinate_mode == "fractional":
        basis_cartesian = [apply_scale(cell_basis_to_cartesian(spec.cell, atom.position), spec.scale) for atom in spec.basis]
    else:
        basis_cartesian = [apply_scale(atom.position, spec.scale) for atom in spec.basis]

    reference_position = basis_cartesian[reference_basis_index]
    neighbor_candidates: list[tuple[float, tuple[float, float, float]]] = []
    seen: set[tuple[int, int, int]] = set()

    for image_x in range(-MAX_IMAGE_RADIUS, MAX_IMAGE_RADIUS + 1):
        for image_y in range(-MAX_IMAGE_RADIUS, MAX_IMAGE_RADIUS + 1):
            for image_z in range(-MAX_IMAGE_RADIUS, MAX_IMAGE_RADIUS + 1):
                image_shift = apply_scale(
                    cell_basis_to_cartesian(spec.cell, (float(image_x), float(image_y), float(image_z))),
                    spec.scale,
                )
                for atom_index, atom_position in enumerate(basis_cartesian):
                    if atom_index == reference_basis_index and image_x == 0 and image_y == 0 and image_z == 0:
                        continue
                    candidate = (
                        atom_position[0] + image_shift[0],
                        atom_position[1] + image_shift[1],
                        atom_position[2] + image_shift[2],
                    )
                    vector = vector_sub(candidate, reference_position)
                    if vector_length_sq(vector) < 1e-16:
                        continue
                    key = (
                        int(round(vector[0] * 1e8)),
                        int(round(vector[1] * 1e8)),
                        int(round(vector[2] * 1e8)),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    neighbor_candidates.append((vector_length_sq(vector), vector))

    neighbor_candidates.sort(key=lambda item: (round(item[0], 12), item[1]))
    if len(neighbor_candidates) < spec.coordination_number:
        raise WrapperError(
            f"No pude derivar {spec.coordination_number} neighbor_vectors para '{spec.name}'. "
            f"Solo encontre {len(neighbor_candidates)}."
        )
    return [vector for _, vector in neighbor_candidates[:spec.coordination_number]]


def materialize_lattice_directories(
    args: argparse.Namespace,
    output_base: str,
) -> tuple[Path, Path, dict[str, object], list[str], str]:
    pattern_definitions = parse_pattern_definitions(args.pattern_definitions)

    generated_root = Path(f"{output_base}_generated_lattices")
    pattern_dir = generated_root / "pattern-structure-matching"
    opendxa_dir = generated_root / "opendxa"

    if generated_root.exists():
        shutil.rmtree(generated_root)

    pattern_dir.mkdir(parents=True, exist_ok=True)
    opendxa_dir.mkdir(parents=True, exist_ok=True)

    selected_patterns = [pattern_definition.spec.name for pattern_definition in pattern_definitions]
    reference_topology = resolve_matrix_pattern_name(
        args.matrix_pattern_name or args.reference_topology,
        pattern_definitions,
    )

    manifest: dict[str, object] = {
        "generated_root": str(generated_root),
        "pattern_lattice_dir": str(pattern_dir),
        "opendxa_lattice_dir": str(opendxa_dir),
        "selected_patterns": selected_patterns,
        "generated_pattern_lattices": [],
        "generated_opendxa_lattices": [],
        "reference_topology": reference_topology,
        "pattern_summaries": [],
    }

    for pattern_definition in pattern_definitions:
        pattern_path = pattern_dir / f"{pattern_definition.spec.name}.yml"
        opendxa_path = opendxa_dir / f"{pattern_definition.spec.name}.yml"

        write_text_file(pattern_path, build_pattern_yaml(pattern_definition.spec))
        write_text_file(
            opendxa_path,
            build_opendxa_yaml(
                pattern_definition.spec.name,
                pattern_definition.spec.coordination_number,
                derive_neighbor_vectors(pattern_definition.spec, pattern_definition.reference_basis_index),
            ),
        )

        manifest["generated_pattern_lattices"].append(str(pattern_path))
        manifest["generated_opendxa_lattices"].append(str(opendxa_path))
        manifest["pattern_summaries"].append({
            "name": pattern_definition.spec.name,
            "coordination_number": pattern_definition.spec.coordination_number,
            "coordinate_mode": pattern_definition.spec.coordinate_mode,
            "basis_atoms": len(pattern_definition.spec.basis),
            "reference_basis_index": pattern_definition.reference_basis_index,
            "is_matrix_phase": pattern_definition.is_matrix_phase,
        })

    manifest_path = Path(f"{output_base}_pattern_structure_matching_manifest.json")
    write_text_file(manifest_path, json.dumps(manifest, indent=2, sort_keys=True))
    manifest["manifest_path"] = str(manifest_path)

    return pattern_dir, opendxa_dir, manifest, selected_patterns, reference_topology


def ensure_outputs(output_base: str) -> dict[str, str]:
    return {
        "annotatedDump": str(require_file(Path(f"{output_base}_annotated.dump"))),
        "clustersTable": str(require_file(Path(f"{output_base}_clusters.table"))),
        "clustersTransitions": str(require_file(Path(f"{output_base}_cluster_transitions.table"))),
        "patternAnalysis": str(require_file(Path(f"{output_base}_pattern_analysis.msgpack"))),
    }


def build_psm_command(
    psm_binary: Path,
    args: argparse.Namespace,
    input_dump: str,
    output_base: str,
    pattern_dir: Path,
    opendxa_dir: Path,
    selected_patterns: list[str],
) -> list[str]:
    command = [
        str(psm_binary),
        input_dump,
        output_base,
        "--lattice_dir",
        str(pattern_dir),
        "--reference_lattice_dir",
        str(opendxa_dir),
        "--patterns",
        ",".join(selected_patterns),
    ]
    if args.dissolve_small_clusters:
        command.append("--dissolve_small_clusters")
    return command


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="pattern_structure_matching_plugin_wrapper.py",
        description="VOLT wrapper for PatternStructureMatching with user-defined lattice generation.",
    )
    parser.add_argument("input_dump")
    parser.add_argument("output_base")
    parser.add_argument("--pattern_definitions", default="")
    parser.add_argument("--matrix_pattern_name", default="")
    parser.add_argument("--reference_topology", default="")
    parser.add_argument("--dissolve_small_clusters", action="store_true")
    parser.add_argument("--only_generate_yaml", action="store_true")

    args, unknown = parser.parse_known_args(argv)
    if unknown:
        raise WrapperError(f"Argumentos no soportados por PatternStructureMatching: {unknown}")
    return args


def filter_reserved_args(args: list[str]) -> list[str]:
    filtered: list[str] = []
    index = 0
    while index < len(args):
        token = str(args[index])
        if token in RESERVED_RUNTIME_FLAGS_WITH_VALUE:
            index += 2
            continue
        filtered.append(token)
        index += 1
    return filtered


def run_pipeline(args: argparse.Namespace) -> dict[str, object]:
    output_base = args.output_base
    ensure_parent_directory(Path(output_base))

    pattern_dir, opendxa_dir, manifest, selected_patterns, reference_topology = materialize_lattice_directories(
        args,
        output_base,
    )

    result: dict[str, object] = {
        "ok": True,
        "outputBase": output_base,
        "patternLatticeDir": str(pattern_dir),
        "opendxaLatticeDir": str(opendxa_dir),
        "selectedPatterns": selected_patterns,
        "reference_topology": reference_topology,
        "manifestPath": manifest["manifest_path"],
    }

    if args.only_generate_yaml:
        result["generatedOnly"] = True
        return result

    psm_binary = resolve_binary()
    ensure_executable(psm_binary)
    if EMBEDDED_LOADER.exists():
        ensure_executable(EMBEDDED_LOADER)

    run(build_psm_command(psm_binary, args, args.input_dump, output_base, pattern_dir, opendxa_dir, selected_patterns))
    result["artifacts"] = ensure_outputs(output_base)
    result["generatedOnly"] = False
    return result


def main(argv: list[str] | None = None) -> int:
    args = parse_args(filter_reserved_args(list(sys.argv[1:] if argv is None else argv)))
    run_pipeline(args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except WrapperError as error:
        log(f"error: {error}")
        raise SystemExit(1)
