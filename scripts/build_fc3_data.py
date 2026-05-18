from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FC3_PATH = ROOT / "FORCE_CONSTANTS_3RD_0.45"
QE_PATH = ROOT / "BASE.si_supper.scf.in"
OUT_PATH = ROOT / "app" / "data" / "fc3_045.json"


@dataclass
class Fc3Block:
    block_id: int
    rj: list[float]
    rk: list[float]
    atoms: list[int]
    values: list[float]

    @property
    def norm(self) -> float:
        return math.sqrt(sum(value * value for value in self.values))


def dot_vec(coeff: list[float], cell: list[list[float]]) -> list[float]:
    return [
        coeff[0] * cell[0][axis] + coeff[1] * cell[1][axis] + coeff[2] * cell[2][axis]
        for axis in range(3)
    ]


def sub_vec(a: list[float], b: list[float]) -> list[float]:
    return [a[i] - b[i] for i in range(3)]


def dist(a: list[float], b: list[float]) -> float:
    d = sub_vec(a, b)
    return math.sqrt(sum(x * x for x in d))


def triangle_perimeter(a: list[float], b: list[float], c: list[float]) -> float:
    return dist(a, b) + dist(a, c) + dist(b, c)


def solve_3x3_transposed(cell: list[list[float]], vec: list[float]) -> list[float]:
    # Solve cell.T * coeff = vec without adding a dependency on numpy.
    a = [[cell[col][row] for col in range(3)] for row in range(3)]
    b = vec[:]
    for pivot in range(3):
        best = max(range(pivot, 3), key=lambda row: abs(a[row][pivot]))
        if abs(a[best][pivot]) < 1e-12:
            raise ValueError("Singular cell matrix")
        if best != pivot:
            a[pivot], a[best] = a[best], a[pivot]
            b[pivot], b[best] = b[best], b[pivot]
        scale = a[pivot][pivot]
        a[pivot] = [x / scale for x in a[pivot]]
        b[pivot] /= scale
        for row in range(3):
            if row == pivot:
                continue
            factor = a[row][pivot]
            a[row] = [a[row][col] - factor * a[pivot][col] for col in range(3)]
            b[row] -= factor * b[pivot]
    return b


def read_qe_supercell(path: Path) -> tuple[list[list[float]], list[dict]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    positions: list[dict] = []
    cell: list[list[float]] = []

    for idx, line in enumerate(lines):
        head = line.strip().lower()
        if head.startswith("atomic_positions"):
            j = idx + 1
            while j < len(lines):
                parts = lines[j].split()
                if len(parts) < 4:
                    break
                try:
                    frac = [float(parts[1]), float(parts[2]), float(parts[3])]
                except ValueError:
                    break
                positions.append({"species": parts[0], "frac": frac})
                j += 1
        if head.startswith("cell_parameters"):
            cell = [[float(x) for x in lines[idx + n].split()[:3]] for n in range(1, 4)]

    if len(cell) != 3 or not positions:
        raise ValueError(f"Failed to parse QE input: {path}")
    return cell, positions


def infer_primitive_from_supercell(
    supercell: list[list[float]], positions: list[dict]
) -> tuple[list[list[float]], list[dict], int]:
    nat = len(positions)
    if nat % 2 != 0:
        raise ValueError("This Si workflow expects a two-atom primitive basis")
    repeat = round((nat / 2) ** (1 / 3))
    if 2 * repeat**3 != nat:
        raise ValueError(f"Cannot infer cubic primitive repeat from nat={nat}")

    primitive = [[value / repeat for value in row] for row in supercell]
    basis: list[dict] = []
    seen: set[tuple[int, int, int]] = set()
    for pos in positions:
        prim_frac_raw = [x * repeat for x in pos["frac"]]
        prim_frac = [x - math.floor(x + 1e-9) for x in prim_frac_raw]
        key = tuple(round(x * 1_000_000) for x in prim_frac)
        if key in seen:
            continue
        seen.add(key)
        basis.append(
            {
                "id": len(basis) + 1,
                "species": pos["species"],
                "frac": prim_frac,
                "cart": dot_vec(prim_frac, primitive),
            }
        )
        if len(basis) == 2:
            break

    if len(basis) != 2:
        raise ValueError("Failed to infer two primitive basis atoms")
    return primitive, basis, repeat


def read_fc3(path: Path) -> list[Fc3Block]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    nblocks = int(lines[0])
    blocks: list[Fc3Block] = []
    idx = 1
    for _ in range(nblocks):
        block_id = int(lines[idx])
        rj = [float(x) for x in lines[idx + 1].split()[:3]]
        rk = [float(x) for x in lines[idx + 2].split()[:3]]
        atoms = [int(x) for x in lines[idx + 3].split()[:3]]
        values = [float(lines[idx + 4 + n].split()[3]) for n in range(27)]
        blocks.append(Fc3Block(block_id, rj, rk, atoms, values))
        idx += 31
    if len(blocks) != nblocks:
        raise ValueError(f"Expected {nblocks} FC3 blocks, parsed {len(blocks)}")
    return blocks


def lattice_index(shift_cart: list[float], primitive: list[list[float]]) -> list[int]:
    coeff = solve_3x3_transposed(primitive, shift_cart)
    rounded = [round(x) for x in coeff]
    if any(abs(coeff[i] - rounded[i]) > 2e-5 for i in range(3)):
        raise ValueError(f"Shift is not an integer primitive lattice vector: {shift_cart} -> {coeff}")
    return [int(x) for x in rounded]


def shell_key(value: float, width: float = 0.25) -> float:
    return round(round(value / width) * width, 3)


def build_payload() -> dict:
    supercell, positions = read_qe_supercell(QE_PATH)
    primitive, basis, repeat = infer_primitive_from_supercell(supercell, positions)
    blocks = read_fc3(FC3_PATH)

    payload_blocks = []
    shell_map: dict[float, dict[str, float]] = {}
    for block in blocks:
        atoms = [basis[idx - 1] for idx in block.atoms]
        rj_idx = lattice_index(block.rj, primitive)
        rk_idx = lattice_index(block.rk, primitive)
        xi = atoms[0]["cart"]
        xj = [atoms[1]["cart"][axis] + dot_vec(rj_idx, primitive)[axis] for axis in range(3)]
        xk = [atoms[2]["cart"][axis] + dot_vec(rk_idx, primitive)[axis] for axis in range(3)]
        dmax = max(dist(xi, xj), dist(xi, xk), dist(xj, xk))
        perimeter = triangle_perimeter(xi, xj, xk)
        norm = block.norm
        key = shell_key(dmax)
        shell = shell_map.setdefault(key, {"dmax": key, "strength": 0.0, "count": 0})
        shell["strength"] += norm
        shell["count"] += 1
        payload_blocks.append(
            {
                "id": block.block_id,
                "atoms": block.atoms,
                "rj": block.rj,
                "rk": block.rk,
                "rjIndex": rj_idx,
                "rkIndex": rk_idx,
                "vertices": [xi, xj, xk],
                "values": block.values,
                "norm": norm,
                "dmax": dmax,
                "perimeter": perimeter,
            }
        )

    norms = [block["norm"] for block in payload_blocks]
    dmax_values = [block["dmax"] for block in payload_blocks]
    perimeter_values = [block["perimeter"] for block in payload_blocks]
    component_values = [value for block in payload_blocks for value in block["values"]]
    component_abs_nonzero = [abs(value) for value in component_values if value != 0]
    component_zero_count = len(component_values) - len(component_abs_nonzero)
    payload_blocks.sort(key=lambda item: item["id"])
    top_blocks = sorted(payload_blocks, key=lambda item: item["norm"], reverse=True)[:40]
    shells = sorted(shell_map.values(), key=lambda item: item["dmax"])

    return {
        "meta": {
            "fc3File": FC3_PATH.name,
            "qeFile": QE_PATH.name,
            "cutoffNm": 0.45,
            "repeat": repeat,
            "blockCount": len(payload_blocks),
            "normMin": min(norms),
            "normMax": max(norms),
            "dmaxMin": min(dmax_values),
            "dmaxMax": max(dmax_values),
            "perimeterMin": min(perimeter_values),
            "perimeterMax": max(perimeter_values),
            "componentCount": len(component_values),
            "componentZeroCount": component_zero_count,
            "componentAbsMin": min(component_abs_nonzero) if component_abs_nonzero else 0.0,
            "componentAbsMax": max(component_abs_nonzero) if component_abs_nonzero else 0.0,
        },
        "cell": primitive,
        "basis": basis,
        "blocks": payload_blocks,
        "shells": shells,
        "topBlocks": top_blocks,
    }


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    OUT_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {OUT_PATH} with {payload['meta']['blockCount']} blocks")


if __name__ == "__main__":
    main()
