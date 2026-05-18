const state = {
  data: null,
  defaultData: null,
  structure: null,
  rotation: { x: -0.55, y: 0.72, z: 0 },
  zoomScale: 1,
  dragging: false,
  lastPointer: null,
  pointerDown: null,
  selectedBlockId: null,
  pickedAtom: null,
  hiddenBlockIds: new Set(),
  hiddenAtomKeys: new Set(),
  hiddenElements: new Set(),
  collapsedTree: new Set(),
  objectTreeInitialized: false,
  scenePickAtoms: [],
  normRange: { min: 0, max: 0 },
  compareStructure: null,
  compareSlots: [],
};

const COMPONENT_PLOT_THRESHOLD = 1e-4;
const BLOCK_PLOT_THRESHOLD = 1e-4;
const SYMMETRY_RTOL = 1e-6;
const SYMMETRY_ATOL = 1e-8;
const COMPARE_COLORS = ["#005eb8", "#c62828", "#2e7d32", "#6a1b9a"];

const $ = (id) => document.getElementById(id);

function add(a, b) {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

function sub(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

function scale(v, s) {
  return [v[0] * s, v[1] * s, v[2] * s];
}

function length(v) {
  return Math.hypot(v[0], v[1], v[2]);
}

function trianglePerimeter(a, b, c) {
  return length(sub(a, b)) + length(sub(a, c)) + length(sub(b, c));
}

function dotCell(index, cell) {
  return [
    index[0] * cell[0][0] + index[1] * cell[1][0] + index[2] * cell[2][0],
    index[0] * cell[0][1] + index[1] * cell[1][1] + index[2] * cell[2][1],
    index[0] * cell[0][2] + index[1] * cell[1][2] + index[2] * cell[2][2],
  ];
}

function solveCellTranspose(cell, vec) {
  const a = [
    [cell[0][0], cell[1][0], cell[2][0]],
    [cell[0][1], cell[1][1], cell[2][1]],
    [cell[0][2], cell[1][2], cell[2][2]],
  ];
  const b = [...vec];
  for (let pivot = 0; pivot < 3; pivot += 1) {
    let best = pivot;
    for (let row = pivot + 1; row < 3; row += 1) {
      if (Math.abs(a[row][pivot]) > Math.abs(a[best][pivot])) best = row;
    }
    if (Math.abs(a[best][pivot]) < 1e-12) throw new Error("Singular cell matrix");
    if (best !== pivot) {
      [a[pivot], a[best]] = [a[best], a[pivot]];
      [b[pivot], b[best]] = [b[best], b[pivot]];
    }
    const denom = a[pivot][pivot];
    a[pivot] = a[pivot].map((value) => value / denom);
    b[pivot] /= denom;
    for (let row = 0; row < 3; row += 1) {
      if (row === pivot) continue;
      const factor = a[row][pivot];
      a[row] = a[row].map((value, col) => value - factor * a[pivot][col]);
      b[row] -= factor * b[pivot];
    }
  }
  return b;
}

function latticeIndex(shiftCart, cell) {
  const coeff = solveCellTranspose(cell, shiftCart);
  const rounded = coeff.map((value) => Math.round(value));
  if (coeff.some((value, axis) => Math.abs(value - rounded[axis]) > 2e-5)) {
    throw new Error(`FC3 translation is not a primitive lattice vector: ${shiftCart.join(", ")}`);
  }
  return rounded;
}

function parseQeStructure(text, name = "uploaded QE structure") {
  const lines = text.split(/\r?\n/);
  const positions = [];
  let cell = null;

  for (let idx = 0; idx < lines.length; idx += 1) {
    const head = lines[idx].trim().toLowerCase();
    if (head.startsWith("atomic_positions")) {
      let row = idx + 1;
      while (row < lines.length) {
        const parts = lines[row].trim().split(/\s+/);
        if (parts.length < 4) break;
        const frac = parts.slice(1, 4).map(Number);
        if (frac.some((value) => Number.isNaN(value))) break;
        positions.push({ species: parts[0], frac });
        row += 1;
      }
    }
    if (head.startsWith("cell_parameters")) {
      cell = [1, 2, 3].map((offset) => lines[idx + offset].trim().split(/\s+/).slice(0, 3).map(Number));
    }
  }

  if (!cell || positions.length === 0) {
    throw new Error(`Failed to parse ${name}. Expected CELL_PARAMETERS and ATOMIC_POSITIONS crystal.`);
  }
  const primitive = inferPrimitive(cell, positions);
  return { ...primitive, source: name };
}

function inferPrimitive(supercell, positions) {
  const nat = positions.length;
  if (nat <= 16) {
    return {
      cell: supercell,
      basis: positions.map((pos, idx) => ({
        id: idx + 1,
        species: normalizeSpecies(pos.species),
        frac: pos.frac.map((value) => value - Math.floor(value)),
        cart: dotCell(pos.frac, supercell),
      })),
      repeat: 1,
    };
  }

  for (let repeat = 8; repeat >= 2; repeat -= 1) {
    const primitiveCount = nat / repeat ** 3;
    if (!Number.isInteger(primitiveCount) || primitiveCount < 1) continue;
    const cell = supercell.map((row) => row.map((value) => value / repeat));
    const basis = [];
    const seen = new Set();
    for (const pos of positions) {
      const primFrac = pos.frac.map((value) => {
        const raw = value * repeat;
        return raw - Math.floor(raw + 1e-9);
      });
      const species = normalizeSpecies(pos.species);
      const key = `${species}:${primFrac.map((value) => Math.round(value * 1_000_000)).join(",")}`;
      if (seen.has(key)) continue;
      seen.add(key);
      basis.push({
        id: basis.length + 1,
        species,
        frac: primFrac,
        cart: dotCell(primFrac, cell),
      });
    }
    if (basis.length === primitiveCount) return { cell, basis, repeat };
  }
  throw new Error(`Cannot infer a primitive cell from nat=${nat}. Upload a primitive QE structure if this is not a cubic repeated supercell.`);
}

function normalizeSpecies(species) {
  const text = String(species || "X").trim();
  return text.charAt(0).toUpperCase() + text.slice(1).toLowerCase();
}

function parseFc3(text, name = "uploaded FORCE_CONSTANTS_3RD") {
  const lines = text.split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
  const nblocks = Number.parseInt(lines[0], 10);
  if (!Number.isFinite(nblocks) || nblocks <= 0) throw new Error(`Invalid FC3 block count in ${name}.`);
  const blocks = [];
  let idx = 1;
  for (let blockIndex = 0; blockIndex < nblocks; blockIndex += 1) {
    const id = Number.parseInt(lines[idx], 10);
    const rj = lines[idx + 1].split(/\s+/).slice(0, 3).map(Number);
    const rk = lines[idx + 2].split(/\s+/).slice(0, 3).map(Number);
    const atoms = lines[idx + 3].split(/\s+/).slice(0, 3).map((value) => Number.parseInt(value, 10));
    const values = [];
    for (let n = 0; n < 27; n += 1) {
      const parts = lines[idx + 4 + n].split(/\s+/);
      values.push(Number(parts[3]));
    }
    if ([id, ...rj, ...rk, ...atoms, ...values].some((value) => !Number.isFinite(value))) {
      throw new Error(`Failed to parse FC3 block near line ${idx + 1} in ${name}.`);
    }
    blocks.push({ id, rj, rk, atoms, values });
    idx += 31;
  }
  if (blocks.length !== nblocks) throw new Error(`Expected ${nblocks} FC3 blocks, parsed ${blocks.length}.`);
  return blocks;
}

function shellKey(value, width = 0.25) {
  return Math.round(Math.round(value / width) * width * 1000) / 1000;
}

function buildPayloadFromFc3Blocks(rawBlocks, structure, fc3File) {
  const shellMap = new Map();
  const blocks = rawBlocks.map((block) => {
    const atomRecords = block.atoms.map((atomId) => {
      const atom = structure.basis[atomId - 1];
      if (!atom) throw new Error(`FC3 block ${block.id} references atom ${atomId}, but structure has ${structure.basis.length} basis atoms.`);
      return atom;
    });
    const rjIndex = latticeIndex(block.rj, structure.cell);
    const rkIndex = latticeIndex(block.rk, structure.cell);
    const xi = atomRecords[0].cart;
    const rjShift = dotCell(rjIndex, structure.cell);
    const rkShift = dotCell(rkIndex, structure.cell);
    const xj = add(atomRecords[1].cart, rjShift);
    const xk = add(atomRecords[2].cart, rkShift);
    const norm = Math.sqrt(block.values.reduce((total, value) => total + value * value, 0));
    const dmax = Math.max(length(sub(xi, xj)), length(sub(xi, xk)), length(sub(xj, xk)));
    const perimeter = trianglePerimeter(xi, xj, xk);
    const key = shellKey(dmax);
    const shell = shellMap.get(key) || { dmax: key, strength: 0, count: 0 };
    shell.strength += norm;
    shell.count += 1;
    shellMap.set(key, shell);
    return {
      id: block.id,
      atoms: block.atoms,
      atomInstances: [
        { atomId: block.atoms[0], species: atomRecords[0].species, cell: [0, 0, 0], key: atomInstanceKey(block.atoms[0], [0, 0, 0]) },
        { atomId: block.atoms[1], species: atomRecords[1].species, cell: rjIndex, key: atomInstanceKey(block.atoms[1], rjIndex) },
        { atomId: block.atoms[2], species: atomRecords[2].species, cell: rkIndex, key: atomInstanceKey(block.atoms[2], rkIndex) },
      ],
      rj: block.rj,
      rk: block.rk,
      rjIndex,
      rkIndex,
      vertices: [xi, xj, xk],
      values: block.values,
      norm,
      dmax,
      perimeter,
    };
  });
  blocks.sort((a, b) => a.id - b.id);
  const norms = blocks.map((block) => block.norm);
  const dmaxValues = blocks.map((block) => block.dmax);
  const perimeterValues = blocks.map((block) => block.perimeter);
  const componentValues = blocks.flatMap((block) => block.values);
  const componentAbsNonzero = componentValues.map(Math.abs).filter((value) => value > 0);
  const componentZeroCount = componentValues.length - componentAbsNonzero.length;
  return {
    meta: {
      fc3File,
      qeFile: structure.source || "current structure",
      cutoffNm: null,
      repeat: structure.repeat,
      blockCount: blocks.length,
      normMin: Math.min(...norms),
      normMax: Math.max(...norms),
      dmaxMin: Math.min(...dmaxValues),
      dmaxMax: Math.max(...dmaxValues),
      perimeterMin: Math.min(...perimeterValues),
      perimeterMax: Math.max(...perimeterValues),
      componentCount: componentValues.length,
      componentZeroCount,
      componentAbsMin: componentAbsNonzero.length ? Math.min(...componentAbsNonzero) : 0,
      componentAbsMax: componentAbsNonzero.length ? Math.max(...componentAbsNonzero) : 0,
    },
    cell: structure.cell,
    basis: structure.basis,
    blocks,
    shells: [...shellMap.values()].sort((a, b) => a.dmax - b.dmax),
    topBlocks: [...blocks].sort((a, b) => b.norm - a.norm).slice(0, 40),
  };
}

function structureFromPayload(payload) {
  return {
    cell: payload.cell,
    basis: payload.basis,
    repeat: payload.meta.repeat,
    source: payload.meta.qeFile,
  };
}

const ELEMENT_COLORS = {
  H: "#f5f5f5",
  C: "#303030",
  N: "#3050f8",
  O: "#ff0d0d",
  Si: "#daa520",
  S: "#ffff30",
  P: "#ff8000",
};

const COVALENT_RADII = {
  H: 0.31,
  C: 0.76,
  N: 0.71,
  O: 0.66,
  Si: 1.11,
  S: 1.05,
  P: 1.07,
};

function speciesColor(species) {
  const key = normalizeSpecies(species);
  if (ELEMENT_COLORS[key]) return ELEMENT_COLORS[key];
  let hash = 0;
  for (const char of key) hash = (hash * 31 + char.charCodeAt(0)) % 360;
  return `hsl(${hash}, 62%, 46%)`;
}

function covalentRadius(species) {
  return COVALENT_RADII[normalizeSpecies(species)] || 0.9;
}

function atomKey(atom) {
  return `${normalizeSpecies(atom.species)}${atom.id}`;
}

function atomInstanceKey(atomId, cellIndex) {
  return `${atomId}@${cellIndex.join(",")}`;
}

function atomInstanceLabel(atomId, species, cellIndex) {
  return `${normalizeSpecies(species)}${atomId} (${cellIndex.join(",")})`;
}

function reset3DInteractionState() {
  state.selectedBlockId = null;
  state.pickedAtom = null;
  state.hiddenBlockIds = new Set();
  state.hiddenAtomKeys = new Set();
  state.hiddenElements = new Set();
  state.collapsedTree = new Set();
  state.objectTreeInitialized = false;
}

function normColor(t) {
  const stops = [
    [0.0, [0, 184, 255]],
    [0.28, [0, 230, 118]],
    [0.52, [255, 238, 0]],
    [0.74, [255, 145, 0]],
    [1.0, [255, 23, 68]],
  ];
  for (let i = 1; i < stops.length; i += 1) {
    if (t <= stops[i][0]) {
      const [aT, a] = stops[i - 1];
      const [bT, b] = stops[i];
      const local = (t - aT) / (bT - aT);
      const rgb = a.map((value, axis) => Math.round(value + (b[axis] - value) * local));
      return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
    }
  }
  return "rgb(255, 23, 68)";
}

function setupTabs() {
  document.querySelectorAll(".tab").forEach((button) => {
    button.addEventListener("click", () => {
      activateTab(button.dataset.tab);
      history.replaceState(null, "", `#${button.dataset.tab}`);
    });
  });
  window.addEventListener("hashchange", () => activateTab(location.hash.slice(1)));
}

function setupIntroTabs() {
  document.querySelectorAll(".intro-tab").forEach((button) => {
    button.addEventListener("click", () => {
      const target = button.dataset.introTab;
      document.querySelectorAll(".intro-tab").forEach((item) => {
        item.classList.toggle("is-active", item === button);
      });
      document.querySelectorAll(".intro-pane").forEach((pane) => {
        pane.classList.toggle("is-active", pane.id === target);
      });
    });
  });
}

function activateTab(tabId) {
  const id = $(tabId) ? tabId : "intro";
  document.querySelectorAll(".tab").forEach((item) => {
    item.classList.toggle("is-active", item.dataset.tab === id);
  });
  document.querySelectorAll(".panel").forEach((item) => {
    item.classList.toggle("is-active", item.id === id);
  });
  $("globalDataPanel").style.display = id === "compare" || id === "intro" ? "none" : "";
  requestAnimationFrame(drawAll);
}

function setViewDirection(value) {
  const map = {
    "100": { x: 0, y: Math.PI / 2, z: 0 },
    "010": { x: -Math.PI / 2, y: 0, z: 0 },
    "001": { x: 0, y: 0, z: 0 },
    "110": { x: -0.65, y: 0.85, z: 0 },
    "111": { x: -0.78, y: 0.62, z: 0.35 },
  };
  if (map[value]) {
    state.rotation = { ...map[value] };
    drawScene();
  }
}

function resetView() {
  state.rotation = { x: -0.55, y: 0.72, z: 0 };
  state.zoomScale = 1;
  $("viewSelect").value = "free";
  drawScene();
}

function normBounds() {
  const min = Math.max(state.data?.meta.normMin || 1e-12, 1e-12);
  const max = Math.max(state.data?.meta.normMax || min * 10, min * 1.0001);
  return { min, max, logMin: Math.log10(min), logMax: Math.log10(max) };
}

function normToSlider(value) {
  const bounds = normBounds();
  const logValue = Math.log10(Math.max(value, bounds.min));
  return Math.round(((logValue - bounds.logMin) / (bounds.logMax - bounds.logMin || 1)) * 1000);
}

function sliderToNorm(value) {
  const bounds = normBounds();
  const t = Number(value) / 1000;
  return 10 ** (bounds.logMin + t * (bounds.logMax - bounds.logMin));
}

function formatNormInput(value) {
  return Number(value).toFixed(3);
}

function clampNormRange(min, max) {
  const bounds = normBounds();
  let nextMin = Math.min(Math.max(min, bounds.min), bounds.max);
  let nextMax = Math.min(Math.max(max, bounds.min), bounds.max);
  if (nextMin > nextMax) [nextMin, nextMax] = [nextMax, nextMin];
  return { min: nextMin, max: nextMax };
}

function setNormRange(min, max, source = "program") {
  if (!state.data) return;
  const next = clampNormRange(min, max);
  state.normRange = next;
  const minSlider = $("normMinSlider");
  const maxSlider = $("normMaxSlider");
  const minInput = $("normMinInput");
  const maxInput = $("normMaxInput");
  if (source !== "slider") {
    minSlider.value = String(normToSlider(next.min));
    maxSlider.value = String(normToSlider(next.max));
  }
  if (source !== "input") {
    minInput.value = formatNormInput(next.min);
    maxInput.value = formatNormInput(next.max);
  }
  $("normRangeLabel").textContent = `${next.min.toFixed(3)} to ${next.max.toFixed(3)} eV/Å³, log slider`;
  drawScene();
}

function initializeNormRange() {
  if (!state.data || !$("normMinInput")) return;
  setNormRange(state.data.meta.normMin, state.data.meta.normMax);
}

function rotatePoint(point) {
  let [x, y, z] = point;
  const rx = state.rotation.x;
  const ry = state.rotation.y;
  const rz = state.rotation.z;
  let cy = Math.cos(ry);
  let sy = Math.sin(ry);
  [x, z] = [x * cy + z * sy, -x * sy + z * cy];
  let cx = Math.cos(rx);
  let sx = Math.sin(rx);
  [y, z] = [y * cx - z * sx, y * sx + z * cx];
  let cz = Math.cos(rz);
  let sz = Math.sin(rz);
  [x, y] = [x * cz - y * sz, x * sz + y * cz];
  return [x, y, z];
}

function project(point, center, zoom, canvas) {
  const rotated = rotatePoint(sub(point, center));
  return {
    x: canvas.width / 2 + rotated[0] * zoom,
    y: canvas.height / 2 - rotated[1] * zoom,
    z: rotated[2],
  };
}

function numericInput(id, fallback) {
  const value = Number($(id).value);
  return Number.isFinite(value) ? Math.trunc(value) : fallback;
}

function readBoundary() {
  const boundary = {
    aMin: numericInput("aMin", -1),
    aMax: numericInput("aMax", 1),
    bMin: numericInput("bMin", -1),
    bMax: numericInput("bMax", 1),
    cMin: numericInput("cMin", -1),
    cMax: numericInput("cMax", 1),
  };
  if (boundary.aMin > boundary.aMax) [boundary.aMin, boundary.aMax] = [boundary.aMax, boundary.aMin];
  if (boundary.bMin > boundary.bMax) [boundary.bMin, boundary.bMax] = [boundary.bMax, boundary.bMin];
  if (boundary.cMin > boundary.cMax) [boundary.cMin, boundary.cMax] = [boundary.cMax, boundary.cMin];
  return boundary;
}

function indexRange(min, max) {
  const values = [];
  for (let value = min; value <= max; value += 1) values.push(value);
  return values;
}

function fc3BoundaryRange() {
  const mins = [0, 0, 0];
  const maxs = [0, 0, 0];
  for (const block of state.data.blocks) {
    for (const index of [[0, 0, 0], block.rjIndex, block.rkIndex]) {
      for (let axis = 0; axis < 3; axis += 1) {
        mins[axis] = Math.min(mins[axis], index[axis]);
        maxs[axis] = Math.max(maxs[axis], index[axis]);
      }
    }
  }
  return { mins, maxs };
}

function initializeBoundaryFromFc3() {
  if (!state.data || !$("aMin")) return;
  const { mins, maxs } = fc3BoundaryRange();
  ["a", "b", "c"].forEach((axis, idx) => {
    $(`${axis}Min`).value = String(mins[idx]);
    $(`${axis}Max`).value = String(maxs[idx]);
  });
  drawScene();
}

function sceneAtoms(boundary) {
  const atoms = [];
  for (const a of indexRange(boundary.aMin, boundary.aMax)) {
    for (const b of indexRange(boundary.bMin, boundary.bMax)) {
      for (const c of indexRange(boundary.cMin, boundary.cMax)) {
        const shift = dotCell([a, b, c], state.data.cell);
        for (const atom of state.data.basis) {
          const cell = [a, b, c];
          atoms.push({
            id: atom.id,
            species: atom.species,
            cell,
            instanceKey: atomInstanceKey(atom.id, cell),
            pos: add(atom.cart, shift),
          });
        }
      }
    }
  }
  return atoms;
}

function sceneBonds(atoms) {
  const bonds = [];
  for (let i = 0; i < atoms.length; i += 1) {
    for (let j = i + 1; j < atoms.length; j += 1) {
      const d = length(sub(atoms[i].pos, atoms[j].pos));
      const cutoff = (covalentRadius(atoms[i].species) + covalentRadius(atoms[j].species)) * 1.25;
      if (d > 0.25 && d <= cutoff) bonds.push([atoms[i], atoms[j]]);
    }
  }
  return bonds;
}

function indexInside(index, boundary) {
  return (
    index[0] >= boundary.aMin &&
    index[0] <= boundary.aMax &&
    index[1] >= boundary.bMin &&
    index[1] <= boundary.bMax &&
    index[2] >= boundary.cMin &&
    index[2] <= boundary.cMax
  );
}

function visibleFc3Blocks(boundary, normRange) {
  const baseBlocks = state.data.blocks.filter((block) => {
    if (block.norm < normRange.min || block.norm > normRange.max) return false;
    if (isBlockHidden(block)) return false;
    return indexInside(block.rjIndex, boundary) && indexInside(block.rkIndex, boundary);
  });
  const selected = state.data.blocks.find((block) => block.id === state.selectedBlockId);
  if (selected && !isBlockHidden(selected) && !baseBlocks.some((block) => block.id === selected.id)) baseBlocks.push(selected);
  return baseBlocks;
}

function isBlockHidden(block) {
  if (state.hiddenBlockIds.has(block.id)) return true;
  return block.atoms.some((atomId) => {
    const atom = state.data.basis[atomId - 1];
    if (!atom) return false;
    return state.hiddenElements.has(atom.species) || state.hiddenAtomKeys.has(atomKey(atom));
  });
}

function sceneBounds(points) {
  const min = [Infinity, Infinity, Infinity];
  const max = [-Infinity, -Infinity, -Infinity];
  for (const point of points) {
    for (let axis = 0; axis < 3; axis += 1) {
      min[axis] = Math.min(min[axis], point[axis]);
      max[axis] = Math.max(max[axis], point[axis]);
    }
  }
  return {
    center: scale(add(min, max), 0.5),
    span: Math.max(max[0] - min[0], max[1] - min[1], max[2] - min[2], 1),
  };
}

function prepareSceneCanvas(canvas) {
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.max(1, Math.min(2.5, window.devicePixelRatio || 1));
  canvas.width = Math.max(900, Math.floor(rect.width * dpr));
  canvas.height = Math.max(560, Math.floor(rect.height * dpr));
  return dpr;
}

function drawLine(ctx, canvas, a, b, center, zoom, style) {
  const pa = project(a, center, zoom, canvas);
  const pb = project(b, center, zoom, canvas);
  ctx.strokeStyle = style.color;
  ctx.lineWidth = style.width;
  ctx.globalAlpha = style.alpha ?? 1;
  ctx.beginPath();
  ctx.moveTo(pa.x, pa.y);
  ctx.lineTo(pb.x, pb.y);
  ctx.stroke();
  ctx.globalAlpha = 1;
}

function drawScene() {
  if (!state.data) return;
  const canvas = $("scene3d");
  const ctx = canvas.getContext("2d");
  const dpr = prepareSceneCanvas(canvas);
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  const boundary = readBoundary();
  if (!state.normRange.max) initializeNormRange();
  const normRange = state.normRange.max ? state.normRange : { min: state.data.meta.normMin, max: state.data.meta.normMax };

  const atoms = sceneAtoms(boundary);
  const blocks = visibleFc3Blocks(boundary, normRange);
  const points = atoms.map((atom) => atom.pos).concat(blocks.flatMap((block) => block.vertices));
  const { center, span } = sceneBounds(points);
  const zoom = ((Math.min(canvas.width, canvas.height) * 0.72) / span) * state.zoomScale;

  if ($("showCell").checked) drawCellBoundaries(ctx, canvas, boundary, center, zoom, dpr);
  if ($("showBonds").checked) {
    for (const [a, b] of sceneBonds(atoms)) {
      drawLine(ctx, canvas, a.pos, b.pos, center, zoom, { color: "rgba(25, 25, 25, 0.32)", width: 1.4 * dpr });
    }
  }
  if ($("showFc3").checked) {
    const selected = blocks.find((block) => block.id === state.selectedBlockId);
    const pickedBlocks = state.pickedAtom
      ? blocks.filter((block) => blockMatchesPickedAtom(block, state.pickedAtom) && block.norm >= normRange.min && block.norm <= normRange.max)
      : [];
    const pickedIds = new Set(pickedBlocks.map((block) => block.id));
    const isolatePicked = $("isolatePickedAtom")?.checked && state.pickedAtom;
    const normalBlocks = isolatePicked ? [] : blocks.filter((block) => block.id !== state.selectedBlockId && !pickedIds.has(block.id)).sort((a, b) => a.norm - b.norm);
    for (const block of normalBlocks) drawFc3Block(ctx, canvas, block, center, zoom, dpr, "normal");
    for (const block of pickedBlocks) drawFc3Block(ctx, canvas, block, center, zoom, dpr, "picked");
    if (selected && (!isolatePicked || pickedIds.has(selected.id))) drawFc3Block(ctx, canvas, selected, center, zoom, dpr, "selected");
  }
  if ($("showAtoms").checked) drawAtoms(ctx, canvas, atoms, center, zoom, dpr);

  $("sceneStats").innerHTML = [
    `FC3 file: ${state.data.meta.fc3File}`,
    `Structure: ${state.data.meta.qeFile}`,
    `Blocks shown: ${blocks.length} / ${state.data.meta.blockCount}`,
    `Norm range: ${normRange.min.toFixed(3)} to ${normRange.max.toFixed(3)} eV/Å³`,
    `Atoms shown: ${atoms.length}`,
    `Zoom: ${state.zoomScale.toFixed(2)}x`,
    state.pickedAtom ? `Picked atom: ${state.pickedAtom.label}; matched blocks: ${pickedBlocks.length}` : "Picked atom: none",
    $("isolatePickedAtom")?.checked && state.pickedAtom ? "Pick isolate: on" : "Pick isolate: off",
    state.selectedBlockId ? `Highlight: block ${state.selectedBlockId}` : "Highlight: none",
  ].join("<br>");
}

function drawFc3Block(ctx, canvas, block, center, zoom, dpr, mode) {
  const t = Math.sqrt(block.norm / state.data.meta.normMax);
  const color = mode === "selected" ? "#ff00d4" : normColor(t);
  const width = mode === "selected" ? 8 * dpr : mode === "picked" ? (6.2 + 7.2 * t) * dpr : (1.2 + 4.2 * t) * dpr;
  const alpha = mode === "normal" ? 0.78 : 1;
  drawLine(ctx, canvas, block.vertices[0], block.vertices[1], center, zoom, { color, width, alpha });
  drawLine(ctx, canvas, block.vertices[1], block.vertices[2], center, zoom, { color, width, alpha });
  drawLine(ctx, canvas, block.vertices[2], block.vertices[0], center, zoom, { color, width, alpha });
}

function drawCellBoundaries(ctx, canvas, boundary, center, zoom, dpr) {
  const corners = [
    [boundary.aMin, boundary.bMin, boundary.cMin],
    [boundary.aMax + 1, boundary.bMin, boundary.cMin],
    [boundary.aMax + 1, boundary.bMax + 1, boundary.cMin],
    [boundary.aMin, boundary.bMax + 1, boundary.cMin],
    [boundary.aMin, boundary.bMin, boundary.cMax + 1],
    [boundary.aMax + 1, boundary.bMin, boundary.cMax + 1],
    [boundary.aMax + 1, boundary.bMax + 1, boundary.cMax + 1],
    [boundary.aMin, boundary.bMax + 1, boundary.cMax + 1],
  ].map((idx) => dotCell(idx, state.data.cell));
  const edges = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]];
  for (const [a, b] of edges) {
    drawLine(ctx, canvas, corners[a], corners[b], center, zoom, { color: "rgba(0,0,0,0.58)", width: 1.2 * dpr });
  }
}

function drawAtoms(ctx, canvas, atoms, center, zoom, dpr) {
  const projected = atoms.map((atom) => ({ atom, p: project(atom.pos, center, zoom, canvas) })).sort((a, b) => a.p.z - b.p.z);
  state.scenePickAtoms = projected;
  for (const item of projected) {
    const radius = Math.max(4.5, 7.5 + item.p.z * 0.03) * dpr;
    const color = speciesColor(item.atom.species);
    const grad = ctx.createRadialGradient(item.p.x - radius * 0.35, item.p.y - radius * 0.35, radius * 0.2, item.p.x, item.p.y, radius);
    grad.addColorStop(0, "#ffffff");
    grad.addColorStop(0.35, color);
    grad.addColorStop(1, shadeColor(color, -42));
    ctx.fillStyle = grad;
    const picked = state.pickedAtom && state.pickedAtom.key === item.atom.instanceKey;
    ctx.strokeStyle = picked ? "#ff00d4" : "#111";
    ctx.lineWidth = (picked ? 2.4 : 0.8) * dpr;
    ctx.beginPath();
    ctx.arc(item.p.x, item.p.y, radius, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
  }
}

function blockMatchesPickedAtom(block, pickedAtom) {
  return Boolean(pickedAtom && block.atomInstances?.some((atom) => atom.key === pickedAtom.key));
}

function shadeColor(color, percent) {
  if (!color.startsWith("#")) return color;
  const num = Number.parseInt(color.slice(1), 16);
  const amt = Math.round(2.55 * percent);
  const r = Math.max(0, Math.min(255, (num >> 16) + amt));
  const g = Math.max(0, Math.min(255, ((num >> 8) & 0xff) + amt));
  const b = Math.max(0, Math.min(255, (num & 0xff) + amt));
  return `rgb(${r}, ${g}, ${b})`;
}

function buildObjectTree() {
  if (!state.data) return;
  const groups = new Map();
  for (const block of state.data.blocks) {
    const uniqueAtoms = [...new Set(block.atoms)];
    for (const atomId of uniqueAtoms) {
      const atom = state.data.basis[atomId - 1];
      if (!atom) continue;
      if (!groups.has(atom.species)) groups.set(atom.species, new Map());
      const sites = groups.get(atom.species);
      const label = `${atom.species}${atom.id}`;
      if (!sites.has(label)) sites.set(label, []);
      sites.get(label).push(block);
    }
  }

  const root = $("objectTree");
  root.innerHTML = "";
  const elementHeader = document.createElement("div");
  elementHeader.className = "tree-section-label";
  elementHeader.textContent = "Element";
  root.appendChild(elementHeader);
  if (!state.objectTreeInitialized) {
    state.collapsedTree = new Set([...groups.keys()].map((element) => `element:${element}`));
    state.objectTreeInitialized = true;
  }
  for (const [element, sites] of groups) {
    const group = document.createElement("div");
    group.className = "tree-group";
    const title = document.createElement("button");
    title.className = "tree-title";
    title.type = "button";
    const elementHidden = state.hiddenElements.has(element);
    const elementCollapsed = state.collapsedTree.has(`element:${element}`);
    title.innerHTML = `<input type="checkbox" ${elementHidden ? "checked" : ""} aria-label="Hide ${element}" /> <span>${element}</span><span class="tree-expander">${elementCollapsed ? "▸" : "▾"}</span>`;
    title.querySelector("input").addEventListener("click", (event) => {
      event.stopPropagation();
      toggleSet(state.hiddenElements, element, event.target.checked);
      drawScene();
    });
    title.addEventListener("click", () => {
      toggleSet(state.collapsedTree, `element:${element}`);
      buildObjectTree();
    });
    group.appendChild(title);
    if (elementCollapsed) {
      root.appendChild(group);
      continue;
    }
    for (const [label, blocks] of sites) {
      blocks.sort((a, b) => b.norm - a.norm);
      const site = document.createElement("div");
      site.className = "tree-site";
      const siteTitle = document.createElement("button");
      siteTitle.className = "tree-site-title";
      siteTitle.type = "button";
      const siteHidden = state.hiddenAtomKeys.has(label);
      const siteCollapsed = state.collapsedTree.has(`site:${label}`);
      siteTitle.innerHTML = `<input type="checkbox" ${siteHidden ? "checked" : ""} aria-label="Hide ${label}" /> <span>${label} (${blocks.length})</span><span class="tree-expander">${siteCollapsed ? "▸" : "▾"}</span>`;
      siteTitle.querySelector("input").addEventListener("click", (event) => {
        event.stopPropagation();
        toggleSet(state.hiddenAtomKeys, label, event.target.checked);
        drawScene();
      });
      siteTitle.addEventListener("click", () => {
        toggleSet(state.collapsedTree, `site:${label}`);
        buildObjectTree();
      });
      site.appendChild(siteTitle);
      if (siteCollapsed) {
        group.appendChild(site);
        continue;
      }
      for (const block of blocks) {
        const row = document.createElement("div");
        row.className = "tree-block-row";
        const checkbox = document.createElement("input");
        checkbox.type = "checkbox";
        checkbox.checked = state.hiddenBlockIds.has(block.id);
        checkbox.title = "Hide this FC3 block";
        checkbox.addEventListener("click", (event) => {
          event.stopPropagation();
          toggleSet(state.hiddenBlockIds, block.id, checkbox.checked);
          drawScene();
        });
        const button = document.createElement("button");
        button.className = "tree-block";
        if (block.id === state.selectedBlockId) button.classList.add("is-selected");
        button.type = "button";
        button.dataset.blockId = String(block.id);
        button.innerHTML = [
          `<strong>block ${block.id}</strong> ||Phi3||F=${block.norm.toFixed(4)}`,
          `(i,j,k)=(${block.atoms.join(",")})`,
          `Rj=(${block.rjIndex.join(",")}) Rk=(${block.rkIndex.join(",")})`,
        ].join("<br>");
        button.addEventListener("click", () => {
          state.selectedBlockId = block.id;
          buildObjectTree();
          drawScene();
        });
        row.appendChild(checkbox);
        row.appendChild(button);
        site.appendChild(row);
      }
      group.appendChild(site);
    }
    root.appendChild(group);
  }
}

function toggleSet(set, key, force = null) {
  const shouldAdd = force === null ? !set.has(key) : force;
  if (shouldAdd) set.add(key);
  else set.delete(key);
}

function setupSceneEvents() {
  const canvas = $("scene3d");
  canvas.addEventListener("pointerdown", (event) => {
    state.dragging = true;
    state.lastPointer = [event.clientX, event.clientY];
    state.pointerDown = [event.clientX, event.clientY];
    canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener("pointermove", (event) => {
    if (!state.dragging) return;
    const [x, y] = state.lastPointer;
    state.rotation.y += (event.clientX - x) * 0.008;
    state.rotation.x += (event.clientY - y) * 0.008;
    state.lastPointer = [event.clientX, event.clientY];
    drawScene();
  });
  canvas.addEventListener("pointerup", (event) => {
    const down = state.pointerDown;
    state.dragging = false;
    state.pointerDown = null;
    if ($("pickAtomMode").checked && down) {
      const moved = Math.hypot(event.clientX - down[0], event.clientY - down[1]);
      if (moved < 5) pickAtomAtEvent(event);
    }
  });
  canvas.addEventListener("wheel", (event) => {
    event.preventDefault();
    const factor = event.deltaY < 0 ? 1.12 : 0.89;
    state.zoomScale = Math.max(0.25, Math.min(6, state.zoomScale * factor));
    drawScene();
  }, { passive: false });

  ["aMin", "aMax", "bMin", "bMax", "cMin", "cMax", "showAtoms", "showBonds", "showFc3", "showCell", "isolatePickedAtom"].forEach((id) => {
    $(id).addEventListener("input", drawScene);
  });
  $("pickAtomMode").addEventListener("input", () => {
    if (!$("pickAtomMode").checked) state.pickedAtom = null;
    drawScene();
  });
  ["normMinSlider", "normMaxSlider"].forEach((id) => {
    $(id).addEventListener("input", () => {
      setNormRange(sliderToNorm($("normMinSlider").value), sliderToNorm($("normMaxSlider").value), "slider");
    });
  });
  ["normMinInput", "normMaxInput"].forEach((id) => {
    $(id).addEventListener("change", () => {
      setNormRange(Number($("normMinInput").value), Number($("normMaxInput").value), "input");
    });
  });
  $("zoomInButton").addEventListener("click", () => {
    state.zoomScale = Math.min(6, state.zoomScale * 1.2);
    drawScene();
  });
  $("zoomOutButton").addEventListener("click", () => {
    state.zoomScale = Math.max(0.25, state.zoomScale / 1.2);
    drawScene();
  });
  $("resetViewButton").addEventListener("click", resetView);
  $("applyViewButton").addEventListener("click", () => setViewDirection($("viewSelect").value));
  $("clearHighlightButton").addEventListener("click", () => {
    state.selectedBlockId = null;
    buildObjectTree();
    drawScene();
  });
  $("viewSelect").addEventListener("input", (event) => setViewDirection(event.target.value));
}

function pickAtomAtEvent(event) {
  const canvas = $("scene3d");
  const rect = canvas.getBoundingClientRect();
  const scaleX = canvas.width / rect.width;
  const scaleY = canvas.height / rect.height;
  const x = (event.clientX - rect.left) * scaleX;
  const y = (event.clientY - rect.top) * scaleY;
  const boundary = readBoundary();
  const normRange = state.normRange.max ? state.normRange : { min: state.data.meta.normMin, max: state.data.meta.normMax };
  const visibleBlocks = visibleFc3Blocks(boundary, normRange);
  const hitRadius = 42 * (window.devicePixelRatio || 1);
  const candidates = [];
  for (const item of state.scenePickAtoms) {
    const distance = Math.hypot(item.p.x - x, item.p.y - y);
    if (distance > hitRadius) continue;
    const key = item.atom.instanceKey;
    const blockCount = visibleBlocks.filter((block) => block.atomInstances?.some((atom) => atom.key === key)).length;
    candidates.push({ ...item, distance, blockCount });
  }
  candidates.sort((a, b) => {
    if (a.blockCount && !b.blockCount) return -1;
    if (!a.blockCount && b.blockCount) return 1;
    if (a.blockCount !== b.blockCount) return b.blockCount - a.blockCount;
    return a.distance - b.distance;
  });
  const best = candidates[0];
  if (best) {
    state.pickedAtom = {
      id: best.atom.id,
      species: best.atom.species,
      cell: [...best.atom.cell],
      key: best.atom.instanceKey,
      label: atomInstanceLabel(best.atom.id, best.atom.species, best.atom.cell),
      visibleBlockCount: best.blockCount,
    };
    drawScene();
  }
}

function setupUploadEvents() {
  $("fc3FileInput").addEventListener("change", () => {
    $("fc3FileDisplay").textContent = $("fc3FileInput").files[0]?.name || "FORCE_CONSTANTS_3RD_Example";
  });
  $("qeFileInput").addEventListener("change", () => {
    $("qeFileDisplay").textContent = $("qeFileInput").files[0]?.name || "BASE.si_supper.scf.in_example";
  });
  $("restoreDefaultButton").addEventListener("click", () => {
    state.data = structuredClone(state.defaultData);
    state.structure = structureFromPayload(state.data);
    reset3DInteractionState();
    $("dataStatus").textContent = "Default example restored: FORCE_CONSTANTS_3RD_Example and BASE.si_supper.scf.in_example.";
    $("fc3FileDisplay").textContent = "FORCE_CONSTANTS_3RD_Example";
    $("qeFileDisplay").textContent = "BASE.si_supper.scf.in_example";
    initializeNormRange();
    initializeBoundaryFromFc3();
    resetView();
    buildObjectTree();
    drawAll();
  });
  $("loadUploadButton").addEventListener("click", async () => {
    try {
      const fc3File = $("fc3FileInput").files[0];
      const qeFile = $("qeFileInput").files[0];
      let structure = state.structure;
      if (qeFile) structure = parseQeStructure(await qeFile.text(), qeFile.name);
      if (!structure) throw new Error("No structure is available. Upload a QE structure file first.");
      if (!fc3File) {
        state.structure = structure;
        $("dataStatus").textContent = `Structure loaded: ${structure.source}. Upload an FC3 file to rebuild FC3 data.`;
        $("qeFileDisplay").textContent = structure.source;
        return;
      }
      const rawBlocks = parseFc3(await fc3File.text(), fc3File.name);
      const payload = buildPayloadFromFc3Blocks(rawBlocks, structure, fc3File.name);
      state.data = payload;
      state.structure = structure;
      reset3DInteractionState();
      $("dataStatus").textContent = `Loaded ${fc3File.name} with ${payload.meta.blockCount} blocks; structure=${payload.meta.qeFile}.`;
      $("fc3FileDisplay").textContent = fc3File.name;
      $("qeFileDisplay").textContent = payload.meta.qeFile;
      initializeNormRange();
      initializeBoundaryFromFc3();
      resetView();
      buildObjectTree();
      updateCompareStatus();
      drawAll();
    } catch (error) {
      $("dataStatus").textContent = error.message;
      console.error(error);
    }
  });
}

function setupCompareEvents() {
  initializeCompareSlots();
  $("loadCompareStructureButton").addEventListener("click", async () => {
    try {
      const file = $("compareStructureInput").files[0];
      if (!file) throw new Error("Select a QE structure file for comparison.");
      state.compareStructure = parseQeStructure(await file.text(), file.name);
      $("compareStructureStatus").textContent = `Compare structure loaded: ${state.compareStructure.source}.`;
      updateCompareStatus();
    } catch (error) {
      $("compareStructureStatus").textContent = error.message;
      console.error(error);
    }
  });
  $("loadCompareBulkButton").addEventListener("click", loadCompareBulkInput);
}

function initializeCompareSlots() {
  state.compareSlots = Array.from({ length: 4 }, (_, idx) => ({
    dataset: null,
    legend: `FC3 ${idx + 1}`,
    color: COMPARE_COLORS[idx],
  }));
  const container = $("compareSlots");
  container.innerHTML = state.compareSlots.map((slot, idx) => `
    <div class="compare-slot" data-slot="${idx}">
      <h4>FC3 slot ${idx + 1}</h4>
      <label>
        FORCE_CONSTANTS_3RD
        <input id="compareFileInput${idx}" type="file" />
      </label>
      <div class="compare-slot-meta">
        <label>
          Legend
          <input id="compareLegend${idx}" type="text" value="${escapeHtml(slot.legend)}" />
        </label>
        <label>
          Color
          <input id="compareColor${idx}" type="color" value="${slot.color}" />
        </label>
      </div>
      <div class="compare-slot-controls">
        <button id="loadCompareButton${idx}" class="mini-command">Load Compare</button>
        <button id="unloadCompareButton${idx}" class="mini-command">Unload Compare</button>
      </div>
      <p id="compareSlotStatus${idx}" class="compare-slot-status"></p>
    </div>
  `).join("");
  state.compareSlots.forEach((_, idx) => {
    $(`compareFileInput${idx}`).addEventListener("change", () => {
      state.compareSlots[idx].pendingFileName = $(`compareFileInput${idx}`).files[0]?.name || null;
    });
    $(`compareLegend${idx}`).addEventListener("input", (event) => {
      state.compareSlots[idx].legend = event.target.value || `FC3 ${idx + 1}`;
      updateCompareStatus();
      drawCompareCharts();
    });
    $(`compareColor${idx}`).addEventListener("input", (event) => {
      state.compareSlots[idx].color = event.target.value;
      drawCompareCharts();
    });
    $(`loadCompareButton${idx}`).addEventListener("click", () => loadCompareSlot(idx));
    $(`unloadCompareButton${idx}`).addEventListener("click", () => unloadCompareSlot(idx));
  });
  updateCompareSlotStatuses();
}

async function loadCompareSlot(idx) {
  try {
    const file = $(`compareFileInput${idx}`).files[0];
    if (!file) throw new Error("Select a FORCE_CONSTANTS_3RD file.");
    const structure = state.compareStructure;
    if (!structure) throw new Error("Load a Compare structure before loading FC3 compare files.");
    const rawBlocks = parseFc3(await file.text(), file.name);
    const dataset = buildPayloadFromFc3Blocks(rawBlocks, structure, file.name);
    const fileChanged = state.compareSlots[idx].fileName !== file.name;
    state.compareSlots[idx].dataset = dataset;
    state.compareSlots[idx].fileName = file.name;
    if (fileChanged) {
      state.compareSlots[idx].legend = file.name;
      $(`compareLegend${idx}`).value = file.name;
    }
    updateCompareStatus();
    drawCompareCharts();
  } catch (error) {
    $(`compareSlotStatus${idx}`).textContent = error.message;
    console.error(error);
  }
}

async function loadCompareBulkInput() {
  try {
    const files = Array.from($("compareBulkFc3Input").files || []).slice(0, 4);
    if (!files.length) throw new Error("Select one to four FORCE_CONSTANTS_3RD files.");
    const structure = state.compareStructure;
    if (!structure) throw new Error("Load a Compare structure before loading batch FC3 files.");

    for (let idx = 0; idx < 4; idx += 1) {
      const slot = state.compareSlots[idx];
      const file = files[idx];
      if (!file) {
        slot.dataset = null;
        slot.fileName = null;
        slot.pendingFileName = null;
        slot.legend = `FC3 ${idx + 1}`;
        const legendInput = $(`compareLegend${idx}`);
        if (legendInput) legendInput.value = slot.legend;
        continue;
      }
      const rawBlocks = parseFc3(await file.text(), file.name);
      const dataset = buildPayloadFromFc3Blocks(rawBlocks, structure, file.name);
      slot.dataset = dataset;
      slot.fileName = file.name;
      slot.pendingFileName = file.name;
      slot.legend = file.name;
      const legendInput = $(`compareLegend${idx}`);
      if (legendInput) legendInput.value = file.name;
    }
    $("compareBulkStatus").textContent = `Batch loaded ${files.length} FC3 file${files.length > 1 ? "s" : ""} into compare slots.`;
    updateCompareStatus();
    updateCompareSlotStatuses();
    drawCompareCharts();
  } catch (error) {
    $("compareBulkStatus").textContent = error.message;
    console.error(error);
  }
}

function unloadCompareSlot(idx) {
  if (idx === 0) {
    state.compareSlots[0].dataset = null;
  } else {
    state.compareSlots[idx].dataset = null;
  }
  updateCompareStatus();
  drawCompareCharts();
}

function escapeHtml(value) {
  return String(value)
    .replace(/&/g, "&amp;")
    .replace(/"/g, "&quot;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

function updateCompareStatus() {
  const datasets = compareDatasets();
  $("compareStatus").innerHTML = datasets.length
    ? datasets
      .map((entry) => `${escapeHtml(entry.legend)}: ${entry.dataset.meta.blockCount} blocks, ${entry.dataset.meta.componentCount || componentValues(entry.dataset).length} tensor components`)
      .join("<br>")
    : "No FC3 compare slots are loaded.";
  updateCompareSlotStatuses();
}

function updateCompareSlotStatuses() {
  state.compareSlots.forEach((slot, idx) => {
    const target = $(`compareSlotStatus${idx}`);
    if (!target) return;
    target.textContent = slot.dataset ? `Loaded: ${slot.dataset.meta.fc3File}` : "Not loaded.";
  });
}

function niceTicks(min, max, count = 5) {
  if (!Number.isFinite(min) || !Number.isFinite(max) || min === max) return [min || 0];
  const span = max - min;
  if (span <= 5 && Math.abs(min - Math.round(min)) < 1e-9 && Math.abs(max - Math.round(max)) < 1e-9) {
    const ticks = [];
    for (let tick = Math.round(min); tick <= Math.round(max); tick += 1) ticks.push(tick);
    return ticks;
  }
  if (span >= 5) {
    const step = Math.max(1, Math.ceil(span / count));
    const start = Math.ceil(min / step) * step;
    const ticks = [];
    for (let tick = start; tick <= max + 1e-9; tick += step) ticks.push(tick);
    return ticks;
  }
  const ticks = [];
  for (let i = 0; i <= count; i += 1) ticks.push(min + ((max - min) * i) / count);
  return ticks;
}

function prepareChartCanvas(canvas) {
  if (!canvas) return null;
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.max(2, Math.min(3, window.devicePixelRatio || 1));
  const width = Math.max(900, Math.floor(rect.width * dpr));
  const height = Math.max(675, Math.floor(width * 0.75));
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }
  return dpr;
}

function chartFrame(canvas, xLabel, yLabel, xMin, xMax, yMin, yMax, options = {}) {
  const dpr = prepareChartCanvas(canvas);
  if (!dpr) return null;
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = "#fff";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  const pad = { l: 100 * dpr, r: 42 * dpr, t: 42 * dpr, b: 78 * dpr };
  if (options.grid) drawChartGrid(ctx, pad, canvas, xMin, xMax, yMin, yMax, options, dpr);
  ctx.strokeStyle = "#000";
  ctx.lineWidth = 1.5 * dpr;
  ctx.strokeRect(pad.l, pad.t, canvas.width - pad.l - pad.r, canvas.height - pad.t - pad.b);

  ctx.fillStyle = "#000";
  ctx.font = `${16 * dpr}px Arial`;
  ctx.textAlign = "center";
  ctx.fillText(xLabel, (pad.l + canvas.width - pad.r) / 2, canvas.height - 26 * dpr);
  ctx.save();
  ctx.translate(38 * dpr, (pad.t + canvas.height - pad.b) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();

  ctx.font = `${13 * dpr}px Arial`;
  ctx.lineWidth = 1 * dpr;
  if (options.logX || options.logY) {
    if (options.logY) drawLogTicks(ctx, "y", pad, canvas, yMin, yMax, dpr);
    else drawLinearTicks(ctx, "y", pad, canvas, yMin, yMax, dpr, options.formatY || formatTick, options.linearYTickMin);
    if (options.logX) drawLogTicks(ctx, "x", pad, canvas, xMin, xMax, dpr);
    else drawLinearTicks(ctx, "x", pad, canvas, xMin, xMax, dpr, options.formatX || formatTick, options.linearXTickMin);
  } else {
    drawLinearTicks(ctx, "y", pad, canvas, yMin, yMax, dpr, options.formatY || formatTick, options.linearYTickMin);
    drawLinearTicks(ctx, "x", pad, canvas, xMin, xMax, dpr, options.formatX || formatTick, options.linearXTickMin);
  }
  return {
    ctx,
    dpr,
    x: (value) => pad.l + ((value - xMin) / (xMax - xMin || 1)) * (canvas.width - pad.l - pad.r),
    y: (value) => canvas.height - pad.b - ((value - yMin) / (yMax - yMin || 1)) * (canvas.height - pad.t - pad.b),
    pad,
  };
}

function logChartFrame(canvas, xLabel, yLabel, xMin, xMax, yMin, yMax) {
  if (xMin <= 0 || yMin <= 0) throw new Error("Log plot values must be positive.");
  const frame = chartFrame(
    canvas,
    xLabel,
    yLabel,
    Math.log10(xMin),
    Math.log10(xMax),
    Math.log10(yMin),
    Math.log10(yMax),
    { logX: true, logY: true, grid: true },
  );
  return {
    ...frame,
    x: (value) => frame.x(Math.log10(value)),
    y: (value) => frame.y(Math.log10(value)),
  };
}

function linearLogYChartFrame(canvas, xLabel, yLabel, xMin, xMax, yMin, yMax, options = {}) {
  if (yMin <= 0) throw new Error("Log plot values must be positive.");
  const frame = chartFrame(
    canvas,
    xLabel,
    yLabel,
    xMin,
    xMax,
    Math.log10(yMin),
    Math.log10(yMax),
    { logY: true, grid: true, ...options },
  );
  return {
    ...frame,
    y: (value) => frame.y(Math.log10(value)),
  };
}

function drawChartGrid(ctx, pad, canvas, xMin, xMax, yMin, yMax, options, dpr) {
  const width = canvas.width - pad.l - pad.r;
  const height = canvas.height - pad.t - pad.b;
  ctx.save();
  ctx.strokeStyle = "#d9d9d9";
  ctx.lineWidth = 0.75 * dpr;
  const drawVertical = (value, major = true) => {
    const x = pad.l + ((value - xMin) / (xMax - xMin || 1)) * width;
    ctx.globalAlpha = major ? 0.75 : 0.36;
    ctx.beginPath();
    ctx.moveTo(x, pad.t);
    ctx.lineTo(x, canvas.height - pad.b);
    ctx.stroke();
  };
  const drawHorizontal = (value, major = true) => {
    const y = canvas.height - pad.b - ((value - yMin) / (yMax - yMin || 1)) * height;
    ctx.globalAlpha = major ? 0.75 : 0.36;
    ctx.beginPath();
    ctx.moveTo(pad.l, y);
    ctx.lineTo(canvas.width - pad.r, y);
    ctx.stroke();
  };
  if (options.logX) {
    for (const tick of logTickValues(xMin, xMax)) drawVertical(tick.value, tick.major);
  } else {
    linearTicksWithFloor(xMin, xMax, options.linearXTickMin).forEach((tick) => drawVertical(tick, true));
  }
  if (options.logY) {
    for (const tick of logTickValues(yMin, yMax)) drawHorizontal(tick.value, tick.major);
  } else {
    niceTicks(yMin, yMax).forEach((tick) => drawHorizontal(tick, true));
  }
  ctx.restore();
}

function logTickValues(logMin, logMax) {
  const ticks = [];
  const start = Math.floor(logMin);
  const end = Math.ceil(logMax);
  for (let exponent = start; exponent <= end; exponent += 1) {
    for (let multiplier = 1; multiplier <= 9; multiplier += 1) {
      const value = Math.log10(multiplier) + exponent;
      if (value < logMin - 1e-9 || value > logMax + 1e-9) continue;
      ticks.push({ value, major: multiplier === 1, exponent, multiplier });
    }
  }
  return ticks;
}

function drawLinearTicks(ctx, axis, pad, canvas, min, max, dpr, formatter, tickMin = null) {
  for (const tick of linearTicksWithFloor(min, max, tickMin)) {
    if (tickMin !== null && tick < tickMin - 1e-12) continue;
    if (axis === "y") {
      const y = canvas.height - pad.b - ((tick - min) / (max - min || 1)) * (canvas.height - pad.t - pad.b);
      ctx.strokeStyle = "#000";
      ctx.beginPath();
      ctx.moveTo(pad.l - 6 * dpr, y);
      ctx.lineTo(pad.l, y);
      ctx.stroke();
      ctx.textAlign = "right";
      ctx.fillText(formatter(tick), pad.l - 11 * dpr, y + 4 * dpr);
    } else {
      const x = pad.l + ((tick - min) / (max - min || 1)) * (canvas.width - pad.l - pad.r);
      ctx.strokeStyle = "#000";
      ctx.beginPath();
      ctx.moveTo(x, canvas.height - pad.b);
      ctx.lineTo(x, canvas.height - pad.b + 6 * dpr);
      ctx.stroke();
      ctx.textAlign = "center";
      ctx.fillText(formatter(tick), x, canvas.height - pad.b + 26 * dpr);
    }
  }
}

function linearTicksWithFloor(min, max, tickMin = null) {
  const ticks = tickMin !== null ? niceTicks(tickMin, max) : niceTicks(min, max);
  if (tickMin !== null && tickMin >= min - 1e-12 && tickMin <= max + 1e-12 && !ticks.some((tick) => Math.abs(tick - tickMin) < 1e-9)) {
    ticks.push(tickMin);
  }
  return ticks.sort((a, b) => a - b);
}

function drawLogTicks(ctx, axis, pad, canvas, logMin, logMax, dpr) {
  for (const tick of logTickValues(logMin, logMax)) {
      const value = tick.value;
      const isMajor = tick.major;
      if (axis === "y") {
        const y = canvas.height - pad.b - ((value - logMin) / (logMax - logMin || 1)) * (canvas.height - pad.t - pad.b);
        ctx.strokeStyle = "#000";
        ctx.lineWidth = isMajor ? 1.4 * dpr : 1 * dpr;
        ctx.beginPath();
        ctx.moveTo(pad.l - (isMajor ? 8 : 4) * dpr, y);
        ctx.lineTo(pad.l, y);
        ctx.stroke();
        if (isMajor) {
          ctx.textAlign = "right";
          ctx.fillText(formatPlainLogTick(10 ** tick.exponent), pad.l - 12 * dpr, y + 4 * dpr);
        }
      } else {
        const x = pad.l + ((value - logMin) / (logMax - logMin || 1)) * (canvas.width - pad.l - pad.r);
        ctx.strokeStyle = "#000";
        ctx.lineWidth = isMajor ? 1.4 * dpr : 1 * dpr;
        ctx.beginPath();
        ctx.moveTo(x, canvas.height - pad.b);
        ctx.lineTo(x, canvas.height - pad.b + (isMajor ? 8 : 4) * dpr);
        ctx.stroke();
        if (isMajor) {
          ctx.textAlign = "center";
          ctx.fillText(formatPlainLogTick(10 ** tick.exponent), x, canvas.height - pad.b + 30 * dpr);
        }
      }
  }
  ctx.lineWidth = 1 * dpr;
}

function formatPlainLogTick(value) {
  if (value >= 1) return formatTick(value).replace(/\.0+$/, "");
  return value.toFixed(Math.max(0, Math.ceil(-Math.log10(value)))).replace(/0+$/, "").replace(/\.$/, "");
}

function formatTick(value) {
  if (Math.abs(value) >= 100) return value.toFixed(0);
  if (Math.abs(value) >= 10) return value.toFixed(1);
  if (Math.abs(value) >= 1) return value.toFixed(2);
  if (Math.abs(value) >= 0.01) return value.toFixed(2);
  return value.toPrecision(2);
}

function formatLinearAxisTick(value) {
  if (Math.abs(value - Math.round(value)) < 1e-9) return String(Math.round(value));
  return value.toFixed(2).replace(/0+$/, "").replace(/\.$/, "");
}

function logRange(values, floorValue = 1e-4) {
  const filtered = values.filter((value) => value > 0);
  const minValue = Math.max(floorValue, filtered.length ? Math.min(...filtered) : floorValue);
  const maxValue = Math.max(floorValue, filtered.length ? Math.max(...filtered) : floorValue * 10);
  return {
    min: Math.max(floorValue, 10 ** Math.floor(Math.log10(minValue))),
    max: 10 ** Math.ceil(Math.log10(maxValue)),
  };
}

function linearRange(values, { includeZeroPadding = false } = {}) {
  const maxValue = values.length ? Math.max(...values) : 1;
  const minValue = values.length ? Math.min(...values) : 0;
  const max = Math.ceil(maxValue);
  const min = includeZeroPadding && minValue <= 0 ? -0.1 : Math.min(0, Math.floor(minValue));
  return { min, max: max === min ? min + 1 : max };
}

function drawHistogram() {
  if (!state.data) return;
  const componentValues = state.data.blocks.flatMap((block) => block.values || []);
  const positive = componentValues.filter((value) => value > 0);
  const negative = componentValues.filter((value) => value < 0).map(Math.abs);
  const zeroCount = componentValues.filter((value) => value === 0).length;
  const belowThreshold = componentValues.filter((value) => Math.abs(value) > 0 && Math.abs(value) < COMPONENT_PLOT_THRESHOLD).length;
  const positiveBelowThreshold = positive.filter((value) => value < COMPONENT_PLOT_THRESHOLD).length;
  const negativeBelowThreshold = negative.filter((value) => value < COMPONENT_PLOT_THRESHOLD).length;
  const positivePlotted = positive.filter((value) => value >= COMPONENT_PLOT_THRESHOLD);
  const negativePlotted = negative.filter((value) => value >= COMPONENT_PLOT_THRESHOLD);
  const symmetry = componentSymmetryStats(positive, negative);
  const plottedSentence = symmetry.matches
    ? `Only the positive component distribution is plotted.`
    : `<span class="warning">The positive and negative component magnitude distributions are not fully paired for this file.</span> Positive and negative component magnitude distributions are both plotted.`;
  $("componentHistogramNote").innerHTML = `The dataset contains ${positive.length} positive, ${negative.length} negative, and ${zeroCount} zero FC3 tensor components; ${belowThreshold} nonzero components fall below the plotting threshold of ${formatScientific(COMPONENT_PLOT_THRESHOLD)} eV/Å³. The sorted magnitudes of positive and negative components ${symmetry.matches ? "are paired within the numerical-noise tolerance used for this check" : "are not fully paired under the numerical-noise tolerance used for this check"}; the maximum absolute mismatch is ${formatScientific(symmetry.maxDiff)}. ${plottedSentence} Zero and below-threshold components are retained in the table but omitted from the logarithmic plot.`;
  if (!positivePlotted.length && !negativePlotted.length) return;
  const allPlotted = symmetry.matches ? positivePlotted : positivePlotted.concat(negativePlotted);
  const range = logRange(allPlotted, COMPONENT_PLOT_THRESHOLD);
  const positiveHistogram = logMinorHistogram(positivePlotted, range.min, range.max);
  const negativeHistogram = symmetry.matches ? null : logMinorHistogram(negativePlotted, range.min, range.max);
  renderComponentDistributionTable(
    {
      positive: positiveHistogram.rows,
      negative: negativeHistogram ? negativeHistogram.rows : [],
      positiveBelow: positiveBelowThreshold,
      negativeBelow: negativeBelowThreshold,
    },
    zeroCount,
    belowThreshold,
    symmetry.matches,
  );
  const canvas = $("histChart");
  const yRange = logRange((negativeHistogram ? positiveHistogram.counts.concat(negativeHistogram.counts) : positiveHistogram.counts).filter((count) => count > 0), 1);
  const frame = logChartFrame(canvas, "FC3_component (eV/Å³)", "count", range.min, range.max, yRange.min, yRange.max);
  drawDistributionCurve(frame, positiveHistogram.points, "#005eb8", false);
  const legendItems = [{ label: symmetry.matches ? "component" : "positive component", color: "#005eb8", dashed: false }];
  if (negativeHistogram) {
    drawDistributionCurve(frame, negativeHistogram.points, "#c62828", false);
    legendItems.push({ label: "negative component", color: "#c62828", dashed: false });
  }
  drawSimpleLegend(frame, legendItems, { boxed: false, corner: "right", vertical: !symmetry.matches });
}

function logHistogramCounts(values, bins, logMin, logMax) {
  const counts = Array.from({ length: bins }, () => 0);
  values.forEach((value) => {
    if (value <= 0) return;
    const idx = Math.min(bins - 1, Math.max(0, Math.floor(((Math.log10(value) - logMin) / (logMax - logMin || 1)) * bins)));
    counts[idx] += 1;
  });
  return counts;
}

function logMinorHistogram(values, min, max) {
  const minExp = Math.floor(Math.log10(min));
  const maxExp = Math.ceil(Math.log10(max));
  const edges = [];
  for (let exp = minExp; exp <= maxExp; exp += 1) {
    for (let multiplier = 1; multiplier <= 9; multiplier += 1) {
      const edge = multiplier * 10 ** exp;
      if (edge >= min * 0.999999 && edge <= max * 1.000001) edges.push(edge);
    }
  }
  if (edges[0] > min) edges.unshift(min);
  if (edges[edges.length - 1] < max) edges.push(max);
  const uniqueEdges = [...new Set(edges)].sort((a, b) => a - b);
  if (uniqueEdges.length < 2) uniqueEdges.push(uniqueEdges[0] * 1.01);
  const counts = Array.from({ length: uniqueEdges.length - 1 }, () => 0);
  values.forEach((value) => {
    let idx = uniqueEdges.findIndex((edge, edgeIdx) => edgeIdx < uniqueEdges.length - 1 && value >= edge && value < uniqueEdges[edgeIdx + 1]);
    if (idx === -1 && value === uniqueEdges[uniqueEdges.length - 1]) idx = counts.length - 1;
    if (idx >= 0) counts[idx] += 1;
  });
  const rows = counts.map((count, idx) => {
    const lower = uniqueEdges[idx];
    const upper = uniqueEdges[idx + 1];
    return { lower, upper, center: Math.sqrt(lower * upper), count };
  });
  return {
    min: uniqueEdges[0],
    max: uniqueEdges[uniqueEdges.length - 1],
    counts,
    rows,
    points: rows.map((row, idx) => ({ ...row, idx })).filter((row) => row.count > 0),
  };
}

function componentSymmetryStats(positiveAbs, negativeAbs) {
  const pos = [...positiveAbs].sort((a, b) => a - b);
  const neg = [...negativeAbs].sort((a, b) => a - b);
  const n = Math.min(pos.length, neg.length);
  let maxDiff = pos.length === neg.length ? 0 : Infinity;
  let matches = pos.length === neg.length;
  for (let idx = 0; idx < n; idx += 1) {
    const diff = Math.abs(pos[idx] - neg[idx]);
    maxDiff = Math.max(maxDiff, diff);
    if (diff > SYMMETRY_ATOL + SYMMETRY_RTOL * Math.abs(neg[idx])) matches = false;
  }
  return { matches, maxDiff };
}

function drawDistributionCurve(frame, points, color, dashed = false) {
  const ctx = frame.ctx;
  ctx.save();
  ctx.strokeStyle = color;
  ctx.fillStyle = color;
  ctx.lineWidth = 2.3 * frame.dpr;
  for (let idx = 1; idx < points.length; idx += 1) {
    const prev = points[idx - 1];
    const current = points[idx];
    ctx.setLineDash(dashed || current.idx !== prev.idx + 1 ? [7 * frame.dpr, 5 * frame.dpr] : []);
    ctx.beginPath();
    ctx.moveTo(frame.x(prev.center), frame.y(prev.count));
    ctx.lineTo(frame.x(current.center), frame.y(current.count));
    ctx.stroke();
  }
  ctx.setLineDash([]);
  points.forEach((point) => {
    ctx.beginPath();
    ctx.arc(frame.x(point.center), frame.y(point.count), 3 * frame.dpr, 0, Math.PI * 2);
    ctx.fill();
  });
  ctx.restore();
}

function drawHistogramCurve(frame, counts, bins, logMin, logMax, color, dashed = false) {
  const ctx = frame.ctx;
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = 2.2 * frame.dpr;
  ctx.setLineDash(dashed ? [8 * frame.dpr, 5 * frame.dpr] : []);
  ctx.beginPath();
  let drawing = false;
  let hasPoint = false;
  counts.forEach((count, idx) => {
    if (count <= 0) {
      drawing = false;
      return;
    }
    const xValue = 10 ** (logMin + ((logMax - logMin) * (idx + 0.5)) / bins);
    const x = frame.x(xValue);
    const y = frame.y(count);
    if (!drawing) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
    drawing = true;
    hasPoint = true;
  });
  if (hasPoint) ctx.stroke();
  ctx.restore();
}

function drawSimpleLegend(frame, items, options = {}) {
  const { ctx, dpr, pad } = frame;
  const boxed = Boolean(options.boxed);
  const estimatedWidth = Math.max(...items.map((item) => item.label.length)) * 7.6 * dpr + 60 * dpr;
  const rowHeight = 20 * dpr;
  const boxHeight = (options.vertical ? items.length * rowHeight : rowHeight) + 12 * dpr;
  let x = options.corner === "right" ? frame.ctx.canvas.width - pad.r - estimatedWidth - 10 * dpr : pad.l + 8 * dpr;
  const y = options.y ? pad.t + options.y * dpr : pad.t + 15 * dpr;
  ctx.save();
  ctx.font = `${13 * dpr}px Arial`;
  ctx.textAlign = "left";
  if (boxed) {
    ctx.fillStyle = "rgba(255,255,255,0.88)";
    ctx.strokeStyle = "#000";
    ctx.lineWidth = 0.9 * dpr;
    ctx.fillRect(x - 8 * dpr, y - 20 * dpr, estimatedWidth + 12 * dpr, boxHeight);
    ctx.strokeRect(x - 8 * dpr, y - 20 * dpr, estimatedWidth + 12 * dpr, boxHeight);
  }
  items.forEach((item, idx) => {
    const itemX = options.vertical ? x : x;
    const itemY = options.vertical ? y + idx * rowHeight : y;
    ctx.strokeStyle = item.color;
    ctx.fillStyle = item.color;
    ctx.lineWidth = 2.2 * dpr;
    ctx.setLineDash(item.dashed ? [8 * dpr, 5 * dpr] : []);
    if (item.marker === "point") {
      ctx.beginPath();
      ctx.arc(itemX + 12 * dpr, itemY - 4 * dpr, 4 * dpr, 0, Math.PI * 2);
      ctx.fill();
    } else {
      ctx.beginPath();
      ctx.moveTo(itemX, itemY - 4 * dpr);
      ctx.lineTo(itemX + 24 * dpr, itemY - 4 * dpr);
      ctx.stroke();
    }
    ctx.setLineDash([]);
    ctx.fillStyle = "#000";
    ctx.fillText(item.label, itemX + 30 * dpr, itemY);
    if (!options.vertical) x += (item.label.length * 7.6 + 62) * dpr;
  });
  ctx.restore();
}

function drawShellChart() {
  if (!state.data) return;
  const data = state.data.blocks.filter((block) => block.norm > 0);
  const plotted = data.filter((block) => block.norm >= BLOCK_PLOT_THRESHOLD);
  const belowThreshold = data.length - plotted.length;
  const norms = data.map((block) => block.norm);
  const perimeters = plotted.map((block) => block.perimeter);
  $("perimeterChartNote").textContent = `The dataset contains ${data.length} FC3 blocks; ${plotted.length} blocks are plotted and ${belowThreshold} blocks fall below the plotting threshold of ${formatScientific(BLOCK_PLOT_THRESHOLD)} eV/Å³. FC3_block ranges from ${formatScientific(Math.min(...norms))} to ${formatScientific(Math.max(...norms))} eV/Å³, and perimeter ranges from ${formatNumber(Math.min(...data.map((block) => block.perimeter)))} to ${formatNumber(Math.max(...data.map((block) => block.perimeter)))} Å.`;
  renderBlockPerimeterTable(data);
  const canvas = $("shellChart");
  const xRange = linearRange(perimeters, { includeZeroPadding: plotted.some((block) => block.perimeter === 0) });
  const yRange = logRange(plotted.map((block) => block.norm), BLOCK_PLOT_THRESHOLD);
  const frame = linearLogYChartFrame(
    canvas,
    "perimeter (Å)",
    "FC3_block (eV/Å³)",
    xRange.min,
    xRange.max,
    yRange.min,
    yRange.max,
    { linearXTickMin: 0, formatX: formatLinearAxisTick },
  );
  plotted.forEach((item) => {
    frame.ctx.fillStyle = "#c62828";
    frame.ctx.beginPath();
    frame.ctx.arc(frame.x(item.perimeter), frame.y(item.norm), 2.7 * frame.dpr, 0, Math.PI * 2);
    frame.ctx.fill();
  });
  drawSimpleLegend(frame, [{ label: "FC3 blocks", color: "#c62828", marker: "point" }], { boxed: true, corner: "right", y: 34 });
}

function drawBlockDistributionChart() {
  if (!state.data) return;
  const values = state.data.blocks.map((block) => block.norm).filter((value) => value > 0);
  const belowThreshold = values.filter((value) => value < BLOCK_PLOT_THRESHOLD).length;
  const plotted = values.filter((value) => value >= BLOCK_PLOT_THRESHOLD);
  $("blockDistributionNote").textContent = `The dataset contains ${values.length} FC3 blocks; ${belowThreshold} blocks fall below the plotting threshold of ${formatScientific(BLOCK_PLOT_THRESHOLD)} eV/Å³ and are retained in the table.`;
  if (!plotted.length) return;
  const range = logRange(plotted, BLOCK_PLOT_THRESHOLD);
  const histogram = logMinorHistogram(plotted, range.min, range.max);
  renderBlockDistributionTable(histogram.rows, belowThreshold);
  const yRange = logRange(histogram.counts.filter((count) => count > 0), 1);
  const frame = logChartFrame($("blockDistChart"), "FC3_block (eV/Å³)", "count", range.min, range.max, yRange.min, yRange.max);
  drawDistributionCurve(frame, histogram.points, "#005eb8", false);
  drawSimpleLegend(frame, [{ label: "FC3_block", color: "#005eb8" }], { corner: "right" });
}

function drawMaxEdgeChart() {
  if (!state.data) return;
  const data = state.data.blocks.filter((block) => block.norm > 0);
  const plotted = data.filter((block) => block.norm >= BLOCK_PLOT_THRESHOLD);
  const belowThreshold = data.length - plotted.length;
  const maxEdges = data.map((block) => blockEdgeLengths(block)[2]);
  $("maxEdgeChartNote").textContent = `The dataset contains ${data.length} FC3 blocks; ${plotted.length} blocks are plotted and ${belowThreshold} blocks fall below the plotting threshold of ${formatScientific(BLOCK_PLOT_THRESHOLD)} eV/Å³. max_edge_length ranges from ${formatNumber(Math.min(...maxEdges))} to ${formatNumber(Math.max(...maxEdges))} Å.`;
  renderBlockMaxEdgeTable(data);
  if (!plotted.length) return;
  const xValues = plotted.map((block) => blockEdgeLengths(block)[2]);
  const xRange = linearRange(xValues, { includeZeroPadding: xValues.some((value) => value === 0) });
  const yRange = logRange(plotted.map((block) => block.norm), BLOCK_PLOT_THRESHOLD);
  const frame = linearLogYChartFrame($("maxEdgeChart"), "max_edge_length (Å)", "FC3_block (eV/Å³)", xRange.min, xRange.max, yRange.min, yRange.max, { linearXTickMin: 0, formatX: formatLinearAxisTick });
  plotted.forEach((block) => {
    frame.ctx.fillStyle = "#c62828";
    frame.ctx.beginPath();
    frame.ctx.arc(frame.x(blockEdgeLengths(block)[2]), frame.y(block.norm), 2.7 * frame.dpr, 0, Math.PI * 2);
    frame.ctx.fill();
  });
  drawSimpleLegend(frame, [{ label: "FC3 blocks", color: "#c62828", marker: "point" }], { boxed: true, corner: "right", y: 34 });
}

function renderComponentDistributionTable(histograms, zeroCount, belowThreshold, symmetric) {
  if (symmetric) {
    const tableRows = [
      [0, 0, 0, zeroCount],
      [0, COMPONENT_PLOT_THRESHOLD, COMPONENT_PLOT_THRESHOLD / 2, belowThreshold],
      ...histograms.positive.filter((row) => row.count > 0).map((row) => [row.lower, row.upper, row.center, row.count]),
    ].map((row) => [formatNumber(row[0]), formatNumber(row[1]), formatNumber(row[2]), String(row[3])]);
    renderTable("componentDistributionTable", ["bin_lower (eV/Å³)", "bin_upper (eV/Å³)", "bin_center (eV/Å³)", "count"], tableRows);
    return;
  }
  const signedRows = [
    ["zero", 0, 0, 0, zeroCount],
    ["positive below threshold", 0, COMPONENT_PLOT_THRESHOLD, COMPONENT_PLOT_THRESHOLD / 2, histograms.positiveBelow || 0],
    ["negative below threshold", 0, COMPONENT_PLOT_THRESHOLD, COMPONENT_PLOT_THRESHOLD / 2, histograms.negativeBelow || 0],
    ...histograms.positive.filter((row) => row.count > 0).map((row) => ["positive", row.lower, row.upper, row.center, row.count]),
    ...histograms.negative.filter((row) => row.count > 0).map((row) => ["negative", row.lower, row.upper, row.center, row.count]),
  ].map((row) => [row[0], formatNumber(row[1]), formatNumber(row[2]), formatNumber(row[3]), String(row[4])]);
  renderTable("componentDistributionTable", ["component_sign", "bin_lower (eV/Å³)", "bin_upper (eV/Å³)", "bin_center (eV/Å³)", "count"], signedRows);
}

function renderBlockPerimeterTable(blocks) {
  const groups = new Map();
  blocks.forEach((block) => {
    const key = Number(block.norm.toPrecision(10)).toString();
    const existing = groups.get(key);
    if (!existing) {
      groups.set(key, { norm: Number(block.norm.toPrecision(10)), count: 1, representative: block });
      return;
    }
    existing.count += 1;
    if (block.id < existing.representative.id) existing.representative = block;
  });
  const rows = [...groups.values()]
    .sort((a, b) => b.norm - a.norm)
    .map((group) => {
      const edges = blockEdgeLengths(group.representative);
      return [
        formatScientific(group.norm),
        String(group.count),
        String(group.representative.id),
        formatNumber(edges[0]),
        formatNumber(edges[1]),
        formatNumber(edges[2]),
        formatNumber(group.representative.perimeter),
        group.norm >= BLOCK_PLOT_THRESHOLD ? "yes" : "no",
      ];
    });
  renderTable("blockPerimeterTable", ["FC3_block (eV/Å³)", "block_count", "min_block_id", "edge_1 (Å)", "edge_2 (Å)", "edge_3 (Å)", "perimeter (Å)", "plotted"], rows);
}

function renderBlockDistributionTable(rows, belowThreshold) {
  const tableRows = [
    [0, BLOCK_PLOT_THRESHOLD, BLOCK_PLOT_THRESHOLD / 2, belowThreshold],
    ...rows.filter((row) => row.count > 0).map((row) => [row.lower, row.upper, row.center, row.count]),
  ].map((row) => [formatNumber(row[0]), formatNumber(row[1]), formatNumber(row[2]), String(row[3])]);
  renderTable("blockDistributionTable", ["bin_lower (eV/Å³)", "bin_upper (eV/Å³)", "bin_center (eV/Å³)", "count"], tableRows);
}


function blockEdgeLengths(block) {
  const [a, b, c] = block.vertices;
  return [length(sub(a, b)), length(sub(a, c)), length(sub(b, c))].sort((x, y) => x - y);
}

function renderBlockMaxEdgeTable(blocks) {
  const groups = new Map();
  blocks.forEach((block) => {
    const key = Number(block.norm.toPrecision(10)).toString();
    const existing = groups.get(key);
    if (!existing) {
      groups.set(key, { norm: Number(block.norm.toPrecision(10)), count: 1, representative: block });
      return;
    }
    existing.count += 1;
    if (block.id < existing.representative.id) existing.representative = block;
  });
  const rows = [...groups.values()]
    .sort((a, b) => b.norm - a.norm)
    .map((group) => {
      const edges = blockEdgeLengths(group.representative);
      return [
        formatScientific(group.norm),
        String(group.count),
        String(group.representative.id),
        formatNumber(edges[2]),
        group.norm >= BLOCK_PLOT_THRESHOLD ? "yes" : "no",
      ];
    });
  renderTable("blockMaxEdgeTable", ["FC3_block (eV/Å³)", "block_count", "min_block_id", "max_edge_length (Å)", "plotted"], rows);
}

function renderTable(id, columns, rows) {
  const table = $(id);
  if (!table) return;
  const header = `<thead><tr>${columns.map((column) => `<th>${column}</th>`).join("")}</tr></thead>`;
  const body = `<tbody>${rows.map((row) => `<tr>${row.map((cell) => `<td>${cell}</td>`).join("")}</tr>`).join("")}</tbody>`;
  table.innerHTML = header + body;
}

function formatNumber(value) {
  if (!Number.isFinite(value)) return "n/a";
  if (value === 0) return "0";
  const abs = Math.abs(value);
  if (abs >= 1000 || abs < 0.001) return value.toExponential(4);
  return Number(value.toPrecision(6)).toString();
}

function formatScientific(value) {
  if (!Number.isFinite(value)) return "n/a";
  if (value === 0) return "0";
  return value.toExponential(4);
}

function componentValues(dataset) {
  return dataset.blocks.flatMap((block) => block.values || []);
}

function componentAbsNonzero(dataset) {
  return componentValues(dataset).map(Math.abs).filter((value) => value > 0);
}

function signedComponentAbs(dataset, sign) {
  return componentValues(dataset)
    .filter((value) => (sign > 0 ? value > 0 : value < 0))
    .map(Math.abs);
}

function compareDatasets() {
  const entries = state.compareSlots
    .map((slot, idx) => ({
      dataset: slot.dataset,
      legend: slot.legend || slot.dataset?.meta.fc3File || `FC3 ${idx + 1}`,
      color: slot.color || COMPARE_COLORS[idx % COMPARE_COLORS.length],
    }))
    .filter((entry) => entry.dataset);
  return entries;
}

function palette(index) {
  return ["#005eb8", "#c62828", "#2e7d32", "#6a1b9a", "#ef6c00", "#00838f"][index % 6];
}

function drawLegend(frame, datasets) {
  const { ctx, dpr, pad } = frame;
  let x = pad.l + 8 * dpr;
  const y = pad.t + 10 * dpr;
  ctx.font = `${11 * dpr}px Arial`;
  ctx.textAlign = "left";
  datasets.forEach((dataset, idx) => {
    ctx.fillStyle = palette(idx);
    ctx.fillRect(x, y - 8 * dpr, 16 * dpr, 3 * dpr);
    ctx.fillStyle = "#000";
    const label = dataset.meta.fc3File;
    ctx.fillText(label, x + 20 * dpr, y);
    x += (label.length * 6.3 + 36) * dpr;
  });
}

function drawCompareComponentChart() {
  const entries = compareDatasets();
  if (!entries.length) {
    clearChartCanvas($("compareComponentChart"));
    clearChartCanvas($("compareNegativeComponentChart"));
    $("compareNegativeFigure").hidden = true;
    return;
  }
  const analyses = entries.map((entry) => {
    const values = componentValues(entry.dataset);
    const positive = values.filter((value) => value > 0);
    const negative = values.filter((value) => value < 0).map(Math.abs);
    return {
      ...entry,
      positive: positive.filter((value) => value >= COMPONENT_PLOT_THRESHOLD),
      negative: negative.filter((value) => value >= COMPONENT_PLOT_THRESHOLD),
      symmetric: componentSymmetryStats(positive, negative).matches,
    };
  });
  const showNegative = analyses.some((entry) => !entry.symmetric);
  $("compareNegativeFigure").hidden = !showNegative;
  $("compareComponentCaption").textContent = showNegative ? "positive FC3_component count distribution" : "FC3_component count distribution";
  const allValues = analyses.flatMap((entry) => showNegative ? entry.positive.concat(entry.negative) : entry.positive);
  if (!allValues.length) return;
  const xRange = logRange(allValues, COMPONENT_PLOT_THRESHOLD);
  const positiveHistograms = analyses.map((entry) => logMinorHistogram(entry.positive, xRange.min, xRange.max));
  const negativeHistograms = analyses.map((entry) => logMinorHistogram(entry.negative, xRange.min, xRange.max));
  const yValues = showNegative
    ? positiveHistograms.flatMap((hist) => hist.counts).concat(negativeHistograms.flatMap((hist) => hist.counts))
    : positiveHistograms.flatMap((hist) => hist.counts);
  const yRange = logRange(yValues.filter((count) => count > 0), 1);
  const drawComponentCompare = (canvas, histograms) => {
    const frame = logChartFrame(canvas, "FC3_component (eV/Å³)", "count", xRange.min, xRange.max, yRange.min, yRange.max);
    histograms.forEach((histogram, idx) => {
      drawDistributionCurve(frame, histogram.points, analyses[idx].color, false);
    });
    drawSimpleLegend(frame, analyses.map((entry) => ({ label: entry.legend, color: entry.color, dashed: false })), { corner: "right", vertical: true });
  };
  drawComponentCompare($("compareComponentChart"), positiveHistograms);
  if (showNegative) drawComponentCompare($("compareNegativeComponentChart"), negativeHistograms);
}

function smoothCounts(counts) {
  return counts.map((count, idx) => {
    const prev = counts[idx - 1] || 0;
    const next = counts[idx + 1] || 0;
    if (prev === 0 && count === 0 && next === 0) return 0;
    return (prev + 2 * count + next) / 4;
  });
}

function drawComparePerimeterChart() {
  const entries = compareDatasets();
  const plottedByEntry = entries.map((entry) => ({
    ...entry,
    blocks: entry.dataset.blocks.filter((block) => block.norm >= BLOCK_PLOT_THRESHOLD),
  }));
  const allBlocks = plottedByEntry.flatMap((entry) => entry.blocks);
  if (!allBlocks.length) {
    clearChartCanvas($("comparePerimeterChart"));
    return;
  }
  const xRange = linearRange(allBlocks.map((block) => block.perimeter), { includeZeroPadding: allBlocks.some((block) => block.perimeter === 0) });
  const yRange = logRange(allBlocks.map((block) => block.norm), BLOCK_PLOT_THRESHOLD);
  const frame = linearLogYChartFrame(
    $("comparePerimeterChart"),
    "perimeter (Å)",
    "FC3_block (eV/Å³)",
    xRange.min,
    xRange.max,
    yRange.min,
    yRange.max,
    { linearXTickMin: 0, formatX: formatLinearAxisTick },
  );
  plottedByEntry.forEach((entry) => {
    frame.ctx.fillStyle = entry.color;
    entry.blocks.forEach((block) => {
      frame.ctx.beginPath();
      frame.ctx.arc(frame.x(block.perimeter), frame.y(block.norm), 2.4 * frame.dpr, 0, Math.PI * 2);
      frame.ctx.fill();
    });
  });
  drawSimpleLegend(frame, plottedByEntry.map((entry) => ({ label: entry.legend, color: entry.color, marker: "point" })), { boxed: true, corner: "right", vertical: true, y: 34 });
}

function drawCompareBlockDistributionChart() {
  const entries = compareDatasets();
  const valuesByEntry = entries.map((entry) => ({ ...entry, values: entry.dataset.blocks.map((block) => block.norm).filter((value) => value >= BLOCK_PLOT_THRESHOLD) }));
  const allValues = valuesByEntry.flatMap((entry) => entry.values);
  if (!allValues.length) return clearChartCanvas($("compareBlockDistChart"));
  const xRange = logRange(allValues, BLOCK_PLOT_THRESHOLD);
  const histograms = valuesByEntry.map((entry) => ({ ...entry, histogram: logMinorHistogram(entry.values, xRange.min, xRange.max) }));
  const yRange = logRange(histograms.flatMap((entry) => entry.histogram.counts).filter((count) => count > 0), 1);
  const frame = logChartFrame($("compareBlockDistChart"), "FC3_block (eV/Å³)", "count", xRange.min, xRange.max, yRange.min, yRange.max);
  histograms.forEach((entry) => drawDistributionCurve(frame, entry.histogram.points, entry.color, false));
  drawSimpleLegend(frame, histograms.map((entry) => ({ label: entry.legend, color: entry.color })), { corner: "right", vertical: true });
}

function drawCompareMaxEdgeChart() {
  const entries = compareDatasets();
  const plottedByEntry = entries.map((entry) => ({ ...entry, blocks: entry.dataset.blocks.filter((block) => block.norm >= BLOCK_PLOT_THRESHOLD) }));
  const allBlocks = plottedByEntry.flatMap((entry) => entry.blocks);
  if (!allBlocks.length) return clearChartCanvas($("compareMaxEdgeChart"));
  const xValues = allBlocks.map((block) => blockEdgeLengths(block)[2]);
  const xRange = linearRange(xValues, { includeZeroPadding: xValues.some((value) => value === 0) });
  const yRange = logRange(allBlocks.map((block) => block.norm), BLOCK_PLOT_THRESHOLD);
  const frame = linearLogYChartFrame($("compareMaxEdgeChart"), "max_edge_length (Å)", "FC3_block (eV/Å³)", xRange.min, xRange.max, yRange.min, yRange.max, { linearXTickMin: 0, formatX: formatLinearAxisTick });
  plottedByEntry.forEach((entry) => {
    frame.ctx.fillStyle = entry.color;
    entry.blocks.forEach((block) => {
      frame.ctx.beginPath();
      frame.ctx.arc(frame.x(blockEdgeLengths(block)[2]), frame.y(block.norm), 2.4 * frame.dpr, 0, Math.PI * 2);
      frame.ctx.fill();
    });
  });
  drawSimpleLegend(frame, plottedByEntry.map((entry) => ({ label: entry.legend, color: entry.color, marker: "point" })), { boxed: true, corner: "right", vertical: true, y: 34 });
}

function clearChartCanvas(canvas) {
  const dpr = prepareChartCanvas(canvas);
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = "#fff";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = "#5f6368";
  ctx.font = `${14 * dpr}px Arial`;
  ctx.textAlign = "center";
  ctx.fillText("No compare data loaded.", canvas.width / 2, canvas.height / 2);
}

function drawCompareCharts() {
  drawCompareComponentChart();
  drawComparePerimeterChart();
  drawCompareBlockDistributionChart();
  drawCompareMaxEdgeChart();
}

function drawTopChart() {
  if (!state.data) return;
  const data = state.data.topBlocks.slice(0, 20);
  const canvas = $("topChart");
  const frame = chartFrame(canvas, "Top block rank", "||Phi3_block||F", 0, data.length, 0, Math.max(...data.map((d) => d.norm)));
  if (!frame) return;
  const barWidth = (canvas.width - frame.pad.l - frame.pad.r) / data.length - 3 * frame.dpr;
  data.forEach((item, idx) => {
    const x = frame.x(idx);
    const y = frame.y(item.norm);
    frame.ctx.fillStyle = idx < 2 ? "#c62828" : "#005eb8";
    frame.ctx.fillRect(x + 2 * frame.dpr, y, barWidth, canvas.height - frame.pad.b - y);
  });
}

function drawScatterChart() {
  if (!state.data) return;
  const data = state.data.blocks;
  const canvas = $("scatterChart");
  const yMin = Math.floor(Math.log10(Math.max(state.data.meta.normMin, 1e-12)));
  const yMax = Math.ceil(Math.log10(state.data.meta.normMax));
  const frame = chartFrame(canvas, "dmax (angstrom)", "log10 ||Phi3_block||F", state.data.meta.dmaxMin, state.data.meta.dmaxMax, yMin, yMax);
  if (!frame) return;
  data.forEach((item) => {
    frame.ctx.fillStyle = "rgba(0, 94, 184, 0.72)";
    frame.ctx.beginPath();
    frame.ctx.arc(frame.x(item.dmax), frame.y(Math.log10(Math.max(item.norm, 1e-12))), 2.7 * frame.dpr, 0, Math.PI * 2);
    frame.ctx.fill();
  });
}

function drawAll() {
  if (!state.data) return;
  drawScene();
  drawHistogram();
  drawBlockDistributionChart();
  drawShellChart();
  drawMaxEdgeChart();
  drawCompareCharts();
}

async function main() {
  setupTabs();
  setupIntroTabs();
  setupSceneEvents();
  setupUploadEvents();
  setupCompareEvents();
  const response = await fetch("./data/fc3_045.json");
  state.defaultData = await response.json();
  state.data = structuredClone(state.defaultData);
  state.structure = structureFromPayload(state.data);
  initializeNormRange();
  initializeBoundaryFromFc3();
  $("dataStatus").textContent = `Default example loaded: FORCE_CONSTANTS_3RD_Example and BASE.si_supper.scf.in_example, ${state.data.meta.blockCount} blocks.`;
  $("fc3FileDisplay").textContent = "FORCE_CONSTANTS_3RD_Example";
  $("qeFileDisplay").textContent = "BASE.si_supper.scf.in_example";
  buildObjectTree();
  updateCompareStatus();
  activateTab(location.hash.slice(1));
  drawAll();
  window.addEventListener("resize", drawAll);
}

main().catch((error) => {
  document.body.innerHTML = `<pre style="padding:24px;color:#b00020">${error.stack || error}</pre>`;
});


