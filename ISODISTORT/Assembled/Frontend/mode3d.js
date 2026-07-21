import * as THREE from "three";
import { TrackballControls } from "three/addons/controls/TrackballControls.js";
import { ConvexGeometry } from "three/addons/geometries/ConvexGeometry.js";
import {
  addLights,
  atomColor,
  atomRadius,
  loadElementStyles,
  makeCellEdges,
} from "./crystal3d.js";

const EPS = 1e-9;
// VESTA models coordination with an A1-A2-specific distance interval. The
// automatic viewer has no user-authored bond table, so widen each pair's
// shortest contact enough to retain a moderately distorted first shell while
// capping the absolute expansion to avoid absorbing a distant second shell.
const COORDINATION_RELATIVE_TOLERANCE = 0.2;
const COORDINATION_ABSOLUTE_TOLERANCE = 0.45;
const REFERENCE_MOMENT = Math.sqrt(3 * 0.569 ** 2);
const REFERENCE_ARROW_LENGTH = 2.72 / 4.17;
const REFERENCE_SHAFT_RADIUS = 0.12 / 4.17;
const REFERENCE_HEAD_RADIUS = 0.33 / 4.17;
const REFERENCE_HEAD_LENGTH = 0.52 / 4.17;
const CUBE_PX = 116;
const CUBE_MARGIN = 6;
const CUBE_VIEW = 1.2;
const Y_AXIS = new THREE.Vector3(0, 1, 0);
const UP_Z = new THREE.Vector3(0, 0, 1);
const UP_Y = new THREE.Vector3(0, 1, 0);
const SHARED_RESOURCE = "modeViewerShared";
const atomGeometryCache = new Map();
const atomMaterialCache = new Map();
const bondMaterialCache = new Map();
const unitCylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, 12);
unitCylinderGeometry.userData[SHARED_RESOURCE] = true;
const unitConeGeometry = new THREE.ConeGeometry(1, 1, 20);
unitConeGeometry.userData[SHARED_RESOURCE] = true;
const magneticArrowMaterial = new THREE.MeshStandardMaterial({
  color: 0xf00000,
  roughness: 0.25,
  metalness: 0.22,
});
magneticArrowMaterial.userData[SHARED_RESOURCE] = true;
const MODE_KINDS = [
  ["displacive", "Atomic displacement", "displacive_definitions"],
  ["magnetic", "Magnetic moment", "magnetic_definitions"],
  ["strain", "Strain", "strain_definitions"],
];

function number(value, fallback = 0) {
  if (typeof value === "number") return Number.isFinite(value) ? value : fallback;
  const text = String(value ?? "").trim();
  if (!text) return fallback;
  if (text.includes("/")) {
    const [top, bottom] = text.split("/", 2).map(Number);
    if (Number.isFinite(top) && Number.isFinite(bottom) && Math.abs(bottom) > EPS) return top / bottom;
  }
  const parsed = Number(text);
  return Number.isFinite(parsed) ? parsed : fallback;
}

function vector3(value) {
  if (!Array.isArray(value) || value.length < 3) return [0, 0, 0];
  return [number(value[0]), number(value[1]), number(value[2])];
}

function matrix3(value) {
  if (!Array.isArray(value)) return null;
  if (value.length === 9) {
    return [vector3(value.slice(0, 3)), vector3(value.slice(3, 6)), vector3(value.slice(6, 9))];
  }
  if (value.length >= 3 && value.slice(0, 3).every(Array.isArray)) {
    return value.slice(0, 3).map(vector3);
  }
  if (value.length === 6) {
    const [e11, e22, e33, e23, e13, e12] = value.map(item => number(item));
    return [[e11, e12, e13], [e12, e22, e23], [e13, e23, e33]];
  }
  return null;
}

function addVector(target, source, scale = 1) {
  target[0] += source[0] * scale;
  target[1] += source[1] * scale;
  target[2] += source[2] * scale;
}

function addMatrix(target, source, scale = 1) {
  for (let row = 0; row < 3; row += 1) {
    for (let col = 0; col < 3; col += 1) target[row][col] += source[row][col] * scale;
  }
}

function multiplyMatrixVector(matrix, value) {
  return matrix.map(row => row[0] * value[0] + row[1] * value[1] + row[2] * value[2]);
}

function cellVectors(cell = {}) {
  const a = Math.max(number(cell.a, 1), EPS);
  const b = Math.max(number(cell.b, 1), EPS);
  const c = Math.max(number(cell.c, 1), EPS);
  const alpha = THREE.MathUtils.degToRad(number(cell.alpha, 90));
  const beta = THREE.MathUtils.degToRad(number(cell.beta, 90));
  const gamma = THREE.MathUtils.degToRad(number(cell.gamma, 90));
  const sinGamma = Math.abs(Math.sin(gamma)) > EPS ? Math.sin(gamma) : EPS;
  const av = new THREE.Vector3(a, 0, 0);
  const bv = new THREE.Vector3(b * Math.cos(gamma), b * Math.sin(gamma), 0);
  const cx = c * Math.cos(beta);
  const cy = c * (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma)) / sinGamma;
  const cz = Math.sqrt(Math.max(c * c - cx * cx - cy * cy, 0));
  return [av, bv, new THREE.Vector3(cx, cy, cz)];
}

function combine(vectors, fractional) {
  return vectors[0].clone().multiplyScalar(fractional[0])
    .add(vectors[1].clone().multiplyScalar(fractional[1]))
    .add(vectors[2].clone().multiplyScalar(fractional[2]));
}

function vectorsFromBasis(parentVectors, basis) {
  const rows = matrix3(basis);
  return rows ? rows.map(row => combine(parentVectors, row)) : null;
}

function fixed(value, digits = 6) {
  return Number(value).toFixed(digits).replace(/0+$/, "").replace(/\.$/, "");
}

function outsideShifts() {
  return [-1, 0, 1].flatMap(i => [-1, 0, 1].flatMap(j => [-1, 0, 1].map(k => [i, j, k])));
}

function isInsideUnit(value) {
  return value.every(component => component >= -1e-8 && component <= 1 + 1e-8);
}

function makeCylinderBetween(left, right, radius, color) {
  const delta = right.clone().sub(left);
  const length = delta.length();
  if (length < EPS) return null;
  const materialKey = Number(color);
  if (!bondMaterialCache.has(materialKey)) {
    const material = new THREE.MeshStandardMaterial({ color, roughness: 0.45, metalness: 0.02 });
    material.userData[SHARED_RESOURCE] = true;
    bondMaterialCache.set(materialKey, material);
  }
  const mesh = new THREE.Mesh(unitCylinderGeometry, bondMaterialCache.get(materialKey));
  mesh.scale.set(radius, length, radius);
  mesh.position.copy(left).add(right).multiplyScalar(0.5);
  mesh.quaternion.setFromUnitVectors(Y_AXIS, delta.normalize());
  return mesh;
}

function makeModeAtomMesh(atom, fallbackRadius) {
  const radius = atomRadius(atom, fallbackRadius);
  const geometryKey = fixed(radius, 5);
  if (!atomGeometryCache.has(geometryKey)) {
    const geometry = new THREE.SphereGeometry(radius, 20, 14);
    geometry.userData[SHARED_RESOURCE] = true;
    atomGeometryCache.set(geometryKey, geometry);
  }
  const color = atomColor(atom);
  if (!atomMaterialCache.has(color)) {
    const base = new THREE.Color(color);
    const material = new THREE.MeshStandardMaterial({
      color: base,
      roughness: 0.42,
      metalness: 0.05,
      emissive: base.clone().multiplyScalar(0.08),
    });
    material.userData[SHARED_RESOURCE] = true;
    atomMaterialCache.set(color, material);
  }
  return new THREE.Mesh(atomGeometryCache.get(geometryKey), atomMaterialCache.get(color));
}

function makeBondMesh(leftPosition, rightPosition, radius, leftAtom, rightAtom) {
  const midpoint = leftPosition.clone().add(rightPosition).multiplyScalar(0.5);
  const group = new THREE.Group();
  const left = makeCylinderBetween(leftPosition, midpoint, radius, atomColor(leftAtom));
  const right = makeCylinderBetween(midpoint, rightPosition, radius, atomColor(rightAtom));
  if (left) group.add(left);
  if (right) group.add(right);
  return group;
}

function makePolyhedronMesh(centerAtom, points) {
  const unique = [];
  const seen = new Set();
  for (const point of points) {
    const key = point.toArray().map(value => fixed(value, 5)).join(",");
    if (seen.has(key)) continue;
    seen.add(key);
    unique.push(point.clone());
  }
  if (unique.length < 4) return null;
  let geometry;
  try {
    geometry = new ConvexGeometry(unique);
  } catch {
    return null;
  }
  const color = atomColor(centerAtom);
  const group = new THREE.Group();
  group.add(new THREE.Mesh(
    geometry,
    new THREE.MeshStandardMaterial({
      color,
      transparent: true,
      opacity: 0.38,
      roughness: 0.55,
      metalness: 0.02,
      side: THREE.DoubleSide,
      depthWrite: false,
    }),
  ));
  group.add(new THREE.LineSegments(
    new THREE.EdgesGeometry(geometry),
    new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.82 }),
  ));
  return group;
}

function minimumDistance(atoms, vectors, leftElement, rightElement) {
  const leftAtoms = atoms.filter(atom => atom.element === leftElement);
  const rightAtoms = atoms.filter(atom => atom.element === rightElement);
  let best = Infinity;
  for (const left of leftAtoms) {
    const leftPosition = combine(vectors, left.frac);
    for (const right of rightAtoms) {
      for (const shift of outsideShifts()) {
        if (left.sourceKey === right.sourceKey && shift.every(value => value === 0)) continue;
        const shifted = right.frac.map((value, axis) => value + shift[axis]);
        const distance = leftPosition.distanceTo(combine(vectors, shifted));
        if (distance > 1e-8 && distance < best) best = distance;
      }
    }
  }
  return Number.isFinite(best) ? best : 2.4;
}

function coordinationMaximum(nearestDistance) {
  const tolerance = Math.min(
    nearestDistance * COORDINATION_RELATIVE_TOLERANCE,
    COORDINATION_ABSOLUTE_TOLERANCE,
  );
  return Math.ceil((nearestDistance + tolerance + 1e-12) * 1e6) / 1e6;
}

function defaultBondRules(atoms, vectors) {
  const counts = new Map();
  for (const atom of atoms) counts.set(atom.element, (counts.get(atom.element) || 0) + 1);
  const entries = [...counts.entries()].sort(([left], [right]) => left.localeCompare(right));
  if (!entries.length) return [];
  const minCount = Math.min(...entries.map(([, count]) => count));
  const maxCount = Math.max(...entries.map(([, count]) => count));
  const centers = entries.filter(([, count]) => count === minCount).map(([element]) => element);
  const ligands = entries.filter(([, count]) => count === maxCount).map(([element]) => element);
  return centers.map(center => {
    const ligand = ligands.find(element => element !== center) || ligands[0] || center;
    const distance = minimumDistance(atoms, vectors, center, ligand);
    return { center, ligand, min: 0, max: coordinationMaximum(distance) };
  });
}

function buildBondTopology(atoms, vectors) {
  const connections = [];
  for (const rule of defaultBondRules(atoms, vectors)) {
    const leftAtoms = atoms.filter(atom => atom.element === rule.center);
    const rightAtoms = atoms.filter(atom => atom.element === rule.ligand);
    const seen = new Set();
    for (const left of leftAtoms) {
      const leftPosition = combine(vectors, left.frac);
      for (const right of rightAtoms) {
        for (const shift of outsideShifts()) {
          if (left.sourceKey === right.sourceKey && shift.every(value => value === 0)) continue;
          const shifted = right.frac.map((value, axis) => value + shift[axis]);
          const distance = leftPosition.distanceTo(combine(vectors, shifted));
          if (distance < rule.min - 1e-8 || distance > rule.max + 1e-8) continue;
          const forward = `${left.sourceKey}>${right.sourceKey}@${shift.join(",")}`;
          const reverse = `${right.sourceKey}>${left.sourceKey}@${shift.map(value => -value).join(",")}`;
          const key = forward < reverse ? forward : reverse;
          if (seen.has(key)) continue;
          seen.add(key);
          connections.push({ leftKey: left.sourceKey, rightKey: right.sourceKey, shift: [...shift] });
        }
      }
    }
  }
  return connections;
}

function buildBondGeometry(atoms, displayedAtoms, vectors, topology, { showBonds, showPolyhedra, visualScale }) {
  const group = new THREE.Group();
  if (!showBonds && !showPolyhedra) return group;
  const referenceLength = Math.max(number(visualScale, 1), 1);
  const radius = referenceLength * 0.012;
  const atomFallback = referenceLength * 0.045;
  const atomByKey = new Map(atoms.map(atom => [atom.sourceKey, atom]));
  const displayedByKey = new Map();
  for (const atom of displayedAtoms) {
    if (!displayedByKey.has(atom.sourceKey)) displayedByKey.set(atom.sourceKey, []);
    displayedByKey.get(atom.sourceKey).push(atom);
  }
  const outsideSeen = new Set();
  const polyhedra = new Map();
  for (const connection of topology) {
    const sourceLeft = atomByKey.get(connection.leftKey);
    const right = atomByKey.get(connection.rightKey);
    if (!sourceLeft || !right) continue;
    for (const left of displayedByKey.get(connection.leftKey) || []) {
      const displayShift = left.displayShift || [0, 0, 0];
      const leftPosition = combine(vectors, left.frac);
      const shifted = right.frac.map((value, axis) => value + connection.shift[axis] + displayShift[axis]);
      const rightPosition = combine(vectors, shifted);
      const entryKey = left.visualKey || `${left.sourceKey}@${displayShift.join(",")}`;
      const entry = polyhedra.get(entryKey) || { atom: sourceLeft, points: [] };
      entry.points.push(rightPosition);
      polyhedra.set(entryKey, entry);
      if (showBonds) group.add(makeBondMesh(leftPosition, rightPosition, radius, sourceLeft, right));
      if (!isInsideUnit(shifted)) {
        const rightKey = rightPosition.toArray().map(value => fixed(value, 5)).join(",");
        const outsideKey = `${right.element}|${rightKey}`;
        if (!outsideSeen.has(outsideKey)) {
          outsideSeen.add(outsideKey);
          const outsideAtom = makeModeAtomMesh(right, atomFallback);
          outsideAtom.position.copy(rightPosition);
          group.add(outsideAtom);
        }
      }
    }
  }
  if (showPolyhedra) {
    for (const entry of polyhedra.values()) {
      const polyhedron = makePolyhedronMesh(entry.atom, entry.points);
      if (polyhedron) group.add(polyhedron);
    }
  }
  return group;
}

function elementFromLabel(value) {
  const match = String(value || "X").match(/[A-Z][a-z]?/);
  return match ? match[0] : "X";
}

function modeRows(definition) {
  return Array.isArray(definition?.rows) ? definition.rows : [];
}

function rowMoment(row) {
  return vector3(row?.moment ?? row?.mxyz ?? row?.magnetic_moment ?? row?.dxyz ?? row?.vector);
}

function strainMatrix(definition) {
  const components = definition?.components;
  if (Array.isArray(components) && components.length === 6) {
    const [e1, e2, e3, e4, e5, e6] = components.map(value => number(value));
    return [[e1, e6 / 2, e5 / 2], [e6 / 2, e2, e4 / 2], [e5 / 2, e4 / 2, e3]];
  }
  return matrix3(definition?.tensor ?? definition?.matrix ?? definition?.strain);
}

function disposeObject(root) {
  const disposedTextures = new Set();
  const disposeTexture = texture => {
    if (!texture || disposedTextures.has(texture)) return;
    disposedTextures.add(texture);
    texture.dispose?.();
  };
  root.traverse(object => {
    if (!object.geometry?.userData?.[SHARED_RESOURCE]) object.geometry?.dispose?.();
    const disposeMaterial = material => {
      if (material?.userData?.[SHARED_RESOURCE]) return;
      disposeTexture(material?.map);
      material?.dispose?.();
    };
    if (Array.isArray(object.material)) object.material.forEach(disposeMaterial);
    else disposeMaterial(object.material);
    disposeTexture(object.userData?.normalTexture);
    disposeTexture(object.userData?.hoverTexture);
  });
}

function arrow(origin, direction, color, scale) {
  const length = direction.length();
  if (length < EPS) return null;
  const headLength = Math.min(Math.max(length * 0.24, scale * 0.08), length * 0.55);
  const headWidth = Math.max(headLength * 0.45, scale * 0.035);
  return new THREE.ArrowHelper(direction.clone().normalize(), origin, length, color, headLength, headWidth);
}

function magneticArrow(origin, direction, scale) {
  const magnitude = direction.length();
  if (magnitude < EPS) return null;
  const amplitudeRatio = magnitude / REFERENCE_MOMENT;
  const totalLength = scale * REFERENCE_ARROW_LENGTH * amplitudeRatio;
  if (totalLength < scale * 0.008) return null;
  const widthScale = Math.min(1, Math.sqrt(amplitudeRatio));
  const headLength = Math.min(scale * REFERENCE_HEAD_LENGTH * widthScale, totalLength * 0.35);
  const shaftLength = Math.max(totalLength - headLength, scale * 0.002);
  const shaftRadius = scale * REFERENCE_SHAFT_RADIUS * widthScale;
  const headRadius = scale * REFERENCE_HEAD_RADIUS * widthScale;
  const group = new THREE.Group();
  const shaft = new THREE.Mesh(unitCylinderGeometry, magneticArrowMaterial);
  shaft.scale.set(shaftRadius, shaftLength, shaftRadius);
  shaft.position.y = -headLength / 2;
  const head = new THREE.Mesh(unitConeGeometry, magneticArrowMaterial);
  head.scale.set(headRadius, headLength, headRadius);
  head.position.y = shaftLength / 2;
  group.add(shaft, head);
  group.position.copy(origin);
  group.quaternion.setFromUnitVectors(Y_AXIS, direction.clone().normalize());
  return group;
}

function matrixFromVectors(vectors) {
  return new THREE.Matrix3().set(
    vectors[0].x, vectors[1].x, vectors[2].x,
    vectors[0].y, vectors[1].y, vectors[2].y,
    vectors[0].z, vectors[1].z, vectors[2].z,
  );
}

function axisLabel(text, color, position) {
  const canvas = document.createElement("canvas");
  canvas.width = canvas.height = 64;
  const context = canvas.getContext("2d");
  context.fillStyle = `#${color.toString(16).padStart(6, "0")}`;
  context.font = "44px sans-serif";
  context.textAlign = "center";
  context.textBaseline = "middle";
  context.fillText(text, 32, 34);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, toneMapped: false, depthTest: false }));
  sprite.position.copy(position);
  sprite.scale.setScalar(0.4);
  return sprite;
}

function cubeFaceTexture(hovered = false) {
  const canvas = document.createElement("canvas");
  canvas.width = canvas.height = 128;
  const context = canvas.getContext("2d");
  context.fillStyle = hovered ? "#c9d4e8" : "#e4e4e1";
  context.fillRect(0, 0, 128, 128);
  context.strokeStyle = "#c2c2bc";
  context.lineWidth = 2;
  context.strokeRect(1, 1, 126, 126);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

function cubeFace(corners, up) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(new Float32Array([
    ...corners[0].toArray(), ...corners[1].toArray(), ...corners[2].toArray(),
    ...corners[0].toArray(), ...corners[2].toArray(), ...corners[3].toArray(),
  ]), 3));
  geometry.setAttribute("uv", new THREE.BufferAttribute(new Float32Array([0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1]), 2));
  const normalTexture = cubeFaceTexture(false);
  const mesh = new THREE.Mesh(geometry, new THREE.MeshBasicMaterial({ map: normalTexture, side: THREE.DoubleSide, toneMapped: false }));
  const edge1 = corners[1].clone().sub(corners[0]);
  const edge2 = corners[2].clone().sub(corners[0]);
  const normal = edge1.cross(edge2).normalize();
  const centroid = corners.reduce((sum, point) => sum.add(point), new THREE.Vector3()).multiplyScalar(0.25);
  if (normal.dot(centroid) < 0) normal.negate();
  mesh.userData = {
    viewDir: normal,
    up: up.clone(),
    normalTexture,
    hoverTexture: cubeFaceTexture(true),
  };
  return mesh;
}

function solidArrow(direction, color, length, shaftRadius, headLength, headRadius) {
  const material = new THREE.MeshBasicMaterial({ color, depthTest: false, depthWrite: false, toneMapped: false });
  const group = new THREE.Group();
  const shaftLength = Math.max(0.01, length - headLength);
  const shaft = new THREE.Mesh(new THREE.CylinderGeometry(shaftRadius, shaftRadius, shaftLength, 20), material);
  shaft.position.y = shaftLength / 2;
  const head = new THREE.Mesh(new THREE.ConeGeometry(headRadius, headLength, 24), material);
  head.position.y = shaftLength + headLength / 2;
  group.add(shaft, head);
  group.quaternion.setFromUnitVectors(Y_AXIS, direction.clone().normalize());
  group.renderOrder = 10;
  return group;
}

export class ModeStructureViewer {
  constructor({ canvas, empty, controls, status, playButton, resetButton, scaleInput, scaleOutput, arrowsInput, momentsInput, bondsInput, polyhedraInput }) {
    this.canvas = canvas;
    this.empty = empty;
    this.controlsNode = controls;
    this.statusNode = status;
    this.playButton = playButton;
    this.resetButton = resetButton;
    this.scaleInput = scaleInput;
    this.scaleOutput = scaleOutput;
    this.arrowsInput = arrowsInput;
    this.momentsInput = momentsInput;
    this.bondsInput = bondsInput;
    this.polyhedraInput = polyhedraInput;
    this.payload = null;
    this.modes = [];
    this.amplitudes = new Map();
    this.structure = null;
    this.playing = false;
    this.masterAmplitude = 1;
    this.animationStart = 0;
    this.animationPhaseOffset = 0;
    this.lastAnimationUpdate = 0;
    this.needsRender = true;
    this.fitted = false;
    this.bondTopology = null;

    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setClearColor(0xffffff, 1);
    this.renderer.autoClear = false;
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.scene = new THREE.Scene();
    addLights(this.scene);
    this.camera = new THREE.OrthographicCamera(-5, 5, 5, -5, -1000, 1000);
    this.camera.position.set(9, 8, 6);
    this.camera.up.set(0, 0, 1);
    this.trackball = new TrackballControls(this.camera, canvas);
    this.trackball.staticMoving = true;
    this.trackball.rotateSpeed = 3.6;
    this.trackball.zoomSpeed = 2.4;
    this.trackball.panSpeed = 0.65;
    this.trackball.addEventListener("change", () => { this.needsRender = true; });
    this.cubeScene = new THREE.Scene();
    this.cubeCamera = new THREE.OrthographicCamera(-CUBE_VIEW, CUBE_VIEW, CUBE_VIEW, -CUBE_VIEW, -10, 10);
    this.cubeRay = new THREE.Raycaster();
    this.cubeDirection = new THREE.Vector3();
    this.cubeFaceMeshes = [];
    this.hoveredCubeFace = null;
    this.cubeDownFace = null;
    this.viewTween = null;
    this.viewCube = null;

    loadElementStyles("/tools/cif_standardizer/vesta_elements.csv").then(() => this.updateScene());
    this.bindControls();
    this.bindViewCube();
    this.resizeObserver = new ResizeObserver(() => {
      this.resize();
      const rect = this.canvas.getBoundingClientRect();
      if (this.structure && rect.width > 10 && rect.height > 10) this.fitView(true);
    });
    this.resizeObserver.observe(canvas.parentElement || canvas);
    this.animate = this.animate.bind(this);
    requestAnimationFrame(this.animate);
  }

  bindControls() {
    this.playButton?.addEventListener("click", () => {
      this.playing = !this.playing;
      this.playButton.textContent = this.playing ? "Pause" : "Animate";
      this.playButton.setAttribute("aria-pressed", String(this.playing));
    });
    this.resetButton?.addEventListener("click", () => this.resetModeAmplitudes());
    [this.scaleInput, this.arrowsInput, this.momentsInput, this.bondsInput, this.polyhedraInput].forEach(input => {
      input?.addEventListener("input", () => {
        if (input === this.scaleInput && this.scaleOutput) {
          this.scaleOutput.textContent = number(this.scaleInput.value, 1).toFixed(2);
        }
        this.updateScene();
      });
      input?.addEventListener("change", () => this.updateScene());
    });
  }

  resetModeAmplitudes() {
    this.playing = false;
    for (const mode of this.modes) this.amplitudes.set(mode.id, 0);
    this.masterAmplitude = 1;
    if (this.playButton) {
      this.playButton.textContent = "Animate";
      this.playButton.setAttribute("aria-pressed", "false");
    }
    this.updateControlValues();
    this.updateScene();
  }

  bindViewCube() {
    this.canvas.addEventListener("pointermove", event => {
      if (event.buttons === 0) this.setCubeHover(this.pickCubeFace(event.clientX, event.clientY));
    });
    this.canvas.addEventListener("pointerleave", () => {
      this.cubeDownFace = null;
      this.setCubeHover(null);
    });
    this.canvas.addEventListener("pointerdown", event => {
      const face = this.pickCubeFace(event.clientX, event.clientY);
      if (!face) return;
      this.cubeDownFace = face;
      event.stopImmediatePropagation();
      event.preventDefault();
    }, true);
    this.canvas.addEventListener("pointerup", event => {
      if (!this.cubeDownFace) return;
      const face = this.pickCubeFace(event.clientX, event.clientY);
      if (face && face === this.cubeDownFace) this.goToCubeView(face.userData.viewDir, face.userData.up);
      this.cubeDownFace = null;
      event.stopImmediatePropagation();
    }, true);
  }

  setData(payload) {
    this.payload = payload && typeof payload === "object" ? payload : null;
    this.modes = [];
    this.amplitudes.clear();
    for (const [kind, groupLabel, key] of MODE_KINDS) {
      const definitions = Array.isArray(this.payload?.[key]) ? this.payload[key] : [];
      definitions.forEach((definition, index) => {
        const mode = {
          id: `${kind}-${index}`,
          kind,
          groupLabel,
          definition,
          label: String(definition?.label || `${groupLabel} ${index + 1}`),
          siteOrder: Number.isInteger(definition?._site_order) ? definition._site_order : 0,
        };
        this.modes.push(mode);
        this.amplitudes.set(mode.id, 0);
      });
    }
    this.playing = false;
    this.masterAmplitude = 1;
    this.fitted = false;
    this.bondTopology = null;
    if (this.scaleInput) this.scaleInput.value = "1";
    if (this.scaleOutput) this.scaleOutput.textContent = "1.00";
    if (this.playButton) {
      this.playButton.textContent = "Animate";
      this.playButton.setAttribute("aria-pressed", "false");
    }
    this.renderModeControls();
    this.updateScene();
  }

  renderModeControls() {
    if (!this.controlsNode) return;
    if (!this.payload || this.payload.status === "unsupported" || this.payload.status === "missing") {
      this.controlsNode.innerHTML = "";
      if (this.statusNode) this.statusNode.textContent = this.payload?.reason || "Mode data is not available.";
      return;
    }
    const groupedModes = new Map();
    for (const mode of this.modes) {
      const siteLabel = mode.kind === "strain" ? "Strain" : this.modeSiteLabel(mode);
      if (!groupedModes.has(siteLabel)) groupedModes.set(siteLabel, []);
      groupedModes.get(siteLabel).push(mode);
    }
    const groups = [...groupedModes.entries()].map(([siteLabel, modes]) => {
      return `<details class="mode-group" open><summary>${this.escape(siteLabel)} Modes</summary>${modes.map(mode => {
        const value = this.amplitudes.get(mode.id) || 0;
        const limit = this.modeAmplitudeLimit(mode);
        const step = Math.max(limit / 500, 0.001);
        return `<label class="mode-slider" title="${this.escape(mode.label)}">
          <input type="range" min="${-limit}" max="${limit}" step="${step}" value="${value}" data-mode-id="${mode.id}">
          <output>${value.toFixed(3)}</output>
          <span>${this.escape(this.modeDisplayLabel(mode))}</span>
        </label>`;
      }).join("")}</details>`;
    }).join("");
    const master = `<section class="mode-group mode-master-group">
      <h4><span>Master Amp</span><span class="mode-master-actions"><button class="primary-action" type="button" data-mode-animate aria-pressed="false">Animate</button><button class="primary-action" type="button" data-mode-reset>Reset</button></span></h4>
      <div class="mode-master">
        <span>Parent</span>
        <input type="range" min="0" max="1" step="0.01" value="1" data-master-amplitude>
        <span>Child</span>
        <output>1.00</output>
      </div>
    </section>`;
    this.controlsNode.innerHTML = `${master}${groups || '<p class="mode-viewer-note">No visualizable mode definitions.</p>'}`;
    const masterInput = this.controlsNode.querySelector("input[data-master-amplitude]");
    masterInput?.addEventListener("input", () => {
      this.playing = false;
      if (this.playButton) {
        this.playButton.textContent = "Animate";
        this.playButton.setAttribute("aria-pressed", "false");
      }
      this.masterAmplitude = Math.min(1, Math.max(0, number(masterInput.value, 1)));
      this.updateControlValues();
      this.updateScene();
    });
    this.playButton = this.controlsNode.querySelector("button[data-mode-animate]");
    this.resetButton = this.controlsNode.querySelector("button[data-mode-reset]");
    this.playButton?.addEventListener("click", () => {
      this.playing = !this.playing;
      if (this.playing) {
        this.animationStart = performance.now();
        this.animationPhaseOffset = Math.asin(Math.sqrt(Math.min(1, Math.max(0, this.masterAmplitude))));
      }
      this.playButton.textContent = this.playing ? "Pause" : "Animate";
      this.playButton.setAttribute("aria-pressed", String(this.playing));
    });
    this.resetButton?.addEventListener("click", () => this.resetModeAmplitudes());
    this.controlsNode.querySelectorAll("input[data-mode-id]").forEach(input => {
      input.addEventListener("input", () => {
        const effective = number(input.value);
        const unscaled = this.masterAmplitude > EPS ? effective / this.masterAmplitude : effective;
        const value = Math.min(number(input.max, unscaled), Math.max(number(input.min, unscaled), unscaled));
        this.amplitudes.set(input.dataset.modeId, value);
        this.updateControlValues();
        this.updateScene();
      });
    });
    const counts = MODE_KINDS.map(([kind, label]) => {
      const count = this.modes.filter(mode => mode.kind === kind).length;
      return count ? `${count} ${label.toLowerCase()}` : "";
    }).filter(Boolean);
    if (this.statusNode) this.statusNode.textContent = counts.length ? counts.join(" / ") : "Undistorted structure";
  }

  updateControlValues() {
    if (!this.controlsNode) return;
    const masterInput = this.controlsNode.querySelector("input[data-master-amplitude]");
    if (masterInput) masterInput.value = String(this.masterAmplitude);
    const masterOutput = masterInput?.closest(".mode-master")?.querySelector("output");
    if (masterOutput) masterOutput.textContent = this.masterAmplitude.toFixed(2);
    this.controlsNode.querySelectorAll("input[data-mode-id]").forEach(input => {
      const amplitude = this.amplitudes.get(input.dataset.modeId) || 0;
      input.value = String(amplitude * this.masterAmplitude);
      const output = input.closest(".mode-slider")?.querySelector("output");
      if (output) output.textContent = (amplitude * this.masterAmplitude).toFixed(3);
    });
  }

  modeAmplitudeLimit(mode) {
    if (mode.kind === "strain") return 0.1;
    const rows = modeRows(mode.definition);
    const norm = Math.abs(number(mode.definition?.normfactor, 1)) || 1;
    const parentCell = this.payload?.viewer_parent?.lattice;
    const parentVectors = cellVectors(parentCell || this.payload?.lattice);
    const basis = this.payload?.subgroup_details?.presentation_basis
      ?? this.payload?.subgroup_details?.display_basis
      ?? this.payload?.subgroup_details?.basis;
    const vectors = vectorsFromBasis(parentVectors, basis) || cellVectors(this.payload?.lattice);
    let largest = 0;
    for (const row of rows) {
      const vector = mode.kind === "magnetic" ? rowMoment(row) : vector3(row?.dxyz);
      largest = Math.max(largest, combine(vectors, vector).length() * norm);
    }
    if (mode.kind === "magnetic") {
      const supplied = number(mode.definition?.max_amplitude ?? mode.definition?.maxAmp, 0);
      return supplied > 0 ? supplied : (largest > EPS ? 1 / largest : 1);
    }
    // ISOVIZ supplies a per-mode max amplitude. Until that value is emitted by
    // the local mode payload, keep the largest single-atom excursion near 1.2 A.
    return largest > EPS ? Math.min(4, Math.max(0.25, 1.2 / largest)) : 1;
  }

  modeDisplayLabel(mode) {
    const label = String(mode.label || "");
    const parentEnd = label.indexOf("]");
    const localLabel = parentEnd >= 0 ? label.slice(parentEnd + 1) : label;
    const siteStart = localLabel.indexOf("[");
    if (siteStart < 0) return localLabel;
    const irrep = localLabel.slice(0, siteStart).replace(/\([^()]*\)$/, "");
    return `${irrep}${localLabel.slice(siteStart)}`;
  }

  modeSiteLabel(mode) {
    const match = String(mode.label || "").match(/\[([^:\]]+):[^\]]+\]/);
    if (match) return match[1];
    const site = this.payload?.sites?.[mode.siteOrder];
    return String(site?.site || site?.label || `Site ${mode.siteOrder + 1}`);
  }

  escape(value) {
    return String(value ?? "").replace(/[&<>"']/g, char => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    }[char]));
  }

  buildAtoms() {
    const sourceAtoms = [];
    const viewerAtoms = Array.isArray(this.payload?.viewer_atoms) ? this.payload.viewer_atoms : [];
    if (viewerAtoms.length) {
      viewerAtoms.forEach((atom, fallbackIndex) => {
        const xyz = atom?.xyz ?? atom?.frac;
        if (!Array.isArray(xyz) || xyz.length < 3) return;
        const siteOrder = Number.isInteger(atom.site_order) ? atom.site_order : 0;
        const atomIndex = Number.isInteger(atom.atom_index) ? atom.atom_index : fallbackIndex;
        const label = String(atom.label || `${atom.element || `X${siteOrder}`}_${atomIndex + 1}`);
        sourceAtoms.push({
          sourceKey: `${siteOrder}:${atomIndex}`,
          siteOrder,
          atomIndex,
          label,
          element: elementFromLabel(atom.element || label),
          frac: vector3(xyz),
          sourceFrac: vector3(xyz),
        });
      });
    } else {
      const siteOrders = [...new Set(this.modes.map(mode => mode.siteOrder))].sort((a, b) => a - b);
      for (const siteOrder of siteOrders) {
        const mode = this.modes.find(item => item.siteOrder === siteOrder && modeRows(item.definition).length)
          || this.modes.find(item => item.siteOrder === siteOrder);
        const siteLabel = this.modeSiteLabel(mode);
        modeRows(mode?.definition).forEach((row, atomIndex) => {
          const xyz = row?.xyz ?? row?.frac;
          if (!Array.isArray(xyz) || xyz.length < 3) return;
          const label = String(row.atom || row.label || `${siteLabel}_${atomIndex + 1}`);
          sourceAtoms.push({
            sourceKey: `${siteOrder}:${atomIndex}`,
            siteOrder,
            atomIndex,
            label,
            element: elementFromLabel(siteLabel || label),
            frac: vector3(xyz),
            sourceFrac: vector3(xyz),
          });
        });
      }
    }
    if (!sourceAtoms.length) {
      const fallback = this.payload?.undistorted_atoms || this.payload?.distorted_atoms || [];
      fallback.forEach((atom, atomIndex) => {
        const xyz = atom.xyz ?? atom.frac;
        if (!Array.isArray(xyz) || xyz.length < 3) return;
        const label = String(atom.label ?? atom.atom ?? atom.type_symbol ?? `X_${atomIndex + 1}`);
        sourceAtoms.push({
          sourceKey: `fallback:${atomIndex}`,
          siteOrder: -1,
          atomIndex,
          label,
          element: elementFromLabel(label),
          frac: vector3(xyz),
          sourceFrac: vector3(xyz),
        });
      });
    }

    // The official atomcoordlist includes periodic copies on cell faces and
    // corners. Recreate those display copies without changing mode identity.
    const displayed = new Map();
    for (const atom of sourceAtoms) {
      const offsets = atom.frac.map(component => {
        if (Math.abs(component) < 1e-7) return [0, 1];
        if (Math.abs(component - 1) < 1e-7) return [0, -1];
        return [0];
      });
      for (const dx of offsets[0]) for (const dy of offsets[1]) for (const dz of offsets[2]) {
        const frac = [atom.frac[0] + dx, atom.frac[1] + dy, atom.frac[2] + dz];
        const visualKey = `${atom.sourceKey}:${frac.map(value => value.toFixed(8)).join(",")}`;
        if (!displayed.has(visualKey)) displayed.set(visualKey, {
          ...atom,
          visualKey,
          displayShift: [dx, dy, dz],
          frac,
        });
      }
    }
    return [...displayed.values()];
  }

  deformationMatrix(scale) {
    const deformation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (const mode of this.modes.filter(item => item.kind === "strain")) {
      const tensor = strainMatrix(mode.definition);
      if (!tensor) continue;
      addMatrix(deformation, tensor, (this.amplitudes.get(mode.id) || 0) * scale * this.masterAmplitude);
    }
    return deformation;
  }

  modeVectors(atoms, scale) {
    const sourceKeys = [...new Set(atoms.map(atom => atom.sourceKey))];
    const displacements = new Map(sourceKeys.map(key => [key, [0, 0, 0]]));
    const moments = new Map(sourceKeys.map(key => [key, [0, 0, 0]]));
    for (const mode of this.modes) {
      const normfactor = Math.abs(number(mode.definition?.normfactor, 1));
      const normalization = normfactor > EPS ? normfactor : 1;
      const amplitude = (this.amplitudes.get(mode.id) || 0) * scale * this.masterAmplitude * normalization;
      if (Math.abs(amplitude) < EPS) continue;
      modeRows(mode.definition).forEach((row, atomIndex) => {
        const key = `${mode.siteOrder}:${atomIndex}`;
        if (!displacements.has(key)) return;
        if (mode.kind === "displacive") {
          addVector(displacements.get(key), vector3(row.dxyz ?? row.displacement ?? row.vector), amplitude);
        } else if (mode.kind === "magnetic") {
          addVector(moments.get(key), rowMoment(row), amplitude);
        }
      });
    }
    return { displacements, moments };
  }

  buildViewCube(vectors) {
    if (this.viewCube) {
      this.setCubeHover(null);
      this.cubeScene.remove(this.viewCube);
      disposeObject(this.viewCube);
    }
    this.cubeFaceMeshes = [];
    const matrix = matrixFromVectors(vectors);
    const center = new THREE.Vector3(0.5, 0.5, 0.5).applyMatrix3(matrix);
    const corner = (i, j, k) => new THREE.Vector3(i, j, k).applyMatrix3(matrix).sub(center);
    let maximum = EPS;
    for (let i = 0; i < 2; i += 1) for (let j = 0; j < 2; j += 1) for (let k = 0; k < 2; k += 1) {
      maximum = Math.max(maximum, corner(i, j, k).length());
    }
    const scale = 0.66 / maximum;
    const c = (i, j, k) => corner(i, j, k).multiplyScalar(scale);
    const group = new THREE.Group();
    [
      { corners: [c(1, 0, 0), c(1, 1, 0), c(1, 1, 1), c(1, 0, 1)], up: UP_Z },
      { corners: [c(0, 0, 0), c(0, 0, 1), c(0, 1, 1), c(0, 1, 0)], up: UP_Z },
      { corners: [c(0, 1, 0), c(0, 1, 1), c(1, 1, 1), c(1, 1, 0)], up: UP_Z },
      { corners: [c(0, 0, 0), c(1, 0, 0), c(1, 0, 1), c(0, 0, 1)], up: UP_Z },
      { corners: [c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)], up: UP_Y },
      { corners: [c(0, 0, 0), c(0, 1, 0), c(1, 1, 0), c(1, 0, 0)], up: UP_Y },
    ].forEach(({ corners, up }) => {
      const face = cubeFace(corners, up);
      this.cubeFaceMeshes.push(face);
      group.add(face);
    });
    const origin = c(0, 0, 0);
    [
      [[1, 0, 0], 0xc0392b, "a"],
      [[0, 1, 0], 0x27ae60, "b"],
      [[0, 0, 1], 0x2980b9, "c"],
    ].forEach(([direction, color, label]) => {
      const edge = new THREE.Vector3(...direction).applyMatrix3(matrix).multiplyScalar(scale);
      const length = edge.length() * 1.28;
      const unit = edge.clone().normalize();
      const axis = solidArrow(unit, color, length, 0.016, 0.12, 0.05);
      axis.position.copy(origin);
      group.add(axis);
      group.add(axisLabel(label, color, origin.clone().add(unit.multiplyScalar(length + 0.22))));
    });
    this.viewCube = group;
    this.cubeScene.add(group);
  }

  cubeNdc(clientX, clientY) {
    const rect = this.canvas.getBoundingClientRect();
    const x = clientX - rect.left;
    const y = clientY - rect.top;
    const left = rect.width - CUBE_PX - CUBE_MARGIN;
    if (x < left || x > left + CUBE_PX || y < CUBE_MARGIN || y > CUBE_MARGIN + CUBE_PX) return null;
    return new THREE.Vector2(
      ((x - left) / CUBE_PX) * 2 - 1,
      -(((y - CUBE_MARGIN) / CUBE_PX) * 2 - 1),
    );
  }

  pickCubeFace(clientX, clientY) {
    const ndc = this.cubeNdc(clientX, clientY);
    if (!ndc || !this.cubeFaceMeshes.length) return null;
    this.cubeRay.setFromCamera(ndc, this.cubeCamera);
    return this.cubeRay.intersectObjects(this.cubeFaceMeshes, false)[0]?.object || null;
  }

  setCubeHover(face) {
    if (this.hoveredCubeFace === face) return;
    if (this.hoveredCubeFace) {
      this.hoveredCubeFace.material.map = this.hoveredCubeFace.userData.normalTexture;
      this.hoveredCubeFace.material.needsUpdate = true;
    }
    this.hoveredCubeFace = face;
    if (face) {
      face.material.map = face.userData.hoverTexture;
      face.material.needsUpdate = true;
    }
    this.canvas.style.cursor = face ? "pointer" : "";
    this.needsRender = true;
  }

  goToCubeView(viewDirection, up) {
    const distance = this.camera.position.distanceTo(this.trackball.target);
    this.viewTween = {
      fromPosition: this.camera.position.clone(),
      toPosition: this.trackball.target.clone().add(viewDirection.clone().normalize().multiplyScalar(distance)),
      fromUp: this.camera.up.clone(),
      toUp: up.clone().normalize(),
      start: performance.now(),
      duration: 320,
    };
    this.needsRender = true;
  }

  applyViewTween() {
    if (!this.viewTween) return;
    const elapsed = Math.min(1, (performance.now() - this.viewTween.start) / this.viewTween.duration);
    const eased = elapsed * elapsed * (3 - 2 * elapsed);
    this.camera.position.lerpVectors(this.viewTween.fromPosition, this.viewTween.toPosition, eased);
    this.camera.up.lerpVectors(this.viewTween.fromUp, this.viewTween.toUp, eased).normalize();
    this.camera.lookAt(this.trackball.target);
    this.needsRender = true;
    if (elapsed >= 1) this.viewTween = null;
  }

  updateScene() {
    if (this.structure) {
      this.scene.remove(this.structure);
      disposeObject(this.structure);
      this.structure = null;
    }
    const lattice = this.payload?.lattice;
    if (!this.payload || !lattice) {
      this.empty.hidden = false;
      this.empty.textContent = this.payload?.reason || "Calculate mode details to open the viewer";
      return;
    }

    const scale = Math.max(number(this.scaleInput?.value, 1), 0);
    const parentCell = this.payload?.viewer_parent?.lattice;
    const parentVectors = cellVectors(parentCell || lattice);
    const childBasis = this.payload?.subgroup_details?.presentation_basis
      ?? this.payload?.subgroup_details?.display_basis
      ?? this.payload?.subgroup_details?.basis;
    const baseVectors = vectorsFromBasis(parentVectors, childBasis) || cellVectors(lattice);
    const deformation = this.deformationMatrix(scale);
    const deformVectors = source => source.map(vector => {
      const transformed = multiplyMatrixVector(deformation, [vector.x, vector.y, vector.z]);
      return new THREE.Vector3(...transformed);
    });
    const distortedParentVectors = deformVectors(parentVectors);
    const vectors = deformVectors(baseVectors);
    const atoms = this.buildAtoms();
    const { displacements, moments } = this.modeVectors(atoms, scale);
    // Viewer primitives are sized from the parent cell, not the child
    // supercell. Otherwise bond/arrow thickness grows with supercell index.
    const visualScale = Math.max(...parentVectors.map(vector => vector.length()), 1);
    const group = new THREE.Group();
    this.buildViewCube(vectors);
    if (parentCell) {
      // Official ISOVIZ receives !parentorigin explicitly in child-cell units.
      // The local payload does not expose it yet; do not substitute OPD origin,
      // which is a different coordinate quantity.
      const viewerOrigin = this.payload?.viewer_parent_origin;
      const parentOffset = Array.isArray(viewerOrigin)
        ? combine(vectors, vector3(viewerOrigin))
        : new THREE.Vector3();
      // ISODISTORT first strains the parent lattice, then applies the subgroup
      // basis transform to obtain the distorted child cell. Both outlines
      // therefore share the same strain tensor.
      group.add(makeCellEdges(distortedParentVectors, 0xc86b6b, parentOffset, 0.7));
    }
    group.add(makeCellEdges(vectors, 0x6879c8));
    const sourceAtoms = new Map();
    for (const atom of atoms) {
      if (sourceAtoms.has(atom.sourceKey)) continue;
      const displacement = displacements.get(atom.sourceKey) || [0, 0, 0];
      sourceAtoms.set(atom.sourceKey, {
        ...atom,
        frac: atom.sourceFrac.map((value, axis) => value + displacement[axis]),
      });
    }
    if (!this.bondTopology) {
      const referenceAtoms = [...sourceAtoms.values()].map(atom => ({ ...atom, frac: [...atom.sourceFrac] }));
      this.bondTopology = buildBondTopology(referenceAtoms, baseVectors);
    }
    const displayedAtoms = atoms.map(atom => {
      const displacement = displacements.get(atom.sourceKey) || [0, 0, 0];
      return {
        ...atom,
        frac: atom.frac.map((value, axis) => value + displacement[axis]),
      };
    });
    group.add(buildBondGeometry([...sourceAtoms.values()], displayedAtoms, vectors, this.bondTopology, {
      showBonds: Boolean(this.bondsInput?.checked),
      showPolyhedra: Boolean(this.polyhedraInput?.checked),
      visualScale,
    }));
    for (const atom of atoms) {
      const basePosition = combine(vectors, atom.frac);
      const displacement = combine(vectors, displacements.get(atom.sourceKey) || [0, 0, 0]);
      const position = basePosition.clone().add(displacement);
      const fallbackRadius = visualScale * 0.045;
      const mesh = makeModeAtomMesh(atom, fallbackRadius);
      mesh.position.copy(position);
      mesh.userData = { label: atom.label, fractional: atom.frac };
      group.add(mesh);
      if (this.arrowsInput?.checked) {
        const displacementArrow = arrow(basePosition, displacement, 0x1677b8, visualScale);
        if (displacementArrow) group.add(displacementArrow);
      }
      if (!this.momentsInput || this.momentsInput.checked) {
        const moment = combine(vectors, moments.get(atom.sourceKey) || [0, 0, 0]);
        const momentArrow = magneticArrow(position, moment, visualScale);
        if (momentArrow) group.add(momentArrow);
      }
    }
    this.structure = group;
    this.scene.add(group);
    this.needsRender = true;
    this.empty.hidden = atoms.length > 0;
    this.empty.textContent = atoms.length ? "" : "No atom positions in mode data";
    if (atoms.length && this.statusNode) {
      this.statusNode.textContent = `${atoms.length} displayed atoms / ${this.modes.length} modes`;
    }
    const canvasRect = this.canvas.getBoundingClientRect();
    if (!this.fitted && canvasRect.width > 10 && canvasRect.height > 10) {
      this.fitView();
      this.fitted = true;
    }
  }

  refreshLayout() {
    this.resize();
    if (!this.structure) return;
    this.fitView();
    this.fitted = true;
  }

  fitView(preserveOrientation = false) {
    if (!this.structure) return;
    const box = new THREE.Box3().setFromObject(this.structure);
    const sphere = box.getBoundingSphere(new THREE.Sphere());
    const radius = Math.max(sphere.radius, 1);
    const center = sphere.center;
    const rect = this.canvas.getBoundingClientRect();
    const aspect = Math.max(rect.width, 1) / Math.max(rect.height, 1);
    const view = radius * 1.2;
    this.camera.left = -view * aspect;
    this.camera.right = view * aspect;
    this.camera.top = view;
    this.camera.bottom = -view;
    const cameraOffset = this.camera.position.clone().sub(this.trackball.target);
    const direction = preserveOrientation && cameraOffset.lengthSq() > EPS
      ? cameraOffset.normalize()
      : new THREE.Vector3(2.4, 2.1, 1.45).normalize();
    this.camera.position.copy(center).add(direction.multiplyScalar(radius * 3.5));
    this.camera.up.set(0, 0, 1);
    this.camera.updateProjectionMatrix();
    this.trackball.target.copy(center);
    this.trackball.update();
  }

  resize() {
    const rect = this.canvas.getBoundingClientRect();
    this.renderer.setSize(Math.max(rect.width, 1), Math.max(rect.height, 1), false);
    if (this.structure) {
      const aspect = Math.max(rect.width, 1) / Math.max(rect.height, 1);
      const height = Math.max(this.camera.top - this.camera.bottom, 1);
      this.camera.left = -height * aspect / 2;
      this.camera.right = height * aspect / 2;
      this.camera.updateProjectionMatrix();
    }
    this.trackball.handleResize();
    this.needsRender = true;
  }

  animate(time) {
    requestAnimationFrame(this.animate);
    this.applyViewTween();
    if (this.playing && time - this.lastAnimationUpdate >= 50) {
      this.lastAnimationUpdate = time;
      const phase = this.animationPhaseOffset + (time - this.animationStart) * (2 * Math.PI / 5000);
      const nextMaster = Math.sin(phase) ** 2;
      if (Math.abs(nextMaster - this.masterAmplitude) > 0.01) {
        this.masterAmplitude = nextMaster;
        this.updateControlValues();
        this.updateScene();
      }
    }
    this.trackball.update();
    if (!this.needsRender) return;
    this.needsRender = false;
    const rect = this.canvas.getBoundingClientRect();
    const width = Math.max(rect.width, 1);
    const height = Math.max(rect.height, 1);
    this.renderer.setViewport(0, 0, width, height);
    this.renderer.setScissorTest(false);
    this.renderer.clear();
    this.renderer.render(this.scene, this.camera);
    if (this.viewCube) {
      this.cubeDirection.copy(this.camera.position).sub(this.trackball.target).normalize().multiplyScalar(4);
      this.cubeCamera.position.copy(this.cubeDirection);
      this.cubeCamera.up.copy(this.camera.up);
      this.cubeCamera.lookAt(0, 0, 0);
      const x = width - CUBE_PX - CUBE_MARGIN;
      const y = height - CUBE_PX - CUBE_MARGIN;
      this.renderer.setViewport(x, y, CUBE_PX, CUBE_PX);
      this.renderer.setScissor(x, y, CUBE_PX, CUBE_PX);
      this.renderer.setScissorTest(true);
      this.renderer.clearDepth();
      this.renderer.render(this.cubeScene, this.cubeCamera);
      this.renderer.setScissorTest(false);
    }
  }
}
