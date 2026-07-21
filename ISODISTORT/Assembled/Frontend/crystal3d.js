// 結晶構造の 3D 描画で共通に使う見た目のプリミティブ。
// cif_standardizer / spacegroup_viewer などが import して同じ見た目を共有する。
// （セル枠の線・原子マテリアル・元素配色・ライト）
import * as THREE from "three";

export const ATOM_RADIUS_SCALE = 0.34;

// 元素 → { radius, color } のスタイル表（VESTA 由来）。
export const elementStyles = new Map();

export function canonicalElement(value) {
  const text = String(value || "").replace(/[^A-Za-z]/g, "");
  if (!text) return "?";
  if (text.toUpperCase() === "D") return "D";
  return text[0].toUpperCase() + text.slice(1).toLowerCase();
}

export function atomElement(atom) {
  return canonicalElement(atom.element || atom.type_symbol || "?");
}

async function fetchTextOrEmpty(path) {
  try {
    const res = await fetch(path);
    return res.ok ? await res.text() : "";
  } catch {
    return "";
  }
}

// 各ツールのディレクトリに置いた vesta_elements.csv を読み込む。
export async function loadElementStyles(path = "vesta_elements.csv") {
  const text = await fetchTextOrEmpty(path);
  for (const line of text.split("\n").slice(1)) {
    const [symbolRaw, radiusRaw, , , colorRaw] = line.split(",");
    const radius = Number(radiusRaw);
    if (!symbolRaw) continue;
    const symbol = canonicalElement(symbolRaw);
    const current = elementStyles.get(symbol) || {};
    if (Number.isFinite(radius)) current.radius = radius;
    if (colorRaw?.startsWith("#")) current.color = parseInt(colorRaw.slice(1), 16);
    elementStyles.set(symbol, current);
  }
}

export function atomColor(atom) {
  return elementStyles.get(atomElement(atom))?.color ?? 0x777777;
}

export function atomRadius(atom, fallback) {
  const radius = elementStyles.get(atomElement(atom))?.radius;
  return Number.isFinite(radius) ? radius * ATOM_RADIUS_SCALE : fallback;
}

export function colorCss(value) {
  return `#${value.toString(16).padStart(6, "0")}`;
}

export function atomMaterial(atom) {
  const color = new THREE.Color(atomColor(atom));
  return new THREE.MeshStandardMaterial({
    color,
    roughness: 0.42,
    metalness: 0.05,
    emissive: color.clone().multiplyScalar(0.08),
  });
}

// 原子球。atom は {element|type_symbol} を持つ。位置は呼び出し側で設定する。
export function makeAtomMesh(atom, fallbackRadius) {
  return new THREE.Mesh(
    new THREE.SphereGeometry(atomRadius(atom, fallbackRadius), 28, 20),
    atomMaterial(atom),
  );
}

// MeshStandardMaterial 用のライト（cif_standardizer と同じ構成）。
export function addLights(scene) {
  scene.add(new THREE.HemisphereLight(0xffffff, 0xd7dde8, 1.8));
  const keyLight = new THREE.DirectionalLight(0xffffff, 2.4);
  keyLight.position.set(2.5, -3.0, 4.5);
  scene.add(keyLight);
  const fillLight = new THREE.DirectionalLight(0xffffff, 0.8);
  fillLight.position.set(-3.0, 2.0, 2.0);
  scene.add(fillLight);
}

export function makeLine(points, color, opacity = 1, dashed = false) {
  const geo = new THREE.BufferGeometry().setFromPoints(points);
  const material = dashed
    ? new THREE.LineDashedMaterial({
        color, opacity, transparent: opacity < 1, depthWrite: opacity >= 1,
        dashSize: 0.18, gapSize: 0.12, linewidth: 2,
      })
    : new THREE.LineBasicMaterial({
        color, opacity, transparent: opacity < 1, depthWrite: opacity >= 1, linewidth: 2,
      });
  const line = new THREE.LineSegments(geo, material);
  if (dashed) line.computeLineDistances();
  return line;
}

// 細い線を 3 本わずかにずらして重ね、太く見せる（線幅非対応環境向け）。
export function makeBoldLine(points, color, opacity = 1, dashed = false) {
  const group = new THREE.Group();
  const offsets = [
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0.012, 0.012, 0),
    new THREE.Vector3(-0.012, 0.012, 0),
  ];
  for (const offset of offsets) {
    group.add(makeLine(points.map((point) => point.clone().add(offset)), color, opacity, dashed));
  }
  return group;
}

// 平行六面体（基底ベクトル a,b,c）の 12 辺を太線で描く。
export function makeCellEdges(vectors, color, offset = new THREE.Vector3(), opacity = 1, dashed = false) {
  const [a, b, c] = vectors;
  const corners = [];
  for (let i = 0; i < 2; i++)
    for (let j = 0; j < 2; j++)
      for (let k = 0; k < 2; k++)
        corners.push(offset.clone().add(a.clone().multiplyScalar(i)).add(b.clone().multiplyScalar(j)).add(c.clone().multiplyScalar(k)));
  const idx = (i, j, k) => i * 4 + j * 2 + k;
  const edges = [];
  for (let i = 0; i < 2; i++)
    for (let j = 0; j < 2; j++)
      for (let k = 0; k < 2; k++) {
        const p = corners[idx(i, j, k)];
        if (i === 0) edges.push(p, corners[idx(1, j, k)]);
        if (j === 0) edges.push(p, corners[idx(i, 1, k)]);
        if (k === 0) edges.push(p, corners[idx(i, j, 1)]);
      }
  const group = new THREE.Group();
  group.add(makeBoldLine(edges, color, opacity, dashed));
  return group;
}
