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

export function softenedMaterialColor(value) {
  const color = Number(value) >>> 0;
  const red = (color >> 16) & 0xff;
  const green = (color >> 8) & 0xff;
  const blue = color & 0xff;
  const neutral = red * 0.2126 + green * 0.7152 + blue * 0.0722;
  const channel = component => Math.max(
    0,
    Math.min(255, Math.round((neutral + (component - neutral) * 0.86) * 0.90)),
  );
  return (channel(red) << 16) | (channel(green) << 8) | channel(blue);
}

export function atomRadius(atom, fallback) {
  const radius = elementStyles.get(atomElement(atom))?.radius;
  return Number.isFinite(radius) ? radius * ATOM_RADIUS_SCALE : fallback;
}

export function colorCss(value) {
  return `#${value.toString(16).padStart(6, "0")}`;
}

export function atomMaterial(atom) {
  const color = new THREE.Color(softenedMaterialColor(atomColor(atom)));
  return new THREE.MeshPhongMaterial({
    color,
    specular: color.clone().lerp(new THREE.Color(0xffffff), 0.27),
    shininess: 28,
    emissive: color.clone().multiplyScalar(0.02),
  });
}

// 原子球。atom は {element|type_symbol} を持つ。位置は呼び出し側で設定する。
export function makeAtomMesh(atom, fallbackRadius) {
  return new THREE.Mesh(
    new THREE.SphereGeometry(atomRadius(atom, fallbackRadius), 32, 24),
    atomMaterial(atom),
  );
}

// Diffuse ambient light keeps element colors legible; one directional key
// produces one readable specular highlight on each atom.
export function addLights(scene) {
  scene.add(new THREE.HemisphereLight(0xffffff, 0xb8c2d0, 1.45));
  scene.add(new THREE.AmbientLight(0xffffff, 0.42));
  const keyLight = new THREE.DirectionalLight(0xffffff, 1.65);
  const keyTarget = new THREE.Object3D();
  keyLight.target = keyTarget;
  scene.add(keyLight, keyTarget);
  return { keyLight, keyTarget };
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
export function makeBoldLine(points, color, opacity = 1, dashed = false, emphasis = 1) {
  const group = new THREE.Group();
  const offsets = [
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0.012, 0.012, 0),
    new THREE.Vector3(-0.012, 0.012, 0),
  ];
  if (emphasis > 1) {
    offsets.push(
      new THREE.Vector3(0.012, -0.012, 0),
      new THREE.Vector3(-0.012, -0.012, 0),
    );
  }
  for (const offset of offsets) {
    group.add(makeLine(points.map((point) => point.clone().add(offset)), color, opacity, dashed));
  }
  return group;
}

// 平行六面体（基底ベクトル a,b,c）の 12 辺を太線で描く。
export function makeCellEdges(
  vectors,
  color,
  offset = new THREE.Vector3(),
  opacity = 1,
  dashed = false,
  emphasis = 1,
) {
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
  group.add(makeBoldLine(edges, color, opacity, dashed, emphasis));
  return group;
}
