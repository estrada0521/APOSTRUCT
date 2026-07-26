import * as THREE from "three";
import { TrackballControls } from "three/addons/controls/TrackballControls.js";
import { ConvexGeometry } from "three/addons/geometries/ConvexGeometry.js";
import {
  addLights,
  atomColor,
  atomRadius,
  loadElementStyles,
  makeCellEdges,
  softenedMaterialColor,
} from "./crystal3d.js?v=7";

const EPS = 1e-9;
// Bond and coordination-polyhedron rules from the active VESTA 3.5.0
// style/default.ini SBOND section (914 entries).
// Columns: [A1, A2, minA, maxA, showPolyhedra(0/1), boundaryMode(0/1/2)]
// boundaryMode 1 = COMPLETE_A1 (default); 0 = INSIDE_ONLY; 2 = RECURSIVE_VISIBLE
const VESTA_SBOND = [
  ["Ac","O",0.0,2.7326,1,1],
  ["Ac","F",0.0,2.58646,1,1],
  ["Ac","Cl",0.0,3.08646,1,1],
  ["Ac","Br",0.0,3.22726,1,1],
  ["Ag","O",0.0,2.81139,1,1],
  ["Ag","S",0.0,3.08839,1,1],
  ["Ag","F",0.0,2.76939,1,1],
  ["Ag","Cl",0.0,3.05939,1,1],
  ["Ag","Br",0.0,2.37642,1,1],
  ["Ag","I",0.0,2.53642,1,1],
  ["Ag","Se",0.0,3.2,1,1],
  ["Ag","Te",0.0,2.66642,1,1],
  ["Ag","N",0.0,2.00642,1,1],
  ["Ag","P",0.0,2.37642,1,1],
  ["Ag","As",0.0,2.45642,1,1],
  ["Ag","H",0.0,1.65642,1,1],
  ["Al","O",0.0,2.1074,1,1],
  ["Al","S",0.0,2.66646,1,1],
  ["Al","Se",0.0,2.72646,1,1],
  ["Al","Te",0.0,2.93646,1,1],
  ["Al","F",0.0,2.00146,1,1],
  ["Al","Cl",0.0,2.48846,1,1],
  ["Al","Br",0.0,2.65646,1,1],
  ["Al","I",0.0,2.86646,1,1],
  ["Al","N",0.0,2.24646,1,1],
  ["Al","P",0.0,2.69646,1,1],
  ["Al","As",0.0,2.75646,1,1],
  ["Al","H",0.0,1.90646,1,1],
  ["Am","O",0.0,2.71649,1,1],
  ["Am","F",0.0,2.61945,1,1],
  ["Am","Cl",0.0,3.08945,1,1],
  ["Am","Br",0.0,3.22944,1,1],
  ["As","S",0.0,2.84649,1,1],
  ["As","Se",0.0,2.98649,1,1],
  ["As","O",0.0,2.24546,1,1],
  ["As","Te",0.0,3.10646,1,1],
  ["As","F",0.0,2.15646,1,1],
  ["As","Cl",0.0,2.61646,1,1],
  ["As","Br",0.0,2.80646,1,1],
  ["As","I",0.0,3.03646,1,1],
  ["As","C",0.0,2.38646,1,1],
  ["Au","Cl",0.0,2.88295,1,1],
  ["Au","I",0.0,3.21295,1,1],
  ["Au","O",0.0,2.34646,1,1],
  ["Au","S",0.0,2.8326,1,1],
  ["Au","F",0.0,2.34646,1,1],
  ["Au","Br",0.0,2.77646,1,1],
  ["Au","N",0.0,2.3826,1,1],
  ["Au","Se",0.0,2.22998,1,1],
  ["Au","Te",0.0,2.45998,1,1],
  ["Au","P",0.0,2.18998,1,1],
  ["Au","As",0.0,2.26998,1,1],
  ["Au","H",0.0,1.41998,1,1],
  ["B","O",0.0,1.82746,1,1],
  ["B","S",0.0,2.27646,1,1],
  ["B","Se",0.0,2.40646,1,1],
  ["B","Te",0.0,2.65646,1,1],
  ["B","F",0.0,1.76646,1,1],
  ["B","Cl",0.0,2.19646,1,1],
  ["B","Br",0.0,2.33646,1,1],
  ["B","I",0.0,2.55646,1,1],
  ["B","N",0.0,1.93846,1,1],
  ["B","P",0.0,2.37646,1,1],
  ["B","As",0.0,2.42646,1,1],
  ["B","H",0.0,1.59646,1,1],
  ["B","B",0.0,1.85846,1,1],
  ["Ba","O",0.0,3.14795,1,1],
  ["Ba","S",0.0,3.66195,1,1],
  ["Ba","Se",0.0,3.74295,1,1],
  ["Ba","Te",0.0,3.94295,1,1],
  ["Ba","F",0.0,3.05095,1,1],
  ["Ba","Cl",0.0,3.55295,1,1],
  ["Ba","Br",0.0,3.74295,1,1],
  ["Ba","I",0.0,3.99295,1,1],
  ["Ba","N",0.0,3.33295,1,1],
  ["Ba","P",0.0,3.48295,1,1],
  ["Ba","As",0.0,3.69,1,1],
  ["Ba","H",0.0,3.08295,1,1],
  ["Be","O",0.0,1.98749,1,1],
  ["Be","S",0.0,2.43649,1,1],
  ["Be","Se",0.0,2.57649,1,1],
  ["Be","Te",0.0,2.81649,1,1],
  ["Be","F",0.0,1.88749,1,1],
  ["Be","Cl",0.0,2.36649,1,1],
  ["Be","Br",0.0,2.50649,1,1],
  ["Be","I",0.0,2.70649,1,1],
  ["Be","N",0.0,2.10649,1,1],
  ["Be","P",0.0,2.55649,1,1],
  ["Be","As",0.0,2.60649,1,1],
  ["Be","H",0.0,1.71649,1,1],
  ["Bi","O",0.0,2.6608,1,1],
  ["Bi","S",0.0,3.13291,1,1],
  ["Bi","Se",0.0,3.5,1,1],
  ["Bi","F",0.0,2.55291,1,1],
  ["Bi","Cl",0.0,3.04291,1,1],
  ["Bi","Br",0.0,3.17993,1,1],
  ["Bi","I",0.0,3.38291,1,1],
  ["Bi","N",0.0,2.56329,1,1],
  ["Bi","Te",0.0,3.02642,1,1],
  ["Bi","P",0.0,2.78642,1,1],
  ["Bi","As",0.0,2.87642,1,1],
  ["Bi","H",0.0,2.12642,1,1],
  ["Bk","O",0.0,2.64329,1,1],
  ["Bk","F",0.0,2.54233,1,1],
  ["Bk","Cl",0.0,3.02291,1,1],
  ["Bk","Br",0.0,3.15233,1,1],
  ["Br","O",0.0,2.35646,1,1],
  ["Br","F",0.0,2.20646,1,1],
  ["Br","Cl",0.0,2.33296,1,1],
  ["C","O",0.0,1.97249,1,2],
  ["C","Cl",0.0,2.11002,1,2],
  ["C","C",0.0,1.89002,0,2],
  ["C","S",0.0,2.15002,0,2],
  ["C","F",0.0,1.76002,1,2],
  ["C","Br",0.0,2.26002,1,2],
  ["C","N",0.0,1.79202,1,2],
  ["C","Se",0.0,2.01998,1,2],
  ["C","I",0.0,2.16998,1,2],
  ["C","Te",0.0,2.25998,1,1],
  ["C","P",0.0,1.93998,1,2],
  ["C","H",0.0,1.2,0,1],
  ["Ca","O",0.0,2.83062,1,1],
  ["Ca","S",0.0,3.31295,1,1],
  ["Ca","Se",0.0,3.42295,1,1],
  ["Ca","Te",0.0,3.62295,1,1],
  ["Ca","F",0.0,2.70495,1,1],
  ["Ca","Cl",0.0,3.23295,1,1],
  ["Ca","Br",0.0,3.36995,1,1],
  ["Ca","I",0.0,3.58295,1,1],
  ["Ca","N",0.0,3.00295,1,1],
  ["Ca","P",0.0,3.41295,1,1],
  ["Ca","As",0.0,3.48295,1,1],
  ["Ca","H",0.0,2.69295,1,1],
  ["Cd","O",0.0,2.76695,1,1],
  ["Cd","S",0.0,3.16695,1,1],
  ["Cd","Se",0.0,3.26295,1,1],
  ["Cd","Te",0.0,3.45295,1,1],
  ["Cd","F",0.0,2.67395,1,1],
  ["Cd","Cl",0.0,3.09295,1,1],
  ["Cd","Br",0.0,3.21295,1,1],
  ["Cd","I",0.0,3.46295,1,1],
  ["Cd","N",0.0,2.82295,1,1],
  ["Cd","P",0.0,3.20295,1,1],
  ["Cd","As",0.0,3.29295,1,1],
  ["Cd","H",0.0,2.52295,1,1],
  ["Ce","O",0.0,2.86393,1,1],
  ["Ce","S",0.0,3.36293,1,1],
  ["Ce","F",0.0,2.75452,1,1],
  ["Ce","Cl",0.0,3.25293,1,1],
  ["Ce","Br",0.0,3.40452,1,1],
  ["Ce","I",0.0,3.62452,1,1],
  ["Ce","N",0.0,2.78549,1,1],
  ["Ce","Se",0.0,3.04644,1,1],
  ["Ce","Te",0.0,3.22644,1,1],
  ["Ce","P",0.0,3.00644,1,1],
  ["Ce","As",0.0,3.08644,1,1],
  ["Ce","H",0.0,2.34644,1,1],
  ["Cf","O",0.0,2.63291,1,1],
  ["Cf","F",0.0,2.53233,1,1],
  ["Cf","Cl",0.0,3.01291,1,1],
  ["Cf","Br",0.0,3.14233,1,1],
  ["Cl","O",0.0,2.16646,1,1],
  ["Cl","F",0.0,2.14646,1,1],
  ["Cl","Cl",0.0,2.14296,1,1],
  ["Cm","O",0.0,2.79291,1,1],
  ["Cm","F",0.0,2.68291,1,1],
  ["Cm","Cl",0.0,3.18291,1,1],
  ["Co","H",0.0,1.9278,1,1],
  ["Co","O",0.0,2.40493,1,1],
  ["Co","S",0.0,2.65293,1,1],
  ["Co","F",0.0,2.35293,1,1],
  ["Co","Cl",0.0,2.74593,1,1],
  ["Co","N",0.0,2.36293,1,1],
  ["Co","C",0.0,2.19691,1,1],
  ["Co","Br",0.0,2.33642,1,1],
  ["Co","I",0.0,2.52878,1,1],
  ["Co","Se",0.0,2.39642,1,1],
  ["Co","Te",0.0,2.61642,1,1],
  ["Co","P",0.0,2.36642,1,1],
  ["Co","As",0.0,2.43642,1,1],
  ["Cr","O",0.0,2.44293,1,1],
  ["Cr","F",0.0,2.45293,1,1],
  ["Cr","Cl",0.0,2.80293,1,1],
  ["Cr","Br",0.0,2.97293,1,1],
  ["Cr","I",0.0,3.19293,1,1],
  ["Cr","N",0.0,2.5152,1,1],
  ["Cr","S",0.0,2.72491,1,1],
  ["Cr","Se",0.0,2.44642,1,1],
  ["Cr","Te",0.0,2.67642,1,1],
  ["Cr","P",0.0,2.42642,1,1],
  ["Cr","As",0.0,2.49642,1,1],
  ["Cr","H",0.0,1.67642,1,1],
  ["Cs","O",0.0,3.53642,1,1],
  ["Cs","S",0.0,4.24942,1,1],
  ["Cs","Se",0.0,4.21655,1,1],
  ["Cs","Te",0.0,4.4631,1,1],
  ["Cs","F",0.0,3.49942,1,1],
  ["Cs","Cl",0.0,3.91042,1,1],
  ["Cs","Br",0.0,4.06942,1,1],
  ["Cs","I",0.0,4.40942,1,1],
  ["Cs","N",0.0,3.94942,1,1],
  ["Cs","P",0.0,3.64194,1,1],
  ["Cs","As",0.0,4.15942,1,1],
  ["Cs","H",0.0,3.55942,1,1],
  ["Cu","O",0.0,2.47295,1,1],
  ["Cu","S",0.0,2.76095,1,1],
  ["Cu","Se",0.0,2.76295,1,1],
  ["Cu","F",0.0,2.46295,1,1],
  ["Cu","Cl",0.0,2.75295,1,1],
  ["Cu","Br",0.0,2.89295,1,1],
  ["Cu","I",0.0,3.01795,1,1],
  ["Cu","N",0.0,2.49295,1,1],
  ["Cu","P",0.0,2.65649,1,1],
  ["Cu","As",0.0,2.71895,1,1],
  ["Cu","C",0.0,2.32649,1,1],
  ["Cu","Te",0.0,2.87649,1,1],
  ["Cu","H",0.0,1.81649,1,1],
  ["Dy","O",0.0,2.65651,1,1],
  ["Dy","F",0.0,2.52944,1,1],
  ["Dy","Cl",0.0,3.01945,1,1],
  ["Dy","Br",0.0,3.16945,1,1],
  ["Dy","I",0.0,3.39945,1,1],
  ["Dy","S",0.0,2.823,1,1],
  ["Dy","Se",0.0,2.81,1,1],
  ["Dy","Te",0.0,3.22,1,1],
  ["Dy","N",0.0,2.38,1,1],
  ["Dy","P",0.0,2.77,1,1],
  ["Dy","As",0.0,3.14,1,1],
  ["Dy","H",0.0,2.09,1,1],
  ["Er","O",0.0,2.63651,1,1],
  ["Er","S",0.0,3.27651,1,1],
  ["Er","Se",0.0,3.18649,1,1],
  ["Er","F",0.0,2.51049,1,1],
  ["Er","Cl",0.0,2.99944,1,1],
  ["Er","Br",0.0,3.14945,1,1],
  ["Er","I",0.0,3.38945,1,1],
  ["Er","Te",0.0,3.22,1,1],
  ["Er","N",0.0,2.36,1,1],
  ["Er","P",0.0,2.75,1,1],
  ["Er","As",0.0,3.18,1,1],
  ["Er","H",0.0,2.06,1,1],
  ["Es","O",0.0,2.70139,1,1],
  ["Eu","O",0.0,2.94249,1,1],
  ["Eu","S",0.0,3.37949,1,1],
  ["Eu","F",0.0,2.83549,1,1],
  ["Eu","Cl",0.0,3.32549,1,1],
  ["Eu","Br",0.0,3.46549,1,1],
  ["Eu","I",0.0,3.69549,1,1],
  ["Eu","N",0.0,2.95549,1,1],
  ["Eu","Se",0.0,2.89898,1,1],
  ["Eu","Te",0.0,3.08898,1,1],
  ["Eu","P",0.0,2.85898,1,1],
  ["Eu","As",0.0,2.93898,1,1],
  ["Eu","H",0.0,2.18898,1,1],
  ["Fe","O",0.0,2.44693,1,1],
  ["Fe","S",0.0,2.83793,1,1],
  ["Fe","F",0.0,2.36293,1,1],
  ["Fe","Cl",0.0,2.86293,1,1],
  ["Fe","Br",0.0,2.8952,1,1],
  ["Fe","I",0.0,3.1552,1,1],
  ["Fe","N",0.0,2.48193,1,1],
  ["Fe","C",0.0,2.25191,1,1],
  ["Fe","Se",0.0,2.43642,1,1],
  ["Fe","Te",0.0,2.68642,1,1],
  ["Fe","P",0.0,2.42642,1,1],
  ["Fe","As",0.0,2.50642,1,1],
  ["Fe","H",0.0,1.68642,1,1],
  ["Ga","Se",0.0,3.41295,1,1],
  ["Ga","O",0.0,2.18646,1,1],
  ["Ga","S",0.0,2.61946,1,1],
  ["Ga","F",0.0,2.14646,1,1],
  ["Ga","Cl",0.0,2.52646,1,1],
  ["Ga","Br",0.0,2.6426,1,1],
  ["Ga","I",0.0,2.91646,1,1],
  ["Ga","Te",0.0,2.58998,1,1],
  ["Ga","N",0.0,1.96,1,1],
  ["Ga","P",0.0,2.46998,1,1],
  ["Ga","As",0.0,2.38998,1,1],
  ["Ga","H",0.0,1.55998,1,1],
  ["Gd","O",0.0,2.76651,1,1],
  ["Gd","F",0.0,3.15651,1,1],
  ["Gd","S",0.0,3.13649,1,1],
  ["Gd","Cl",0.0,3.06349,1,1],
  ["Gd","Br",0.0,3.19945,1,1],
  ["Gd","I",0.0,3.41945,1,1],
  ["Gd","Se",0.0,2.85,1,1],
  ["Gd","Te",0.0,3.24,1,1],
  ["Gd","N",0.0,2.42,1,1],
  ["Gd","P",0.0,2.81,1,1],
  ["Gd","As",0.0,3.18,1,1],
  ["Gd","H",0.0,2.13,1,1],
  ["Ge","O",0.0,2.09802,1,1],
  ["Ge","S",0.0,2.56702,1,1],
  ["Ge","Se",0.0,2.70002,1,1],
  ["Ge","F",0.0,2.01002,1,1],
  ["Ge","Cl",0.0,2.49002,1,1],
  ["Ge","Br",0.0,2.34998,1,1],
  ["Ge","I",0.0,2.54998,1,1],
  ["Ge","Te",0.0,2.60998,1,1],
  ["Ge","N",0.0,1.92998,1,1],
  ["Ge","P",0.0,2.36998,1,1],
  ["Ge","As",0.0,2.47998,1,1],
  ["Ge","H",0.0,1.59998,1,1],
  ["Ge","Ge",0.0,2.6,1,1],
  ["O","H",0.0,1.2,0,1],
  ["H","O",1.2,2.1,0,0],
  ["H","F",0.0,1.1,0,1],
  ["H","Cl",0.0,1.5,0,1],
  ["H","N",0.0,1.2,0,1],
  ["O","D",0.0,1.2,0,1],
  ["D","O",1.2,2.1,0,0],
  ["D","F",0.0,1.1,0,1],
  ["D","Cl",0.0,1.5,0,1],
  ["D","N",0.0,1.2,0,1],
  ["Hf","F",0.0,3.18291,1,1],
  ["Hf","O",0.0,2.37946,1,1],
  ["Hf","Cl",0.0,2.75646,1,1],
  ["Hf","Br",0.0,2.62642,1,1],
  ["Hf","S",0.0,2.64642,1,1],
  ["Hf","Se",0.0,2.67642,1,1],
  ["Hf","Te",0.0,2.87642,1,1],
  ["Hf","I",0.0,2.83642,1,1],
  ["Hf","N",0.0,2.24642,1,1],
  ["Hf","P",0.0,2.63642,1,1],
  ["Hf","As",0.0,2.71642,1,1],
  ["Hf","H",0.0,1.93642,1,1],
  ["Hg","O",0.0,2.86939,1,1],
  ["Hg","F",0.0,2.88293,1,1],
  ["Hg","Cl",0.0,3.24939,1,1],
  ["Hg","S",0.0,3.42093,1,1],
  ["Hg","Br",0.0,3.09293,1,1],
  ["Hg","I",0.0,3.33293,1,1],
  ["Hg","Se",0.0,2.82642,1,1],
  ["Hg","Te",0.0,2.76642,1,1],
  ["Hg","N",0.0,2.17642,1,1],
  ["Hg","P",0.0,2.57642,1,1],
  ["Hg","As",0.0,2.65642,1,1],
  ["Hg","H",0.0,1.86642,1,1],
  ["Hg","Hg",0.0,3.1952,1,1],
  ["Ho","O",0.0,2.67047,1,1],
  ["Ho","S",0.0,3.13547,1,1],
  ["Ho","F",0.0,2.56159,1,1],
  ["Ho","Cl",0.0,3.05159,1,1],
  ["Ho","Br",0.0,3.20159,1,1],
  ["Ho","I",0.0,3.44159,1,1],
  ["Ho","Se",0.0,2.84898,1,1],
  ["Ho","Te",0.0,3.03898,1,1],
  ["Ho","N",0.0,2.41898,1,1],
  ["Ho","P",0.0,2.79898,1,1],
  ["Ho","As",0.0,2.87898,1,1],
  ["Ho","H",0.0,2.11898,1,1],
  ["I","I",0.0,2.4,1,1],
  ["I","F",0.0,3.18295,1,1],
  ["I","Cl",0.0,3.33295,1,1],
  ["I","O",0.0,2.47646,1,1],
  ["In","Cl",0.0,3.52939,1,1],
  ["In","O",0.0,2.46491,1,1],
  ["In","S",0.0,2.93291,1,1],
  ["In","F",0.0,2.35491,1,1],
  ["In","Br",0.0,3.05329,1,1],
  ["In","I",0.0,3.19291,1,1],
  ["In","Co",0.0,3.13629,1,1],
  ["In","Mn",0.0,3.14729,1,1],
  ["In","Se",0.0,2.62642,1,1],
  ["In","Te",0.0,2.84642,1,1],
  ["In","N",0.0,2.18642,1,1],
  ["In","P",0.0,2.68642,1,1],
  ["In","As",0.0,2.66642,1,1],
  ["In","H",0.0,1.87642,1,1],
  ["Ir","O",0.0,2.27746,1,1],
  ["Ir","F",0.0,2.15002,1,1],
  ["Ir","Cl",0.0,2.56746,1,1],
  ["Ir","S",0.0,2.42998,1,1],
  ["Ir","Se",0.0,2.55998,1,1],
  ["Ir","Te",0.0,2.75998,1,1],
  ["Ir","Br",0.0,2.49998,1,1],
  ["Ir","I",0.0,2.70998,1,1],
  ["Ir","N",0.0,2.10998,1,1],
  ["Ir","P",0.0,2.50998,1,1],
  ["Ir","As",0.0,2.58998,1,1],
  ["Ir","H",0.0,1.80998,1,1],
  ["K","O",0.0,3.25142,1,1],
  ["K","S",0.0,3.79285,1,1],
  ["K","Se",0.0,4.00186,1,1],
  ["K","Te",0.0,4.23284,1,1],
  ["K","F",0.0,3.11142,1,1],
  ["K","Cl",0.0,3.65976,1,1],
  ["K","Br",0.0,3.8513,1,1],
  ["K","I",0.0,4.11717,1,1],
  ["K","N",0.0,3.41942,1,1],
  ["K","P",0.0,3.44942,1,1],
  ["K","As",0.0,3.94942,1,1],
  ["K","H",0.0,3.21942,1,1],
  ["Kr","F",0.0,2.67549,1,1],
  ["La","O",0.0,2.90983,1,1],
  ["La","S",0.0,3.35593,1,1],
  ["La","Se",0.0,3.45293,1,1],
  ["La","Te",0.0,3.65293,1,1],
  ["La","F",0.0,2.79293,1,1],
  ["La","Cl",0.0,3.33452,1,1],
  ["La","Br",0.0,3.43293,1,1],
  ["La","I",0.0,3.64293,1,1],
  ["La","N",0.0,3.05293,1,1],
  ["La","P",0.0,3.32293,1,1],
  ["La","As",0.0,3.51293,1,1],
  ["La","H",0.0,2.77293,1,1],
  ["Li","O",0.0,2.60087,1,1],
  ["Li","S",0.0,3.02481,1,1],
  ["Li","Se",0.0,3.2433,1,1],
  ["Li","Te",0.0,3.42496,1,1],
  ["Li","F",0.0,2.34276,1,1],
  ["Li","Cl",0.0,2.91814,1,1],
  ["Li","Br",0.0,3.11654,1,1],
  ["Li","I",0.0,3.37676,1,1],
  ["Li","N",0.0,2.66213,1,1],
  ["Lu","O",0.0,2.57749,1,1],
  ["Lu","S",0.0,3.03649,1,1],
  ["Lu","Se",0.0,3.16649,1,1],
  ["Lu","Te",0.0,3.35649,1,1],
  ["Lu","F",0.0,2.48249,1,1],
  ["Lu","Cl",0.0,2.96944,1,1],
  ["Lu","Br",0.0,3.11945,1,1],
  ["Lu","I",0.0,3.36945,1,1],
  ["Lu","N",0.0,2.71649,1,1],
  ["Lu","P",0.0,3.11649,1,1],
  ["Lu","As",0.0,3.19649,1,1],
  ["Lu","H",0.0,2.42649,1,1],
  ["Mg","O",0.0,2.41824,1,1],
  ["Mg","S",0.0,2.89293,1,1],
  ["Mg","Se",0.0,3.03293,1,1],
  ["Mg","Te",0.0,3.24293,1,1],
  ["Mg","F",0.0,2.29093,1,1],
  ["Mg","Cl",0.0,2.79293,1,1],
  ["Mg","Br",0.0,2.99293,1,1],
  ["Mg","I",0.0,3.17293,1,1],
  ["Mg","N",0.0,2.56293,1,1],
  ["Mg","P",0.0,3.00293,1,1],
  ["Mg","As",0.0,3.09293,1,1],
  ["Mg","H",0.0,2.24293,1,1],
  ["Mn","O",0.0,2.51652,1,1],
  ["Mn","S",0.0,2.93293,1,1],
  ["Mn","F",0.0,2.41093,1,1],
  ["Mn","Cl",0.0,2.84593,1,1],
  ["Mn","Br",0.0,3.05293,1,1],
  ["Mn","I",0.0,3.23293,1,1],
  ["Mn","N",0.0,2.56193,1,1],
  ["Mn","Se",0.0,2.47642,1,1],
  ["Mn","Te",0.0,2.70642,1,1],
  ["Mn","P",0.0,2.39642,1,1],
  ["Mn","As",0.0,2.51642,1,1],
  ["Mn","H",0.0,1.70642,1,1],
  ["Mo","S",0.0,2.80067,1,1],
  ["Mo","Cl",0.0,2.80447,1,1],
  ["Mo","O",0.0,2.3475,1,1],
  ["Mo","F",0.0,2.2998,1,1],
  ["Mo","Br",0.0,2.8535,1,1],
  ["Mo","N",0.0,2.4735,1,1],
  ["Mo","I",0.0,2.74701,1,1],
  ["Mo","Se",0.0,2.59701,1,1],
  ["Mo","Te",0.0,2.79701,1,1],
  ["Mo","P",0.0,2.54701,1,1],
  ["Mo","As",0.0,2.62701,1,1],
  ["Mo","H",0.0,1.83701,1,1],
  ["N","O",0.0,1.81746,1,1],
  ["N","F",0.0,1.82646,1,1],
  ["N","Cl",0.0,2.20646,1,1],
  ["N","N",0.0,1.8826,1,1],
  ["Na","O",0.0,2.95693,1,1],
  ["Na","S",0.0,3.57685,1,1],
  ["Na","Se",0.0,3.71593,1,1],
  ["Na","Te",0.0,3.95459,1,1],
  ["Na","F",0.0,2.80398,1,1],
  ["Na","Cl",0.0,3.39412,1,1],
  ["Na","Br",0.0,3.57715,1,1],
  ["Na","I",0.0,3.88251,1,1],
  ["Na","N",0.0,3.12942,1,1],
  ["Na","P",0.0,3.47942,1,1],
  ["Na","As",0.0,3.64942,1,1],
  ["Na","H",0.0,2.79942,1,1],
  ["Nb","O",0.0,2.45329,1,1],
  ["Nb","F",0.0,2.35646,1,1],
  ["Nb","Cl",0.0,2.76291,1,1],
  ["Nb","Br",0.0,3.07646,1,1],
  ["Nb","N",0.0,2.46046,1,1],
  ["Nb","I",0.0,3.1439,1,1],
  ["Nb","S",0.0,2.74,1,1],
  ["Nb","Se",0.0,2.66642,1,1],
  ["Nb","Te",0.0,2.85642,1,1],
  ["Nb","P",0.0,2.61642,1,1],
  ["Nb","As",0.0,2.69642,1,1],
  ["Nb","H",0.0,1.90642,1,1],
  ["Nd","O",0.0,2.8587,1,1],
  ["Nd","S",0.0,3.42712,1,1],
  ["Nd","Se",0.0,3.42293,1,1],
  ["Nd","Te",0.0,3.60293,1,1],
  ["Nd","F",0.0,2.73452,1,1],
  ["Nd","Cl",0.0,3.22493,1,1],
  ["Nd","Br",0.0,3.37293,1,1],
  ["Nd","I",0.0,3.59452,1,1],
  ["Nd","N",0.0,3.01293,1,1],
  ["NH","O",0.0,3.08895,1,1],
  ["NH","F",0.0,2.99195,1,1],
  ["NH","Cl",0.0,3.48195,1,1],
  ["Ni","O",0.0,2.28149,1,1],
  ["Ni","S",0.0,2.58649,1,1],
  ["Ni","F",0.0,2.20249,1,1],
  ["Ni","Cl",0.0,2.62649,1,1],
  ["Ni","Br",0.0,2.80649,1,1],
  ["Ni","I",0.0,3.00649,1,1],
  ["Ni","N",0.0,2.30649,1,1],
  ["Ni","Se",0.0,2.18998,1,1],
  ["Ni","Te",0.0,2.47998,1,1],
  ["Ni","P",0.0,2.31998,1,1],
  ["Ni","As",0.0,2.28998,1,1],
  ["Ni","H",0.0,1.44998,1,1],
  ["Np","F",0.0,2.59233,1,1],
  ["Np","Cl",0.0,3.07233,1,1],
  ["Np","S",0.0,3.2,1,1],
  ["Np","Br",0.0,3.21233,1,1],
  ["Np","I",0.0,3.44233,1,1],
  ["Np","O",0.0,2.63646,1,1],
  ["O","O",0.0,1.7,0,1],
  ["Os","O",0.0,2.23,1,1],
  ["Os","S",0.0,2.56002,1,1],
  ["Os","F",0.0,2.07746,1,1],
  ["Os","Cl",0.0,2.54002,1,1],
  ["Os","Br",0.0,2.72002,1,1],
  ["P","O",0.0,2.08646,1,1],
  ["P","S",0.0,2.57646,1,1],
  ["P","Se",0.0,2.69646,1,1],
  ["P","F",0.0,2.01002,1,1],
  ["P","Cl",0.0,2.28746,1,1],
  ["P","Br",0.0,2.44293,1,1],
  ["P","N",0.0,1.97146,1,1],
  ["P","I",0.0,2.44998,1,1],
  ["P","P",0.0,2.48381,1,1],
  ["P","As",0.0,2.29998,1,1],
  ["P","H",0.0,1.45998,1,1],
  ["Pa","O",0.0,2.63383,1,1],
  ["Pa","F",0.0,2.54437,1,1],
  ["Pa","Cl",0.0,3.01437,1,1],
  ["Pa","Br",0.0,3.18437,1,1],
  ["Pb","O",0.0,3.04096,1,1],
  ["Pb","S",0.0,3.40395,1,1],
  ["Pb","Se",0.0,3.55295,1,1],
  ["Pb","F",0.0,2.92045,1,1],
  ["Pb","Cl",0.0,3.39295,1,1],
  ["Pb","Br",0.0,3.62451,1,1],
  ["Pb","I",0.0,3.69562,1,1],
  ["Pb","N",0.0,3.0967,1,1],
  ["Pb","Te",0.0,3.14644,1,1],
  ["Pb","P",0.0,2.94644,1,1],
  ["Pb","As",0.0,3.02644,1,1],
  ["Pb","H",0.0,2.27644,1,1],
  ["Pd","O",0.0,2.39849,1,1],
  ["Pd","S",0.0,2.69649,1,1],
  ["Pd","F",0.0,2.34649,1,1],
  ["Pd","Cl",0.0,2.65649,1,1],
  ["Pd","Br",0.0,2.80649,1,1],
  ["Pd","I",0.0,2.96649,1,1],
  ["Pd","N",0.0,2.40451,1,1],
  ["Pd","C",0.0,2.33649,1,1],
  ["Pd","Se",0.0,2.26998,1,1],
  ["Pd","Te",0.0,2.52998,1,1],
  ["Pd","P",0.0,2.46998,1,1],
  ["Pd","As",0.0,2.34998,1,1],
  ["Pd","H",0.0,1.51998,1,1],
  ["Pm","F",0.0,2.55233,1,1],
  ["Pm","Cl",0.0,3.41233,1,1],
  ["Pm","Br",0.0,3.18233,1,1],
  ["Po","O",0.0,2.64646,1,1],
  ["Po","F",0.0,2.83646,1,1],
  ["Pr","O",0.0,2.74449,1,1],
  ["Pr","S",0.0,3.20649,1,1],
  ["Pr","Se",0.0,3.32649,1,1],
  ["Pr","Te",0.0,3.50649,1,1],
  ["Pr","F",0.0,2.62945,1,1],
  ["Pr","Cl",0.0,3.12749,1,1],
  ["Pr","Br",0.0,3.27649,1,1],
  ["Pr","I",0.0,3.49649,1,1],
  ["Pr","N",0.0,2.90649,1,1],
  ["Pr","P",0.0,3.28649,1,1],
  ["Pr","As",0.0,3.35649,1,1],
  ["Pr","H",0.0,2.62649,1,1],
  ["Pt","O",0.0,2.40649,1,1],
  ["Pt","S",0.0,2.76649,1,1],
  ["Pt","F",0.0,2.54002,1,1],
  ["Pt","Cl",0.0,2.75646,1,1],
  ["Pt","Br",0.0,2.94191,1,1],
  ["Pt","C",0.0,2.36649,1,1],
  ["Pt","N",0.0,2.41649,1,1],
  ["Pt","I",0.0,2.41998,1,1],
  ["Pt","Se",0.0,2.23998,1,1],
  ["Pt","Te",0.0,2.49998,1,1],
  ["Pt","P",0.0,2.23998,1,1],
  ["Pt","As",0.0,2.30998,1,1],
  ["Pt","H",0.0,1.44998,1,1],
  ["Pu","O",0.0,2.68329,1,1],
  ["Pu","F",0.0,2.58233,1,1],
  ["Pu","Cl",0.0,3.05233,1,1],
  ["Pu","S",0.0,3.3,1,1],
  ["Pu","Br",0.0,3.19233,1,1],
  ["Pu","I",0.0,3.43233,1,1],
  ["Rb","O",0.0,3.43945,1,1],
  ["Rb","S",0.0,3.97645,1,1],
  ["Rb","Se",0.0,4.13773,1,1],
  ["Rb","Te",0.0,4.28802,1,1],
  ["Rb","F",0.0,3.37645,1,1],
  ["Rb","Cl",0.0,3.86664,1,1],
  ["Rb","Br",0.0,4.05498,1,1],
  ["Rb","I",0.0,4.33462,1,1],
  ["Rb","N",0.0,3.79645,1,1],
  ["Rb","P",0.0,3.50645,1,1],
  ["Rb","As",0.0,4.04645,1,1],
  ["Rb","H",0.0,3.43645,1,1],
  ["Re","Cl",0.0,3.44712,1,1],
  ["Re","O",0.0,2.3426,1,1],
  ["Re","F",0.0,2.16002,1,1],
  ["Re","Br",0.0,2.70002,1,1],
  ["Re","I",0.0,2.65998,1,1],
  ["Re","S",0.0,2.61998,1,1],
  ["Re","Se",0.0,2.54998,1,1],
  ["Re","Te",0.0,2.74998,1,1],
  ["Re","N",0.0,2.10998,1,1],
  ["Re","P",0.0,2.50998,1,1],
  ["Re","As",0.0,2.61998,1,1],
  ["Re","H",0.0,1.79998,1,1],
  ["Rh","O",0.0,2.24946,1,1],
  ["Rh","F",0.0,2.16646,1,1],
  ["Rh","Cl",0.0,2.62646,1,1],
  ["Rh","Br",0.0,2.7126,1,1],
  ["Rh","N",0.0,2.2626,1,1],
  ["Rh","I",0.0,2.52998,1,1],
  ["Rh","S",0.0,2.19998,1,1],
  ["Rh","Se",0.0,2.37998,1,1],
  ["Rh","Te",0.0,2.59998,1,1],
  ["Rh","P",0.0,2.43998,1,1],
  ["Rh","As",0.0,2.41998,1,1],
  ["Rh","H",0.0,1.59998,1,1],
  ["Ru","Se",0.0,2.69451,1,1],
  ["Ru","F",0.0,2.57646,1,1],
  ["Ru","O",0.0,2.22646,1,1],
  ["Ru","S",0.0,2.6426,1,1],
  ["Ru","Cl",0.0,2.70646,1,1],
  ["Ru","N",0.0,2.2626,1,1],
  ["Ru","Br",0.0,2.30998,1,1],
  ["Ru","I",0.0,2.52998,1,1],
  ["Ru","Te",0.0,2.58998,1,1],
  ["Ru","P",0.0,2.33998,1,1],
  ["Ru","As",0.0,2.40998,1,1],
  ["Ru","H",0.0,1.65998,1,1],
  ["S","O",0.0,1.8,1,1],
  ["S","S",0.0,2.2,1,1],
  ["S","N",0.0,2.28849,1,1],
  ["S","F",0.0,1.95002,1,1],
  ["S","Cl",0.0,2.37002,1,1],
  ["S","Br",0.0,2.21998,1,1],
  ["S","I",0.0,2.40998,1,1],
  ["S","H",0.0,1.42998,1,1],
  ["Sb","O",0.0,2.45237,1,1],
  ["Sb","S",0.0,3.05,1,1],
  ["Sb","Se",0.0,3.05646,1,1],
  ["Sb","F",0.0,2.35646,1,1],
  ["Sb","Cl",0.0,2.80646,1,1],
  ["Sb","Br",0.0,2.96646,1,1],
  ["Sb","I",0.0,3.21646,1,1],
  ["Sb","N",0.0,2.56446,1,1],
  ["Sb","Te",0.0,2.82998,1,1],
  ["Sb","P",0.0,2.56998,1,1],
  ["Sb","As",0.0,2.64998,1,1],
  ["Sb","H",0.0,2.81998,1,1],
  ["Sc","O",0.0,2.42029,1,1],
  ["Sc","S",0.0,2.88391,1,1],
  ["Sc","Se",0.0,3.00291,1,1],
  ["Sc","Te",0.0,3.20291,1,1],
  ["Sc","F",0.0,2.32291,1,1],
  ["Sc","Cl",0.0,2.92291,1,1],
  ["Sc","Br",0.0,2.94291,1,1],
  ["Sc","I",0.0,3.15291,1,1],
  ["Sc","N",0.0,2.54291,1,1],
  ["Sc","P",0.0,2.96291,1,1],
  ["Sc","As",0.0,3.04291,1,1],
  ["Sc","H",0.0,2.24291,1,1],
  ["Se","S",0.0,2.81649,1,1],
  ["Se","Se",0.0,2.93649,1,1],
  ["Se","O",0.0,2.16102,1,1],
  ["Se","F",0.0,2.08002,1,1],
  ["Se","Cl",0.0,2.57002,1,1],
  ["Se","Br",0.0,2.78002,1,1],
  ["Se","N",0.0,2.1,1,1],
  ["Se","I",0.0,2.58998,1,1],
  ["Se","H",0.0,1.58998,1,1],
  ["Si","O",0.0,1.99002,1,1],
  ["Si","S",0.0,2.47602,1,1],
  ["Si","Se",0.0,2.61002,1,1],
  ["Si","Te",0.0,2.84002,1,1],
  ["Si","F",0.0,1.93002,1,1],
  ["Si","Cl",0.0,2.38002,1,1],
  ["Si","Br",0.0,2.55002,1,1],
  ["Si","I",0.0,2.76002,1,1],
  ["Si","C",0.0,2.23302,1,1],
  ["Si","N",0.0,2.12002,1,1],
  ["Si","P",0.0,2.58002,1,1],
  ["Si","As",0.0,2.66002,1,1],
  ["Si","H",0.0,1.82002,1,1],
  ["Si","Si",0.0,2.6,1,1],
  ["Sm","O",0.0,2.88251,1,1],
  ["Sm","N",0.0,3.02351,1,1],
  ["Sm","S",0.0,3.15649,1,1],
  ["Sm","Se",0.0,3.27649,1,1],
  ["Sm","Te",0.0,3.46649,1,1],
  ["Sm","F",0.0,2.60649,1,1],
  ["Sm","Cl",0.0,3.08749,1,1],
  ["Sm","Br",0.0,3.26649,1,1],
  ["Sm","I",0.0,3.44649,1,1],
  ["Sm","P",0.0,3.23649,1,1],
  ["Sm","As",0.0,3.30649,1,1],
  ["Sm","H",0.0,2.56649,1,1],
  ["Sn","O",0.0,2.82146,1,1],
  ["Sn","S",0.0,3.15293,1,1],
  ["Sn","F",0.0,2.63793,1,1],
  ["Sn","Cl",0.0,3.13111,1,1],
  ["Sn","Br",0.0,3.2152,1,1],
  ["Sn","I",0.0,3.52293,1,1],
  ["Sn","N",0.0,2.7152,1,1],
  ["Sn","Se",0.0,2.96646,1,1],
  ["Sn","Te",0.0,2.91642,1,1],
  ["Sn","P",0.0,2.60642,1,1],
  ["Sn","As",0.0,2.77642,1,1],
  ["Sn","H",0.0,2.00642,1,1],
  ["Sr","O",0.0,2.98095,1,1],
  ["Sr","S",0.0,3.51295,1,1],
  ["Sr","Se",0.0,3.58295,1,1],
  ["Sr","Te",0.0,3.73295,1,1],
  ["Sr","F",0.0,2.88195,1,1],
  ["Sr","Cl",0.0,3.37295,1,1],
  ["Sr","Br",0.0,3.54295,1,1],
  ["Sr","I",0.0,3.74295,1,1],
  ["Sr","N",0.0,3.09295,1,1],
  ["Sr","P",0.0,3.43295,1,1],
  ["Sr","As",0.0,3.62295,1,1],
  ["Sr","H",0.0,2.87295,1,1],
  ["Ta","O",0.0,2.74646,1,1],
  ["Ta","S",0.0,2.8439,1,1],
  ["Ta","F",0.0,2.2539,1,1],
  ["Ta","Cl",0.0,2.6739,1,1],
  ["Ta","Br",0.0,2.60642,1,1],
  ["Ta","I",0.0,2.81642,1,1],
  ["Ta","Se",0.0,2.66642,1,1],
  ["Ta","Te",0.0,2.85642,1,1],
  ["Ta","N",0.0,2.16642,1,1],
  ["Ta","P",0.0,2.62642,1,1],
  ["Ta","As",0.0,2.70642,1,1],
  ["Ta","H",0.0,1.91642,1,1],
  ["Tb","O",0.0,2.65549,1,1],
  ["Tb","S",0.0,3.11649,1,1],
  ["Tb","Se",0.0,3.23649,1,1],
  ["Tb","Te",0.0,3.42649,1,1],
  ["Tb","F",0.0,2.54249,1,1],
  ["Tb","Cl",0.0,3.04349,1,1],
  ["Tb","Br",0.0,3.18649,1,1],
  ["Tb","I",0.0,3.40945,1,1],
  ["Tb","N",0.0,2.80649,1,1],
  ["Tb","P",0.0,3.19649,1,1],
  ["Tb","As",0.0,3.26649,1,1],
  ["Tb","H",0.0,2.51649,1,1],
  ["Tc","O",0.0,2.22446,1,1],
  ["Tc","F",0.0,2.24219,1,1],
  ["Tc","Cl",0.0,2.56002,1,1],
  ["Te","O",0.0,2.3334,1,1],
  ["Te","S",0.0,2.79002,1,1],
  ["Te","F",0.0,2.22002,1,1],
  ["Te","Cl",0.0,2.73906,1,1],
  ["Te","Br",0.0,2.90002,1,1],
  ["Te","I",0.0,3.13702,1,1],
  ["Te","Se",0.0,2.57998,1,1],
  ["Te","Te",0.0,2.80998,1,1],
  ["Te","N",0.0,2.16998,1,1],
  ["Te","P",0.0,2.56998,1,1],
  ["Te","H",0.0,1.87998,1,1],
  ["Th","O",0.0,2.77349,1,1],
  ["Th","S",0.0,3.24649,1,1],
  ["Th","Se",0.0,3.36649,1,1],
  ["Th","Te",0.0,3.54649,1,1],
  ["Th","F",0.0,2.68945,1,1],
  ["Th","Cl",0.0,3.15945,1,1],
  ["Th","Br",0.0,3.31945,1,1],
  ["Th","I",0.0,3.56649,1,1],
  ["Th","N",0.0,2.94649,1,1],
  ["Th","P",0.0,3.33649,1,1],
  ["Th","As",0.0,3.40649,1,1],
  ["Th","H",0.0,2.67649,1,1],
  ["Ti","F",0.0,2.86293,1,1],
  ["Ti","Cl",0.0,3.02293,1,1],
  ["Ti","Br",0.0,3.20293,1,1],
  ["Ti","O",0.0,2.35391,1,1],
  ["Ti","S",0.0,2.74646,1,1],
  ["Ti","I",0.0,3.08291,1,1],
  ["Ti","Se",0.0,2.53642,1,1],
  ["Ti","Te",0.0,2.95642,1,1],
  ["Ti","N",0.0,2.08642,1,1],
  ["Ti","P",0.0,2.51642,1,1],
  ["Ti","As",0.0,2.87642,1,1],
  ["Ti","H",0.0,1.76642,1,1],
  ["Tl","O",0.0,3.36945,1,1],
  ["Tl","S",0.0,3.66442,1,1],
  ["Tl","F",0.0,3.26942,1,1],
  ["Tl","Cl",0.0,3.72942,1,1],
  ["Tl","Br",0.0,3.80942,1,1],
  ["Tl","I",0.0,3.94142,1,1],
  ["Tl","Se",0.0,3.00644,1,1],
  ["Tl","Te",0.0,3.23644,1,1],
  ["Tl","N",0.0,2.59644,1,1],
  ["Tl","P",0.0,3.01644,1,1],
  ["Tl","As",0.0,3.09644,1,1],
  ["Tl","H",0.0,2.35644,1,1],
  ["Tm","O",0.0,2.60649,1,1],
  ["Tm","S",0.0,3.05649,1,1],
  ["Tm","Se",0.0,3.18649,1,1],
  ["Tm","Te",0.0,3.37649,1,1],
  ["Tm","F",0.0,2.51649,1,1],
  ["Tm","Cl",0.0,2.98944,1,1],
  ["Tm","Br",0.0,3.13945,1,1],
  ["Tm","I",0.0,3.37945,1,1],
  ["Tm","N",0.0,2.74649,1,1],
  ["Tm","P",0.0,3.13649,1,1],
  ["Tm","As",0.0,3.22649,1,1],
  ["Tm","H",0.0,2.45649,1,1],
  ["U","O",0.0,2.94295,1,1],
  ["U","S",0.0,3.25293,1,1],
  ["U","F",0.0,2.80293,1,1],
  ["U","Cl",0.0,3.24452,1,1],
  ["U","Br",0.0,3.39452,1,1],
  ["U","I",0.0,3.62452,1,1],
  ["U","N",0.0,2.78649,1,1],
  ["U","Se",0.0,3.00644,1,1],
  ["U","Te",0.0,3.16644,1,1],
  ["U","P",0.0,2.94644,1,1],
  ["U","As",0.0,3.02644,1,1],
  ["U","H",0.0,2.27644,1,1],
  ["V","O",0.0,2.84939,1,1],
  ["V","Cl",0.0,3.15293,1,1],
  ["V","S",0.0,2.82293,1,1],
  ["V","F",0.0,2.87293,1,1],
  ["V","Br",0.0,2.87329,1,1],
  ["V","N",0.0,2.38329,1,1],
  ["V","I",0.0,2.66642,1,1],
  ["V","Se",0.0,2.48642,1,1],
  ["V","Te",0.0,2.72642,1,1],
  ["V","P",0.0,2.46642,1,1],
  ["V","As",0.0,2.54642,1,1],
  ["V","H",0.0,1.73642,1,1],
  ["W","O",0.0,2.20102,1,1],
  ["W","F",0.0,2.03,1,1],
  ["W","Cl",0.0,2.47,1,1],
  ["W","Br",0.0,2.49998,1,1],
  ["W","I",0.0,2.70998,1,1],
  ["W","S",0.0,2.43998,1,1],
  ["W","Se",0.0,2.55998,1,1],
  ["W","Te",0.0,2.75998,1,1],
  ["W","N",0.0,2.10998,1,1],
  ["W","P",0.0,2.50998,1,1],
  ["W","As",0.0,2.58998,1,1],
  ["W","H",0.0,1.80998,1,1],
  ["Xe","O",0.0,2.63451,1,1],
  ["Xe","F",0.0,2.62649,1,1],
  ["Y","O",0.0,2.62549,1,1],
  ["Y","S",0.0,3.08649,1,1],
  ["Y","Se",0.0,3.21649,1,1],
  ["Y","Te",0.0,3.40649,1,1],
  ["Y","F",0.0,2.51049,1,1],
  ["Y","Cl",0.0,3.00649,1,1],
  ["Y","Br",0.0,3.15649,1,1],
  ["Y","I",0.0,3.37649,1,1],
  ["Y","N",0.0,2.77649,1,1],
  ["Y","P",0.0,3.17649,1,1],
  ["Y","As",0.0,3.24649,1,1],
  ["Y","H",0.0,2.46649,1,1],
  ["Yb","O",0.0,2.74551,1,1],
  ["Yb","N",0.0,2.84851,1,1],
  ["Yb","S",0.0,3.03649,1,1],
  ["Yb","Se",0.0,3.16649,1,1],
  ["Yb","Te",0.0,3.36649,1,1],
  ["Yb","F",0.0,2.50649,1,1],
  ["Yb","Cl",0.0,2.98249,1,1],
  ["Yb","Br",0.0,3.12945,1,1],
  ["Yb","I",0.0,3.37945,1,1],
  ["Yb","P",0.0,3.13649,1,1],
  ["Yb","As",0.0,3.19649,1,1],
  ["Yb","H",0.0,2.42649,1,1],
  ["Zn","O",0.0,2.41693,1,1],
  ["Zn","S",0.0,2.80293,1,1],
  ["Zn","Se",0.0,2.93293,1,1],
  ["Zn","Te",0.0,3.16293,1,1],
  ["Zn","F",0.0,2.38293,1,1],
  ["Zn","Cl",0.0,2.72293,1,1],
  ["Zn","Br",0.0,2.86293,1,1],
  ["Zn","I",0.0,3.07293,1,1],
  ["Zn","N",0.0,2.43293,1,1],
  ["Zn","P",0.0,2.86293,1,1],
  ["Zn","As",0.0,2.95293,1,1],
  ["Zn","H",0.0,2.13293,1,1],
  ["Zr","O",0.0,3.09651,1,1],
  ["Zr","F",0.0,2.99651,1,1],
  ["Zr","Cl",0.0,3.33651,1,1],
  ["Zr","S",0.0,2.91004,1,1],
  ["Zr","Se",0.0,3.03004,1,1],
  ["Zr","Te",0.0,3.17004,1,1],
  ["Zr","Br",0.0,2.98004,1,1],
  ["Zr","I",0.0,3.19004,1,1],
  ["Zr","N",0.0,2.65004,1,1],
  ["Zr","P",0.0,3.02004,1,1],
  ["Zr","As",0.0,3.07004,1,1],
  ["Zr","H",0.0,2.29004,1,1],
];
// Build lookup: Map<a1Element, [{a2, minA, maxA, showPolyhedra, boundaryMode}]>
const VESTA_SBOND_MAP = new Map();
for (const [a1, a2, minA, maxA, showPoly, bMode] of VESTA_SBOND) {
  if (!VESTA_SBOND_MAP.has(a1)) VESTA_SBOND_MAP.set(a1, []);
  VESTA_SBOND_MAP.get(a1).push({ a2, minA, maxA, showPolyhedra: showPoly === 1, boundaryMode: bMode });
}
const VESTA_MAX_CUTOFF = 4.4631; // largest maxA in the table
const PARENT_CELL_COLOR = 0xe53935;
const CHILD_CELL_COLOR = 0x8c8c8c;
const MAX_BOND_RADIUS = 0.14;
// VESTA reference geometry is expressed in scene-space angstroms. Scaling it
// by a lattice axis makes identical moments grow with supercell dimensions.
const REFERENCE_MOMENT = Math.sqrt(3 * 0.569 ** 2);
const REFERENCE_ARROW_LENGTH = 2.72;
const REFERENCE_SHAFT_RADIUS = 0.12;
const REFERENCE_HEAD_RADIUS = 0.33;
const REFERENCE_HEAD_LENGTH = 0.52;
const MIN_MAGNETIC_ARROW_LENGTH = 0.03;
const MIN_MAGNETIC_SHAFT_LENGTH = 0.008;
const CUBE_PX = 116;
const CUBE_MARGIN = 6;
const CUBE_VIEW = 1.2;
const DEPTH_CUE_NEAR_RADIUS = 3.6;
const DEPTH_CUE_FAR_RADIUS = 9.0;
const Y_AXIS = new THREE.Vector3(0, 1, 0);
const UP_Z = new THREE.Vector3(0, 0, 1);
const UP_Y = new THREE.Vector3(0, 1, 0);
const SHARED_RESOURCE = "modeViewerShared";
const atomGeometryCache = new Map();
const atomMaterialCache = new Map();
const bondMaterialCache = new Map();
const unitCylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, 20);
unitCylinderGeometry.userData[SHARED_RESOURCE] = true;
const unitConeGeometry = new THREE.ConeGeometry(1, 1, 20);
unitConeGeometry.userData[SHARED_RESOURCE] = true;
const magneticArrowMaterial = new THREE.MeshPhongMaterial({
  color: 0xf00000,
  specular: 0xff8080,
  shininess: 38,
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

// Compute safe periodic translation bounds from the reciprocal lattice.
// For a Cartesian cutoff R, any neighbor vector of length <= R satisfies
//   |fractional_i| <= R * ||row_i(A^-1)||.
// Add 1 for the canonical atom-pair fractional difference in (-1,1).
function periodicShifts(vectors, cutoffA) {
  // vectors is [THREE.Vector3 a, b, c] (columns of the Cartesian lattice matrix)
  // Build 3x3 matrix A whose columns are a,b,c, then invert.
  const [av, bv, cv] = vectors;
  // A columns: a,b,c. Rows of A:
  const A = [
    [av.x, bv.x, cv.x],
    [av.y, bv.y, cv.y],
    [av.z, bv.z, cv.z],
  ];
  // 3x3 inverse via cofactors
  const det = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
             -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
             +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
  if (Math.abs(det) < 1e-15) return [-1,0,1].flatMap(i=>[-1,0,1].flatMap(j=>[-1,0,1].map(k=>[i,j,k])));
  const invDet = 1 / det;
  // Rows of G = A^-1
  const G = [
    [(A[1][1]*A[2][2]-A[1][2]*A[2][1])*invDet, -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invDet, (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invDet],
    [-(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invDet, (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invDet, -(A[0][0]*A[1][2]-A[0][2]*A[1][0])*invDet],
    [(A[1][0]*A[2][1]-A[1][1]*A[2][0])*invDet, -(A[0][0]*A[2][1]-A[0][1]*A[2][0])*invDet, (A[0][0]*A[1][1]-A[0][1]*A[1][0])*invDet],
  ];
  const bounds = G.map(row => Math.ceil(cutoffA * Math.hypot(...row) + 1));
  const [nx, ny, nz] = bounds;
  const shifts = [];
  for (let i = -nx; i <= nx; i++)
    for (let j = -ny; j <= ny; j++)
      for (let k = -nz; k <= nz; k++)
        shifts.push([i, j, k]);
  return shifts;
}

function isInsideUnit(value) {
  return value.every(component => component >= -1e-8 && component <= 1 + 1e-8);
}

function makeCylinderBetween(left, right, radius, color) {
  const delta = right.clone().sub(left);
  const length = delta.length();
  if (length < EPS) return null;
  const materialKey = softenedMaterialColor(color);
  if (!bondMaterialCache.has(materialKey)) {
    const base = new THREE.Color(materialKey);
    const material = new THREE.MeshPhongMaterial({
      color: base,
      specular: base.clone().lerp(new THREE.Color(0xffffff), 0.25),
      shininess: 28,
    });
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
    const geometry = new THREE.SphereGeometry(radius, 32, 24);
    geometry.userData[SHARED_RESOURCE] = true;
    atomGeometryCache.set(geometryKey, geometry);
  }
  const color = softenedMaterialColor(atomColor(atom));
  if (!atomMaterialCache.has(color)) {
    const base = new THREE.Color(color);
    const material = new THREE.MeshPhongMaterial({
      color: base,
      specular: base.clone().lerp(new THREE.Color(0xffffff), 0.27),
      shininess: 28,
      emissive: base.clone().multiplyScalar(0.02),
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

function makePolyhedronEdgeMesh(edgeGeometry, color) {
  // WebGL ignores linewidth > 1, so draw each edge as a thin cylinder.
  const pos = edgeGeometry.attributes.position;
  const group = new THREE.Group();
  const radius = 0.013;
  const edgeMat = new THREE.MeshPhongMaterial({
    color: 0xd7d9dc,
    emissive: new THREE.Color(0xd7d9dc).multiplyScalar(0.32),
    specular: 0xe4e6e9,
    shininess: 60,
    transparent: true,
    opacity: 0.86,
    depthWrite: false,
  });
  for (let i = 0; i < pos.count; i += 2) {
    const a = new THREE.Vector3().fromBufferAttribute(pos, i);
    const b = new THREE.Vector3().fromBufferAttribute(pos, i + 1);
    const delta = b.clone().sub(a);
    const length = delta.length();
    if (length < 1e-8) continue;
    const mesh = new THREE.Mesh(unitCylinderGeometry, edgeMat);
    mesh.scale.set(radius, length, radius);
    mesh.position.copy(a).add(b).multiplyScalar(0.5);
    mesh.quaternion.setFromUnitVectors(Y_AXIS, delta.normalize());
    group.add(mesh);
  }
  edgeMat.userData[SHARED_RESOURCE] = true;
  return group;
}

function makePolyhedronMesh(centerAtom, points, opacity) {
  const unique = [];
  const seen = new Set();
  for (const point of points) {
    const key = point.toArray().map(value => fixed(value, 5)).join(",");
    if (seen.has(key)) continue;
    seen.add(key);
    unique.push(point.clone());
  }
  // No chemistry-based vertex cap: ConvexGeometry requires >= 4 non-coplanar points.
  // A hard upper renderer safety limit (200) protects only against pathological cases.
  if (unique.length < 4 || unique.length > 200) return null;
  let geometry;
  try {
    geometry = new ConvexGeometry(unique);
  } catch {
    return null;
  }
  const color = softenedMaterialColor(atomColor(centerAtom));
  const base = new THREE.Color(color);
  const group = new THREE.Group();
  group.add(new THREE.Mesh(
    geometry,
    new THREE.MeshPhongMaterial({
      color: base,
      emissive: new THREE.Color(0x000000),   // no self-emission: face angle drives brightness
      specular: base.clone().lerp(new THREE.Color(0xffffff), 0.22),
      shininess: 30,
      flatShading: true,
      transparent: true,
      opacity,
      side: THREE.DoubleSide,
      depthWrite: false,
      polygonOffset: true,
      polygonOffsetFactor: 1,
      polygonOffsetUnits: 1,
    }),
  ));
  group.add(makePolyhedronEdgeMesh(new THREE.EdgesGeometry(geometry), color));
  return group;
}

// Build bond rules from the VESTA SBOND table.
// Only rules whose A1 and A2 are both present in the structure are active.
// A1 is always the center (polyhedron center); A2 supplies vertices.
function vestaStyleBondRules(atoms, vectors) {
  const presentElements = new Set(atoms.map(a => a.element));
  const rules = [];
  for (const [a1Elem, specs] of VESTA_SBOND_MAP) {
    if (!presentElements.has(a1Elem)) continue;
    for (const spec of specs) {
      if (!presentElements.has(spec.a2)) continue;
      rules.push({
        center: a1Elem,
        ligand: spec.a2,
        min: spec.minA,
        max: spec.maxA,
        polyhedron: spec.showPolyhedra,
        boundaryMode: spec.boundaryMode,
      });
    }
  }
  return rules;
}

function buildBondTopology(atoms, vectors) {
  const connections = [];
  const shifts = periodicShifts(vectors, VESTA_MAX_CUTOFF);
  for (const rule of vestaStyleBondRules(atoms, vectors)) {
    const leftAtoms = atoms.filter(atom => atom.element === rule.center);
    const rightAtoms = atoms.filter(atom => atom.element === rule.ligand);
    const seen = new Set();
    for (const left of leftAtoms) {
      const leftPosition = combine(vectors, left.frac);
      for (const right of rightAtoms) {
        for (const shift of shifts) {
          if (left.sourceKey === right.sourceKey && shift.every(value => value === 0)) continue;
          const shifted = right.frac.map((value, axis) => value + shift[axis]);
          const distance = leftPosition.distanceTo(combine(vectors, shifted));
          if (distance < rule.min - 1e-8 || distance > rule.max + 1e-8) continue;
          const forward = `${left.sourceKey}>${right.sourceKey}@${shift.join(",")}`;
          const reverse = `${right.sourceKey}>${left.sourceKey}@${shift.map(value => -value).join(",")}`;
          const key = forward < reverse ? forward : reverse;
          if (seen.has(key)) continue;
          seen.add(key);
          connections.push({
            leftKey: left.sourceKey,
            rightKey: right.sourceKey,
            shift: [...shift],
            polyhedron: rule.polyhedron,
          });
        }
      }
    }
  }
  return connections;
}

function buildBondGeometry(
  atoms,
  displayedAtoms,
  vectors,
  topology,
  { showBonds, showPolyhedra, polyhedronOpacity, visualScale },
) {
  const group = new THREE.Group();
  if (!showBonds && !showPolyhedra) return group;
  const referenceLength = Math.max(number(visualScale, 1), 1);
  const radius = Math.min(referenceLength * 0.012, MAX_BOND_RADIUS);
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
      if (connection.polyhedron) {
        const entryKey = left.visualKey || `${left.sourceKey}@${displayShift.join(",")}`;
        const entry = polyhedra.get(entryKey) || { atom: sourceLeft, points: [] };
        entry.points.push(rightPosition);
        polyhedra.set(entryKey, entry);
      }
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
      const polyhedron = makePolyhedronMesh(entry.atom, entry.points, polyhedronOpacity);
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

function magneticArrow(origin, direction) {
  const magnitude = direction.length();
  if (magnitude < EPS) return null;
  const amplitudeRatio = magnitude / REFERENCE_MOMENT;
  const totalLength = REFERENCE_ARROW_LENGTH * amplitudeRatio;
  if (totalLength < MIN_MAGNETIC_ARROW_LENGTH) return null;
  const widthScale = Math.min(1, Math.sqrt(amplitudeRatio));
  const headLength = Math.min(REFERENCE_HEAD_LENGTH * widthScale, totalLength * 0.35);
  const shaftLength = Math.max(totalLength - headLength, MIN_MAGNETIC_SHAFT_LENGTH);
  const shaftRadius = REFERENCE_SHAFT_RADIUS * widthScale;
  const headRadius = REFERENCE_HEAD_RADIUS * widthScale;
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
  constructor({
    canvas,
    empty,
    controls,
    status,
    playButton,
    resetButton,
    scaleInput,
    scaleOutput,
    arrowsInput,
    momentsInput,
    arrowOnlyInput,
    bondsInput,
    polyhedraInput,
    polyhedronOpacityInput,
    polyhedronOpacityOutput,
  }) {
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
    this.arrowOnlyInput = arrowOnlyInput;
    this.bondsInput = bondsInput;
    this.polyhedraInput = polyhedraInput;
    this.polyhedronOpacityInput = polyhedronOpacityInput;
    this.polyhedronOpacityOutput = polyhedronOpacityOutput;
    this.payload = null;
    this.modes = [];
    this.amplitudes = new Map();
    this.hiddenChildSites = new Set();
    this.structure = null;
    this.playing = false;
    this.masterAmplitude = 1;
    this.animationStart = 0;
    this.animationPhaseOffset = 0;
    this.lastAnimationUpdate = 0;
    this.needsRender = true;
    this.fitted = false;
    this.bondTopology = null;
    this.liveLayoutResize = false;
    this.keyLightRight = new THREE.Vector3();
    this.keyLightUp = new THREE.Vector3();

    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setClearColor(0xffffff, 1);
    this.renderer.autoClear = false;
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.Fog(0xffffff, 1e6, 1e6 + 1);
    this.lighting = addLights(this.scene);
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

    loadElementStyles("./vesta_elements.csv").then(() => this.updateScene());
    this.bindControls();
    this.bindViewCube();
    this.resizeObserver = new ResizeObserver(() => {
      this.resize();
      const rect = this.canvas.getBoundingClientRect();
      if (this.structure && !this.liveLayoutResize && rect.width > 10 && rect.height > 10) {
        this.fitView(true);
      }
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
    [
      this.scaleInput,
      this.arrowsInput,
      this.momentsInput,
      this.arrowOnlyInput,
      this.bondsInput,
      this.polyhedraInput,
      this.polyhedronOpacityInput,
    ].forEach(input => {
      input?.addEventListener("input", () => {
        if (input === this.scaleInput && this.scaleOutput) {
          this.scaleOutput.textContent = number(this.scaleInput.value, 1).toFixed(2);
        }
        if (input === this.polyhedronOpacityInput && this.polyhedronOpacityOutput) {
          this.polyhedronOpacityOutput.textContent = `${Math.round(number(input.value, 70))}%`;
        }
        this.syncDisplayControls();
        this.updateScene();
      });
      input?.addEventListener("change", () => {
        this.syncDisplayControls();
        this.updateScene();
      });
    });
    this.syncDisplayControls();
  }

  syncDisplayControls() {
    const arrowsOnly = Boolean(this.arrowOnlyInput?.checked);
    if (this.bondsInput) this.bondsInput.disabled = arrowsOnly;
    if (this.polyhedraInput) this.polyhedraInput.disabled = arrowsOnly;
    if (this.polyhedronOpacityInput) {
      this.polyhedronOpacityInput.disabled = arrowsOnly || !this.polyhedraInput?.checked;
    }
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
    this.hiddenChildSites.clear();
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
    const childSites = this.undistortedChildSites();
    const childSitesByElement = new Map();
    for (const site of childSites) {
      if (!childSitesByElement.has(site.element)) childSitesByElement.set(site.element, []);
      childSitesByElement.get(site.element).push(site);
    }
    const atomVisibility = childSites.length ? `<fieldset class="mode-atom-visibility">
      <legend>Undistorted structure</legend>
      <div class="mode-element-list">${[...childSitesByElement.entries()].map(([element, sites]) => `<label class="mode-atom-toggle mode-element-toggle">
        <input type="checkbox" data-mode-element="${this.escape(element)}" checked>
        <span>${this.escape(`${element} (${sites.length})`)}</span>
      </label>`).join("")}</div>
      <div class="mode-atom-list">${childSites.map(site => `<label class="mode-atom-toggle">
        <input type="checkbox" data-mode-child-site="${this.escape(site.id)}" data-mode-site-element="${this.escape(site.element)}" checked>
        <span>${this.escape(`${site.id} (${site.wyckoff})`)}</span>
      </label>`).join("")}</div>
    </fieldset>` : "";
    this.controlsNode.innerHTML = `${master}${atomVisibility}${groups || '<p class="mode-viewer-note">No visualizable mode definitions.</p>'}`;
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
    const siteVisibilityInputs = [...this.controlsNode.querySelectorAll("input[data-mode-child-site]")];
    const elementVisibilityInputs = [...this.controlsNode.querySelectorAll("input[data-mode-element]")];
    const syncElementVisibility = element => {
      const elementInput = elementVisibilityInputs.find(input => input.dataset.modeElement === element);
      const siteInputs = siteVisibilityInputs.filter(input => input.dataset.modeSiteElement === element);
      if (!elementInput || !siteInputs.length) return;
      const checkedCount = siteInputs.filter(input => input.checked).length;
      elementInput.checked = checkedCount === siteInputs.length;
      elementInput.indeterminate = checkedCount > 0 && checkedCount < siteInputs.length;
    };
    siteVisibilityInputs.forEach(input => {
      input.addEventListener("change", () => {
        const childSite = input.dataset.modeChildSite;
        if (!childSite) return;
        if (input.checked) this.hiddenChildSites.delete(childSite);
        else this.hiddenChildSites.add(childSite);
        syncElementVisibility(input.dataset.modeSiteElement);
        this.updateScene();
      });
    });
    elementVisibilityInputs.forEach(input => {
      input.addEventListener("change", () => {
        const element = input.dataset.modeElement;
        if (!element) return;
        input.indeterminate = false;
        siteVisibilityInputs.filter(siteInput => siteInput.dataset.modeSiteElement === element).forEach(siteInput => {
          siteInput.checked = input.checked;
          const childSite = siteInput.dataset.modeChildSite;
          if (input.checked) this.hiddenChildSites.delete(childSite);
          else this.hiddenChildSites.add(childSite);
        });
        this.updateScene();
      });
    });
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

  undistortedChildSites() {
    const atoms = Array.isArray(this.payload?.undistorted_atoms) ? this.payload.undistorted_atoms : [];
    const sites = [];
    const seenSites = new Set();
    const seenAtomIds = new Set();
    for (const atom of atoms) {
      const id = atom?.child_site;
      const wyckoff = atom?.site;
      const atomIds = atom?.atom_ids;
      const multiplicity = atom?.multiplicity;
      if (
        typeof id !== "string"
        || !id
        || typeof wyckoff !== "string"
        || !Array.isArray(atomIds)
        || atomIds.length !== multiplicity
        || seenSites.has(id)
        || atomIds.some(atomId => (
          typeof atomId !== "string" || !atomId || seenAtomIds.has(atomId)
        ))
      ) return [];
      seenSites.add(id);
      atomIds.forEach(atomId => seenAtomIds.add(atomId));
      sites.push({ id, wyckoff, element: elementFromLabel(id), atomIds: [...atomIds] });
    }
    return sites;
  }

  escape(value) {
    return String(value ?? "").replace(/[&<>"']/g, char => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    }[char]));
  }

  buildAtoms() {
    const sourceAtoms = [];
    const childSiteByAtomId = new Map();
    for (const site of this.undistortedChildSites()) {
      site.atomIds.forEach(atomId => childSiteByAtomId.set(atomId, site.id));
    }
    const viewerAtoms = Array.isArray(this.payload?.viewer_atoms) ? this.payload.viewer_atoms : [];
    const seenAtomIds = new Set();
    let invalid = viewerAtoms.length !== childSiteByAtomId.size;
    viewerAtoms.forEach(atom => {
      const xyz = atom?.xyz;
      if (
        !Array.isArray(xyz)
        || xyz.length !== 3
        || typeof atom.atom_id !== "string"
        || !atom.atom_id
        || typeof atom.child_site !== "string"
        || !atom.child_site
        || typeof atom.element !== "string"
        || !atom.element
        || !Number.isInteger(atom.site_order)
        || !Number.isInteger(atom.atom_index)
        || childSiteByAtomId.get(atom.atom_id) !== atom.child_site
        || seenAtomIds.has(atom.atom_id)
      ) {
        invalid = true;
        return;
      }
      seenAtomIds.add(atom.atom_id);
      const label = atom.child_site;
      sourceAtoms.push({
        sourceKey: atom.atom_id,
        siteOrder: atom.site_order,
        atomIndex: atom.atom_index,
        childSite: atom.child_site,
        label,
        element: elementFromLabel(atom.element),
        frac: vector3(xyz),
        sourceFrac: vector3(xyz),
      });
    });
    if (invalid || seenAtomIds.size !== childSiteByAtomId.size) return [];

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
      modeRows(mode.definition).forEach(row => {
        if (typeof row?.atom_id !== "string" || !row.atom_id) return;
        const key = row.atom_id;
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
    const atoms = this.buildAtoms().filter(
      atom => !this.hiddenChildSites.has(atom.childSite),
    );
    const { displacements, moments } = this.modeVectors(atoms, scale);
    const arrowsOnly = Boolean(this.arrowOnlyInput?.checked);
    const polyhedronOpacity = Math.min(
      1,
      Math.max(0, number(this.polyhedronOpacityInput?.value, 70) / 100),
    );
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
      group.add(makeCellEdges(distortedParentVectors, PARENT_CELL_COLOR, parentOffset, 1, false, 2));
    }
    group.add(makeCellEdges(vectors, CHILD_CELL_COLOR, new THREE.Vector3(), 0.62));
    if (!arrowsOnly) {
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
        const referenceAtoms = [...sourceAtoms.values()].map(atom => ({
          ...atom,
          frac: [...atom.sourceFrac],
        }));
        this.bondTopology = buildBondTopology(referenceAtoms, baseVectors);
      }
      const displayedAtoms = atoms.map(atom => {
        const displacement = displacements.get(atom.sourceKey) || [0, 0, 0];
        return {
          ...atom,
          frac: atom.frac.map((value, axis) => value + displacement[axis]),
        };
      });
      group.add(buildBondGeometry(
        [...sourceAtoms.values()],
        displayedAtoms,
        vectors,
        this.bondTopology,
        {
          showBonds: Boolean(this.bondsInput?.checked),
          showPolyhedra: Boolean(this.polyhedraInput?.checked),
          polyhedronOpacity,
          visualScale,
        },
      ));
    }
    for (const atom of atoms) {
      const basePosition = combine(vectors, atom.frac);
      const displacement = combine(vectors, displacements.get(atom.sourceKey) || [0, 0, 0]);
      const position = basePosition.clone().add(displacement);
      const fallbackRadius = visualScale * 0.045;
      if (!arrowsOnly) {
        const mesh = makeModeAtomMesh(atom, fallbackRadius);
        mesh.position.copy(position);
        mesh.userData = { label: atom.label, fractional: atom.frac };
        group.add(mesh);
      }
      if (arrowsOnly || this.arrowsInput?.checked) {
        const displacementArrow = arrow(basePosition, displacement, 0x1677b8, visualScale);
        if (displacementArrow) group.add(displacementArrow);
      }
      if (!this.momentsInput || this.momentsInput.checked) {
        const moment = combine(vectors, moments.get(atom.sourceKey) || [0, 0, 0]);
        const momentArrow = magneticArrow(position, moment);
        if (momentArrow) group.add(momentArrow);
      }
    }
    this.structure = group;
    this.scene.add(group);
    this.needsRender = true;
    this.empty.hidden = atoms.length > 0;
    this.empty.textContent = atoms.length ? "" : "All atoms are hidden";
    if (atoms.length && this.statusNode) {
      this.statusNode.textContent = `${atoms.length} displayed atoms / ${this.modes.length} modes`;
    }
    const canvasRect = this.canvas.getBoundingClientRect();
    if (!this.fitted && canvasRect.width > 10 && canvasRect.height > 10) {
      this.fitView();
      this.fitted = true;
    }
  }

  setLiveLayoutResize(active) {
    const next = Boolean(active);
    if (next === this.liveLayoutResize) return;
    this.liveLayoutResize = next;
    if (!next) this.refreshLayout(true);
  }

  refreshLayout(preserveOrientation = false) {
    this.resize();
    if (!this.structure) return;
    this.fitView(preserveOrientation);
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
    this.scene.fog.near = radius * DEPTH_CUE_NEAR_RADIUS;
    this.scene.fog.far = radius * DEPTH_CUE_FAR_RADIUS;
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
    const lightDistance = this.camera.position.distanceTo(this.trackball.target);
    this.keyLightRight.set(1, 0, 0).applyQuaternion(this.camera.quaternion);
    this.keyLightUp.set(0, 1, 0).applyQuaternion(this.camera.quaternion);
    this.lighting.keyLight.position.copy(this.camera.position)
      .addScaledVector(this.keyLightRight, -lightDistance * 0.14)
      .addScaledVector(this.keyLightUp, lightDistance * 0.1);
    this.lighting.keyTarget.position.copy(this.trackball.target);
    this.lighting.keyTarget.updateMatrixWorld();
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
