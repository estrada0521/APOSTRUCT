import { OutputTabs } from "./output-tabs.js?v=1";
import { favoriteStorageKey, isFavorite, toggleFavorite } from "./favorites.js?v=1";
import { ModeStructureViewer } from "./mode3d.js?v=81";

const $ = id => document.getElementById(id);
function syncFavorite() {
  const entry = tabs.active?.snapshot.savedInput;
  const saved = isFavorite(entry);
  const button = $("output-favorite");
  button.disabled = !entry;
  button.setAttribute("aria-pressed", String(saved));
  button.setAttribute("aria-label", saved ? "Remove from favorites" : "Add to favorites");
}
$("output-favorite").addEventListener("click", () => {
  try {
    toggleFavorite(tabs.active?.snapshot.savedInput);
    syncFavorite();
    outputError(null);
  } catch {
    outputError("Favorites could not be saved in this browser.");
  }
});
window.addEventListener("storage", event => {
  if (event.key === favoriteStorageKey || event.key === null) syncFavorite();
});
const modeViewer = new ModeStructureViewer({
  canvas: $("mode-viewer-canvas"),
  empty: $("mode-viewer-empty"),
  controls: $("mode-viewer-controls"),
  homeButton: $("mode-viewer-home"),
  bondsInput: $("mode-viewer-bonds"),
  polyhedraInput: $("mode-viewer-polyhedra"),
});

const modeWorkspace = document.querySelector(".mode-viewer-workspace");
const modeSidebar = document.querySelector(".mode-viewer-sidebar");
const modeResizer = document.querySelector(".mode-viewer-resizer");
function sliderLabelWidth() {
  const rows = [...document.querySelectorAll("#mode-viewer-controls .mode-slider")];
  if (!rows.length) return 0;
  const row = getComputedStyle(rows[0]);
  const card = rows[0].closest(".mode-site-card");
  const padding = card ? parseFloat(getComputedStyle(card).paddingLeft) + parseFloat(getComputedStyle(card).paddingRight) : 36;
  const label = Math.max(...rows.map(item => item.querySelector(".mode-slider-label")?.scrollWidth || 0));
  return padding + 80 + 44 + 2 * parseFloat(row.columnGap || "6") + label + 24;
}

function setModeSidebarWidth(width) {
  if (window.matchMedia && window.matchMedia("(max-width: 600px)").matches) return;
  if (!modeWorkspace || !modeSidebar || !modeResizer) return;
  const bounds = modeWorkspace.getBoundingClientRect();
  const handleWidth = modeResizer.getBoundingClientRect().width;
  const availableWidth = bounds.width - handleWidth;
  if (availableWidth <= 0) return;
  const minimumSidebar = Math.max(Math.min(220, Math.max(100, availableWidth * 0.3)), sliderLabelWidth());
  const minimumViewer = Math.min(260, Math.max(120, availableWidth * 0.3));
  const maximumSidebar = Math.max(minimumSidebar, availableWidth - minimumViewer);
  const nextWidth = Math.min(maximumSidebar, Math.max(minimumSidebar, width));
  modeWorkspace.style.setProperty("--mode-sidebar-width", `${nextWidth}px`);
  modeResizer.setAttribute("aria-valuemin", String(Math.round(minimumSidebar)));
  modeResizer.setAttribute("aria-valuemax", String(Math.round(maximumSidebar)));
  modeResizer.setAttribute("aria-valuenow", String(Math.round(nextWidth)));
}

modeResizer?.addEventListener("pointerdown", event => {
  event.preventDefault();
  modeResizer.focus();
  modeViewer.setLiveLayoutResize(true);
  modeResizer.setPointerCapture(event.pointerId);
  modeResizer.classList.add("is-dragging");
});
modeResizer?.addEventListener("pointermove", event => {
  if (!modeResizer.hasPointerCapture(event.pointerId) || !modeWorkspace) return;
  setModeSidebarWidth(event.clientX - modeWorkspace.getBoundingClientRect().left);
});
const stopModeResize = event => {
  if (modeResizer?.hasPointerCapture(event.pointerId)) {
    modeResizer.releasePointerCapture(event.pointerId);
  }
  modeResizer?.classList.remove("is-dragging");
  modeViewer.setLiveLayoutResize(false);
};
modeResizer?.addEventListener("pointerup", stopModeResize);
modeResizer?.addEventListener("pointercancel", stopModeResize);
modeResizer?.addEventListener("lostpointercapture", () => {
  modeResizer.classList.remove("is-dragging");
  modeViewer.setLiveLayoutResize(false);
});
modeResizer?.addEventListener("keydown", event => {
  if (!modeSidebar || !["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) return;
  event.preventDefault();
  const current = modeSidebar.getBoundingClientRect().width;
  const target = event.key === "Home"
    ? 0
    : event.key === "End"
      ? Number.POSITIVE_INFINITY
      : current + (event.key === "ArrowLeft" ? -20 : 20);
  setModeSidebarWidth(target);
});
window.addEventListener("resize", () => {
  if (modeSidebar) setModeSidebarWidth(modeSidebar.getBoundingClientRect().width);
});
new MutationObserver(() => {
  if (modeSidebar && !modeSidebar.closest("[hidden]")) setModeSidebarWidth(modeSidebar.getBoundingClientRect().width);
}).observe($("mode-viewer-controls"), { childList: true });

function showViewer() {
  modeViewer.setVisible(!document.hidden);
  requestAnimationFrame(() => {
    modeViewer.resize();
    if (!modeViewer.fitted) { modeViewer.fitView(); modeViewer.fitted = true; }
  });
}

const channel = location.hash.slice(1);
const source = window.opener;
const tabs = new OutputTabs();
let displayedRevision = null;

function send(type, detail = {}) {
  if (source && !source.closed) source.postMessage({ app: "apostruct-output", channel, type, ...detail }, location.origin);
}
function outputError(error) {
  $("output-error").hidden = !error;
  $("output-error").textContent = error?.message || String(error || "");
}
function saveActive() {
  const entry = tabs.entries.get(displayedRevision);
  if (!entry) return;
  const scroll = $("mode-viewer-controls").querySelector(".mode-controls-scroll");
  entry.ui = {
    viewer: modeViewer.captureState(),
    controlsScroll: scroll?.scrollTop || 0,
    pageScroll: [window.scrollX, window.scrollY],
    sidebarWidth: modeWorkspace.style.getPropertyValue("--mode-sidebar-width"),
  };
  entry.amplitudes = entry.ui.viewer.amplitudes;
  entry.masterAmplitude = entry.ui.viewer.masterAmplitude;
  entry.view = "viewer";
  send("result-state", { revision: entry.revision, state: {
    ui: entry.ui, amplitudes: entry.amplitudes, masterAmplitude: entry.masterAmplitude, view: entry.view,
  } });
}
let tabStripVisible = false;
const tabReflowAnimations = new WeakMap();
function renderTabs() {
  const bar = $("output-tabs");
  const closing = [...bar.children].some(item => !tabs.entries.has(Number(item.dataset.revision)));
  const previousPositions = new Map();
  if (closing && !window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
    for (const item of bar.children) {
      previousPositions.set(item, item.getBoundingClientRect().left);
      tabReflowAnimations.get(item)?.cancel();
    }
  }
  if (tabs.entries.size > 1) tabStripVisible = true;
  else if (tabs.entries.size === 0) tabStripVisible = false;
  bar.hidden = !tabStripVisible;
  document.querySelector(".output-page").classList.toggle("single-result", !tabStripVisible);
  for (const item of [...bar.children]) {
    if (!tabs.entries.has(Number(item.dataset.revision))) item.remove();
  }
  for (const entry of tabs.entries.values()) {
    let item = [...bar.children].find(node => Number(node.dataset.revision) === entry.revision);
    if (!item) {
      item = document.createElement("div");
      item.className = "output-tab";
      item.dataset.revision = entry.revision;
      const tab = document.createElement("button");
      tab.type = "button"; tab.id = "result-tab-" + entry.revision;
      tab.setAttribute("role", "tab");
      tab.setAttribute("aria-controls", "mode-details-section");
      tab.addEventListener("click", () => activate(entry.revision));
      const close = document.createElement("button");
      close.type = "button"; close.className = "output-tab-close";
      const icon = document.createElementNS("http://www.w3.org/2000/svg", "svg");
      icon.setAttribute("viewBox", "0 0 16 16"); icon.setAttribute("aria-hidden", "true");
      const cross = document.createElementNS("http://www.w3.org/2000/svg", "path");
      cross.setAttribute("d", "M4 4 L12 12 M12 4 L4 12");
      cross.setAttribute("fill", "none"); cross.setAttribute("stroke", "currentColor");
      cross.setAttribute("stroke-width", "1.5"); cross.setAttribute("stroke-linecap", "round");
      icon.appendChild(cross); close.appendChild(icon);
      close.addEventListener("click", () => {
        saveActive();
        tabs.close(entry.revision);
        send("close", { revision: entry.revision });
        if (displayedRevision !== tabs.activeRevision) displayActive();
        else renderTabs();
        $("result-tab-" + tabs.activeRevision)?.focus({ preventScroll: true });
      });
      item.append(tab, close); bar.append(item);
    }
    const [tab, close] = item.children;
    const active = entry.revision === tabs.activeRevision;
    item.classList.toggle("is-active", active);
    tab.setAttribute("aria-selected", String(active));
    tab.tabIndex = active ? 0 : -1;
    const title = entry.snapshot.title || "Result";
    if (tab.textContent !== title) tab.textContent = title;
    tab.title = title;
    close.setAttribute("aria-label", "Close " + title);
  }
  for (const item of bar.children) {
    const previousLeft = previousPositions.get(item);
    if (previousLeft === undefined) continue;
    const offset = previousLeft - item.getBoundingClientRect().left;
    if (Math.abs(offset) < 0.5) continue;
    const animation = item.animate([
      { transform: `translateX(${offset}px)` },
      { transform: "translateX(0)" },
    ], { duration: 180, easing: "cubic-bezier(0.2, 0, 0, 1)" });
    tabReflowAnimations.set(item, animation);
  }
}
function activate(revision) {
  if (revision === tabs.activeRevision) return;
  saveActive();
  if (tabs.select(revision)) displayActive();
}
function displayActive() {
  renderTabs();
  const entry = tabs.active;
  syncFavorite();
  $("mode-details-section").hidden = !entry;
  $("output-empty").hidden = Boolean(entry);
  displayedRevision = entry?.revision ?? null;
  send("select", { revision: displayedRevision });
  if (!entry) {
    modeViewer.setVisible(false);
    modeViewer.setData(null);
    document.title = "APOSTRUCT — Output";
    return;
  }
  outputError(null);
  modeWorkspace.style.removeProperty("--mode-sidebar-width");
  modeViewer.setVisible(false);
  modeViewer.setData(entry.snapshot.viewer);
  if (!modeViewer.bondOptionsSuppressed) {
    $("mode-viewer-bonds").checked = true;
    $("mode-viewer-polyhedra").checked = true;
  }
  if (entry.ui) modeViewer.restoreState(entry.ui.viewer);
  else {
    if (entry.amplitudes?.length) modeViewer.setAmplitudes(entry.amplitudes);
    modeViewer.masterAmplitude = entry.masterAmplitude ?? 1;
    modeViewer.updateControlValues();
    modeViewer.updateScene();
  }
  if (entry.ui?.sidebarWidth) modeWorkspace.style.setProperty("--mode-sidebar-width", entry.ui.sidebarWidth);
  let heading = entry.snapshot.heading || [entry.snapshot.title || "APOSTRUCT — Output"];
  if (entry.snapshot.material && entry.snapshot.headingLayout !== 2) {
    const groups = (entry.snapshot.groupTransition || "").split(" → ");
    heading = [entry.snapshot.material, [groups[0], (heading[1] || "").replace(" | ", "・"), groups[1]].filter(Boolean).join(" → ")];
  }
  $("result-title").replaceChildren(...heading.map((text, index) => {
    const line = document.createElement("span");
    if (index) {
      line.textContent = text;
      line.className = "result-heading-detail";
    } else {
      const material = document.createElement("span");
      material.className = "result-material";
      material.textContent = text;
      line.append(material);
    }
    return line;
  }));
  document.title = (entry.snapshot.title || "APOSTRUCT") + " — Output";
  $("mode-details-section").setAttribute("aria-labelledby", "result-tab-" + entry.revision);
  $("output-status").hidden = true;
  showViewer();
  requestAnimationFrame(() => {
    if (displayedRevision !== entry.revision) return;
    if (entry.ui) {
      // Layout may have changed while the tab was hidden; preserve its camera and scale.
      modeViewer.restoreState(entry.ui.viewer);
      const scroll = $("mode-viewer-controls").querySelector(".mode-controls-scroll");
      if (scroll) scroll.scrollTop = entry.ui.controlsScroll;
      window.scrollTo(...entry.ui.pageScroll);
    }
    saveActive();
  });
}
$("output-tabs").addEventListener("keydown", event => {
  if (!event.target.matches('[role="tab"]')) return;
  const keys = [...tabs.entries.keys()], i = keys.indexOf(tabs.activeRevision);
  const index = event.key === "ArrowRight" ? (i + 1) % keys.length
    : event.key === "ArrowLeft" ? (i + keys.length - 1) % keys.length
    : event.key === "Home" ? 0 : event.key === "End" ? keys.length - 1 : -1;
  if (index < 0) return;
  event.preventDefault(); activate(keys[index]);
  $("result-tab-" + keys[index])?.focus();
});
window.addEventListener("message", event => {
  if (event.origin !== location.origin || event.source !== source) return;
  const message = event.data;
  if (!message || message.app !== "apostruct-output" || message.channel !== channel) return;
  try {
    if (message.type === "results" || message.type === "result") {
      const empty = !tabs.active;
      const fresh = message.type === "result" && message.entry?.snapshot && !tabs.entries.has(message.entry.revision);
      if (fresh) saveActive();
      for (const entry of message.entries || [message.entry]) if (entry?.snapshot) tabs.add(entry);
      if (fresh) tabs.select(message.entry.revision);
      else if (empty && message.activeRevision !== undefined) tabs.select(message.activeRevision);
      if (displayedRevision !== tabs.activeRevision) displayActive();
      else renderTabs();
      if (fresh || (empty && tabs.activeRevision === [...tabs.entries.keys()].at(-1))) {
        const bar = $("output-tabs");
        bar.scrollLeft = bar.scrollWidth;
      }
    } else if (message.type === "amplitudes") {
      const entry = tabs.entries.get(message.revision);
      if (!entry) return;
      entry.amplitudes = message.amplitudes;
      if (entry.ui) {
        const values = new Map(entry.ui.viewer.amplitudes.map(row => [row.definition_id, row.value]));
        for (const row of message.amplitudes) values.set(row.definition_id, row.value);
        entry.ui.viewer.amplitudes = [...values].map(([definition_id, value]) => ({ definition_id, value }));
      }
      if (displayedRevision === message.revision) {
        modeViewer.setAmplitudes(message.amplitudes);
        saveActive();
      }
    } else if (message.type === "input-unavailable") disconnected();
  } catch (error) { outputError(error); }
});
$("mode-viewer-controls").addEventListener("input", saveActive);
$("mode-viewer-controls").addEventListener("click", () => queueMicrotask(saveActive));
$("mode-viewer-canvas").addEventListener("pointerup", () => requestAnimationFrame(saveActive));
let wheelSave;
$("mode-viewer-canvas").addEventListener("wheel", () => {
  clearTimeout(wheelSave); wheelSave = setTimeout(saveActive, 150);
}, { passive: true });
for (const id of ["mode-viewer-bonds", "mode-viewer-polyhedra"]) {
  $(id).addEventListener("change", saveActive);
}
window.addEventListener("pagehide", saveActive);
document.addEventListener("visibilitychange", () => {
  modeViewer.setVisible(!document.hidden && Boolean(tabs.active));
});
function disconnected() {
  $("output-status").hidden = false;
  $("output-status").textContent = tabs.entries.size
    ? "Input window disconnected — saved results remain available."
    : "Input window disconnected. Open output from the input window.";
}
if (!source || !channel) $("output-empty").textContent = "Open this page using “Open output” in the input window.";
else {
  send("ready");
  setInterval(() => { if (source.closed) disconnected(); }, 2000);
}
modeViewer.setVisible(false);
