import { OutputTabs } from "./output-tabs.js?v=1";
import { ModeStructureViewer } from "./mode3d.js?v=64";

const $ = id => document.getElementById(id);
const state = { modeDetailText: null, modeDetailTextRendered: false };
function renderModeDetailText(output, text) {
  const fragment = document.createDocumentFragment();
  const source = String(text || "");
  const headingPattern = /^(Parent structure|Subgroup details|Undistorted superstructure|Displacive mode definitions|Magnetic mode definitions|Parent-cell strain mode definitions)$/gm;
  let cursor = 0;
  for (const match of source.matchAll(headingPattern)) {
    appendModeText(fragment, source.slice(cursor, match.index));
    const heading = document.createElement("span");
    heading.className = "mode-details-section-heading";
    heading.textContent = match[0];
    fragment.append(heading);
    cursor = match.index + match[0].length;
  }
  appendModeText(fragment, source.slice(cursor));
  output.replaceChildren(fragment);
}

function ensureModeDetailText() {
  if (state.modeDetailTextRendered) return;
  const output = $("mode-details-output");
  if (state.modeDetailText !== null) {
    renderModeDetailText(output, state.modeDetailText);
    state.modeDetailTextRendered = true;
    return;
  }
  renderModeDetailText(output, "");
  state.modeDetailTextRendered = true;
}

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

function setModeDetailsView(view) {
  const showViewer = view === "viewer";
  if (!showViewer) ensureModeDetailText();
  $("mode-details-text-panel").hidden = showViewer;
  $("mode-details-viewer-panel").hidden = !showViewer;
  $("mode-viewer-controls").hidden = !showViewer;
  $("mode-details-view-select").value = showViewer ? "viewer" : "text";
  modeViewer.setVisible(showViewer && !document.hidden);
  if (showViewer) requestAnimationFrame(() => {
    modeViewer.resize();
    if (!modeViewer.fitted) { modeViewer.fitView(); modeViewer.fitted = true; }
  });
}

$("mode-details-view-select").addEventListener("change", event => setModeDetailsView(event.target.value));

const CHOICE_TITLE = "Depends on a choice (irrep basis, domain, subgroup frame or star arm), not on the subgroup itself";
function choiceNode(text) {
  const span = document.createElement("span");
  span.className = "choice";
  span.title = CHOICE_TITLE;
  span.textContent = text;
  return span;
}

// the subgroup line of the mode details text, with its basis, origin and k-active parts marked as choices
function appendModeText(fragment, text) {
  const line = /^(.*?)(, basis=\{[^}\n]*\})(, origin=\([^)\n]*\))(, s=\d+, i=\d+, )(k-active= [^\n]*)$/m;
  const match = text.match(line);
  if (!match) {
    fragment.append(document.createTextNode(text));
    return;
  }
  fragment.append(document.createTextNode(text.slice(0, match.index) + match[1]));
  fragment.append(document.createTextNode(", "), choiceNode(match[2].slice(2)));
  fragment.append(document.createTextNode(", "), choiceNode(match[3].slice(2)));
  fragment.append(document.createTextNode(match[4]), choiceNode(match[5]));
  fragment.append(document.createTextNode(text.slice(match.index + match[0].length)));
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
    textScroll: [$("mode-details-output").scrollLeft, $("mode-details-output").scrollTop],
    pageScroll: [window.scrollX, window.scrollY],
    sidebarWidth: modeWorkspace.style.getPropertyValue("--mode-sidebar-width"),
  };
  entry.amplitudes = entry.ui.viewer.amplitudes;
  entry.masterAmplitude = entry.ui.viewer.masterAmplitude;
  entry.view = $("mode-details-view-select").value;
  send("result-state", { revision: entry.revision, state: {
    ui: entry.ui, amplitudes: entry.amplitudes, masterAmplitude: entry.masterAmplitude, view: entry.view,
  } });
}
function renderTabs() {
  const bar = $("output-tabs");
  bar.replaceChildren();
  for (const entry of tabs.entries.values()) {
    const item = document.createElement("div");
    item.className = "output-tab";
    const tab = document.createElement("button");
    tab.type = "button"; tab.id = "result-tab-" + entry.revision;
    tab.setAttribute("role", "tab");
    tab.setAttribute("aria-controls", "mode-details-section");
    tab.setAttribute("aria-selected", String(entry.revision === tabs.activeRevision));
    tab.tabIndex = entry.revision === tabs.activeRevision ? 0 : -1;
    tab.textContent = entry.snapshot.title || "Result";
    tab.title = tab.textContent;
    tab.addEventListener("click", () => activate(entry.revision));
    const close = document.createElement("button");
    close.type = "button"; close.className = "output-tab-close"; close.textContent = "×";
    close.setAttribute("aria-label", "Close " + tab.textContent);
    close.addEventListener("click", () => {
      saveActive();
      tabs.close(entry.revision);
      send("close", { revision: entry.revision });
      if (displayedRevision !== tabs.activeRevision) displayActive();
      else renderTabs();
    });
    item.append(tab, close); bar.append(item);
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
  $("mode-details-section").hidden = !entry;
  $("output-empty").hidden = Boolean(entry);
  displayedRevision = entry?.revision ?? null;
  send("select", { revision: displayedRevision });
  if (!entry) {
    state.modeDetailText = null; state.modeDetailTextRendered = false;
    $("mode-details-output").replaceChildren();
    modeViewer.setVisible(false);
    modeViewer.setData(null);
    document.title = "APOSTRUCT — Output";
    return;
  }
  outputError(null);
  state.modeDetailText = entry.snapshot.text || "";
  state.modeDetailTextRendered = false;
  $("mode-details-output").replaceChildren();
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
  $("result-title").textContent = entry.snapshot.title || "APOSTRUCT — Output";
  document.title = (entry.snapshot.title || "APOSTRUCT") + " — Output";
  $("mode-details-section").setAttribute("aria-labelledby", "result-tab-" + entry.revision);
  $("output-status").hidden = true;
  setModeDetailsView(entry.view || "viewer");
  requestAnimationFrame(() => {
    if (displayedRevision !== entry.revision) return;
    if (entry.ui) {
      // Layout may have changed while the tab was hidden; preserve its camera and scale.
      modeViewer.restoreState(entry.ui.viewer);
      const scroll = $("mode-viewer-controls").querySelector(".mode-controls-scroll");
      if (scroll) scroll.scrollTop = entry.ui.controlsScroll;
      [$("mode-details-output").scrollLeft, $("mode-details-output").scrollTop] = entry.ui.textScroll;
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
      for (const entry of message.entries || [message.entry]) if (entry?.snapshot) tabs.add(entry);
      if (empty && message.activeRevision !== undefined) tabs.select(message.activeRevision);
      if (displayedRevision !== tabs.activeRevision) displayActive();
      else renderTabs();
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
for (const id of ["mode-viewer-bonds", "mode-viewer-polyhedra", "mode-details-view-select"]) {
  $(id).addEventListener("change", saveActive);
}
window.addEventListener("pagehide", saveActive);
document.addEventListener("visibilitychange", () => {
  modeViewer.setVisible(!document.hidden && Boolean(tabs.active) && $("mode-details-view-select").value === "viewer");
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
