import { ModeStructureViewer } from "./mode3d.js?v=62";

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
  if (showViewer) requestAnimationFrame(() => modeViewer.refreshLayout());
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
let revision = -1;

let hasResult = false;
function send(type, detail = {}) {
  if (source && !source.closed) source.postMessage({ app: "apostruct-output", channel, type, ...detail }, location.origin);
}
function outputError(error) {
  $("output-error").hidden = !error;
  $("output-error").textContent = error?.message || String(error || "");
}
function resultStatus(stale = false) {
  $("output-status").hidden = !stale;
  $("output-status").textContent = stale
    ? "Input changed — showing the previous calculated result."
    : "Showing the last calculated result.";
}
window.addEventListener("message", event => {
  if (event.origin !== location.origin || event.source !== source) return;
  const message = event.data;
  if (!message || message.app !== "apostruct-output" || message.channel !== channel) return;
  if (message.type === "snapshot") {
    if (!message.snapshot || message.revision < revision) return;
    try {
      const isNew = message.revision !== revision;
      const snapshot = message.snapshot;
      if (isNew) {
        outputError(null);
        state.modeDetailText = snapshot.text || "";
        state.modeDetailTextRendered = false;
        $("mode-details-output").replaceChildren();
        modeViewer.setData(snapshot.viewer);
        if (message.amplitudes?.length) modeViewer.setAmplitudes(message.amplitudes);
        if (Number.isFinite(message.masterAmplitude)) modeViewer.masterAmplitude = message.masterAmplitude;
        modeViewer.updateControlValues();
        modeViewer.updateScene();
        $("result-title").textContent = snapshot.title || "APOSTRUCT — Output";
        document.title = (snapshot.title || "APOSTRUCT") + " — Output";
        setModeDetailsView(message.view || "viewer");
        revision = message.revision;
        hasResult = true;
      }
      resultStatus(message.stale);
    } catch (error) { outputError(error); }
  } else if (message.type === "amplitudes" && message.revision === revision) {
    try { modeViewer.setAmplitudes(message.amplitudes); } catch (error) { outputError(error); }
  } else if (message.type === "stale") {
    if (hasResult) resultStatus(true);
  } else if (message.type === "input-unavailable") {
    disconnected();
  }
});
function captureControls() {
  queueMicrotask(() => send("controls", {
    revision,
    amplitudes: [...modeViewer.amplitudes].map(([definition_id, value]) => ({ definition_id, value })),
    masterAmplitude: modeViewer.masterAmplitude,
  }));
}
$("mode-viewer-controls").addEventListener("input", captureControls);
$("mode-viewer-controls").addEventListener("click", event => {
  if (event.target.closest("button[data-mode-reset]")) captureControls();
});
window.addEventListener("pagehide", captureControls);
$("mode-details-view-select").addEventListener("change", event => send("view", { view: event.target.value }));
function disconnected() {
  $("output-status").hidden = false;
  $("output-status").textContent = hasResult
    ? "Input window disconnected — this result remains available. Reopen output from the input window to reconnect."
    : "Input window disconnected. Open output from the input window.";
}
if (!source || !channel) {
  $("output-status").hidden = false;
  $("output-status").textContent = "Open this page using “Open output” in the input window.";
} else {
  send("ready");
  setInterval(() => { if (source.closed) disconnected(); }, 2000);
}
