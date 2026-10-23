// One output window per input page. Only serializable result data crosses the boundary.
export class OutputWindowHost {
  constructor(status) {
    this.status = status;
    this.channel = crypto.randomUUID();
    this.popup = null;
    this.snapshot = null;
    this.revision = 0;
    this.amplitudes = [];
    this.masterAmplitude = 1;
    this.view = "viewer";
    this.stale = false;
    window.addEventListener("message", event => {
      if (event.origin !== location.origin || event.source !== this.popup) return;
      const message = event.data;
      if (!message || message.app !== "apostruct-output" || message.channel !== this.channel) return;
      if (message.type === "ready") this.publish();
      else if (message.type === "controls" && message.revision === this.revision) {
        this.amplitudes = message.amplitudes;
        this.masterAmplitude = message.masterAmplitude;
      } else if (message.type === "view") this.view = message.view;
    });
    window.addEventListener("pagehide", () => this.send("input-unavailable"));
    window.addEventListener("pageshow", () => this.publish());
  }
  send(type, detail = {}) {
    if (this.popup && !this.popup.closed) {
      this.popup.postMessage({ app: "apostruct-output", channel: this.channel, type, ...detail }, location.origin);
    }
  }
  publish() {
    if (this.snapshot) this.send("snapshot", {
      snapshot: this.snapshot, revision: this.revision, amplitudes: this.amplitudes,
      masterAmplitude: this.masterAmplitude, view: this.view, stale: this.stale,
    });
  }
  setData(snapshot, view = this.view) {
    this.snapshot = snapshot;
    this.view = view;
    this.revision++;
    this.amplitudes = [];
    this.masterAmplitude = 1;
    this.stale = false;
    this.publish();
  }
  setAmplitudes(amplitudes) {
    this.amplitudes = amplitudes;
    this.send("amplitudes", { amplitudes, revision: this.revision });
  }
  markStale() {
    this.stale = true;
    this.send("stale");
  }
  open() {
    if (!this.snapshot) return;
    if (!this.popup || this.popup.closed) {
      const url = new URL("./output.html", import.meta.url);
      url.hash = this.channel;
      this.popup = window.open(url, "apostruct-output-" + this.channel, "popup,width=1200,height=850");
      if (!this.popup) {
        this.status.textContent = "The browser blocked the output window. Allow popups for this page, then try again.";
        return;
      }
    }
    this.popup.focus();
    this.publish();
    this.status.textContent = "Output opens in a separate window or tab.";
  }
}
