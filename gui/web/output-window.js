import { OutputTabs } from "./output-tabs.js?v=1";

// One output window, many independent results. The input page keeps them across popup reloads.
export class OutputWindowHost {
  constructor(status) {
    this.status = status;
    this.channel = crypto.randomUUID();
    this.popup = null;
    this.tabs = new OutputTabs();
    this.revision = 0;
    window.addEventListener("message", event => {
      if (event.origin !== location.origin || event.source !== this.popup) return;
      const message = event.data;
      if (!message || message.app !== "apostruct-output" || message.channel !== this.channel) return;
      if (message.type === "ready") this.publish();
      else if (message.type === "result-state") {
        const entry = this.tabs.entries.get(message.revision);
        if (entry) Object.assign(entry, message.state);
      } else if (message.type === "select") this.tabs.select(message.revision);
      else if (message.type === "close") this.tabs.close(message.revision);
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
    this.send("results", { entries: [...this.tabs.entries.values()], activeRevision: this.tabs.activeRevision });
  }
  setData(snapshot, view = "viewer") {
    const entry = this.tabs.add({ revision: ++this.revision, snapshot, view: view || "viewer" });
    this.tabs.select(entry.revision);
    this.send("result", { entry });
  }
  setAmplitudes(amplitudes) {
    const entry = this.tabs.entries.get(this.revision);
    if (!entry) return;
    entry.amplitudes = amplitudes;
    this.send("amplitudes", { amplitudes, revision: entry.revision });
  }
  // Editing input does not invalidate an independently retained result.
  markStale() {}
  open() {
    if (!this.tabs.entries.size && (!this.popup || this.popup.closed)) return;
    if (!this.popup || this.popup.closed) {
      const url = new URL("./output.html", import.meta.url);
      url.hash = this.channel;
      this.popup = window.open(url, "apostruct-output-" + this.channel, "popup,width=800,height=567");
      if (!this.popup) {
        this.status.textContent = "The browser blocked the output window. Allow popups for this page, then try again.";
        return;
      }
    }
    this.popup.focus();
    this.publish();
    this.status.textContent = "Results remain in output tabs until closed.";
  }
}
