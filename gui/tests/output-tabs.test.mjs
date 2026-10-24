import test from "node:test";
import assert from "node:assert/strict";
import { OutputTabs } from "../web/output-tabs.js";
import { OutputWindowHost } from "../web/output-window.js";

test("appending results preserves selection and independent UI state", () => {
  const tabs = new OutputTabs();
  tabs.add({ revision: 1, snapshot: { title: "A" } });
  tabs.active.ui = { camera: [1, 2, 3], amplitudes: [0.5] };
  tabs.add({ revision: 2, snapshot: { title: "B" } });
  assert.equal(tabs.activeRevision, 1);
  tabs.select(2);
  assert.equal(tabs.active.ui, null);
  tabs.active.ui = { camera: [4, 5, 6], amplitudes: [0.9] };
  tabs.select(1);
  assert.deepEqual(tabs.active.ui, { camera: [1, 2, 3], amplitudes: [0.5] });
  tabs.add({ revision: 1, snapshot: { title: "stale duplicate" } });
  assert.equal(tabs.entries.size, 2);
  assert.equal(tabs.active.snapshot.title, "A");
});

test("closing an inactive tab does not change selection; closing active chooses a neighbor", () => {
  const tabs = new OutputTabs();
  for (let revision = 1; revision <= 4; revision++) tabs.add({ revision });
  tabs.select(2);
  tabs.close(4); assert.equal(tabs.activeRevision, 2);
  tabs.close(2); assert.equal(tabs.activeRevision, 3);
  tabs.close(3); assert.equal(tabs.activeRevision, 1);
  tabs.close(1); assert.equal(tabs.active, null);
  tabs.add({ revision: 5 }); assert.equal(tabs.activeRevision, 5);
  assert.equal(tabs.entries.has(2), false);
});

test("host retains all results and routes CLI amplitudes to their own result", () => {
  const listeners = new Map(), messages = [];
  const popup = { closed: false, focus() {}, postMessage(message) { messages.push(message); } };
  globalThis.window = { addEventListener(type, handler) { listeners.set(type, handler); }, open() { return popup; } };
  globalThis.location = { origin: "http://localhost" };
  const host = new OutputWindowHost({ textContent: "" });
  host.setData({ title: "A" });
  host.setAmplitudes([{ definition_id: "displacive-1", value: "1/2" }]);
  host.open();
  host.setData({ title: "B" });
  host.setAmplitudes([{ definition_id: "displacive-1", value: "3/4" }]);
  assert.equal(host.tabs.activeRevision, 1);
  assert.equal(host.tabs.entries.get(1).amplitudes[0].value, "1/2");
  assert.equal(host.tabs.entries.get(2).amplitudes[0].value, "3/4");
  const receive = (type, detail = {}, origin = location.origin) => listeners.get("message")({
    origin, source: popup, data: { app: "apostruct-output", channel: host.channel, type, ...detail },
  });
  receive("close", { revision: 1 }, "https://untrusted.invalid");
  assert.equal(host.tabs.entries.size, 2);
  receive("select", { revision: 2 });
  receive("result-state", { revision: 2, state: { ui: { camera: [4, 5, 6] } } });
  host.markStale();
  assert.equal(host.tabs.entries.get(2).ui.camera[0], 4);
  receive("close", { revision: 1 });
  receive("ready");
  const restored = messages.at(-1);
  assert.equal(restored.type, "results");
  assert.deepEqual(restored.entries.map(e => e.revision), [2]);
  assert.equal(restored.activeRevision, 2);
  assert.deepEqual(restored.entries[0].ui.camera, [4, 5, 6]);
});
