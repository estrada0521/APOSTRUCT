// Results live until explicitly closed. A new result never steals the current selection.
export class OutputTabs {
  constructor() { this.entries = new Map(); this.activeRevision = null; }
  get active() { return this.entries.get(this.activeRevision) || null; }
  add(entry) {
    if (this.entries.has(entry.revision)) return this.entries.get(entry.revision);
    const result = { amplitudes: [], masterAmplitude: 1, view: "viewer", ui: null, ...entry };
    this.entries.set(result.revision, result);
    if (this.activeRevision === null) this.activeRevision = result.revision;
    return result;
  }
  select(revision) {
    if (!this.entries.has(revision)) return false;
    this.activeRevision = revision;
    return true;
  }
  close(revision) {
    const keys = [...this.entries.keys()], index = keys.indexOf(revision);
    if (index < 0) return false;
    this.entries.delete(revision);
    if (this.activeRevision === revision) this.activeRevision = keys[index + 1] ?? keys[index - 1] ?? null;
    return true;
  }
}
