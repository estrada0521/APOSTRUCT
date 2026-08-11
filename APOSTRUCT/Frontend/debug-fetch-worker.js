self.onmessage = async event => {
  if (event.data?.type === "mode-detail-text") {
    self.postMessage({ type: "mode-detail-text", text: self.modeDetailText || "" });
    return;
  }
  const { url, payload } = event.data || {};
  try {
    const response = await fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const data = await response.json();
    if (!response.ok) throw new Error(data?.error || response.statusText);
    self.modeDetailText = String(data?.state?.selected?.mode_detail_text || "");
    if (data?.state?.selected) delete data.state.selected.mode_detail_text;
    self.postMessage({ type: "result", ok: true, data });
  } catch (error) {
    self.postMessage({
      ok: false,
      error: {
        name: error?.name || "Error",
        message: error?.message || String(error),
      },
    });
  }
};
