export const favoriteStorageKey = "apostruct.favorite-inputs.v1";

export function loadFavorites() {
  try {
    const entries = JSON.parse(localStorage.getItem(favoriteStorageKey) || "[]");
    return Array.isArray(entries) ? entries.filter(entry =>
      typeof entry?.key === "string" && typeof entry?.title === "string" &&
      typeof entry?.launch?.state_request?.cif === "string"
    ) : [];
  } catch { return []; }
}

export function isFavorite(entry) {
  return !!entry && loadFavorites().some(item => item.key === entry.key);
}

export function toggleFavorite(entry) {
  if (!entry) return;
  const entries = loadFavorites();
  const next = entries.some(item => item.key === entry.key)
    ? entries.filter(item => item.key !== entry.key)
    : [entry, ...entries];
  localStorage.setItem(favoriteStorageKey, JSON.stringify(next));
}
