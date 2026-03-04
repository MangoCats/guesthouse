/* ========== Catalog & Placement Tool (TL-18, TL-19, TL-20, TL-21) ========== */
"use strict";

/** Hardcoded catalog items (mirrors floorplan constants, per NF-4 pattern). */
const CATALOG = {
  furniture: [
    { key: "bed",      label: "King Bed",     width: 76.0/12, depth: 80.0/12 },
    { key: "dresser",  label: "Dresser",      width: 34.0/12, depth: 18.0/12 },
    { key: "shelves",  label: "Shelves",      width: 16.0/12, depth: 16.0/12 },
    { key: "loveseat", label: "Loveseat",     width: 35.0/12, depth: 72.0/12 },
    { key: "desk",     label: "Desk",         width: 60.0/12, depth: 24.0/12 },
    { key: "chair",    label: "Chair",        width: 32.0/12, depth: 33.0/12 },
    { key: "sofa",     label: "Sofa",         width: 97.2/12, depth: 38.0/12 },
    { key: "rocker",   label: "Rocking Chair",width: 26.75/12,depth: 32.0/12 },
  ],
  appliance: [
    { key: "washer",       label: "Washer",       width: 35.0/12, depth: 31.0/12 },
    { key: "dryer",        label: "Dryer",        width: 35.0/12, depth: 31.0/12 },
    { key: "stove",        label: "Stove",        width: 30.0/12, depth: 25.0/12 },
    { key: "dishwasher",   label: "Dishwasher",   width: 28.0/12, depth: 24.0/12 },
    { key: "ice_maker",    label: "Ice Maker",    width: 17.7/12, depth: 23.0/12 },
    { key: "kitchen_sink", label: "Kitchen Sink",  width: 45.0/12, depth: 22.0/12 },
  ],
  fixture: [
    { key: "toilet",   label: "Toilet",   width: 15.0/12, depth: 28.0/12 },
    { key: "tub",      label: "Bathtub",  width: 30.0/12, depth: 60.0/12 },
    { key: "vanity",   label: "Vanity",   width: 24.0/12, depth: 21.0/12 },
  ],
};


/** Placement tool state. */
const PlaceTool = {
  active: false,
  itemTemplate: null,
  itemType: null,
  previewEl: null, // reserved for future hover-preview during placement
};


/** Auto-generate a unique name for a placed item. */
function nextPlacedName(type, key) {
  const prefix = `custom_${key}`;
  const elements = App.state.elements || [];
  let max = 0;
  for (const e of elements) {
    if (e.name && e.name.startsWith(prefix)) {
      const suffix = e.name.slice(prefix.length);
      const n = suffix ? parseInt(suffix, 10) : 1;
      if (n > max) max = n;
    }
  }
  return max === 0 ? prefix : `${prefix}${max + 1}`;
}


/** Compute a centered rect polygon from position, width, depth. */
function rectPoly(cx, cy, w, d) {
  const hw = w / 2, hd = d / 2;
  return [
    [cx - hw, cy + hd],
    [cx + hw, cy + hd],
    [cx + hw, cy - hd],
    [cx - hw, cy - hd],
  ];
}


/** Show the catalog dialog for a given item type. */
function showCatalog(itemType) {
  const items = CATALOG[itemType];
  if (!items || items.length === 0) {
    showToast(`No ${itemType} items in catalog`, "warning");
    return;
  }

  const grid = document.createElement("div");
  grid.className = "catalog-grid";

  for (const item of items) {
    const card = document.createElement("div");
    card.className = "catalog-item";
    card.textContent = item.label;
    const dims = document.createElement("small");
    dims.textContent = `${fmtFtIn(item.width)} x ${fmtFtIn(item.depth)}`;
    card.appendChild(document.createElement("br"));
    card.appendChild(dims);
    card.addEventListener("click", () => {
      Dialog.close();
      startPlacement(item, itemType);
    });
    grid.appendChild(card);
  }

  Dialog.show({
    title: `Add ${itemType.charAt(0).toUpperCase() + itemType.slice(1)}`,
    customContent: grid,
    fields: [],
    onSubmit() {},
  });
}


/** Enter placement mode: next canvas click places the item. */
function startPlacement(item, type) {
  PlaceTool.active = true;
  PlaceTool.itemTemplate = item;
  PlaceTool.itemType = type;
  App.els["viewport"].style.cursor = "copy";
  showToast(`Click to place ${item.label}`, "info");
}


/** Handle canvas click during placement mode. */
function placeToolMouseDown(e) {
  if (e.button !== 0 || !PlaceTool.active) return;

  const rect = App.els["viewport"].getBoundingClientRect();
  let [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  // Grid snap
  if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    wx = Math.round(wx / snap) * snap;
    wy = Math.round(wy / snap) * snap;
  }

  placeElement(wx, wy);
  PlaceTool.active = false;
  App.els["viewport"].style.cursor = "crosshair";
}


/** Create the placed element via API. */
async function placeElement(wx, wy) {
  const item = PlaceTool.itemTemplate;
  const type = PlaceTool.itemType;
  const name = nextPlacedName(type, item.key);
  const poly = rectPoly(wx, wy, item.width, item.depth);

  const resp = await fetch("/api/elements", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      type: type,
      name: name,
      variant: App.state.variant,
      properties: {
        source: "placed",
        center: [wx, wy],
        width: item.width,
        depth: item.depth,
        rotation: 0,
        poly: poly,
        shape: "rect",
        catalog_key: item.key,
      },
    }),
  });

  if (resp.ok) {
    showToast(`Placed ${item.label}`, "success");
    await loadElements();
    await loadGeometry();
  } else {
    const data = await resp.json();
    showToast(data.error || "Placement failed", "error");
  }
}


/** Cancel placement mode. */
function cancelPlacement() {
  PlaceTool.active = false;
  PlaceTool.itemTemplate = null;
  PlaceTool.itemType = null;
  if (PlaceTool.previewEl && PlaceTool.previewEl.parentNode) {
    PlaceTool.previewEl.remove();
  }
  PlaceTool.previewEl = null;
}
