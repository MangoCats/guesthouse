/* ========== Catalog & Placement Tool (TL-18, TL-19, TL-20, TL-21) ========== */
"use strict";

/** Cached catalog data from server (populated by loadCatalog). */
let _catalogData = null;

/** Placement tool state. */
const PlaceTool = {
  active: false,
  itemTemplate: null,
  itemType: null,
  previewEl: null, // reserved for future hover-preview during placement
};


/** Fetch catalog from server, caching result. */
async function loadCatalog(force) {
  if (_catalogData && !force) return _catalogData;
  try {
    const resp = await apiFetch("/api/catalog");
    _catalogData = await resp.json();
  } catch (e) {
    showToast("Failed to load catalog: " + e.message, "error");
    _catalogData = { furniture: [], appliance: [], fixture: [] };
  }
  return _catalogData;
}


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


/** Format item dimensions for catalog display. */
function catalogDimText(item) {
  if (item.shape === "circle" && item.radius) {
    return `${fmtFtIn(item.radius * 2)} dia`;
  }
  return `${fmtFtIn(item.width)} × ${fmtFtIn(item.depth)}`;
}


/** Build a catalog card element for an item or group. */
function buildCatalogCard(item, itemType, onSelect) {
  const card = document.createElement("div");
  card.className = "catalog-item";

  const nameEl = document.createElement("span");
  nameEl.className = "catalog-item-label";
  nameEl.textContent = item.label;
  card.appendChild(nameEl);

  if (item.product_url) {
    const linkIcon = document.createElement("span");
    linkIcon.className = "catalog-link-icon";
    linkIcon.textContent = "\u2197"; // ↗
    linkIcon.title = "Has product link";
    card.appendChild(linkIcon);
  }

  const dims = document.createElement("small");
  dims.textContent = catalogDimText(item);
  card.appendChild(document.createElement("br"));
  card.appendChild(dims);

  if (item.shape && item.shape !== "rect") {
    const shapeTag = document.createElement("span");
    shapeTag.className = "catalog-shape-tag";
    shapeTag.textContent = item.shape;
    card.appendChild(shapeTag);
  }

  card.addEventListener("click", () => onSelect(item));
  return card;
}


/** Show the catalog dialog for a given item type. */
async function showCatalog(itemType) {
  const catalog = await loadCatalog();
  const items = catalog[itemType];
  if (!items || items.length === 0) {
    showToast(`No ${itemType} items in catalog`, "warning");
    return;
  }

  // Show all items from all variants (not filtered by current variant)
  const available = items;

  // Group items by label to create sub-menus for items like COUNTER
  const groups = new Map();
  for (const item of available) {
    const label = item.label;
    if (!groups.has(label)) groups.set(label, []);
    groups.get(label).push(item);
  }

  const grid = document.createElement("div");
  grid.className = "catalog-grid";

  for (const [label, groupItems] of groups) {
    if (groupItems.length === 1) {
      // Single item — direct placement
      const card = buildCatalogCard(groupItems[0], itemType, (item) => {
        Dialog.close();
        startPlacement(item, itemType);
      });
      grid.appendChild(card);
    } else {
      // Multiple items sharing a label — show sub-menu on click
      const card = document.createElement("div");
      card.className = "catalog-item catalog-group";

      const nameEl = document.createElement("span");
      nameEl.className = "catalog-item-label";
      nameEl.textContent = label;
      card.appendChild(nameEl);

      const countEl = document.createElement("small");
      countEl.textContent = `${groupItems.length} variants \u25B6`;
      card.appendChild(document.createElement("br"));
      card.appendChild(countEl);

      card.addEventListener("click", () => {
        Dialog.close();
        showCatalogSubMenu(label, groupItems, itemType);
      });
      grid.appendChild(card);
    }
  }

  Dialog.show({
    title: `Add ${itemType.charAt(0).toUpperCase() + itemType.slice(1)}`,
    customContent: grid,
    fields: [],
    onSubmit() {},
  });
}


/** Show sub-menu for items sharing a label (e.g. COUNTER variants). */
function showCatalogSubMenu(label, items, itemType) {
  const grid = document.createElement("div");
  grid.className = "catalog-grid";

  for (const item of items) {
    const card = document.createElement("div");
    card.className = "catalog-item";

    // Show the specific key name for disambiguation
    const nameEl = document.createElement("span");
    nameEl.className = "catalog-item-label";
    nameEl.textContent = item.key.replace(/_/g, " ");
    card.appendChild(nameEl);

    if (item.product_url) {
      const linkIcon = document.createElement("span");
      linkIcon.className = "catalog-link-icon";
      linkIcon.textContent = "\u2197";
      linkIcon.title = "Has product link";
      card.appendChild(linkIcon);
    }

    const dims = document.createElement("small");
    dims.textContent = catalogDimText(item);
    card.appendChild(document.createElement("br"));
    card.appendChild(dims);

    if (item.variants && item.variants.length < 3) {
      const varTag = document.createElement("span");
      varTag.className = "catalog-shape-tag";
      varTag.textContent = item.variants.join(", ");
      card.appendChild(varTag);
    }

    card.addEventListener("click", () => {
      Dialog.close();
      startPlacement(item, itemType);
    });
    grid.appendChild(card);
  }

  Dialog.show({
    title: `${label} Variants`,
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
  showToast(`Click to place ${item.label} (${item.key})`, "info");
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

  const properties = {
    source: "placed",
    center: [wx, wy],
    width: item.width,
    depth: item.depth,
    rotation: 0,
    poly: poly,
    shape: item.shape || "rect",
    catalog_key: item.key,
    label: item.label,
  };

  // Carry over metadata from catalog — but NOT variants: the catalog's
  // variants list controls ADD-menu visibility, not where placed instances
  // appear.  Placed instances always show in the current variant (set via
  // the element record's variant column, handled by the server).
  if (item.product_url) properties.product_url = item.product_url;
  if (item.door) properties.door = item.door;
  if (item.clearance) properties.clearance = item.clearance;
  if (item.stacked) properties.stacked = true;
  if (item.clip_to_inner) properties.clip_to_inner = true;

  const resp = await fetch("/api/elements", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      type: type,
      name: name,
      variant: App.state.variant,
      properties: properties,
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
