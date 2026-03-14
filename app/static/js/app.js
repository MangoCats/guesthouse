/* ========== ADU Editor — Main Application ========== */
"use strict";

const App = {
  state: {
    constants: [],
    geometry: null,
    views: [],
    activeView: "interactive",
    activeTool: "select",
    selection: null,
    selections: [],
    pan: { x: 0, y: 0 },
    zoom: 1,
    showPoints: true,
    showLabels: true,
    showDims: false,
    showUserDims: true,
    showGrid: false,
    showOpenings: true,
    showFurniture: true,
    showRooms: true,
    showDoors: true,
    showClearance: false,
    openLinks: true,
    showAreas: false,
    variant: "standard",
    measureStart: null,
    rubberBand: null,
    isDragging: false,
    dragStart: null,
    lastPan: null,
    constSortKey: "name",
    constSortAsc: true,
    svgView: { pan: { x: 0, y: 0 }, zoom: 1 },
    outlineChain: [],
    outlineSelectedSeq: null,
    innerWallOverrides: {},
    plumbingElements: [],
    variants: [],
    configName: null,
    configDirty: false,
  },
  els: {},
  sse: null,
};


/* ========== INITIALISATION ========== */

document.addEventListener("DOMContentLoaded", async () => {
  cacheElements();
  document.querySelectorAll(".plumb-swatch").forEach(el => {
    const c = PLUMBING_COLORS[el.dataset.colorKey];
    if (c) { el.setAttribute("stroke", c); el.setAttribute("fill", c); }
  });
  setupEventListeners();
  setupOutlineToolbar();
  connectSSE();
  loadViews();
  loadConstants();
  await loadElements();  // must complete before first render
  loadPlumbingElements();
  await loadVariants();  // populate variant dropdown before geometry
  updateDeleteVariantBtn();
  loadGeometry();
  loadShapes();
  loadBuildLabel();
  loadRoofStyle();
  loadConfigName();
  loadCatalog();
});

function cacheElements() {
  const ids = [
    "canvas", "canvas-transform", "viewport",
    "layer-outline", "layer-inner", "layer-walls",
    "layer-openings", "layer-doors", "layer-furniture", "layer-clearance",
    "layer-rooms", "layer-points",
    "layer-labels", "layer-dims",
    "layer-plumbing-pipes", "layer-plumbing-fittings", "layer-plumbing-fixtures",
    "layer-measure",
    "grid-rect", "svg-view-container",
    "coord-display", "connection-status", "zoom-level",
    "selection-info", "measure-info",
    "view-tabs", "const-category-filter", "const-search",
    "constants-table", "outline-table", "outline-closure-indicator",
    "outline-add-btn", "outline-remove-btn", "openings-table",
    "rough-openings-table", "interior-walls-table", "furniture-table",
    "props-empty", "props-detail", "props-title", "props-table",
    "show-points", "show-labels", "show-dims", "show-user-dims", "show-grid",
    "show-openings", "show-furniture", "show-rooms",
    "show-doors", "show-clearance", "open-links", "show-areas", "roof-style",
    "plumbing-tools", "plumbing-fixtures-table", "plumbing-pipes-table",
    "variant-select", "variant-selector",
    "error-banner", "error-banner-text", "error-banner-action", "error-banner-dismiss",
    "config-name",
  ];
  for (const id of ids) {
    App.els[id] = document.getElementById(id);
  }
}


/* ========== SSE (Server-Sent Events) ========== */

function connectSSE() {
  if (App.sse) App.sse.close();
  App.sse = new EventSource("/api/events");

  App.sse.addEventListener("connected", () => {
    App.els["connection-status"].textContent = "Connected";
    App.els["connection-status"].className = "connected";
  });

  App.sse.addEventListener("constants_changed", () => {
    loadConstants();
    markConfigDirty();
  });

  App.sse.addEventListener("element_changed", async () => {
    markConfigDirty();
    await loadElements();
    if (App.state.geometry) {
      if (App.state.activeView === "interactive") {
        if (App.state.variant === "plumbing") renderPlumbingCanvas();
        else renderCanvas();
      }
    }
  });

  App.sse.addEventListener("formula_locked", () => {
    loadGeometry();
  });

  App.sse.addEventListener("plumbing_changed", () => {
    markConfigDirty();
    loadPlumbingElements();
  });

  App.sse.addEventListener("geometry_changed", () => {
    _spanDataCache = null;
    _spanRotationCache = null;
    loadGeometry();
  });

  App.sse.addEventListener("undo_status", (e) => {
    const data = JSON.parse(e.data);
    updateUndoButtons(data.can_undo, data.can_redo);
  });

  App.sse.addEventListener("outline_changed", () => {
    markConfigDirty();
    loadOutlineTable();
  });

  App.sse.addEventListener("variant_changed", async () => {
    markConfigDirty();
    await loadVariants();
    ensureActiveVariantValid();
  });

  App.sse.addEventListener("svg_updated", (e) => {
    const data = JSON.parse(e.data);
    if (data.view && App.state.activeView === data.view) {
      loadSVGView(data.view);
    }
    showToast(`View "${data.view}" updated`, "success");
  });

  App.sse.onerror = () => {
    App.els["connection-status"].textContent = "Disconnected";
    App.els["connection-status"].className = "disconnected";
    setTimeout(connectSSE, 3000);
  };
}


/* ========== API HELPERS ========== */

/** Parse element properties from string or object, always returns an object. */
function parseProps(elemOrProps) {
  if (!elemOrProps) return {};
  const p = elemOrProps.properties !== undefined ? elemOrProps.properties : elemOrProps;
  return typeof p === "string" ? JSON.parse(p) : (p || {});
}

async function apiFetch(url, opts) {
  const resp = await fetch(url, opts);
  if (!resp.ok) {
    const body = await resp.text();
    let msg;
    try { msg = JSON.parse(body).error; } catch (_) { msg = body; }
    throw new Error(msg || `HTTP ${resp.status}`);
  }
  return resp;
}


/* ========== DATA LOADING ========== */

async function loadBuildLabel() {
  try {
    const resp = await fetch("/api/version");
    const data = await resp.json();
    const el = document.getElementById("build-label");
    if (!el) return;
    const text = `${data.git} | ${data.started}`;
    el.textContent = text;
    const btn = document.getElementById("build-copy");
    if (btn) {
      btn.addEventListener("click", () => {
        navigator.clipboard.writeText(text).then(
          () => { btn.textContent = "\u2713"; setTimeout(() => { btn.textContent = "\u2398"; }, 1200); },
          () => { btn.textContent = "!"; setTimeout(() => { btn.textContent = "\u2398"; }, 1200); }
        );
      });
    }
  } catch (_) {}
}

async function loadRoofStyle() {
  try {
    const resp = await apiFetch("/api/config/roof_style");
    const data = await resp.json();
    if (App.els["roof-style"]) App.els["roof-style"].value = data.value;
  } catch (_) {}
}

async function loadConstants() {
  try {
    const resp = await apiFetch("/api/constants");
    App.state.constants = await resp.json();
    renderConstantsTable();
    loadCategories();
  } catch (e) {
    console.error("Constants load failed:", e);
  }
}

async function loadCategories() {
  const resp = await apiFetch("/api/constants/categories");
  const cats = await resp.json();
  const sel = App.els["const-category-filter"];
  sel.innerHTML = '<option value="">All Categories</option>';
  for (const cat of cats) {
    const opt = document.createElement("option");
    opt.value = cat;
    opt.textContent = cat;
    sel.appendChild(opt);
  }
}

async function loadGeometry() {
  try {
    const resp = await fetch(`/api/geometry?variant=${App.state.variant}`);
    if (!resp.ok) {
      let data;
      try { data = await resp.json(); } catch (_) { data = {}; }
      if (data.db_issue) {
        showErrorBanner(
          data.hint || "Database error — geometry computation failed.",
          "Reset Database", resetDatabase
        );
      }
      throw new Error(data.error || `HTTP ${resp.status}`);
    }
    App.state.geometry = await resp.json();
    hideErrorBanner();
    if (App.state.activeView === "interactive") {
      if (App.state.variant === "plumbing") renderPlumbingCanvas();
      else renderCanvas();
    }
    updateOpeningsTable();
    updateElementsTable();
  } catch (e) {
    console.error("Geometry load failed:", e);
    showToast("Geometry computation error: " + e.message, "error");
  }
}

async function reloadAfterChange() {
  await reloadAfterChange();
}

async function loadViews() {
  try {
    const resp = await apiFetch("/api/views");
    App.state.views = await resp.json();
    renderViewTabs();
  } catch (e) {
    console.error("Views load failed:", e);
  }
}

async function loadVariants() {
  try {
    const resp = await fetch("/api/variants");
    App.state.variants = await resp.json();
    populateVariantSelect();
  } catch (e) {
    console.error("Variants load failed:", e);
  }
}

function populateVariantSelect() {
  const sel = App.els["variant-select"];
  if (!sel) return;
  sel.innerHTML = "";
  for (const v of App.state.variants) {
    const opt = document.createElement("option");
    opt.value = v.name;
    opt.textContent = v.label;
    sel.appendChild(opt);
  }
  sel.value = App.state.variant;
}

function variantLabel(name) {
  const v = (App.state.variants || []).find(v => v.name === name);
  return v ? v.label : name;
}

async function saveCurrentLayerConfig() {
  const v = (App.state.variants || []).find(v => v.name === App.state.variant);
  if (!v) return;
  const cfg = {
    points: App.state.showPoints, labels: App.state.showLabels,
    dims: App.state.showDims, userDims: App.state.showUserDims,
    grid: App.state.showGrid,
    openings: App.state.showOpenings, furniture: App.state.showFurniture,
    rooms: App.state.showRooms, doors: App.state.showDoors,
    clearance: App.state.showClearance, areas: App.state.showAreas,
  };
  try {
    await apiFetch(`/api/variants/${v.id}`, {
      method: "PUT", headers: {"Content-Type": "application/json"},
      body: JSON.stringify({layer_config: cfg}),
    });
  } catch (e) {
    console.error("Failed to save layer config:", e);
  }
}

function restoreLayerConfig() {
  const v = (App.state.variants || []).find(v => v.name === App.state.variant);
  const cfg = v?.layer_config || {};
  const map = {
    points: "showPoints", labels: "showLabels", dims: "showDims",
    userDims: "showUserDims", grid: "showGrid",
    openings: "showOpenings", furniture: "showFurniture",
    rooms: "showRooms", doors: "showDoors",
    clearance: "showClearance", areas: "showAreas",
  };
  for (const [cfgKey, stateKey] of Object.entries(map)) {
    const val = cfg[cfgKey] !== undefined ? cfg[cfgKey] : App.state[stateKey];
    App.state[stateKey] = val;
    const elKeyMap = { userDims: "show-user-dims" };
    const elKey = elKeyMap[cfgKey] || (cfgKey === "points" ? "show-points" : `show-${cfgKey}`);
    const el = App.els[elKey];
    if (el) el.checked = val;
  }
}

function ensureActiveVariantValid() {
  const names = (App.state.variants || []).map(v => v.name);
  if (!names.includes(App.state.variant)) {
    App.state.variant = "standard";
    populateVariantSelect();
    restoreLayerConfig();
    loadGeometry();
  }
  updateDeleteVariantBtn();
}

function updateDeleteVariantBtn() {
  const btn = document.getElementById("delete-variant-btn");
  if (!btn) return;
  const v = (App.state.variants || []).find(v => v.name === App.state.variant);
  btn.style.display = (v && !v.is_builtin) ? "inline-block" : "none";
}

async function createNewVariant() {
  const name = prompt("New variant name:");
  if (!name || !name.trim()) return;
  const slug = name.trim().toLowerCase().replace(/[^a-z0-9_]/g, "_");
  try {
    const resp = await apiFetch("/api/variants", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({
        name: slug,
        label: name.trim(),
        source_variant: App.state.variant,
      }),
    });
    await loadVariants();
    App.state.variant = slug;
    populateVariantSelect();
    restoreLayerConfig();
    loadGeometry();
    updateDeleteVariantBtn();
    showToast(`Created variant "${name.trim()}"`, "success");
  } catch (e) {
    showToast("Failed to create variant", "error");
  }
}

async function deleteActiveVariant() {
  const v = (App.state.variants || []).find(v => v.name === App.state.variant);
  if (!v || v.is_builtin) return;
  if (!confirm(`Delete variant "${v.label}"?`)) return;
  try {
    await apiFetch(`/api/variants/${v.id}`, { method: "DELETE" });
    await loadVariants();
    App.state.variant = "standard";
    populateVariantSelect();
    restoreLayerConfig();
    loadGeometry();
    updateDeleteVariantBtn();
    showToast(`Deleted variant "${v.label}"`, "success");
  } catch (e) {
    showToast("Failed to delete variant", "error");
  }
}


/* ========== VIEW TABS ========== */

function renderViewTabs() {
  const container = App.els["view-tabs"];
  // Preserve variant selector before clearing (it's a static child of view-tabs)
  const vs = App.els["variant-selector"];
  container.innerHTML = "";

  // Interactive view tab (always first)
  const interTab = document.createElement("button");
  interTab.className = "view-tab active";
  interTab.textContent = "Interactive";
  interTab.dataset.view = "interactive";
  interTab.onclick = () => switchView("interactive");
  container.appendChild(interTab);

  // Floorplan tab (always second, right after Interactive)
  const fp = App.state.views.find(v => v.name === "floorplan");
  if (fp) {
    const fpTab = document.createElement("button");
    fpTab.className = "view-tab";
    fpTab.textContent = fp.label;
    fpTab.dataset.view = fp.name;
    fpTab.onclick = () => switchView(fp.name);
    container.appendChild(fpTab);
  }

  // Variant selector right after Floorplan
  container.appendChild(vs);

  // Remaining generated SVG view tabs
  for (const v of App.state.views) {
    if (v.name === "floorplan" || v.name === "plumbing") continue;
    const tab = document.createElement("button");
    tab.className = "view-tab";
    tab.textContent = v.label;
    tab.dataset.view = v.name;
    tab.onclick = () => switchView(v.name);
    container.appendChild(tab);
  }

}

function isCanvasView(name) {
  const v = name || App.state.activeView;
  return v === "interactive";
}

function switchView(viewName) {
  App.state.activeView = viewName;
  document.querySelectorAll(".view-tab").forEach(t => {
    t.classList.toggle("active", t.dataset.view === viewName);
  });

  // Show variant selector for interactive and floorplan views
  const showVariant = viewName === "interactive" || viewName === "floorplan";
  App.els["variant-selector"].style.display = showVariant ? "inline-block" : "none";

  // Show/hide plumbing tools (visible when plumbing variant is active)
  if (App.els["plumbing-tools"]) {
    const isPlumbing = App.state.variant === "plumbing" &&
                       (viewName === "interactive" || viewName === "floorplan");
    App.els["plumbing-tools"].style.display = isPlumbing ? "" : "none";
  }

  if (isCanvasView(viewName)) {
    App.els["canvas"].style.display = "block";
    App.els["svg-view-container"].style.display = "none";
    App.els["viewport"].style.cursor = App.state.activeTool === "pan" ? "grab" : "crosshair";
    if (App.state.variant === "plumbing") renderPlumbingCanvas();
    else renderCanvas();
  } else {
    App.els["canvas"].style.display = "none";
    App.els["svg-view-container"].style.display = "block";
    App.els["viewport"].style.cursor = "grab";
    const viewDef = App.state.views.find(v => v.name === viewName);
    if (viewDef && viewDef.svg_path.endsWith(".pdf")) {
      loadPDFView(viewName);
    } else if (viewDef && viewDef.svg_path.endsWith(".png")) {
      loadImageView(viewName);
    } else {
      loadSVGView(viewName);
    }
  }
}

async function loadSVGView(viewName) {
  const container = App.els["svg-view-container"];
  container.innerHTML = "<p style='padding:20px;color:#888'>Loading...</p>";
  try {
    let url = `/api/svg/${viewName}`;
    if (viewName === "floorplan") url += `?variant=${App.state.variant}`;
    const resp = await fetch(url);
    if (!resp.ok) {
      container.innerHTML = `<p style='padding:20px;color:#f88'>SVG not available. Click File > Regenerate to generate it.</p>`;
      return;
    }
    const svg = await resp.text();
    container.innerHTML = svg;
    // Initialise pan/zoom to fit the SVG in the viewport
    const svgEl = container.querySelector("svg");
    if (svgEl) {
      svgEl.style.transformOrigin = "0 0";
      svgViewFit();
    }
    // Set up span hover tooltip for analysis views
    if (["span", "span_min", "span_minmax"].includes(viewName)) {
      setupSpanTooltip(viewName);
    }
    // Show span rotation properties for span_minmax view
    if (viewName === "span_minmax") {
      showSpanRotationProperties();
    }
  } catch (e) {
    const p = document.createElement("p");
    p.style.cssText = "padding:20px;color:#f88";
    p.textContent = "Error loading SVG: " + e.message;
    container.innerHTML = "";
    container.appendChild(p);
  }
}

function loadPDFView(viewName) {
  const container = App.els["svg-view-container"];
  container.innerHTML = "";
  const embed = document.createElement("embed");
  embed.src = `/api/svg/${viewName}/file`;
  embed.type = "application/pdf";
  embed.style.cssText = "width:100%;height:100%;border:none";
  container.appendChild(embed);
  if (viewName.startsWith("site_plan")) {
    showSitePlanProperties();
  }
}

function loadImageView(viewName) {
  const container = App.els["svg-view-container"];
  container.innerHTML = "";
  const img = document.createElement("img");
  img.src = `/api/svg/${viewName}/file`;
  img.alt = viewName;
  img.style.cssText = "max-width:100%;max-height:100%;object-fit:contain;display:block;margin:auto";
  img.onerror = () => {
    container.innerHTML = '<p style="padding:20px;color:#f88">Image not available. Click File &gt; Generate 3D Model to render it.</p>';
  };
  container.appendChild(img);
}

function svgViewFit() {
  const container = App.els["svg-view-container"];
  const svgEl = container.querySelector("svg");
  if (!svgEl) return;
  const cw = container.clientWidth;
  const ch = container.clientHeight;
  // Use the SVG's intrinsic size from viewBox or width/height attributes
  const vb = svgEl.viewBox.baseVal;
  let sw, sh;
  if (vb && vb.width > 0) {
    sw = vb.width;
    sh = vb.height;
  } else {
    sw = svgEl.width.baseVal.value || cw;
    sh = svgEl.height.baseVal.value || ch;
  }
  const scale = Math.min(cw / sw, ch / sh) * 0.95;
  const px = (cw - sw * scale) / 2;
  const py = (ch - sh * scale) / 2;
  App.state.svgView = { pan: { x: px, y: py }, zoom: scale };
  svgViewApplyTransform();
}

function svgViewApplyTransform() {
  const container = App.els["svg-view-container"];
  const svgEl = container.querySelector("svg");
  if (!svgEl) return;
  const { pan, zoom } = App.state.svgView;
  svgEl.style.transform = `translate(${pan.x}px, ${pan.y}px) scale(${zoom})`;
  // Update zoom display
  App.els["zoom-level"].textContent = `${Math.round(zoom * 100)}%`;
}


/* ========== SVG VIEW TOOLTIP (ANALYSIS-1, ANALYSIS-2) ========== */

/** Cached span data (fetched once per session, invalidated on geometry change). */
let _spanDataCache = null;
let _spanRotationCache = null;

/** Ensure the tooltip div exists. */
function getOrCreateSvgTooltip() {
  let tip = document.getElementById("svg-tooltip");
  if (!tip) {
    tip = document.createElement("div");
    tip.id = "svg-tooltip";
    tip.className = "svg-tooltip";
    document.getElementById("viewport-area").appendChild(tip);
  }
  return tip;
}

/** Show span hover tooltip for the span SVG view. */
function setupSpanTooltip(viewName) {
  const container = App.els["svg-view-container"];
  const tip = getOrCreateSvgTooltip();

  // Fetch span data if not cached
  if (viewName === "span" && !_spanDataCache) {
    fetch("/api/span-data").then(r => r.json()).then(d => { _spanDataCache = d; });
  }

  container.addEventListener("mousemove", function spanHover(e) {
    if (App.state.activeView !== viewName) {
      tip.style.display = "none";
      container.removeEventListener("mousemove", spanHover);
      return;
    }
    const svgEl = container.querySelector("svg");
    if (!svgEl) return;
    const vb = svgEl.viewBox.baseVal;
    if (!vb || vb.width <= 0) return;

    // Translate mouse to SVG viewBox coords
    const { pan, zoom } = App.state.svgView;
    const rect = container.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const svgX = (mx - pan.x) / zoom;
    const svgY = (my - pan.y) / zoom;

    // Position tooltip near cursor
    tip.style.left = (mx + 12) + "px";
    tip.style.top = (my - 8) + "px";

    if (viewName === "span" && _spanDataCache) {
      // Map svgX to easting (span graph has a left margin and x-scale)
      // Approximate: find nearest easting by fraction of viewBox width
      const data = _spanDataCache;
      const frac = svgX / vb.width;
      const idx = Math.round(frac * (data.eastings.length - 1));
      if (idx >= 0 && idx < data.eastings.length) {
        const east = data.eastings[idx];
        const span = data.spans[idx];
        tip.textContent = `E: ${fmtFtIn(east)} | Span: ${fmtFtIn(span)}`;
        tip.style.display = "block";
      } else {
        tip.style.display = "none";
      }
    } else {
      // Generic SVG coordinate readout for other span views
      tip.textContent = `X: ${svgX.toFixed(0)} Y: ${svgY.toFixed(0)}`;
      tip.style.display = "block";
    }
  });

  container.addEventListener("mouseleave", () => {
    tip.style.display = "none";
  });
}

/** Show span rotation summary in properties panel. */
async function showSpanRotationProperties() {
  if (!_spanRotationCache) {
    const resp = await fetch("/api/span-rotation");
    _spanRotationCache = await resp.json();
  }
  const data = _spanRotationCache;
  const tbody = App.els["props-table"];
  tbody.innerHTML = "";
  App.els["props-title"].textContent = "Span vs Rotation";
  App.els["props-detail"].style.display = "block";
  App.els["props-empty"].style.display = "none";

  addPropRow(tbody, "Min Span", fmtFtIn(data.min_span));
  addPropRow(tbody, "At Rotation", data.min_angle + "\u00B0");
  addPropRow(tbody, "Max Span", fmtFtIn(data.max_span));
  addPropRow(tbody, "At Rotation", data.max_angle + "\u00B0");
  addPropRow(tbody, "Data Points", data.data.length + " angles");
}

/** Show site plan properties panel with setbacks and survey points. */
async function showSitePlanProperties() {
  const tbody = App.els["props-table"];
  tbody.innerHTML = "";
  App.els["props-title"].textContent = "Site Plan";
  App.els["props-detail"].style.display = "block";
  App.els["props-empty"].style.display = "none";

  // Setback config values (SITE-1)
  try {
    const r216 = await fetch("/api/config/setback_216");
    const d216 = await r216.json();
    const r275 = await fetch("/api/config/setback_275");
    const d275 = await r275.json();
    addConfigRow(tbody, "Setback from 216.73\u2032", "setback_216", d216.value);
    addConfigRow(tbody, "Setback from 275.08\u2032", "setback_275", d275.value);
  } catch (_) {}

  // Survey points (SITE-4)
  try {
    const resp = await fetch("/api/survey-points");
    const data = await resp.json();
    addPropRow(tbody, "\u2014", "Survey Points");
    for (const name of ["POB", "P2", "P3", "P4", "P5"]) {
      const pt = data.points[name];
      if (pt) addPropRow(tbody, name, `E: ${fmtFtIn(pt[0])}, N: ${fmtFtIn(pt[1])}`);
    }
    addPropRow(tbody, "\u2014", "Traverse Distances");
    for (const d of data.distances) {
      addPropRow(tbody, `${d.from} \u2192 ${d.to}`, fmtFtIn(d.distance));
    }
  } catch (_) {}
}

/** Add an editable config row with save-on-change. */
function addConfigRow(tbody, label, configKey, value) {
  const tr = document.createElement("tr");
  const tdLabel = document.createElement("td");
  tdLabel.textContent = label;
  tr.appendChild(tdLabel);
  const tdVal = document.createElement("td");
  const inp = document.createElement("input");
  inp.type = "text";
  inp.value = value;
  inp.addEventListener("change", async () => {
    await apiFetch(`/api/config/${configKey}`, {
      method: "PUT",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({value: inp.value}),
    });
    showToast(`${label} updated`, "success");
  });
  tdVal.appendChild(inp);
  tr.appendChild(tdVal);
  tbody.appendChild(tr);
}


/* ========== CANVAS RENDERING ========== */

function renderCanvas() {
  const g = App.state.geometry;
  if (!g) return;

  clearLayers();
  const overrides = itemOverrides();
  renderOutline(g);
  renderInnerWalls(g);
  renderInteriorWalls(g);
  renderOpenings(g);
  renderDoors(g, overrides);
  renderFurniture(g, overrides);
  renderClearanceZones(g, overrides);
  renderRoomLabels(g);
  renderPoints(g);
  renderUserDimensions(g);
  renderUserLabels(g);
  renderLockIcons(g);

  if (App.state.zoom === 1 && App.state.pan.x === 0 && App.state.pan.y === 0) {
    fitToWindow();
  } else {
    applyTransform();
  }
}

function clearLayers() {
  const layers = ["layer-outline", "layer-inner", "layer-walls",
    "layer-openings", "layer-doors", "layer-furniture", "layer-clearance",
    "layer-rooms", "layer-points",
    "layer-labels", "layer-dims",
    "layer-plumbing-pipes", "layer-plumbing-fittings", "layer-plumbing-fixtures",
    "layer-measure"];
  for (const id of layers) {
    App.els[id].innerHTML = "";
  }
  // Reset ghost opacity
  App.els["layer-outline"].style.opacity = "";
  App.els["layer-inner"].style.opacity = "";
  App.els["layer-walls"].style.opacity = "";
  App.els["layer-furniture"].style.opacity = "";
}

function polyToStr(poly) {
  return poly.map(p => `${p[0]},${-p[1]}`).join(" ");
}

/** Compute bounding box {w, e, s, n} from a polygon. */
function bboxFromPoly(poly) {
  let w = Infinity, e = -Infinity, s = Infinity, n = -Infinity;
  for (const p of poly) {
    if (p[0] < w) w = p[0];
    if (p[0] > e) e = p[0];
    if (p[1] < s) s = p[1];
    if (p[1] > n) n = p[1];
  }
  return { w, e, s, n };
}

/** Build override offset lookup from App.state.elements. */
function itemOverrides() {
  const ov = {};
  for (const e of (App.state.elements || [])) {
    if (e.type === "furniture" || e.type === "appliance" || e.type === "fixture") {
      let props = e.properties;
      if (typeof props === "string") props = JSON.parse(props);
      if (props && props.source === "override") {
        ov[e.name] = props;
      }
    }
  }
  return ov;
}

/** Get (ox, oy) offset for a named item from an overrides map. */
function itemOffset(overrides, name) {
  const ov = overrides[name];
  return [ov ? (ov.offset_x || 0) : 0, ov ? (ov.offset_y || 0) : 0];
}

function fmtFtIn(ft) {
  const totalIn = Math.round(Math.abs(ft) * 1200) / 100;
  const wholeFt = Math.floor(totalIn / 12);
  const remainIn = totalIn - wholeFt * 12;
  const inStr = remainIn.toFixed(2).replace(/0+$/, '').replace(/\.$/, '');
  const sign = ft < 0 ? '-' : '';
  return `${sign}${wholeFt}' ${inStr}"`;
}

function fmtDeg(deg) {
  const s = deg.toFixed(4).replace(/0+$/, '').replace(/\.$/, '');
  return `${s}°`;
}

function renderOutline(g) {
  const layer = App.els["layer-outline"];
  if (g.outline_poly && g.outline_poly.length > 0) {
    const fill = svgEl("polygon", {
      points: polyToStr(g.outline_poly),
      class: "outline-fill",
    });
    layer.appendChild(fill);

    const stroke = svgEl("polygon", {
      points: polyToStr(g.outline_poly),
      class: "outline-stroke",
    });
    layer.appendChild(stroke);
  }
}

function renderInnerWalls(g) {
  if (!App.state.showPoints) return; // inner line visible when points shown
  const layer = App.els["layer-inner"];
  if (g.inner_poly && g.inner_poly.length > 0) {
    const el = svgEl("polygon", {
      points: polyToStr(g.inner_poly),
      class: "inner-stroke",
    });
    layer.appendChild(el);
  }
}

function renderInteriorWalls(g) {
  const layer = App.els["layer-walls"];
  for (const [name, wall] of Object.entries(g.interior_walls || {})) {
    const el = svgEl("polygon", {
      points: polyToStr(wall.poly),
      class: "wall-fill selectable",
      "data-type": "wall",
      "data-name": name,
    });
    el.addEventListener("click", (e) => selectElement("wall", name, wall, e));
    layer.appendChild(el);
  }
  // Render custom drawn walls from elements
  for (const elem of (App.state.elements || [])) {
    const props = typeof elem.properties === "string"
      ? JSON.parse(elem.properties) : elem.properties;
    if (elem.type === "wall" && props && props.source === "drawn" && props.poly) {
      const el = svgEl("polygon", {
        points: polyToStr(props.poly),
        class: "wall-fill selectable",
        "data-type": "wall",
        "data-name": elem.name,
      });
      el.addEventListener("click", (e) => selectElement("wall", elem.name, { poly: props.poly, ...props }, e));
      layer.appendChild(el);
    }
  }
}

function renderOpenings(g) {
  if (!App.state.showOpenings) return;
  const layer = App.els["layer-openings"];
  for (const op of (g.outer_openings || [])) {
    const el = svgEl("polygon", {
      points: polyToStr(op.poly),
      class: "opening-fill selectable",
      "data-type": "opening",
      "data-name": op.name,
    });
    el.addEventListener("click", (e) => selectElement("opening", op.name, op, e));
    layer.appendChild(el);
  }
  for (const ro of (g.rough_openings || [])) {
    if (ro.poly) {
      const el = svgEl("polygon", {
        points: polyToStr(ro.poly),
        class: "opening-fill selectable",
        "data-type": "rough_opening",
        "data-name": ro.name,
      });
      el.addEventListener("click", (e) => selectElement("rough_opening", ro.name, ro, e));
      layer.appendChild(el);
    } else {
      const b = ro.bbox;
      const poly = [[b.w, b.s], [b.e, b.s], [b.e, b.n], [b.w, b.n]];
      const el = svgEl("polygon", {
        points: polyToStr(poly),
        class: "opening-fill selectable",
        "data-type": "rough_opening",
        "data-name": ro.name,
      });
      el.addEventListener("click", (e) => selectElement("rough_opening", ro.name, ro, e));
      layer.appendChild(el);
    }
  }
}

/** Render a single appliance door arc (line + polyline) into a layer. */
function renderApplDoor(layer, ad, ox, oy) {
  layer.appendChild(svgEl("line", {
    x1: ad.hinge[0] + ox, y1: -(ad.hinge[1] + oy),
    x2: ad.tip[0] + ox, y2: -(ad.tip[1] + oy),
    class: "appl-door-line",
  }));
  const pts = ad.arc_pts.map(p => `${p[0] + ox},${-(p[1] + oy)}`).join(" ");
  layer.appendChild(svgEl("polyline", {
    points: pts,
    class: "appl-door-arc",
  }));
}

function renderDoors(g, overrides) {
  if (!App.state.showDoors) return;
  const layer = App.els["layer-doors"];
  // Structural door arcs (openings)
  for (const door of (g.door_arcs || [])) {
    for (const leaf of door.leaves) {
      layer.appendChild(svgEl("line", {
        x1: leaf.hinge[0], y1: -leaf.hinge[1],
        x2: leaf.tip[0], y2: -leaf.tip[1],
        class: "door-line",
      }));
      const pts = leaf.arc_pts.map(p => `${p[0]},${-p[1]}`).join(" ");
      layer.appendChild(svgEl("polyline", {
        points: pts,
        class: "door-arc",
      }));
      layer.appendChild(svgEl("circle", {
        cx: leaf.hinge[0], cy: -leaf.hinge[1], r: 0.04,
        class: "door-hinge",
      }));
    }
  }
  // Appliance door arcs (fridge, washer, dryer, microwave)
  // Stacked appliance doors are rendered in renderFurniture() so they
  // appear above the counter they sit on (layer-furniture paints after
  // layer-doors in SVG paint order).
  // Appliance doors only show when both Doors and Furniture toggles are on.
  if (App.state.showFurniture) {
    for (const ad of (g.appliance_doors || [])) {
      if (ad.stacked) continue;
      const [ox, oy] = itemOffset(overrides, ad.item_name);
      renderApplDoor(layer, ad, ox, oy);
    }
  }
}

function renderClearanceZones(g, overrides) {
  if (!App.state.showClearance) return;
  const layer = App.els["layer-clearance"];
  for (const cz of (g.clearance_zones || [])) {
    // Extract parent item name from "itemname_clearance"
    const parentName = cz.name.replace(/_clearance$/, "");
    const [ox, oy] = itemOffset(overrides, parentName);
    let poly = cz.poly;
    if (ox !== 0 || oy !== 0) {
      poly = poly.map(p => [p[0] + ox, p[1] + oy]);
    }
    layer.appendChild(svgEl("polygon", {
      points: polyToStr(poly),
      class: "clearance-zone",
    }));
  }
}

function renderFurniture(g, overrides) {
  if (!App.state.showFurniture) return;
  const layer = App.els["layer-furniture"];

  // Render variant items (comprehensive set) if available; empty dict = no items
  if (g.variant_items !== undefined) {

    // Sort: non-stacked items first, stacked items on top (SVG paint order)
    const entries = Object.entries(g.variant_items);
    entries.sort((a, b) => (a[1].stacked ? 1 : 0) - (b[1].stacked ? 1 : 0));
    for (const [name, item] of entries) {
      const cssClass = `item-${item.type} selectable` + (item.stacked ? " item-stacked" : "");
      const [ox, oy] = itemOffset(overrides, name);

      const itemStyle = resolveItemStyle(overrides, name, App.state.variant, item);

      if (item.shape === "circle") {
        const c = item.center;
        const el = svgEl("circle", {
          cx: c[0] + ox, cy: -(c[1] + oy), r: item.radius,
          class: cssClass,
          "data-type": item.type,
          "data-name": name,
        });
        applyElementStyle(el, itemStyle);
        el.addEventListener("click", (e) => selectElement(item.type, name, item, e));
        // Wrap in link if product URL exists and links enabled
        if (App.state.openLinks && itemStyle && itemStyle.product_url && /^https?:\/\//.test(itemStyle.product_url)) {
          const a = document.createElementNS("http://www.w3.org/2000/svg", "a");
          a.setAttributeNS("http://www.w3.org/1999/xlink", "xlink:href", itemStyle.product_url);
          a.setAttribute("target", "_blank");
          a.appendChild(el);
          layer.appendChild(a);
        } else {
          layer.appendChild(el);
        }
      } else {
        // Shift polygon points by override offset
        let poly = item.poly;
        if (ox !== 0 || oy !== 0) {
          poly = poly.map(p => [p[0] + ox, p[1] + oy]);
        }
        const el = svgEl("polygon", {
          points: polyToStr(poly),
          class: cssClass,
          "data-type": item.type,
          "data-name": name,
        });
        applyElementStyle(el, itemStyle);
        el.addEventListener("click", (e) => selectElement(item.type, name, item, e));
        if (App.state.openLinks && itemStyle && itemStyle.product_url && /^https?:\/\//.test(itemStyle.product_url)) {
          const a = document.createElementNS("http://www.w3.org/2000/svg", "a");
          a.setAttributeNS("http://www.w3.org/1999/xlink", "xlink:href", itemStyle.product_url);
          a.setAttribute("target", "_blank");
          a.appendChild(el);
          layer.appendChild(a);
        } else {
          layer.appendChild(el);
        }
      }

      // Link icon overlay (CV-12) — anchor to actual shape, not bbox
      if (App.state.openLinks && itemStyle && itemStyle.product_url && /^https?:\/\//.test(itemStyle.product_url)) {
        let iconX, iconY;
        if (item.shape === "circle") {
          const c = item.center;
          // Top-right of circle at 45 degrees
          const r = item.radius || 0;
          iconX = c[0] + ox + r * 0.707;
          iconY = -(c[1] + oy + r * 0.707);
        } else if (item.poly) {
          // Find the topmost-rightmost vertex of the actual polygon
          let best = null;
          const pts = item.poly;
          for (const p of pts) {
            const px = p[0] + ox, py = p[1] + oy;
            // Prefer high N (large py) then high E (large px)
            if (!best || py > best[1] || (py === best[1] && px > best[0]))
              best = [px, py];
          }
          iconX = best[0];
          iconY = -best[1];
        } else if (item.bbox) {
          const bx = item.bbox;
          iconX = bx.e + ox - 0.1;
          iconY = -(bx.n + oy) + 0.15;
        } else {
          iconX = null;
        }
        if (iconX == null) continue;
        const icon = svgEl("text", {
          x: iconX, y: iconY,
          class: "link-icon",
        });
        icon.textContent = "\u{1F517}";
        icon.addEventListener("click", (e) => {
          e.stopPropagation();
          window.open(itemStyle.product_url, "_blank");
        });
        layer.appendChild(icon);
      }

      // Add text label at bbox center (shifted by override)
      if (item.label && item.bbox) {
        const bx = item.bbox;
        const cx = (bx.w + bx.e) / 2 + ox;
        const cy = -((bx.s + bx.n) / 2 + oy);
        const lbl = svgEl("text", {
          x: cx, y: cy,
          class: "item-label",
          "text-anchor": "middle",
          "pointer-events": "none",
        });
        lbl.textContent = item.label;
        layer.appendChild(lbl);
      }
    }

    // Stacked appliance door arcs: rendered here (above all items) so they
    // are not hidden by counter polygons in this same layer.
    if (App.state.showDoors) {
      for (const ad of (g.appliance_doors || [])) {
        if (!ad.stacked) continue;
        const [ox, oy] = itemOffset(overrides, ad.item_name);
        renderApplDoor(layer, ad, ox, oy);
      }
    }
  }

  // Placed items now go through the formula evaluator and appear in
  // variant_items above — no fallback rendering loop needed.
}

function renderRoomLabels(g) {
  if (!App.state.showRooms) return;
  const layer = App.els["layer-rooms"];
  if (!g.room_labels) return;

  // SF dashed partition lines (render first, behind labels)
  if (g.sf_lines) {
    for (const line of g.sf_lines) {
      const el = svgEl("line", {
        x1: line.start[0], y1: -line.start[1],
        x2: line.end[0], y2: -line.end[1],
        class: "sf-partition",
      });
      layer.appendChild(el);
    }
  }

  for (const lbl of g.room_labels) {
    const [e, n] = lbl.pos;
    const hasArea = lbl.area !== undefined;
    const showArea = hasArea && (App.state.variant === "sf" || App.state.showAreas);

    // Clickable highlight polygon (SF variant)
    if (lbl.poly) {
      const hlPoly = svgEl("polygon", {
        points: polyToStr(lbl.poly),
        class: "room-highlight",
        "data-room": lbl.name,
      });
      hlPoly.addEventListener("click", () => {
        // Toggle: remove if already highlighted, else clear others and show
        if (hlPoly.classList.contains("room-highlight-active")) {
          hlPoly.classList.remove("room-highlight-active");
        } else {
          layer.querySelectorAll(".room-highlight-active").forEach(
            el => el.classList.remove("room-highlight-active"));
          hlPoly.classList.add("room-highlight-active");
        }
      });
      layer.appendChild(hlPoly);
    }

    // Look up font_size from label_elements
    const leRec = (g.label_elements || []).find(le => le.name === lbl.name);
    const fontSize = leRec && leRec.properties.font_size ? leRec.properties.font_size : null;

    // Group for name + area text (clickable when SF)
    const nameEl = svgEl("text", {
      x: e, y: -n + (showArea ? -0.15 : 0),
      class: "room-label selectable" + (showArea ? " room-label-sf" : ""),
      "text-anchor": "middle",
      "dominant-baseline": "middle",
      "data-type": "label", "data-name": lbl.name,
    });
    if (fontSize) nameEl.style.fontSize = fontSize + "px";
    nameEl.textContent = lbl.name;
    if (lbl.poly) {
      nameEl.style.cursor = "pointer";
      nameEl.addEventListener("click", () => {
        const hl = layer.querySelector(`.room-highlight[data-room="${lbl.name}"]`);
        if (hl) hl.dispatchEvent(new Event("click"));
      });
    }
    layer.appendChild(nameEl);

    if (showArea) {
      const areaEl = svgEl("text", {
        x: e, y: -n + 0.15,
        class: "room-label room-area",
        "text-anchor": "middle",
        "dominant-baseline": "middle",
      });
      areaEl.textContent = lbl.area + " sf";
      if (lbl.poly) {
        areaEl.style.cursor = "pointer";
        areaEl.addEventListener("click", () => {
          const hl = layer.querySelector(`.room-highlight[data-room="${lbl.name}"]`);
          if (hl) hl.dispatchEvent(new Event("click"));
        });
      }
      layer.appendChild(areaEl);
    }
  }
}

function renderPoints(g) {
  if (!App.state.showPoints) return;
  const pointLayer = App.els["layer-points"];
  const labelLayer = App.els["layer-labels"];

  for (const [name, pt] of Object.entries(g.points || {})) {
    // Filter to F/W/C series and survey points
    const isFSeries = name.startsWith("F") && !name.startsWith("FC");
    const isWSeries = name.startsWith("W");
    const isCSeries = name.startsWith("C") && name !== "CTR";
    if (!isFSeries && !isWSeries && !isCSeries) continue;

    const r = isFSeries ? 0.12 : 0.08;
    const cls = isWSeries ? "point-marker inner" : "point-marker";
    const circle = svgEl("circle", {
      cx: pt[0], cy: -pt[1], r: r, class: cls,
      "data-type": "point", "data-name": name,
    });
    circle.addEventListener("click", (e) => selectElement("point", name, { pos: pt }, e));
    pointLayer.appendChild(circle);

    if (App.state.showLabels) {
      let labelClass = "point-label";
      if (isFSeries) labelClass += " f-series";
      else if (isWSeries) labelClass += " w-series";
      else if (isCSeries) labelClass += " c-series";

      const label = svgEl("text", {
        x: pt[0] + 0.2, y: -pt[1] - 0.15,
        class: labelClass,
      });
      label.textContent = name;
      labelLayer.appendChild(label);
    }
  }
}


/* ========== Phase 12f: Lock Icons & Formula Dependency Highlighting ========== */

function renderLockIcons(g) {
  const locked = g.locked_elements || [];
  if (locked.length === 0) return;
  const layer = App.els["layer-labels"];
  const lockedSet = new Set(locked.map(n => n.toLowerCase()));
  // Place lock icon at bbox center of locked elements
  const sources = [
    [g.interior_walls, true],
    [g.variant_items, false],
    [g.furniture, false],
    [g.appliances, false],
  ];
  for (const [dict, upperKey] of sources) {
    if (!dict) continue;
    for (const [name, item] of Object.entries(dict)) {
      const key = upperKey ? name.toUpperCase() : name.toLowerCase();
      if (!lockedSet.has(key.toLowerCase())) continue;
      const b = item.bbox;
      if (!b) continue;
      const cx = (b.w + b.e) / 2;
      const cy = -((b.s + b.n) / 2);
      const icon = svgEl("text", {
        x: cx, y: cy - 0.4,
        class: "lock-icon",
      });
      icon.textContent = "\u{1F512}";
      layer.appendChild(icon);
    }
  }
}

async function highlightFormulaDeps(elemName) {
  clearFormulaDeps();
  try {
    const [depsResp, deptsResp] = await Promise.all([
      fetch(`/api/formulas/${encodeURIComponent(elemName)}/deps`),
      fetch(`/api/formulas/${encodeURIComponent(elemName)}/dependents`),
    ]);
    const upstream = new Set();
    const downstream = new Set();
    if (depsResp.ok) {
      const deps = await depsResp.json();
      for (const d of deps) {
        if (d.dep_type === "element") upstream.add(d.dep_name);
      }
    }
    if (deptsResp.ok) {
      const depts = await deptsResp.json();
      for (const d of depts) downstream.add(d.element_name);
    }
    if (upstream.size === 0 && downstream.size === 0) return;
    document.querySelectorAll("[data-name]").forEach(el => {
      const name = el.getAttribute("data-name");
      if (!name) return;
      if (upstream.has(name) || upstream.has(name.toUpperCase()) || upstream.has(name.toLowerCase())) {
        el.classList.add("formula-dep-upstream");
      } else if (downstream.has(name) || downstream.has(name.toUpperCase()) || downstream.has(name.toLowerCase())) {
        el.classList.add("formula-dep-downstream");
      }
    });
  } catch (e) {
    // Non-critical — skip silently
  }
}

function clearFormulaDeps() {
  document.querySelectorAll(".formula-dep-upstream, .formula-dep-downstream").forEach(el => {
    el.classList.remove("formula-dep-upstream", "formula-dep-downstream");
  });
}


/* ========== PLUMBING DATA ========== */

async function loadPlumbingElements() {
  try {
    const resp = await apiFetch("/api/plumbing");
    App.state.plumbingElements = await resp.json();
    if (App.state.activeView === "interactive" && App.state.variant === "plumbing") renderPlumbingCanvas();
    updatePlumbingTable();
  } catch (e) {
    console.error("Plumbing load failed:", e);
  }
}

function updatePlumbingTable() {
  const elems = App.state.plumbingElements || [];
  // Fixtures table
  const ftBody = App.els["plumbing-fixtures-table"]
    ? App.els["plumbing-fixtures-table"].querySelector("tbody") : null;
  if (ftBody) {
    ftBody.innerHTML = "";
    for (const e of elems.filter(x => x.type === "fixture_connection")) {
      const p = e.properties || {};
      const tr = document.createElement("tr");
      tr.className = "selectable";
      if (p.show_in_table === false) tr.style.opacity = "0.4";
      const label = p.table_label || e.name;
      const fmtConn = v => { const n = (v === true ? 1 : v === false ? 0 : Number(v) || 0); return n > 1 ? "\u2713".repeat(n) : n === 1 ? "\u2713" : ""; };
      tr.innerHTML = `<td>${label}</td><td>${fmtConn(p.cold)}</td><td>${fmtConn(p.hot)}</td><td>${fmtConn(p.drain)}</td>`;
      tr.addEventListener("click", () => selectElement("plumbing", e.name, e));
      ftBody.appendChild(tr);
    }
  }
  // Pipes & fittings table
  const ptBody = App.els["plumbing-pipes-table"]
    ? App.els["plumbing-pipes-table"].querySelector("tbody") : null;
  if (ptBody) {
    ptBody.innerHTML = "";
    for (const e of elems.filter(x => x.type !== "fixture_connection")) {
      const tr = document.createElement("tr");
      tr.className = "selectable";
      tr.innerHTML = `<td>${e.name}</td><td>${e.type}</td>`;
      tr.addEventListener("click", () => selectElement("plumbing", e.name, e));
      ptBody.appendChild(tr);
    }
  }
}


/* ========== PLUMBING CANVAS ========== */

function renderPlumbingCanvas() {
  const g = App.state.geometry;
  if (!g) return;

  clearLayers();

  // Render building outline as ghost context (reduced opacity)
  renderOutline(g);
  renderInnerWalls(g);
  renderInteriorWalls(g);
  App.els["layer-outline"].style.opacity = "0.25";
  App.els["layer-inner"].style.opacity = "0.25";
  App.els["layer-walls"].style.opacity = "0.25";

  // Render plumbing-relevant fixtures/appliances (ghosted, matching SVG)
  const savedShowFurn = App.state.showFurniture;
  App.state.showFurniture = true;
  renderFurniture(g, {});
  App.state.showFurniture = savedShowFurn;
  App.els["layer-furniture"].style.opacity = "0.6";

  // Render plumbing elements (on top of ghosted fixtures)
  renderPlumbingPipes();
  renderPlumbingFittings();
  renderPlumbingFixtures();

  // Preserve plumbing draw preview if active
  if (PlumbingDraw.points.length > 0) {
    renderPlumbingDrawPreview();
  }

  if (App.state.zoom === 1 && App.state.pan.x === 0 && App.state.pan.y === 0) {
    fitToWindow();
  } else {
    applyTransform();
  }
}

const PLUMBING_COLORS = {
  supply_cold: "#1E90FF",
  supply_hot: "#FF4444",
  drain_pipe: "#C4956A",
};

function renderPlumbingPipes() {
  const layer = App.els["layer-plumbing-pipes"];
  const elems = App.state.plumbingElements || [];
  for (const e of elems) {
    if (e.type !== "supply_pipe" && e.type !== "drain_pipe") continue;
    const path = e.path || [];
    if (path.length < 2) continue;
    const p = e.properties || {};
    let color;
    if (e.type === "supply_pipe") {
      color = (p.hot_cold === "hot") ? PLUMBING_COLORS.supply_hot : PLUMBING_COLORS.supply_cold;
    } else {
      color = PLUMBING_COLORS.drain_pipe;
    }
    const pts = path.map(pt => `${pt[0]},${-pt[1]}`).join(" ");
    const attrs = {
      points: pts,
      class: "plumbing-pipe selectable",
      "data-type": "plumbing",
      "data-name": e.name,
      fill: "none",
      stroke: color,
      "stroke-width": "0.06",
      "stroke-linecap": "round",
      "stroke-linejoin": "round",
    };
    if (p.buried) attrs["stroke-dasharray"] = "0.15 0.08";
    const el = svgEl("polyline", attrs);
    el.addEventListener("click", (ev) => selectElement("plumbing", e.name, e, ev));
    layer.appendChild(el);

    // Slope annotation for drain pipes
    if (e.type === "drain_pipe" && p.slope && path.length >= 2) {
      const mid = [
        (path[0][0] + path[path.length - 1][0]) / 2,
        (path[0][1] + path[path.length - 1][1]) / 2,
      ];
      const label = svgEl("text", {
        x: mid[0], y: -mid[1] - 0.15,
        class: "plumbing-slope-label",
      });
      label.textContent = p.slope;
      layer.appendChild(label);
    }
  }
}

function renderPlumbingFittings() {
  const layer = App.els["layer-plumbing-fittings"];
  const elems = App.state.plumbingElements || [];
  for (const e of elems) {
    if (e.type !== "fitting") continue;
    const path = e.path || [];
    if (path.length < 1) continue;
    const pt = path[0];
    const p = e.properties || {};
    const ft = p.fitting_type || "tee";
    const rot = p.rotation || 0;
    const g = svgEl("g", {
      transform: `translate(${pt[0]},${-pt[1]}) rotate(${-rot})`,
      class: "plumbing-fitting selectable",
      "data-type": "plumbing",
      "data-name": e.name,
    });
    // Draw fitting symbol
    if (ft === "elbow90" || ft === "elbow45") {
      g.appendChild(svgEl("path", {
        d: ft === "elbow90"
          ? "M -0.15,0 L 0,0 L 0,-0.15"
          : "M -0.15,0 L 0,0 L -0.05,-0.14",
        fill: "none", stroke: "#ccc", "stroke-width": "0.04",
      }));
    } else if (ft === "valve") {
      g.appendChild(svgEl("path", {
        d: "M -0.1,-0.08 L 0.1,0.08 M -0.1,0.08 L 0.1,-0.08",
        fill: "none", stroke: "#f0c040", "stroke-width": "0.04",
      }));
    } else {
      // Default: tee
      g.appendChild(svgEl("line", {
        x1: -0.15, y1: 0, x2: 0.15, y2: 0,
        stroke: "#ccc", "stroke-width": "0.04",
      }));
      g.appendChild(svgEl("line", {
        x1: 0, y1: 0, x2: 0, y2: -0.15,
        stroke: "#ccc", "stroke-width": "0.04",
      }));
    }
    g.appendChild(svgEl("circle", {
      cx: 0, cy: 0, r: 0.06, fill: "#555", stroke: "#ccc", "stroke-width": "0.02",
    }));
    g.addEventListener("click", (ev) => selectElement("plumbing", e.name, e, ev));
    layer.appendChild(g);
  }
}

function renderPlumbingFixtures() {
  const layer = App.els["layer-plumbing-fixtures"];
  const elems = App.state.plumbingElements || [];
  for (const e of elems) {
    if (e.type !== "fixture_connection") continue;
    const path = e.path || [];
    if (path.length < 1) continue;
    const pt = path[0];
    const p = e.properties || {};
    // Determine marker color based on connections
    let color = "#888";
    if (p.hot && p.cold) color = "#b060d0";
    else if (p.hot) color = PLUMBING_COLORS.supply_hot;
    else if (p.cold) color = PLUMBING_COLORS.supply_cold;
    const g = svgEl("g", {
      transform: `translate(${pt[0]},${-pt[1]})`,
      class: "plumbing-fixture-marker selectable",
      "data-type": "plumbing",
      "data-name": e.name,
    });
    g.appendChild(svgEl("circle", {
      cx: 0, cy: 0, r: 0.15, fill: "none", stroke: color, "stroke-width": "0.03",
    }));
    g.appendChild(svgEl("circle", {
      cx: 0, cy: 0, r: 0.05, fill: color,
    }));
    // Label
    const label = svgEl("text", {
      x: 0.2, y: 0.05, class: "plumbing-fixture-label",
    });
    label.textContent = e.name;
    g.appendChild(label);
    g.addEventListener("click", (ev) => selectElement("plumbing", e.name, e, ev));
    layer.appendChild(g);
  }
}


/* ========== PLUMBING DRAWING TOOLS ========== */

const PlumbingDraw = {
  points: [],
  type: null,       // "supply_pipe" or "drain_pipe"
  hotCold: null,    // "cold" or "hot"
};

function isPlumbingDrawTool(tool) {
  return tool === "supply-cold" || tool === "supply-hot" || tool === "drain";
}

function plumbingDrawClick(wx, wy) {
  PlumbingDraw.points.push([Math.round(wx * 100) / 100, Math.round(wy * 100) / 100]);
  renderPlumbingDrawPreview();
}

function plumbingDrawDblClick(wx, wy) {
  // Double-click finishes path (the point was already added by click)
  finishPlumbingDraw();
}

async function finishPlumbingDraw() {
  if (PlumbingDraw.points.length < 2) {
    cancelPlumbingDraw();
    return;
  }
  const props = {};
  if (PlumbingDraw.type === "supply_pipe") {
    props.hot_cold = PlumbingDraw.hotCold;
    props.pipe_size = "0.75";
  } else {
    props.slope = "0.25 in/ft";
  }
  const name = `${PlumbingDraw.type}_${Date.now()}`;
  try {
    await apiFetch("/api/plumbing", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        type: PlumbingDraw.type,
        name: name,
        path: PlumbingDraw.points,
        properties: props,
      }),
    });
    showToast(`${PlumbingDraw.type} created`, "success");
  } catch (e) {
    showToast(`Failed: ${e.message}`, "error");
  }
  PlumbingDraw.points = [];
  PlumbingDraw.type = null;
  PlumbingDraw.hotCold = null;
}

function cancelPlumbingDraw() {
  PlumbingDraw.points = [];
  PlumbingDraw.type = null;
  PlumbingDraw.hotCold = null;
  if (App.state.activeView === "interactive" && App.state.variant === "plumbing") renderPlumbingCanvas();
}

function renderPlumbingDrawPreview() {
  const layer = App.els["layer-plumbing-pipes"];
  // Remove existing preview
  const old = layer.querySelector(".plumbing-draw-preview");
  if (old) old.remove();
  if (PlumbingDraw.points.length < 1) return;
  let color = PLUMBING_COLORS.drain_pipe;
  if (PlumbingDraw.type === "supply_pipe") {
    color = PlumbingDraw.hotCold === "hot" ? PLUMBING_COLORS.supply_hot : PLUMBING_COLORS.supply_cold;
  }
  const pts = PlumbingDraw.points.map(p => `${p[0]},${-p[1]}`).join(" ");
  const el = svgEl("polyline", {
    points: pts,
    class: "plumbing-draw-preview",
    fill: "none",
    stroke: color,
    "stroke-width": "0.06",
    "stroke-dasharray": "0.15 0.08",
    "stroke-linecap": "round",
    "stroke-linejoin": "round",
    opacity: "0.7",
  });
  layer.appendChild(el);
  // Vertex markers
  for (const p of PlumbingDraw.points) {
    layer.appendChild(svgEl("circle", {
      cx: p[0], cy: -p[1], r: 0.08,
      fill: color, opacity: "0.5",
      class: "plumbing-draw-preview",
    }));
  }
}

async function placeFixture(wx, wy) {
  const pt = [Math.round(wx * 100) / 100, Math.round(wy * 100) / 100];
  const pendingId = App.state._pendingFixturePlacement;
  try {
    if (pendingId) {
      // Update existing fixture's position
      await apiFetch(`/api/plumbing/${pendingId}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ path: [pt] }),
      });
      delete App.state._pendingFixturePlacement;
      showToast("Fixture placed", "success");
    } else {
      // Create new fixture at click location
      const elems = App.state.plumbingElements || [];
      const n = elems.filter(e => e.type === "fixture_connection").length + 1;
      const name = `Fixture ${n}`;
      await apiFetch("/api/plumbing", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          type: "fixture_connection", name,
          path: [pt],
          properties: { cold: 1, hot: 0, drain: 0 },
        }),
      });
      showToast("Fixture placed", "success");
    }
  } catch (e) {
    showToast(`Failed: ${e.message}`, "error");
  }
}

async function addFixtureConnection() {
  const elems = App.state.plumbingElements || [];
  const n = elems.filter(e => e.type === "fixture_connection").length + 1;
  const name = prompt("Fixture name:", `Fixture ${n}`);
  if (!name) return;
  try {
    const resp = await apiFetch("/api/plumbing", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        type: "fixture_connection", name,
        path: [],
        properties: { cold: 0, hot: 0, drain: 0 },
      }),
    });
    const created = await resp.json();
    await loadPlumbingElements();
    // Select the new fixture and prompt placement
    selectElement("plumbing", name, created);
    showToast("Fixture added — click canvas to place it", "success");
    // Activate place-fixture tool with pending update for this fixture
    App.state._pendingFixturePlacement = created.id;
    setTool("place-fixture");
  } catch (e) {
    showToast(`Failed: ${e.message}`, "error");
  }
}

async function placeFitting(wx, wy) {
  const name = `fitting_${Date.now()}`;
  try {
    await apiFetch("/api/plumbing", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        type: "fitting",
        name: name,
        path: [[Math.round(wx * 100) / 100, Math.round(wy * 100) / 100]],
        properties: { fitting_type: "tee", rotation: 0, size: "0.75" },
      }),
    });
    showToast("Fitting placed", "success");
  } catch (e) {
    showToast(`Failed: ${e.message}`, "error");
  }
}


/* ========== SVG HELPERS ========== */

function svgEl(tag, attrs) {
  const el = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) {
    el.setAttribute(k, v);
  }
  return el;
}


/* ========== PAN / ZOOM ========== */

function applyTransform() {
  const s = App.state;
  App.els["canvas-transform"].setAttribute("transform",
    `translate(${s.pan.x}, ${s.pan.y}) scale(${s.zoom})`
  );
  App.els["zoom-level"].textContent = `${Math.round(s.zoom * 100)}%`;

  // Update grid
  const grid = App.els["grid-rect"];
  grid.setAttribute("visibility", s.showGrid ? "visible" : "hidden");
}

function fitToWindow() {
  if (!isCanvasView()) {
    svgViewFit();
    return;
  }
  const g = App.state.geometry;
  if (!g || !g.bbox) return;

  const vp = App.els["viewport"];
  const vw = vp.clientWidth;
  const vh = vp.clientHeight;
  const margin = 40;

  const bw = g.bbox.e - g.bbox.w;
  const bh = g.bbox.n - g.bbox.s;
  const cx = (g.bbox.w + g.bbox.e) / 2;
  const cy = -(g.bbox.s + g.bbox.n) / 2; // negated for SVG Y

  const scaleX = (vw - 2 * margin) / bw;
  const scaleY = (vh - 2 * margin) / bh;
  App.state.zoom = Math.min(scaleX, scaleY);

  App.state.pan.x = vw / 2 - cx * App.state.zoom;
  App.state.pan.y = vh / 2 - cy * App.state.zoom;

  applyTransform();
}

function screenToWorld(sx, sy) {
  const s = App.state;
  const wx = (sx - s.pan.x) / s.zoom;
  const wy = -(sy - s.pan.y) / s.zoom; // negate for E/N
  return [wx, wy];
}

/** Convert a mouse event to world coordinates. */
function mouseToWorld(e) {
  const rect = App.els["viewport"].getBoundingClientRect();
  return screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
}


/* ========== SELECTION ========== */

function selectElement(type, name, data, event) {
  if (event) event.stopPropagation();

  // After a move-tool drag, suppress the click that fires on mouseup
  if (MoveTool._suppressClick) {
    MoveTool._suppressClick = false;
    return;
  }

  // Ctrl+Click or Shift+Click: toggle in/out of multi-selection
  if (event && (event.ctrlKey || event.metaKey || event.shiftKey)) {
    const idx = App.state.selections.findIndex(s => s.name === name && s.type === type);
    if (idx >= 0) {
      App.state.selections.splice(idx, 1);
    } else {
      App.state.selections.push({ type, name, data });
    }
    // Update highlight for multi-select
    document.querySelectorAll(".multi-selected").forEach(el => {
      el.classList.remove("multi-selected");
    });
    for (const s of App.state.selections) {
      const el = document.querySelector(`[data-name="${s.name}"][data-type="${s.type}"]`);
      if (el) el.classList.add("multi-selected");
    }
    // Keep last-clicked as primary selection
    App.state.selection = { type, name, data };
    App.els["selection-info"].textContent =
      App.state.selections.length > 1
        ? `${App.state.selections.length} selected`
        : `${type}: ${name}`;
    showProperties(type, name, data);
    return;
  }

  // Regular click: clear multi-selection, set single
  App.state.selections = [];
  document.querySelectorAll(".multi-selected").forEach(el => {
    el.classList.remove("multi-selected");
  });

  // Remove previous selection highlight
  document.querySelectorAll(".selected-highlight").forEach(el => {
    el.classList.remove("selected-highlight");
  });

  App.state.selection = { type, name, data };

  // Highlight selected element
  const el = document.querySelector(`[data-name="${name}"][data-type="${type}"]`);
  if (el) el.classList.add("selected-highlight");

  App.els["selection-info"].textContent = `${type}: ${name}`;
  showProperties(type, name, data);

  // Phase 12f: highlight formula dependencies
  highlightFormulaDeps(name);
}

function clearSelection() {
  document.querySelectorAll(".selected-highlight").forEach(el => {
    el.classList.remove("selected-highlight");
  });
  document.querySelectorAll(".multi-selected").forEach(el => {
    el.classList.remove("multi-selected");
  });
  clearWallHandles();
  clearConstantHighlights();
  clearFormulaDeps();
  App.state.selection = null;
  App.state.selections = [];
  App.els["selection-info"].textContent = "No selection";
  App.els["props-empty"].style.display = "block";
  App.els["props-detail"].style.display = "none";
}

function switchToPanel(panelName) {
  document.querySelectorAll(".panel-tab").forEach(t => t.classList.remove("active"));
  document.querySelectorAll(".panel-content").forEach(p => p.classList.remove("active"));
  const tab = document.querySelector(`.panel-tab[data-panel="${panelName}"]`);
  const panel = document.getElementById(`panel-${panelName}`);
  if (tab) tab.classList.add("active");
  if (panel) panel.classList.add("active");
}

function showProperties(type, name, data) {
  switchToPanel("properties");
  App.els["props-empty"].style.display = "none";
  App.els["props-detail"].style.display = "block";
  App.els["props-title"].textContent = `${type.toUpperCase()}: ${name}`;

  const tbody = App.els["props-table"].querySelector("tbody");
  tbody.innerHTML = "";

  if (type === "point") {
    addPropRow(tbody, "Easting", fmtFtIn(data.pos[0]));
    addPropRow(tbody, "Northing", fmtFtIn(data.pos[1]));
  } else if (type === "wall") {
    // Check if this is a drawn wall
    const elemRec = (App.state.elements || []).find(e => e.name === name);
    const props = elemRec ? parseProps(elemRec) : null;
    if (props && props.source === "drawn") {
      // TL-16: Drawn wall — show editable thickness
      addPropRow(tbody, "Source", "drawn");
      if (props.start) addPropRow(tbody, "Start", `${fmtFtIn(props.start[0])}, ${fmtFtIn(props.start[1])}`);
      if (props.end) addPropRow(tbody, "End", `${fmtFtIn(props.end[0])}, ${fmtFtIn(props.end[1])}`);
      // Editable thickness
      const thickTr = document.createElement("tr");
      const thickTd1 = document.createElement("td");
      thickTd1.textContent = "Thickness";
      thickTr.appendChild(thickTd1);
      const thickTd2 = document.createElement("td");
      const thickInp = document.createElement("input");
      thickInp.type = "text";
      thickInp.className = "prop-edit-input";
      thickInp.value = fmtFtIn(props.thickness || 4.0 / 12.0);
      thickInp.addEventListener("change", async () => {
        const newThick = parseDimension(thickInp.value);
        if (newThick === null || newThick <= 0) {
          showToast("Invalid thickness", "error");
          return;
        }
        const newPoly = wallPoly(props.start, props.end, newThick);
        if (!newPoly) return;
        const newProps = { ...props, thickness: newThick, poly: newPoly };
        await fetch(`/api/elements/${elemRec.id}`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ properties: newProps }),
        });
        await reloadAfterChange();
      });
      thickTd2.appendChild(thickInp);
      thickTr.appendChild(thickTd2);
      tbody.appendChild(thickTr);
      // TL-17: render endpoint drag handles
      renderWallHandles(elemRec.id, props);
    } else {
      const b = data.bbox;
      addPropRow(tbody, "Width", fmtFtIn(b.e - b.w));
      addPropRow(tbody, "Height", fmtFtIn(b.n - b.s));
      addPropRow(tbody, "West", fmtFtIn(b.w));
      addPropRow(tbody, "South", fmtFtIn(b.s));
      addPropRow(tbody, "East", fmtFtIn(b.e));
      addPropRow(tbody, "North", fmtFtIn(b.n));
    }
    // Show related constants
    const related = findRelatedConstants(name);
    if (related.length > 0) {
      addPropRow(tbody, "—", "Related Constants");
      for (const c of related) {
        addPropRow(tbody, c.name, formatConstValue(c), true, c.name);
      }
    }
    // Phase 12b: formula section
    addFormulaSection(tbody, name);
    // Layout checkboxes + delete (like furniture addElementActions)
    addWallActions(tbody, name, elemRec);
  } else if (type === "opening" || type === "rough_opening") {
    addPropRow(tbody, "Name", data.name);
    if (data.seg_start) addPropRow(tbody, "Segment", `${data.seg_start}–${data.seg_end}`);
    if (data.wall_name) addPropRow(tbody, "Wall", data.wall_name);
    // SEL-10: Editable width via controlling constant
    const widthConstName = findWidthConstant(data.name);
    const widthConst = widthConstName && App.state.constants.find(c => c.name === widthConstName);
    if (widthConst) {
      const label = widthConstName.includes("HALF") ? "Half-width" : "Width";
      addPropRow(tbody, label, formatConstValue(widthConst), true, widthConstName);
    } else if (data.width) {
      addPropRow(tbody, "Width", fmtFtIn(data.width));
    }
    if (data.orientation) addPropRow(tbody, "Orient", data.orientation);
    if (data.poly) {
      const dx = data.poly[1][0] - data.poly[0][0];
      const dy = data.poly[1][1] - data.poly[0][1];
      const w = Math.sqrt(dx * dx + dy * dy);
      addPropRow(tbody, "Actual width", fmtFtIn(w));
    }
    // SEL-7: Door properties for openings
    if ((type === "rough_opening" || type === "opening") && data.name) {
      showDoorProperties(tbody, data.name);
    }
    const related = findRelatedConstants(data.name);
    if (related.length > 0) {
      addPropRow(tbody, "—", "Related Constants");
      for (const c of related) {
        addPropRow(tbody, c.name, formatConstValue(c), true, c.name);
      }
    }
  } else if (type === "appliance" || type === "furniture" || type === "fixture") {
    // SEL-8: Enhanced furniture/appliance properties
    const b = data.bbox;
    const w = b.e - b.w;
    const d = b.n - b.s;
    addPropRow(tbody, "Type", type);
    // Editable Width/Depth for all items with DB records (unified system)
    const elemRec = (App.state.elements || []).find(e => e.name === name);
    const elemProps = elemRec ? parseProps(elemRec) : null;
    // Use geometry-computed width/depth if available, fall back to props then bbox
    const dispW = data.width || (elemProps && elemProps.width) || w;
    const dispD = data.depth || (elemProps && elemProps.depth) || d;
    if (elemRec) {
      const tr1 = document.createElement("tr");
      tr1.innerHTML = `<td>Width</td><td><input type="text" value="${fmtFtIn(dispW)}" /></td>`;
      tr1.querySelector("input").addEventListener("change", (e) => handleElementPropEdit(elemRec.id, "width", e.target.value));
      tbody.appendChild(tr1);
      const tr2 = document.createElement("tr");
      tr2.innerHTML = `<td>Depth</td><td><input type="text" value="${fmtFtIn(dispD)}" /></td>`;
      tr2.querySelector("input").addEventListener("change", (e) => handleElementPropEdit(elemRec.id, "depth", e.target.value));
      tbody.appendChild(tr2);
    } else {
      addPropRow(tbody, "Width", fmtFtIn(dispW));
      addPropRow(tbody, "Depth", fmtFtIn(dispD));
    }
    if (data.center) {
      addPropRow(tbody, "Center E", fmtFtIn(data.center[0]));
      addPropRow(tbody, "Center N", fmtFtIn(data.center[1]));
    } else {
      addPropRow(tbody, "Center E", fmtFtIn((b.w + b.e) / 2));
      addPropRow(tbody, "Center N", fmtFtIn((b.s + b.n) / 2));
    }
    if (data.rotation !== undefined) {
      if (elemRec) {
        const rotVal = Math.round(data.rotation * 10) / 10;
        const trR = document.createElement("tr");
        trR.innerHTML = `<td>Rotation</td><td><input type="text" value="${rotVal}°" /></td>`;
        trR.querySelector("input").addEventListener("change", (e) => {
          const v = parseFloat(e.target.value.replace("°", ""));
          if (!isNaN(v)) handleElementPropEdit(elemRec.id, "rotation", v);
        });
        tbody.appendChild(trR);
      } else {
        addPropRow(tbody, "Rotation", data.rotation.toFixed(1) + "°");
      }
    }
    if (elemRec) {
      const props = typeof elemRec.properties === "string"
        ? JSON.parse(elemRec.properties) : elemRec.properties;
      if (props && typeof addShapePicker === "function") {
        addShapePicker(tbody, elemRec, props);
      }
    }
    const related = findRelatedConstants(name);
    if (related.length > 0) {
      addPropRow(tbody, "—", "Related Constants");
      for (const c of related) {
        addPropRow(tbody, c.name, formatConstValue(c), true, c.name);
      }
    }
    // Phase 12b: formula section
    addFormulaSection(tbody, name);
    // Style controls and product URL
    const _elemRec = elemRec || (App.state.elements || []).find(e => e.name === name);
    if (_elemRec) {
      const _props = parseProps(_elemRec);
      // Use variant item's product_url as fallback if not in DB props
      if (!_props.product_url && data && data.product_url) _props.product_url = data.product_url;
      addStyleControls(tbody, _elemRec, _props, type);
      addViewOverrideControls(tbody, _elemRec, _props);
      addProductUrlField(tbody, _elemRec, _props);
      addElementActions(tbody, _elemRec);
    } else {
      // No DB record — show product URL (if any) + formula-aware delete
      if (data && data.product_url) {
        addProductUrlField(tbody, null, { product_url: data.product_url });
      }
      addFormulaDeleteButton(tbody, name);
    }
  } else if (type === "dimension") {
    const elemRec = (App.state.elements || []).find(e => e.name === name && e.type === "dimension");
    if (elemRec) {
      const props = parseProps(elemRec);
      addPropRow(tbody, "Source", props.source || "user");
      if (props.start) addPropRow(tbody, "Start", `${fmtFtIn(props.start[0])}, ${fmtFtIn(props.start[1])}`);
      if (props.end) addPropRow(tbody, "End", `${fmtFtIn(props.end[0])}, ${fmtFtIn(props.end[1])}`);
      const dist = props.start && props.end ?
        Math.sqrt((props.end[0] - props.start[0]) ** 2 + (props.end[1] - props.start[1]) ** 2) : 0;
      addPropRow(tbody, "Distance", fmtFtIn(dist));
      addPropRow(tbody, "Offset", fmtFtIn(props.offset || 0));
      addPropRow(tbody, "Label Rot.", props.label_rotation || "parallel");
      // Style selector (solid = builtin look, dashed = user look)
      const styleTr = document.createElement("tr");
      const styleTd1 = document.createElement("td"); styleTd1.textContent = "Style";
      styleTr.appendChild(styleTd1);
      const styleTd2 = document.createElement("td");
      const styleSel = document.createElement("select");
      styleSel.className = "prop-edit-input";
      const curStyle = props.dim_style || (props.source === "builtin" ? "solid" : "dashed");
      for (const [val, lbl] of [["solid", "Solid"], ["dashed", "Dashed"]]) {
        const opt = document.createElement("option");
        opt.value = val; opt.textContent = lbl;
        if (val === curStyle) opt.selected = true;
        styleSel.appendChild(opt);
      }
      styleSel.addEventListener("change", async () => {
        const cur = (App.state.elements || []).find(e => e.id === elemRec.id);
        const curProps = cur ? parseProps(cur) : props;
        const newProps = { ...curProps, dim_style: styleSel.value };
        await fetch(`/api/elements/${elemRec.id}`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ properties: newProps }),
        });
        await reloadAfterChange();
      });
      styleTd2.appendChild(styleSel);
      styleTr.appendChild(styleTd2);
      tbody.appendChild(styleTr);
      // Anchor info
      const fmtAnchor = (a) => {
        if (!a) return "(absolute)";
        if (a.type === "point") return `${a.target} (point)`;
        if (a.type === "wall_face") return `${a.target} ${a.face} face`;
        if (a.type === "opening_face") return `${a.target} ${a.face} face`;
        if (a.type === "line_intersection") return "line intersection";
        if (a.type === "computed") return "computed";
        return JSON.stringify(a);
      };
      addPropRow(tbody, "Start Anchor", fmtAnchor(props.start_anchor));
      addPropRow(tbody, "End Anchor", fmtAnchor(props.end_anchor));
      addStyleControls(tbody, elemRec, props, "dimension");
      addViewOverrideControls(tbody, elemRec, props);
      addElementActions(tbody, elemRec);
    }
  } else if (type === "label") {
    const elemRec = (App.state.elements || []).find(e => e.name === name && e.type === "label");
    if (elemRec) {
      const props = parseProps(elemRec);
      addPropRow(tbody, "Source", props.source || "unknown");
      addPropRow(tbody, "Text", props.text || name);
      // Editable font size
      const fsTr = document.createElement("tr");
      const fsTd1 = document.createElement("td"); fsTd1.textContent = "Font Size";
      fsTr.appendChild(fsTd1);
      const fsTd2 = document.createElement("td");
      const fsInp = document.createElement("input");
      fsInp.type = "text"; fsInp.className = "prop-edit-input";
      fsInp.value = (props.font_size || 0.25).toString();
      fsInp.addEventListener("change", async () => {
        const newFs = parseFloat(fsInp.value);
        if (isNaN(newFs) || newFs <= 0) { showToast("Invalid font size", "error"); return; }
        const newProps = { ...props, font_size: newFs };
        await fetch(`/api/elements/${elemRec.id}`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ properties: newProps }),
        });
        await reloadAfterChange();
      });
      fsTd2.appendChild(fsInp);
      fsTr.appendChild(fsTd2);
      tbody.appendChild(fsTr);
      addStyleControls(tbody, elemRec, props, "label");
      addViewOverrideControls(tbody, elemRec, props);
      addElementActions(tbody, elemRec);
    }
  } else if (type === "plumbing") {
    const pe = data;  // full plumbing element object
    if (!pe || !pe.id) return;
    const p = pe.properties || {};

    if (pe.type === "fixture_connection") {
      // ── Fixture connection properties ──
      App.els["props-title"].textContent = pe.name;

      // Editable name
      const nameTr = document.createElement("tr");
      const nameTd1 = document.createElement("td"); nameTd1.textContent = "Name";
      nameTr.appendChild(nameTd1);
      const nameTd2 = document.createElement("td");
      const nameInp = document.createElement("input");
      nameInp.type = "text"; nameInp.className = "prop-edit-input";
      nameInp.value = pe.name;
      nameInp.addEventListener("change", async () => {
        const v = nameInp.value.trim();
        if (!v) { showToast("Name cannot be empty", "error"); return; }
        await apiFetch(`/api/plumbing/${pe.id}`, {
          method: "PUT", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ name: v }),
        });
        loadPlumbingElements();
      });
      nameTd2.appendChild(nameInp);
      nameTr.appendChild(nameTd2);
      tbody.appendChild(nameTr);

      // Read-only position
      if (pe.path && pe.path.length > 0) {
        addPropRow(tbody, "E", fmtFtIn(pe.path[0][0]));
        addPropRow(tbody, "N", fmtFtIn(pe.path[0][1]));
      }

      // Service connection counts (0 = none, 1+ = number of connections)
      for (const flag of ["cold", "hot", "drain"]) {
        const tr = document.createElement("tr");
        const td1 = document.createElement("td"); td1.textContent = flag.charAt(0).toUpperCase() + flag.slice(1);
        tr.appendChild(td1);
        const td2 = document.createElement("td");
        const inp = document.createElement("input");
        inp.type = "number"; inp.min = "0"; inp.max = "9"; inp.step = "1";
        inp.style.width = "3em";
        // Normalize: true→1, false→0, number→number
        const cur = p[flag];
        inp.value = (cur === true ? 1 : cur === false ? 0 : Number(cur) || 0);
        inp.addEventListener("change", async () => {
          const v = Math.max(0, parseInt(inp.value) || 0);
          inp.value = v;
          const newProps = { ...p, [flag]: v };
          await apiFetch(`/api/plumbing/${pe.id}`, {
            method: "PUT", headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ properties: newProps }),
          });
          loadPlumbingElements();
        });
        td2.appendChild(inp);
        tr.appendChild(td2);
        tbody.appendChild(tr);
      }

      // Show in table checkbox
      {
        const tr = document.createElement("tr");
        const td1 = document.createElement("td"); td1.textContent = "In Table";
        tr.appendChild(td1);
        const td2 = document.createElement("td");
        const cb = document.createElement("input");
        cb.type = "checkbox"; cb.checked = p.show_in_table !== false;
        cb.addEventListener("change", async () => {
          const newProps = { ...p, show_in_table: cb.checked };
          await apiFetch(`/api/plumbing/${pe.id}`, {
            method: "PUT", headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ properties: newProps }),
          });
          loadPlumbingElements();
        });
        td2.appendChild(cb);
        tr.appendChild(td2);
        tbody.appendChild(tr);
      }

      // Table label override
      {
        const tr = document.createElement("tr");
        const td1 = document.createElement("td"); td1.textContent = "Table Label";
        tr.appendChild(td1);
        const td2 = document.createElement("td");
        const inp = document.createElement("input");
        inp.type = "text"; inp.className = "prop-edit-input";
        inp.value = p.table_label || "";
        inp.placeholder = pe.name.toUpperCase();
        inp.addEventListener("change", async () => {
          const v = inp.value.trim();
          const newProps = { ...p };
          if (v) newProps.table_label = v;
          else delete newProps.table_label;
          await apiFetch(`/api/plumbing/${pe.id}`, {
            method: "PUT", headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ properties: newProps }),
          });
          loadPlumbingElements();
        });
        td2.appendChild(inp);
        tr.appendChild(td2);
        tbody.appendChild(tr);
      }

      // Delete button
      const delTr = document.createElement("tr");
      const delTd = document.createElement("td"); delTd.colSpan = 2;
      const delBtn = document.createElement("button");
      delBtn.textContent = "Delete"; delBtn.className = "prop-delete-btn";
      delBtn.addEventListener("click", async () => {
        await apiFetch(`/api/plumbing/${pe.id}`, { method: "DELETE" });
        clearSelection();
        loadPlumbingElements();
      });
      delTd.appendChild(delBtn);
      delTr.appendChild(delTd);
      tbody.appendChild(delTr);

    } else {
      // ── Pipe / fitting properties ──
      App.els["props-title"].textContent = pe.name;
      addPropRow(tbody, "Type", pe.type);
      if (pe.path) addPropRow(tbody, "Points", pe.path.length.toString());
      for (const [k, v] of Object.entries(p)) {
        if (k === "buried") continue;  // handled below
        addPropRow(tbody, k, String(v));
      }
      // Always show buried checkbox for pipes
      if (pe.type && pe.type.includes("pipe")) {
        const tr = document.createElement("tr");
        const td1 = document.createElement("td"); td1.textContent = "Buried";
        tr.appendChild(td1);
        const td2 = document.createElement("td");
        const cb = document.createElement("input");
        cb.type = "checkbox"; cb.checked = !!p.buried;
        cb.addEventListener("change", async () => {
          const newProps = { ...p, buried: cb.checked };
          await apiFetch(`/api/plumbing/${pe.id}`, {
            method: "PUT", headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ properties: newProps }),
          });
          loadPlumbingElements();
        });
        td2.appendChild(cb);
        tr.appendChild(td2);
        tbody.appendChild(tr);
      }

      // Delete button
      const delTr = document.createElement("tr");
      const delTd = document.createElement("td"); delTd.colSpan = 2;
      const delBtn = document.createElement("button");
      delBtn.textContent = "Delete"; delBtn.className = "prop-delete-btn";
      delBtn.addEventListener("click", async () => {
        await apiFetch(`/api/plumbing/${pe.id}`, { method: "DELETE" });
        clearSelection();
        loadPlumbingElements();
      });
      delTd.appendChild(delBtn);
      delTr.appendChild(delTd);
      tbody.appendChild(delTr);
    }
  }
}

// VARIANT_LABELS removed — use variantLabel(name) for dynamic lookup

// ── Style defaults matching CSS classes ─────────────────────────
const STYLE_DEFAULTS = {
  appliance: { fill_color: "#2a3a4a", stroke_color: "#4682B4", stroke_width: 0.02, stroke_style: "solid", opacity: 100 },
  furniture: { fill_color: "#2a3a2a", stroke_color: "#5a8a5a", stroke_width: 0.02, stroke_style: "solid", opacity: 100 },
  fixture:   { fill_color: "#3a2a3a", stroke_color: "#8a5a8a", stroke_width: 0.02, stroke_style: "solid", opacity: 100 },
  wall:      { fill_color: "#334",    stroke_color: "#556",    stroke_width: 0.02, stroke_style: "solid", opacity: 100 },
  opening:   { fill_color: "#2a4a6a", stroke_color: "#4488cc", stroke_width: 0.02, stroke_style: "solid", opacity: 100 },
  dimension: { fill_color: null,      stroke_color: null,      stroke_width: 0.02, stroke_style: "solid", opacity: 100 },
  label:     { fill_color: null,      stroke_color: null,      stroke_width: null,  stroke_style: null,    opacity: 100 },
};

/** Apply inline style overrides to an SVG element from resolved properties. */
function applyElementStyle(el, props) {
  if (!props) return;
  if (props.fill_color) el.style.fill = props.fill_color;
  if (props.stroke_color) el.style.stroke = props.stroke_color;
  if (props.stroke_width != null) el.style.strokeWidth = props.stroke_width;
  if (props.stroke_style === "dashed") el.style.strokeDasharray = "0.06 0.03";
  else if (props.stroke_style === "dotted") el.style.strokeDasharray = "0.02 0.02";
  else if (props.stroke_style === "solid") el.style.strokeDasharray = "none";
  if (props.opacity != null && props.opacity !== 100) el.style.opacity = props.opacity / 100;
}

/** Resolve style properties with per-view override merging. */
function resolveItemStyle(overrides, name, variant, variantItem) {
  const ov = overrides[name];
  // Start with variant item's product_url (from engine) as base
  const base = {};
  if (variantItem && variantItem.product_url) base.product_url = variantItem.product_url;
  if (!ov) return Object.keys(base).length ? base : null;
  const resolved = { ...base, ...ov };
  const viewOv = (ov.view_overrides || {})[variant];
  if (viewOv) Object.assign(resolved, viewOv);
  return resolved;
}

/** Add style controls (fill, stroke, opacity) to properties panel. */
function addStyleControls(tbody, elemRec, props, elementType) {
  if (!elemRec) return;
  const defaults = STYLE_DEFAULTS[elementType] || STYLE_DEFAULTS.furniture;

  // Helper to save a style property
  async function saveStyleProp(key, value) {
    const cur = (App.state.elements || []).find(e => e.id === elemRec.id);
    const curProps = cur ? parseProps(cur) : props;
    const newProps = { ...curProps, [key]: value };
    await fetch(`/api/elements/${elemRec.id}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ properties: newProps }),
    });
    await reloadAfterChange();
  }

  // Fill colour (not for dimension/label)
  if (defaults.fill_color !== null) {
    const tr = document.createElement("tr");
    const td1 = document.createElement("td"); td1.textContent = "Fill";
    tr.appendChild(td1);
    const td2 = document.createElement("td");
    const inp = document.createElement("input");
    inp.type = "color";
    inp.className = "prop-color-input";
    inp.value = props.fill_color || defaults.fill_color || "#333333";
    inp.addEventListener("change", () => saveStyleProp("fill_color", inp.value));
    td2.appendChild(inp);
    if (props.fill_color) {
      const rst = document.createElement("button");
      rst.textContent = "Reset";
      rst.className = "prop-reset-btn";
      rst.addEventListener("click", () => saveStyleProp("fill_color", null));
      td2.appendChild(rst);
    }
    tr.appendChild(td2);
    tbody.appendChild(tr);
  }

  // Stroke colour (not for label)
  if (defaults.stroke_color !== null || elementType === "dimension") {
    const tr = document.createElement("tr");
    const td1 = document.createElement("td"); td1.textContent = "Stroke";
    tr.appendChild(td1);
    const td2 = document.createElement("td");
    const inp = document.createElement("input");
    inp.type = "color";
    inp.className = "prop-color-input";
    inp.value = props.stroke_color || defaults.stroke_color || "#666666";
    inp.addEventListener("change", () => saveStyleProp("stroke_color", inp.value));
    td2.appendChild(inp);
    if (props.stroke_color) {
      const rst = document.createElement("button");
      rst.textContent = "Reset";
      rst.className = "prop-reset-btn";
      rst.addEventListener("click", () => saveStyleProp("stroke_color", null));
      td2.appendChild(rst);
    }
    tr.appendChild(td2);
    tbody.appendChild(tr);
  }

  // Stroke style (not for label)
  if (defaults.stroke_style !== null) {
    const tr = document.createElement("tr");
    const td1 = document.createElement("td"); td1.textContent = "Line Style";
    tr.appendChild(td1);
    const td2 = document.createElement("td");
    const sel = document.createElement("select");
    sel.className = "prop-edit-input";
    const curSS = props.stroke_style || defaults.stroke_style;
    for (const [val, lbl] of [["solid", "Solid"], ["dashed", "Dashed"], ["dotted", "Dotted"]]) {
      const opt = document.createElement("option");
      opt.value = val; opt.textContent = lbl;
      if (val === curSS) opt.selected = true;
      sel.appendChild(opt);
    }
    sel.addEventListener("change", () => saveStyleProp("stroke_style", sel.value));
    td2.appendChild(sel);
    tr.appendChild(td2);
    tbody.appendChild(tr);
  }

  // Opacity (all types)
  {
    const tr = document.createElement("tr");
    const td1 = document.createElement("td"); td1.textContent = "Opacity";
    tr.appendChild(td1);
    const td2 = document.createElement("td");
    td2.style.display = "flex";
    td2.style.alignItems = "center";
    td2.style.gap = "6px";
    const slider = document.createElement("input");
    slider.type = "range";
    slider.className = "prop-opacity-input";
    slider.min = "0"; slider.max = "100"; slider.step = "1";
    slider.value = props.opacity != null ? props.opacity : defaults.opacity;
    const readout = document.createElement("span");
    readout.style.fontSize = "11px";
    readout.style.minWidth = "30px";
    readout.textContent = slider.value + "%";
    slider.addEventListener("input", () => { readout.textContent = slider.value + "%"; });
    slider.addEventListener("change", () => saveStyleProp("opacity", parseInt(slider.value)));
    td2.appendChild(slider);
    td2.appendChild(readout);
    tr.appendChild(td2);
    tbody.appendChild(tr);
  }
}

/** Add per-view style override controls for the current variant. */
function addViewOverrideControls(tbody, elemRec, props) {
  if (!elemRec) return;
  const variant = App.state.variant;
  const vLabel = variantLabel(variant);
  const viewOv = (props.view_overrides || {})[variant] || {};

  // Header
  const hdr = document.createElement("tr");
  const hTd = document.createElement("td");
  hTd.colSpan = 2;
  hTd.style.paddingTop = "8px";
  hTd.style.fontStyle = "italic";
  hTd.style.fontSize = "11px";
  hTd.textContent = `Override for ${vLabel}`;
  hdr.appendChild(hTd);
  tbody.appendChild(hdr);

  async function saveViewOverride(key, value) {
    const cur = (App.state.elements || []).find(e => e.id === elemRec.id);
    const curProps = cur ? parseProps(cur) : props;
    const allOverrides = { ...(curProps.view_overrides || {}) };
    const curView = { ...(allOverrides[variant] || {}) };
    if (value === null) {
      delete curView[key];
    } else {
      curView[key] = value;
    }
    if (Object.keys(curView).length === 0) {
      delete allOverrides[variant];
    } else {
      allOverrides[variant] = curView;
    }
    const newProps = { ...curProps, view_overrides: allOverrides };
    await fetch(`/api/elements/${elemRec.id}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ properties: newProps }),
    });
    await reloadAfterChange();
  }

  // Opacity override for current view
  {
    const tr = document.createElement("tr");
    const td1 = document.createElement("td"); td1.textContent = "Opacity";
    tr.appendChild(td1);
    const td2 = document.createElement("td");
    td2.style.display = "flex";
    td2.style.alignItems = "center";
    td2.style.gap = "6px";
    const slider = document.createElement("input");
    slider.type = "range";
    slider.className = "prop-opacity-input";
    slider.min = "0"; slider.max = "100"; slider.step = "1";
    slider.value = viewOv.opacity != null ? viewOv.opacity : "";
    const readout = document.createElement("span");
    readout.style.fontSize = "11px";
    readout.style.minWidth = "30px";
    readout.textContent = viewOv.opacity != null ? viewOv.opacity + "%" : "—";
    slider.addEventListener("input", () => { readout.textContent = slider.value + "%"; });
    slider.addEventListener("change", () => saveViewOverride("opacity", parseInt(slider.value)));
    td2.appendChild(slider);
    td2.appendChild(readout);
    if (viewOv.opacity != null) {
      const clr = document.createElement("button");
      clr.textContent = "Clear";
      clr.className = "prop-reset-btn";
      clr.addEventListener("click", () => saveViewOverride("opacity", null));
      td2.appendChild(clr);
    }
    tr.appendChild(td2);
    tbody.appendChild(tr);
  }

  // Fill override for current view
  {
    const tr = document.createElement("tr");
    const td1 = document.createElement("td"); td1.textContent = "Fill";
    tr.appendChild(td1);
    const td2 = document.createElement("td");
    const inp = document.createElement("input");
    inp.type = "color";
    inp.className = "prop-color-input";
    inp.value = viewOv.fill_color || props.fill_color || "#333333";
    inp.addEventListener("change", () => saveViewOverride("fill_color", inp.value));
    td2.appendChild(inp);
    if (viewOv.fill_color) {
      const clr = document.createElement("button");
      clr.textContent = "Clear";
      clr.className = "prop-reset-btn";
      clr.addEventListener("click", () => saveViewOverride("fill_color", null));
      td2.appendChild(clr);
    }
    tr.appendChild(td2);
    tbody.appendChild(tr);
  }
}

/** Add product URL field to properties panel. */
function addProductUrlField(tbody, elemRec, props) {
  const url = (props && props.product_url) || "";
  if (!elemRec && !url) return;
  const tr = document.createElement("tr");
  const td1 = document.createElement("td"); td1.textContent = "Link";
  tr.appendChild(td1);
  const td2 = document.createElement("td");
  td2.style.display = "flex";
  td2.style.gap = "4px";
  td2.style.alignItems = "center";
  const inp = document.createElement("input");
  inp.type = "text";
  inp.className = "prop-edit-input";
  inp.placeholder = "https://...";
  inp.style.flex = "1";
  inp.value = url;
  if (!elemRec) inp.readOnly = true;
  if (elemRec) {
    inp.addEventListener("change", async () => {
      const cur = (App.state.elements || []).find(e => e.id === elemRec.id);
      const curProps = cur ? parseProps(cur) : props;
      const newProps = { ...curProps, product_url: inp.value || null };
      await fetch(`/api/elements/${elemRec.id}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ properties: newProps }),
      });
      await reloadAfterChange();
    });
  }
  td2.appendChild(inp);
  if (url && /^https?:\/\//.test(url)) {
    const btn = document.createElement("button");
    btn.textContent = "Open";
    btn.className = "prop-url-open";
    btn.addEventListener("click", () => window.open(url, "_blank"));
    td2.appendChild(btn);
  }
  tr.appendChild(td2);
  tbody.appendChild(tr);
}

function addElementActions(tbody, elemRec) {
  if (!elemRec) return;
  const props = typeof elemRec.properties === "string"
    ? JSON.parse(elemRec.properties) : (elemRec.properties || {});
  // Determine which layouts are currently checked
  const allVariants = (App.state.variants || []).map(v => v.name);
  let checked;
  if (Array.isArray(props.variants)) {
    checked = new Set(props.variants);
  } else if (elemRec.variant) {
    checked = new Set([elemRec.variant]);
  } else {
    checked = new Set(allVariants); // NULL variant = all
  }
  // Layout checkboxes
  const vTr = document.createElement("tr");
  const vTd1 = document.createElement("td");
  vTd1.textContent = "Layouts";
  vTd1.style.verticalAlign = "top";
  vTr.appendChild(vTd1);
  const vTd2 = document.createElement("td");
  const boxes = [];
  for (const v of (App.state.variants || [])) {
    const val = v.name, lbl = v.label;
    const label = document.createElement("label");
    label.className = "prop-checkbox-label";
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.value = val;
    cb.checked = checked.has(val);
    cb.addEventListener("change", async () => {
      const selected = boxes.filter(b => b.checked).map(b => b.value);
      // Re-read current properties from App.state to avoid stale closure
      const cur = (App.state.elements || []).find(e => e.id === elemRec.id);
      const curProps = cur ? parseProps(cur) : props;
      const newProps = { ...curProps, variants: selected };
      await fetch(`/api/elements/${elemRec.id}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ properties: newProps }),
      });
      await reloadAfterChange();
    });
    boxes.push(cb);
    label.appendChild(cb);
    label.appendChild(document.createTextNode(" " + lbl));
    vTd2.appendChild(label);
  }
  vTr.appendChild(vTd2);
  tbody.appendChild(vTr);
  // Delete button (uses shared helper)
  addFormulaDeleteButton(tbody, elemRec.name, elemRec.type);
}

/** Add layout checkboxes + delete for IW walls using variant_exclusions. */
function addWallActions(tbody, wallName, elemRec) {
  const allExclusions = App.state.exclusions || {};
  const allVariants = App.state.variants || [];

  // Layout checkboxes — one per variant, all interactive
  const vTr = document.createElement("tr");
  const vTd1 = document.createElement("td");
  vTd1.textContent = "Layouts";
  vTd1.style.verticalAlign = "top";
  vTr.appendChild(vTd1);
  const vTd2 = document.createElement("td");
  for (const v of allVariants) {
    const varExcl = allExclusions[v.name] || {};
    const excludedWalls = new Set(varExcl.wall || []);
    const label = document.createElement("label");
    label.className = "prop-checkbox-label";
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.checked = !excludedWalls.has(wallName);
    cb.addEventListener("change", async () => {
      try {
        const resp = await fetch("/api/exclusions", {
          method: "PUT",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify({
            variant: v.name,
            element_type: "wall",
            element_name: wallName,
            excluded: !cb.checked,
          }),
        });
        if (!resp.ok) {
          const data = await resp.json();
          showToast(data.error || "Failed to update", "error");
          cb.checked = !cb.checked;
          return;
        }
        await reloadAfterChange();
      } catch (e) {
        showToast("Failed to toggle wall: " + e.message, "error");
        cb.checked = !cb.checked;
      }
    });
    label.appendChild(cb);
    label.appendChild(document.createTextNode(" " + v.label));
    vTd2.appendChild(label);
  }
  vTr.appendChild(vTd2);
  tbody.appendChild(vTr);

  // Delete button (uses shared helper with dependents check)
  if (elemRec) {
    addFormulaDeleteButton(tbody, wallName);
  }
}

/** Show a modal dialog with custom buttons. Returns a Promise that resolves
 *  with the button's data-action value, or null if cancelled/dismissed. */
function showDeleteDialog(title, message, buttons) {
  return new Promise(resolve => {
    const overlay = document.createElement("div");
    overlay.className = "dialog-overlay";
    const dialog = document.createElement("div");
    dialog.className = "dialog";
    const h3 = document.createElement("h3");
    h3.textContent = title;
    dialog.appendChild(h3);
    if (message) {
      const p = document.createElement("p");
      p.style.fontSize = "12px";
      p.style.marginBottom = "8px";
      p.style.whiteSpace = "pre-line";
      p.textContent = message;
      dialog.appendChild(p);
    }
    const btnRow = document.createElement("div");
    btnRow.className = "dialog-buttons";
    btnRow.style.flexWrap = "wrap";
    for (const b of buttons) {
      const btn = document.createElement("button");
      btn.textContent = b.label;
      btn.className = b.className || "dialog-btn-cancel";
      btn.addEventListener("click", () => {
        overlay.remove();
        resolve(b.action);
      });
      btnRow.appendChild(btn);
    }
    dialog.appendChild(btnRow);
    overlay.appendChild(dialog);
    overlay.addEventListener("click", (e) => {
      if (e.target === overlay) { overlay.remove(); resolve(null); }
    });
    document.body.appendChild(overlay);
  });
}

/** Add delete button for formula-driven elements.  Disabled when other
 *  elements depend on this one (dependents check via API). */
async function addFormulaDeleteButton(tbody, elemName, elemType) {
  const tr = document.createElement("tr");
  const td = document.createElement("td");
  td.colSpan = 2;
  const btn = document.createElement("button");
  btn.textContent = "Delete";
  btn.className = "prop-delete-btn";

  const isCatalogType = elemType === "furniture" || elemType === "appliance" || elemType === "fixture";

  // Check for dependents — include in confirmation message
  let depNames = [];
  try {
    const depResp = await fetch(
      `/api/formulas/${encodeURIComponent(elemName)}/dependents?type=element`);
    if (depResp.ok) {
      const deps = await depResp.json();
      depNames = [...new Set(deps.map(d => d.element_name))];
    }
  } catch (_) { /* proceed without dep info */ }

  btn.addEventListener("click", async () => {
    let depMsg = "";
    if (depNames.length > 0) {
      depMsg = `${depNames.join(", ")} will be rebased to ${elemName}'s parent references.`;
    }

    let action;
    if (isCatalogType) {
      // Three-button dialog for catalog items
      action = await showDeleteDialog(
        `Delete ${elemName}?`,
        depMsg || null,
        [
          { label: "Delete", action: "delete", className: "dialog-btn-primary" },
          { label: "Also delete from ADD menu", action: "delete_catalog", className: "dialog-btn-danger" },
          { label: "Cancel", action: null, className: "dialog-btn-cancel" },
        ]
      );
    } else {
      let msg = `Delete ${elemName}?`;
      if (depMsg) msg += "\n\n" + depMsg;
      action = confirm(msg) ? "delete" : null;
    }

    if (!action) return;

    const v = App.state.variant || "standard";
    const resp = await fetch(
      `/api/formulas/${encodeURIComponent(elemName)}/element?variant=${encodeURIComponent(v)}`,
      { method: "DELETE" }
    );
    if (!resp.ok) {
      const data = await resp.json().catch(() => ({}));
      showToast(data.error || "Delete failed", "error");
      return;
    }
    const result = await resp.json();

    // If user chose to also delete from catalog, remove catalog entry
    if (action === "delete_catalog") {
      const catalogKey = elemName;
      await fetch(`/api/catalog/${encodeURIComponent(catalogKey)}`, { method: "DELETE" });
    }

    clearSelection();
    await reloadAfterChange();
    const rebased = result.rebased || [];
    let toast = rebased.length > 0
      ? `Deleted ${elemName} (re-based ${rebased.join(", ")})`
      : `Deleted ${elemName}`;
    if (action === "delete_catalog") toast += " and removed from ADD menu";
    showToast(toast, "success");
  });
  td.appendChild(btn);
  tr.appendChild(td);
  tbody.appendChild(tr);
}

/** Phase 12b: Show formula section in properties panel for elements with formulas. */
async function addFormulaSection(tbody, elemName) {
  try {
    const resp = await fetch(`/api/formulas/${encodeURIComponent(elemName)}`);
    if (!resp.ok) return;
    const formulas = await resp.json();
    if (!formulas || formulas.length === 0) return;

    addPropRow(tbody, "\u2014", "Formula");

    for (const f of formulas) {
      const fj = typeof f.formula_json === "string" ? JSON.parse(f.formula_json) : f.formula_json;
      addPropRow(tbody, "Type", fj.type || "unknown");
      addPropRow(tbody, "Param", f.param_name);

      // Lock toggle
      const lockTr = document.createElement("tr");
      const lockTd1 = document.createElement("td");
      lockTd1.textContent = f.locked ? "\u{1F512} Locked" : "\u{1F513} Unlocked";
      lockTr.appendChild(lockTd1);
      const lockTd2 = document.createElement("td");
      const lockBtn = document.createElement("button");
      lockBtn.className = "prop-btn";
      lockBtn.textContent = f.locked ? "Unlock" : "Lock";
      lockBtn.title = f.locked
        ? "Unlock: allow upstream changes to propagate"
        : "Lock: freeze at current computed value";
      lockBtn.addEventListener("click", async () => {
        const newLocked = !f.locked;
        await fetch(`/api/formulas/${encodeURIComponent(elemName)}/${f.param_name}/lock`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ locked: newLocked, locked_value: null }),
        });
        await loadGeometry();
        // Re-show properties to refresh lock state
        if (App.state.selection && App.state.selection.name === elemName) {
          showProperties(App.state.selection.type, elemName, App.state.selection.data);
        }
      });
      lockTd2.appendChild(lockBtn);
      lockTr.appendChild(lockTd2);
      tbody.appendChild(lockTr);
    }

    // Dependencies
    const depsResp = await fetch(`/api/formulas/${encodeURIComponent(elemName)}/deps`);
    if (depsResp.ok) {
      const deps = await depsResp.json();
      if (deps.length > 0) {
        const depNames = deps.map(d => `${d.dep_name} (${d.dep_type})`).join(", ");
        addPropRow(tbody, "Depends on", depNames);
      }
    }

    // Dependents
    const deptsResp = await fetch(`/api/formulas/${encodeURIComponent(elemName)}/dependents`);
    if (deptsResp.ok) {
      const depts = await deptsResp.json();
      if (depts.length > 0) {
        const deptNames = [...new Set(depts.map(d => d.element_name))].join(", ");
        addPropRow(tbody, "Used by", deptNames);
      }
    }
  } catch (e) {
    // Formula fetch failed — not critical, skip silently
  }
}

function addPropRow(tbody, label, value, editable = false, constName = null) {
  const tr = document.createElement("tr");
  const tdLabel = document.createElement("td");
  tdLabel.textContent = label;
  tr.appendChild(tdLabel);

  const tdVal = document.createElement("td");
  if (editable && constName) {
    const inp = document.createElement("input");
    inp.type = "text";
    inp.value = value;
    inp.dataset.constName = constName;
    inp.addEventListener("change", () => handleConstantEdit(constName, inp.value));
    // SEL-15: highlight dependent elements on focus
    inp.addEventListener("focus", () => { inp.select(); highlightConstantDeps(constName); });
    inp.addEventListener("blur", () => { setTimeout(clearConstantHighlights, 200); });
    tdVal.appendChild(inp);
  } else {
    tdVal.textContent = value;
  }
  tr.appendChild(tdVal);
  tbody.appendChild(tr);
}

/** Return the width-controlling constant name for an opening, or null. */
function findWidthConstant(openingName) {
  // Outer openings: O1-O6 use {name}_WIDTH; O4,O7-O11 use {name}_HALF_WIDTH
  // Rough openings: RO1→IW1_RO_WIDTH, RO2→IW2_RO_WIDTH, RO3→RO3_WIDTH,
  //   RO4→IW2_RO_WIDTH, RO5→IW4_RO_WIDTH, RO6→IW11_RO_WIDTH, RO7→IW9_RO_WIDTH
  const RO_MAP = {
    RO1: "IW1_RO_WIDTH", RO2: "IW2_RO_WIDTH", RO3: "RO3_WIDTH",
    RO4: "IW2_RO_WIDTH", RO5: "IW4_RO_WIDTH", RO6: "IW11_RO_WIDTH",
    RO7: "IW9_RO_WIDTH",
  };
  const n = openingName.toUpperCase();
  if (RO_MAP[n]) return RO_MAP[n];
  // Try direct WIDTH, then HALF_WIDTH
  for (const suffix of ["_WIDTH", "_HALF_WIDTH"]) {
    if (App.state.constants.find(c => c.name === n + suffix)) return n + suffix;
  }
  return null;
}

/* ========== SEL-15: Constant Dependency Highlighting ========== */

/** Reverse map: constant name → IW wall names that depend on it. */
const CONSTANT_TO_IW = {};
(function() {
  const map = {
    IW1: "IW1_OFFSET_FROM_W9", IW2: "IW2_DIST_W2W5",
    IW2S: "IW2S_W2REF_OFFSET", IW3: "IW3_DIST_W2W5",
    IW4: "IW4_GAP_IW11", IW5: "IW5_S_OFFSET_FROM_IW1",
    IW6: "IW6_OFFSET_FROM_W6", IW7: "IW7_OFFSET_FROM_W18W1",
    IW9: "IW3_OFFSET_IW9", IW11: "IW9_IW11_GAP",
    IW12: "IW12_S_OFFSET_W18W1",
  };
  for (const [iw, cname] of Object.entries(map)) {
    if (!CONSTANT_TO_IW[cname]) CONSTANT_TO_IW[cname] = [];
    CONSTANT_TO_IW[cname].push(iw);
  }
})();

/** Constant prefix → furniture/appliance item names. */
const CONSTANT_PREFIX_TO_ITEMS = {
  BED_: ["bed"], DRESSER_: ["dresser"], SHELVES_: ["shelves"],
  LOVESEAT_: ["loveseat"], DESK_: ["desk"], CHAIR_: ["chair"],
  SOFA_: ["sofa"], ROCKER_: ["rocker"],
  APPLIANCE_: ["washer", "dryer"], COUNTER_: ["counter"],
  STOVE_: ["stove"], FRIDGE_: ["fridge"], DW_: ["dishwasher"],
  MINIK_: ["minik_counter"], ICE_: ["ice_maker"],
  TOILET_: ["toilet"], SINK_: ["bath_sink"],
  WH_: ["water_heater"],
};

function highlightConstantDeps(constName) {
  clearConstantHighlights();
  const firstOrder = new Set();
  const secondOrder = new Set();
  // 1. Reverse IW map
  const iws = CONSTANT_TO_IW[constName] || [];
  for (const iw of iws) {
    firstOrder.add(iw);
    // Hosted openings are second-order
    if (typeof IW_HOSTED_OPENINGS !== "undefined") {
      for (const ro of (IW_HOSTED_OPENINGS[iw] || [])) secondOrder.add(ro);
    }
  }
  // 2. Prefix-based furniture/appliance mapping
  for (const [prefix, items] of Object.entries(CONSTANT_PREFIX_TO_ITEMS)) {
    if (constName.startsWith(prefix)) {
      for (const item of items) firstOrder.add(item);
    }
  }
  // 3. Opening width constants (O1_WIDTH → O1, RO3_WIDTH → RO3, etc.)
  const m = constName.match(/^((?:RO?\d+|O\d+))_/i);
  if (m) firstOrder.add(m[1].toUpperCase());
  // 4. Apply CSS classes
  document.querySelectorAll("[data-name]").forEach(el => {
    const name = el.getAttribute("data-name");
    if (!name) return;
    const upper = name.toUpperCase();
    if (firstOrder.has(name) || firstOrder.has(upper)) {
      el.classList.add("const-dep-first");
    } else if (secondOrder.has(name) || secondOrder.has(upper)) {
      el.classList.add("const-dep-second");
    }
  });
}

function clearConstantHighlights() {
  document.querySelectorAll(".const-dep-first, .const-dep-second").forEach(el => {
    el.classList.remove("const-dep-first", "const-dep-second");
  });
}

function findRelatedConstants(elementName) {
  const name = elementName.toUpperCase().replace(/^IW/, "IW");
  return App.state.constants.filter(c => {
    const cn = c.name.toUpperCase();
    if (name.startsWith("O") && name.length <= 3) {
      return cn.startsWith(name + "_") || cn.includes("_" + name + "_") || cn === name + "_WIDTH";
    }
    if (name.startsWith("IW")) {
      return cn.includes(name) || cn.startsWith(name + "_");
    }
    if (name.startsWith("RO")) {
      return cn.includes(name) || cn.startsWith(name + "_");
    }
    // furniture/appliance: match by keyword
    const keywords = {
      bed: "BED_", dresser: "DRESSER_", shelves: "SHELVES_",
      washer: "APPLIANCE_", dryer: "APPLIANCE_", counter: "COUNTER_",
    };
    const prefix = keywords[name.toLowerCase()];
    if (prefix) return cn.startsWith(prefix);
    return false;
  });
}


/* ========== CONSTANTS TABLE ========== */

function renderConstantsTable() {
  const tbody = App.els["constants-table"].querySelector("tbody");
  tbody.innerHTML = "";

  let filtered = App.state.constants;

  // Category filter
  const catFilter = App.els["const-category-filter"].value;
  if (catFilter) {
    filtered = filtered.filter(c => c.category === catFilter);
  }

  // Search filter
  const search = App.els["const-search"].value.toLowerCase();
  if (search) {
    filtered = filtered.filter(c =>
      c.name.toLowerCase().includes(search) ||
      c.description.toLowerCase().includes(search)
    );
  }

  // Sort
  const key = App.state.constSortKey;
  const asc = App.state.constSortAsc ? 1 : -1;
  filtered.sort((a, b) => {
    if (a[key] < b[key]) return -asc;
    if (a[key] > b[key]) return asc;
    return 0;
  });

  for (const c of filtered) {
    const tr = document.createElement("tr");
    tr.dataset.name = c.name;

    const tdName = document.createElement("td");
    tdName.textContent = c.name;
    tdName.title = c.description || c.name;
    tr.appendChild(tdName);

    const tdVal = document.createElement("td");
    const inp = document.createElement("input");
    inp.type = "text";
    inp.value = formatConstValue(c);
    inp.dataset.origValue = c.value;
    inp.dataset.constName = c.name;
    inp.addEventListener("focus", () => inp.select());
    inp.addEventListener("change", () => handleConstantEdit(c.name, inp.value));
    inp.addEventListener("input", () => {
      inp.classList.toggle("changed", inp.value !== inp.dataset.origValue);
    });
    tdVal.appendChild(inp);
    tr.appendChild(tdVal);

    const tdCat = document.createElement("td");
    tdCat.textContent = c.category;
    tdCat.style.color = categoryColor(c.category);
    tr.appendChild(tdCat);

    tbody.appendChild(tr);
  }
}

function formatConstValue(c) {
  if (c.unit === "ft") {
    const inches = Math.round(c.value * 1200) / 100;
    const inStr = inches.toFixed(2).replace(/0+$/, '').replace(/\.$/, '');
    return `${inStr}"`;
  }
  const suffix = {cm: " cm", mm: " mm", m: " m"}[c.unit];
  if (suffix) {
    return c.value.toFixed(6).replace(/0+$/, '').replace(/\.$/, '') + suffix;
  }
  return c.value.toFixed(6).replace(/0+$/, '').replace(/\.$/, '');
}

function categoryColor(cat) {
  const colors = {
    wall: "#89b4fa",
    interior_wall: "#74c7ec",
    opening: "#f9e2af",
    appliance: "#fab387",
    furniture: "#a6e3a1",
    fixture: "#cba6f7",
    geometry: "#f38ba8",
    construction: "#94e2d5",
    misc: "#8888aa",
  };
  return colors[cat] || "#8888aa";
}

/**
 * Unit-aware dimension parser (CT-7a through CT-7j, CT-8).
 * Parses user input with optional unit suffixes and returns feet.
 * Returns NaN on invalid input.
 */
function parseDimension(text) {
  if (text == null) return NaN;
  const trimmed = text.trim();
  if (!trimmed) return NaN;

  // CT-8: fraction shortcut — "1/3", "3/4", etc.
  if (trimmed.includes("/") && !/['"\s]/.test(trimmed)) {
    const parts = trimmed.split("/");
    if (parts.length === 2) {
      const num = parseFloat(parts[0]), den = parseFloat(parts[1]);
      if (!isNaN(num) && !isNaN(den) && den !== 0) return num / den;
      return NaN;
    }
  }

  const conversions = {
    ft: 1, feet: 1,
    "in": 1 / 12, inches: 1 / 12,
    cm: 1 / 30.48, centimeters: 1 / 30.48,
    mm: 1 / 304.8, millimeters: 1 / 304.8,
    m: 1 / 0.3048, meters: 1 / 0.3048,
  };

  const tokenRe = /(-?(?:\d+\.?\d*|\.\d+))\s*(?:(['\u2032])|(["\u2033])|(feet|ft|inches|in|centimeters|cm|millimeters|mm|meters|m)(?![a-z]))?/gi;
  const tokens = [];
  let match;
  while ((match = tokenRe.exec(trimmed)) !== null) {
    tokens.push(match);
    if (tokens.length >= 4) break;  // CT-7i: stop after 4th (only use 3)
  }
  if (tokens.length === 0) return NaN;

  let result = 0, bareCount = 0;
  for (let i = 0; i < Math.min(tokens.length, 3); i++) {
    const [, numStr, footMark, inchMark, wordUnit] = tokens[i];
    const value = parseFloat(numStr);
    if (isNaN(value)) return NaN;

    if (footMark) {
      result += value;
    } else if (inchMark) {
      result += value / 12;
    } else if (wordUnit) {
      const factor = conversions[wordUnit.toLowerCase()];
      if (factor == null) return NaN;
      result += value * factor;
    } else {
      bareCount++;
      result += bareCount === 1 ? value : value / 12;
    }
  }
  return result;
}

async function handleElementPropEdit(elemId, propKey, rawValue) {
  let value;
  if (propKey === "rotation") {
    // Rotation is in degrees, not a dimension
    value = parseFloat(String(rawValue).replace("°", ""));
  } else {
    value = parseDimension(rawValue);
  }
  if (isNaN(value)) {
    showToast(`Invalid value: ${rawValue}`, "error");
    return;
  }
  try {
    // Use formula update endpoint — replaces symbolic constant with literal
    const body = {};
    if (propKey === "rotation") {
      body.rotation_deg = value;
      body.world_rotation = value;
    } else {
      body[propKey] = value;
    }
    const resp = await apiFetch(`/api/elements/${elemId}/update-formula`, {
      method: "PUT",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify(body),
    });
    if (!resp.ok) throw new Error((await resp.json()).error);
    const msg = propKey === "rotation" ? `${propKey} = ${value}°` : `${propKey} = ${fmtFtIn(value)}`;
    showToast(msg, "success");
    loadElements();
    loadGeometry();
  } catch (e) {
    showToast(`Error: ${e.message}`, "error");
  }
}

async function handleConstantEdit(name, rawValue) {
  const value = parseDimension(rawValue);

  if (isNaN(value)) {
    showToast(`Invalid value: ${rawValue}`, "error");
    return;
  }

  try {
    const resp = await fetch(`/api/constants/${name}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ value }),
    });
    if (!resp.ok) throw new Error((await resp.json()).error);
    // Update local state
    const c = App.state.constants.find(x => x.name === name);
    if (c) c.value = value;
    const disp = c ? formatConstValue(c) : value;
    showToast(`${name} = ${disp}`, "success");
  } catch (e) {
    showToast(`Error: ${e.message}`, "error");
  }
}


/* ========== OUTLINE TABLE ========== */

async function loadOutlineTable() {
  const [resp, ovResp] = await Promise.all([
    apiFetch("/api/outline?_=" + Date.now()),
    apiFetch("/api/inner-wall-overrides"),
  ]);
  const chain = await resp.json();
  const overrides = await ovResp.json();
  App.state.outlineChain = chain;
  App.state.innerWallOverrides = overrides;
  const n = chain.length;

  // Update closure indicator
  updateClosureIndicator();

  const tbody = App.els["outline-table"].querySelector("tbody");
  tbody.innerHTML = "";
  for (const seg of chain) {
    const tr = document.createElement("tr");
    tr.dataset.seq = seg.seq;
    if (seg.seq === App.state.outlineSelectedSeq) {
      tr.classList.add("outline-selected");
    }
    tr.addEventListener("click", () => selectOutlineRow(seg.seq));

    // Seq column
    const tdSeq = document.createElement("td");
    tdSeq.textContent = seg.seq;
    tr.appendChild(tdSeq);

    // Type column
    const tdType = document.createElement("td");
    tdType.textContent = seg.seg_type;
    tr.appendChild(tdType);

    // Dist/R column — editable for non-solved segments
    const tdDist = document.createElement("td");
    const isSolvedDist = seg.seq === 0 || seg.seq === n - 2;
    const distVal = seg.seg_type === "L" ? (seg.distance || 0) : (seg.radius || 0);
    if (isSolvedDist && seg.seg_type === "L") {
      tdDist.textContent = fmtFtIn(distVal);
      tdDist.classList.add("solved");
      tdDist.title = "Solved by closure";
    } else {
      const inp = document.createElement("input");
      inp.type = "text";
      inp.value = fmtFtIn(distVal);
      inp.dataset.origValue = String(distVal);
      inp.dataset.seq = seg.seq;
      inp.addEventListener("focus", () => inp.select());
      inp.addEventListener("click", (e) => e.stopPropagation());
      inp.addEventListener("change", () => handleOutlineEdit(seg.seq, "dist_or_radius", inp.value));
      tdDist.appendChild(inp);
    }
    tr.appendChild(tdDist);

    // Sweep column — editable for arcs, except closure arc (last seg)
    const tdSweep = document.createElement("td");
    const isSolvedSweep = seg.seq === n - 1;
    if (seg.seg_type !== "L") {
      const sweepDeg = (seg.sweep || 0) * 180 / Math.PI;
      if (isSolvedSweep) {
        tdSweep.textContent = fmtDeg(sweepDeg);
        tdSweep.classList.add("solved");
        tdSweep.title = "Solved by closure";
      } else {
        const inp = document.createElement("input");
        inp.type = "text";
        inp.value = fmtDeg(sweepDeg);
        inp.dataset.origValue = String(seg.sweep);
        inp.dataset.seq = seg.seq;
        inp.addEventListener("focus", () => inp.select());
        inp.addEventListener("click", (e) => e.stopPropagation());
        inp.addEventListener("change", () => handleOutlineSweepEdit(seg.seq, inp.value));
        tdSweep.appendChild(inp);
      }
    } else {
      tdSweep.textContent = "\u2014";
    }
    tr.appendChild(tdSweep);

    // End name column
    const tdEnd = document.createElement("td");
    tdEnd.textContent = seg.end_name;
    tr.appendChild(tdEnd);

    // W Override column — shows button if override exists or to add one
    const tdOv = document.createElement("td");
    tdOv.className = "td-override";
    const innerIdx = _outlineToInnerIndex(seg.seq, chain);
    if (innerIdx !== null) {
      // Check if this segment is a span start, mid-span, or standalone
      const spanInfo = _getSpanInfo(innerIdx, overrides);
      if (spanInfo.type === "start") {
        const btn = document.createElement("button");
        btn.className = "override-btn has-override";
        btn.textContent = "\u2726";
        btn.title = `Edit W override (span ${spanInfo.start}\u2013${spanInfo.end})`;
        btn.addEventListener("click", (e) => {
          e.stopPropagation();
          openOverrideEditor(spanInfo.start, seg.end_name, null, spanInfo.end);
        });
        tdOv.appendChild(btn);
      } else if (spanInfo.type === "mid") {
        const btn = document.createElement("button");
        btn.className = "override-btn has-override";
        btn.textContent = "\u2502";
        btn.title = `Part of span ${spanInfo.start}\u2013${spanInfo.end}`;
        btn.addEventListener("click", (e) => {
          e.stopPropagation();
          openOverrideEditor(spanInfo.start, seg.end_name, null, spanInfo.end);
        });
        tdOv.appendChild(btn);
      } else {
        // Standalone: either has single-segment override or none
        const hasOverride = String(innerIdx) in overrides;
        const btn = document.createElement("button");
        btn.className = "override-btn" + (hasOverride ? " has-override" : "");
        btn.textContent = hasOverride ? "\u2726" : "\u2727";
        btn.title = hasOverride ? "Edit W override" : "Add W override";
        btn.addEventListener("click", (e) => {
          e.stopPropagation();
          openOverrideEditor(innerIdx, seg.end_name);
        });
        tdOv.appendChild(btn);
      }
    }
    tr.appendChild(tdOv);

    tbody.appendChild(tr);
  }
}

/**
 * Map outline table seq to inner_segs index.
 *
 * The outline chain (DB) starts at F2 (seq 0 = F2→F5), while
 * inner_segs starts at W1 (index 0 = W1→W2).  inner_segs mirrors
 * outline_segs which starts one segment earlier (F1→F2), so
 * inner_segs index = (outline seq + 1) % n.
 */
function _outlineToInnerIndex(seq, chain) {
  if (seq < 0 || seq >= chain.length) return null;
  return (seq + 1) % chain.length;
}

/**
 * Determine span info for a given inner_seg index relative to overrides.
 * Returns { type: "start"|"mid"|"none", start, end } where:
 *   - "start": this index is the start of a span override
 *   - "mid": this index is inside a span override (not the start)
 *   - "none": no span covers this index
 */
function _getSpanInfo(innerIdx, overrides) {
  // Check if this index is a span start
  const key = String(innerIdx);
  if (key in overrides) {
    const chain = overrides[key];
    const spanEnd = chain[0]?.span_end;
    if (spanEnd != null && spanEnd > innerIdx) {
      return { type: "start", start: innerIdx, end: spanEnd };
    }
    // Single-segment override — handled by caller
    return { type: "none" };
  }
  // Check if this index is mid-span of another override
  for (const [startStr, chain] of Object.entries(overrides)) {
    const start = parseInt(startStr, 10);
    const spanEnd = chain[0]?.span_end;
    if (spanEnd != null && innerIdx > start && innerIdx <= spanEnd) {
      return { type: "mid", start, end: spanEnd };
    }
  }
  return { type: "none" };
}

// ---------------------------------------------------------------------------
// Override geometry helpers (JS port of gen_provider.py)
// ---------------------------------------------------------------------------

/**
 * Compute compass bearing (degrees) at the start of an inner segment.
 * Mirrors _seg_start_bearing() in gen_provider.py.
 */
function segStartBearing(seg, points) {
  let dx, dy;
  if (seg.type === "line") {
    const s = points[seg.start], e = points[seg.end];
    dx = e[0] - s[0];
    dy = e[1] - s[1];
  } else {
    const c = points[seg.center], s = points[seg.start];
    const rx = s[0] - c[0], ry = s[1] - c[1];
    if (seg.direction === "CW") {
      dx = ry; dy = -rx;
    } else {
      dx = -ry; dy = rx;
    }
  }
  return ((Math.atan2(dx, dy) * 180 / Math.PI) % 360 + 360) % 360;
}

/**
 * Walk a parametric override chain, returning a polyline.
 * Mirrors walk_override_chain() in gen_provider.py.
 * Bearings use compass convention (0=N, 90=E, 180=S, 270=W).
 */
function walkOverrideChain(chain, startPt, startBearingDeg) {
  const polyline = [startPt.slice()];
  let cx = startPt[0], cy = startPt[1];
  let bearing = startBearingDeg;

  for (const sub of chain) {
    if (sub.seg_type === "L") {
      const brg = sub.bearing || 0;
      const dist = sub.distance || 0;
      const rad = brg * Math.PI / 180;
      cx += dist * Math.sin(rad);
      cy += dist * Math.cos(rad);
      polyline.push([cx, cy]);
      bearing = brg;
    } else {
      const radius = sub.radius || 0;
      const sweepDeg = sub.sweep || 0;
      const nPts = sub.n_pts || 20;
      const sweepRad = sweepDeg * Math.PI / 180;
      const dirX = Math.sin(bearing * Math.PI / 180);
      const dirY = Math.cos(bearing * Math.PI / 180);

      let acx, acy, entryAngle;
      if (sub.seg_type === "CCW") {
        acx = cx - dirY * radius;
        acy = cy + dirX * radius;
        entryAngle = Math.atan2(cy - acy, cx - acx);
        for (let i = 1; i <= nPts; i++) {
          const a = entryAngle + i * sweepRad / nPts;
          polyline.push([acx + radius * Math.cos(a), acy + radius * Math.sin(a)]);
        }
        bearing = ((bearing - sweepDeg) % 360 + 360) % 360;
      } else {
        acx = cx + dirY * radius;
        acy = cy - dirX * radius;
        entryAngle = Math.atan2(cy - acy, cx - acx);
        for (let i = 1; i <= nPts; i++) {
          const a = entryAngle - i * sweepRad / nPts;
          polyline.push([acx + radius * Math.cos(a), acy + radius * Math.sin(a)]);
        }
        bearing = ((bearing + sweepDeg) % 360 + 360) % 360;
      }
      const last = polyline[polyline.length - 1];
      cx = last[0]; cy = last[1];
    }
  }
  return polyline;
}

// ---------------------------------------------------------------------------
// Override preview rendering
// ---------------------------------------------------------------------------

const OverridePreview = {
  segIndex: null,
  chain: null,
  clickDefineMode: false,
  waypoints: [],

  /** Render override polyline on canvas from current chain state. */
  update(segIndex, chain, spanEnd) {
    this.segIndex = segIndex;
    this.chain = chain;
    const layer = document.getElementById("layer-measure");
    // Remove old preview
    layer.querySelectorAll(".override-preview").forEach(el => el.remove());

    const g = App.state.geometry;
    if (!g || !g.inner_segments || !g.inner_segments[segIndex]) return;
    const seg = g.inner_segments[segIndex];
    const startPt = g.points[seg.start];
    if (!startPt) return;

    const bearing = segStartBearing(seg, g.points);

    // Walk the chain to get a polyline
    const validChain = chain.filter(s => {
      if (s.seg_type === "L") return s.distance > 0;
      return s.radius > 0 && s.sweep > 0;
    });
    if (validChain.length === 0) return;

    const poly = walkOverrideChain(validChain, startPt, bearing);
    if (poly.length < 2) return;

    // Draw polyline
    const pts = poly.map(p => `${p[0]},${-p[1]}`).join(" ");
    const el = svgEl("polyline", {
      points: pts,
      class: "override-preview",
      fill: "none",
      stroke: "#f9e2af",
      "stroke-width": "0.04",
      "stroke-dasharray": "0.12 0.06",
      "stroke-linecap": "round",
      "stroke-linejoin": "round",
      opacity: "0.85",
    });
    layer.appendChild(el);

    // Draw vertex dots
    for (const p of poly) {
      layer.appendChild(svgEl("circle", {
        cx: p[0], cy: -p[1], r: 0.06,
        fill: "#f9e2af", opacity: "0.7",
        class: "override-preview",
      }));
    }

    // Draw endpoint marker (target: where the chain should end)
    // For spans, use the end of the last segment in the span
    const endSegIdx = spanEnd != null ? spanEnd : segIndex;
    const endSeg = g.inner_segments[endSegIdx];
    const endPt = endSeg ? g.points[endSeg.end] : null;
    if (endPt) {
      layer.appendChild(svgEl("circle", {
        cx: endPt[0], cy: -endPt[1], r: 0.08,
        fill: "none", stroke: "#a6e3a1", "stroke-width": "0.02",
        "stroke-dasharray": "0.04 0.04",
        class: "override-preview",
      }));
    }
  },

  /** Clear preview from canvas. */
  clear() {
    this.segIndex = null;
    this.chain = null;
    const layer = document.getElementById("layer-measure");
    layer.querySelectorAll(".override-preview").forEach(el => el.remove());
  },
};

/**
 * Open the inner wall override editor dialog.
 * Shows a dynamic table of sub-segments (line/arc chains).
 * @param {number} segIndex - Inner segment start index
 * @param {string} endName - End point name for display
 * @param {Array} [initialChain] - Pre-filled chain (from click-to-define)
 * @param {number|null} [spanEnd] - Span end index (null = single segment)
 */
function openOverrideEditor(segIndex, endName, initialChain, spanEnd) {
  const existing = App.state.innerWallOverrides[String(segIndex)] || [];
  // Determine initial span_end from existing data or parameter
  let currentSpanEnd = spanEnd != null ? spanEnd
    : (existing.length && existing[0].span_end != null) ? existing[0].span_end
    : null;
  const chain = initialChain
    ? initialChain.map(s => ({ ...s }))
    : existing.length
      ? existing.map(s => ({ ...s }))
      : [{ seg_type: "L", bearing: 0, distance: 0, radius: null, sweep: null, n_pts: 20 }];

  const container = document.createElement("div");
  container.className = "override-editor";

  function updatePreview() {
    OverridePreview.update(segIndex, chain, currentSpanEnd);
  }

  function render() {
    container.innerHTML = "";

    // Span range picker
    const g = App.state.geometry;
    const nSegs = g && g.inner_segments ? g.inner_segments.length : 20;
    const spanRow = document.createElement("div");
    spanRow.className = "override-btn-row";
    spanRow.style.marginBottom = "6px";
    const spanLabel = document.createElement("span");
    spanLabel.textContent = "Span: ";
    spanLabel.style.marginRight = "4px";
    spanRow.appendChild(spanLabel);

    const startLabel = document.createElement("span");
    startLabel.textContent = `seg ${segIndex}`;
    startLabel.style.marginRight = "8px";
    spanRow.appendChild(startLabel);

    const toLabel = document.createElement("span");
    toLabel.textContent = " \u2192 ";
    toLabel.style.marginRight = "4px";
    spanRow.appendChild(toLabel);

    const endSel = document.createElement("select");
    endSel.style.width = "auto";
    for (let i = segIndex; i < nSegs; i++) {
      const opt = document.createElement("option");
      opt.value = i === segIndex ? "" : String(i);
      const segName = g && g.inner_segments[i]
        ? g.inner_segments[i].end : `seg ${i}`;
      opt.textContent = i === segIndex ? `seg ${i} (single)` : `seg ${i} (${segName})`;
      if ((currentSpanEnd == null && i === segIndex) ||
          (currentSpanEnd != null && i === currentSpanEnd)) {
        opt.selected = true;
      }
      endSel.appendChild(opt);
    }
    endSel.onchange = () => {
      const val = endSel.value;
      currentSpanEnd = val ? parseInt(val, 10) : null;
      updatePreview();
    };
    spanRow.appendChild(endSel);
    container.appendChild(spanRow);

    const tbl = document.createElement("table");
    tbl.className = "override-chain-table";
    const thead = document.createElement("thead");
    thead.innerHTML = "<tr><th>#</th><th>Type</th><th>Bearing</th><th>Distance</th><th>Radius</th><th>Sweep</th><th></th></tr>";
    tbl.appendChild(thead);
    const tbody = document.createElement("tbody");

    chain.forEach((sub, i) => {
      const tr = document.createElement("tr");
      const td0 = document.createElement("td");
      td0.textContent = i;
      tr.appendChild(td0);
      // Type selector
      const td1 = document.createElement("td");
      const sel = document.createElement("select");
      for (const t of ["L", "CW", "CCW"]) {
        const opt = document.createElement("option");
        opt.value = t; opt.textContent = t;
        if (t === sub.seg_type) opt.selected = true;
        sel.appendChild(opt);
      }
      sel.onchange = () => { sub.seg_type = sel.value; render(); };
      td1.appendChild(sel);
      tr.appendChild(td1);
      // Bearing (lines only)
      const td2 = document.createElement("td");
      if (sub.seg_type === "L") {
        const inp = document.createElement("input");
        inp.type = "number"; inp.step = "0.01";
        inp.value = sub.bearing != null ? sub.bearing : "";
        inp.onchange = () => { sub.bearing = parseFloat(inp.value); updatePreview(); };
        td2.appendChild(inp);
      } else { td2.textContent = "\u2014"; }
      tr.appendChild(td2);
      // Distance (lines only)
      const td3 = document.createElement("td");
      if (sub.seg_type === "L") {
        const inp = document.createElement("input");
        inp.type = "text";
        inp.value = sub.distance != null ? fmtFtIn(sub.distance) : "";
        inp.onchange = () => { sub.distance = parseDimension(inp.value); updatePreview(); };
        td3.appendChild(inp);
      } else { td3.textContent = "\u2014"; }
      tr.appendChild(td3);
      // Radius (arcs only)
      const td4 = document.createElement("td");
      if (sub.seg_type !== "L") {
        const inp = document.createElement("input");
        inp.type = "text";
        inp.value = sub.radius != null ? fmtFtIn(sub.radius) : "";
        inp.onchange = () => { sub.radius = parseDimension(inp.value); updatePreview(); };
        td4.appendChild(inp);
      } else { td4.textContent = "\u2014"; }
      tr.appendChild(td4);
      // Sweep (arcs only)
      const td5 = document.createElement("td");
      if (sub.seg_type !== "L") {
        const inp = document.createElement("input");
        inp.type = "number"; inp.step = "0.01";
        inp.value = sub.sweep != null ? sub.sweep : "";
        inp.onchange = () => { sub.sweep = parseFloat(inp.value); updatePreview(); };
        td5.appendChild(inp);
      } else { td5.textContent = "\u2014"; }
      tr.appendChild(td5);
      // Remove button
      const td6 = document.createElement("td");
      const rmBtn = document.createElement("button");
      rmBtn.textContent = "\u2212";
      rmBtn.className = "override-rm-btn";
      rmBtn.onclick = () => { chain.splice(i, 1); render(); };
      td6.appendChild(rmBtn);
      tr.appendChild(td6);

      tbody.appendChild(tr);
    });
    tbl.appendChild(tbody);
    container.appendChild(tbl);

    // Button row
    const btnRow = document.createElement("div");
    btnRow.className = "override-btn-row";

    // Add sub-segment
    const addBtn = document.createElement("button");
    addBtn.className = "override-add-btn";
    addBtn.textContent = "+ Add";
    addBtn.onclick = () => {
      chain.push({ seg_type: "L", bearing: 0, distance: 0, radius: null, sweep: null, n_pts: 20 });
      render();
    };
    btnRow.appendChild(addBtn);

    // Compute from default
    const defaultBtn = document.createElement("button");
    defaultBtn.className = "override-add-btn";
    defaultBtn.textContent = "Compute Default";
    defaultBtn.title = "Compute parametric chain from current geometry";
    defaultBtn.onclick = async () => {
      let url = `/api/inner-wall-overrides/${segIndex}/compute-default`;
      if (currentSpanEnd != null) url += `?span_end=${currentSpanEnd}`;
      const resp = await apiFetch(url);
      if (!resp.ok) {
        showToast("No default available for this segment", "error");
        return;
      }
      const data = await resp.json();
      chain.length = 0;
      chain.push(...data.chain);
      render();
    };
    btnRow.appendChild(defaultBtn);

    // Click-to-define
    const clickBtn = document.createElement("button");
    clickBtn.className = "override-add-btn";
    clickBtn.textContent = "Click to Define";
    clickBtn.title = "Click waypoints on the canvas to define the chain";
    clickBtn.onclick = () => {
      Dialog.close();
      OverridePreview.clear();
      startClickDefineMode(segIndex, endName, chain, currentSpanEnd);
    };
    btnRow.appendChild(clickBtn);

    // Remove Override (existing overrides only)
    if (existing.length > 0) {
      const delBtn = document.createElement("button");
      delBtn.className = "override-add-btn";
      delBtn.style.color = "var(--red)";
      delBtn.textContent = "Remove";
      delBtn.onclick = async () => {
        Dialog.close();
        OverridePreview.clear();
        await apiFetch(`/api/inner-wall-overrides/${segIndex}`, { method: "DELETE" });
        await loadOutlineTable();
        await loadGeometry();
      };
      btnRow.appendChild(delBtn);
    }

    container.appendChild(btnRow);

    updatePreview();
  }

  render();

  const spanTitle = currentSpanEnd != null
    ? `W Override \u2014 Seg ${segIndex}\u2013${currentSpanEnd}`
    : `W Override \u2014 Seg ${segIndex} (${endName})`;
  Dialog.show({
    title: spanTitle,
    customContent: container,
    fields: [],
    onCancel() { OverridePreview.clear(); },
    async onSubmit() {
      OverridePreview.clear();
      if (chain.length === 0) {
        showToast("Chain is empty \u2014 use Remove to delete", "error");
        return;
      }
      const payload = { chain };
      if (currentSpanEnd != null) payload.span_end = currentSpanEnd;
      const resp = await apiFetch(`/api/inner-wall-overrides/${segIndex}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      if (!resp.ok) {
        const err = await resp.json();
        showToast(err.error || "Failed to save override", "error");
        return;
      }
      const result = await resp.json();
      if (result.warnings && result.warnings.length > 0) {
        showToast("Saved with warnings: " + result.warnings.join("; "), "warning");
      }
      await loadOutlineTable();
      await loadGeometry();
    },
  });
}

// ---------------------------------------------------------------------------
// Click-to-define mode
// ---------------------------------------------------------------------------

function startClickDefineMode(segIndex, endName, existingChain, spanEnd) {
  const g = App.state.geometry;
  const seg = g.inner_segments[segIndex];
  const startPt = g.points[seg.start];
  const startBearing = segStartBearing(seg, g.points);

  const waypoints = [startPt.slice()];
  OverridePreview.clickDefineMode = true;
  OverridePreview.waypoints = waypoints;

  showToast("Click waypoints on canvas. Double-click or Enter to finish.", "info");

  // Draw current waypoints as preview
  function drawClickPreview() {
    const layer = document.getElementById("layer-measure");
    layer.querySelectorAll(".override-preview").forEach(el => el.remove());

    if (waypoints.length < 1) return;
    const pts = waypoints.map(p => `${p[0]},${-p[1]}`).join(" ");
    layer.appendChild(svgEl("polyline", {
      points: pts,
      class: "override-preview",
      fill: "none",
      stroke: "#f9e2af",
      "stroke-width": "0.04",
      "stroke-dasharray": "0.12 0.06",
      "stroke-linecap": "round",
      "stroke-linejoin": "round",
      opacity: "0.85",
    }));
    for (const p of waypoints) {
      layer.appendChild(svgEl("circle", {
        cx: p[0], cy: -p[1], r: 0.06,
        fill: "#f9e2af", opacity: "0.7",
        class: "override-preview",
      }));
    }
    // Target endpoint (use span end if set)
    const endSegIdx = spanEnd != null ? spanEnd : segIndex;
    const endSeg = g.inner_segments[endSegIdx];
    const endPt = endSeg ? g.points[endSeg.end] : null;
    if (endPt) {
      layer.appendChild(svgEl("circle", {
        cx: endPt[0], cy: -endPt[1], r: 0.08,
        fill: "none", stroke: "#a6e3a1", "stroke-width": "0.02",
        "stroke-dasharray": "0.04 0.04",
        class: "override-preview",
      }));
    }
  }
  drawClickPreview();

  function onClick(e) {
    if (e.button !== 0) return;
    const [wx, wy] = mouseToWorld(e);
    waypoints.push([wx, wy]);
    drawClickPreview();
    e.stopPropagation();
    e.preventDefault();
  }

  function onDblClick(e) {
    e.stopPropagation();
    e.preventDefault();
    finish();
  }

  function onKeyDown(e) {
    if (e.key === "Enter") { finish(); e.preventDefault(); }
    if (e.key === "Escape") { cancel(); e.preventDefault(); }
  }

  function finish() {
    cleanup();
    if (waypoints.length < 2) {
      openOverrideEditor(segIndex, endName, null, spanEnd);
      return;
    }
    // Convert waypoints to parametric chain (line segments)
    const newChain = [];
    for (let i = 1; i < waypoints.length; i++) {
      const dx = waypoints[i][0] - waypoints[i - 1][0];
      const dy = waypoints[i][1] - waypoints[i - 1][1];
      const dist = Math.sqrt(dx * dx + dy * dy);
      const bearing = ((Math.atan2(dx, dy) * 180 / Math.PI) % 360 + 360) % 360;
      newChain.push({
        seg_type: "L", bearing: Math.round(bearing * 100) / 100,
        distance: Math.round(dist * 1e6) / 1e6,
        radius: null, sweep: null, n_pts: 20,
      });
    }
    openOverrideEditor(segIndex, endName, newChain, spanEnd);
  }

  function cancel() {
    cleanup();
    openOverrideEditor(segIndex, endName, null, spanEnd);
  }

  function cleanup() {
    OverridePreview.clickDefineMode = false;
    OverridePreview.waypoints = [];
    OverridePreview.clear();
    const vp = App.els["viewport"];
    vp.removeEventListener("click", onClick, true);
    vp.removeEventListener("dblclick", onDblClick, true);
    document.removeEventListener("keydown", onKeyDown, true);
  }

  const vp = App.els["viewport"];
  vp.addEventListener("click", onClick, true);
  vp.addEventListener("dblclick", onDblClick, true);
  document.addEventListener("keydown", onKeyDown, true);
}

function selectOutlineRow(seq) {
  App.state.outlineSelectedSeq = seq;
  const tbody = App.els["outline-table"].querySelector("tbody");
  for (const tr of tbody.querySelectorAll("tr")) {
    tr.classList.toggle("outline-selected",
      String(tr.dataset.seq) === String(seq));
  }
}

async function updateClosureIndicator() {
  const indicator = App.els["outline-closure-indicator"];
  if (!indicator) return;
  try {
    const resp = await apiFetch("/api/outline/validate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ changes: {} }),
    });
    const data = await resp.json();
    if (data.valid) {
      indicator.className = "closed";
      indicator.innerHTML = '<span class="dot"></span>Closed';
    } else {
      indicator.className = "open";
      indicator.innerHTML = `<span class="dot"></span>Open (error: ${data.closure_error.toFixed(6)})`;
    }
  } catch (e) {
    indicator.className = "open";
    indicator.innerHTML = '<span class="dot"></span>Unknown';
  }
}

async function handleOutlineEdit(seq, field, rawValue) {
  const value = parseDimension(rawValue);
  if (isNaN(value)) {
    showToast(`Invalid value: ${rawValue}`, "error");
    return;
  }
  try {
    const body = {};
    body[field] = value;
    const resp = await apiFetch(`/api/outline/${seq}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    showToast(`Seg ${seq}: ${fmtFtIn(value)}`, "success");
    await loadOutlineTable();
  } catch (e) {
    showToast(`Error: ${e.message}`, "error");
    await loadOutlineTable();
  }
}

async function handleOutlineSweepEdit(seq, rawValue) {
  // Parse degrees — strip trailing degree sign
  const cleaned = rawValue.replace(/[°\u00b0]/g, "").trim();
  const deg = parseFloat(cleaned);
  if (isNaN(deg)) {
    showToast(`Invalid angle: ${rawValue}`, "error");
    return;
  }
  const rad = deg * Math.PI / 180;
  try {
    const resp = await apiFetch(`/api/outline/${seq}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sweep: rad }),
    });
    showToast(`Seg ${seq}: ${fmtDeg(deg)}`, "success");
    await loadOutlineTable();
  } catch (e) {
    showToast(`Error: ${e.message}`, "error");
    await loadOutlineTable();
  }
}

function setupOutlineToolbar() {
  const addBtn = App.els["outline-add-btn"];
  const removeBtn = App.els["outline-remove-btn"];

  if (addBtn) {
    addBtn.addEventListener("click", async () => {
      const seq = App.state.outlineSelectedSeq;
      if (seq == null) {
        showToast("Select a segment first", "error");
        return;
      }
      const name = prompt("New point name (e.g., F9a):");
      if (!name) return;
      try {
        await apiFetch("/api/outline/add-point", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ after_seq: seq, end_name: name }),
        });
        showToast(`Added point ${name}`, "success");
        loadOutlineTable();
      } catch (e) {
        showToast(`Error: ${e.message}`, "error");
      }
    });
  }

  if (removeBtn) {
    removeBtn.addEventListener("click", async () => {
      const seq = App.state.outlineSelectedSeq;
      if (seq == null) {
        showToast("Select a segment first", "error");
        return;
      }
      if (!confirm(`Delete segment at seq ${seq}?`)) return;
      try {
        await apiFetch(`/api/outline/${seq}`, { method: "DELETE" });
        App.state.outlineSelectedSeq = null;
        showToast(`Removed segment ${seq}`, "success");
        loadOutlineTable();
      } catch (e) {
        showToast(`Error: ${e.message}`, "error");
      }
    });
  }
}


/* ========== OPENINGS TABLE ========== */

function updateOpeningsTable() {
  const g = App.state.geometry;
  if (!g) return;

  // Outer openings
  const tbody1 = App.els["openings-table"].querySelector("tbody");
  tbody1.innerHTML = "";
  for (const op of (g.outer_openings || [])) {
    const dx = op.poly[1][0] - op.poly[0][0];
    const dy = op.poly[1][1] - op.poly[0][1];
    const w = Math.sqrt(dx * dx + dy * dy);
    const tr = document.createElement("tr");
    const tdName = document.createElement("td");
    tdName.textContent = op.name;
    tr.appendChild(tdName);
    const tdSeg = document.createElement("td");
    tdSeg.textContent = `${op.seg_start}–${op.seg_end}`;
    tr.appendChild(tdSeg);
    // DT-7: Editable width cell
    const tdW = document.createElement("td");
    const wConstName = findWidthConstant(op.name);
    const wConst = wConstName && App.state.constants.find(c => c.name === wConstName);
    if (wConst) {
      const inp = document.createElement("input");
      inp.type = "text";
      inp.className = "table-edit-input";
      inp.value = formatConstValue(wConst);
      inp.addEventListener("change", () => handleConstantEdit(wConstName, inp.value));
      inp.addEventListener("click", (e) => e.stopPropagation());
      tdW.appendChild(inp);
    } else {
      tdW.textContent = fmtFtIn(w);
    }
    tr.appendChild(tdW);
    tr.classList.add("selectable");
    tr.addEventListener("click", () => selectElement("opening", op.name, op));
    tbody1.appendChild(tr);
  }

  // Rough openings
  const tbody2 = App.els["rough-openings-table"].querySelector("tbody");
  tbody2.innerHTML = "";
  for (const ro of (g.rough_openings || [])) {
    const tr = document.createElement("tr");
    const tdName = document.createElement("td");
    tdName.textContent = ro.name;
    tr.appendChild(tdName);
    const tdWall = document.createElement("td");
    tdWall.textContent = ro.wall_name;
    tr.appendChild(tdWall);
    // DT-7: Editable width cell
    const tdW = document.createElement("td");
    const wConstName = findWidthConstant(ro.name);
    const wConst = wConstName && App.state.constants.find(c => c.name === wConstName);
    if (wConst) {
      const inp = document.createElement("input");
      inp.type = "text";
      inp.className = "table-edit-input";
      inp.value = formatConstValue(wConst);
      inp.addEventListener("change", () => handleConstantEdit(wConstName, inp.value));
      inp.addEventListener("click", (e) => e.stopPropagation());
      tdW.appendChild(inp);
    } else {
      tdW.textContent = fmtFtIn(ro.width);
    }
    tr.appendChild(tdW);
    const tdOri = document.createElement("td");
    tdOri.textContent = ro.orientation;
    tr.appendChild(tdOri);
    tr.classList.add("selectable");
    tr.addEventListener("click", () => selectElement("rough_opening", ro.name, ro));
    tbody2.appendChild(tr);
  }
}


/* ========== ELEMENTS TABLE ========== */

async function loadElements() {
  try {
    const [elemResp, doorResp, exclResp] = await Promise.all([
      fetch("/api/elements"),
      fetch("/api/doors"),
      fetch("/api/exclusions"),
    ]);
    App.state.elements = await elemResp.json();
    App.state.doors = await doorResp.json();
    App.state.exclusions = await exclResp.json();
    updateElementsTable();
  } catch (e) {
    console.error("Elements load failed:", e);
  }
}

function updateElementsTable() {
  const g = App.state.geometry;
  const elements = App.state.elements || [];
  const doors = App.state.doors || [];

  // Interior walls table
  const tbody1 = App.els["interior-walls-table"]
    ? App.els["interior-walls-table"].querySelector("tbody") : null;
  if (tbody1) {
    tbody1.innerHTML = "";
    const walls = elements.filter(e => e.type === "wall");
    for (const wall of walls) {
      const tr = document.createElement("tr");
      const iwData = g ? g.interior_walls[wall.name] : null;
      let thickness = "—";
      let length = "—";
      let orientation = "—";
      let openings = "—";
      if (iwData) {
        const b = iwData.bbox;
        const ew = b.e - b.w;
        const ns = b.n - b.s;
        const thick = Math.min(ew, ns);
        const len = Math.max(ew, ns);
        thickness = (thick * 12).toFixed(1).replace(/\.0$/, '') + '"';
        length = fmtFtIn(len);
        orientation = ew > ns ? "H" : "V";
      }
      // Find hosted openings from rough_openings data
      const hosted = g ? (g.rough_openings || [])
        .filter(ro => ro.wall_name === wall.name)
        .map(ro => ro.name) : [];
      openings = hosted.length > 0 ? hosted.join(", ") : "—";

      tr.innerHTML = `
        <td>${wall.name}</td>
        <td>${thickness}</td>
        <td>${length}</td>
        <td>${orientation}</td>
        <td>${openings}</td>
      `;
      tr.classList.add("selectable");
      tr.addEventListener("click", (e) => {
        // Use geometry data if available, otherwise a minimal stub so
        // excluded walls can still be selected to access properties/delete
        const data = iwData || {bbox: {w: 0, e: 0, s: 0, n: 0}};
        selectElement("wall", wall.name, data, e);
      });
      tbody1.appendChild(tr);
    }
  }

  // Furniture & Appliances table
  const tbody2 = App.els["furniture-table"]
    ? App.els["furniture-table"].querySelector("tbody") : null;
  if (tbody2) {
    tbody2.innerHTML = "";
    const items = g ? g.variant_items || {} : {};
    for (const [name, item] of Object.entries(items)) {
      const tr = document.createElement("tr");
      const b = item.bbox;
      const w = b.e - b.w;
      const d = b.n - b.s;
      tr.innerHTML = `
        <td>${item.label || name.toUpperCase()}</td>
        <td>${item.type || "—"}</td>
        <td>${fmtFtIn(w)}</td>
        <td>${fmtFtIn(d)}</td>
      `;
      tr.classList.add("selectable");
      tr.addEventListener("click", (e) => {
        selectElement(item.type || "furniture", name, item, e);
      });
      tbody2.appendChild(tr);
    }
  }
}

function showDoorProperties(tbody, roName) {
  const doors = App.state.doors || [];
  const door = doors.find(d => d.opening_name === roName);
  if (!door) {
    // No door — show "Add Door" button
    const tr = document.createElement("tr");
    const td = document.createElement("td");
    td.colSpan = 2;
    const btn = document.createElement("button");
    btn.textContent = "Add Door";
    btn.className = "prop-btn";
    btn.addEventListener("click", () => showAddDoorDialog(roName));
    td.appendChild(btn);
    tr.appendChild(td);
    tbody.appendChild(tr);
    return;
  }
  addPropRow(tbody, "—", "Door");
  addPropRow(tbody, "Width", door.width + '"');
  addDoorDropdownRow(tbody, "Hinge", door.hinge_side, roName, "hinge_side", null, door);
  addDoorDropdownRow(tbody, "Swing", door.swing_direction, roName, "swing_direction", null, door);
  addDoorDropdownRow(tbody, "Type", door.door_type, roName, "door_type",
    ["single", "double"], door);

  // Flip buttons
  const tr = document.createElement("tr");
  const td = document.createElement("td");
  td.colSpan = 2;
  td.style.display = "flex";
  td.style.gap = "4px";
  for (const [label, field] of [["Flip Hinge", "hinge_side"], ["Flip Swing", "swing_direction"]]) {
    const btn = document.createElement("button");
    btn.textContent = label;
    btn.className = "prop-btn";
    btn.addEventListener("click", () => flipDoorProperty(roName, field, door));
    td.appendChild(btn);
  }
  tr.appendChild(td);
  tbody.appendChild(tr);
}

function addDoorDropdownRow(tbody, label, value, openingName, field, options, door) {
  const dirs = options || ["east", "west", "north", "south"];
  const tr = document.createElement("tr");
  const tdLabel = document.createElement("td");
  tdLabel.textContent = label;
  tr.appendChild(tdLabel);
  const tdVal = document.createElement("td");
  const sel = document.createElement("select");
  sel.className = "prop-select";
  for (const d of dirs) {
    const opt = document.createElement("option");
    opt.value = d;
    opt.textContent = d;
    if (d === value) opt.selected = true;
    sel.appendChild(opt);
  }
  sel.addEventListener("change", async () => {
    try {
      await apiFetch(`/api/doors/${openingName}`, {
        method: "PUT",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify({[field]: sel.value}),
      });
      if (door) door[field] = sel.value;
    } catch (e) {
      showToast("Door update failed: " + e.message, "error");
    }
  });
  tdVal.appendChild(sel);
  tr.appendChild(tdVal);
  tbody.appendChild(tr);
}

const FLIP_MAP = {
  east: "west", west: "east", north: "south", south: "north",
};

async function flipDoorProperty(openingName, field, door) {
  const cur = door[field];
  const flipped = FLIP_MAP[cur] || cur;
  try {
    await apiFetch(`/api/doors/${openingName}`, {
      method: "PUT",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({[field]: flipped}),
    });
    door[field] = flipped;
  } catch (e) {
    showToast("Flip failed: " + e.message, "error");
  }
}

function showAddDoorDialog(openingName) {
  Dialog.show({
    title: "Add Door — " + openingName,
    fields: [
      {name: "width", label: "Width (inches)", type: "number", value: "36"},
      {name: "hinge_side", label: "Hinge side", type: "select",
       options: ["east", "west", "north", "south"], value: "east"},
      {name: "swing_direction", label: "Swing direction", type: "select",
       options: ["east", "west", "north", "south"], value: "south"},
      {name: "door_type", label: "Type", type: "select",
       options: ["single", "double"], value: "single"},
    ],
    onSubmit: async (vals) => {
      try {
        await apiFetch("/api/doors", {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify({
            opening_name: openingName,
            width: parseFloat(vals.width) || 36,
            hinge_side: vals.hinge_side,
            swing_direction: vals.swing_direction,
            door_type: vals.door_type,
          }),
        });
        Dialog.close();
      } catch (e) {
        showToast("Add door failed: " + e.message, "error");
      }
    },
  });
}


/* ========== EVENT LISTENERS ========== */

function setupEventListeners() {
  const vp = App.els["viewport"];

  // Mouse events on viewport
  vp.addEventListener("mousedown", onMouseDown);
  vp.addEventListener("mousemove", onMouseMove);
  vp.addEventListener("mouseup", onMouseUp);
  vp.addEventListener("wheel", onWheel, { passive: false });
  vp.addEventListener("click", onViewportClick);

  // Context menu for dimensions
  vp.addEventListener("contextmenu", onContextMenu);

  // Double-click for inline label editing (LABEL-3)
  vp.addEventListener("dblclick", onLabelDblClick);

  // Keyboard shortcuts
  document.addEventListener("keydown", onKeyDown);

  // Tool buttons
  document.querySelectorAll(".tool-btn").forEach(btn => {
    btn.addEventListener("click", () => setTool(btn.dataset.tool));
  });

  // Display toggles
  function rerender() {
    if (App.state.variant === "plumbing") renderPlumbingCanvas();
    else renderCanvas();
  }
  // Data-driven display toggle listeners
  const toggleMap = [
    ["show-points", "showPoints"], ["show-labels", "showLabels"],
    ["show-dims", "showDims"], ["show-user-dims", "showUserDims"],
    ["show-openings", "showOpenings"], ["show-furniture", "showFurniture"],
    ["show-rooms", "showRooms"], ["show-doors", "showDoors"],
    ["show-clearance", "showClearance"], ["open-links", "openLinks"],
    ["show-areas", "showAreas"],
  ];
  for (const [elId, stateKey] of toggleMap) {
    App.els[elId].addEventListener("change", (e) => {
      App.state[stateKey] = e.target.checked;
      rerender();
    });
  }
  App.els["show-grid"].addEventListener("change", (e) => {
    App.state.showGrid = e.target.checked;
    applyTransform();
  });

  // Roof style selector (SCAD-2)
  App.els["roof-style"].addEventListener("change", async (e) => {
    try {
      await apiFetch("/api/config/roof_style", {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ value: e.target.value }),
      });
    } catch (err) { console.error("Failed to save roof style:", err); }
  });

  // Variant selector
  App.els["variant-select"].addEventListener("change", async (e) => {
    await saveCurrentLayerConfig();
    App.state.variant = e.target.value;
    restoreLayerConfig();
    loadElements();
    loadGeometry();
    updateDeleteVariantBtn();
    // Show/hide plumbing tools based on variant
    if (App.els["plumbing-tools"]) {
      const isPlumbing = App.state.variant === "plumbing" &&
                         (App.state.activeView === "interactive" || App.state.activeView === "floorplan");
      App.els["plumbing-tools"].style.display = isPlumbing ? "" : "none";
    }
    if (App.state.activeView === "floorplan") loadSVGView("floorplan");
  });

  // New/Delete variant buttons
  document.getElementById("new-variant-btn").addEventListener("click", createNewVariant);
  document.getElementById("delete-variant-btn").addEventListener("click", deleteActiveVariant);

  // Constants filters
  App.els["const-category-filter"].addEventListener("change", renderConstantsTable);
  App.els["const-search"].addEventListener("input", renderConstantsTable);

  // Sort headers
  document.querySelectorAll("#constants-table th.sortable").forEach(th => {
    th.addEventListener("click", () => {
      const key = th.dataset.sort;
      if (App.state.constSortKey === key) {
        App.state.constSortAsc = !App.state.constSortAsc;
      } else {
        App.state.constSortKey = key;
        App.state.constSortAsc = true;
      }
      renderConstantsTable();
    });
  });

  // Panel tabs
  document.querySelectorAll(".panel-tab").forEach(tab => {
    tab.addEventListener("click", () => {
      document.querySelectorAll(".panel-tab").forEach(t => t.classList.remove("active"));
      document.querySelectorAll(".panel-content").forEach(p => p.classList.remove("active"));
      tab.classList.add("active");
      const panel = document.getElementById(`panel-${tab.dataset.panel}`);
      if (panel) panel.classList.add("active");
      // Load data for panel if needed
      if (tab.dataset.panel === "outline") loadOutlineTable();
      if (tab.dataset.panel === "elements") updateElementsTable();
      if (tab.dataset.panel === "plumbing") updatePlumbingTable();
    });
  });

  // Menu actions
  document.querySelectorAll("[data-action]").forEach(btn => {
    btn.addEventListener("click", () => handleMenuAction(btn.dataset.action));
  });

  // Window resize
  window.addEventListener("resize", () => {
    if (isCanvasView() && App.state.geometry) {
      fitToWindow();
    }
  });
}


/* ========== MOUSE HANDLERS ========== */

function onMouseDown(e) {
  // SVG view: any left-click starts pan drag
  if (!isCanvasView()) {
    if (e.button === 0 || e.button === 1) {
      App.state.isDragging = true;
      App.state.dragStart = { x: e.clientX, y: e.clientY };
      App.state.lastPan = { ...App.state.svgView.pan };
      App.els["viewport"].style.cursor = "grabbing";
      e.preventDefault();
    }
    return;
  }
  if (e.button === 1 || (e.button === 0 && App.state.activeTool === "pan")) {
    // Pan mode
    App.state.isDragging = true;
    App.state.dragStart = { x: e.clientX, y: e.clientY };
    App.state.lastPan = { ...App.state.pan };
    App.els["viewport"].style.cursor = "grabbing";
    e.preventDefault();
  } else if (e.button === 0 && OpeningTool.active) {
    openingToolMouseDown(e);
  } else if (e.button === 0 && PlaceTool.active) {
    placeToolMouseDown(e);
  } else if (e.button === 0 && App.state.activeTool === "select") {
    // OE-1: outline F-point drag (only if select tool + points visible)
    if (typeof outlineEditorMouseDown === "function") {
      outlineEditorMouseDown(e);
    }
    // SEL-4: Start rubber-band drag-select if clicking empty space
    if (!OutlineEditor || !OutlineEditor.pending) {
      const rect = App.els["viewport"].getBoundingClientRect();
      const sx = e.clientX - rect.left;
      const sy = e.clientY - rect.top;
      App.state.rubberBand = { sx, sy, active: false };
    }
  } else if (e.button === 0 && App.state.activeTool === "move") {
    moveToolMouseDown(e);
  } else if (e.button === 0 && App.state.activeTool === "draw-wall") {
    drawWallMouseDown(e);
  } else if (e.button === 0 && isPlumbingDrawTool(App.state.activeTool)) {
    const [wx, wy] = mouseToWorld(e);
    // Initialise draw state if starting
    if (PlumbingDraw.points.length === 0) {
      if (App.state.activeTool === "supply-cold") {
        PlumbingDraw.type = "supply_pipe"; PlumbingDraw.hotCold = "cold";
      } else if (App.state.activeTool === "supply-hot") {
        PlumbingDraw.type = "supply_pipe"; PlumbingDraw.hotCold = "hot";
      } else {
        PlumbingDraw.type = "drain_pipe"; PlumbingDraw.hotCold = null;
      }
    }
    plumbingDrawClick(wx, wy);
  } else if (e.button === 0 && App.state.activeTool === "place-fitting") {
    const [wx, wy] = mouseToWorld(e);
    placeFitting(wx, wy);
  } else if (e.button === 0 && App.state.activeTool === "place-fixture") {
    const [wx, wy] = mouseToWorld(e);
    placeFixture(wx, wy);
  } else if (e.button === 0 && App.state.activeTool === "dimension") {
    dimToolMouseDown(e);
  } else if (e.button === 0 && App.state.activeTool === "label") {
    labelToolMouseDown(e);
  } else if (e.button === 0 && App.state.activeTool === "measure") {
    const [wx, wy] = mouseToWorld(e);
    App.state.measureStart = [wx, wy];
  }
}

function onMouseMove(e) {
  // SVG view pan drag
  if (!isCanvasView()) {
    if (App.state.isDragging) {
      const dx = e.clientX - App.state.dragStart.x;
      const dy = e.clientY - App.state.dragStart.y;
      App.state.svgView.pan.x = App.state.lastPan.x + dx;
      App.state.svgView.pan.y = App.state.lastPan.y + dy;
      svgViewApplyTransform();
    }
    return;
  }
  const rect = App.els["viewport"].getBoundingClientRect();
  const sx = e.clientX - rect.left;
  const sy = e.clientY - rect.top;
  const [wx, wy] = screenToWorld(sx, sy);

  // Update coordinate display
  App.els["coord-display"].textContent =
    `E: ${fmtFtIn(wx)}  N: ${fmtFtIn(wy)}`;

  // TL-17: Endpoint drag
  if (EndpointDragTool.active) {
    endpointDragMouseMove(e);
    return;
  }

  // Outline editor drag
  if (typeof OutlineEditor !== "undefined" && (OutlineEditor.active || OutlineEditor.pending)) {
    outlineEditorMouseMove(e);
    return;
  }

  // Move tool drag (active or pending threshold check)
  if (MoveTool.active || MoveTool.pending) {
    moveToolMouseMove(e);
    return;
  }

  // Draw wall preview
  if (App.state.activeTool === "draw-wall" && DrawWallTool.start) {
    drawWallMouseMove(e);
  }

  // Plumbing draw rubber-band
  if (isPlumbingDrawTool(App.state.activeTool) && PlumbingDraw.points.length > 0) {
    const layer = App.els["layer-plumbing-pipes"];
    const oldRB = layer.querySelector(".plumbing-rubber-band");
    if (oldRB) oldRB.remove();
    const last = PlumbingDraw.points[PlumbingDraw.points.length - 1];
    let color = PLUMBING_COLORS.drain_pipe;
    if (PlumbingDraw.type === "supply_pipe") {
      color = PlumbingDraw.hotCold === "hot" ? PLUMBING_COLORS.supply_hot : PLUMBING_COLORS.supply_cold;
    }
    const rb = svgEl("line", {
      x1: last[0], y1: -last[1], x2: wx, y2: -wy,
      stroke: color, "stroke-width": "0.04",
      "stroke-dasharray": "0.1 0.05",
      opacity: "0.5",
      class: "plumbing-rubber-band",
    });
    layer.appendChild(rb);
  }

  // Dimension tool preview + snap indicator
  if (App.state.activeTool === "dimension") {
    dimToolMouseMove(e);
  }

  // SEL-4: Rubber-band drag-select
  if (App.state.rubberBand) {
    const rb = App.state.rubberBand;
    const dx = sx - rb.sx;
    const dy = sy - rb.sy;
    if (!rb.active && Math.abs(dx) + Math.abs(dy) > 4) {
      rb.active = true;
      // Create rubber-band SVG rect in screen coords projected to world
      const el = svgEl("rect", { class: "rubber-band" });
      document.getElementById("layer-measure").appendChild(el);
      rb.el = el;
    }
    if (rb.active && rb.el) {
      // Convert screen corners to world SVG coords
      const [wx1, wy1] = screenToWorld(Math.min(rb.sx, sx), Math.min(rb.sy, sy));
      const [wx2, wy2] = screenToWorld(Math.max(rb.sx, sx), Math.max(rb.sy, sy));
      // SVG uses (E, -N): x = world E, y = -world N
      const svgX = Math.min(wx1, wx2);
      const svgY = Math.min(-wy1, -wy2);
      const svgW = Math.abs(wx2 - wx1);
      const svgH = Math.abs(wy2 - wy1);
      rb.el.setAttribute("x", svgX);
      rb.el.setAttribute("y", svgY);
      rb.el.setAttribute("width", svgW);
      rb.el.setAttribute("height", svgH);
    }
  }

  if (App.state.isDragging) {
    const dx = e.clientX - App.state.dragStart.x;
    const dy = e.clientY - App.state.dragStart.y;
    App.state.pan.x = App.state.lastPan.x + dx;
    App.state.pan.y = App.state.lastPan.y + dy;
    applyTransform();
  }

  if (App.state.activeTool === "measure" && App.state.measureStart) {
    drawMeasureLine(App.state.measureStart, [wx, wy]);
  }
}

function onMouseUp(e) {
  // TL-17: Endpoint drag end
  if (EndpointDragTool.active) {
    endpointDragMouseUp(e);
    return;
  }

  // Outline editor drag end
  if (typeof OutlineEditor !== "undefined" && (OutlineEditor.active || OutlineEditor.pending)) {
    outlineEditorMouseUp(e);
    return;
  }

  // Move tool drag end (or pending click release)
  if (MoveTool.active || MoveTool.pending) {
    moveToolMouseUp(e);
    return;
  }

  // SEL-4: Complete rubber-band drag-select
  if (App.state.rubberBand) {
    const rb = App.state.rubberBand;
    if (rb.active) {
      // Compute world-coordinate bbox
      const rect = App.els["viewport"].getBoundingClientRect();
      const sx = e.clientX - rect.left;
      const sy = e.clientY - rect.top;
      const [wx1, wy1] = screenToWorld(Math.min(rb.sx, sx), Math.min(rb.sy, sy));
      const [wx2, wy2] = screenToWorld(Math.max(rb.sx, sx), Math.max(rb.sy, sy));
      const minE = Math.min(wx1, wx2), maxE = Math.max(wx1, wx2);
      const minN = Math.min(wy1, wy2), maxN = Math.max(wy1, wy2);

      // Find all elements that intersect the rubber-band bbox
      clearSelection();
      const hits = document.querySelectorAll("[data-name][data-type]");
      for (const el of hits) {
        const elName = el.getAttribute("data-name");
        const elType = el.getAttribute("data-type");
        // Use SVG bounding box for hit test
        try {
          const bb = el.getBBox();
          // SVG coords: x = E, y = -N
          const elMinE = bb.x, elMaxE = bb.x + bb.width;
          const elMinN = -bb.y - bb.height, elMaxN = -bb.y;
          // Check overlap
          if (elMaxE >= minE && elMinE <= maxE &&
              elMaxN >= minN && elMinN <= maxN) {
            App.state.selections.push({ type: elType, name: elName, data: {} });
            el.classList.add("multi-selected");
          }
        } catch (_) { /* ignore elements without bbox */ }
      }
      if (App.state.selections.length > 0) {
        App.state.selection = App.state.selections[App.state.selections.length - 1];
        App.els["selection-info"].textContent =
          `${App.state.selections.length} selected`;
      }
      // Remove rubber-band rect
      if (rb.el && rb.el.parentNode) rb.el.remove();
    }
    App.state.rubberBand = null;
  }

  if (App.state.isDragging) {
    App.state.isDragging = false;
    if (!isCanvasView()) {
      App.els["viewport"].style.cursor = "grab";
    } else {
      App.els["viewport"].style.cursor = App.state.activeTool === "pan" ? "grab" : "crosshair";
    }
  }

  if (isCanvasView() &&
      App.state.activeTool === "measure" && App.state.measureStart) {
    const [wx, wy] = mouseToWorld(e);
    const [sx, sy] = App.state.measureStart;
    const dx = wx - sx;
    const dy = wy - sy;
    const dist = Math.sqrt(dx * dx + dy * dy);
    App.els["measure-info"].textContent = `Dist: ${fmtFtIn(dist)}`;
    App.state.measureStart = null;
  }
}

function onWheel(e) {
  e.preventDefault();
  const rect = App.els["viewport"].getBoundingClientRect();
  const mx = e.clientX - rect.left;
  const my = e.clientY - rect.top;
  const factor = e.deltaY > 0 ? 0.9 : 1.1;

  if (!isCanvasView()) {
    // SVG view zoom toward cursor
    const sv = App.state.svgView;
    sv.pan.x = mx - (mx - sv.pan.x) * factor;
    sv.pan.y = my - (my - sv.pan.y) * factor;
    sv.zoom *= factor;
    svgViewApplyTransform();
    return;
  }

  const newZoom = App.state.zoom * factor;
  // Zoom toward cursor
  App.state.pan.x = mx - (mx - App.state.pan.x) * factor;
  App.state.pan.y = my - (my - App.state.pan.y) * factor;
  App.state.zoom = newZoom;

  applyTransform();
}

function onViewportClick(e) {
  // After a move-tool drag, suppress this click
  if (MoveTool._suppressClick) {
    MoveTool._suppressClick = false;
    return;
  }
  // If clicking on the viewport background (not on an element), clear selection
  if (e.target === App.els["canvas"] || e.target === App.els["viewport"]) {
    clearSelection();
  }
}

function onLabelDblClick(e) {
  // Finish plumbing drawing on double-click
  if (isPlumbingDrawTool(App.state.activeTool) && PlumbingDraw.points.length >= 2) {
    finishPlumbingDraw();
    return;
  }
  const target = e.target.closest("[data-type='label']");
  if (!target) return;
  const name = target.getAttribute("data-name");
  const elemRec = (App.state.elements || []).find(el => el.name === name && el.type === "label");
  if (!elemRec) return;
  const props = parseProps(elemRec);

  Dialog.show({
    title: `Edit Label: ${name}`,
    fields: [{ label: "Text", name: "text", value: props.text || name }],
    async onSubmit(vals) {
      props.text = vals.text;
      await fetch(`/api/elements/${elemRec.id}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ properties: props }),
      });
      Dialog.close();
      await reloadAfterChange();
    },
  });
}

function onContextMenu(e) {
  // Context menu for dimensions
  const dimTarget = e.target.closest("[data-type='dimension']");
  if (dimTarget) {
    e.preventDefault();
    const name = dimTarget.getAttribute("data-name");
    const items = [
      { label: "Horizontal", action: () => setDimRotation(name, "horizontal") },
      { label: "Vertical", action: () => setDimRotation(name, "vertical") },
      { label: "Parallel", action: () => setDimRotation(name, "parallel") },
      { label: "Perpendicular", action: () => setDimRotation(name, "perpendicular") },
    ];
    const elemRec = (App.state.elements || []).find(el => el.name === name && el.type === "dimension");
    if (elemRec) {
      const props = parseProps(elemRec);
      const hasStart = !!props.start_anchor;
      const hasEnd = !!props.end_anchor;
      if (hasStart || hasEnd) {
        items.push({ label: "---" });
        if (hasStart) items.push({ label: "Detach Start", action: () => detachAnchor(name, "start") });
        if (hasEnd) items.push({ label: "Detach End", action: () => detachAnchor(name, "end") });
        if (hasStart && hasEnd) items.push({ label: "Detach Both", action: () => detachAnchor(name, "both") });
      }
    }
    showContextMenu(e.clientX, e.clientY, items);
    return;
  }

  // Context menu for furniture/appliance/fixture items
  const itemTarget = e.target.closest("[data-type='furniture'], [data-type='appliance'], [data-type='fixture']");
  if (itemTarget) {
    e.preventDefault();
    const name = itemTarget.getAttribute("data-name");
    const type = itemTarget.getAttribute("data-type");
    const elemRec = (App.state.elements || []).find(el => el.name === name);
    // Select the item
    const geom = App.state.geometry;
    const vi = geom ? (geom.variant_items || {}) : {};
    const itemData = vi[name];
    if (itemData) selectElement(type, name, itemData);

    const menuItems = [
      { label: `Rotate...`, action: () => showRotationDialog() },
    ];
    showContextMenu(e.clientX, e.clientY, menuItems);
    return;
  }
}

function showContextMenu(x, y, items) {
  hideContextMenu();
  const menu = document.createElement("div");
  menu.className = "context-menu";
  menu.style.left = x + "px";
  menu.style.top = y + "px";
  for (const item of items) {
    if (item.label === "---") {
      const hr = document.createElement("hr");
      hr.style.margin = "4px 0"; hr.style.border = "none"; hr.style.borderTop = "1px solid var(--border)";
      menu.appendChild(hr);
      continue;
    }
    const btn = document.createElement("div");
    btn.className = "context-menu-item";
    btn.textContent = item.label;
    btn.addEventListener("click", () => { hideContextMenu(); item.action(); });
    menu.appendChild(btn);
  }
  document.body.appendChild(menu);
  // Close on click outside
  setTimeout(() => document.addEventListener("click", hideContextMenu, { once: true }), 0);
}

function hideContextMenu() {
  const existing = document.querySelector(".context-menu");
  if (existing) existing.remove();
}

async function setDimRotation(name, rotation) {
  const elemRec = (App.state.elements || []).find(e => e.name === name && e.type === "dimension");
  if (!elemRec) return;
  const props = parseProps(elemRec);
  props.label_rotation = rotation;
  await fetch(`/api/elements/${elemRec.id}`, {
    method: "PUT",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ properties: props }),
  });
  await reloadAfterChange();
}

async function detachAnchor(name, which) {
  const elemRec = (App.state.elements || []).find(e => e.name === name && e.type === "dimension");
  if (!elemRec) return;
  const props = parseProps(elemRec);
  if (which === "start" || which === "both") delete props.start_anchor;
  if (which === "end" || which === "both") delete props.end_anchor;
  await fetch(`/api/elements/${elemRec.id}`, {
    method: "PUT",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ properties: props }),
  });
  await reloadAfterChange();
}

function onKeyDown(e) {
  if (e.target.tagName === "INPUT" || e.target.tagName === "SELECT") return;

  // Config Save / Save As / Load
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "s") {
    e.preventDefault();
    if (e.shiftKey) handleMenuAction("config-save-as");
    else handleMenuAction("config-save");
    return;
  }
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "o") {
    e.preventDefault();
    handleMenuAction("config-load");
    return;
  }

  // Undo / Redo
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "z") {
    e.preventDefault();
    if (e.shiftKey) { doRedo(); } else { doUndo(); }
    return;
  }

  switch (e.key) {
    case "v": case "V": setTool("select"); break;
    case "h": case "H": setTool("pan"); break;
    case "m": case "M": setTool("measure"); break;
    case "g": case "G": setTool("move"); break;
    case "w": case "W": setTool("draw-wall"); break;
    case "d": setTool("dimension"); break;
    case "l": setTool("label"); break;
    case "r": case "R": showRotationDialog(); break;
    case "f": case "F": fitToWindow(); break;
    case "Delete": case "Backspace":
      e.preventDefault();
      deleteSelectedElements();
      break;
    case "Enter":
      if (PlumbingDraw.points.length >= 2) {
        finishPlumbingDraw();
        break;
      }
      if (App.state.activeTool === "move" && App.state.selection) {
        showOffsetDialog();
      }
      break;
    case "Escape":
      if (PlumbingDraw.points.length > 0) {
        cancelPlumbingDraw();
        break;
      }
      if (MoveTool.active) {
        // Cancel drag: remove ghosts, reset state
        for (const g of MoveTool.origTransforms) {
          if (g.ghost && g.ghost.parentNode) g.ghost.remove();
        }
        MoveTool.active = false;
        MoveTool.origTransforms = [];
        break;
      }
      if (DrawWallTool.start) {
        cancelDrawWall();
        break;
      }
      if (DimTool.start) {
        cancelDimTool();
        break;
      }
      clearSelection();
      clearMeasure();
      break;
    case "+": case "=":
      App.state.zoom *= 1.2;
      applyTransform();
      break;
    case "-":
      App.state.zoom *= 0.8;
      applyTransform();
      break;
  }
}

/** Delete selected elements after confirmation (TL-22, TL-23). */
async function deleteSelectedElements() {
  // Gather targets from multi-selection or single selection
  const targets = App.state.selections.length > 0
    ? App.state.selections
    : (App.state.selection ? [App.state.selection] : []);
  if (targets.length === 0) return;

  // Find DB element records for each target
  const elements = App.state.elements || [];
  const toDelete = [];
  for (const t of targets) {
    const el = elements.find(e => e.name === t.name);
    if (el) toDelete.push(el);
  }
  if (toDelete.length === 0) {
    showToast("Selected items cannot be deleted", "warning");
    return;
  }

  // Check if any targets are catalog items (furniture/appliance/fixture)
  const catalogTypes = new Set(["furniture", "appliance", "fixture"]);
  const hasCatalogItems = toDelete.some(el => catalogTypes.has(el.type));

  // Build confirmation message with cascade warnings
  const lines = toDelete.map(el => {
    let line = el.name;
    if (el.type === "wall") {
      const hosted = IW_HOSTED_OPENINGS[el.name] || [];
      if (hosted.length > 0) {
        line += ` (will also delete ${hosted.join(", ")})`;
      }
    }
    return line;
  });

  let action;
  if (hasCatalogItems) {
    action = await showDeleteDialog(
      `Delete ${lines.join(", ")}?`,
      null,
      [
        { label: "Delete", action: "delete", className: "dialog-btn-primary" },
        { label: "Also delete from ADD menu", action: "delete_catalog", className: "dialog-btn-danger" },
        { label: "Cancel", action: null, className: "dialog-btn-cancel" },
      ]
    );
  } else {
    action = confirm(`Delete ${lines.join(", ")}?`) ? "delete" : null;
  }
  if (!action) return;

  // Delete each element via API
  for (const el of toDelete) {
    const resp = await fetch(`/api/elements/${el.id}`, { method: "DELETE" });
    if (!resp.ok) {
      const data = await resp.json();
      showToast(data.error || "Delete failed", "error");
      return;
    }
  }

  // If user chose to also delete from catalog, remove catalog entries
  if (action === "delete_catalog") {
    for (const el of toDelete) {
      if (catalogTypes.has(el.type)) {
        const catalogKey = el.name;
        await fetch(`/api/catalog/${encodeURIComponent(catalogKey)}`, { method: "DELETE" });
      }
    }
  }

  clearSelection();
  await reloadAfterChange();
  let toast = `Deleted ${toDelete.map(e => e.name).join(", ")}`;
  if (action === "delete_catalog") toast += " and removed from ADD menu";
  showToast(toast, "success");
}

/** Compute a rotated rectangle polygon from center, width, depth, angle (degrees). */
function rotatedRectPoly(cx, cy, w, d, angleDeg) {
  const rad = angleDeg * Math.PI / 180;
  const cos = Math.cos(rad), sin = Math.sin(rad);
  const hw = w / 2, hd = d / 2;
  // Local corners before rotation
  const corners = [[-hw, hd], [hw, hd], [hw, -hd], [-hw, -hd]];
  return corners.map(([lx, ly]) => [
    cx + lx * cos - ly * sin,
    cy + lx * sin + ly * cos,
  ]);
}

/* ========== TL-17: Endpoint Drag Handles ========== */

const EndpointDragTool = {
  active: false,
  elementId: null,
  props: null,
  handle: null,  // "start" or "end"
  origWorld: null,
};

function clearWallHandles() {
  document.querySelectorAll(".wall-handle").forEach(el => el.remove());
}

function renderWallHandles(elemId, props) {
  clearWallHandles();
  if (!props.start || !props.end) return;
  const layer = App.els["layer-measure"];
  const hs = 0.15; // handle size in world units
  for (const which of ["start", "end"]) {
    const pt = props[which];
    const rect = svgEl("rect", {
      x: pt[0] - hs / 2, y: -pt[1] - hs / 2,
      width: hs, height: hs,
      class: "wall-handle",
      "data-handle": which,
    });
    rect.addEventListener("mousedown", (e) => {
      e.stopPropagation();
      e.preventDefault();
      EndpointDragTool.active = true;
      EndpointDragTool.elementId = elemId;
      EndpointDragTool.props = { ...props };
      EndpointDragTool.handle = which;
      EndpointDragTool.origWorld = [...pt];
    });
    layer.appendChild(rect);
  }
}

function endpointDragMouseMove(e) {
  let [wx, wy] = mouseToWorld(e);
  // Shift-constrain to horizontal or vertical
  if (e.shiftKey && EndpointDragTool.origWorld) {
    const ox = EndpointDragTool.origWorld[0], oy = EndpointDragTool.origWorld[1];
    if (Math.abs(wx - ox) > Math.abs(wy - oy)) wy = oy;
    else wx = ox;
  }
  // Grid snap
  if (App.state.showGrid) {
    const gs = 1.0 / 12;
    wx = Math.round(wx / gs) * gs;
    wy = Math.round(wy / gs) * gs;
  }
  // Move the handle visually
  const hs = 0.15;
  const handleEl = document.querySelector(`.wall-handle[data-handle="${EndpointDragTool.handle}"]`);
  if (handleEl) {
    handleEl.setAttribute("x", wx - hs / 2);
    handleEl.setAttribute("y", -wy - hs / 2);
  }
  // Preview the wall outline
  const p = { ...EndpointDragTool.props };
  p[EndpointDragTool.handle] = [wx, wy];
  const poly = wallPoly(p.start, p.end, p.thickness || 4.0 / 12);
  if (poly) {
    let preview = document.querySelector(".wall-drag-preview");
    if (!preview) {
      preview = svgEl("polygon", { class: "wall-drag-preview", style: "fill:none;stroke:var(--accent);stroke-width:0.02;stroke-dasharray:0.06 0.03;opacity:0.6" });
      App.els["layer-measure"].appendChild(preview);
    }
    preview.setAttribute("points", poly.map(([x, y]) => `${x},${-y}`).join(" "));
  }
}

async function endpointDragMouseUp(e) {
  if (!EndpointDragTool.active) return;
  EndpointDragTool.active = false;
  // Remove preview
  const preview = document.querySelector(".wall-drag-preview");
  if (preview) preview.remove();
  // Compute final position
  let [wx, wy] = mouseToWorld(e);
  if (e.shiftKey && EndpointDragTool.origWorld) {
    const ox = EndpointDragTool.origWorld[0], oy = EndpointDragTool.origWorld[1];
    if (Math.abs(wx - ox) > Math.abs(wy - oy)) wy = oy;
    else wx = ox;
  }
  if (App.state.showGrid) {
    const gs = 1.0 / 12;
    wx = Math.round(wx / gs) * gs;
    wy = Math.round(wy / gs) * gs;
  }
  const p = { ...EndpointDragTool.props };
  p[EndpointDragTool.handle] = [wx, wy];
  p.poly = wallPoly(p.start, p.end, p.thickness || 4.0 / 12);
  if (!p.poly) return;
  try {
    const resp = await apiFetch(`/api/elements/${EndpointDragTool.elementId}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ properties: p }),
    });
    if (!resp.ok) throw new Error((await resp.json()).error);
    await reloadAfterChange();
  } catch (err) {
    showToast(`Error: ${err.message}`, "error");
  }
}

/* ========== TL-21: Add Opening Tool ========== */

const OpeningTool = { active: false, defaultWidth: 32.0 / 12 };

function startOpeningPlacement() {
  OpeningTool.active = true;
  App.els["viewport"].style.cursor = "crosshair";
  showToast("Click on a wall to place an opening", "info");
}

function cancelOpeningPlacement() {
  OpeningTool.active = false;
}

/** Distance from point (px, py) to line segment (x1,y1)-(x2,y2). */
function distToSeg(px, py, x1, y1, x2, y2) {
  const dx = x2 - x1, dy = y2 - y1;
  const len2 = dx * dx + dy * dy;
  if (len2 < 1e-12) return Math.sqrt((px - x1) ** 2 + (py - y1) ** 2);
  let t = ((px - x1) * dx + (py - y1) * dy) / len2;
  t = Math.max(0, Math.min(1, t));
  const cx = x1 + t * dx, cy = y1 + t * dy;
  return Math.sqrt((px - cx) ** 2 + (py - cy) ** 2);
}

/** Find nearest wall polygon edge to (wx, wy). Returns { name, type } or null. */
function findNearestWall(wx, wy) {
  let bestDist = Infinity, bestWall = null;
  const g = App.state.geometry;
  if (!g) return null;
  // Interior walls
  for (const [name, wall] of Object.entries(g.interior_walls || {})) {
    const poly = wall.poly;
    if (!poly || poly.length < 2) continue;
    for (let i = 0; i < poly.length; i++) {
      const j = (i + 1) % poly.length;
      const d = distToSeg(wx, wy, poly[i][0], poly[i][1], poly[j][0], poly[j][1]);
      if (d < bestDist) { bestDist = d; bestWall = { name, type: "iw" }; }
    }
  }
  // Drawn walls from elements
  for (const elem of (App.state.elements || [])) {
    const props = parseProps(elem);
    if (elem.type === "wall" && props?.source === "drawn" && props.poly) {
      for (let i = 0; i < props.poly.length; i++) {
        const j = (i + 1) % props.poly.length;
        const d = distToSeg(wx, wy, props.poly[i][0], props.poly[i][1], props.poly[j][0], props.poly[j][1]);
        if (d < bestDist) { bestDist = d; bestWall = { name: elem.name, type: "drawn" }; }
      }
    }
  }
  return bestDist < 0.5 ? bestWall : null;
}

function nextOpeningName() {
  const elements = App.state.elements || [];
  let max = 0;
  for (const e of elements) {
    if (e.type === "opening" && e.name && e.name.startsWith("UO")) {
      const n = parseInt(e.name.slice(2), 10);
      if (n > max) max = n;
    }
  }
  return `UO${max + 1}`;
}

function openingToolMouseDown(e) {
  const [wx, wy] = mouseToWorld(e);
  const wall = findNearestWall(wx, wy);
  if (!wall) {
    showToast("No wall found near click", "warning");
    return;
  }
  const defName = nextOpeningName();
  const defWidth = fmtFtIn(OpeningTool.defaultWidth);
  Dialog.show({
    title: "Add Opening",
    fields: [
      { label: "Name", name: "name", value: defName },
      { label: "Width", name: "width", value: defWidth },
    ],
    async onSubmit(values) {
      const width = parseDimension(values.width);
      if (!width || width <= 0) { showToast("Invalid width", "error"); return; }
      const name = values.name.trim() || defName;
      try {
        const resp = await apiFetch("/api/openings", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            name,
            segment: wall.name,
            width,
            offset: 0,
            host_wall: wall.name,
            variant: App.state.variant,
          }),
        });
        if (!resp.ok) throw new Error((await resp.json()).error);
        showToast(`Created opening ${name} on ${wall.name}`, "success");
        OpeningTool.active = false;
        App.els["viewport"].style.cursor = "crosshair";
        await reloadAfterChange();
      } catch (err) {
        showToast(`Error: ${err.message}`, "error");
      }
    },
  });
}

/** Show rotation dialog for selected element (TL-24, unified). */
function showRotationDialog() {
  const sel = App.state.selection;
  if (!sel) {
    showToast("Select an element first", "warning");
    return;
  }
  const elements = App.state.elements || [];
  const elemRec = elements.find(e => e.name === sel.name);
  if (!elemRec) {
    showToast("No editable element selected", "warning");
    return;
  }
  const props = typeof elemRec.properties === "string"
    ? JSON.parse(elemRec.properties) : elemRec.properties;
  const elemType = elemRec.type;
  // Allow rotation for furniture, appliance, fixture, and drawn walls
  if (!["furniture", "appliance", "fixture"].includes(elemType) &&
      !(props && props.source === "drawn")) {
    showToast("Rotation is available for furniture, appliances, and fixtures", "warning");
    return;
  }

  // Get current rotation from geometry if available, fall back to props
  const geom = App.state.geometry;
  const vi = geom ? (geom.variant_items || {}) : {};
  const itemGeom = vi[sel.name] || {};
  const currentAngle = itemGeom.rotation || (props && props.rotation) || 0;

  Dialog.show({
    title: `Rotate ${sel.name}`,
    fields: [
      { label: "Angle (degrees)", name: "angle", value: String(Math.round(currentAngle * 10) / 10) },
    ],
    presetButtons: {
      target: "angle",
      values: [
        { label: "0", value: "0" },
        { label: "90", value: "90" },
        { label: "180", value: "180" },
        { label: "270", value: "270" },
      ],
    },
    async onSubmit(values) {
      const angle = parseFloat(values.angle);
      if (isNaN(angle)) {
        showToast("Invalid angle", "error");
        return;
      }
      // For drawn walls, use legacy property update
      if (props && props.source === "drawn") {
        const newProps = { ...props, rotation: angle };
        await fetch(`/api/elements/${elemRec.id}`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ properties: newProps }),
        });
      } else {
        // Use formula update endpoint for all items
        await fetch(`/api/elements/${elemRec.id}/update-formula`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ rotation_deg: angle, world_rotation: angle }),
        });
      }
      await reloadAfterChange();
      showToast(`Rotated ${sel.name} to ${angle}°`, "success");
    },
  });
}

function updateUndoButtons(canUndo, canRedo) {
  const undoBtn = document.querySelector('[data-action="undo"]');
  const redoBtn = document.querySelector('[data-action="redo"]');
  if (undoBtn) undoBtn.disabled = !canUndo;
  if (redoBtn) redoBtn.disabled = !canRedo;
}

async function doUndoRedo(endpoint, label) {
  const resp = await fetch(`/api/${endpoint}`, { method: "POST" });
  const data = await resp.json();
  if (resp.ok) {
    showToast(`${label}: ` + (data.description || data.action), "success");
    updateUndoButtons(data.can_undo, data.can_redo);
    loadConstants();
    loadGeometry();
    loadElements();
    loadPlumbingElements();
    if (["variant_update", "variant_create", "variant_delete"].includes(data.action)) {
      await loadVariants();
      ensureActiveVariantValid();
    }
  } else {
    showToast(data.error || `Nothing to ${endpoint}`, "warning");
  }
}

async function doUndo() { return doUndoRedo("undo", "Undo"); }
async function doRedo() { return doUndoRedo("redo", "Redo"); }


/* ========== TOOLS ========== */

function setTool(tool) {
  App.state.activeTool = tool;
  document.querySelectorAll(".tool-btn").forEach(btn => {
    btn.classList.toggle("active", btn.dataset.tool === tool);
  });
  const cursors = {
    select: "crosshair", pan: "grab", measure: "crosshair",
    move: "move", "draw-wall": "crosshair", dimension: "crosshair", label: "crosshair",
    "supply-cold": "crosshair", "supply-hot": "crosshair", drain: "crosshair",
    "place-fitting": "crosshair", "place-fixture": "crosshair",
  };
  App.els["viewport"].style.cursor = cursors[tool] || "crosshair";

  if (tool !== "measure") clearMeasure();
  if (tool !== "draw-wall") cancelDrawWall();
  if (tool !== "dimension" && typeof cancelDimTool === "function") cancelDimTool();
  if (!isPlumbingDrawTool(tool)) cancelPlumbingDraw();
  cancelOpeningPlacement();
}


/* ========== DIMENSION LINES ========== */

function renderUserDimensions(g) {
  if (!g.user_dimensions) return;
  if (!App.state.showDims && !App.state.showUserDims) return;
  const layer = App.els["layer-labels"];

  for (const ud of g.user_dimensions) {
    const p = ud.properties;
    if (!p.start || !p.end) continue;
    const isBuiltin = p.source === "builtin";
    // DIS-7: gate on appropriate toggle
    if (isBuiltin && !App.state.showDims) continue;
    if (!isBuiltin && !App.state.showUserDims) continue;
    const offset = p.offset || 0;
    const ax = p.start[0], ay = p.start[1];
    const bx = p.end[0], by = p.end[1];
    const dimStyle = p.dim_style || (isBuiltin ? "solid" : "dashed");
    const lineCls = dimStyle === "solid" ? "dim-line" : "user-dim-line";
    const labelCls = dimStyle === "solid" ? "dim-label" : "user-dim-label";

    // Direction and perpendicular
    const dx = bx - ax, dy = by - ay;
    const len = Math.sqrt(dx * dx + dy * dy);
    if (len < 0.01) continue;
    const ux = dx / len, uy = dy / len;
    const px = -uy, py = ux;

    // Offset endpoints
    const ax2 = ax + offset * px, ay2 = ay + offset * py;
    const bx2 = bx + offset * px, by2 = by + offset * py;

    // Extension lines (A→A', B→B')
    if (Math.abs(offset) > 0.01) {
      layer.appendChild(svgEl("line", {
        x1: ax, y1: -ay, x2: ax2, y2: -ay2, class: "user-dim-ext",
      }));
      layer.appendChild(svgEl("line", {
        x1: bx, y1: -by, x2: bx2, y2: -by2, class: "user-dim-ext",
      }));
    }

    // Main dimension line (A'→B')
    const mainLine = svgEl("line", {
      x1: ax2, y1: -ay2, x2: bx2, y2: -by2,
      class: `${lineCls} selectable`,
      "data-type": "dimension", "data-name": ud.name,
    });
    // Apply inline style overrides (STYLE-1..3)
    if (p.stroke_color) mainLine.style.stroke = p.stroke_color;
    if (p.stroke_style === "dashed") mainLine.style.strokeDasharray = "0.06 0.03";
    else if (p.stroke_style === "dotted") mainLine.style.strokeDasharray = "0.02 0.02";
    if (p.opacity != null && p.opacity !== 100) mainLine.style.opacity = p.opacity / 100;
    mainLine.addEventListener("click", (e) => selectElement("dimension", ud.name, ud, e));
    layer.appendChild(mainLine);

    // Tick marks at endpoints
    const tw = 0.15;
    for (const [tx, ty] of [[ax2, -ay2], [bx2, -by2]]) {
      layer.appendChild(svgEl("line", {
        x1: tx - px * tw, y1: ty + py * tw,
        x2: tx + px * tw, y2: ty - py * tw,
        class: lineCls,
      }));
    }

    // Label at midpoint
    const mx = (ax2 + bx2) / 2, my = (-ay2 + -by2) / 2;
    const dist = Math.sqrt((bx2 - ax2) ** 2 + (by2 - ay2) ** 2);
    let ang = Math.atan2(-by2 - (-ay2), bx2 - ax2) * 180 / Math.PI;
    // Determine rotation based on label_rotation property
    const rot = p.label_rotation || "parallel";
    if (rot === "horizontal") ang = 0;
    else if (rot === "vertical") ang = -90;
    else if (rot === "perpendicular") ang += 90;
    // Keep readable (not upside-down)
    if (ang >= 90) ang -= 180;
    else if (ang < -90) ang += 180;

    const angRad = ang * Math.PI / 180;
    const lx = mx + 0.15 * Math.sin(angRad);
    const ly = my - 0.15 * Math.cos(angRad);

    const label = svgEl("text", {
      x: lx, y: ly, class: `${labelCls} selectable`,
      "data-type": "dimension", "data-name": ud.name,
      transform: `rotate(${ang},${lx},${ly})`,
    });
    if (p.stroke_color) label.style.fill = p.stroke_color;
    if (p.opacity != null && p.opacity !== 100) label.style.opacity = p.opacity / 100;
    label.textContent = fmtFtIn(dist);
    label.addEventListener("click", (e) => selectElement("dimension", ud.name, ud, e));
    layer.appendChild(label);
  }
}

function renderUserLabels(g) {
  if (!g.label_elements) return;
  const layer = App.els["layer-labels"];

  for (const le of g.label_elements) {
    const p = le.properties;
    if (p.source === "room") continue; // rendered by renderRoomLabels
    if (!p.position) continue;

    const [ex, en] = p.position;
    const fontSize = p.font_size || 0.25;
    const rotation = p.rotation || 0;
    const text = p.text || le.name;

    const el = svgEl("text", {
      x: ex, y: -en,
      class: "user-label selectable",
      "text-anchor": "middle",
      "dominant-baseline": "middle",
      "data-type": "label", "data-name": le.name,
    });
    el.style.fontSize = fontSize + "px";
    if (p.fill_color) el.style.fill = p.fill_color;
    if (p.opacity != null && p.opacity !== 100) el.style.opacity = p.opacity / 100;
    if (rotation) el.setAttribute("transform", `rotate(${-rotation},${ex},${-en})`);
    el.textContent = text;
    el.addEventListener("click", (e) => selectElement("label", le.name, le, e));
    layer.appendChild(el);
  }
}


/* ========== MEASURE ========== */

function drawMeasureLine(p1, p2) {
  const layer = App.els["layer-measure"];
  layer.innerHTML = "";

  const line = svgEl("line", {
    x1: p1[0], y1: -p1[1],
    x2: p2[0], y2: -p2[1],
    class: "measure-line",
  });
  layer.appendChild(line);

  const dx = p2[0] - p1[0];
  const dy = p2[1] - p1[1];
  const dist = Math.sqrt(dx * dx + dy * dy);
  const mx = (p1[0] + p2[0]) / 2;
  const my = -(p1[1] + p2[1]) / 2;

  const label = svgEl("text", {
    x: mx, y: my - 0.3,
    class: "measure-label",
  });
  label.textContent = fmtFtIn(dist);
  layer.appendChild(label);
}

function clearMeasure() {
  App.els["layer-measure"].innerHTML = "";
  App.els["measure-info"].textContent = "";
  App.state.measureStart = null;
}


/* ========== MENU ACTIONS ========== */

async function handleMenuAction(action) {
  switch (action) {
    case "regenerate-all": {
      showToast("Regenerating all views...");
      try {
        await apiFetch("/api/regenerate", { method: "POST" });
        showToast("All views regenerated", "success");
      } catch (e) { showToast(`Regenerate failed: ${e.message}`, "error"); }
      break;
    }
    case "regenerate-view": {
      if (isCanvasView()) {
        showToast("Cannot regenerate interactive view — edit constants instead");
        return;
      }
      showToast(`Regenerating ${App.state.activeView}...`);
      try {
        await apiFetch("/api/regenerate", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ view: App.state.activeView }),
        });
        showToast(`${App.state.activeView} regenerated`, "success");
      } catch (e) { showToast(`Regenerate failed: ${e.message}`, "error"); }
      break;
    }
    case "generate-3d": {
      showToast("Generating 3D model...");
      try {
        const resp = await apiFetch("/api/generate-3d", { method: "POST" });
        const data = await resp.json();
        if (data.ok) showToast(`3D model generated (${data.roof_style})`, "success");
        else showToast(`Generation failed: ${data.error || "unknown"}`, "error");
      } catch (e) { showToast(`Generation failed: ${e.message}`, "error"); }
      break;
    }
    case "generate-views": {
      showToast("Generating views (this may take a moment)...");
      try {
        const resp = await apiFetch("/api/generate-views", { method: "POST" });
        const data = await resp.json();
        if (data.ok) showToast("Views generated (3views.pdf)", "success");
        else showToast(`Generation failed: ${data.error || "unknown"}`, "error");
      } catch (e) { showToast(`Generation failed: ${e.message}`, "error"); }
      break;
    }
    case "generate-site-plan": {
      showToast("Generating site plan...");
      try {
        const resp = await apiFetch("/api/generate-site-plan", { method: "POST" });
        const data = await resp.json();
        if (data.ok) {
          showToast("Site plan regenerated", "success");
          if (App.state.activeView.startsWith("site_plan")) {
            loadPDFView(App.state.activeView);
          }
        } else showToast(`Generation failed: ${data.error || "unknown"}`, "error");
      } catch (e) { showToast(`Generation failed: ${e.message}`, "error"); }
      break;
    }
    case "add-drainfield": {
      try {
        const resp = await apiFetch("/api/elements", {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify({
            type: "site_element",
            name: `drainfield_${Date.now()}`,
            properties: { subtype: "drainfield", width: 25.0, height: 10.0,
                          label: "NEW DRAINFIELD", x: 0, y: 0, rotation: 0 },
          }),
        });
        if (resp.ok) showToast("Drainfield added", "success");
        else showToast("Failed to add drainfield", "error");
      } catch (e) { showToast(`Failed: ${e.message}`, "error"); }
      break;
    }
    case "add-site-note": {
      const text = prompt("Annotation text:");
      if (!text) break;
      try {
        const resp = await apiFetch("/api/elements", {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify({
            type: "site_annotation",
            name: `note_${Date.now()}`,
            properties: { text, style: "text", x: 0, y: 0, font_size: 8.0, rotation: 0 },
          }),
        });
        if (resp.ok) showToast("Annotation added", "success");
        else showToast("Failed to add annotation", "error");
      } catch (e) { showToast(`Failed: ${e.message}`, "error"); }
      break;
    }
    case "export-all": {
      showToast("Exporting all products...");
      try {
        await apiFetch("/api/regenerate", { method: "POST" });
        showToast("Export complete", "success");
      } catch (e) { showToast(`Export failed: ${e.message}`, "error"); }
      break;
    }
    case "config-save": {
      if (App.state.configName) {
        await saveConfig(App.state.configName, true);
      } else {
        await promptSaveConfigAs();
      }
      break;
    }
    case "config-save-as": {
      await promptSaveConfigAs();
      break;
    }
    case "config-load": {
      await showLoadConfigDialog();
      break;
    }
    case "reset-constants": {
      if (!confirm("Reset all constants, outline, and elements to original values?")) return;
      try {
        await apiFetch("/api/constants/reset", { method: "POST" });
        await loadConstants();
        await loadOutlineTable();
        await reloadAfterChange();
        showToast("Reset to defaults", "success");
      } catch (e) { showToast(`Reset failed: ${e.message}`, "error"); }
      break;
    }
    case "reset-database": {
      await resetDatabase();
      break;
    }
    case "zoom-fit": fitToWindow(); break;
    case "zoom-in": App.state.zoom *= 1.3; applyTransform(); break;
    case "zoom-out": App.state.zoom *= 0.7; applyTransform(); break;
    case "toggle-points":
      App.state.showPoints = !App.state.showPoints;
      App.els["show-points"].checked = App.state.showPoints;
      renderCanvas();
      break;
    case "toggle-grid":
      App.state.showGrid = !App.state.showGrid;
      App.els["show-grid"].checked = App.state.showGrid;
      applyTransform();
      break;
    case "toggle-dims":
      App.state.showDims = !App.state.showDims;
      App.els["show-dims"].checked = App.state.showDims;
      renderCanvas();
      break;
    case "tool-select": setTool("select"); break;
    case "tool-pan": setTool("pan"); break;
    case "tool-measure": setTool("measure"); break;
    case "tool-move": setTool("move"); break;
  }
}


/* ========== TOAST ========== */

function showToast(msg, type = "") {
  const existing = document.querySelector(".toast");
  if (existing) existing.remove();

  const toast = document.createElement("div");
  toast.className = `toast ${type}`;
  toast.textContent = msg;
  document.body.appendChild(toast);
  setTimeout(() => toast.remove(), 3000);
}


/* ========== ERROR BANNER ========== */

function showErrorBanner(text, actionLabel, actionFn) {
  const banner = App.els["error-banner"];
  if (!banner) return;
  App.els["error-banner-text"].textContent = text;
  const btn = App.els["error-banner-action"];
  btn.textContent = actionLabel || "Reset Database";
  btn.onclick = actionFn || resetDatabase;
  App.els["error-banner-dismiss"].onclick = () => { banner.style.display = "none"; };
  banner.style.display = "flex";
}

function hideErrorBanner() {
  const banner = App.els["error-banner"];
  if (banner) banner.style.display = "none";
}

/* ========== CONFIGURATION FILES ========== */

function markConfigDirty() {
  App.state.configDirty = true;
  updateConfigIndicator();
}

function updateConfigIndicator() {
  const el = App.els["config-name"];
  if (!el) return;
  const name = App.state.configName || "Unsaved";
  el.textContent = App.state.configDirty ? name + " *" : name;
}

async function loadConfigName() {
  try {
    const resp = await apiFetch("/api/config/config_name");
    const data = await resp.json();
    if (data.value) {
      App.state.configName = data.value;
      App.state.configDirty = false;
      updateConfigIndicator();
    }
  } catch (_) { /* no config name set yet */ }
}

async function saveConfig(name, overwrite = false) {
  try {
    const resp = await apiFetch("/api/configs/save", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name, overwrite }),
    });
    if (resp.status === 409) {
      if (confirm(`Configuration "${name}" already exists. Overwrite?`)) {
        return saveConfig(name, true);
      }
      return;
    }
    if (!resp.ok) {
      const err = await resp.json();
      showToast(`Save failed: ${err.error}`, "error");
      return;
    }
    App.state.configName = name;
    App.state.configDirty = false;
    updateConfigIndicator();
    showToast(`Saved "${name}"`, "success");
  } catch (e) { showToast(`Save failed: ${e.message}`, "error"); }
}

async function promptSaveConfigAs() {
  const name = prompt("Configuration name:", App.state.configName || "");
  if (!name || !name.trim()) return;
  await saveConfig(name.trim());
}

async function showLoadConfigDialog() {
  if (App.state.configDirty) {
    const ans = confirm("You have unsaved changes. Discard and load a different configuration?");
    if (!ans) return;
  }
  try {
    const resp = await apiFetch("/api/configs");
    const configs = await resp.json();
    if (configs.length === 0) {
      showToast("No saved configurations", "error");
      return;
    }
    // Build modal
    const overlay = document.createElement("div");
    overlay.className = "config-modal-overlay";
    const modal = document.createElement("div");
    modal.className = "config-modal";
    modal.innerHTML = `<h3>Load Configuration</h3><div class="config-list"></div><button class="config-modal-close">Cancel</button>`;
    const list = modal.querySelector(".config-list");
    for (const cfg of configs) {
      const row = document.createElement("div");
      row.className = "config-item";
      const modified = new Date(cfg.modified).toLocaleString();
      const sizeKB = (cfg.size / 1024).toFixed(0);
      row.innerHTML = `<span class="config-item-name">${escapeHtml(cfg.name)}</span><span class="config-item-info">${modified} (${sizeKB} KB)</span><button class="config-item-delete" title="Delete">&#x2715;</button>`;
      row.querySelector(".config-item-name").addEventListener("click", async () => {
        overlay.remove();
        await loadConfig(cfg.name);
      });
      row.querySelector(".config-item-delete").addEventListener("click", async (e) => {
        e.stopPropagation();
        if (!confirm(`Delete configuration "${cfg.name}"?`)) return;
        try {
          await apiFetch(`/api/configs/${encodeURIComponent(cfg.name)}`, { method: "DELETE" });
          row.remove();
          if (!list.children.length) {
            overlay.remove();
            showToast("No saved configurations", "error");
          }
        } catch (err) { showToast(`Delete failed: ${err.message}`, "error"); }
      });
      list.appendChild(row);
    }
    modal.querySelector(".config-modal-close").addEventListener("click", () => overlay.remove());
    overlay.addEventListener("click", (e) => { if (e.target === overlay) overlay.remove(); });
    overlay.appendChild(modal);
    document.body.appendChild(overlay);
  } catch (e) { showToast(`Failed to list configs: ${e.message}`, "error"); }
}

async function loadConfig(name) {
  showToast(`Loading "${name}"...`);
  try {
    const resp = await apiFetch("/api/configs/load", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name }),
    });
    if (!resp.ok) {
      const err = await resp.json();
      showToast(`Load failed: ${err.error}`, "error");
      return;
    }
    App.state.configName = name;
    App.state.configDirty = false;
    updateConfigIndicator();
    // Full reload (same pattern as resetDatabase)
    loadConstants();
    await loadElements();
    loadPlumbingElements();
    loadGeometry();
    loadViews();
    loadOutlineTable();
    await loadVariants();
    ensureActiveVariantValid();
    showToast(`Loaded "${name}"`, "success");
  } catch (e) { showToast(`Load failed: ${e.message}`, "error"); }
}

function escapeHtml(str) {
  const d = document.createElement("div");
  d.textContent = str;
  return d.innerHTML;
}

async function resetDatabase() {
  if (!confirm("This will delete the current database and recreate it with default seed data. All changes will be lost.\n\nContinue?")) return;
  showToast("Resetting database\u2026");
  try {
    const resp = await apiFetch("/api/reset-database", { method: "POST" });
    const data = await resp.json();
    if (data.ok) {
      showToast("Database reset \u2014 reloading", "success");
      hideErrorBanner();
      loadConstants();
      await loadElements();
      loadGeometry();
      loadViews();
    } else {
      showToast("Reset failed: " + (data.issues || []).join(", "), "error");
    }
  } catch (e) {
    showToast("Reset failed: " + e.message, "error");
  }
}
