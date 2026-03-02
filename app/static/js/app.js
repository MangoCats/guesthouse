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
    pan: { x: 0, y: 0 },
    zoom: 1,
    showPoints: true,
    showLabels: true,
    showDims: false,
    showGrid: false,
    showOpenings: true,
    showFurniture: true,
    variant: "standard",
    measureStart: null,
    isDragging: false,
    dragStart: null,
    lastPan: null,
    constSortKey: "name",
    constSortAsc: true,
  },
  els: {},
  sse: null,
};


/* ========== INITIALISATION ========== */

document.addEventListener("DOMContentLoaded", () => {
  cacheElements();
  setupEventListeners();
  connectSSE();
  loadViews();
  loadConstants();
  loadGeometry();
});

function cacheElements() {
  const ids = [
    "canvas", "canvas-transform", "viewport",
    "layer-outline", "layer-inner", "layer-walls",
    "layer-openings", "layer-furniture", "layer-points",
    "layer-labels", "layer-dims", "layer-measure",
    "grid-rect", "svg-view-container",
    "coord-display", "connection-status", "zoom-level",
    "selection-info", "measure-info",
    "view-tabs", "const-category-filter", "const-search",
    "constants-table", "outline-table", "openings-table",
    "rough-openings-table",
    "props-empty", "props-detail", "props-title", "props-table",
    "show-points", "show-labels", "show-dims", "show-grid",
    "show-openings", "show-furniture",
    "variant-select", "variant-selector",
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

  App.sse.addEventListener("geometry_changed", () => {
    loadGeometry();
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


/* ========== DATA LOADING ========== */

async function loadConstants() {
  const resp = await fetch("/api/constants");
  App.state.constants = await resp.json();
  renderConstantsTable();
  loadCategories();
}

async function loadCategories() {
  const resp = await fetch("/api/constants/categories");
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
    if (!resp.ok) throw new Error(await resp.text());
    App.state.geometry = await resp.json();
    if (App.state.activeView === "interactive") renderCanvas();
    updateOpeningsTable();
  } catch (e) {
    console.error("Geometry load failed:", e);
    showToast("Geometry computation error", "error");
  }
}

async function loadViews() {
  const resp = await fetch("/api/views");
  App.state.views = await resp.json();
  renderViewTabs();
}


/* ========== VIEW TABS ========== */

function renderViewTabs() {
  const container = App.els["view-tabs"];
  container.innerHTML = "";

  // Interactive view tab (always first)
  const interTab = document.createElement("button");
  interTab.className = "view-tab active";
  interTab.textContent = "Interactive";
  interTab.dataset.view = "interactive";
  interTab.onclick = () => switchView("interactive");
  container.appendChild(interTab);

  // Generated SVG view tabs
  for (const v of App.state.views) {
    if (v.svg_path.endsWith(".pdf")) continue; // skip PDFs for now
    const tab = document.createElement("button");
    tab.className = "view-tab";
    tab.textContent = v.label;
    tab.dataset.view = v.name;
    tab.onclick = () => switchView(v.name);
    container.appendChild(tab);
  }
}

function switchView(viewName) {
  App.state.activeView = viewName;
  document.querySelectorAll(".view-tab").forEach(t => {
    t.classList.toggle("active", t.dataset.view === viewName);
  });

  // Show variant selector only for interactive view
  App.els["variant-selector"].style.display = viewName === "interactive" ? "inline-block" : "none";

  if (viewName === "interactive") {
    App.els["canvas"].style.display = "block";
    App.els["svg-view-container"].style.display = "none";
    renderCanvas();
  } else {
    App.els["canvas"].style.display = "none";
    App.els["svg-view-container"].style.display = "block";
    loadSVGView(viewName);
  }
}

async function loadSVGView(viewName) {
  const container = App.els["svg-view-container"];
  container.innerHTML = "<p style='padding:20px;color:#888'>Loading...</p>";
  try {
    const resp = await fetch(`/api/svg/${viewName}`);
    if (!resp.ok) {
      container.innerHTML = `<p style='padding:20px;color:#f88'>SVG not available. Click File > Regenerate to generate it.</p>`;
      return;
    }
    const svg = await resp.text();
    container.innerHTML = svg;
  } catch (e) {
    container.innerHTML = `<p style='padding:20px;color:#f88'>Error loading SVG: ${e.message}</p>`;
  }
}


/* ========== CANVAS RENDERING ========== */

function renderCanvas() {
  const g = App.state.geometry;
  if (!g) return;

  clearLayers();
  renderOutline(g);
  renderInnerWalls(g);
  renderInteriorWalls(g);
  renderOpenings(g);
  renderFurniture(g);
  renderPoints(g);

  if (App.state.zoom === 1 && App.state.pan.x === 0 && App.state.pan.y === 0) {
    fitToWindow();
  } else {
    applyTransform();
  }
}

function clearLayers() {
  const layers = ["layer-outline", "layer-inner", "layer-walls",
    "layer-openings", "layer-furniture", "layer-points",
    "layer-labels", "layer-dims", "layer-measure"];
  for (const id of layers) {
    App.els[id].innerHTML = "";
  }
}

function polyToStr(poly) {
  return poly.map(p => `${p[0]},${-p[1]}`).join(" ");
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

function renderFurniture(g) {
  if (!App.state.showFurniture) return;
  const layer = App.els["layer-furniture"];

  // Render variant items (comprehensive set) if available
  if (g.variant_items && Object.keys(g.variant_items).length > 0) {
    for (const [name, item] of Object.entries(g.variant_items)) {
      const cssClass = `item-${item.type} selectable`;
      if (item.shape === "circle") {
        const c = item.center;
        const el = svgEl("circle", {
          cx: c[0], cy: -c[1], r: item.radius,
          class: cssClass,
          "data-type": item.type,
          "data-name": name,
        });
        el.addEventListener("click", (e) => selectElement(item.type, name, item, e));
        layer.appendChild(el);
      } else {
        const el = svgEl("polygon", {
          points: polyToStr(item.poly),
          class: cssClass,
          "data-type": item.type,
          "data-name": name,
        });
        el.addEventListener("click", (e) => selectElement(item.type, name, item, e));
        layer.appendChild(el);
      }
    }
    return;
  }

  // Fallback: legacy appliances/furniture (if variant_items not present)
  for (const [name, item] of Object.entries(g.appliances || {})) {
    const poly = item.clip || item.poly;
    const el = svgEl("polygon", {
      points: polyToStr(poly),
      class: "appliance-fill selectable",
      "data-type": "appliance",
      "data-name": name,
    });
    el.addEventListener("click", (e) => selectElement("appliance", name, item, e));
    layer.appendChild(el);
  }

  for (const [name, item] of Object.entries(g.furniture || {})) {
    const el = svgEl("polygon", {
      points: polyToStr(item.poly),
      class: "furniture-fill selectable",
      "data-type": "furniture",
      "data-name": name,
    });
    el.addEventListener("click", (e) => selectElement("furniture", name, item, e));
    layer.appendChild(el);
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


/* ========== SELECTION ========== */

function selectElement(type, name, data, event) {
  if (event) event.stopPropagation();

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
}

function clearSelection() {
  document.querySelectorAll(".selected-highlight").forEach(el => {
    el.classList.remove("selected-highlight");
  });
  App.state.selection = null;
  App.els["selection-info"].textContent = "No selection";
  App.els["props-empty"].style.display = "block";
  App.els["props-detail"].style.display = "none";
}

function showProperties(type, name, data) {
  App.els["props-empty"].style.display = "none";
  App.els["props-detail"].style.display = "block";
  App.els["props-title"].textContent = `${type.toUpperCase()}: ${name}`;

  const tbody = App.els["props-table"].querySelector("tbody");
  tbody.innerHTML = "";

  if (type === "point") {
    addPropRow(tbody, "Easting", `${data.pos[0].toFixed(4)} ft`);
    addPropRow(tbody, "Northing", `${data.pos[1].toFixed(4)} ft`);
    addPropRow(tbody, "E (inches)", `${(data.pos[0] * 12).toFixed(2)}"`);
    addPropRow(tbody, "N (inches)", `${(data.pos[1] * 12).toFixed(2)}"`);
  } else if (type === "wall") {
    const b = data.bbox;
    addPropRow(tbody, "Width", `${((b.e - b.w) * 12).toFixed(1)}"`);
    addPropRow(tbody, "Height", `${((b.n - b.s) * 12).toFixed(1)}"`);
    addPropRow(tbody, "West", `${b.w.toFixed(4)} ft`);
    addPropRow(tbody, "South", `${b.s.toFixed(4)} ft`);
    addPropRow(tbody, "East", `${b.e.toFixed(4)} ft`);
    addPropRow(tbody, "North", `${b.n.toFixed(4)} ft`);
    // Show related constants
    const related = findRelatedConstants(name);
    if (related.length > 0) {
      addPropRow(tbody, "—", "Related Constants");
      for (const c of related) {
        addPropRow(tbody, c.name, `${c.value.toFixed(6)} ${c.unit}`, true, c.name);
      }
    }
  } else if (type === "opening" || type === "rough_opening") {
    addPropRow(tbody, "Name", data.name);
    if (data.seg_start) addPropRow(tbody, "Segment", `${data.seg_start}–${data.seg_end}`);
    if (data.wall_name) addPropRow(tbody, "Wall", data.wall_name);
    if (data.width) addPropRow(tbody, "Width", `${(data.width * 12).toFixed(1)}"`);
    if (data.orientation) addPropRow(tbody, "Orient", data.orientation);
    if (data.poly) {
      const dx = data.poly[1][0] - data.poly[0][0];
      const dy = data.poly[1][1] - data.poly[0][1];
      const w = Math.sqrt(dx * dx + dy * dy);
      addPropRow(tbody, "Actual width", `${(w * 12).toFixed(1)}"`);
    }
    const related = findRelatedConstants(data.name);
    if (related.length > 0) {
      addPropRow(tbody, "—", "Related Constants");
      for (const c of related) {
        addPropRow(tbody, c.name, `${c.value.toFixed(6)} ${c.unit}`, true, c.name);
      }
    }
  } else if (type === "appliance" || type === "furniture") {
    const b = data.bbox;
    addPropRow(tbody, "Width", `${((b.e - b.w) * 12).toFixed(1)}"`);
    addPropRow(tbody, "Depth", `${((b.n - b.s) * 12).toFixed(1)}"`);
    addPropRow(tbody, "Center E", `${((b.w + b.e) / 2).toFixed(4)} ft`);
    addPropRow(tbody, "Center N", `${((b.s + b.n) / 2).toFixed(4)} ft`);
    const related = findRelatedConstants(name);
    if (related.length > 0) {
      addPropRow(tbody, "—", "Related Constants");
      for (const c of related) {
        addPropRow(tbody, c.name, `${c.value.toFixed(6)} ${c.unit}`, true, c.name);
      }
    }
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
    tdVal.appendChild(inp);
  } else {
    tdVal.textContent = value;
  }
  tr.appendChild(tdVal);
  tbody.appendChild(tr);
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

    const tdUnit = document.createElement("td");
    tdUnit.textContent = c.unit;
    tr.appendChild(tdUnit);

    const tdCat = document.createElement("td");
    tdCat.textContent = c.category;
    tdCat.style.color = categoryColor(c.category);
    tr.appendChild(tdCat);

    tbody.appendChild(tr);
  }
}

function formatConstValue(c) {
  if (c.unit === "ft") {
    const inches = c.value * 12;
    if (Math.abs(inches - Math.round(inches)) < 0.001) {
      return `${Math.round(inches)}"`;
    }
    return `${inches.toFixed(3)}"`;
  }
  return c.value.toFixed(6);
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

async function handleConstantEdit(name, rawValue) {
  // Parse value: handle inches (e.g. "35"") and fractions
  let value;
  const trimmed = rawValue.trim();
  if (trimmed.endsWith('"')) {
    value = parseFloat(trimmed.slice(0, -1)) / 12.0;
  } else if (trimmed.includes("/")) {
    const parts = trimmed.split("/");
    value = parseFloat(parts[0]) / parseFloat(parts[1]);
  } else {
    value = parseFloat(trimmed);
  }

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
    showToast(`${name} = ${value.toFixed(6)}`, "success");

    // Update local state
    const c = App.state.constants.find(c => c.name === name);
    if (c) c.value = value;
  } catch (e) {
    showToast(`Error: ${e.message}`, "error");
  }
}


/* ========== OUTLINE TABLE ========== */

async function loadOutlineTable() {
  const resp = await fetch("/api/outline");
  const chain = await resp.json();
  const tbody = App.els["outline-table"].querySelector("tbody");
  tbody.innerHTML = "";
  for (const seg of chain) {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${seg.seq}</td>
      <td>${seg.seg_type}</td>
      <td>${seg.seg_type === "L" ? (seg.distance || 0).toFixed(4) : (seg.radius || 0).toFixed(6)}</td>
      <td>${seg.seg_type === "L" ? "—" : (seg.sweep || 0).toFixed(6)}</td>
      <td>${seg.end_name}</td>
    `;
    tbody.appendChild(tr);
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
    const w = Math.sqrt(dx * dx + dy * dy) * 12;
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${op.name}</td>
      <td>${op.seg_start}–${op.seg_end}</td>
      <td>${w.toFixed(1)}"</td>
    `;
    tr.classList.add("selectable");
    tr.addEventListener("click", () => selectElement("opening", op.name, op));
    tbody1.appendChild(tr);
  }

  // Rough openings
  const tbody2 = App.els["rough-openings-table"].querySelector("tbody");
  tbody2.innerHTML = "";
  for (const ro of (g.rough_openings || [])) {
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${ro.name}</td>
      <td>${ro.wall_name}</td>
      <td>${(ro.width * 12).toFixed(1)}"</td>
      <td>${ro.orientation}</td>
    `;
    tr.classList.add("selectable");
    tr.addEventListener("click", () => selectElement("rough_opening", ro.name, ro));
    tbody2.appendChild(tr);
  }
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

  // Keyboard shortcuts
  document.addEventListener("keydown", onKeyDown);

  // Tool buttons
  document.querySelectorAll(".tool-btn").forEach(btn => {
    btn.addEventListener("click", () => setTool(btn.dataset.tool));
  });

  // Display toggles
  App.els["show-points"].addEventListener("change", (e) => {
    App.state.showPoints = e.target.checked;
    renderCanvas();
  });
  App.els["show-labels"].addEventListener("change", (e) => {
    App.state.showLabels = e.target.checked;
    renderCanvas();
  });
  App.els["show-dims"].addEventListener("change", (e) => {
    App.state.showDims = e.target.checked;
    renderCanvas();
  });
  App.els["show-grid"].addEventListener("change", (e) => {
    App.state.showGrid = e.target.checked;
    applyTransform();
  });
  App.els["show-openings"].addEventListener("change", (e) => {
    App.state.showOpenings = e.target.checked;
    renderCanvas();
  });
  App.els["show-furniture"].addEventListener("change", (e) => {
    App.state.showFurniture = e.target.checked;
    renderCanvas();
  });

  // Variant selector
  App.els["variant-select"].addEventListener("change", (e) => {
    App.state.variant = e.target.value;
    loadGeometry();
  });

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
    });
  });

  // Menu actions
  document.querySelectorAll("[data-action]").forEach(btn => {
    btn.addEventListener("click", () => handleMenuAction(btn.dataset.action));
  });

  // Window resize
  window.addEventListener("resize", () => {
    if (App.state.activeView === "interactive" && App.state.geometry) {
      fitToWindow();
    }
  });
}


/* ========== MOUSE HANDLERS ========== */

function onMouseDown(e) {
  if (e.button === 1 || (e.button === 0 && App.state.activeTool === "pan")) {
    // Pan mode
    App.state.isDragging = true;
    App.state.dragStart = { x: e.clientX, y: e.clientY };
    App.state.lastPan = { ...App.state.pan };
    App.els["viewport"].style.cursor = "grabbing";
    e.preventDefault();
  } else if (e.button === 0 && App.state.activeTool === "measure") {
    const rect = App.els["viewport"].getBoundingClientRect();
    const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
    App.state.measureStart = [wx, wy];
  }
}

function onMouseMove(e) {
  const rect = App.els["viewport"].getBoundingClientRect();
  const sx = e.clientX - rect.left;
  const sy = e.clientY - rect.top;
  const [wx, wy] = screenToWorld(sx, sy);

  // Update coordinate display
  App.els["coord-display"].textContent =
    `E: ${wx.toFixed(3)}  N: ${wy.toFixed(3)}  (${(wx*12).toFixed(1)}" / ${(wy*12).toFixed(1)}")`;

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
  if (App.state.isDragging) {
    App.state.isDragging = false;
    App.els["viewport"].style.cursor = App.state.activeTool === "pan" ? "grab" : "crosshair";
  }

  if (App.state.activeTool === "measure" && App.state.measureStart) {
    const rect = App.els["viewport"].getBoundingClientRect();
    const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
    const [sx, sy] = App.state.measureStart;
    const dx = wx - sx;
    const dy = wy - sy;
    const dist = Math.sqrt(dx * dx + dy * dy);
    App.els["measure-info"].textContent =
      `Dist: ${dist.toFixed(4)} ft (${(dist * 12).toFixed(2)}")`;
    App.state.measureStart = null;
  }
}

function onWheel(e) {
  e.preventDefault();
  const rect = App.els["viewport"].getBoundingClientRect();
  const mx = e.clientX - rect.left;
  const my = e.clientY - rect.top;

  const factor = e.deltaY > 0 ? 0.9 : 1.1;
  const newZoom = App.state.zoom * factor;

  // Zoom toward cursor
  App.state.pan.x = mx - (mx - App.state.pan.x) * factor;
  App.state.pan.y = my - (my - App.state.pan.y) * factor;
  App.state.zoom = newZoom;

  applyTransform();
}

function onViewportClick(e) {
  // If clicking on the viewport background (not on an element), clear selection
  if (e.target === App.els["canvas"] || e.target === App.els["viewport"]) {
    clearSelection();
  }
}

function onKeyDown(e) {
  if (e.target.tagName === "INPUT" || e.target.tagName === "SELECT") return;

  switch (e.key) {
    case "v": case "V": setTool("select"); break;
    case "h": case "H": setTool("pan"); break;
    case "m": case "M": setTool("measure"); break;
    case "g": case "G": setTool("move"); break;
    case "f": case "F": fitToWindow(); break;
    case "Escape":
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


/* ========== TOOLS ========== */

function setTool(tool) {
  App.state.activeTool = tool;
  document.querySelectorAll(".tool-btn").forEach(btn => {
    btn.classList.toggle("active", btn.dataset.tool === tool);
  });
  const cursors = { select: "crosshair", pan: "grab", measure: "crosshair", move: "move" };
  App.els["viewport"].style.cursor = cursors[tool] || "crosshair";

  if (tool !== "measure") clearMeasure();
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
  label.textContent = `${(dist * 12).toFixed(1)}" (${dist.toFixed(3)}')`;
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
      const resp = await fetch("/api/regenerate", { method: "POST" });
      const data = await resp.json();
      showToast("All views regenerated", "success");
      break;
    }
    case "regenerate-view": {
      if (App.state.activeView === "interactive") {
        showToast("Cannot regenerate interactive view — edit constants instead");
        return;
      }
      showToast(`Regenerating ${App.state.activeView}...`);
      const resp = await fetch("/api/regenerate", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ view: App.state.activeView }),
      });
      break;
    }
    case "export-all": {
      showToast("Exporting all products...");
      const resp = await fetch("/api/regenerate", { method: "POST" });
      showToast("Export complete", "success");
      break;
    }
    case "reset-constants": {
      if (!confirm("Reset all constants to original values?")) return;
      await fetch("/api/constants/reset", { method: "POST" });
      await loadConstants();
      showToast("Constants reset to defaults", "success");
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
