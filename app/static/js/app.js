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
    showGrid: false,
    showOpenings: true,
    showFurniture: true,
    showRooms: true,
    showDoors: true,
    showClearance: false,
    variant: "standard",
    measureStart: null,
    isDragging: false,
    dragStart: null,
    lastPan: null,
    constSortKey: "name",
    constSortAsc: true,
    svgView: { pan: { x: 0, y: 0 }, zoom: 1 },
    outlineChain: [],
    outlineSelectedSeq: null,
  },
  els: {},
  sse: null,
};


/* ========== INITIALISATION ========== */

document.addEventListener("DOMContentLoaded", () => {
  cacheElements();
  setupEventListeners();
  setupOutlineToolbar();
  connectSSE();
  loadViews();
  loadConstants();
  loadGeometry();
  loadElements();
  loadBuildLabel();
});

function cacheElements() {
  const ids = [
    "canvas", "canvas-transform", "viewport",
    "layer-outline", "layer-inner", "layer-walls",
    "layer-openings", "layer-doors", "layer-furniture", "layer-clearance",
    "layer-rooms", "layer-points",
    "layer-labels", "layer-dims", "layer-measure",
    "grid-rect", "svg-view-container",
    "coord-display", "connection-status", "zoom-level",
    "selection-info", "measure-info",
    "view-tabs", "const-category-filter", "const-search",
    "constants-table", "outline-table", "outline-closure-indicator",
    "outline-add-btn", "outline-remove-btn", "openings-table",
    "rough-openings-table", "interior-walls-table", "furniture-table",
    "props-empty", "props-detail", "props-title", "props-table",
    "show-points", "show-labels", "show-dims", "show-grid",
    "show-openings", "show-furniture", "show-rooms",
    "show-doors", "show-clearance",
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

  App.sse.addEventListener("constants_changed", () => {
    loadConstants();
  });

  App.sse.addEventListener("element_changed", () => {
    loadElements();
  });

  App.sse.addEventListener("geometry_changed", () => {
    loadGeometry();
  });

  App.sse.addEventListener("undo_status", (e) => {
    const data = JSON.parse(e.data);
    updateUndoButtons(data.can_undo, data.can_redo);
  });

  App.sse.addEventListener("outline_changed", () => {
    loadOutlineTable();
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
    if (!resp.ok) throw new Error(await resp.text());
    App.state.geometry = await resp.json();
    if (App.state.activeView === "interactive") renderCanvas();
    updateOpeningsTable();
    updateElementsTable();
  } catch (e) {
    console.error("Geometry load failed:", e);
    showToast("Geometry computation error", "error");
  }
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

  // Remaining generated SVG view tabs
  for (const v of App.state.views) {
    if (v.name === "floorplan") continue; // already added above
    if (v.svg_path.endsWith(".pdf")) continue; // skip PDFs for now
    const tab = document.createElement("button");
    tab.className = "view-tab";
    tab.textContent = v.label;
    tab.dataset.view = v.name;
    tab.onclick = () => switchView(v.name);
    container.appendChild(tab);
  }

  // Re-append variant selector
  container.appendChild(vs);
}

function switchView(viewName) {
  App.state.activeView = viewName;
  document.querySelectorAll(".view-tab").forEach(t => {
    t.classList.toggle("active", t.dataset.view === viewName);
  });

  // Show variant selector for interactive and floorplan views
  const showVariant = viewName === "interactive" || viewName === "floorplan";
  App.els["variant-selector"].style.display = showVariant ? "inline-block" : "none";

  if (viewName === "interactive") {
    App.els["canvas"].style.display = "block";
    App.els["svg-view-container"].style.display = "none";
    App.els["viewport"].style.cursor = App.state.activeTool === "pan" ? "grab" : "crosshair";
    renderCanvas();
  } else {
    App.els["canvas"].style.display = "none";
    App.els["svg-view-container"].style.display = "block";
    App.els["viewport"].style.cursor = "grab";
    loadSVGView(viewName);
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
  } catch (e) {
    const p = document.createElement("p");
    p.style.cssText = "padding:20px;color:#f88";
    p.textContent = "Error loading SVG: " + e.message;
    container.innerHTML = "";
    container.appendChild(p);
  }
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


/* ========== CANVAS RENDERING ========== */

function renderCanvas() {
  const g = App.state.geometry;
  if (!g) return;

  clearLayers();
  renderOutline(g);
  renderInnerWalls(g);
  renderInteriorWalls(g);
  renderOpenings(g);
  renderDoors(g);
  renderFurniture(g);
  renderClearanceZones(g);
  renderRoomLabels(g);
  renderPoints(g);
  renderDimensions(g);

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
    "layer-labels", "layer-dims", "layer-measure"];
  for (const id of layers) {
    App.els[id].innerHTML = "";
  }
}

function polyToStr(poly) {
  return poly.map(p => `${p[0]},${-p[1]}`).join(" ");
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

function renderDoors(g) {
  if (!App.state.showDoors) return;
  const layer = App.els["layer-doors"];
  for (const door of (g.door_arcs || [])) {
    for (const leaf of door.leaves) {
      // Door line (hinge → open tip)
      layer.appendChild(svgEl("line", {
        x1: leaf.hinge[0], y1: -leaf.hinge[1],
        x2: leaf.tip[0], y2: -leaf.tip[1],
        class: "door-line",
      }));
      // Swing arc polyline
      const pts = leaf.arc_pts.map(p => `${p[0]},${-p[1]}`).join(" ");
      layer.appendChild(svgEl("polyline", {
        points: pts,
        class: "door-arc",
      }));
      // Hinge point circle
      layer.appendChild(svgEl("circle", {
        cx: leaf.hinge[0], cy: -leaf.hinge[1], r: 0.04,
        class: "door-hinge",
      }));
    }
  }
}

function renderClearanceZones(g) {
  if (!App.state.showClearance) return;
  const layer = App.els["layer-clearance"];
  for (const cz of (g.clearance_zones || [])) {
    layer.appendChild(svgEl("polygon", {
      points: polyToStr(cz.poly),
      class: "clearance-zone",
    }));
  }
}

function renderFurniture(g) {
  if (!App.state.showFurniture) return;
  const layer = App.els["layer-furniture"];

  // Render variant items (comprehensive set) if available; empty dict = no items
  if (g.variant_items !== undefined) {
    // Build override lookup from DB elements
    const overrides = {};
    for (const e of (App.state.elements || [])) {
      if (e.type === "furniture" || e.type === "appliance" || e.type === "fixture") {
        let props = e.properties;
        if (typeof props === "string") props = JSON.parse(props);
        if (props && props.source === "override") {
          overrides[e.name] = props;
        }
      }
    }

    // Sort: non-stacked items first, stacked items on top (SVG paint order)
    const entries = Object.entries(g.variant_items);
    entries.sort((a, b) => (a[1].stacked ? 1 : 0) - (b[1].stacked ? 1 : 0));
    for (const [name, item] of entries) {
      const cssClass = `item-${item.type} selectable` + (item.stacked ? " item-stacked" : "");
      // Apply override offset if present
      const ov = overrides[name];
      const ox = ov ? (ov.offset_x || 0) : 0;
      const oy = ov ? (ov.offset_y || 0) : 0;

      if (item.shape === "circle") {
        const c = item.center;
        const el = svgEl("circle", {
          cx: c[0] + ox, cy: -(c[1] + oy), r: item.radius,
          class: cssClass,
          "data-type": item.type,
          "data-name": name,
        });
        el.addEventListener("click", (e) => selectElement(item.type, name, item, e));
        layer.appendChild(el);
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
        el.addEventListener("click", (e) => selectElement(item.type, name, item, e));
        layer.appendChild(el);
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

    // Group for name + area text (clickable when SF)
    const nameEl = svgEl("text", {
      x: e, y: -n + (hasArea ? -0.15 : 0),
      class: "room-label" + (hasArea ? " room-label-sf" : ""),
      "text-anchor": "middle",
      "dominant-baseline": "middle",
    });
    nameEl.textContent = lbl.name;
    if (lbl.poly) {
      nameEl.style.cursor = "pointer";
      nameEl.addEventListener("click", () => {
        const hl = layer.querySelector(`.room-highlight[data-room="${lbl.name}"]`);
        if (hl) hl.dispatchEvent(new Event("click"));
      });
    }
    layer.appendChild(nameEl);

    if (hasArea) {
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
  if (App.state.activeView !== "interactive") {
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


/* ========== SELECTION ========== */

function selectElement(type, name, data, event) {
  if (event) event.stopPropagation();

  // After a move-tool drag, suppress the click that fires on mouseup
  if (MoveTool._suppressClick) {
    MoveTool._suppressClick = false;
    return;
  }

  // Ctrl+Click: toggle in/out of multi-selection
  if (event && (event.ctrlKey || event.metaKey)) {
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
}

function clearSelection() {
  document.querySelectorAll(".selected-highlight").forEach(el => {
    el.classList.remove("selected-highlight");
  });
  document.querySelectorAll(".multi-selected").forEach(el => {
    el.classList.remove("multi-selected");
  });
  App.state.selection = null;
  App.state.selections = [];
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
    addPropRow(tbody, "Easting", fmtFtIn(data.pos[0]));
    addPropRow(tbody, "Northing", fmtFtIn(data.pos[1]));
  } else if (type === "wall") {
    const b = data.bbox;
    addPropRow(tbody, "Width", fmtFtIn(b.e - b.w));
    addPropRow(tbody, "Height", fmtFtIn(b.n - b.s));
    addPropRow(tbody, "West", fmtFtIn(b.w));
    addPropRow(tbody, "South", fmtFtIn(b.s));
    addPropRow(tbody, "East", fmtFtIn(b.e));
    addPropRow(tbody, "North", fmtFtIn(b.n));
    // Show related constants
    const related = findRelatedConstants(name);
    if (related.length > 0) {
      addPropRow(tbody, "—", "Related Constants");
      for (const c of related) {
        addPropRow(tbody, c.name, fmtConstProp(c), true, c.name);
      }
    }
  } else if (type === "opening" || type === "rough_opening") {
    addPropRow(tbody, "Name", data.name);
    if (data.seg_start) addPropRow(tbody, "Segment", `${data.seg_start}–${data.seg_end}`);
    if (data.wall_name) addPropRow(tbody, "Wall", data.wall_name);
    if (data.width) addPropRow(tbody, "Width", fmtFtIn(data.width));
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
        addPropRow(tbody, c.name, fmtConstProp(c), true, c.name);
      }
    }
  } else if (type === "appliance" || type === "furniture") {
    // SEL-8: Enhanced furniture/appliance properties
    const b = data.bbox;
    const w = b.e - b.w;
    const d = b.n - b.s;
    addPropRow(tbody, "Type", type);
    addPropRow(tbody, "Width", fmtFtIn(w));
    addPropRow(tbody, "Depth", fmtFtIn(d));
    if (data.center) {
      addPropRow(tbody, "Center E", fmtFtIn(data.center[0]));
      addPropRow(tbody, "Center N", fmtFtIn(data.center[1]));
    } else {
      addPropRow(tbody, "Center E", fmtFtIn((b.w + b.e) / 2));
      addPropRow(tbody, "Center N", fmtFtIn((b.s + b.n) / 2));
    }
    if (data.rotation !== undefined) addPropRow(tbody, "Rotation", data.rotation.toFixed(1) + "°");
    const related = findRelatedConstants(name);
    if (related.length > 0) {
      addPropRow(tbody, "—", "Related Constants");
      for (const c of related) {
        addPropRow(tbody, c.name, fmtConstProp(c), true, c.name);
      }
    }
  }
}

function fmtConstProp(c) {
  return formatConstValue(c);
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
  const resp = await apiFetch("/api/outline?_=" + Date.now());
  const chain = await resp.json();
  App.state.outlineChain = chain;
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

    tbody.appendChild(tr);
  }
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
    const data = await resp.json();
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
    const data = await resp.json();
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
    tr.innerHTML = `
      <td>${op.name}</td>
      <td>${op.seg_start}–${op.seg_end}</td>
      <td>${fmtFtIn(w)}</td>
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
      <td>${fmtFtIn(ro.width)}</td>
      <td>${ro.orientation}</td>
    `;
    tr.classList.add("selectable");
    tr.addEventListener("click", () => selectElement("rough_opening", ro.name, ro));
    tbody2.appendChild(tr);
  }
}


/* ========== ELEMENTS TABLE ========== */

async function loadElements() {
  try {
    const [elemResp, doorResp] = await Promise.all([
      fetch("/api/elements"),
      fetch("/api/doors"),
    ]);
    App.state.elements = await elemResp.json();
    App.state.doors = await doorResp.json();
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
        // Thickness is the smaller dimension, length is the larger
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
        if (iwData) selectElement("wall", wall.name, iwData, e);
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
  addDoorDropdownRow(tbody, "Hinge", door.hinge_side, roName, "hinge_side");
  addDoorDropdownRow(tbody, "Swing", door.swing_direction, roName, "swing_direction");
  addDoorDropdownRow(tbody, "Type", door.door_type, roName, "door_type",
    ["single", "double"]);

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

function addDoorDropdownRow(tbody, label, value, openingName, field, options) {
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
  App.els["show-rooms"].addEventListener("change", (e) => {
    App.state.showRooms = e.target.checked;
    renderCanvas();
  });
  App.els["show-doors"].addEventListener("change", (e) => {
    App.state.showDoors = e.target.checked;
    renderCanvas();
  });
  App.els["show-clearance"].addEventListener("change", (e) => {
    App.state.showClearance = e.target.checked;
    renderCanvas();
  });

  // Variant selector
  App.els["variant-select"].addEventListener("change", (e) => {
    App.state.variant = e.target.value;
    loadGeometry();
    if (App.state.activeView === "floorplan") loadSVGView("floorplan");
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
      if (tab.dataset.panel === "elements") updateElementsTable();
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
  // SVG view: any left-click starts pan drag
  if (App.state.activeView !== "interactive") {
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
  } else if (e.button === 0 && App.state.activeTool === "select") {
    // OE-1: outline F-point drag (only if select tool + points visible)
    if (typeof outlineEditorMouseDown === "function") {
      outlineEditorMouseDown(e);
    }
  } else if (e.button === 0 && App.state.activeTool === "move") {
    moveToolMouseDown(e);
  } else if (e.button === 0 && App.state.activeTool === "measure") {
    const rect = App.els["viewport"].getBoundingClientRect();
    const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
    App.state.measureStart = [wx, wy];
  }
}

function onMouseMove(e) {
  // SVG view pan drag
  if (App.state.activeView !== "interactive") {
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

  if (App.state.isDragging) {
    App.state.isDragging = false;
    if (App.state.activeView !== "interactive") {
      App.els["viewport"].style.cursor = "grab";
    } else {
      App.els["viewport"].style.cursor = App.state.activeTool === "pan" ? "grab" : "crosshair";
    }
  }

  if (App.state.activeView === "interactive" &&
      App.state.activeTool === "measure" && App.state.measureStart) {
    const rect = App.els["viewport"].getBoundingClientRect();
    const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
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

  if (App.state.activeView !== "interactive") {
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

function onKeyDown(e) {
  if (e.target.tagName === "INPUT" || e.target.tagName === "SELECT") return;

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
    case "f": case "F": fitToWindow(); break;
    case "Enter":
      if (App.state.activeTool === "move" && App.state.selection) {
        showOffsetDialog();
      }
      break;
    case "Escape":
      if (MoveTool.active) {
        // Cancel drag: remove ghosts, reset state
        for (const g of MoveTool.origTransforms) {
          if (g.ghost && g.ghost.parentNode) g.ghost.remove();
        }
        MoveTool.active = false;
        MoveTool.origTransforms = [];
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

function updateUndoButtons(canUndo, canRedo) {
  const undoBtn = document.querySelector('[data-action="undo"]');
  const redoBtn = document.querySelector('[data-action="redo"]');
  if (undoBtn) undoBtn.disabled = !canUndo;
  if (redoBtn) redoBtn.disabled = !canRedo;
}

async function doUndo() {
  const resp = await fetch("/api/undo", { method: "POST" });
  const data = await resp.json();
  if (resp.ok) {
    showToast("Undo: " + (data.description || data.action), "success");
    updateUndoButtons(data.can_undo, data.can_redo);
    loadConstants();
    loadGeometry();
    loadElements();
  } else {
    showToast(data.error || "Nothing to undo", "warning");
  }
}

async function doRedo() {
  const resp = await fetch("/api/redo", { method: "POST" });
  const data = await resp.json();
  if (resp.ok) {
    showToast("Redo: " + (data.description || data.action), "success");
    updateUndoButtons(data.can_undo, data.can_redo);
    loadConstants();
    loadGeometry();
    loadElements();
  } else {
    showToast(data.error || "Nothing to redo", "warning");
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


/* ========== DIMENSION LINES ========== */

function renderDimensions(g) {
  const layer = App.els["layer-dims"];
  if (!App.state.showDims || !g.dimensions) return;

  for (const [name, dim] of Object.entries(g.dimensions)) {
    const x1 = dim.A[0], y1 = -dim.A[1];
    const x2 = dim.B[0], y2 = -dim.B[1];

    // Dimension line
    layer.appendChild(svgEl("line", {
      x1, y1, x2, y2, class: "dim-line"
    }));

    // Tick marks (perpendicular to line)
    const dx = x2 - x1, dy = y2 - y1;
    const len = Math.sqrt(dx * dx + dy * dy);
    if (len < 0.01) continue;
    const px = -dy / len * 0.15, py = dx / len * 0.15;
    for (const [tx, ty] of [[x1, y1], [x2, y2]]) {
      layer.appendChild(svgEl("line", {
        x1: tx - px, y1: ty - py, x2: tx + px, y2: ty + py,
        class: "dim-line"
      }));
    }

    // Label at midpoint, rotated for readability
    const mx = (x1 + x2) / 2, my = (y1 + y2) / 2;
    let ang = Math.atan2(dy, dx) * 180 / Math.PI;
    if (ang >= 90) ang -= 180;
    else if (ang < -90) ang += 180;
    const angRad = ang * Math.PI / 180;
    const lx = mx + 0.15 * Math.sin(angRad);
    const ly = my - 0.15 * Math.cos(angRad);

    const label = svgEl("text", {
      x: lx, y: ly, class: "dim-label",
      transform: `rotate(${ang},${lx},${ly})`
    });
    label.textContent = fmtFtIn(dim.dist);
    layer.appendChild(label);
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
      if (App.state.activeView === "interactive") {
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
    case "export-all": {
      showToast("Exporting all products...");
      try {
        await apiFetch("/api/regenerate", { method: "POST" });
        showToast("Export complete", "success");
      } catch (e) { showToast(`Export failed: ${e.message}`, "error"); }
      break;
    }
    case "reset-constants": {
      if (!confirm("Reset all constants and outline to original values?")) return;
      try {
        await apiFetch("/api/constants/reset", { method: "POST" });
        await loadConstants();
        await loadOutlineTable();
        showToast("Constants and outline reset to defaults", "success");
      } catch (e) { showToast(`Reset failed: ${e.message}`, "error"); }
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
