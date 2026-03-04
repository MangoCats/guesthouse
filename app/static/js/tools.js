/* ========== Tool Utilities ========== */
"use strict";

/**
 * Client-side hosted-openings map (mirrors app/elements.py IW_HOSTED_OPENINGS).
 * Used for cascade warnings when deleting walls.
 */
const IW_HOSTED_OPENINGS = {
  IW1: ["RO1"], IW2: [], IW2O: ["RO4"], IW2S: [], IW3: [], IW4: [],
  IW5: [], IW6: ["RO5"], IW7: [], IW8: [], IW9: ["RO3", "RO7"],
  IW11: ["RO2", "RO6"], IW12: [],
};

/* ========== Move Tool ========== */

/**
 * Client-side IW move axis mapping (mirrors app/elements.py IW_MOVE_AXIS).
 * axis: "x" or "y" — the world-coordinate axis the constant affects.
 * sign: +1 or -1 — delta_constant = drag_delta_on_axis * sign.
 */
const IW_MOVE_AXIS = {
  IW1:  { axis: "y", sign: -1 },
  IW2:  { axis: "x", sign: +1 },
  IW2O: null,  // not movable
  IW2S: { axis: "x", sign: +1 },
  IW3:  { axis: "x", sign: +1 },
  IW4:  { axis: "x", sign: +1 },
  IW5:  { axis: "y", sign: -1 },
  IW6:  { axis: "y", sign: -1 },
  IW7:  { axis: "y", sign: +1 },
  IW8:  null,  // not movable
  IW9:  { axis: "x", sign: +1 },
  IW11: { axis: "x", sign: +1 },
  IW12: { axis: "y", sign: +1 },
};

/** Minimum screen-pixel drag distance before a drag actually starts. */
const DRAG_THRESHOLD = 4;

const MoveTool = {
  /** True once drag threshold is exceeded and ghosts are created. */
  active: false,
  /** True from mousedown until mouseup (pending drag). */
  pending: false,
  /** True while commitMove is in progress (prevents new drags during reload). */
  _committing: false,
  startScreen: null,   // {x, y} screen coords at mousedown
  startWorld: null,     // [wx, wy] at mousedown
  targets: [],          // [{type, name, elementId, svgEl}]
  origTransforms: [],   // original transform data for ghosts
  /** Suppress the next selectElement() call after a committed drag. */
  _suppressClick: false,
};


/**
 * Find the DB element record matching a rendered SVG element.
 * Returns {id, type, name, properties, variant} or null.
 */
function findElementRecord(type, name) {
  const elements = App.state.elements || [];
  if (type === "wall") {
    return elements.find(e => e.type === "wall" && e.name === name) || null;
  }
  // For furniture/appliance, check if there's an override element
  return elements.find(e => e.name === name) || null;
}


/**
 * Get the SVG element(s) matching a selection.
 */
function getSelectionSvgEl(sel) {
  if (!sel) return null;
  return document.querySelector(
    `[data-name="${sel.name}"][data-type="${sel.type}"]`
  );
}


/**
 * Create a ghost (semi-transparent clone) of an SVG element.
 */
function createGhost(svgEl) {
  if (!svgEl) return null;
  const clone = svgEl.cloneNode(true);
  clone.classList.add("move-ghost");
  clone.classList.remove("selectable", "selected-highlight", "multi-selected");
  clone.removeAttribute("data-name");
  clone.removeAttribute("data-type");
  // Insert ghost into measure layer (top)
  App.els["layer-measure"].appendChild(clone);
  return clone;
}


/**
 * Shift polygon points string by (dx, dy) in SVG coordinates.
 */
function shiftPolygonPoints(originalPoints, dxSvg, dySvg) {
  return originalPoints.split(/\s+/).map(pair => {
    const [x, y] = pair.split(",").map(Number);
    return `${x + dxSvg},${y + dySvg}`;
  }).join(" ");
}


/**
 * Apply axis constraint, shift-constrain, and grid snap to a delta.
 * Returns [dx, dy] in world coordinates.
 */
function applyMoveConstraints(dx, dy, shiftKey) {
  // Shift key: constrain to dominant axis (TL-7)
  if (shiftKey) {
    if (Math.abs(dx) > Math.abs(dy)) dy = 0;
    else dx = 0;
  }

  // Grid snap (TL-8): snap to 1 inch = 1/12 ft
  if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    dx = Math.round(dx / snap) * snap;
    dy = Math.round(dy / snap) * snap;
  }

  // For single IW walls, constrain to move axis
  if (MoveTool.targets.length === 1) {
    const t = MoveTool.targets[0];
    if (t.type === "wall" && IW_MOVE_AXIS[t.name]) {
      const axis = IW_MOVE_AXIS[t.name].axis;
      if (axis === "x") dy = 0;
      else dx = 0;
    }
  }

  return [dx, dy];
}


/** Element types that can be dragged with the move tool. */
const MOVABLE_TYPES = new Set(["wall", "furniture", "appliance", "fixture", "dimension", "label"]);


/**
 * Find the topmost movable SVG element at a screen coordinate.
 * Uses elementsFromPoint to look through overlapping layers
 * (e.g. point markers above walls) and find the first movable element.
 */
function findMovableAtPoint(clientX, clientY) {
  const hits = document.elementsFromPoint(clientX, clientY);
  for (const hit of hits) {
    const el = hit.closest("[data-name][data-type]");
    if (el && MOVABLE_TYPES.has(el.getAttribute("data-type"))) {
      return el;
    }
  }
  return null;
}


/* ========== Move Tool Mouse Handlers ========== */

function moveToolMouseDown(e) {
  if (e.button !== 0) return;
  if (App.state.activeView !== "interactive") return;

  // Don't start a new drag while a previous move is being committed/reloaded
  if (MoveTool._committing) return;

  // Reset stuck _suppressClick from a previous drag where no click fired
  // (happens when mousedown and mouseup land on different elements).
  MoveTool._suppressClick = false;

  // Find the topmost movable element at the click point, looking through
  // non-movable elements like point markers, labels, and openings.
  const targetEl = findMovableAtPoint(e.clientX, e.clientY);
  if (!targetEl) return; // clicked empty space or non-movable element

  const targetName = targetEl.getAttribute("data-name");
  const targetType = targetEl.getAttribute("data-type");

  // If this element is part of a multi-selection, move the whole group.
  // Otherwise, move just the clicked element.
  const inMultiSel = App.state.selections.length > 0 &&
    App.state.selections.some(s => s.name === targetName && s.type === targetType);

  let targets;
  if (inMultiSel) {
    targets = [];
    for (const s of App.state.selections) {
      const rec = findElementRecord(s.type, s.name);
      const el = getSelectionSvgEl(s);
      if (el) {
        targets.push({ type: s.type, name: s.name, elementId: rec ? rec.id : null, svgEl: el });
      }
    }
  } else {
    const rec = findElementRecord(targetType, targetName);
    targets = [{ type: targetType, name: targetName, elementId: rec ? rec.id : null, svgEl: targetEl }];
  }

  if (targets.length === 0) return;

  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  // Enter "pending" state — don't create ghosts yet.
  // Ghosts are created only after the drag threshold is exceeded (in mousemove).
  MoveTool.pending = true;
  MoveTool.active = false;
  MoveTool.startScreen = { x: e.clientX, y: e.clientY };
  MoveTool.startWorld = [wx, wy];
  MoveTool.targets = targets;
  MoveTool.origTransforms = [];

  e.preventDefault();
}

function moveToolMouseMove(e) {
  if (!MoveTool.pending && !MoveTool.active) return;

  // Check drag threshold before starting actual drag
  if (MoveTool.pending && !MoveTool.active) {
    const dx = e.clientX - MoveTool.startScreen.x;
    const dy = e.clientY - MoveTool.startScreen.y;
    if (Math.sqrt(dx * dx + dy * dy) < DRAG_THRESHOLD) return;

    // Threshold exceeded — promote to active drag, create ghosts
    MoveTool.active = true;
    MoveTool.pending = false;
    for (const t of MoveTool.targets) {
      const ghost = createGhost(t.svgEl);
      if (ghost) {
        if (ghost.tagName === "polygon") {
          MoveTool.origTransforms.push({
            ghost,
            origPoints: ghost.getAttribute("points"),
            isCircle: false,
          });
        } else if (ghost.tagName === "circle") {
          MoveTool.origTransforms.push({
            ghost,
            origCx: parseFloat(ghost.getAttribute("cx")),
            origCy: parseFloat(ghost.getAttribute("cy")),
            isCircle: true,
          });
        }
      }
    }
  }

  if (!MoveTool.active) return;

  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  let [dx, dy] = applyMoveConstraints(
    wx - MoveTool.startWorld[0],
    wy - MoveTool.startWorld[1],
    e.shiftKey,
  );

  // SVG coordinates: dx stays same, dy is negated
  const dxSvg = dx;
  const dySvg = -dy;

  // Update ghost positions
  for (const g of MoveTool.origTransforms) {
    if (g.isCircle) {
      g.ghost.setAttribute("cx", g.origCx + dxSvg);
      g.ghost.setAttribute("cy", g.origCy + dySvg);
    } else {
      g.ghost.setAttribute("points", shiftPolygonPoints(g.origPoints, dxSvg, dySvg));
    }
  }
}

function moveToolMouseUp(e) {
  if (!MoveTool.pending && !MoveTool.active) return;

  const wasActive = MoveTool.active;

  // Remove ghosts
  for (const g of MoveTool.origTransforms) {
    if (g.ghost && g.ghost.parentNode) g.ghost.remove();
  }

  // If the drag never started (below threshold), treat as a click-to-select.
  if (!wasActive) {
    MoveTool.pending = false;
    MoveTool.active = false;
    MoveTool.origTransforms = [];
    // Let the click event fire normally to select the element
    return;
  }

  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  let [dx, dy] = applyMoveConstraints(
    wx - MoveTool.startWorld[0],
    wy - MoveTool.startWorld[1],
    e.shiftKey,
  );

  MoveTool.active = false;
  MoveTool.pending = false;
  MoveTool.origTransforms = [];

  // Only commit if there was actual movement
  if (Math.abs(dx) < 1e-6 && Math.abs(dy) < 1e-6) return;

  // Suppress the click event that will fire after this mouseup,
  // so it doesn't re-select/change selection after the move.
  MoveTool._suppressClick = true;

  commitMove(MoveTool.targets, dx, dy);
}


/**
 * Commit a move to the server for one or more targets.
 */
async function commitMove(targets, dx, dy) {
  MoveTool._committing = true;
  try {
    for (const t of targets) {
      let elementId = t.elementId;

      // For furniture/appliance with no DB record, create override first
      // Dimension move: update perpendicular offset (TL-12)
      if (t.type === "dimension" && elementId) {
        const elemRec = (App.state.elements || []).find(e => e.id === elementId);
        if (elemRec) {
          const props = typeof elemRec.properties === "string"
            ? JSON.parse(elemRec.properties) : elemRec.properties;
          if (props.start && props.end) {
            // Compute perpendicular component of drag
            const sdx = props.end[0] - props.start[0];
            const sdy = props.end[1] - props.start[1];
            const slen = Math.sqrt(sdx * sdx + sdy * sdy);
            if (slen > 1e-9) {
              const perpX = -sdy / slen, perpY = sdx / slen;
              const perpDelta = dx * perpX + dy * perpY;
              props.offset = (props.offset || 0) + perpDelta;
              await fetch(`/api/elements/${elementId}`, {
                method: "PUT",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({ properties: props }),
              });
            }
          }
        }
        continue;
      }

      // Label move: update offset or position (LABEL-2)
      if (t.type === "label" && elementId) {
        const elemRec = (App.state.elements || []).find(e => e.id === elementId);
        if (elemRec) {
          const props = typeof elemRec.properties === "string"
            ? JSON.parse(elemRec.properties) : elemRec.properties;
          if (props.source === "room") {
            props.offset_e = (props.offset_e || 0) + dx;
            props.offset_n = (props.offset_n || 0) + dy;
          } else if (props.position) {
            props.position = [props.position[0] + dx, props.position[1] + dy];
          }
          await fetch(`/api/elements/${elementId}`, {
            method: "PUT",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ properties: props }),
          });
        }
        continue;
      }

      if (!elementId && (t.type === "furniture" || t.type === "appliance" || t.type === "fixture")) {
        try {
          const resp = await fetch("/api/elements", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
              type: t.type,
              name: t.name,
              variant: App.state.variant,
              properties: { offset_x: 0, offset_y: 0, source: "override" },
            }),
          });
          if (resp.ok) {
            const created = await resp.json();
            elementId = created.id;
          } else {
            showToast(`Failed to create override for ${t.name}`, "error");
            continue;
          }
        } catch (err) {
          showToast(`Error creating override: ${err.message}`, "error");
          continue;
        }
      }

      if (!elementId) {
        showToast(`Cannot move ${t.name}: no element record`, "error");
        continue;
      }

      try {
        const resp = await fetch(`/api/elements/${elementId}/move`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ dx, dy }),
        });
        const data = await resp.json();
        if (resp.ok) {
          updateUndoButtons(data.can_undo, data.can_redo);
        } else {
          showToast(data.error || `Move failed for ${t.name}`, "error");
        }
      } catch (err) {
        showToast(`Move error: ${err.message}`, "error");
      }
    }

    // Reload state: elements/constants FIRST (updates App.state data),
    // then geometry LAST (calls renderCanvas which reads that data).
    await loadElements();
    await loadConstants();
    await loadGeometry();
  } finally {
    MoveTool._committing = false;
  }
}


/* ========== Draw Wall Tool (TL-15, TL-16, TL-17) ========== */

const DrawWallTool = {
  /** First click position [wx, wy], or null if not started. */
  start: null,
  /** Preview SVG elements. */
  previewLine: null,
  /** Default wall thickness in feet (4 inches). */
  defaultThickness: 4.0 / 12.0,
};

/** Auto-generate the next custom wall name (CW1, CW2, ...). */
function nextCustomWallName() {
  const elements = App.state.elements || [];
  let max = 0;
  for (const e of elements) {
    if (e.type === "wall" && e.name && e.name.startsWith("CW")) {
      const n = parseInt(e.name.slice(2), 10);
      if (n > max) max = n;
    }
  }
  return `CW${max + 1}`;
}

/** Compute a rectangle polygon from start/end/thickness. */
function wallPoly(start, end, thickness) {
  const dx = end[0] - start[0];
  const dy = end[1] - start[1];
  const len = Math.sqrt(dx * dx + dy * dy);
  if (len < 1e-9) return null;
  // Perpendicular unit vector
  const px = -dy / len * (thickness / 2);
  const py = dx / len * (thickness / 2);
  return [
    [start[0] + px, start[1] + py],
    [end[0] + px, end[1] + py],
    [end[0] - px, end[1] - py],
    [start[0] - px, start[1] - py],
  ];
}

function drawWallMouseDown(e) {
  if (e.button !== 0) return;
  const rect = App.els["viewport"].getBoundingClientRect();
  let [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  // Grid snap
  if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    wx = Math.round(wx / snap) * snap;
    wy = Math.round(wy / snap) * snap;
  }

  if (!DrawWallTool.start) {
    // First click — set start point
    DrawWallTool.start = [wx, wy];
  } else {
    // Second click — create the wall
    createDrawnWall(DrawWallTool.start, [wx, wy]);
    // Clear preview
    if (DrawWallTool.previewLine && DrawWallTool.previewLine.parentNode) {
      DrawWallTool.previewLine.remove();
    }
    DrawWallTool.previewLine = null;
    DrawWallTool.start = null;
  }
}

function drawWallMouseMove(e) {
  if (!DrawWallTool.start) return;
  const rect = App.els["viewport"].getBoundingClientRect();
  let [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    wx = Math.round(wx / snap) * snap;
    wy = Math.round(wy / snap) * snap;
  }

  // Update preview line from start to cursor
  const layer = document.getElementById("layer-measure");
  if (!DrawWallTool.previewLine) {
    DrawWallTool.previewLine = svgEl("line", { class: "draw-preview" });
    layer.appendChild(DrawWallTool.previewLine);
  }
  const s = DrawWallTool.start;
  DrawWallTool.previewLine.setAttribute("x1", s[0]);
  DrawWallTool.previewLine.setAttribute("y1", -s[1]);
  DrawWallTool.previewLine.setAttribute("x2", wx);
  DrawWallTool.previewLine.setAttribute("y2", -wy);
}

function cancelDrawWall() {
  DrawWallTool.start = null;
  if (DrawWallTool.previewLine && DrawWallTool.previewLine.parentNode) {
    DrawWallTool.previewLine.remove();
  }
  DrawWallTool.previewLine = null;
}

async function createDrawnWall(start, end) {
  const thickness = DrawWallTool.defaultThickness;
  const poly = wallPoly(start, end, thickness);
  if (!poly) return;

  const name = nextCustomWallName();
  const resp = await fetch("/api/elements", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      type: "wall",
      name: name,
      properties: {
        source: "drawn",
        start: start,
        end: end,
        thickness: thickness,
        poly: poly,
      },
    }),
  });
  if (resp.ok) {
    showToast(`Created wall ${name}`, "success");
    await loadElements();
    await loadGeometry();
  } else {
    const data = await resp.json();
    showToast(data.error || "Failed to create wall", "error");
  }
}


/**
 * Show offset dialog for precise move entry (TL-6).
 */
function showOffsetDialog() {
  const sel = App.state.selection;
  if (!sel) {
    showToast("Select an element first", "warning");
    return;
  }

  Dialog.show({
    title: `Move ${sel.name}`,
    fields: [
      {
        label: 'Offset (e.g. "6in east", "2ft north", or "1.5 west")',
        name: "offset",
        placeholder: "6in east",
      },
    ],
    onSubmit(values) {
      const parsed = parseOffsetString(values.offset);
      if (!parsed) {
        showToast("Could not parse offset. Use format: '6in east'", "error");
        return;
      }

      const selections = App.state.selections && App.state.selections.length > 0
        ? App.state.selections
        : [sel];

      const targets = [];
      for (const s of selections) {
        const rec = findElementRecord(s.type, s.name);
        targets.push({
          type: s.type,
          name: s.name,
          data: s.data,
          elementId: rec ? rec.id : null,
        });
      }

      commitMove(targets, parsed.dx, parsed.dy);
    },
  });
}


/* ========== Dimension Tool (TL-11) ========== */

const DimTool = {
  start: null,
  startAnchor: null,
  previewLine: null,
  snapIndicator: null,
  snapTargets: null,
};

/**
 * Build flat list of snap targets from current geometry.
 * Each entry: { type, target, face?, pos: [E, N] }
 */
function buildSnapTargets(g) {
  const targets = [];
  if (!g) return targets;

  // Points (F/W/C-series)
  if (g.points) {
    for (const [name, pt] of Object.entries(g.points)) {
      targets.push({ type: "point", target: name, pos: pt });
    }
  }

  // Wall faces
  if (g.interior_walls) {
    const faces = [
      ["south", 0, 1], ["east", 1, 2], ["north", 2, 3], ["west", 3, 0],
    ];
    for (const [wname, wall] of Object.entries(g.interior_walls)) {
      const poly = wall.poly;
      if (!poly || poly.length < 4) continue;
      for (const [face, i, j] of faces) {
        const mid = [(poly[i][0] + poly[j][0]) / 2, (poly[i][1] + poly[j][1]) / 2];
        targets.push({ type: "wall_face", target: wname, face, pos: mid });
      }
    }
  }

  // Outer opening faces
  if (g.outer_openings) {
    const faces = [
      ["south", 0, 1], ["east", 1, 2], ["north", 2, 3], ["west", 3, 0],
    ];
    for (const op of g.outer_openings) {
      if (!op.poly || op.poly.length < 4) continue;
      for (const [face, i, j] of faces) {
        const mid = [(op.poly[i][0] + op.poly[j][0]) / 2, (op.poly[i][1] + op.poly[j][1]) / 2];
        targets.push({ type: "opening_face", target: op.name, face, pos: mid });
      }
    }
  }

  // Rough opening faces
  if (g.rough_openings) {
    const faces = [
      ["south", 0, 1], ["east", 1, 2], ["north", 2, 3], ["west", 3, 0],
    ];
    for (const ro of g.rough_openings) {
      if (!ro.poly || ro.poly.length < 4) continue;
      for (const [face, i, j] of faces) {
        const mid = [(ro.poly[i][0] + ro.poly[j][0]) / 2, (ro.poly[i][1] + ro.poly[j][1]) / 2];
        targets.push({ type: "opening_face", target: ro.name, face, pos: mid });
      }
    }
  }

  return targets;
}

/**
 * Find nearest snap target within pixel threshold.
 * Returns { anchor: {type, target, face?}, pos: [E, N] } or null.
 */
function findNearestSnap(wx, wy, snapTargets, thresholdPx) {
  if (!snapTargets || snapTargets.length === 0) return null;
  const worldThreshold = (thresholdPx || 12) / (App.state.zoom || 1);
  let best = null;
  let bestDist = worldThreshold;
  for (const t of snapTargets) {
    const dx = t.pos[0] - wx;
    const dy = t.pos[1] - wy;
    const d = Math.sqrt(dx * dx + dy * dy);
    if (d < bestDist) {
      bestDist = d;
      const anchor = { type: t.type, target: t.target };
      if (t.face) anchor.face = t.face;
      best = { anchor, pos: t.pos };
    }
  }
  return best;
}

/**
 * Show or update the snap indicator circle.
 */
function updateSnapIndicator(snapResult) {
  const layer = document.getElementById("layer-measure");
  if (!snapResult) {
    if (DimTool.snapIndicator && DimTool.snapIndicator.parentNode) {
      DimTool.snapIndicator.remove();
    }
    DimTool.snapIndicator = null;
    return;
  }
  if (!DimTool.snapIndicator) {
    DimTool.snapIndicator = svgEl("circle", { class: "snap-indicator", r: 0.1 });
    layer.appendChild(DimTool.snapIndicator);
  }
  DimTool.snapIndicator.setAttribute("cx", snapResult.pos[0]);
  DimTool.snapIndicator.setAttribute("cy", -snapResult.pos[1]);
}

function nextDimensionName() {
  const elements = App.state.elements || [];
  let max = 0;
  for (const e of elements) {
    if (e.type === "dimension" && e.name && e.name.startsWith("UD")) {
      const n = parseInt(e.name.slice(2), 10);
      if (n > max) max = n;
    }
  }
  return `UD${max + 1}`;
}

function dimToolMouseDown(e) {
  if (e.button !== 0) return;
  const rect = App.els["viewport"].getBoundingClientRect();
  let [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  // Build snap targets on first use
  if (!DimTool.snapTargets) {
    DimTool.snapTargets = buildSnapTargets(App.state.geometry);
  }

  // Try geometry snap first, then grid snap
  let anchor = null;
  const snapResult = findNearestSnap(wx, wy, DimTool.snapTargets, 12);
  if (snapResult) {
    wx = snapResult.pos[0];
    wy = snapResult.pos[1];
    anchor = snapResult.anchor;
  } else if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    wx = Math.round(wx / snap) * snap;
    wy = Math.round(wy / snap) * snap;
  }

  if (!DimTool.start) {
    DimTool.start = [wx, wy];
    DimTool.startAnchor = anchor;
  } else {
    createDimension(DimTool.start, [wx, wy], DimTool.startAnchor, anchor);
    if (DimTool.previewLine && DimTool.previewLine.parentNode) {
      DimTool.previewLine.remove();
    }
    DimTool.previewLine = null;
    DimTool.start = null;
    DimTool.startAnchor = null;
    updateSnapIndicator(null);
  }
}

function dimToolMouseMove(e) {
  const rect = App.els["viewport"].getBoundingClientRect();
  let [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  // Build snap targets on first use
  if (!DimTool.snapTargets) {
    DimTool.snapTargets = buildSnapTargets(App.state.geometry);
  }

  // Show snap indicator even before first click
  const snapResult = findNearestSnap(wx, wy, DimTool.snapTargets, 12);
  updateSnapIndicator(snapResult);
  if (snapResult) {
    wx = snapResult.pos[0];
    wy = snapResult.pos[1];
  } else if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    wx = Math.round(wx / snap) * snap;
    wy = Math.round(wy / snap) * snap;
  }

  if (!DimTool.start) return;

  const layer = document.getElementById("layer-measure");
  if (!DimTool.previewLine) {
    DimTool.previewLine = svgEl("line", { class: "draw-preview" });
    layer.appendChild(DimTool.previewLine);
  }
  const s = DimTool.start;
  DimTool.previewLine.setAttribute("x1", s[0]);
  DimTool.previewLine.setAttribute("y1", -s[1]);
  DimTool.previewLine.setAttribute("x2", wx);
  DimTool.previewLine.setAttribute("y2", -wy);
}

function cancelDimTool() {
  DimTool.start = null;
  DimTool.startAnchor = null;
  DimTool.snapTargets = null;
  if (DimTool.previewLine && DimTool.previewLine.parentNode) {
    DimTool.previewLine.remove();
  }
  DimTool.previewLine = null;
  updateSnapIndicator(null);
}

async function createDimension(start, end, startAnchor, endAnchor) {
  const name = nextDimensionName();
  const props = {
    source: "user",
    start: start,
    end: end,
    offset: 0.0,
    label_rotation: "parallel",
  };
  if (startAnchor) props.start_anchor = startAnchor;
  if (endAnchor) props.end_anchor = endAnchor;
  const resp = await fetch("/api/elements", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      type: "dimension",
      name: name,
      properties: props,
    }),
  });
  if (resp.ok) {
    showToast(`Created dimension ${name}`, "success");
    // Auto-enable User Dims if not already
    if (!App.state.showUserDims) {
      App.state.showUserDims = true;
      App.els["show-user-dims"].checked = true;
    }
    await loadElements();
    await loadGeometry();
  } else {
    const data = await resp.json();
    showToast(data.error || "Failed to create dimension", "error");
  }
}


/* ========== Label Tool (LABEL-1) ========== */

function labelToolMouseDown(e) {
  if (e.button !== 0) return;
  const rect = App.els["viewport"].getBoundingClientRect();
  let [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    wx = Math.round(wx / snap) * snap;
    wy = Math.round(wy / snap) * snap;
  }

  Dialog.show({
    title: "New Label",
    fields: [{ label: "Text", name: "text", placeholder: "Label text" }],
    async onSubmit(vals) {
      if (!vals.text || !vals.text.trim()) {
        showToast("Label text required", "error");
        return;
      }
      const name = nextLabelName();
      const resp = await fetch("/api/elements", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          type: "label",
          name: name,
          properties: {
            source: "user",
            text: vals.text.trim(),
            position: [wx, wy],
            rotation: 0,
            font_size: 0.25,
          },
        }),
      });
      Dialog.close();
      if (resp.ok) {
        showToast(`Created label ${name}`, "success");
        await loadElements();
        await loadGeometry();
      } else {
        const data = await resp.json();
        showToast(data.error || "Failed to create label", "error");
      }
    },
  });
}

function nextLabelName() {
  const elements = App.state.elements || [];
  let max = 0;
  for (const e of elements) {
    if (e.type === "label" && e.name && e.name.startsWith("UL")) {
      const n = parseInt(e.name.slice(2), 10);
      if (n > max) max = n;
    }
  }
  return `UL${max + 1}`;
}
