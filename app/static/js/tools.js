/* ========== Tool Utilities ========== */
"use strict";

/* ========== Move Tool ==========
 * IW_MOVE_AXIS and IW_HOSTED_OPENINGS are no longer hardcoded here.
 * They are fetched from /api/iw-config at startup and stored in
 * App.state.iwConfig.  Use App.state.iwConfig.iw_move_axis[name] and
 * App.state.iwConfig.iw_hosted_openings[name] everywhere below.
 */

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
    const axisInfo = App.state.iwConfig.iw_move_axis[t.name];
    if (t.type === "wall" && axisInfo) {
      if (axisInfo.axis === "x") dy = 0;
      else dx = 0;
    }
  }

  return [dx, dy];
}


/** Element types that can be dragged with the move tool. */
const MOVABLE_TYPES = new Set(["wall", "furniture", "appliance", "fixture", "dimension", "label", "plumbing"]);


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
  if (!isCanvasView()) return;

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
        const tag = ghost.tagName;
        if (tag === "polygon" || tag === "polyline") {
          MoveTool.origTransforms.push({
            ghost,
            origPoints: ghost.getAttribute("points"),
            kind: "points",
          });
        } else if (tag === "circle") {
          MoveTool.origTransforms.push({
            ghost,
            origCx: parseFloat(ghost.getAttribute("cx")),
            origCy: parseFloat(ghost.getAttribute("cy")),
            kind: "circle",
          });
        } else if (tag === "g") {
          // Groups (plumbing fixtures/fittings): extract translate
          const tf = ghost.getAttribute("transform") || "";
          const m = tf.match(/translate\(\s*([\d.e+-]+)\s*,\s*([\d.e+-]+)\s*\)/);
          MoveTool.origTransforms.push({
            ghost,
            origTx: m ? parseFloat(m[1]) : 0,
            origTy: m ? parseFloat(m[2]) : 0,
            origTransform: tf,
            kind: "group",
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
    if (g.kind === "circle") {
      g.ghost.setAttribute("cx", g.origCx + dxSvg);
      g.ghost.setAttribute("cy", g.origCy + dySvg);
    } else if (g.kind === "group") {
      const newTf = g.origTransform.replace(
        /translate\([^)]*\)/,
        `translate(${g.origTx + dxSvg},${g.origTy + dySvg})`
      );
      g.ghost.setAttribute("transform", newTf);
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

      // Plumbing element move: shift all path waypoints
      if (t.type === "plumbing") {
        const pe = (App.state.plumbingElements || []).find(e => e.name === t.name);
        if (pe && pe.path) {
          const newPath = pe.path.map(pt => [pt[0] + dx, pt[1] + dy]);
          try {
            await fetch(`/api/plumbing/${pe.id}`, {
              method: "PUT",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ path: newPath }),
            });
          } catch (err) {
            showToast(`Plumbing move error: ${err.message}`, "error");
          }
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
    if (targets.some(t => t.type === "plumbing")) {
      await loadPlumbingElements();
      renderPlumbingCanvas();
    } else {
      await loadElements();
      await loadConstants();
      await loadGeometry();
    }
  } finally {
    MoveTool._committing = false;
  }
}


/* ========== (DrawWallTool removed — use Add > Wall wizard instead) ========== */



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

  // Wall faces: 2" interval snap points along each face
  if (g.interior_walls) {
    const FACE_IDX = [
      ["south", 0, 1], ["east", 1, 2], ["north", 2, 3], ["west", 3, 0],
    ];
    const SNAP_IN = 2;
    for (const [wname, wall] of Object.entries(g.interior_walls)) {
      const poly = wall.poly;
      if (!poly || poly.length < 4) continue;
      for (const [face, i, j] of FACE_IDX) {
        const p1 = poly[i], p2 = poly[j];
        const dx = p2[0] - p1[0], dy = p2[1] - p1[1];
        const lenFt = Math.hypot(dx, dy);
        if (lenFt < 1e-9) continue;
        const ux = dx / lenFt, uy = dy / lenFt;
        for (let distIn = 0; distIn <= lenFt * 12 + 1e-6; distIn += SNAP_IN) {
          const clamped = Math.min(distIn, Math.floor(lenFt * 12 / SNAP_IN) * SNAP_IN);
          const d = clamped / 12;
          targets.push({ type: "wall_face", target: wname, face, distIn: clamped,
                         pos: [p1[0] + d * ux, p1[1] + d * uy] });
          if (clamped >= lenFt * 12 - 1e-6) break;
        }
      }
    }
  }

  // Inner wall segments (W-series straight sections): 2" interval snap points,
  // excluding endpoints which are already covered as named points above.
  if (g.inner_segments) {
    const pts = g.points || {};
    const SNAP_IN = 2;
    for (const seg of g.inner_segments) {
      if (seg.type !== "line") continue;
      const p1 = pts[seg.start], p2 = pts[seg.end];
      if (!p1 || !p2) continue;
      const dx = p2[0] - p1[0], dy = p2[1] - p1[1];
      const lenFt = Math.hypot(dx, dy);
      const lenIn = lenFt * 12;
      if (lenIn < SNAP_IN * 2 - 1e-6) continue;
      const ux = dx / lenFt, uy = dy / lenFt;
      for (let distIn = SNAP_IN; distIn <= lenIn - SNAP_IN + 1e-6; distIn += SNAP_IN) {
        const d = distIn / 12;
        targets.push({ type: "inner_seg", target: seg.start,
                       pos: [p1[0] + d * ux, p1[1] + d * uy] });
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
      if (t.face !== undefined) anchor.face = t.face;
      if (t.distIn !== undefined) anchor.distIn = t.distIn;
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
    const anc = snapResult.anchor;
    if (anc && (anc.type === "wall_face" || anc.type === "inner_seg")) {
      const label = anc.type === "wall_face"
        ? `${anc.target} ${anc.face}  (${fmtFtIn(wx)}, ${fmtFtIn(wy)})`
        : `${anc.target}  (${fmtFtIn(wx)}, ${fmtFtIn(wy)})`;
      App.els["coord-display"].textContent = label;
    }
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
    // Auto-enable Dims if not already
    if (!App.state.showDims) {
      App.state.showDims = true;
      App.els["show-dims"].checked = true;
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

/* ========== Area Edit Tool ========== */
/**
 * Vertex drag editor for area_poly elements.
 * Activated from the area properties panel "Edit Vertices" button.
 */
const AreaEditTool = {
  active: false,
  elementId: null,
  name: null,
  vertices: [],        // formula spec list (may be strings or dicts)
  resolvedPoly: [],    // [[E, N], ...] from g.room_labels
  arcAdjustments: [],
  subtractElements: [],
  handles: [],         // SVG circle elements
  ghostLine: null,     // SVG polygon outline ghost
  dragIdx: null,       // index of handle being dragged
  dragHandle: null,    // SVG element being dragged
  startScreen: null,
  startWorld: null,
  _layer: null,
};

const AREA_SNAP_R = 0.25; // feet, snap threshold

/** Return candidate snap targets: [{spec, pt: [E,N]}] */
function _areaSnapTargets() {
  const g = App.state.geometry;
  if (!g) return [];
  const targets = [];
  // W-series and other named points
  for (const [pname, pt] of Object.entries(g.points || {})) {
    targets.push({ spec: pname, pt: [pt[0], pt[1]] });
  }
  // Interior wall corners
  for (const [wname, wdata] of Object.entries(g.interior_walls || {})) {
    const poly = wdata.poly || [];
    for (let i = 0; i < poly.length; i++) {
      targets.push({ spec: { element: wname, corner: i }, pt: [poly[i][0], poly[i][1]] });
    }
  }
  // Opening corners
  for (const op of [...(g.outer_openings || []), ...(g.rough_openings || [])]) {
    if (!op.poly) continue;
    for (let i = 0; i < op.poly.length; i++) {
      targets.push({ spec: { element: op.name, corner: i }, pt: [op.poly[i][0], op.poly[i][1]] });
    }
  }
  return targets;
}

function _areaFindSnap(wx, wy) {
  const targets = _areaSnapTargets();
  let best = null, bestDist = AREA_SNAP_R;
  for (const t of targets) {
    const d = Math.hypot(wx - t.pt[0], wy - t.pt[1]);
    if (d < bestDist) { bestDist = d; best = t; }
  }
  return best;
}

function _areaRenderHandles() {
  const tool = AreaEditTool;
  if (!tool.active || !tool._layer) return;
  // Remove old handles
  for (const h of tool.handles) h.remove();
  tool.handles = [];
  if (tool.ghostLine) { tool.ghostLine.remove(); tool.ghostLine = null; }

  // Ghost polygon outline
  const ptStr = tool.resolvedPoly.map(p => `${p[0]},${-p[1]}`).join(" ");
  tool.ghostLine = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
  tool.ghostLine.setAttribute("points", ptStr);
  tool.ghostLine.setAttribute("fill", "none");
  tool.ghostLine.setAttribute("stroke", "#4af");
  tool.ghostLine.setAttribute("stroke-width", "0.04");
  tool.ghostLine.setAttribute("stroke-dasharray", "0.1 0.05");
  tool._layer.appendChild(tool.ghostLine);

  // Vertex handles
  for (let i = 0; i < tool.resolvedPoly.length; i++) {
    const [e, n] = tool.resolvedPoly[i];
    const h = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    h.setAttribute("cx", e);
    h.setAttribute("cy", -n);
    h.setAttribute("r", "0.12");
    h.setAttribute("fill", "#4af");
    h.setAttribute("stroke", "#fff");
    h.setAttribute("stroke-width", "0.025");
    h.setAttribute("cursor", "grab");
    h.dataset.idx = i;
    h.addEventListener("mousedown", _areaHandleMouseDown);
    tool._layer.appendChild(h);
    tool.handles.push(h);
  }
}

function _areaHandleMouseDown(e) {
  const tool = AreaEditTool;
  if (!tool.active) return;
  e.preventDefault();
  e.stopPropagation();
  const idx = parseInt(e.currentTarget.dataset.idx, 10);
  tool.dragIdx = idx;
  tool.dragHandle = e.currentTarget;
  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
  tool.startScreen = { x: e.clientX, y: e.clientY };
  tool.startWorld = [wx, wy];
  tool.dragHandle.setAttribute("fill", "#fa0");
  tool.dragHandle.setAttribute("cursor", "grabbing");
}

async function _areaHandleMouseUp(e) {
  const tool = AreaEditTool;
  if (!tool.active || tool.dragIdx === null) return;
  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
  const snap = _areaFindSnap(wx, wy);
  let newSpec;
  if (snap) {
    newSpec = snap.spec;
    tool.resolvedPoly[tool.dragIdx] = [...snap.pt];
  } else {
    newSpec = [wx, wy];
    tool.resolvedPoly[tool.dragIdx] = [wx, wy];
  }
  tool.vertices[tool.dragIdx] = newSpec;
  tool.dragIdx = null;
  tool.dragHandle = null;
  _areaRenderHandles();
  // Commit to server
  await _areaCommitVertices();
}

function _areaHandleMouseMove(e) {
  const tool = AreaEditTool;
  if (!tool.active || tool.dragIdx === null || !tool.dragHandle) return;
  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
  // Move handle and ghost
  tool.dragHandle.setAttribute("cx", wx);
  tool.dragHandle.setAttribute("cy", -wy);
  // Update ghost polygon
  const pts = tool.resolvedPoly.map((p, i) =>
    i === tool.dragIdx ? `${wx},${-wy}` : `${p[0]},${-p[1]}`
  ).join(" ");
  if (tool.ghostLine) tool.ghostLine.setAttribute("points", pts);
  // Snap highlight
  const snap = _areaFindSnap(wx, wy);
  tool.dragHandle.setAttribute("fill", snap ? "#0f0" : "#fa0");
}

async function _areaCommitVertices() {
  const tool = AreaEditTool;
  try {
    await fetch(`/api/elements/${tool.elementId}/vertices`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        vertices: tool.vertices,
        arc_adjustments: tool.arcAdjustments,
        subtract_elements: tool.subtractElements,
      }),
    });
    await reloadAfterChange();
    // Update resolved poly from fresh geometry
    const g = App.state.geometry;
    if (g && g.room_labels) {
      const rl = g.room_labels.find(r => r.name === tool.name);
      if (rl && rl.poly) {
        tool.resolvedPoly = rl.poly.map(p => [...p]);
      }
    }
    _areaRenderHandles();
  } catch (err) {
    showToast("Failed to save vertices", "error");
  }
}

AreaEditTool.activate = async function(elementId, name, data) {
  const tool = AreaEditTool;
  tool.elementId = elementId;
  tool.name = name;
  tool.resolvedPoly = (data.poly || []).map(p => [...p]);
  tool.arcAdjustments = data.arc_adjustments || [];
  tool.subtractElements = [];
  tool._layer = App.els["layer-measure"];

  // Fetch the formula to get the spec vertices
  try {
    const resp = await fetch(`/api/formulas/${encodeURIComponent(name)}`);
    if (resp.ok) {
      const formulas = await resp.json();
      const polyRec = formulas.find(f => f.param_name === "poly");
      if (polyRec) {
        const fj = typeof polyRec.formula_json === "string"
          ? JSON.parse(polyRec.formula_json) : polyRec.formula_json;
        tool.vertices = fj.vertices || tool.resolvedPoly.map(p => [...p]);
        tool.arcAdjustments = fj.arc_adjustments || [];
        tool.subtractElements = fj.subtract_elements || [];
      } else {
        tool.vertices = tool.resolvedPoly.map(p => [...p]);
      }
    }
  } catch (_) {
    tool.vertices = tool.resolvedPoly.map(p => [...p]);
  }

  tool.active = true;
  _areaRenderHandles();

  // Register mouse event listeners on viewport
  const vp = App.els["viewport"];
  vp.addEventListener("mousemove", _areaHandleMouseMove);
  vp.addEventListener("mouseup", _areaHandleMouseUp);

  showToast(`Editing ${name} vertices — drag handles, click Done when finished`, "info");
};

AreaEditTool.deactivate = function() {
  const tool = AreaEditTool;
  tool.active = false;
  for (const h of tool.handles) h.remove();
  tool.handles = [];
  if (tool.ghostLine) { tool.ghostLine.remove(); tool.ghostLine = null; }
  const vp = App.els["viewport"];
  vp.removeEventListener("mousemove", _areaHandleMouseMove);
  vp.removeEventListener("mouseup", _areaHandleMouseUp);
  tool.elementId = null;
  tool.dragIdx = null;
  tool.dragHandle = null;
};


/* ========== Area Draw Tool ========== */
/**
 * Draw a new area polygon by clicking canvas vertices.
 * Finish: click near first vertex (≤ AREA_SNAP_R) or press Enter.
 */
const AreaDrawTool = {
  active: false,
  vertices: [],         // [[E, N], ...]
  _layer: null,
  _polyEl: null,        // live polygon preview
  _rubberLine: null,    // line from last vertex to cursor
  _handles: [],         // vertex circles
};

function _areaDrawRender() {
  const tool = AreaDrawTool;
  const layer = tool._layer;
  if (!layer) return;
  if (tool._polyEl) tool._polyEl.remove();
  if (tool._rubberLine) tool._rubberLine.remove();
  for (const h of tool._handles) h.remove();
  tool._handles = [];
  tool._polyEl = null;
  tool._rubberLine = null;

  if (tool.vertices.length === 0) return;

  // Polygon preview
  tool._polyEl = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
  tool._polyEl.setAttribute("points",
    tool.vertices.map(p => `${p[0]},${-p[1]}`).join(" "));
  tool._polyEl.setAttribute("fill", "#88f");
  tool._polyEl.setAttribute("fill-opacity", "0.15");
  tool._polyEl.setAttribute("stroke", "#88f");
  tool._polyEl.setAttribute("stroke-width", "0.04");
  layer.appendChild(tool._polyEl);

  // Vertex handles
  for (const [e, n] of tool.vertices) {
    const h = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    h.setAttribute("cx", e); h.setAttribute("cy", -n); h.setAttribute("r", "0.1");
    h.setAttribute("fill", "#88f"); h.setAttribute("stroke", "#fff");
    h.setAttribute("stroke-width", "0.02");
    layer.appendChild(h);
    tool._handles.push(h);
  }
}

function _areaDrawMouseMove(e) {
  const tool = AreaDrawTool;
  if (!tool.active || tool.vertices.length === 0) return;
  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
  if (tool._rubberLine) tool._rubberLine.remove();
  const last = tool.vertices[tool.vertices.length - 1];
  tool._rubberLine = document.createElementNS("http://www.w3.org/2000/svg", "line");
  tool._rubberLine.setAttribute("x1", last[0]); tool._rubberLine.setAttribute("y1", -last[1]);
  tool._rubberLine.setAttribute("x2", wx); tool._rubberLine.setAttribute("y2", -wy);
  tool._rubberLine.setAttribute("stroke", "#88f"); tool._rubberLine.setAttribute("stroke-width", "0.03");
  tool._rubberLine.setAttribute("stroke-dasharray", "0.1 0.05");
  tool._layer.appendChild(tool._rubberLine);
}

function _areaDrawClick(e) {
  const tool = AreaDrawTool;
  if (!tool.active) return;
  e.preventDefault(); e.stopPropagation();
  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  // Close if near first vertex
  if (tool.vertices.length >= 3) {
    const first = tool.vertices[0];
    if (Math.hypot(wx - first[0], wy - first[1]) <= AREA_SNAP_R) {
      _areaDrawFinish();
      return;
    }
  }

  tool.vertices.push([wx, wy]);
  _areaDrawRender();
}

async function _areaDrawFinish() {
  const tool = AreaDrawTool;
  if (tool.vertices.length < 3) {
    showToast("Need at least 3 vertices", "error");
    return;
  }
  // Prompt for name/label
  Dialog.open({
    title: "New Area",
    fields: [
      { id: "label", label: "Label", type: "text", value: "New Area" },
    ],
    onConfirm: async (vals) => {
      const label = vals.label.trim() || "New Area";
      const autoName = "AREA" + Date.now().toString().slice(-5);
      // Create element
      const resp1 = await fetch("/api/elements", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          type: "area",
          name: autoName,
          properties: { label, color: "#dddddd" },
        }),
      });
      if (!resp1.ok) { showToast("Failed to create area", "error"); return; }
      const created = await resp1.json();
      // Create formula
      const resp2 = await fetch(`/api/elements/${created.id}/vertices`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          vertices: tool.vertices.map(p => [...p]),
          arc_adjustments: [],
        }),
      });
      if (resp2.ok) {
        showToast(`Created area ${label}`, "success");
        AreaDrawTool.deactivate();
        await reloadAfterChange();
      } else {
        showToast("Failed to save area vertices", "error");
      }
    },
    onCancel: () => {},
  });
}

function _areaDrawKeyDown(e) {
  if (!AreaDrawTool.active) return;
  if (e.key === "Enter") { e.preventDefault(); _areaDrawFinish(); }
  if (e.key === "Escape") { e.preventDefault(); AreaDrawTool.deactivate(); }
}

AreaDrawTool.activate = function() {
  const tool = AreaDrawTool;
  tool.active = true;
  tool.vertices = [];
  tool._layer = App.els["layer-measure"];
  const vp = App.els["viewport"];
  vp.addEventListener("click", _areaDrawClick);
  vp.addEventListener("mousemove", _areaDrawMouseMove);
  document.addEventListener("keydown", _areaDrawKeyDown);
  showToast("Click to place vertices. Click near first vertex or press Enter to close.", "info");
};

AreaDrawTool.deactivate = function() {
  const tool = AreaDrawTool;
  tool.active = false;
  if (tool._polyEl) { tool._polyEl.remove(); tool._polyEl = null; }
  if (tool._rubberLine) { tool._rubberLine.remove(); tool._rubberLine = null; }
  for (const h of tool._handles) h.remove();
  tool._handles = [];
  tool.vertices = [];
  const vp = App.els["viewport"];
  vp.removeEventListener("click", _areaDrawClick);
  vp.removeEventListener("mousemove", _areaDrawMouseMove);
  document.removeEventListener("keydown", _areaDrawKeyDown);
};


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
