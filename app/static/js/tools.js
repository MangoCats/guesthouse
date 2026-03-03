/* ========== Move Tool ========== */
"use strict";

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
const MOVABLE_TYPES = new Set(["wall", "furniture", "appliance", "fixture"]);


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
  for (const t of targets) {
    let elementId = t.elementId;

    // For furniture/appliance with no DB record, create override first
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

  // Reload state — await so canvas isn't rebuilt mid-interaction
  await loadGeometry();
  await loadElements();
  await loadConstants();
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
