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

const MoveTool = {
  active: false,
  startScreen: null,   // {x, y} screen coords at mousedown
  startWorld: null,     // [wx, wy] at mousedown
  ghost: null,          // cloned SVG group
  targets: [],          // [{type, name, data, elementId, svgEl}]
  origTransforms: [],   // original transform data for ghosts
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
 * SVG Y is negated from world Y.
 */
function shiftPolygonPoints(originalPoints, dxSvg, dySvg) {
  return originalPoints.split(/\s+/).map(pair => {
    const [x, y] = pair.split(",").map(Number);
    return `${x + dxSvg},${y + dySvg}`;
  }).join(" ");
}


/* ========== Move Tool Mouse Handlers ========== */

function moveToolMouseDown(e) {
  if (e.button !== 0) return;
  if (App.state.activeView !== "interactive") return;

  // Need a selection to move
  const sel = App.state.selection;
  if (!sel) return;

  const svgEl = getSelectionSvgEl(sel);
  if (!svgEl) return;

  // Build targets from selections (multi-select) or single selection
  const selections = App.state.selections && App.state.selections.length > 0
    ? App.state.selections
    : [sel];

  MoveTool.targets = [];
  for (const s of selections) {
    const rec = findElementRecord(s.type, s.name);
    const el = getSelectionSvgEl(s);
    if (el) {
      MoveTool.targets.push({
        type: s.type,
        name: s.name,
        data: s.data,
        elementId: rec ? rec.id : null,
        svgEl: el,
      });
    }
  }

  if (MoveTool.targets.length === 0) return;

  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  MoveTool.active = true;
  MoveTool.startScreen = { x: e.clientX, y: e.clientY };
  MoveTool.startWorld = [wx, wy];

  // Create ghosts
  MoveTool.origTransforms = [];
  for (const t of MoveTool.targets) {
    const ghost = createGhost(t.svgEl);
    if (ghost) {
      // Store original points/position for shifting
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

  e.preventDefault();
  e.stopPropagation();
}

function moveToolMouseMove(e) {
  if (!MoveTool.active) return;

  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  let dx = wx - MoveTool.startWorld[0];
  let dy = wy - MoveTool.startWorld[1];

  // Shift key: constrain to dominant axis (TL-7)
  if (e.shiftKey) {
    if (Math.abs(dx) > Math.abs(dy)) {
      dy = 0;
    } else {
      dx = 0;
    }
  }

  // Grid snap (TL-8): snap to 1 inch = 1/12 ft
  if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    dx = Math.round(dx / snap) * snap;
    dy = Math.round(dy / snap) * snap;
  }

  // For IW walls, constrain to move axis
  if (MoveTool.targets.length === 1) {
    const t = MoveTool.targets[0];
    if (t.type === "wall" && IW_MOVE_AXIS[t.name]) {
      const axis = IW_MOVE_AXIS[t.name].axis;
      if (axis === "x") dy = 0;
      else dx = 0;
    }
  }

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
  if (!MoveTool.active) return;

  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  let dx = wx - MoveTool.startWorld[0];
  let dy = wy - MoveTool.startWorld[1];

  // Apply same constraints as mousemove
  if (e.shiftKey) {
    if (Math.abs(dx) > Math.abs(dy)) dy = 0;
    else dx = 0;
  }
  if (App.state.showGrid) {
    const snap = 1.0 / 12.0;
    dx = Math.round(dx / snap) * snap;
    dy = Math.round(dy / snap) * snap;
  }
  if (MoveTool.targets.length === 1) {
    const t = MoveTool.targets[0];
    if (t.type === "wall" && IW_MOVE_AXIS[t.name]) {
      const axis = IW_MOVE_AXIS[t.name].axis;
      if (axis === "x") dy = 0;
      else dx = 0;
    }
  }

  // Remove ghosts
  for (const g of MoveTool.origTransforms) {
    if (g.ghost && g.ghost.parentNode) {
      g.ghost.remove();
    }
  }

  MoveTool.active = false;
  MoveTool.origTransforms = [];

  // Only commit if there was actual movement
  if (Math.abs(dx) < 1e-6 && Math.abs(dy) < 1e-6) return;

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

  // Reload state
  loadGeometry();
  loadElements();
  loadConstants();
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
