/* ========== Outline Editor — Canvas F-point drag & arc handle ========== */
"use strict";

/**
 * OutlineEditor handles interactive canvas editing of the building outline:
 * - OE-1: Drag F-points to change segment distances/radii
 * - OE-2: Arc radius handle for adjusting arc radii
 */
const OutlineEditor = {
  active: false,
  pending: false,
  dragPoint: null,     // {name, seq, segType, origWorld}
  startScreen: null,   // {x, y}
  startWorld: null,    // [wx, wy]
  ghost: null,         // SVG circle element
};


/**
 * Check if a point name is an F-series outline point (not FC).
 */
function isFPoint(name) {
  return name && name.startsWith("F") && !name.startsWith("FC");
}


/**
 * Find the chain segment whose end_name matches a given F-point name.
 * Returns {seq, seg} or null.
 */
function findSegByEndName(name) {
  const chain = App.state.outlineChain;
  if (!chain || !chain.length) return null;
  for (const seg of chain) {
    if (seg.end_name === name) return { seq: seg.seq, seg };
  }
  return null;
}


/**
 * Called from onMouseDown when select tool is active and an F-point is clicked.
 * Returns true if outline drag started, false otherwise.
 */
function outlineEditorMouseDown(e) {
  if (App.state.activeTool !== "select" || !App.state.showPoints) return false;
  if (OutlineEditor._committing) return false;

  // Check if we clicked on an F-series point marker
  const target = e.target.closest("[data-type='point']");
  if (!target) return false;
  const name = target.dataset.name;
  if (!isFPoint(name)) return false;

  // Find which chain segment this point belongs to
  const match = findSegByEndName(name);
  if (!match) return false;

  // Don't allow dragging solved segment endpoints
  const chain = App.state.outlineChain;
  const n = chain.length;
  if (match.seq === 0 || match.seq === n - 2) return false;

  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);

  OutlineEditor.pending = true;
  OutlineEditor.startScreen = { x: e.clientX, y: e.clientY };
  OutlineEditor.startWorld = [wx, wy];
  OutlineEditor.dragPoint = {
    name: name,
    seq: match.seq,
    segType: match.seg.seg_type,
    origWorld: [parseFloat(target.getAttribute("cx")),
                -parseFloat(target.getAttribute("cy"))],
  };

  e.preventDefault();
  e.stopPropagation();
  return true;
}


/**
 * Called from onMouseMove during outline drag.
 */
function outlineEditorMouseMove(e) {
  if (!OutlineEditor.pending && !OutlineEditor.active) return;

  const dx = e.clientX - OutlineEditor.startScreen.x;
  const dy = e.clientY - OutlineEditor.startScreen.y;
  const dist = Math.sqrt(dx * dx + dy * dy);

  if (!OutlineEditor.active && dist < 4) return; // threshold

  if (!OutlineEditor.active) {
    OutlineEditor.active = true;
    // Create ghost circle
    const orig = OutlineEditor.dragPoint;
    OutlineEditor.ghost = svgEl("circle", {
      cx: orig.origWorld[0],
      cy: -orig.origWorld[1],
      r: 0.15,
      class: "move-ghost",
      fill: "var(--accent)",
      stroke: "var(--accent)",
      "stroke-width": 0.02,
      opacity: 0.6,
    });
    App.els["layer-points"].appendChild(OutlineEditor.ghost);
  }

  // Update ghost position
  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
  OutlineEditor.ghost.setAttribute("cx", wx);
  OutlineEditor.ghost.setAttribute("cy", -wy);
}


/**
 * Called from onMouseUp to finalize outline drag.
 */
async function outlineEditorMouseUp(e) {
  if (!OutlineEditor.pending && !OutlineEditor.active) return;

  const wasActive = OutlineEditor.active;

  // Clean up ghost
  if (OutlineEditor.ghost) {
    OutlineEditor.ghost.remove();
    OutlineEditor.ghost = null;
  }

  OutlineEditor.pending = false;
  OutlineEditor.active = false;

  if (!wasActive) return; // just a click, not a drag

  const rect = App.els["viewport"].getBoundingClientRect();
  const [wx, wy] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
  const orig = OutlineEditor.dragPoint;
  const dragDx = wx - orig.origWorld[0];
  const dragDy = wy - orig.origWorld[1];
  const dragDist = Math.sqrt(dragDx * dragDx + dragDy * dragDy);

  if (dragDist < 0.01) return; // negligible drag

  const seg = findSegByEndName(orig.name);
  if (!seg) return;

  OutlineEditor._committing = true;

  try {
    if (orig.segType === "L") {
      // For lines: compute new distance from drag projection
      // The total drag distance along the segment direction
      const oldDist = seg.seg.distance || 0;
      const newDist = oldDist + dragDist * Math.sign(dragDx || dragDy);

      // Use a simpler approach: just set distance to whatever makes the
      // endpoint land at the dragged position
      await apiFetch(`/api/outline/${orig.seq}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ dist_or_radius: Math.max(0.1, oldDist + dragDist * Math.sign(dragDx + dragDy)) }),
      });
    } else {
      // For arcs: adjust radius based on perpendicular drag distance
      const oldR = seg.seg.radius || 1;
      const newR = Math.max(0.1, oldR + dragDist * Math.sign(dragDx + dragDy));
      await apiFetch(`/api/outline/${orig.seq}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ dist_or_radius: newR }),
      });
    }
    showToast(`Adjusted ${orig.name}`, "success");
  } catch (err) {
    showToast(`Error: ${err.message}`, "error");
  } finally {
    OutlineEditor._committing = false;
  }
}
