/* ========== Shape Editor (TL-25, TL-26, TL-27) ========== */
"use strict";

const ShapeEditor = {
  shapes: [],
};


/** Load shapes from API into ShapeEditor.shapes. */
async function loadShapes() {
  try {
    const resp = await fetch("/api/shapes");
    ShapeEditor.shapes = await resp.json();
  } catch (e) {
    console.error("Failed to load shapes:", e);
  }
}


/**
 * Show the shape editor dialog (TL-25).
 * @param {string} [shapeName] - If given, edit this shape; otherwise create new.
 * @param {Function} [onSave]  - Callback when shape is saved.
 */
function showShapeEditor(shapeName, onSave) {
  let vertices = [];
  let name = shapeName || "";
  let description = "";

  // Load existing shape data
  if (shapeName) {
    const existing = ShapeEditor.shapes.find(s => s.name === shapeName);
    if (existing) {
      const parsed = typeof existing.poly_json === "string"
        ? JSON.parse(existing.poly_json) : existing.poly_json;
      vertices = parsed.map(v => [...v]);
      name = existing.name;
      description = existing.description || "";
    }
  }

  // Default triangle if no vertices
  if (vertices.length === 0) {
    vertices = [[0, 1], [1, -1], [-1, -1]];
  }

  // Build dialog content
  const content = document.createElement("div");
  content.className = "shape-editor-dialog";

  // SVG canvas for interactive editing
  const svgNS = "http://www.w3.org/2000/svg";
  const svgCanvas = document.createElementNS(svgNS, "svg");
  svgCanvas.setAttribute("viewBox", "-3 -3 6 6");
  svgCanvas.setAttribute("width", "280");
  svgCanvas.setAttribute("height", "280");
  svgCanvas.style.background = "var(--surface0)";
  svgCanvas.style.borderRadius = "6px";
  svgCanvas.style.display = "block";
  svgCanvas.style.margin = "0 auto 8px";
  content.appendChild(svgCanvas);

  // Polygon element
  const polyEl = document.createElementNS(svgNS, "polygon");
  polyEl.setAttribute("fill", "rgba(137,180,250,0.3)");
  polyEl.setAttribute("stroke", "var(--accent)");
  polyEl.setAttribute("stroke-width", "0.06");
  svgCanvas.appendChild(polyEl);

  // Vertex handles container
  const handlesGroup = document.createElementNS(svgNS, "g");
  svgCanvas.appendChild(handlesGroup);

  // Vertex list
  const vertList = document.createElement("div");
  vertList.style.maxHeight = "120px";
  vertList.style.overflowY = "auto";
  vertList.style.fontSize = "0.8rem";
  vertList.style.marginBottom = "6px";
  content.appendChild(vertList);

  // Buttons row
  const btnRow = document.createElement("div");
  btnRow.style.display = "flex";
  btnRow.style.gap = "6px";
  btnRow.style.marginBottom = "6px";

  const addBtn = document.createElement("button");
  addBtn.className = "dialog-preset-btn";
  addBtn.textContent = "+ Vertex";
  addBtn.onclick = () => {
    // Add vertex at centroid
    const cx = vertices.reduce((s, v) => s + v[0], 0) / vertices.length;
    const cy = vertices.reduce((s, v) => s + v[1], 0) / vertices.length;
    vertices.push([cx, cy]);
    updateEditor();
  };
  btnRow.appendChild(addBtn);

  const removeBtn = document.createElement("button");
  removeBtn.className = "dialog-preset-btn";
  removeBtn.textContent = "- Vertex";
  removeBtn.onclick = () => {
    if (vertices.length > 3) {
      vertices.pop();
      updateEditor();
    }
  };
  btnRow.appendChild(removeBtn);

  // TL-27: SVG import button
  const importBtn = document.createElement("button");
  importBtn.className = "dialog-preset-btn";
  importBtn.textContent = "Import SVG";
  importBtn.onclick = () => {
    const input = document.createElement("textarea");
    input.placeholder = '<polygon points="..."/> or <path d="M...L..."/>';
    input.style.width = "100%";
    input.style.height = "60px";
    input.style.fontFamily = "monospace";
    input.style.fontSize = "0.75rem";

    Dialog.show({
      title: "Import SVG",
      customContent: input,
      fields: [],
      onSubmit() {
        const parsed = parseSvgPolygon(input.value);
        if (parsed && parsed.length >= 3) {
          vertices = parsed;
          // Reopen the shape editor with imported vertices
          Dialog.close();
          showShapeEditor(shapeName, onSave);
          // Vertices will be set from the re-opened dialog
        } else {
          showToast("Could not parse SVG polygon", "error");
        }
      },
    });
  };
  btnRow.appendChild(importBtn);
  content.appendChild(btnRow);

  function updateEditor() {
    // Update polygon
    const pts = vertices.map(v => `${v[0]},${-v[1]}`).join(" ");
    polyEl.setAttribute("points", pts);

    // Update handles
    while (handlesGroup.firstChild) handlesGroup.removeChild(handlesGroup.firstChild);
    vertices.forEach((v, i) => {
      const circle = document.createElementNS(svgNS, "circle");
      circle.setAttribute("cx", v[0]);
      circle.setAttribute("cy", -v[1]);
      circle.setAttribute("r", "0.15");
      circle.setAttribute("fill", "var(--accent)");
      circle.setAttribute("stroke", "var(--surface0)");
      circle.setAttribute("stroke-width", "0.03");
      circle.style.cursor = "grab";

      // Drag handling
      let dragging = false;
      circle.addEventListener("mousedown", (e) => {
        e.preventDefault();
        dragging = true;
        const onMove = (me) => {
          if (!dragging) return;
          const pt = svgCanvas.createSVGPoint();
          pt.x = me.clientX;
          pt.y = me.clientY;
          const svgPt = pt.matrixTransform(svgCanvas.getScreenCTM().inverse());
          vertices[i] = [Math.round(svgPt.x * 100) / 100, Math.round(-svgPt.y * 100) / 100];
          updateEditor();
        };
        const onUp = () => {
          dragging = false;
          document.removeEventListener("mousemove", onMove);
          document.removeEventListener("mouseup", onUp);
        };
        document.addEventListener("mousemove", onMove);
        document.addEventListener("mouseup", onUp);
      });
      handlesGroup.appendChild(circle);
    });

    // Update vertex list
    vertList.innerHTML = "";
    vertices.forEach((v, i) => {
      const row = document.createElement("div");
      row.textContent = `V${i + 1}: (${v[0].toFixed(2)}, ${v[1].toFixed(2)})`;
      row.style.padding = "1px 4px";
      vertList.appendChild(row);
    });
  }

  updateEditor();

  Dialog.show({
    title: shapeName ? `Edit Shape: ${shapeName}` : "New Shape",
    customContent: content,
    fields: [
      { label: "Name", name: "name", value: name },
      { label: "Description", name: "description", value: description },
    ],
    async onSubmit(values) {
      const sName = values.name.trim();
      if (!sName) {
        showToast("Shape name is required", "error");
        return;
      }
      const body = {
        name: sName,
        poly_json: vertices,
        description: values.description || "",
      };

      let resp;
      if (shapeName && ShapeEditor.shapes.find(s => s.name === shapeName)) {
        // Update existing
        resp = await fetch(`/api/shapes/${shapeName}`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(body),
        });
      } else {
        // Create new
        resp = await fetch("/api/shapes", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(body),
        });
      }

      if (resp.ok) {
        showToast(`Shape "${sName}" saved`, "success");
        await loadShapes();
        if (onSave) onSave(sName);
      } else {
        const data = await resp.json();
        showToast(data.error || "Failed to save shape", "error");
      }
    },
  });
}


/**
 * Parse SVG polygon or simple path data into vertices (TL-27).
 * Supports: <polygon points="x1,y1 x2,y2 ..."/>
 *           <path d="M x1 y1 L x2 y2 ... Z"/>
 */
function parseSvgPolygon(text) {
  text = text.trim();

  // Try <polygon points="...">
  const polyMatch = text.match(/points\s*=\s*"([^"]+)"/);
  if (polyMatch) {
    return polyMatch[1].trim().split(/\s+/).map(pair => {
      const parts = pair.split(",").map(Number);
      return [parts[0] || 0, parts[1] || 0];
    }).filter(p => !isNaN(p[0]) && !isNaN(p[1]));
  }

  // Try <path d="M ... L ... Z">
  const pathMatch = text.match(/d\s*=\s*"([^"]+)"/);
  if (pathMatch) {
    const d = pathMatch[1];
    const points = [];
    // Match M/L followed by coordinates
    const re = /[ML]\s*(-?[\d.]+)[,\s]+(-?[\d.]+)/gi;
    let m;
    while ((m = re.exec(d)) !== null) {
      points.push([parseFloat(m[1]), parseFloat(m[2])]);
    }
    if (points.length >= 3) return points;
  }

  return null;
}


/**
 * Show a shape picker dropdown for element properties (TL-26).
 * @param {HTMLElement} tbody - The properties table body.
 * @param {Object} elemRec - The element DB record.
 * @param {Object} props - The element's parsed properties.
 */
function addShapePicker(tbody, elemRec, props) {
  const tr = document.createElement("tr");
  const td1 = document.createElement("td");
  td1.textContent = "Shape";
  tr.appendChild(td1);

  const td2 = document.createElement("td");
  const select = document.createElement("select");
  select.className = "prop-edit-input";

  // Add "rect" option
  const rectOpt = document.createElement("option");
  rectOpt.value = "rect";
  rectOpt.textContent = "Rectangle";
  select.appendChild(rectOpt);

  // Add shapes from DB
  for (const s of ShapeEditor.shapes) {
    const opt = document.createElement("option");
    opt.value = s.name;
    opt.textContent = s.name;
    select.appendChild(opt);
  }

  // Select current value
  select.value = props.shape || "rect";

  select.addEventListener("change", async () => {
    const shapeName = select.value;
    const newProps = { ...props, shape: shapeName };

    if (shapeName !== "rect") {
      // Load shape poly and transform to element position/rotation
      const shape = ShapeEditor.shapes.find(s => s.name === shapeName);
      if (shape) {
        const shapePoly = typeof shape.poly_json === "string"
          ? JSON.parse(shape.poly_json) : shape.poly_json;
        const cx = props.center ? props.center[0] : 0;
        const cy = props.center ? props.center[1] : 0;
        const w = props.width || 1;
        const d = props.depth || 1;
        const angle = props.rotation || 0;
        const rad = angle * Math.PI / 180;
        const cos = Math.cos(rad), sin = Math.sin(rad);
        // Scale shape poly to element dimensions and transform
        newProps.poly = shapePoly.map(([sx, sy]) => {
          const lx = sx * w / 2, ly = sy * d / 2;
          return [
            cx + lx * cos - ly * sin,
            cy + lx * sin + ly * cos,
          ];
        });
      }
    } else if (props.center && props.width && props.depth) {
      newProps.poly = rotatedRectPoly(
        props.center[0], props.center[1],
        props.width, props.depth,
        props.rotation || 0
      );
    }

    await fetch(`/api/elements/${elemRec.id}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ properties: newProps }),
    });
    await loadElements();
    await loadGeometry();
  });

  // Edit shape button
  const editBtn = document.createElement("button");
  editBtn.className = "dialog-preset-btn";
  editBtn.textContent = "Edit";
  editBtn.style.marginLeft = "4px";
  editBtn.onclick = () => {
    const currentShape = select.value !== "rect" ? select.value : null;
    showShapeEditor(currentShape);
  };

  td2.appendChild(select);
  td2.appendChild(editBtn);
  tr.appendChild(td2);
  tbody.appendChild(tr);
}
