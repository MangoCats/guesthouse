/* ========== Dialog Framework ========== */
"use strict";

const Dialog = {
  _overlay: null,

  /**
   * Show a modal dialog.
   * @param {Object} opts
   * @param {string} opts.title - Dialog title
   * @param {Array<{label: string, name: string, value?: string, placeholder?: string}>} opts.fields
   * @param {Function} opts.onSubmit - Called with {name: value} dict
   * @param {Function} [opts.onCancel]
   */
  show(opts) {
    this.close();

    const overlay = document.createElement("div");
    overlay.className = "dialog-overlay";

    const dialog = document.createElement("div");
    dialog.className = "dialog";

    const title = document.createElement("h3");
    title.textContent = opts.title;
    dialog.appendChild(title);

    // Custom content (e.g. catalog grid)
    if (opts.customContent) {
      dialog.appendChild(opts.customContent);
    }

    const inputs = {};
    for (const f of (opts.fields || [])) {
      const field = document.createElement("div");
      field.className = "dialog-field";
      const label = document.createElement("label");
      label.textContent = f.label;
      field.appendChild(label);
      let input;
      if (f.type === "select" && f.options) {
        input = document.createElement("select");
        input.name = f.name;
        for (const opt of f.options) {
          const o = document.createElement("option");
          o.value = typeof opt === "object" ? opt.value : opt;
          o.textContent = typeof opt === "object" ? opt.label : opt;
          if (o.value === f.value) o.selected = true;
          input.appendChild(o);
        }
      } else {
        input = document.createElement("input");
        input.type = "text";
        input.name = f.name;
        input.value = f.value || "";
        input.placeholder = f.placeholder || "";
      }
      field.appendChild(input);
      dialog.appendChild(field);
      inputs[f.name] = input;
    }

    // Preset buttons (e.g. rotation presets)
    if (opts.presetButtons && inputs[opts.presetButtons.target]) {
      const presetDiv = document.createElement("div");
      presetDiv.className = "dialog-presets";
      for (const preset of opts.presetButtons.values) {
        const btn = document.createElement("button");
        btn.className = "dialog-preset-btn";
        btn.textContent = preset.label;
        btn.onclick = () => { inputs[opts.presetButtons.target].value = preset.value; };
        presetDiv.appendChild(btn);
      }
      dialog.appendChild(presetDiv);
    }

    const buttons = document.createElement("div");
    buttons.className = "dialog-buttons";

    const cancelBtn = document.createElement("button");
    cancelBtn.className = "dialog-btn-cancel";
    cancelBtn.textContent = "Cancel";
    cancelBtn.onclick = () => {
      this.close();
      if (opts.onCancel) opts.onCancel();
    };
    buttons.appendChild(cancelBtn);

    const okBtn = document.createElement("button");
    okBtn.className = "dialog-btn-primary";
    okBtn.textContent = "OK";
    okBtn.onclick = () => {
      const values = {};
      for (const [name, input] of Object.entries(inputs)) {
        values[name] = input.value;
      }
      this.close();
      opts.onSubmit(values);
    };
    buttons.appendChild(okBtn);

    dialog.appendChild(buttons);
    overlay.appendChild(dialog);
    document.body.appendChild(overlay);
    this._overlay = overlay;

    // Focus first input
    const firstInput = dialog.querySelector("input");
    if (firstInput) firstInput.focus();

    // Keyboard: Enter submits, Escape closes
    overlay.addEventListener("keydown", (e) => {
      if (e.key === "Enter") {
        e.preventDefault();
        okBtn.click();
      } else if (e.key === "Escape") {
        e.preventDefault();
        cancelBtn.click();
      }
    });

    // Click overlay background to cancel
    overlay.addEventListener("click", (e) => {
      if (e.target === overlay) cancelBtn.click();
    });
  },

  close() {
    if (this._overlay) {
      this._overlay.remove();
      this._overlay = null;
    }
  },
};


/**
 * Parse an offset string like "6in east", "2ft north", "3.5 south".
 * Returns {dx, dy} in feet, or null if unparseable.
 */
function parseOffsetString(str) {
  str = str.trim().toLowerCase();
  if (!str) return null;

  // Match number + optional unit + direction
  const m = str.match(/^(-?[\d.]+)\s*(in|inch|inches|ft|feet|foot|')?\s*(north|south|east|west|n|s|e|w)$/);
  if (!m) return null;

  let val = parseFloat(m[1]);
  if (isNaN(val)) return null;

  const unit = m[2];
  if (unit && (unit === "in" || unit === "inch" || unit === "inches")) {
    val /= 12.0;
  }
  // ft/feet/foot/' or no unit → already in feet

  const dir = m[3];
  switch (dir) {
    case "north": case "n": return { dx: 0, dy: val };
    case "south": case "s": return { dx: 0, dy: -val };
    case "east":  case "e": return { dx: val, dy: 0 };
    case "west":  case "w": return { dx: -val, dy: 0 };
    default: return null;
  }
}
