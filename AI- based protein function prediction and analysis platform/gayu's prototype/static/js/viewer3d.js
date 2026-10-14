/**
 * 3D Molecular Structure Viewer (Powered by 3Dmol.js)
 */

let glViewer = null;
let currentStructureId = null;
let isSpinning = false;

function init3DViewer(containerId, structureId, isPdb = true) {
  const container = document.getElementById(containerId);
  if (!container) return;

  if (typeof $3Dmol === 'undefined') {
    container.innerHTML = `
      <div style="display:flex; align-items:center; justify-content:center; height:100%; color:#94a3b8; text-align:center; padding:20px;">
        <div>
          <p style="font-weight:600; margin-bottom:8px;">3Dmol.js library loading or offline.</p>
          <a href="https://www.rcsb.org/structure/${structureId}" target="_blank" class="btn btn-outline btn-sm">View on RCSB PDB</a>
        </div>
      </div>
    `;
    return;
  }

  container.innerHTML = "";
  currentStructureId = structureId;

  const config = { backgroundColor: '#0f172a' };
  glViewer = $3Dmol.createViewer(container, config);

  let uri = "";
  if (isPdb) {
    uri = `https://files.rcsb.org/download/${structureId}.pdb`;
  } else {
    // AlphaFold DB model URI
    uri = `https://alphafold.ebi.ac.uk/files/AF-${structureId}-F1-model_v4.pdb`;
  }

  // Load structure via fetch
  fetch(uri)
    .then(res => {
      if (!res.ok) throw new Error("Could not load PDB structure file");
      return res.text();
    })
    .then(pdbData => {
      glViewer.clear();
      glViewer.addModel(pdbData, "pdb");
      glViewer.setStyle({}, { cartoon: { color: 'spectrum' } });
      glViewer.zoomTo();
      glViewer.render();
      glViewer.spin("y", 0.5);
      isSpinning = true;
    })
    .catch(err => {
      console.warn("Could not load remote 3D structure coordinates:", err);
      container.innerHTML = `
        <div style="display:flex; align-items:center; justify-content:center; height:100%; color:#cbd5e1; text-align:center; padding:30px;">
          <div>
            <div style="font-size: 32px; margin-bottom: 12px;">🧬</div>
            <h4 style="font-size:15px; margin-bottom:6px; color:#ffffff;">Direct 3D Coordinates Not Downloaded</h4>
            <p style="font-size:12px; color:#94a3b8; margin-bottom:16px; max-width:320px;">
              Predicted AlphaFold structure or experimental PDB record is available via official repositories.
            </p>
            <div style="display:flex; gap:8px; justify-content:center;">
              <a href="https://alphafold.ebi.ac.uk/entry/${structureId}" target="_blank" class="btn btn-primary btn-sm">AlphaFold DB Model</a>
              <a href="https://www.rcsb.org/structure/${structureId}" target="_blank" class="btn btn-secondary btn-sm">RCSB PDB</a>
            </div>
          </div>
        </div>
      `;
    });
}

function toggleViewerSpin() {
  if (!glViewer) return;
  isSpinning = !isSpinning;
  if (isSpinning) {
    glViewer.spin("y", 0.5);
  } else {
    glViewer.spin(false);
  }
}

function setViewerStyle(styleType) {
  if (!glViewer) return;
  glViewer.setStyle({}, {});
  if (styleType === 'cartoon') {
    glViewer.setStyle({}, { cartoon: { color: 'spectrum' } });
  } else if (styleType === 'sphere') {
    glViewer.setStyle({}, { sphere: { scale: 0.3 } });
  } else if (styleType === 'stick') {
    glViewer.setStyle({}, { stick: {} });
  }
  glViewer.render();
}

function resetViewerView() {
  if (!glViewer) return;
  glViewer.zoomTo();
  glViewer.render();
}
