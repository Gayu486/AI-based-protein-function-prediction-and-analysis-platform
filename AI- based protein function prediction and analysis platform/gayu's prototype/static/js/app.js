/**
 * Main Application Frontend Logic
 * AI-BASED PROTEIN FUNCTION PREDICTION AND ANALYSIS PLATFORM
 */

let currentAnalysisData = null;
let sampleProteinsList = [];
let validateDebounceTimer = null;

document.addEventListener("DOMContentLoaded", () => {
  setupNavigation();
  loadSampleProteins();
  setupInputListeners();
  setupFileUpload();
  setupReportDownloads();
});

// ----------------------------------------------------
// TAB NAVIGATION
// ----------------------------------------------------
function setupNavigation() {
  const navLinks = document.querySelectorAll(".nav-link[data-tab]");
  navLinks.forEach(link => {
    link.addEventListener("click", (e) => {
      e.preventDefault();
      const tabId = link.getAttribute("data-tab");
      switchToTab(tabId);
    });
  });

  // Action buttons that switch tabs
  document.querySelectorAll("[data-action='start-analysis']").forEach(btn => {
    btn.addEventListener("click", () => switchToTab("analysis"));
  });

  document.querySelectorAll("[data-action='try-example']").forEach(btn => {
    btn.addEventListener("click", () => {
      switchToTab("analysis");
      if (sampleProteinsList.length > 0) {
        loadSampleIntoInput(sampleProteinsList[0]);
      }
    });
  });
}

function switchToTab(tabId) {
  // Update nav links
  document.querySelectorAll(".nav-link").forEach(link => {
    link.classList.toggle("active", link.getAttribute("data-tab") === tabId);
  });

  // Update tab panes
  document.querySelectorAll(".tab-pane").forEach(pane => {
    pane.classList.toggle("active", pane.id === tabId);
  });

  window.scrollTo({ top: 0, behavior: "smooth" });
}

// ----------------------------------------------------
// SAMPLE PROTEINS PRESET LOADING
// ----------------------------------------------------
async function loadSampleProteins() {
  try {
    const res = await fetch("/api/samples");
    const data = await res.json();
    if (data.success && data.samples) {
      sampleProteinsList = data.samples;
      renderSampleList(sampleProteinsList);
    }
  } catch (err) {
    console.warn("Could not fetch preset sample proteins:", err);
  }
}

function renderSampleList(samples) {
  const container = document.getElementById("samples-container");
  if (!container) return;

  container.innerHTML = "";
  samples.forEach(sample => {
    const item = document.createElement("div");
    item.className = "sample-item";
    item.innerHTML = `
      <div class="sample-item-header">
        <span class="sample-name">${sample.name}</span>
        <span class="sample-acc">${sample.accession}</span>
      </div>
      <div class="sample-cat">${sample.category}</div>
      <div class="sample-desc">${sample.description}</div>
    `;
    item.addEventListener("click", () => {
      loadSampleIntoInput(sample);
    });
    container.appendChild(item);
  });
}

function loadSampleIntoInput(sample) {
  const textarea = document.getElementById("sequence-input");
  const accInput = document.getElementById("accession-input");
  
  if (textarea) {
    textarea.value = `>${sample.accession} ${sample.name}\n${sample.sequence}`;
    updateCharCounter();
    triggerLiveValidation(textarea.value);
  }
  if (accInput) {
    accInput.value = sample.accession;
  }
}

// ----------------------------------------------------
// INPUT & VALIDATION LISTENERS
// ----------------------------------------------------
function setupInputListeners() {
  const textarea = document.getElementById("sequence-input");
  const charCounter = document.getElementById("char-counter");
  const clearBtn = document.getElementById("clear-btn");
  const analyzeBtn = document.getElementById("analyze-btn");
  const accInput = document.getElementById("accession-input");
  const accSearchBtn = document.getElementById("acc-search-btn");

  const tabBtnSeq = document.getElementById("tab-btn-seq");
  const tabBtnAcc = document.getElementById("tab-btn-acc");
  const accBlock = document.getElementById("accession-lookup-block");
  const seqWrap = document.getElementById("seq-textarea-wrap");

  if (tabBtnSeq && tabBtnAcc) {
    tabBtnSeq.addEventListener("click", () => {
      tabBtnSeq.classList.add("active");
      tabBtnAcc.classList.remove("active");
      if (accBlock) accBlock.style.display = "none";
      if (seqWrap) seqWrap.style.display = "block";
    });

    tabBtnAcc.addEventListener("click", () => {
      tabBtnAcc.classList.add("active");
      tabBtnSeq.classList.remove("active");
      if (accBlock) accBlock.style.display = "block";
      if (seqWrap) seqWrap.style.display = "none";
    });
  }

  if (textarea) {
    textarea.addEventListener("input", () => {
      updateCharCounter();
      clearTimeout(validateDebounceTimer);
      validateDebounceTimer = setTimeout(() => {
        triggerLiveValidation(textarea.value);
      }, 350);
    });
  }

  if (clearBtn) {
    clearBtn.addEventListener("click", () => {
      if (textarea) textarea.value = "";
      if (accInput) accInput.value = "";
      updateCharCounter();
      hideValidationBanner();
    });
  }

  if (analyzeBtn) {
    analyzeBtn.addEventListener("click", () => {
      const isAccActive = tabBtnAcc && tabBtnAcc.classList.contains("active");
      const seq = textarea ? textarea.value.trim() : "";
      const acc = accInput ? accInput.value.trim() : "";

      if (isAccActive) {
        if (!acc) {
          showValidationBanner("invalid", "Please specify a UniProt accession code (e.g. P68871).");
          return;
        }
        runProteinAnalysis({ accession: acc });
      } else {
        if (!seq) {
          showValidationBanner("invalid", "Please enter an amino acid sequence to analyze.");
          return;
        }
        runProteinAnalysis({ sequence: seq });
      }
    });
  }

  if (accSearchBtn) {
    accSearchBtn.addEventListener("click", () => {
      const acc = accInput ? accInput.value.trim() : "";
      if (!acc) {
        showValidationBanner("invalid", "Please specify a UniProt accession code (e.g. P68871).");
        return;
      }
      runProteinAnalysis({ accession: acc });
    });
  }
}

function updateCharCounter() {
  const textarea = document.getElementById("sequence-input");
  const counter = document.getElementById("char-counter");
  if (!textarea || !counter) return;

  const raw = textarea.value;
  // Estimate length without fasta header and whitespace
  const lines = raw.split("\n");
  let len = 0;
  for (const l of lines) {
    if (!l.startsWith(">")) {
      len += l.replace(/[^A-Za-z]/g, "").length;
    }
  }
  counter.textContent = `${len} AA`;
}

async function triggerLiveValidation(rawText) {
  if (!rawText.trim()) {
    hideValidationBanner();
    return;
  }

  try {
    const res = await fetch("/api/validate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sequence: rawText })
    });
    const val = await res.json();
    renderValidationFeedback(val);
  } catch (err) {
    console.error("Validation request error:", err);
  }
}

function renderValidationFeedback(val) {
  if (!val.is_valid) {
    const msg = val.errors.join(" ");
    showValidationBanner("invalid", `Validation Error: ${msg}`);
  } else if (val.warnings && val.warnings.length > 0) {
    const msg = val.warnings.join(" ");
    showValidationBanner("warning", `Sequence Valid with Notices (${val.length} AA): ${msg}`);
  } else {
    showValidationBanner("valid", `Sequence format is valid. Ready for analysis (${val.length} standard amino acids).`);
  }
}

function showValidationBanner(type, message) {
  const banner = document.getElementById("validation-banner");
  if (!banner) return;
  banner.className = `validation-banner ${type}`;
  banner.textContent = message;
  banner.style.display = "block";
}

function hideValidationBanner() {
  const banner = document.getElementById("validation-banner");
  if (banner) banner.style.display = "none";
}

// ----------------------------------------------------
// FASTA FILE UPLOAD
// ----------------------------------------------------
function setupFileUpload() {
  const fileInput = document.getElementById("fasta-file-input");
  const fileDropZone = document.getElementById("file-drop-zone");

  if (!fileInput) return;

  fileInput.addEventListener("change", (e) => {
    const file = e.target.files[0];
    if (file) handleFastaFile(file);
  });

  if (fileDropZone) {
    fileDropZone.addEventListener("dragover", (e) => {
      e.preventDefault();
      fileDropZone.classList.add("dragover");
    });
    fileDropZone.addEventListener("dragleave", () => {
      fileDropZone.classList.remove("dragover");
    });
    fileDropZone.addEventListener("drop", (e) => {
      e.preventDefault();
      fileDropZone.classList.remove("dragover");
      if (e.dataTransfer.files && e.dataTransfer.files[0]) {
        handleFastaFile(e.dataTransfer.files[0]);
      }
    });
  }
}

function handleFastaFile(file) {
  const reader = new FileReader();
  reader.onload = (e) => {
    const content = e.target.result;
    const textarea = document.getElementById("sequence-input");
    if (textarea) {
      textarea.value = content;
      updateCharCounter();
      triggerLiveValidation(content);
    }
  };
  reader.readAsText(file);
}

// ----------------------------------------------------
// ANALYSIS EXECUTION & LOADING SIMULATION
// ----------------------------------------------------
async function runProteinAnalysis(payload) {
  showLoadingModal();

  const stages = [
    { pct: 15, msg: "Validating sequence syntax and headers..." },
    { pct: 35, msg: "Computing physicochemical properties (MW, pI, GRAVY)..." },
    { pct: 60, msg: "Querying biological databases & pairwise sequence alignment..." },
    { pct: 80, msg: "Scanning transmembrane domains & signal peptides..." },
    { pct: 95, msg: "Running Random Forest function & localization classifiers..." },
  ];

  let currentStage = 0;
  const stageInterval = setInterval(() => {
    if (currentStage < stages.length) {
      updateLoadingStatus(stages[currentStage].pct, stages[currentStage].msg);
      currentStage++;
    }
  }, 400);

  try {
    const response = await fetch("/api/analyze", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    });

    clearInterval(stageInterval);
    const data = await response.json();

    if (!data.success) {
      hideLoadingModal();
      showValidationBanner("invalid", data.error || "Analysis failed.");
      alert(`Analysis notice: ${data.error || "Unknown server response."}`);
      return;
    }

    updateLoadingStatus(100, "Rendering results dashboard...");
    currentAnalysisData = data;

    setTimeout(() => {
      hideLoadingModal();
      populateResultsDashboard(data);
      switchToTab("results");
    }, 450);

  } catch (err) {
    clearInterval(stageInterval);
    hideLoadingModal();
    console.error("Network or execution failure:", err);
    alert("Connection error: The analysis service is temporarily unavailable or timed out. Please verify your connection and retry.");
  }
}

function showLoadingModal() {
  const modal = document.getElementById("loading-modal");
  if (modal) modal.classList.add("active");
  updateLoadingStatus(5, "Initiating biological analysis engine...");
}

function hideLoadingModal() {
  const modal = document.getElementById("loading-modal");
  if (modal) modal.classList.remove("active");
}

function updateLoadingStatus(pct, msg) {
  const fill = document.getElementById("loading-progress-fill");
  const text = document.getElementById("loading-status-text");
  if (fill) fill.style.width = `${pct}%`;
  if (text) text.textContent = msg;
}

// ----------------------------------------------------
// POPULATE RESULTS DASHBOARD
// ----------------------------------------------------
function populateResultsDashboard(data) {
  const seq = data.sequence || "";
  const phys = data.physicochemical || {};
  const ident = data.identification || {};
  const rec = ident.record;
  const aln = ident.alignment;
  const pred = data.predictions || {};

  // 1. Header Title & Badges
  const titleElem = document.getElementById("results-protein-title");
  const badgesContainer = document.getElementById("results-header-badges");
  
  if (titleElem) {
    titleElem.textContent = rec ? `${rec.protein_name}` : "Analyzed Protein Sequence";
  }

  if (badgesContainer) {
    badgesContainer.innerHTML = "";
    if (rec) {
      const revBadge = rec.reviewed ? 
        `<span class="badge badge-reviewed">Reviewed (Swiss-Prot)</span>` : 
        `<span class="badge badge-unreviewed">Unreviewed (TrEMBL)</span>`;
      badgesContainer.innerHTML += revBadge;
      badgesContainer.innerHTML += `<span class="badge badge-blue">${rec.accession}</span>`;
      badgesContainer.innerHTML += `<span class="badge badge-teal">${rec.organism}</span>`;
    } else {
      badgesContainer.innerHTML += `<span class="badge badge-purple">Custom / Unannotated Sequence</span>`;
      badgesContainer.innerHTML += `<span class="badge badge-blue">${seq.length} AA</span>`;
    }
  }

  // 2. Metrics Bar
  document.getElementById("res-metric-length").textContent = `${seq.length} AA`;
  document.getElementById("res-metric-mw").textContent = `${phys.molecular_weight_kda || '0'} kDa`;
  document.getElementById("res-metric-pi").textContent = `${phys.isoelectric_point || '0'}`;
  document.getElementById("res-metric-charge").textContent = `${phys.net_charge_ph7 > 0 ? '+' : ''}${phys.net_charge_ph7}`;
  document.getElementById("res-metric-gravy").textContent = `${phys.gravy || '0'}`;
  document.getElementById("res-metric-stability").textContent = `${phys.instability_index || '0'}`;
  document.getElementById("res-metric-stability-desc").textContent = phys.stability_class || '';

  // 3. Identification & Alignment Card
  const identCard = document.getElementById("identification-section");
  if (identCard) {
    if (rec) {
      identCard.style.display = "block";
      document.getElementById("ident-protein-name").textContent = rec.protein_name;
      document.getElementById("ident-entry-name").textContent = rec.entry_name;
      document.getElementById("ident-accession").innerHTML = `<a href="https://www.uniprot.org/uniprotkb/${rec.accession}" target="_blank" class="nav-link" style="padding:0; color:#2563eb; text-decoration:underline;">${rec.accession}</a>`;
      document.getElementById("ident-gene").textContent = rec.gene_name || 'Not specified';
      document.getElementById("ident-organism").textContent = rec.organism;
      document.getElementById("ident-source").textContent = ident.database_source || 'UniProtKB';

      // Alignment statistics
      if (aln) {
        document.getElementById("aln-identity-val").textContent = `${aln.identity_pct}%`;
        document.getElementById("aln-identity-sub").textContent = `${aln.identical_count} / ${aln.alignment_length} positions`;
        document.getElementById("aln-coverage-val").textContent = `${aln.coverage_pct}%`;
        document.getElementById("aln-evalue-val").textContent = aln.e_value_str;

        // Render alignment snippet
        renderAlignmentSnippet("alignment-view-container", aln.blocks);
      }
    } else {
      identCard.style.display = "block";
      document.getElementById("ident-protein-name").textContent = "No identical database match";
      document.getElementById("ident-accession").textContent = "N/A";
      document.getElementById("ident-organism").textContent = "Unknown";
      document.getElementById("aln-identity-val").textContent = "N/A";
      document.getElementById("aln-coverage-val").textContent = "N/A";
      document.getElementById("aln-evalue-val").textContent = "N/A";
      document.getElementById("alignment-view-container").innerHTML = `<p style="padding:15px; color:#94a3b8;">No alignment available (sequence identity with Swiss-Prot reference set < 25%).</p>`;
    }
  }

  // 4. Functional Annotations & GO
  const funcCard = document.getElementById("function-section");
  if (funcCard) {
    if (rec && rec.function) {
      funcCard.style.display = "block";
      document.getElementById("annotated-function-text").textContent = rec.function;

      const goMF = document.getElementById("go-mf-list");
      const goBP = document.getElementById("go-bp-list");
      const goCC = document.getElementById("go-cc-list");

      if (goMF) renderGOList(goMF, rec.go_terms ? rec.go_terms.molecular_function : []);
      if (goBP) renderGOList(goBP, rec.go_terms ? rec.go_terms.biological_process : []);
      if (goCC) renderGOList(goCC, rec.go_terms ? rec.go_terms.cellular_component : []);
    } else {
      funcCard.style.display = "block";
      document.getElementById("annotated-function-text").textContent = "No direct functional description retrieved from UniProtKB for this exact sequence.";
      document.getElementById("go-mf-list").innerHTML = "<li>None annotated</li>";
      document.getElementById("go-bp-list").innerHTML = "<li>None annotated</li>";
      document.getElementById("go-cc-list").innerHTML = "<li>None annotated</li>";
    }
  }

  // 5. Protein Domains & Families
  const domainCard = document.getElementById("domain-section");
  if (domainCard) {
    const domainTableBody = document.getElementById("domains-table-body");
    const sitesContainer = document.getElementById("functional-sites-container");

    if (rec && (rec.domains.length > 0 || rec.functional_sites.length > 0)) {
      domainCard.style.display = "block";

      if (domainTableBody) {
        domainTableBody.innerHTML = "";
        rec.domains.forEach(d => {
          const row = document.createElement("tr");
          row.innerHTML = `
            <td><strong>${d.name}</strong></td>
            <td><a href="https://www.ebi.ac.uk/interpro/" target="_blank" style="color:#2563eb;">${d.database}</a></td>
            <td>${d.start > 0 ? `${d.start} - ${d.end}` : 'N/A'}</td>
          `;
          domainTableBody.appendChild(row);
        });
      }

      if (sitesContainer) {
        sitesContainer.innerHTML = "";
        rec.functional_sites.forEach(s => {
          const badge = document.createElement("div");
          badge.className = "sample-item";
          badge.style.cursor = "default";
          badge.innerHTML = `
            <div style="font-weight:700; color:#0f172a; font-size:12px;">${s.type} (Pos: ${s.position})</div>
            <div style="font-size:11px; color:#64748b;">${s.description}</div>
          `;
          sitesContainer.appendChild(badge);
        });
      }
    } else {
      domainCard.style.display = "block";
      if (domainTableBody) {
        domainTableBody.innerHTML = `<tr><td colspan="3" style="text-align:center; color:#94a3b8; padding:16px;">No domain annotations retrieved. <a href="https://www.ebi.ac.uk/interpro/search/sequence/" target="_blank" style="color:#2563eb;">Search InterPro</a></td></tr>`;
      }
      if (sitesContainer) sitesContainer.innerHTML = "<p style='color:#94a3b8; font-size:12px;'>No active or functional sites cataloged.</p>";
    }
  }

  // 6. Machine Learning Function Predictions
  renderPredictionsSection(pred);

  // 7. Composition Table & Charts
  renderCompositionTable(phys.composition || {});
  renderAminoAcidCompositionChart("aa-chart-canvas", phys.composition || {});
  renderBiochemicalGroupChart("biogroup-chart-canvas", phys.biochemical_groups || {});
  renderHydropathyProfileChart("hydropathy-chart-canvas", phys.hydropathy_profile || {});

  // 8. 3D Structure Viewer
  if (rec && (rec.pdb_ids && rec.pdb_ids.length > 0)) {
    init3DViewer("viewer3d-container", rec.pdb_ids[0], true);
  } else if (rec && rec.alphafold_id) {
    init3DViewer("viewer3d-container", rec.alphafold_id, false);
  } else {
    init3DViewer("viewer3d-container", "1HHO", true);
  }
}

function renderAlignmentSnippet(containerId, blocks) {
  const container = document.getElementById(containerId);
  if (!container) return;

  if (!blocks || blocks.length === 0) {
    container.innerHTML = `<p style="padding:10px; color:#94a3b8;">No alignment blocks available.</p>`;
    return;
  }

  let text = "";
  blocks.slice(0, 5).forEach(b => {
    text += `Query   ${String(b.query_start).padStart(4)}  ${b.query_str}  ${b.query_end}\n`;
    text += `              ${b.match_str}\n`;
    text += `Sbjct   ${String(b.subject_start).padStart(4)}  ${b.subject_str}  ${b.subject_end}\n\n`;
  });

  container.textContent = text;
}

function renderGOList(container, items) {
  container.innerHTML = "";
  if (!items || items.length === 0) {
    container.innerHTML = "<li>None annotated</li>";
    return;
  }
  items.slice(0, 6).forEach(item => {
    const li = document.createElement("li");
    li.innerHTML = `
      <a href="https://www.ebi.ac.uk/QuickGO/term/${item.id}" target="_blank" style="color:#2563eb; text-decoration:none;">${item.name}</a>
      <span style="font-size:10px; color:#94a3b8; margin-left:4px;">(${item.evidence})</span>
    `;
    container.appendChild(li);
  });
}

function renderCompositionTable(composition) {
  const tbody = document.getElementById("aa-table-body");
  if (!tbody) return;

  tbody.innerHTML = "";
  const sorted = Object.keys(composition).sort();
  sorted.forEach(aa => {
    const c = composition[aa];
    const row = document.createElement("tr");
    row.innerHTML = `
      <td><strong>${c.code1}</strong></td>
      <td>${c.code3}</td>
      <td>${c.name}</td>
      <td>${c.count}</td>
      <td>${c.percentage}%</td>
    `;
    tbody.appendChild(row);
  });
}

function renderPredictionsSection(pred) {
  const enz = pred.enzyme_prediction || {};
  const loc = pred.subcellular_localization || {};
  const struct = pred.structural_motifs || {};

  // Enzyme Box
  const enzTitle = document.getElementById("pred-enzyme-title");
  const enzConfidence = document.getElementById("pred-enzyme-confidence");
  const enzBar = document.getElementById("pred-enzyme-bar");
  const topEcText = document.getElementById("pred-top-ec");

  if (enzTitle) enzTitle.textContent = enz.label || "N/A";
  if (enzConfidence) enzConfidence.textContent = `${enz.confidence || 0}%`;
  if (enzBar) enzBar.style.width = `${enz.confidence || 0}%`;
  if (topEcText) {
    topEcText.textContent = enz.top_ec_class ? `${enz.top_ec_class.name}` : "Non-enzymatic activity";
  }

  // Localization Box
  const locTitle = document.getElementById("pred-loc-title");
  const locConfidence = document.getElementById("pred-loc-confidence");
  const locBar = document.getElementById("pred-loc-bar");
  const locDesc = document.getElementById("pred-loc-desc");

  const topLoc = loc.top_compartment;
  if (locTitle) locTitle.textContent = topLoc ? topLoc.name : "N/A";
  if (locConfidence) locConfidence.textContent = topLoc ? `${topLoc.probability}%` : "0%";
  if (locBar) locBar.style.width = topLoc ? `${topLoc.probability}%` : "0%";
  if (locDesc) locDesc.textContent = topLoc ? topLoc.description : "";

  // Transmembrane & Signal Peptide
  const tmCountElem = document.getElementById("pred-tm-count");
  const signalPeptideElem = document.getElementById("pred-signal-peptide");

  if (tmCountElem) {
    const tmCount = struct.tm_count || 0;
    tmCountElem.textContent = tmCount > 0 ? 
      `${tmCount} Putative Helix/Helices detected` : 
      "No transmembrane segments detected (Soluble)";
  }

  if (signalPeptideElem) {
    const sp = struct.signal_peptide;
    signalPeptideElem.textContent = sp && sp.predicted ? 
      `Predicted (${Math.round(sp.probability * 100)}% likelihood; ${sp.cleavage_site})` : 
      "Negative / Not detected";
  }

  // Explanatory Insights
  const insightsList = document.getElementById("pred-insights-list");
  if (insightsList && pred.insights) {
    insightsList.innerHTML = "";
    pred.insights.forEach(ins => {
      const li = document.createElement("li");
      li.textContent = ins;
      insightsList.appendChild(li);
    });
  }
}

// ----------------------------------------------------
// REPORT DOWNLOAD HANDLERS
// ----------------------------------------------------
function setupReportDownloads() {
  const btnHtml = document.getElementById("download-html-report-btn");
  const btnJson = document.getElementById("download-json-report-btn");
  const btnFasta = document.getElementById("download-fasta-btn");

  if (btnHtml) {
    btnHtml.addEventListener("click", () => downloadReport("/api/report/html", "report.html"));
  }
  if (btnJson) {
    btnJson.addEventListener("click", () => downloadReport("/api/report/json", "analysis.json"));
  }
  if (btnFasta) {
    btnFasta.addEventListener("click", () => downloadReport("/api/report/fasta", "sequence.fasta"));
  }
}

async function downloadReport(endpoint, fallbackFilename) {
  if (!currentAnalysisData) {
    alert("Please perform an analysis before downloading a report.");
    return;
  }

  try {
    const res = await fetch(endpoint, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(currentAnalysisData)
    });

    if (!res.ok) throw new Error("Report export request failed.");

    const blob = await res.blob();
    const disposition = res.headers.get("Content-Disposition");
    let filename = fallbackFilename;
    if (disposition && disposition.indexOf("filename=") !== -1) {
      const match = disposition.match(/filename="?([^"]+)"?/);
      if (match && match[1]) filename = match[1];
    }

    const downloadUrl = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = downloadUrl;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    window.URL.revokeObjectURL(downloadUrl);

  } catch (err) {
    console.error("Download error:", err);
    alert("Failed to generate download file. Please retry.");
  }
}
