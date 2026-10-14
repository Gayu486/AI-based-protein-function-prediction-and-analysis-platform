/**
 * Chart.js Visualization Helpers for Protein Analysis
 */

let aaChartInstance = null;
let bioGroupChartInstance = null;
let hydropathyChartInstance = null;

// Color mapping for 20 amino acids based on physicochemical properties
const AA_COLORS = {
  'A': '#3b82f6', 'V': '#3b82f6', 'I': '#3b82f6', 'L': '#3b82f6', 'M': '#3b82f6', // Aliphatic/Hydrophobic
  'F': '#8b5cf6', 'W': '#8b5cf6', 'Y': '#8b5cf6',                                 // Aromatic
  'S': '#10b981', 'T': '#10b981', 'N': '#10b981', 'Q': '#10b981',                 // Polar uncharged
  'K': '#ef4444', 'R': '#ef4444', 'H': '#ef4444',                                 // Basic (+)
  'D': '#f59e0b', 'E': '#f59e0b',                                                 // Acidic (-)
  'G': '#6b7280', 'P': '#6b7280', 'C': '#6b7280'                                  // Special conformation
};

function renderAminoAcidCompositionChart(canvasId, compositionData) {
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;

  if (aaChartInstance) {
    aaChartInstance.destroy();
  }

  const sortedKeys = Object.keys(compositionData).sort();
  const labels = sortedKeys.map(k => `${k} (${compositionData[k].code3})`);
  const counts = sortedKeys.map(k => compositionData[k].count);
  const percentages = sortedKeys.map(k => compositionData[k].percentage);
  const bgColors = sortedKeys.map(k => AA_COLORS[k] || '#94a3b8');

  aaChartInstance = new Chart(ctx, {
    type: 'bar',
    data: {
      labels: labels,
      datasets: [{
        label: 'Residue Percentage (%)',
        data: percentages,
        backgroundColor: bgColors,
        borderRadius: 4,
        borderWidth: 0
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: function(context) {
              const idx = context.dataIndex;
              const aa = sortedKeys[idx];
              return `${compositionData[aa].name}: ${percentages[idx]}% (${counts[idx]} residues)`;
            }
          }
        }
      },
      scales: {
        y: {
          beginAtZero: true,
          title: {
            display: true,
            text: 'Frequency (%)',
            font: { size: 12, weight: 'bold' }
          },
          grid: { color: '#f1f5f9' }
        },
        x: {
          grid: { display: false }
        }
      }
    }
  });
}

function renderBiochemicalGroupChart(canvasId, groupsData) {
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;

  if (bioGroupChartInstance) {
    bioGroupChartInstance.destroy();
  }

  const labels = Object.keys(groupsData);
  const dataVals = labels.map(k => groupsData[k].percentage);
  const colors = labels.map(k => groupsData[k].color);

  bioGroupChartInstance = new Chart(ctx, {
    type: 'doughnut',
    data: {
      labels: labels,
      datasets: [{
        data: dataVals,
        backgroundColor: colors,
        borderWidth: 2,
        borderColor: '#ffffff'
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          position: 'bottom',
          labels: {
            boxWidth: 12,
            font: { size: 11 }
          }
        },
        tooltip: {
          callbacks: {
            label: function(context) {
              return ` ${context.label}: ${context.raw}%`;
            }
          }
        }
      },
      cutout: '65%'
    }
  });
}

function renderHydropathyProfileChart(canvasId, profileData) {
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;

  if (hydropathyChartInstance) {
    hydropathyChartInstance.destroy();
  }

  const positions = profileData.positions || [];
  const scores = profileData.scores || [];

  hydropathyChartInstance = new Chart(ctx, {
    type: 'line',
    data: {
      labels: positions,
      datasets: [{
        label: `Kyte-Doolittle Hydropathy (Window=${profileData.window_size || 9})`,
        data: scores,
        borderColor: '#2563eb',
        backgroundColor: 'rgba(37, 99, 235, 0.08)',
        fill: true,
        tension: 0.2,
        borderWidth: 2,
        pointRadius: scores.length > 100 ? 0 : 2,
        pointHoverRadius: 4
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: true },
        tooltip: {
          callbacks: {
            title: function(items) {
              return `Residue Position: ${items[0].label}`;
            },
            label: function(context) {
              const val = context.raw;
              const type = val >= 0 ? 'Hydrophobic' : 'Hydrophilic';
              return ` Score: ${val} (${type})`;
            }
          }
        }
      },
      scales: {
        x: {
          title: {
            display: true,
            text: 'Residue Position (N -> C terminus)',
            font: { size: 11 }
          },
          grid: { color: '#f8fafc' }
        },
        y: {
          title: {
            display: true,
            text: 'Hydropathy Index (Kyte-Doolittle)',
            font: { size: 11 }
          },
          grid: {
            color: function(context) {
              if (context.tick.value === 0) return '#475569';
              return '#f1f5f9';
            },
            lineWidth: function(context) {
              if (context.tick.value === 0) return 1.5;
              return 1;
            }
          }
        }
      }
    }
  });
}
