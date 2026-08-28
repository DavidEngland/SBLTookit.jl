function wrapLabel(label) {
  if (typeof label !== "string") return label;
  if (label.length <= 16) return label;
  const words = label.split(" ");
  const lines = [];
  let currentLine = "";
  words.forEach((word) => {
    if ((currentLine + word).length > 16) {
      if (currentLine.trim()) lines.push(currentLine.trim());
      currentLine = word + " ";
    } else {
      currentLine += word + " ";
    }
  });
  if (currentLine.trim()) lines.push(currentLine.trim());
  return lines;
}

const sharedTooltipConfig = {
  plugins: {
    tooltip: {
      callbacks: {
        title: function (tooltipItems) {
          const item = tooltipItems[0];
          let label = item.chart.data.labels
            ? item.chart.data.labels[item.dataIndex]
            : undefined;
          if (Array.isArray(label)) {
            return label.join(" ");
          } else if (label) {
            return label;
          } else if (item.parsed && item.parsed.x !== undefined) {
            return "Shear S: " + item.parsed.x.toFixed(4) + " s⁻¹";
          }
          return "";
        },
      },
    },
  },
};

function simulateGSPT(G0, epsilon, l0) {
  const beta = 5.0;
  const N2 = 0.1;
  const B0_max = 0.05;
  const delta_reg = 0.01;
  const gamma_s = 1.8;
  const r_s = 0.15;
  const dt = 0.002;
  const t_max = 80.0;
  const n_steps = Math.round(t_max / dt);

  const target = (27.0 * B0_max * B0_max) / (4.0 * Math.pow(l0, 4));
  const a = beta * N2;
  const b = 0.0;
  const c = -target;

  const p_c = b - (a * a) / 3.0;
  const q_c = (2.0 * Math.pow(a, 3)) / 27.0 - (a * b) / 3.0 + c;
  const delta_c = Math.pow(q_c / 2.0, 2) + Math.pow(p_c / 3.0, 3);

  let S_fold = 0.39989;
  if (delta_c <= 0) {
    const arg = ((3.0 * q_c) / (2.0 * p_c)) * Math.sqrt(-3.0 / p_c);
    const clampedArg = Math.max(-1.0, Math.min(1.0, arg));
    const phi = Math.acos(clampedArg);
    const R = 2.0 * Math.sqrt(-p_c / 3.0);
    const x1 = R * Math.cos(phi / 3.0) - a / 3.0;
    const x2 = R * Math.cos((phi + 2.0 * Math.PI) / 3.0) - a / 3.0;
    const x3 = R * Math.cos((phi + 4.0 * Math.PI) / 3.0) - a / 3.0;
    const validX = [x1, x2, x3].filter((x) => x > 0);
    if (validX.length > 0) {
      const maxX = Math.max(...validX);
      S_fold = Math.sqrt(maxX);
    }
  } else {
    const term1 = -q_c / 2.0 + Math.sqrt(delta_c);
    const term2 = -q_c / 2.0 - Math.sqrt(delta_c);
    const y = Math.cbrt(term1) + Math.cbrt(term2);
    const x = y - a / 3.0;
    if (x > 0) S_fold = Math.sqrt(x);
  }
  const e_fold = (3.0 * B0_max) / (2.0 * S_fold * S_fold);

  let e = 0.5;
  let S = 1.0;
  const times = [];
  const e_arr = [];
  const S_arr = [];
  const trajectory = [];

  const sampleInterval = Math.round(n_steps / 80);

  for (let i = 0; i < n_steps; i++) {
    const t = i * dt;
    const B0 = B0_max * ((e * e) / (e * e + delta_reg * delta_reg));
    const denom = 1.0 + (beta * N2) / Math.max(S * S, 1e-8);
    const de_dt =
      (1.0 / epsilon) * (l0 * e * S * S - B0 - Math.pow(e, 3) / (l0 * denom));
    const dS_dt = G0 - gamma_s * e * S - r_s * S;

    e = Math.max(e + de_dt * dt, 1e-4);
    S = Math.max(S + dS_dt * dt, 0.0);

    if (i % sampleInterval === 0) {
      times.push(t.toFixed(1) + "s");
      e_arr.push(e);
      S_arr.push(S);
      trajectory.push({ x: S, y: e });
    }
  }

  return { S_fold, e_fold, times, e_arr, S_arr, trajectory };
}

let bifurcationChartInstance = null;
let phasePortraitChartInstance = null;

document.addEventListener("DOMContentLoaded", () => {
  const campaignLabelsRaw = [
    "CASES-99 (Kansas Prairie)",
    "GABLS3 (Cabauw Mast)",
    "SHEBA (Arctic Pack)",
  ];
  const campaignLabelsWrapped = campaignLabelsRaw.map(wrapLabel);

  new Chart(document.getElementById("campaignChart"), {
    type: "bar",
    data: {
      labels: campaignLabelsWrapped,
      datasets: [
        {
          label: "Vertical Levels (Nz)",
          data: [7, 20, 2],
          backgroundColor: "#14b8a6",
          borderRadius: 6,
        },
        {
          label: "Max Height z (m)",
          data: [55, 200, 10],
          backgroundColor: "#38bdf8",
          borderRadius: 6,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: {
          ticks: { color: "#94a3b8", font: { size: 11 } },
          grid: { display: false },
        },
        y: { ticks: { color: "#94a3b8" }, grid: { color: "#334155" } },
      },
      plugins: {
        legend: { labels: { color: "#cbd5e1", font: { size: 12 } } },
        ...sharedTooltipConfig.plugins,
      },
    },
  });

  const ceofVarianceLabelsRaw = [
    "Mode 1: Standing LLJ",
    "Mode 2: Propagating IGW",
  ];
  const ceofVarianceLabelsWrapped = ceofVarianceLabelsRaw.map(wrapLabel);

  new Chart(document.getElementById("ceofVarianceChart"), {
    type: "doughnut",
    data: {
      labels: ceofVarianceLabelsWrapped,
      datasets: [
        {
          data: [91.37, 8.62],
          backgroundColor: ["#38bdf8", "#6366f1"],
          borderWidth: 0,
          hoverOffset: 6,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          position: "bottom",
          labels: { color: "#cbd5e1", font: { size: 12 } },
        },
        ...sharedTooltipConfig.plugins,
      },
    },
  });

  const heightLabels = [
    "1.5m",
    "5.0m",
    "10.0m",
    "20.0m",
    "30.0m",
    "45.0m",
    "55.0m",
  ];

  new Chart(document.getElementById("ceofProfileChart"), {
    type: "line",
    data: {
      labels: heightLabels,
      datasets: [
        {
          label: "Mode 1 (Standing Jet)",
          data: [0.0908, 0.0439, 0.0062, -0.002, -0.0013, 0.0099, 0.0237],
          borderColor: "#38bdf8",
          backgroundColor: "rgba(56, 189, 248, 0.1)",
          borderWidth: 3,
          tension: 0.3,
          fill: true,
        },
        {
          label: "Mode 2 (Propagating Wave)",
          data: [-0.1476, -0.1436, -0.1365, -0.1713, -0.1799, -0.137, -0.1493],
          borderColor: "#818cf8",
          backgroundColor: "rgba(129, 140, 248, 0.1)",
          borderWidth: 3,
          tension: 0.3,
          fill: true,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: { ticks: { color: "#94a3b8" }, grid: { color: "#334155" } },
        y: {
          title: {
            display: true,
            text: "Phase Gradient (rad/m)",
            color: "#94a3b8",
          },
          ticks: { color: "#94a3b8" },
          grid: { color: "#334155" },
        },
      },
      plugins: {
        legend: { labels: { color: "#cbd5e1" } },
        ...sharedTooltipConfig.plugins,
      },
    },
  });

  const regHeights = ["5m", "15m", "25m", "30m (Nose)", "35m", "45m", "55m"];

  new Chart(document.getElementById("regularizationChart"), {
    type: "line",
    data: {
      labels: regHeights,
      datasets: [
        {
          label: "Unregularized Ri (Singular Spike)",
          data: [0.15, 0.22, 0.85, 12.5, 0.9, 0.35, 0.28],
          borderColor: "#f43f5e",
          borderDash: [5, 5],
          borderWidth: 2,
          pointRadius: 4,
          pointBackgroundColor: "#f43f5e",
        },
        {
          label: "Z0HR Regularized Ri (Smooth)",
          data: [0.15, 0.22, 0.48, 0.62, 0.51, 0.35, 0.28],
          borderColor: "#14b8a6",
          borderWidth: 3,
          pointRadius: 5,
          pointBackgroundColor: "#14b8a6",
          tension: 0.3,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: { ticks: { color: "#94a3b8" }, grid: { color: "#334155" } },
        y: {
          title: {
            display: true,
            text: "Richardson Number (Ri)",
            color: "#94a3b8",
          },
          ticks: { color: "#94a3b8" },
          grid: { color: "#334155" },
        },
      },
      plugins: {
        legend: { labels: { color: "#cbd5e1" } },
        ...sharedTooltipConfig.plugins,
      },
    },
  });

  function updateSimulationUI() {
    const G0 = parseFloat(document.getElementById("slider-g0").value);
    const epsilon = parseFloat(document.getElementById("slider-eps").value);
    const l0 = parseFloat(document.getElementById("slider-l0").value);

    document.getElementById("val-g0").textContent = G0.toFixed(2);
    document.getElementById("val-eps").textContent = epsilon.toFixed(2);
    document.getElementById("val-l0").textContent = l0.toFixed(2);

    const sim = simulateGSPT(G0, epsilon, l0);

    document.getElementById("kpi-s-fold").textContent = sim.S_fold.toFixed(4);
    document.getElementById("kpi-e-fold").textContent = sim.e_fold.toFixed(4);
    document.getElementById("ctrl-s-fold").textContent =
      sim.S_fold.toFixed(5) + " s⁻¹";
    document.getElementById("ctrl-e-fold").textContent =
      sim.e_fold.toFixed(5) + " m²s⁻²";
    document.getElementById("info-s-fold").textContent =
      sim.S_fold.toFixed(5) + " s⁻¹";
    document.getElementById("info-e-fold").textContent =
      sim.e_fold.toFixed(5) + " m²s⁻²";
    document.getElementById("fold-summary-text").innerHTML =
      "S<sub>fold</sub> = " + sim.S_fold.toFixed(5) + " s<sup>-1</sup>";

    if (bifurcationChartInstance) {
      bifurcationChartInstance.data.labels = sim.times;
      bifurcationChartInstance.data.datasets[0].data = sim.e_arr;
      bifurcationChartInstance.data.datasets[1].data = sim.S_arr;
      bifurcationChartInstance.update();
    } else {
      bifurcationChartInstance = new Chart(
        document.getElementById("bifurcationChart"),
        {
          type: "line",
          data: {
            labels: sim.times,
            datasets: [
              {
                label: "Turbulent Kinetic Energy e (m²/s²)",
                data: sim.e_arr,
                borderColor: "#f59e0b",
                backgroundColor: "rgba(245, 158, 11, 0.15)",
                borderWidth: 2,
                tension: 0.1,
                fill: true,
                yAxisID: "y",
              },
              {
                label: "Shear S (s⁻¹)",
                data: sim.S_arr,
                borderColor: "#38bdf8",
                borderDash: [4, 4],
                borderWidth: 2,
                pointRadius: 0,
                yAxisID: "y1",
              },
            ],
          },
          options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
              x: {
                ticks: { color: "#94a3b8", maxTicksLimit: 12 },
                grid: { color: "#334155" },
              },
              y: {
                type: "linear",
                display: true,
                position: "left",
                title: {
                  display: true,
                  text: "TKE e (m²/s²)",
                  color: "#f59e0b",
                },
                ticks: { color: "#f59e0b" },
                grid: { color: "#334155" },
              },
              y1: {
                type: "linear",
                display: true,
                position: "right",
                title: {
                  display: true,
                  text: "Shear S (s⁻¹)",
                  color: "#38bdf8",
                },
                ticks: { color: "#38bdf8" },
                grid: { drawOnChartArea: false },
              },
            },
            plugins: {
              legend: { labels: { color: "#cbd5e1" } },
              ...sharedTooltipConfig.plugins,
            },
          },
        },
      );
    }

    if (phasePortraitChartInstance) {
      phasePortraitChartInstance.data.datasets[0].data = sim.trajectory;
      phasePortraitChartInstance.update();
    } else {
      phasePortraitChartInstance = new Chart(
        document.getElementById("phasePortraitChart"),
        {
          type: "scatter",
          data: {
            datasets: [
              {
                label: "GSPT Trajectory Orbit (S vs e)",
                data: sim.trajectory,
                showLine: true,
                borderColor: "#14b8a6",
                backgroundColor: "rgba(20, 184, 166, 0.2)",
                borderWidth: 2,
                pointRadius: 2,
                tension: 0.2,
              },
            ],
          },
          options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
              x: {
                type: "linear",
                position: "bottom",
                title: {
                  display: true,
                  text: "Shear S (s⁻¹)",
                  color: "#38bdf8",
                },
                ticks: { color: "#94a3b8" },
                grid: { color: "#334155" },
              },
              y: {
                title: {
                  display: true,
                  text: "TKE e (m²/s²)",
                  color: "#f59e0b",
                },
                ticks: { color: "#f59e0b" },
                grid: { color: "#334155" },
              },
            },
            plugins: {
              legend: { labels: { color: "#cbd5e1" } },
              ...sharedTooltipConfig.plugins,
            },
          },
        },
      );
    }
  }

  document
    .getElementById("slider-g0")
    .addEventListener("input", updateSimulationUI);
  document
    .getElementById("slider-eps")
    .addEventListener("input", updateSimulationUI);
  document
    .getElementById("slider-l0")
    .addEventListener("input", updateSimulationUI);

  updateSimulationUI();
});
