// =============================================
// 1. Background particle field (About page)
// =============================================
(function () {
  var canvas = document.getElementById("bg");
  if (!canvas) return;
  var ctx = canvas.getContext("2d");
  var particles = [];
  var N = 100;
  var connectDist = 130;
  var mouse = { x: -9999, y: -9999 };

  function resize() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
  }

  function init() {
    resize();
    particles = [];
    for (var i = 0; i < N; i++) {
      particles.push({
        x: Math.random() * canvas.width,
        y: Math.random() * canvas.height,
        vx: (Math.random() - 0.5) * 0.25,
        vy: (Math.random() - 0.5) * 0.25,
        r: Math.random() * 1.1 + 0.4,
        br: Math.random() * 0.4 + 0.15
      });
    }
  }

  function tick() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (var i = 0; i < particles.length; i++) {
      for (var j = i + 1; j < particles.length; j++) {
        var dx = particles[i].x - particles[j].x;
        var dy = particles[i].y - particles[j].y;
        var d = Math.sqrt(dx * dx + dy * dy);
        if (d < connectDist) {
          ctx.beginPath();
          ctx.moveTo(particles[i].x, particles[i].y);
          ctx.lineTo(particles[j].x, particles[j].y);
          ctx.strokeStyle = "rgba(130,150,130," + ((1 - d / connectDist) * 0.07) + ")";
          ctx.lineWidth = 0.5;
          ctx.stroke();
        }
      }
    }
    for (var k = 0; k < particles.length; k++) {
      var p = particles[k];
      var dx2 = p.x - mouse.x, dy2 = p.y - mouse.y;
      var dm = Math.sqrt(dx2 * dx2 + dy2 * dy2);
      if (dm < 140 && dm > 0) {
        var f = (140 - dm) / 140 * 0.12;
        p.vx += (dx2 / dm) * f;
        p.vy += (dy2 / dm) * f;
      }
      p.vx *= 0.99; p.vy *= 0.99;
      p.x += p.vx; p.y += p.vy;
      if (p.x < -10) p.x = canvas.width + 10;
      if (p.x > canvas.width + 10) p.x = -10;
      if (p.y < -10) p.y = canvas.height + 10;
      if (p.y > canvas.height + 10) p.y = -10;
      ctx.beginPath();
      ctx.arc(p.x, p.y, p.r, 0, 6.2832);
      ctx.fillStyle = "rgba(190,200,190," + p.br + ")";
      ctx.fill();
    }
    requestAnimationFrame(tick);
  }

  window.addEventListener("resize", resize);
  document.addEventListener("mousemove", function (e) { mouse.x = e.clientX; mouse.y = e.clientY; });
  document.addEventListener("mouseleave", function () { mouse.x = -9999; mouse.y = -9999; });
  init();
  tick();
})();

// =============================================
// 2. Tab switching
// =============================================
(function () {
  var zelStarted = false;

  document.addEventListener("DOMContentLoaded", function () {
    var tabs = document.querySelectorAll(".tab");
    var pages = { about: document.getElementById("about"), research: document.getElementById("research") };

    tabs.forEach(function (tab) {
      tab.addEventListener("click", function (e) {
        e.preventDefault();
        var target = this.getAttribute("data-tab");
        tabs.forEach(function (t) { t.classList.remove("active"); });
        this.classList.add("active");
        Object.keys(pages).forEach(function (key) {
          if (key === target) pages[key].classList.remove("hidden");
          else pages[key].classList.add("hidden");
        });
        if (target === "research" && !zelStarted) {
          zelStarted = true;
          ZEL.init();
        }
      });
    });
  });
})();

// =============================================
// 3. Zel'dovich approximation simulation
// =============================================
var ZEL = (function () {

  // --- Config ---
  var NG = 150;            // particle grid dimension
  var KMAX = 12;           // max wavenumber
  var NP = NG * NG;        // total particles
  var L = 2 * Math.PI;     // box size
  var RES = 500;            // canvas pixel resolution

  // --- State ---
  var canvas, ctx, imgData, px;
  var qx, qy, psiX, psiY;
  var dens;
  var colorR, colorG, colorB; // LUT
  var growth = 0, target = 0;
  var autoPlaying = false;
  var autoT0 = 0;
  var running = false;

  // --- Helpers ---
  function grand() {
    return Math.sqrt(-2 * Math.log(Math.random())) * Math.cos(6.2832 * Math.random());
  }

  // --- Color LUT (inferno-inspired) ---
  function buildLUT() {
    colorR = new Uint8Array(256);
    colorG = new Uint8Array(256);
    colorB = new Uint8Array(256);
    var stops = [
      [0.00, 8, 6, 16],
      [0.08, 22, 12, 55],
      [0.18, 55, 18, 95],
      [0.30, 100, 28, 100],
      [0.42, 145, 40, 82],
      [0.55, 190, 70, 48],
      [0.68, 225, 120, 30],
      [0.80, 245, 180, 50],
      [0.90, 252, 225, 140],
      [1.00, 255, 252, 240]
    ];
    for (var i = 0; i < 256; i++) {
      var t = i / 255;
      var lo = 0, hi = 1;
      for (var s = 1; s < stops.length; s++) {
        if (stops[s][0] >= t) { hi = s; lo = s - 1; break; }
      }
      var range = stops[hi][0] - stops[lo][0];
      var f = range > 0 ? (t - stops[lo][0]) / range : 0;
      colorR[i] = Math.round(stops[lo][1] + f * (stops[hi][1] - stops[lo][1]));
      colorG[i] = Math.round(stops[lo][2] + f * (stops[hi][2] - stops[lo][2]));
      colorB[i] = Math.round(stops[lo][3] + f * (stops[hi][3] - stops[lo][3]));
    }
  }

  // --- Generate displacement field ---
  function generate() {
    qx = new Float64Array(NP);
    qy = new Float64Array(NP);
    psiX = new Float64Array(NP);
    psiY = new Float64Array(NP);

    // Lagrangian grid
    var idx = 0;
    for (var iy = 0; iy < NG; iy++) {
      for (var ix = 0; ix < NG; ix++) {
        qx[idx] = (ix + 0.5) / NG * L;
        qy[idx] = (iy + 0.5) / NG * L;
        idx++;
      }
    }

    // Collect Fourier modes (flat arrays for speed)
    var mkx = [], mky = [], mk2 = [], mamp = [], mph = [];
    for (var a = -KMAX; a <= KMAX; a++) {
      for (var b = -KMAX; b <= KMAX; b++) {
        if (a === 0 && b === 0) continue;
        var k2 = a * a + b * b;
        var k = Math.sqrt(k2);
        // Scale-invariant primordial spectrum P(k) ~ k^n_s
        // with Gaussian damping near Nyquist to suppress aliasing
        var ns = 0.96;
        var kNy = KMAX * 0.8;
        var Pk = Math.pow(k, ns) * Math.exp(-k * k / (2 * kNy * kNy));
        mkx.push(a); mky.push(b); mk2.push(k2);
        mamp.push(grand() * Math.sqrt(Pk));
        mph.push(Math.random() * 6.2832);
      }
    }
    var nm = mkx.length;

    // Compute Psi(q) = sum_k (amplitude / k^2) * k_hat * sin(k.q + phase)
    for (var i = 0; i < NP; i++) {
      var sx = 0, sy = 0;
      var px_ = qx[i], py_ = qy[i];
      for (var m = 0; m < nm; m++) {
        var v = mamp[m] / mk2[m] * Math.sin(mkx[m] * px_ + mky[m] * py_ + mph[m]);
        sx += mkx[m] * v;
        sy += mky[m] * v;
      }
      psiX[i] = sx;
      psiY[i] = sy;
    }

    // Normalize: RMS displacement = 1.8 grid spacings at D=1
    var ss = 0;
    for (var i = 0; i < NP; i++) ss += psiX[i] * psiX[i] + psiY[i] * psiY[i];
    var rms = Math.sqrt(ss / NP);
    var tgt = (L / NG) * 1.8;
    var sc = tgt / rms;
    for (var i = 0; i < NP; i++) { psiX[i] *= sc; psiY[i] *= sc; }
  }

  // --- Render density field ---
  function render(D) {
    dens.fill(0);
    var scale = RES / L;

    // CIC deposit
    for (var i = 0; i < NP; i++) {
      var x = ((qx[i] + D * psiX[i]) % L + L) % L;
      var y = ((qy[i] + D * psiY[i]) % L + L) % L;
      var fx = x * scale, fy = y * scale;
      var ix = Math.floor(fx), iy = Math.floor(fy);
      var dx = fx - ix, dy = fy - iy;
      var ix1 = (ix + 1) % RES, iy1 = (iy + 1) % RES;

      dens[iy  * RES + ix]  += (1 - dx) * (1 - dy);
      dens[iy  * RES + ix1] += dx       * (1 - dy);
      dens[iy1 * RES + ix]  += (1 - dx) * dy;
      dens[iy1 * RES + ix1] += dx       * dy;
    }

    // Find max (skip top 0.01% to avoid outlier domination)
    var sorted = new Float64Array(dens);
    sorted.sort();
    var maxD = sorted[Math.floor(sorted.length * 0.9997)] || 1;
    if (maxD < 0.5) maxD = 0.5;
    var logMax = Math.log(1 + maxD);

    // Write pixels
    for (var i = 0; i < RES * RES; i++) {
      var val = Math.log(1 + dens[i]) / logMax;
      if (val > 1) val = 1;
      var ci = Math.floor(val * 255);
      var pi = i * 4;
      px[pi]     = colorR[ci];
      px[pi + 1] = colorG[ci];
      px[pi + 2] = colorB[ci];
      px[pi + 3] = 255;
    }
    ctx.putImageData(imgData, 0, 0);
  }

  // --- Animation loop ---
  function frame() {
    if (!running) return;

    if (autoPlaying) {
      var elapsed = (performance.now() - autoT0) / 1000;
      // Animate D from 0 to 3.5 over 4 seconds, ease-out
      var progress = Math.min(elapsed / 4.0, 1.0);
      var eased = 1 - Math.pow(1 - progress, 2.5);
      var D = eased * 3.5;
      growth = D;
      target = D;

      var slider = document.getElementById("growthSlider");
      var label = document.getElementById("growthLabel");
      if (slider) slider.value = Math.round(D * 100);
      if (label) label.textContent = "D = " + D.toFixed(2);

      if (progress >= 1) autoPlaying = false;
    } else {
      growth += (target - growth) * 0.1;
    }

    render(growth);
    requestAnimationFrame(frame);
  }

  // --- Public init ---
  function initSim() {
    canvas = document.getElementById("zelCanvas");
    if (!canvas) return;
    canvas.width = RES;
    canvas.height = RES;
    ctx = canvas.getContext("2d");
    imgData = ctx.createImageData(RES, RES);
    px = imgData.data;
    dens = new Float64Array(RES * RES);

    buildLUT();

    // Show loading message
    ctx.fillStyle = "#08080c";
    ctx.fillRect(0, 0, RES, RES);
    ctx.fillStyle = "#444";
    ctx.font = "14px Inter, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Generating initial conditions…", RES / 2, RES / 2);

    setTimeout(function () {
      generate();
      render(0);

      // Wire controls
      var slider = document.getElementById("growthSlider");
      var label = document.getElementById("growthLabel");
      if (slider) {
        slider.addEventListener("input", function () {
          autoPlaying = false;
          target = parseFloat(this.value) / 100;
          if (label) label.textContent = "D = " + target.toFixed(2);
        });
      }

      var btn = document.getElementById("regenBtn");
      if (btn) {
        btn.addEventListener("click", function () {
          autoPlaying = false;
          ctx.fillStyle = "#08080c";
          ctx.fillRect(0, 0, RES, RES);
          ctx.fillStyle = "#444";
          ctx.font = "14px Inter, sans-serif";
          ctx.textAlign = "center";
          ctx.fillText("Generating…", RES / 2, RES / 2);
          setTimeout(function () {
            generate();
            growth = 0;
            target = 0;
            if (slider) slider.value = 0;
            if (label) label.textContent = "D = 0.00";
            render(0);
            // Auto-play again
            autoPlaying = true;
            autoT0 = performance.now();
          }, 30);
        });
      }

      // Auto-play on first load
      running = true;
      autoPlaying = true;
      autoT0 = performance.now();
      frame();
    }, 50);
  }

  return { init: initSim };
})();
