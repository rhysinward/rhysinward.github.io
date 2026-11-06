// Maths for the in-page explorers in renewal-rt.qmd and genotype-replacement.qmd.
// Plain functions with no dependencies; random samplers are passed in from d3.
// They mirror the R code in those tutorials, so the defaults give the same estimator.

// Log gamma function (Lanczos approximation).
export function lgamma(z) {
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
  z -= 1;
  let x = c[0];
  for (let i = 1; i < 9; i++) x += c[i] / (z + i);
  const t = z + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

// Regularised lower incomplete gamma function P(a, x).
export function pgamma(x, a) {
  if (x <= 0) return 0;
  const gln = lgamma(a);
  if (x < a + 1) {
    let ap = a, sum = 1 / a, del = sum;
    for (let n = 0; n < 1000; n++) {
      ap++; del *= x / ap; sum += del;
      if (Math.abs(del) < Math.abs(sum) * 1e-15) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - gln);
  }
  let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
  for (let i = 1; i < 1000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
    c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-15) break;
  }
  return 1 - Math.exp(-x + a * Math.log(x) - gln) * h;
}

// Quantile of a gamma distribution with the given shape and rate, by bisection.
export function qgamma(p, shape, rate) {
  let lo = 0, hi = Math.max(1, shape) * 10;
  while (pgamma(hi, shape) < p) hi *= 2;
  for (let i = 0; i < 100; i++) {
    const mid = (lo + hi) / 2;
    if (pgamma(mid, shape) < p) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2 / rate;
}

// Generation interval: a gamma density with shape = mean and rate = 1 (so the
// tutorial's mean of 6 days is dgamma(1:21, shape = 6, rate = 1)), evaluated at
// whole days and normalised. The cut-off grows with the mean so long intervals
// are not truncated.
export function genInterval(mean) {
  const shape = mean, maxDay = Math.max(21, Math.ceil(mean + 4 * Math.sqrt(mean)));
  const w = [];
  for (let s = 1; s <= maxDay; s++) w.push(Math.exp((shape - 1) * Math.log(s) - s - lgamma(shape)));
  const total = w.reduce((x, y) => x + y, 0);
  return w.map(v => v / total);
}

// Total infectiousness on day t (0-based): sum over s of cases[t - s] * w[s].
export function infectiousness(cases, t, w) {
  let total = 0;
  for (let s = 1; s <= Math.min(t, w.length); s++) total += cases[t - s] * w[s - 1];
  return total;
}

// Simulate local cases from the renewal equation. Imported cases arrive each day
// and infect others, but are not themselves caused by local transmission.
// rpois is a seeded Poisson sampler factory, e.g. d3.randomPoisson.source(d3.randomLcg(seed)).
export function simulate({ days = 120, seedCases = 5, seedDays = 5, rt, w, rpois, imported = 0 }) {
  const local = new Array(days).fill(0), imp = new Array(days).fill(0), all = new Array(days).fill(0);
  for (let t = 0; t < days; t++) {
    imp[t] = imported;
    local[t] = t < seedDays ? seedCases : rpois(rt(t) * infectiousness(all, t, w))();
    all[t] = local[t] + imp[t];
  }
  return { local, imp };
}

// Cori et al. (2013): gamma posterior for Rt over a sliding window of tau days,
// with a gamma prior of shape a and scale b. If imp is given, imported cases add
// to infectiousness but are not counted as offspring. Each row records how many
// cases fell in the window, so thin windows can be left out of the plot.
export function cori({ local, imp = null, w, tau = 7, a = 1, b = 5, firstDay = 15 }) {
  const n = local.length;
  const all = imp ? local.map((v, i) => v + imp[i]) : local;
  const lambda = all.map((_, t) => (t === 0 ? 0 : infectiousness(all, t, w)));
  const out = [];
  for (let t = Math.max(firstDay - 1, tau - 1); t < n; t++) {
    let cases = 0, lam = 0;
    for (let k = t - tau + 1; k <= t; k++) { cases += local[k]; lam += lambda[k]; }
    const shape = a + cases, rate = 1 / b + lam;
    out.push({ day: t + 1, cases, mean: shape / rate, lo: qgamma(0.025, shape, rate), hi: qgamma(0.975, shape, rate) });
  }
  return out;
}

// Binomial GLM with a logit link, fitted by iteratively reweighted least squares,
// as glm(cbind(k, n - k) ~ x, family = binomial) does in R.
export function logisticFit(x, k, n) {
  let b0 = 0, b1 = 0, cov = [[NaN, NaN], [NaN, NaN]];
  for (let it = 0; it < 50; it++) {
    let s00 = 0, s01 = 0, s11 = 0, g0 = 0, g1 = 0;
    for (let i = 0; i < x.length; i++) {
      const p = 1 / (1 + Math.exp(-(b0 + b1 * x[i]))), W = n[i] * p * (1 - p), r = k[i] - n[i] * p;
      s00 += W; s01 += W * x[i]; s11 += W * x[i] * x[i]; g0 += r; g1 += r * x[i];
    }
    const det = s00 * s11 - s01 * s01;
    if (!(det > 0)) break;
    cov = [[s11 / det, -s01 / det], [-s01 / det, s00 / det]];
    const d0 = cov[0][0] * g0 + cov[0][1] * g1, d1 = cov[1][0] * g0 + cov[1][1] * g1;
    b0 += d0; b1 += d1;
    if (Math.abs(d0) + Math.abs(d1) < 1e-12) break;
  }
  return { b0, b1, se0: Math.sqrt(cov[0][0]), se1: Math.sqrt(cov[1][1]), cov };
}

// Monthly sequence counts for a genotype taking over along a logistic curve.
// bias(m) is how many times more likely a case of the new genotype is to be
// sequenced in month m than a case of the old one (1 = random sampling).
// rpois draws the monthly totals and runif one uniform per sequence, from separate
// streams, so changing the sampling scheme changes only which sequences count as
// the new genotype, not how many were sequenced or the random draws behind them.
export function simulateSequences({ months = 36, perMonth = 40, growth = 0.25, midpoint = 18, bias = () => 1, rpois, runif }) {
  const rows = [];
  const draw = rpois(perMonth);
  for (let m = 0; m < months; m++) {
    const p = 1 / (1 + Math.exp(-growth * (m - midpoint)));
    const b = bias(m), shown = b * p / (b * p + 1 - p);
    const n = draw();
    let newer = 0;
    for (let i = 0; i < n; i++) if (runif() < shown) newer++;
    rows.push({ month: m, n, newer, truth: p });
  }
  return rows;
}

// Mean of the shown estimates between two days (inclusive), or NaN if none are shown.
export function meanBetween(est, from, to) {
  const xs = est.filter(r => r.day >= from && r.day <= to && Number.isFinite(r.mean)).map(r => r.mean);
  return xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : NaN;
}

// Accessibility helpers for Observable Inputs and Plot.

// Inputs.range labels only its number box; give the slider the same name.
export function labelRange(view) {
  const text = view.querySelector("label")?.textContent?.trim();
  if (text) view.querySelectorAll("input[type=range]").forEach(r => r.setAttribute("aria-label", text));
  return view;
}

// Plot labels every group inside the chart, which is not allowed on a plain <g>.
// Expose the chart as one image (its ariaLabel and ariaDescription stay on the
// <svg>), and keep each mark's ariaDescription as data-mark for styling.
export function tidyPlot(node) {
  const svg = node.tagName.toLowerCase() === "svg" ? node : node.querySelector("svg");
  if (!svg) return node;
  svg.setAttribute("role", "img");
  svg.querySelectorAll("g[aria-label]").forEach(g => g.removeAttribute("aria-label"));
  svg.querySelectorAll("g[aria-description]").forEach(g => {
    g.dataset.mark = g.getAttribute("aria-description");
    g.removeAttribute("aria-description");
  });
  return node;
}
