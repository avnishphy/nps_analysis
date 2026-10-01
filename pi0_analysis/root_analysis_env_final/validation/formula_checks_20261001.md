# Synthetic formula checks — 2026-10-01

Status: **synthetic-check-only**. No experimental or SIMC production file was
read. These calculations test algebra and expose estimator assumptions; they do
not measure a correction or certify publication readiness.

## Source definitions checked first

```bash
nl -ba src/analysis/nps_time_bg.h | sed -n '1,270p'
nl -ba src/analysis/nps_analysis_main.C | sed -n '2585,2720p'
nl -ba src/analysis/nps_analysis_main.C | sed -n '2960,3080p'
sed -n '140,178p' /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/docs/pi0_uncertainty_audit_2026-10-01.md
```

The first three commands confirm the region widths, the histogram range, and
that pair selection precedes timing-region filling. The fourth recovers the
audit's precise toy definitions for an independent rerun.

## Corrected executable check

Run from the `FINAL` workspace:

```bash
python3 -s - <<'PY_CHECK'
import numpy as np

rng = np.random.default_rng(20261001)
on, off, alpha = 100.0, 250.0, 0.2
signal = on - alpha * off
replicas = rng.poisson(on, 200000) - alpha * rng.poisson(off, 200000)
print("frozen-purity variance:", signal**2 / (on + off))
print("subtraction variance:", on + alpha**2 * off)
print("toy variance:", np.var(replicas, ddof=1))

d = np.array([10.0, 20.0])
m = np.array([12.0, 18.0])
v = np.array([5.0, 8.0])
for c in [1.0, 0.001]:
    y, mu, variance = c*d, c*m, c*c*v
    terms = 2*(mu-y+y*np.log(y/mu))
    print(c, "ordinary:", terms.sum(),
          "scaled:", (terms/(variance/y)).sum())

print("clone/event variance:", 80*(2/80)**2, 2**2)

windows = [(141.,143.), (143.,145.), (145.,147.),
           (153.,155.), (155.,157.), (157.,159.)]
print("areas:", 2*2,
      sum((hi-lo)**2 for lo,hi in windows),
      sum(2*(hi-lo) for lo,hi in windows), 6*6)

def accepted_square(base_difference, pair_cut, width=6.0):
    side = max(0.0, min(width, width + pair_cut - base_difference))
    return 0.5*side*side

print("full-box accepted areas:",
      accepted_square(12.0, 10.0),
      accepted_square(16.0, 13.0),
      36.0)
PY_CHECK
```

Observed values:

```text
frozen-purity variance: 7.142857142857143
subtraction variance: 110.0
toy variance: 110.96253747547135
1.0 ordinary: 0.5679894904339626 scaled: 1.2431892940244524
0.001 ordinary: 0.0005679894904339583 scaled: 1.2431892940244418
clone/event variance: 0.05000000000000001 4
areas: 4 24.0 24.0 36
full-box accepted areas: 8.0 4.5 36.0
```

Interpretation:

- The on/off example has `S=100-0.2*250=50`. Treating the fitted purity
  `S/(on+off)` as fixed gives the displayed 7.14 variance proxy, whereas direct
  independent-Poisson subtraction gives `100+0.2^2*250=110`. The fixed-seed
  ensemble is consistent with the latter. These counts are illustrative only.
- Multiplying the same weighted histogram by `0.001` multiplies ordinary
  Poisson deviance by `0.001`; it leaves the shown scaled-Poisson expression
  invariant. This demonstrates unit behavior, not validity for fitted weights.
- Eighty repeated smears with weight `2/80` give copy-wise squared-weight sum
  `0.05`, while one original generator event carrying total weight 2 has square
  4. Repeated smears are correlated integration samples, not 80 independent
  generated events.
- The nominal prompt, diagonal, horizontal/vertical, and each full-box areas are
  4, 24, 24, and 36 ns2. Intersecting a 6x6 ns2 full box with the pair-time
  selection leaves a triangle of area 8 ns2 in the unshifted >=3-cluster mode
  and 4.5 ns2 in the shifted >=3-cluster mode. The exactly-two-cluster bypass
  retains all 36 ns2.

## Failed and superseded checks

The first attempted Python one-liner failed before execution with:

```text
SyntaxError: EOL while scanning string literal
```

Cause: an accidental newline inside the `ordinary_deviance_scale0p001` f-string.
No file or state changed.

After fixing the syntax, a hand-written preliminary script printed timing areas
of 48 and 72 ns2 and zero accepted gated area. Review against
`nps_time_bg.h:168-218` showed that it had incorrectly summed distinct
diagonal/horizontal/full components and used the wrong triangle side. Those
numbers were discarded. The source-matched command above reports each aggregate
separately and uses `side = width + pair_cut - base_difference`, yielding the
verified 24/36 and 8/4.5 ns2 values.

The audit check and this rerun use the same fixed NumPy generator and seed, so
the toy variance matches exactly. This is a reproducibility check, not an
independent random ensemble.
