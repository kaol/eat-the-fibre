# Discordianism: Open Questions Explored (v2.28) — draft

*Predictions recorded before the censuses and analyses run (this file is completed when they finish).*

The post hoc reading of v2.27 §5: Folds(kQ) has infinitely many local-torus orbits iff Q has consecutive parallel
classes (size ≥ 2), or converging/diverging ones whose outer ends are not joined by an arrow.

- **P1.** 1⇉2⇇3 + 1→3 over 𝔽₇: 366 local-torus orbits (as over 𝔽₃ and 𝔽₅).
- **P2.** 1⇛2⇇3 + 1→3 (a larger converging class, outer ends adjacent): the local-torus count is the same over 𝔽₃ and 𝔽₅.
- **P3.** 1⇉2⇉3 + 1→3: its moduli are in the radical, as for 1⇉2⇉3, so the number of vertex blocks is the same over 𝔽₃
  and 𝔽₅.

**Outcomes so far (appended in order):**

- **P1 holds.** 1⇉2⇇3 + 1→3 over 𝔽₇: N = 502 306, 651 G-orbits, 392 under G × T, **366** under the local torus.
- **P2: not yet decided.** Over 𝔽₃: 752 local-torus orbits (N = 316 612). The 𝔽₅ run exceeds the container: first the
  symmetry plan (H has ~10⁶ elements; the engine is now rewritten, see `tests/test_symmetry_plan.py`), then memory
  (OOM at 6.1 GB with ~5.7·10⁵ fibres). Handed over.
- **P3 holds.** 1⇉2⇉3 + 1→3 has 21 vertex blocks over 𝔽₃ and over 𝔽₅ (`analysis/loc_blocks.py`). Its local-torus
  growth (1668 → 1840) sits inside 13 of them, in orbit types of size q²(q − 1)^k, k = 2, 3, 4.

**P4 (recorded before the census).** 1⇉2⇇3 + 1→3 + 2→4 (converging classes, adjacent outer ends, a tail at v): the
local-torus count is the same over 𝔽₃ and 𝔽₅.
