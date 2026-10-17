# Discordianism: Open Questions Explored (v2.27) — draft

*Predictions recorded before the censuses run (this file is completed when they finish).*

- **1⇉2⇇3** (two parallel classes, not consecutive): Conjecture 3ᶜ predicts finitely many classes. Prediction: the
  orbit count under the local class torus (`orbits_GxTloc`) is the same over 𝔽₃ and 𝔽₅. (G-orbit counts may grow: the
  quiver has a same-orientation fork at vertex 2, which v2.22 showed gives grading/torus-absorbed families.)
- **1⇉2⇉3 + 1→3** (consecutive classes, but Theorem 2's hypothesis fails): no prediction.

**Outcome (recorded before any analysis of it): the 1⇉2⇇3 prediction failed.** Local-torus orbit counts over
𝔽₃, 𝔽₅, 𝔽₇, 𝔽₁₁: 453, 461, 469, 485 = 445 + 4(q − 1). G-orbits 485, 557, 629, 773. G × grading = G × T = local torus
at every prime. The only orbit type whose count grows is the one of size q²(q − 1)² (63, 71, 79, 95). Raw lines:
`logs/results_container.jsonl`.

**Second set of predictions (recorded after Theorem 3 below was found and checked, before these censuses run):**

- **1⇇2⇉3** (opposite quiver): identical numbers to 1⇉2⇇3 over 𝔽₃ and 𝔽₅ (N, and the orbit count under every group).
  This follows from op-duality, so it tests the engine, not the mathematics.
- **1⇉2⇇3, 2→4**: the four families of Theorem 3 persist (symbolic check), and no other mechanism is present. Prediction:
  the local-torus orbit count grows by exactly 8 from 𝔽₃ to 𝔽₅, and every vertex block that carries exactly one orbit is a
  generic member of one of the four families.

**1⇇2⇉3 outcome:** N = 6376 / 59 972 and 485 / 557 G-orbits, 453 / 461 under every larger group, over 𝔽₃ / 𝔽₅: identical
to 1⇉2⇇3, as op-duality requires.

**Third prediction (recorded before the census):** Conjecture 3ᵛ (below) says the local-torus orbit count is finite iff no
two parallel classes of size ≥ 2 share a vertex. **1⇉2⇇3 + 1→3** violates the hypothesis of Theorem 3 (the families
P_s are not folds there), but its classes share vertex 2. Prediction under 3ᵛ: its local-torus count grows from 𝔽₃ to
𝔽₅. A constant count refutes 3ᵛ as stated.
