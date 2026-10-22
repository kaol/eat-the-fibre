# v2.40 predictions (committed before writing the proof)

Goal: remove the generation hypothesis from v2.39 on 1⇉2⇉3.

1. Over every reduced commutative ring k, every Φ that is multiplicative (Φ(vw) = Φv·Φw) with a two-sided matrix
   inverse Ψ is `Filt`. Planned proof: b² = 0 for rad basis elements gives vertex coordinates squaring to 0; products
   give rad → rad²; the sink annihilates rad on the left and the source on the right; (1 − Φ₀₀)·Φ(e₀) and
   (1 − Φ₂₂)·Φ(e₂) are idempotents in rad, hence 0; row sums from Φ(1) = 1 finish. No domain or connectedness
   hypothesis on k is needed beyond reducedness.
2. Corollary, on the standard axioms: over a reduced nontrivial commutative ring, an algebra automorphism intertwining
   members of families i and j forces i = j. So the 19 families are 19 behaviours with no appeal to v2.33's
   generation statement.

## Outcomes

1. **Held** (`Mult.filt`, `lean/EatTheFibre/AutFilt.lean`), by the planned proof. The first build had mechanical
   failures, none mathematical:
   - `rw` on `fin_cases`-shaped terms;
   - `col` rewriting the wrong occurrence;
   - one sign in a `linear_combination`;
   - a `simp_all` timeout.
2. **Held** (`fam_behaviours_distinct_aut`), on the standard axioms.
