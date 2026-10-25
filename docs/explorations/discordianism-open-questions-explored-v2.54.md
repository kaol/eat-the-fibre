# Discordianism: Open Questions Explored (v2.54): finitely many census orbits on 1⇉2⇇3 + 1→3 (characteristic 0)

*v2.53 found the 118 irreducible components of the fold variety of Q₀ = 1⇉2⇇3 + 1→3, each with a dense census orbit.
This round handles the boundaries. On every component, the locus where census orbits have dimension ≤ r has dimension
≤ r, for every r. By the standard orbit-dimension criterion, each component then contains finitely many census
orbits. So Folds(k̄Q₀) is a finite union of census orbits. This is lemma (L) of v2.34 in its strong form, for the
joined-ends configuration, and ChatGPT's Program 2 asked for it first. Prediction 53 held; 54 was not computed.*

---

## 1. The theorem

Let k be an algebraically closed field of characteristic 0 and Q₀ = 1⇉2⇇3 + 1→3.

**Theorem.** The folds of kQ₀ form finitely many orbits of the census symmetries. Precisely, Folds(kQ₀) is the union of
the 118 components Z of v2.53. On each Z, the connected group H_Z = G°·T_{W_Z} has finitely many orbits. Here G° is
the identity component of the substitution group, which equals Aut_Δ° (their Lie algebras agree: dim 4, `logs/aut_audit.log`),
and T_{W_Z} is the torus of class weightings compatible on all of Z.

**The criterion.** A connected group H acts on an irreducible variety Z with finitely many orbits if, for every
r < dim Z, the set Z_r = {P : dim H·P ≤ r} has dimension ≤ r.

*Proof.* Z_r is closed and H-stable. Let Y be a component of the closure of Z_r \ Z_{r−1} with dim Y = r. The orbits
of dimension r in Y are open in Y, so there is exactly one of them, and the rest of Y lies in Z_{r−1}. Components of
smaller dimension are covered by lower r. Induction on r ends at Z_0, which is a finite set by hypothesis. Above all
these sits the open orbit, of dimension dim Z (v2.53). In characteristic 0, dim H·P = rank [Lie H, P]. ∎

**Applying it** (`analysis/v2_54_orbit_strata.py`, all exact over ℚ):
- **Lie H_Z.** It is 𝔥_Z = Der_Δ ⊕ E_{W_Z}, a direct sum: substitutions within a class commute with the scaling of
  that class.
- **H_Z preserves Z.** It preserves folds compatible with W_Z and fixes vertices, so it preserves the vertex-block
  stratum. Being connected, it preserves each component.
- **The bound on Z_r.** M(P) is the matrix of the [h_i, P]. The method keeps a set F of (r+1)-minors of M, each
  reduced modulo the prime of Z. Starting with F = ∅:
  1. Decompose Z ∩ V(F) with Singular (`minAssGTZ`).
  2. For each component of dimension > r, take rational points. At a point where rank M > r, add to F the minor on M's
     pivot rows and columns at that point. It is nonzero there, so it cuts that component.
  3. Repeat until no component of dimension > r remains. Since every minor vanishes on Z_r, Z_r ⊆ Z ∩ V(F), and so
     dim Z_r ≤ r.

**Result.** 7 components are points. For the other 111, every r < dim Z succeeded: no component stayed at rank ≤ r
on all sampled points, and no Singular call timed out (`logs/v2_54/orbit_strata.txt`). ∎

## 2. Scope and caveats

- **What is proved.** Finiteness over k̄ in characteristic 0, for the connected census groups. It does not give the
  number of orbits. Stabilisers may be disconnected, so one k̄-orbit can split into several 𝔽_q-orbits, but only
  boundedly many: at most the number of components of the stabiliser, by Lang's theorem.
- **From k̄ to 𝔽_q (Conjecture 3ᶜ for Q₀).** The decomposition is defined over ℚ: finitely many components, rank
  loci, and minors. By spreading out it holds over 𝔽_p for all p outside a finite set, with the same number of orbits
  over 𝔽̄_p. That gives a uniform bound on the census count for all large p. This standard argument is sketched, not
  written out, and the finite set of exceptional primes is not determined. The census itself gives 366 at
  q = 3, 5, 7, 11.
- **Why the joining arrow matters**, in two steps:
  1. v2.52: the arrow c freezes the vertex block, which was Theorem 3's modulus.
  2. v2.53–54: once the block is frozen, the census group acts with finitely many orbits on everything that remains.

  On the unjoined 1⇉2⇇3, step 1 already fails.
- **Methods.** Exact over ℚ throughout: Singular for decompositions and dimensions, sympy for the minors and the
  ranks at rational points. The Lean files cover the realisations (v2.52), not this finiteness.

## 3. Program 2: what has been done

| ChatGPT's deliverable | status |
|---|---|
| A generic twist-fibre theorem for the joined-ends configuration, with all substitutions and stabilisers | done in a stronger form: every fold lies in one of finitely many census orbits (v2.53–54); the generic part was v2.53 |
| How the connecting arrow handles the dangerous direction 1 + s(δ_b − δ_a) | it forbids it: vertex rigidity (v2.52) |
| A reduction toward Conjecture 3ᶜ, or a minimal counterexample to the criterion | no counterexample. The method (vertex-block strata, components, orbit-dimension strata) is algorithmic and quiver-independent, so it is the candidate reduction. |

## 4. Next

1. The same pipeline on the next unobstructed quivers in the census: 1⇉2⇇3 + 1→3 + 2→4 (job 14's quiver, 1220
   orbits) and 1⇛2⇇3 + 1→3.
2. A control: 1⇉2⇉3 + 1→3. It is obstructed (consecutive classes), and its merged count 873 + 22q grows. The pipeline
   must fail there, at components with no open orbit or with dim Z_r > r, and it should fail where v2.44–v2.50 located
   the moduli.
3. Write out the spreading-out step and find the exceptional primes (compare with q = 3).

## Files

- analysis/v2_54_orbit_strata.py (uses v2_53_orbits.py for the components, W_Z and 𝔥_Z)
- logs/v2_54/: PREDICTIONS.md (with the method changes recorded before the full run), orbit_strata.txt
