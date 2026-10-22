# Discordianism: Open Questions Explored (v2.44): families that change behaviour

*Since v2.33 every family in this project has been a single behaviour, an orbit of a one-parameter group of algebra
automorphisms whose parameter moves the splitting relative to Δ. This round finds the four families of
1⇉2⇉3 + 1→3 that v2.42 could not locate, and they are different in kind. Along them the splitting (im P, ker P)
itself deforms: the image stays fixed while the kernel tilts, and to first order the deformation is not an
automorphism. In characteristic 0 the members realise infinitely many behaviours. One of them is checked exactly over
ℚ and as a fold in Lean. Six predictions were committed in advance, plus one written before but committed after its
run (logs/v2_42, logs/v2_44). Every outcome held as stated, though for 11 and 12 the counting of families is not
settled.*

---

## 0. Summary

| # | Statement | Status |
|---|---|---|
| 1 | With the shear merged in, 1⇉2⇉3 + 1→3 has 939, 983, 1027 = 873 + 22q orbits at q = 3, 5, 7. The 18 inner families meet q orbits at each prime. The two image-dimension-7 groups grow as 17 + 2q each | computed (v2.43; predictions 6, 7) |
| 2 | v2.34's invariants separate the two kinds of family. Inner-family points have b = 0, e ≥ 1. In the dim-7 groups only the classes (b, e) = (1, 0) and (2, 0) grow, by q each | computed (predictions 8, 9 held) |
| 3 | At a (1, 0) point, T_P Folds has dimension 8 and the derivation directions [Der, P] ∩ T have dimension 7. There is exactly one transversal tangent direction, and no derivation produces it | computed (mod a 7-digit prime) |
| 4 | P(t) = P₀ + t·u·wᵀ is a fold over ℚ[t] for every t, with u = e₃ + b₁ − b₂ ∈ im P₀ and w = a₁* + a₂* − c* + (a₁b₁)* − (a₂b₂)*, which kills im P₀. Its points meet exactly q orbits at q = 3, 5, 7 | **exact over ℚ** (prediction 10a); orbit counts computed |
| 5 | N = u·wᵀ ∉ [Der(kQ), P₀] over ℚ (Der has dimension 22; rank 22 → 23). In characteristic 0 each Aut(kQ)-orbit meets the line in finitely many points | **exact over ℚ** (prediction 10b) + the orbit-tangent argument below |
| 6 | `line_isFold`: the line is a fold over every commutative ring, for every t | **proved (Lean)** (prediction 13) |
| 7 | Kernel-moving lines through the growing orbits are abundant: 32 distinct lines from 16 representatives, in both groups | computed (predictions 11, 12 held as stated) |

## 1. Why the behaviour changes

Over a field of characteristic 0, Aut(kQ) is an algebraic group with Lie algebra Der(kQ). The orbit Aut(kQ)·P₀ is
smooth and locally closed, with tangent space [Der(kQ), P₀] at P₀.

Suppose infinitely many points of the line P(t) were in that orbit. Then the line would lie in the orbit's closure,
and an open part of it would lie in the orbit. So its tangent N would be in [Der, P₀]. The exact rank computation says
it is not. Hence the line meets each behaviour in finitely many points, and over ℚ it carries infinitely many
behaviours.

This is the first such family in the project:
- on 1⇉2⇉3 the 19 families are 19 behaviours, one each (v2.36–v2.40);
- on 1⇉2⇇3 the four families are three behaviours (v2.41);
- here, a family is a curve of behaviours.

## 2. The shape of the deformation

Every line found moves the kernel and fixes the image. P(t) = P₀ + t·u·wᵀ keeps im P(t) = im P₀, a 7-dimensional
subalgebra, and tilts the 5-dimensional kernel along w. The algebra im P is the same for all t. What varies is the
complement ker P(t), a subalgebra and coideal, relative to it.

The vector u = e₃ + b₁ − b₂ involves the sink vertex. The functional w involves the joining arrow c and the paths
a₁b₁, a₂b₂ parallel to it. This is consistent with the joining arrow being the source of the phenomenon: 1⇉2⇉3 has
no such families.

## 3. What is not settled

- **Counting.** The growth says 4 families: 2 per group, one each in classes (1, 0) and (2, 0). The lines are not in
  one-to-one correspondence with them. 32 different lines overlap in orbits, and two lines through the same point
  can share 6 of 7 orbits (`lineclass_5107_q7.txt`).
  - In the L group at q = 3, 5, 7, the lines cover q + 1 of the (1, 0) orbits and q + 3 of the (2, 0) orbits.
  - The R group search used 6 representatives, and its coverage (6, 7, 9 and 8, 8, 9) is incomplete.
  - The natural next object is the moduli curve of behaviours, not individual lines.
- **The behaviour statement** is proved for one line in characteristic 0, from the exact rank. Over 𝔽_q it is
  plausible but not shown.
- **Lean** has the fold property, not yet the non-derivation property. That needs Der(kQ) for a quiver with paths
  parallel to an arrow.

## 4. Consequences for the finiteness question

The finiteness criterion (Conjecture 3ᶜ, the open half) is about orbit counts, so it is not affected. The behaviour
picture, though, is now two-layered:
- Moduli of position: inner families, one behaviour each. These come from Inn(1 + t·x) not being Δ-compatible.
- Moduli of behaviour: kernel-tilting lines. These come from complementary subalgebra-coideals to a fixed image
  forming a positive-dimensional family modulo automorphisms.

1⇉2⇉3 + 1→3 has both kinds.

## 5. Next

1. The moduli of behaviours of the dim-7 groups: classify the kernels complementary to the fixed image, modulo the
   stabiliser of the image in Aut(kQ), and match their dimension to the 4 families.
2. Lean: closed forms for the 18 inner families on 1⇉2⇉3 + 1→3, and the non-derivation certificate for the line.
3. Job 14 on the GPU machine.

Files:
- `analysis/v2_44_{be_orbits,lines,exact,line_driver,lineclass,coverage}.py`;
- `lean/gen/gen_line.py`, `lean/EatTheFibre/Gen/two_pairs_short_line.lean`;
- `logs/v2_44/`.
