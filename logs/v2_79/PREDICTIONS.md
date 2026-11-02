# v2.79 predictions (committed before the analysis below; after reading the user's v2.78 heavy-run results)

Observed (v2.78 runs on the user's machine): 1⇛2→3⇉4 (K3_after_pair) is NOT flat: 2236 (q = 2) → 2240 (q = 3), all
in the identity block (338 → 342). Prediction 140 failed, and with it Conjecture 3ᶜ″'s exception for branchings
through a parallel class, although v2.77 found no modulus in the src path family there.

143. The vertex-part analysis (v2_76_vertexpart.py, paths b·c₁ and b·c₂) on 1⇛2→3⇉4 finds families with m·m′ ≠ 0;
     among them some on which the infinitesimal test (v2_77_tangent-style, per family) finds exactly 1 modulus. The
     number of such families is 4, matching the +4.
144. Conjecture 3ᶜ‴ (path count): an acyclic quiver has finitely many census orbits iff it has no consecutive classes,
     no converging/diverging pair of classes with unjoined ends, and for every class of size ≥ 3 at most one path
     (in kQ, parallel arrows giving different paths) runs from its target to a sink and at most one from a source to
     its source.
