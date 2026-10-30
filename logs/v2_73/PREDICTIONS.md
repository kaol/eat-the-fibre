# v2.73 predictions (committed before the runs)

The three case-split blocks: Q₁ 0000,0000,0000,0001 (53 leaf files, 17 orbit closures), Q₂ 1000,010m,010m,0000
(15), 1⇉2→3→4→5 00000,00000,0010m,00000,00000 (3). Lemma L at p for a case-split block, from ℤ-certificates:
  A  the case split (v2.57/59, decision tree) is valid at p: p divides no recorded branch constant;
  B  each leaf's points (rem = 0, nz ≠ 0, in its chart) lie in its leaf prime mod p (Rabinowitsch integer per generator);
  C  each leaf prime that is not H-stable lies, mod p, in an H-stable root (the root's full-coordinate ideal, pulled
     back to the leaf's chart, reduces to 0 modulo a strong ℤ-basis of the leaf prime);
  D  the trees of the H-stable roots, in their leaf charts: C1–C4 (v2_68_uniform, leaf mode) and the node
     certificates (v2_71_gcert, leaf mode);
  E  nodes of different trees: the same full-coordinate key means the same node; different keys of equal dimension
     are separated at a point (full coordinates).

112. Q₂'s and 1⇉2→3→4→5's blocks: A holds at every p (branch constants ±1). Q₁'s block: A fails at p = 2 only
     (constants ±2, ±1/2).
113. B, C, D, E have candidate primes in {2, 3} after direct checks.
114. Hence Conjecture G on these blocks for all p ≥ 5 (and p = 3 where 113's direct checks at 3 pass); the counts are
     17, 15 and 3.
