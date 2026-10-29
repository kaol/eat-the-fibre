# v2.71 predictions (committed before the runs)

Erratum found before this round: v2_70_gcert.py combined the separation primes as (∩_x cand(x)) ∪ (∩_x sep(x)); the
sound set is ∩_x (cand(x) ∪ sep(x)) (one point must certify every claim at once). The former can miss a prime.
v2.71 recomputes every node with separation per point, and extends the certificates to p = 3:
  * adaptive sampling targets 3 as well (only 2 is exempt); direct checks in characteristic q (G1, stabiliser) for
    every residual q ≥ 3 at a point where q comes only from the G1 / S basis leading coefficients;
  * the Lemma C claims (v2.68 C1–C4) checked directly over 𝔽_3 on the 34 trees whose ℤ-candidates contain 3.

105. The corrected combination changes no conclusion of v2.70: every node is certified for every p ≥ 5.
106. At p = 3 every dense node and every point node is certified (directly where needed).
107. The direct C1–C4 check at p = 3 passes on all but a few of the 34 trees; it fails on 1⇛2⇇3 + 1→3 block
     100,001,001 (as in v2.68), which then needs the char-3 tree.
108. Hence Conjecture G holds at every q = 3^k for every quiver whose trees all pass 107 (expected: all but
     1⇛2⇇3 + 1→3).
