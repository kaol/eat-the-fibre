"""Generate the behaviour-distinctness modules (v2.39):
  lean/EatTheFibre/Gen/two_pairs_behaviour.lean  — the 19 families of 1⇉2⇉3 are 19 distinct behaviours;
  lean/EatTheFibre/Gen/pairs_in_behaviour.lean   — Theorem 3's F1, F2, F3, F4 on 1⇉2⇇3 are 4 distinct behaviours
                                                  modulo Inn ⋊ Sub (for t, t' ≠ 0).

Filt Φ: the rows of vertices are unit rows (ε∘Φ = ε, ε : kQ → kQ/rad) and Φ l c = 0 when len l < len c (Φ preserves
rad^n). Filt is closed under products. Every class substitution and every inner automorphism Inn(1 + x), x ∈ rad
(an explicit matrix in the coordinates of x, proved to be y ↦ u y u⁻¹), is Filt. On 1⇉2⇉3 these generate Aut(kQ)
(v2.33: Aut(kQ) = Inn ⋊ (Sub ⋊ Aut(Q)), and Aut(Q) only permutes arrows within classes there). On 1⇉2⇇3 they generate
Inn ⋊ Sub; Aut(Q) adds the source swap, which v2.38 found to exchange F1 and F3.

If Filt Φ, Ψ are mutually inverse and Φ·P = P'·Φ, then for a vertex v: row v of P vanishes on every column of length
≥ n whenever row v of P' does. For each ordered pair, a vertex where the depths differ, and an entry ±1 (or ±t,
using t ≠ 0) at full depth in the deeper row, refute the intertwining with Φ or with Ψ.
Usage: python3 lean/gen/gen_behaviour.py CERT.pkl"""
import os, sys, re, ast, itertools
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
cert_file = sys.argv[1]
sys.argv = ['x', cert_file, '/nonexistent']
exec(open(os.path.join(HERE, '..', '..', 'analysis', 'v2_38_disjoint.py')).read().split("os.makedirs(outdir")[0])
OUT = os.path.join(HERE, '..', 'EatTheFibre', 'Gen')
sys.argv = ['x']
GEN = {'__file__': os.path.join(HERE, 'gen.py'), '__name__': 'gen'}
exec(open(os.path.join(HERE, 'gen.py')).read().split("\nif __name__")[0], GEN)
tsym = sp.Symbol('t')

def lean(e):
    s = str(sp.expand(e)).replace('**', '^')
    return re.sub(r'\b([vwx])(\d+)\b', r'\1 \2', s)

def generic(qname):
    """the quiver-independent part, instantiated for the path algebra of `qname`"""
    basis, M, DL, nms = GEN['names_of'](qname); d = len(basis)
    ln = [0 if b[0] == 'v' else len(b[1]) for b in basis]
    nv = ln.count(0)
    V = sp.symbols(f'v0:{d}'); W = sp.symbols(f'w0:{d}'); X = sp.symbols(f'x0:{d}')
    def mul(v, w):
        return [sp.expand(sum(int(M[i][j][m]) * v[i] * w[j] for i in range(d) for j in range(d) if int(M[i][j][m])))
                for m in range(d)]
    one = [1] * nv + [0] * (d - nv)
    x = [0] * nv + list(X[nv:])
    u = [one[i] + x[i] for i in range(d)]
    pw, ui, sign = one[:], one[:], 1
    while True:                                  # u⁻¹ = Σ (−x)^n, finite since rad is nilpotent
        pw = mul(pw, x); sign = -sign
        if all(e == 0 for e in pw): break
        ui = [sp.expand(a + sign * b) for a, b in zip(ui, pw)]
    assert mul(u, ui) == one and mul(ui, u) == one
    E = lambda j: [1 if i == j else 0 for i in range(d)]
    Inn = sp.Matrix([mul(mul(u, E(j)), ui) for j in range(d)]).T
    for i in range(d):
        for j in range(d):
            if ln[i] == 0: assert Inn[i, j] == (1 if i == j else 0)
            if ln[i] < ln[j]: assert Inn[i, j] == 0
    mulv = mul(V, W)
    vec = lambda xs: "![" + ", ".join(str(a) for a in xs) + "]"
    txt = f'''variable {{k : Type*}} [CommRing k]

/-- path length of a basis element -/
def len : Fin {d} → ℕ := {vec(ln)}

/-- `ε ∘ Φ = ε` (vertex rows are unit rows) and `Φ (rad^n) ⊆ rad^n` -/
structure Filt (Φ : Matrix (Fin {d}) (Fin {d}) k) : Prop where
  vert : ∀ i j, len i = 0 → Φ i j = if i = j then 1 else 0
  tri : ∀ i j, len i < len j → Φ i j = 0

lemma Filt.one : Filt (1 : Matrix (Fin {d}) (Fin {d}) k) :=
  ⟨fun i j _ => Matrix.one_apply, fun i j h => Matrix.one_apply_ne (fun e => by rw [e] at h; exact lt_irrefl _ h)⟩

lemma Filt.mul {{A B : Matrix (Fin {d}) (Fin {d}) k}} (hA : Filt A) (hB : Filt B) : Filt (A * B) := by
  constructor
  · intro i j hi
    rw [Matrix.mul_apply, Finset.sum_eq_single i]
    · rw [hA.vert i i hi, if_pos rfl, one_mul, hB.vert i j hi]
    · intro l _ hl; rw [hA.vert i l hi, if_neg (fun e => hl e.symm), zero_mul]
    · intro h; exact absurd (Finset.mem_univ i) h
  · intro i j hij
    rw [Matrix.mul_apply]
    apply Finset.sum_eq_zero
    intro l _
    by_cases h : len l < len j
    · rw [hB.tri l j h, mul_zero]
    · rw [hA.tri i l (by omega), zero_mul]

lemma subst_filt (ga gb : Matrix (Fin 2) (Fin 2) k) : Filt (subst ga gb) := by
  constructor <;> intro i j h <;> fin_cases i <;> fin_cases j <;> simp_all [subst, len]

/-- the product of the path algebra, coordinatewise -/
lemma mul_explicit (v w : Fin {d} → k) : (alg k).mul v w = ![{", ".join(lean(e) for e in mulv)}] := by
  funext m; fin_cases m <;> simp [CoordAlg.mul, CoordAlg.ofInt, c, Fin.sum_univ_succ] <;> ring

/-- `u = 1 + x` and its inverse `Σ (−x)^n`, for `x ∈ rad` (coordinates of the arrows and paths) -/
def uU (x : Fin {d} → k) : Fin {d} → k := ![{", ".join(lean(e) for e in u)}]
def uI (x : Fin {d} → k) : Fin {d} → k := ![{", ".join(lean(e) for e in ui)}]

lemma uU_uI (x : Fin {d} → k) : (alg k).mul (uU x) (uI x) = {vec(one)} ∧
    (alg k).mul (uI x) (uU x) = {vec(one)} := by
  constructor <;> rw [mul_explicit] <;> funext m <;> fin_cases m <;> simp [uU, uI] <;> ring

/-- the matrix of `Inn(1 + x) : y ↦ u y u⁻¹` -/
def innM (x : Fin {d} → k) : Matrix (Fin {d}) (Fin {d}) k := !![{";\n    ".join(", ".join(lean(Inn[i, j]) for j in range(d)) for i in range(d))}]

lemma innM_conj (x y : Fin {d} → k) : innM x *ᵥ y = (alg k).mul ((alg k).mul (uU x) y) (uI x) := by
  rw [mul_explicit, mul_explicit]; funext m; fin_cases m <;>
    simp [innM, uU, uI, Matrix.mulVec, dotProduct, Fin.sum_univ_succ] <;> ring

lemma innM_filt (x : Fin {d} → k) : Filt (innM x) := by
  constructor <;> intro i j h <;> fin_cases i <;> fin_cases j <;> simp_all [innM, len]

/-- if `Φ P = P' Φ` with `Φ` filtered and row `v` of `P'` vanishes from length `len c` on, then `P v c = 0` -/
lemma row_zero {{Φ A B : Matrix (Fin {d}) (Fin {d}) k}} (hΦ : Filt Φ) (hc : Φ * A = B * Φ) {{v c : Fin {d}}}
    (hv : len v = 0) (hB : ∀ l, len c ≤ len l → B v l = 0) : A v c = 0 := by
  have h := congrFun (congrFun hc v) c
  rw [Matrix.mul_apply, Matrix.mul_apply, Finset.sum_eq_single v
    (fun l _ hl => by rw [hΦ.vert v l hv, if_neg (fun e => hl e.symm), zero_mul])
    (fun h => absurd (Finset.mem_univ v) h), hΦ.vert v v hv, if_pos rfl, one_mul] at h
  rw [h]; apply Finset.sum_eq_zero; intro l _
  by_cases hl : len l < len c
  · rw [hΦ.tri l c hl, mul_zero]
  · rw [hB l (by omega), zero_mul]

variable [Nontrivial k]
'''
    return txt, ln, nv, d

def depth_info(P, ln, nv, d):
    res = {}
    for v in range(nv):
        row = [sp.expand(P[v, j]) for j in range(d)]
        dv = max([ln[j] for j in range(d) if row[j] != 0], default=-1)
        full = [j for j in range(d) if dv >= 0 and ln[j] == dv]
        cs = [j for j in full if row[j] in (1, -1)] or [j for j in full if row[j] in (tsym, -tsym)]
        assert dv < 0 or cs, (v, row)
        res[v] = (dv, cs[0] if cs else None)
    return res

def pair_lemmas(names, mats, ln, nv, d, sfx, tvar_hyps):
    """names[i]: Lean name of the Poly table; mats[i]: sympy matrix in t"""
    DI = [depth_info(P, ln, nv, d) for P in mats]
    out = []
    for a, b in itertools.permutations(range(len(names)), 2):
        vs = [v for v in range(nv) if DI[a][v][0] != DI[b][v][0]]
        assert vs, (names[a], names[b])
        v = vs[0]
        if DI[a][v][0] > DI[b][v][0]: A, B, c, hyp = a, b, DI[a][v][1], 'hΦ hc'
        else: A, B, c, hyp = b, a, DI[b][v][1], "hΨ hc'"
        out.append(f"lemma {names[a]}_{names[b]}_{sfx} {{Φ Ψ : Matrix (Fin {d}) (Fin {d}) k}} (hΦ : Filt Φ) (hΨ : Filt Ψ) {{t t' : k}}{tvar_hyps}\n"
                   f"    (hc : Φ * evalM {names[a]} t = evalM {names[b]} t' * Φ) (hc' : Ψ * evalM {names[b]} t' = evalM {names[a]} t * Ψ) :\n"
                   f"    False := by\n"
                   f"  have h := row_zero {hyp} (v := {v}) (c := {c}) rfl\n"
                   f"    (fun l hl => by fin_cases l <;> simp [len] at hl ⊢ <;> simp [evalM, {names[B]}])\n"
                   f"  simp [evalM, {names[A]}] at h" + ("\n  all_goals first | exact ht h | exact ht' h\n" if tvar_hyps else "\n"))
    return out

def packaged(thm, doc, idx_name, names, d, n, sfx, tvar_hyps, targs):
    pre = f"/-- {doc} -/\ndef {idx_name} : Fin {n} → Fin {d} → Fin {d} → Poly :=\n  ![" + ", ".join(names) + "]\n\n" \
        if doc else ""
    return pre + f'''theorem {thm} {{Φ Ψ : Matrix (Fin {d}) (Fin {d}) k}} (hΦ : Filt Φ) (hΨ : Filt Ψ)
    (h1 : Φ * Ψ = 1) (h2 : Ψ * Φ = 1) {{t t' : k}}{tvar_hyps} (i j : Fin {n})
    (hc : Φ * evalM ({idx_name} i) t = evalM ({idx_name} j) t' * Φ) : i = j := by
  have hc' : Ψ * evalM ({idx_name} j) t' = evalM ({idx_name} i) t * Ψ := by
    calc Ψ * evalM ({idx_name} j) t' = Ψ * evalM ({idx_name} j) t' * (Φ * Ψ) := by rw [h1, Matrix.mul_one]
      _ = Ψ * (evalM ({idx_name} j) t' * Φ) * Ψ := by simp only [Matrix.mul_assoc]
      _ = Ψ * (Φ * evalM ({idx_name} i) t) * Ψ := by rw [hc]
      _ = (Ψ * Φ) * evalM ({idx_name} i) t * Ψ := by simp only [Matrix.mul_assoc]
      _ = evalM ({idx_name} i) t * Ψ := by rw [h2, Matrix.one_mul]
  fin_cases i <;> fin_cases j <;> simp only [{idx_name}, Fin.isValue, Matrix.cons_val, Fin.reduceFinMk] at hc hc'
''' + "\n".join(("  · rfl" if a == b else f"  · exact ({names[a]}_{names[b]}_{sfx} hΦ hΨ{targs} hc hc').elim")
                for a in range(n) for b in range(n)) + "\n"

# ---------------- 1⇉2⇉3
txt, ln, nv, d = generic('two_pairs')
ks = sorted(R)
names = [f"fam{k}P" for k in ks]
L = ['''/-
# The nineteen families of 1⇉2⇉3 are nineteen behaviours (v2.39)

`Filt Φ`: vertex rows of `Φ` are unit rows and `Φ` is triangular for path length (`ε ∘ Φ = ε`, `Φ (rad^n) ⊆ rad^n`).
Class substitutions and inner automorphisms `Inn(1 + x)`, `x ∈ rad`, are `Filt`, and `Filt` is closed under products,
so the group they generate, which is Aut(kQ) for 1⇉2⇉3 (v2.33), consists of `Filt` matrices. For mutually inverse
`Filt` matrices `Φ, Ψ`, `Φ · P^{(i)}_t = P^{(j)}_{t'} · Φ` forces `i = j`: the depths of the vertex rows differ.

Generated by lean/gen/gen_behaviour.py; do not edit by hand.
-/
import EatTheFibre.Gen.two_pairs_disjoint

open Matrix EatTheFibre.Reflect

namespace EatTheFibre.Q_two_pairs
''', txt]
L += pair_lemmas(names, [famP(k, tsym) for k in ks], ln, nv, d, 'behav', '')
L.append(packaged('fam_behaviours_distinct',
                  None, 'famP', names, d, len(ks), 'behav', '', '').replace(
    "theorem fam_behaviours_distinct",
    "/-- **The nineteen families of 1⇉2⇉3 are nineteen behaviours.** If mutually inverse filtered matrices (for example\n"
    "products of class substitutions and inner automorphisms) intertwine member `t` of family `i` with member `t'` of\n"
    "family `j`, then `i = j`. -/\ntheorem fam_behaviours_distinct"))
L.append("end EatTheFibre.Q_two_pairs\n")
open(os.path.join(OUT, 'two_pairs_behaviour.lean'), 'w').write("\n".join(L))
print("wrote two_pairs_behaviour.lean")

# ---------------- 1⇉2⇇3
txt, ln, nv, d = generic('pairs_in')
src = open(os.path.join(OUT, 'pairs_in.lean')).read()
def table(name):
    m = re.search(rf"def {name} : Fin {d} → Fin {d} → Poly :=\s*(.*?)\n\s*\n", src, re.S)
    T = ast.literal_eval(m.group(1).replace("![", "["))
    return sp.Matrix(d, d, lambda i, j: sum(c_ * tsym ** e for e, c_ in enumerate(T[i][j])))
names = [f"thm3F{F}P" for F in (1, 2, 3, 4)]
L = ['''/-
# Theorem 3's four families on 1⇉2⇇3 are four behaviours modulo Inn ⋊ Sub (v2.39)

Same method as `two_pairs_behaviour.lean`. F1 and F3 degenerate at `t = 0` (v2.33), where a vertex row loses depth, so
the parameters are assumed nonzero. The source swap (an automorphism of Q, not in Inn ⋊ Sub) exchanges F1 and F3
(v2.38, computed), so modulo all of Aut(kQ) there are three behaviour types.

Generated by lean/gen/gen_behaviour.py; do not edit by hand.
-/
import EatTheFibre.Invariants

open Matrix EatTheFibre.Reflect

namespace EatTheFibre.Q_pairs_in
''', txt]
L += pair_lemmas(names, [table(n) for n in names], ln, nv, d, 'behav', " (ht : t ≠ 0) (ht' : t' ≠ 0)")
L.append(packaged('thm3_behaviours_distinct', "Theorem 3's families F1–F4, indexed by `Fin 4`", 'thm3FP', names, d, 4,
                  'behav', " (ht : t ≠ 0) (ht' : t' ≠ 0)", " ht ht'").replace(
    "theorem thm3_behaviours_distinct",
    "/-- **F1, F2, F3, F4 are four behaviours modulo Inn ⋊ Sub** (nonzero parameters). -/\ntheorem thm3_behaviours_distinct"))
L.append("end EatTheFibre.Q_pairs_in\n")
open(os.path.join(OUT, 'pairs_in_behaviour.lean'), 'w').write("\n".join(L))
print("wrote pairs_in_behaviour.lean")
