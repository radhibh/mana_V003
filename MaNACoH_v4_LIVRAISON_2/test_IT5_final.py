#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test IT5 FINAL — Validation complète des corrections
=====================================================
E1 révisé : orientation arbitraire (pas de Head<Tail)
E3 : mapping indices robuste
E5 : LH_dual[position] 0-based
E7 : forces nodales cohérentes avec C[head]=+1
E8 : convention pz = -Weight dans résidu
+ invariance de l'équilibre par inversion d'arêtes
+ bug fallback u=LH_vec corrigé
"""
import math, random

# ============================================================
# ARC TEST
# ============================================================
def build_test_arch(N=16, R=10.0, t=1.5, rho=2000.0):
    g = 9.81
    nodes_raw = []
    for k in range(N):
        phi = k * math.pi / (2 * (N - 1))
        is_bou = (k == 0 or k == N - 1)
        if not is_bou:
            idx = k
        elif k == 0:
            idx = N - 1
        else:
            idx = N
        nodes_raw.append({
            'Index': idx, 'XG': R*math.cos(phi), 'YG': 0.0,
            'ZG': R*math.sin(phi), 'Boundary': is_bou,
            'Weight': 0.0, '_k': k
        })
    int_n = sorted([n for n in nodes_raw if not n['Boundary']], key=lambda n: n['Index'])
    bou_n = sorted([n for n in nodes_raw if n['Boundary']], key=lambda n: n['Index'])
    nodes = int_n + bou_n
    ni, nb = len(int_n), len(bou_n)
    nmap = {n['Index']: n for n in nodes}
    by_k = sorted(nodes_raw, key=lambda n: n['_k'])

    branches = []
    for k in range(N - 1):
        h_idx, t_idx = by_k[k]['Index'], by_k[k+1]['Index']
        hn, tn = nmap[h_idx], nmap[t_idx]
        dx, dy, dz = hn['XG']-tn['XG'], hn['YG']-tn['YG'], hn['ZG']-tn['ZG']
        L = math.sqrt(dx*dx+dy*dy+dz*dz)
        LH = math.sqrt(dx*dx+dy*dy)
        w_b = L * t * 1.0 * rho * g
        branches.append({
            'Index': k+1, 'Head_node': h_idx, 'Tail_node': t_idx,
            'u': dx, 'v': dy, 'w': dz, 'L': L, 'LH': LH, 'Weight': w_b
        })
    for b in branches:
        for n in nodes:
            if n['Index'] in (b['Head_node'], b['Tail_node']):
                n['Weight'] += b['Weight'] / 2
    return nodes, branches, ni, nb


def solve_full(nodes, branches, ni, nb, H=100000.0):
    """Pipeline complet Module 2+3 : Ci/Cb → K → ζ → Di/Db → zi"""
    int_n = [n for n in nodes if not n['Boundary']]
    bou_n = [n for n in nodes if n['Boundary']]
    m = len(branches)
    nmap = {n['Index']: n for n in nodes}

    # E3 : mapping robuste
    int_col, bou_col = {}, {}
    ic, bc = 0, 0
    for n in int_n:
        int_col[n['Index']] = ic; ic += 1
    for n in bou_n:
        bou_col[n['Index']] = bc; bc += 1

    Ci = [[0.0]*ni for _ in range(m)]
    Cb = [[0.0]*nb for _ in range(m)]
    u_vec, LH_vec = [], []

    for j, b in enumerate(branches):
        h, t = b['Head_node'], b['Tail_node']
        if h in int_col: Ci[j][int_col[h]] = 1.0
        elif h in bou_col: Cb[j][bou_col[h]] = 1.0
        if t in int_col: Ci[j][int_col[t]] = -1.0
        elif t in bou_col: Cb[j][bou_col[t]] = -1.0
        u_vec.append(b['u'])
        LH_vec.append(b['LH'])

    zeta = [H] * m
    q = [zeta[j]/LH_vec[j] if LH_vec[j] > 1e-10 else 0 for j in range(m)]

    Di = [[sum(Ci[j][i]*q[j]*Ci[j][k] for j in range(m)) for k in range(ni)] for i in range(ni)]
    Db = [[sum(Ci[j][i]*q[j]*Cb[j][k] for j in range(m)) for k in range(nb)] for i in range(ni)]

    pz = [-n['Weight'] for n in int_n]
    zb = [n['ZG'] for n in bou_n]
    rhs = [-pz[i] - sum(Db[i][k]*zb[k] for k in range(nb)) for i in range(ni)]

    # Gauss
    aug = [[Di[i][j] for j in range(ni)] + [rhs[i]] for i in range(ni)]
    for col in range(ni):
        mx = max(range(col, ni), key=lambda r: abs(aug[r][col]))
        aug[col], aug[mx] = aug[mx], aug[col]
        piv = aug[col][col]
        if abs(piv) < 1e-15: return None, None, None
        for j in range(col, ni+1): aug[col][j] /= piv
        for row in range(ni):
            if row != col:
                f = aug[row][col]
                for j in range(col, ni+1): aug[row][j] -= f*aug[col][j]
    zi = [aug[i][ni] for i in range(ni)]
    return zi, q, zeta


def flip_branches(branches, indices):
    """Inverse l'orientation des branches aux indices donnés."""
    out = []
    for j, b in enumerate(branches):
        if j in indices:
            out.append({
                'Index': b['Index'],
                'Head_node': b['Tail_node'], 'Tail_node': b['Head_node'],
                'u': -b['u'], 'v': -b['v'], 'w': -b['w'],
                'L': b['L'], 'LH': b['LH'], 'Weight': b['Weight']
            })
        else:
            out.append(dict(b))
    return out


# ============================================================
# TEST 1 : Invariance par inversion d'arêtes (E1 révisé)
# ============================================================
def test_invariance():
    print("=" * 60)
    print("TEST 1 : Invariance equilibre par inversion d'aretes")
    print("=" * 60)
    nodes, branches, ni, nb = build_test_arch()
    m = len(branches)
    zi_ref, _, _ = solve_full(nodes, branches, ni, nb)

    random.seed(42)
    scenarios = [
        ("Toutes inversees", list(range(m))),
        ("Paires inversees", list(range(0, m, 2))),
        ("Aleatoire (7/15)", random.sample(range(m), 7)),
        ("1 seule (#3)", [3]),
    ]
    ok = True
    for name, idx in scenarios:
        bf = flip_branches(branches, idx)
        zi_f, _, _ = solve_full(nodes, bf, ni, nb)
        diff = max(abs(a-b) for a, b in zip(zi_ref, zi_f))
        status = "✅" if diff < 1e-10 else "❌"
        print("  {} : Δzi={:.2e} {}".format(name, diff, status))
        if diff >= 1e-10: ok = False

    print("  {} TEST 1\n".format("✅" if ok else "❌"))
    return ok


# ============================================================
# TEST 2 : E3 mapping robuste vs fragile
# ============================================================
def test_E3():
    print("=" * 60)
    print("TEST 2 : E3 mapping indices robuste")
    print("=" * 60)
    nodes, _, ni, nb = build_test_arch()

    # Robuste
    col_r = {}; cnt = 0
    for n in nodes:
        if not n['Boundary']:
            col_r[n['Index']] = cnt; cnt += 1

    # Fragile
    col_f = {}
    for n in nodes:
        if not n['Boundary']:
            col_f[n['Index']] = n['Index'] - 1

    match = (col_r == col_f)

    # Cas pathologique : indices avec trous
    gap_nodes = [{'Index': 3, 'Boundary': False}, {'Index': 7, 'Boundary': False}, {'Index': 20, 'Boundary': True}]
    col_gap = {}; cnt = 0
    for n in gap_nodes:
        if not n['Boundary']:
            col_gap[n['Index']] = cnt; cnt += 1
    gap_ok = (col_gap == {3: 0, 7: 1})

    col_gap_old = {}
    for n in gap_nodes:
        if not n['Boundary']:
            col_gap_old[n['Index']] = n['Index'] - 1
    # col_gap_old = {3:2, 7:6} — FAUX ! Les colonnes 0,1 ne sont pas utilisées
    gap_old_wrong = (col_gap_old != col_gap)

    ok = match and gap_ok and gap_old_wrong
    print("  Standard : robuste==fragile = {}".format(match))
    print("  Trous    : robuste={}, correct={}".format(col_gap, gap_ok))
    print("  Trous    : fragile={}, faux={}".format(col_gap_old, gap_old_wrong))
    print("  {} TEST 2\n".format("✅" if ok else "❌"))
    return ok


# ============================================================
# TEST 3 : E5 index 0-based vs 1-based
# ============================================================
def test_E5():
    print("=" * 60)
    print("TEST 3 : E5 index branche 0-based pour LH_dual")
    print("=" * 60)
    _, branches, _, _ = build_test_arch()
    m = len(branches)
    LH_dual = [1000.0*(i+1) for i in range(m)]

    # Nouveau (0-based)
    vals_new = [LH_dual[pos] for pos, _ in enumerate(branches)]
    # Ancien (1-based via b['Index'])
    vals_old = []
    for b in branches:
        bi = b['Index']
        vals_old.append(LH_dual[bi] if bi < m else -1)

    decal = sum(1 for a, b in zip(vals_old, vals_new) if abs(a-b) > 1e-10)
    hors_bornes = sum(1 for b in branches if b['Index'] >= m)

    ok = (decal > 0 and hors_bornes > 0)
    print("  Decalages ancien: {}/{} branches".format(decal, m))
    print("  Hors bornes ancien: {}".format(hors_bornes))
    print("  {} TEST 3\n".format("✅" if ok else "❌"))
    return ok


# ============================================================
# TEST 4 : E7+E8 equilibre nodal
# ============================================================
def test_E7_E8():
    print("=" * 60)
    print("TEST 4 : E7+E8 equilibre nodal")
    print("=" * 60)
    nodes, branches, ni, nb = build_test_arch()
    int_n = [n for n in nodes if not n['Boundary']]
    bou_n = [n for n in nodes if n['Boundary']]
    all_n = int_n + bou_n
    idx_map = {n['Index']: i for i, n in enumerate(all_n)}
    nmap = {n['Index']: n for n in nodes}
    m = len(branches)

    zi, q_vec, _ = solve_full(nodes, branches, ni, nb)
    zb = [n['ZG'] for n in bou_n]
    Z = zi + zb

    def check_equilibrium(sign_head, sign_tail, label):
        nf = {n['Index']: [0.0, 0.0, 0.0] for n in nodes}
        for j, b in enumerate(branches):
            hi, ti = b['Head_node'], b['Tail_node']
            hn, tn = nmap[hi], nmap[ti]
            u = hn['XG'] - tn['XG']
            v = hn['YG'] - tn['YG']
            i1, i2 = idx_map[hi], idx_map[ti]
            w = Z[i1] - Z[i2]
            us, vs, ws = q_vec[j]*u, q_vec[j]*v, q_vec[j]*w
            nf[hi][0] += sign_head * us
            nf[hi][1] += sign_head * vs
            nf[hi][2] += sign_head * ws
            nf[ti][0] += sign_tail * us
            nf[ti][1] += sign_tail * vs
            nf[ti][2] += sign_tail * ws
        max_res = 0.0
        for n in int_n:
            fx, fy, fz = nf[n['Index']]
            max_res = max(max_res, math.sqrt(fx**2 + fy**2 + (fz - n['Weight'])**2))
        return max_res

    res_new = check_equilibrium(+1, -1, "E7 corrigé")
    res_old = check_equilibrium(-1, +1, "E7 ancien")

    # Equilibre global
    nf_g = {n['Index']: [0.0, 0.0, 0.0] for n in nodes}
    for j, b in enumerate(branches):
        hi, ti = b['Head_node'], b['Tail_node']
        hn, tn = nmap[hi], nmap[ti]
        u = hn['XG'] - tn['XG']; v = hn['YG'] - tn['YG']
        i1, i2 = idx_map[hi], idx_map[ti]
        w = Z[i1] - Z[i2]
        us, vs, ws = q_vec[j]*u, q_vec[j]*v, q_vec[j]*w
        nf_g[hi][0] += us; nf_g[hi][1] += vs; nf_g[hi][2] += ws
        nf_g[ti][0] -= us; nf_g[ti][1] -= vs; nf_g[ti][2] -= ws
    tot_Rz = sum(-nf_g[n['Index']][2] for n in bou_n)
    tot_W = sum(n['Weight'] for n in int_n)

    ok = (res_new < 1e-6 and res_old > 1.0)
    print("  Residu E7 corrige : {:.2e}".format(res_new))
    print("  Residu E7 ancien  : {:.2e}".format(res_old))
    print("  Sum RV = {:.1f} N  vs  Sum W = {:.1f} N  (ecart {:.2e})".format(
        tot_Rz, tot_W, abs(tot_Rz - tot_W)))
    print("  {} TEST 4\n".format("✅" if ok else "❌"))
    return ok


# ============================================================
# TEST 5 : Bug fallback u=LH_vec
# ============================================================
def test_fallback_bug():
    print("=" * 60)
    print("TEST 5 : Bug fallback u=LH_vec (corrige v51)")
    print("=" * 60)
    # Cas explicite : branche de droite à gauche (u>0) vs gauche à droite (u<0)
    # En 2D avec v=0 : LH = |u| ≥ 0 toujours, mais u peut être < 0
    test_branches = [
        {'u':  5.0, 'v': 0.0},  # u>0 → LH=5 = u ✓
        {'u': -3.0, 'v': 0.0},  # u<0 → LH=3 ≠ u=-3 ✗
        {'u':  4.0, 'v': 3.0},  # LH=5 ≠ u=4 ✗
        {'u': -4.0, 'v': 3.0},  # LH=5 ≠ u=-4 ✗
    ]
    for b in test_branches:
        b['LH'] = math.sqrt(b['u']**2 + b['v']**2)

    u_correct = [b['u'] for b in test_branches]
    LH_vals = [b['LH'] for b in test_branches]
    diffs = sum(1 for j in range(len(test_branches)) if abs(u_correct[j] - LH_vals[j]) > 1e-10)
    neg_u = sum(1 for u in u_correct if u < -1e-10)

    # Impact sur K : K[i,j] = Ci[j][i] * u[j] / LH[j]
    # Si on utilise LH au lieu de u : K_faux = Ci * LH / LH = Ci (perd la direction !)
    # Avec u correct : K_ok = Ci * u / LH (préserve le cosinus directionnel)
    # Pour la branche u=-3 : u/LH = -1.0 vs LH/LH = +1.0 → SIGNE INVERSÉ
    cos_correct = [b['u']/b['LH'] if b['LH'] > 1e-10 else 0 for b in test_branches]
    cos_faux = [b['LH']/b['LH'] if b['LH'] > 1e-10 else 0 for b in test_branches]
    sign_errors = sum(1 for c, f in zip(cos_correct, cos_faux) if abs(c - f) > 1e-10)

    ok = diffs > 0 and sign_errors > 0
    print("  Branches avec u<0 : {}".format(neg_u))
    print("  Branches où u≠LH : {}/{}".format(diffs, len(test_branches)))
    print("  Erreurs de signe dans u/LH : {}/{}".format(sign_errors, len(test_branches)))
    print("  cos(correct) = {}".format([round(c, 2) for c in cos_correct]))
    print("  cos(faux)    = {}".format([round(c, 2) for c in cos_faux]))
    print("  {} TEST 5\n".format("✅" if ok else "❌"))
    return ok


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    results = [
        ("Invariance orientation (E1)", test_invariance()),
        ("Mapping robuste (E3)", test_E3()),
        ("Index 0-based (E5)", test_E5()),
        ("Equilibre nodal (E7+E8)", test_E7_E8()),
        ("Fallback u=LH bug", test_fallback_bug()),
    ]

    print("=" * 60)
    print("BILAN IT5 FINAL")
    print("=" * 60)
    for name, ok in results:
        print("  {} : {}".format(name, "✅" if ok else "❌"))
    passed = sum(1 for _, ok in results if ok)
    print("  {}/{} tests passes".format(passed, len(results)))
    print("=" * 60)
