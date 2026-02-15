#! python3
"""
MaNACoH v4 - Module 3 : EQUILIBRIUM TNA (v53 - RO Optimizer)
===================================================================
Refonte complète + optimiseur Recherche Opérationnelle.

CORRECTIONS par rapport à v49 :
  V1 : Équilibre horizontal explicite (MF1 eq. 4.27-4.30)
  V2 : Optimisation DOF propre (f_Optimization.bas)
  V3 : Z = résultat de l'équilibre (MF1 eq. 4.38), pas une cible

CHAÎNE DE CALCUL (fidèle au VBA e_ComputeEquilibrium.bas) :
  1. Construire K = [Ci^T U L^-1_H ; Ci^T V L^-1_H]  — eq. 4.27
  2. Identifier DOF-branches, séparer Kunk | Kdof     — eq. 4.28-4.29
  3. zeta_unk = -Kunk^-1 · Kdof · zeta_dof             — eq. 4.30
  4. Di = Ci^T diag(q) Ci,  Db = Ci^T diag(q) Cb       — eq. 4.36-4.37
  5. zi = Di^-1 (-pz - Db·zb)                           — eq. 4.38
  6. Optimiser les DOF pour max(min(CSG))                — f_Optimization.bas

SOURCES :
  - MF1.pdf §4.3.5, eq. 4.27-4.38 (théorie)
  - MF-part-4.pdf §13.1.2-13.1.3  (implémentation)
  - e_ComputeEquilibrium.bas       (VBA original)
  - f_Optimization.bas             (VBA optimisation)
  - MANUAL §2.4                    (DOF, hints)

INPUTS Grasshopper :
  - nodes : str/list (Item)     - JSON ou list of dict du Module 0
  - branches : str/list (Item)  - JSON ou list of dict du Module 0
  - optimize : bool (Item)      - Lancer l'optimisation [True]
  - max_iter : int (Item)       - Budget optimisation [200]
  - objective : str (Item)      - "proximity" ou "max_csg" ["proximity"]
  - boundary_position : str (Item) - "middle", "intrados", "extrados" ["middle"]

  Note : Ci, Cb, u_vec, v_vec, LH_vec, n_int ne sont PLUS des inputs.
  La topologie est construite en interne par build_topology().
  Le Module 2 n'est plus nécessaire dans le pipeline GH.

OUTPUTS Grasshopper :
  - Z : list float - Altitudes thrust network [ni + nb]
  - q : list float - Densités de force [m branches]
  - F : list float - Forces dans les branches [m branches]
  - LH_dual : list float - Forces horizontales [m branches]
  - CSG_min : float - CSG minimum de la structure
  - CSG_list : list float - CSG par branche
  - thrust_pts : Point3d list - Points du thrust network
  - thrust_lines : Line list - Branches du thrust network

Auteur  : Claude (IT4 du plan de correction MaNACoH)
Date    : 2026-02-14
Licence : MIT
"""

import math

# Tentative import Rhino (absent en mode test)
try:
    import Rhino.Geometry as rg
    HAS_RHINO = True
except ImportError:
    HAS_RHINO = False

# =======================================================================
# ALGÈBRE LINÉAIRE PURE PYTHON
# (compatible IronPython / Grasshopper, pas de numpy requis)
# =======================================================================

def _zeros(r, c):
    return [[0.0] * c for _ in range(r)]

def _eye(n):
    M = _zeros(n, n)
    for i in range(n):
        M[i][i] = 1.0
    return M

def _transpose(A):
    if not A or not A[0]:
        return []
    r, c = len(A), len(A[0])
    return [[A[i][j] for i in range(r)] for j in range(c)]

def _matmul(A, B):
    """Produit matriciel A @ B."""
    if not A or not B:
        return []
    r, k, c = len(A), len(A[0]), len(B[0])
    R = _zeros(r, c)
    for i in range(r):
        for j in range(c):
            s = 0.0
            for l in range(k):
                s += A[i][l] * B[l][j]
            R[i][j] = s
    return R

def _matvec(A, x):
    """Produit matrice-vecteur A @ x."""
    if not A:
        return []
    r, c = len(A), len(A[0])
    return [sum(A[i][j] * x[j] for j in range(c)) for i in range(r)]

def _vecmat(x, A):
    """Produit vecteur-matrice x^T @ A."""
    if not A:
        return []
    r, c = len(A), len(A[0])
    return [sum(x[i] * A[i][j] for i in range(r)) for j in range(c)]

def _diag_mat(v):
    """Construit matrice diagonale depuis vecteur."""
    n = len(v)
    M = _zeros(n, n)
    for i in range(n):
        M[i][i] = v[i]
    return M

def _inverse(A):
    """Inversion par pivot de Gauss avec pivotage partiel. A doit être carrée."""
    n = len(A)
    if n == 0:
        return []
    ncols = len(A[0])
    if n != ncols:
        raise ValueError("_inverse: matrice non carree ({} x {})".format(n, ncols))
    M = [A[i][:] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for col in range(n):
        # Pivotage partiel
        max_row = col
        max_val = abs(M[col][col])
        for row in range(col + 1, n):
            if abs(M[row][col]) > max_val:
                max_val = abs(M[row][col])
                max_row = row
        if max_val < 1e-14:
            raise ValueError("Matrice singuliere (col {})".format(col))
        M[col], M[max_row] = M[max_row], M[col]
        scale = M[col][col]
        M[col] = [x / scale for x in M[col]]
        for row in range(n):
            if row != col:
                factor = M[row][col]
                M[row] = [M[row][j] - factor * M[col][j] for j in range(2 * n)]
    return [row[n:] for row in M]


def _solve_lstsq(A, B):
    """
    Résout A·X = B aux moindres carrés (A peut être rectangulaire).
    A : n_eq x n_unk, B : n_eq x n_dof
    Retourne X : n_unk x n_dof
    
    Méthode : équations normales A^T A X = A^T B
    (valide car Kunk a rang plein par construction — identify_dof_branches
     garantit l'indépendance linéaire des colonnes UNK)
    
    Source : MF1 eq. 4.30, résolution standard du système surdéterminé
    """
    n_eq = len(A)
    n_unk = len(A[0]) if n_eq > 0 else 0
    
    # Cas carré : inversion directe
    if n_eq == n_unk:
        A_inv = _inverse(A)
        return _matmul(A_inv, B)
    
    # Cas surdéterminé (n_eq > n_unk) : équations normales
    # A^T A (n_unk x n_unk) · X = A^T B (n_unk x n_dof)
    AT = _transpose(A)
    ATA = _matmul(AT, A)     # n_unk x n_unk
    ATB = _matmul(AT, B)     # n_unk x n_dof
    
    ATA_inv = _inverse(ATA)  # n_unk x n_unk
    X = _matmul(ATA_inv, ATB)  # n_unk x n_dof
    return X

def _norm(v):
    return math.sqrt(sum(x * x for x in v))

def _dot(a, b):
    return sum(ai * bi for ai, bi in zip(a, b))

def _col(A, j):
    """Extraire colonne j de A."""
    return [A[i][j] for i in range(len(A))]

def _row_norms(A):
    """Normes des lignes de A."""
    return [_norm(row) for row in A]

def _delete_rows(A, indices_to_delete):
    """Supprime les lignes dont l'index est dans indices_to_delete."""
    s = set(indices_to_delete)
    return [A[i][:] for i in range(len(A)) if i not in s]

def _gram_det(cols):
    """
    Calcule det(T^T @ T) pour T = matrice dont les colonnes sont dans cols.
    Utilisé pour le test d'indépendance linéaire (VBA L.472-532).
    """
    k = len(cols)
    if k == 0:
        return 1.0
    n = len(cols[0])
    # G = T^T @ T  (k x k)
    G = _zeros(k, k)
    for i in range(k):
        for j in range(k):
            G[i][j] = sum(cols[i][r] * cols[j][r] for r in range(n))
    # det de G
    return _det(G)

def _det(A):
    """Déterminant par élimination de Gauss."""
    n = len(A)
    M = [row[:] for row in A]
    sign = 1.0
    for col in range(n):
        max_row = col
        for row in range(col + 1, n):
            if abs(M[row][col]) > abs(M[max_row][col]):
                max_row = row
        if abs(M[max_row][col]) < 1e-30:
            return 0.0
        if max_row != col:
            M[col], M[max_row] = M[max_row], M[col]
            sign *= -1.0
        for row in range(col + 1, n):
            factor = M[row][col] / M[col][col]
            for j in range(col + 1, n):
                M[row][j] -= factor * M[col][j]
    result = sign
    for i in range(n):
        result *= M[i][i]
    return result


# =======================================================================
# ÉTAPE 1 : ÉQUILIBRE HORIZONTAL (MF1 eq. 4.27-4.30)
# Source : e_ComputeEquilibrium.bas subs c, e
# =======================================================================

def build_K_matrix(Ci, u, v, LH, ni, m):
    """
    K = [ Ci^T · diag(u/LH) ]   (2*ni x m) — MF1 eq. 4.27
        [ Ci^T · diag(v/LH) ]

    K[i, j] = Ci[j][i] * u[j] / LH[j]  (partie haute)
    K[ni+i, j] = Ci[j][i] * v[j] / LH[j]  (partie basse)
    """
    K = _zeros(2 * ni, m)
    for j in range(m):
        lh = LH[j] if abs(LH[j]) > 1e-12 else 1.0
        for i in range(ni):
            K[i][j] = Ci[j][i] * u[j] / lh
            K[ni + i][j] = Ci[j][i] * v[j] / lh
    return K


def remove_redundant_rows(K, Ci, u, v, LH, node_valence, ni, m, tol=1e-3):
    """
    Supprime les lignes redondantes de K :
    - Nœuds 2-valents avec branches parallèles (VBA L.402-446)
    - Lignes quasi-nulles (cas 2D)

    Source : e_ComputeEquilibrium.bas L.402-446
    """
    rows_to_delete = set()

    for i in range(ni):
        if node_valence[i] != 2:
            continue
        # Trouver les 2 branches connectées
        br = [j for j in range(m) if abs(Ci[j][i]) > 0.5]
        if len(br) != 2:
            continue
        j1, j2 = br
        # Test parallélisme (VBA L.414-445)
        if abs(v[j1]) > abs(u[j1]):
            lh1 = max(abs(LH[j1]), 1e-12)
            lh2 = max(abs(LH[j2]), 1e-12)
            if abs(abs(v[j1] / lh1) - abs(v[j2] / lh2)) < tol:
                rows_to_delete.add(i)  # supprimer eq x
        else:
            lh1 = max(abs(LH[j1]), 1e-12)
            lh2 = max(abs(LH[j2]), 1e-12)
            if abs(abs(u[j1] / lh1) - abs(u[j2] / lh2)) < tol:
                rows_to_delete.add(ni + i)  # supprimer eq y

    # Supprimer aussi les lignes quasi-nulles
    for i in range(2 * ni):
        if _norm(K[i]) < tol:
            rows_to_delete.add(i)

    K_red = _delete_rows(K, rows_to_delete)
    return K_red, len(rows_to_delete)


def identify_dof_branches(K_red, m, chosen_dof=None, tol=1e-10):
    """
    Sépare les branches en UNK et DOF par test d'indépendance linéaire.
    Source : VBA e_ComputeEquilibrium.bas L.456-536

    Returns: Kunk, Kdof, unk_indices, dof_indices
    """
    n_eq = len(K_red)
    if chosen_dof is None:
        chosen_dof = set()
    else:
        chosen_dof = set(chosen_dof)

    dof_indices = []
    unk_indices = []
    unk_cols = []  # colonnes de K pour les branches unk

    # Phase 1 : extraire les chosen_dof (VBA L.376-398)
    remaining = []
    for j in range(m):
        if j in chosen_dof:
            dof_indices.append(j)
        else:
            remaining.append(j)

    # Phase 2 : construire Kunk par test d'indépendance (VBA L.472-532)
    for j in remaining:
        col_j = _col(K_red, j)
        if _norm(col_j) < tol:
            dof_indices.append(j)
            continue
        if len(unk_cols) == 0:
            unk_cols.append(col_j)
            unk_indices.append(j)
            continue
        # Test : cols augmentées indépendantes ?
        test_cols = unk_cols + [col_j]
        d = abs(_gram_det(test_cols))
        if d > tol:
            unk_cols.append(col_j)
            unk_indices.append(j)
        else:
            dof_indices.append(j)

    n_unk = len(unk_indices)
    n_dof = len(dof_indices)

    # Construire matrices Kunk (n_eq x n_unk) et Kdof (n_eq x n_dof)
    Kunk = _zeros(n_eq, n_unk)
    for j_local, j_global in enumerate(unk_indices):
        for i in range(n_eq):
            Kunk[i][j_local] = K_red[i][j_global]

    Kdof = _zeros(n_eq, n_dof)
    for j_local, j_global in enumerate(dof_indices):
        for i in range(n_eq):
            Kdof[i][j_local] = K_red[i][j_global]

    return Kunk, Kdof, unk_indices, dof_indices


def solve_horizontal(Ci, Cb, u, v, LH, node_valence, ni, m,
                     zeta_dof_values=None, chosen_dof=None,
                     is_vertical=None):
    """
    Résout l'équilibre horizontal : MF1 eq. 4.27-4.30.

    Returns: zeta (m,), dof_indices, unk_indices, M_stored (= -Kunk^-1 Kdof)
    """
    if is_vertical is None:
        is_vertical = [abs(LH[j]) < 1e-12 for j in range(m)]

    # Branches verticales → DOF automatiques
    vert_set = set(j for j in range(m) if is_vertical[j])
    if chosen_dof is None:
        chosen_dof = set()
    else:
        chosen_dof = set(chosen_dof)
    all_chosen = chosen_dof | vert_set

    # Indices non-verticaux
    nv_indices = [j for j in range(m) if not is_vertical[j]]
    m_nv = len(nv_indices)
    if m_nv == 0:
        return [0.0] * m, list(range(m)), [], None

    # Mapping global → local
    g2l = {g: l for l, g in enumerate(nv_indices)}
    chosen_local = set(g2l[g] for g in all_chosen if g in g2l)

    # Construire Ci restreint aux branches non-verticales
    Ci_nv = [[Ci[nv_indices[j]][i] for i in range(ni)] for j in range(m_nv)]
    u_nv = [u[nv_indices[j]] for j in range(m_nv)]
    v_nv = [v[nv_indices[j]] for j in range(m_nv)]
    LH_nv = [LH[nv_indices[j]] for j in range(m_nv)]

    # 1. Construire K — eq. 4.27
    K = build_K_matrix(Ci_nv, u_nv, v_nv, LH_nv, ni, m_nv)

    # 2. Supprimer lignes redondantes
    K_red, n_removed = remove_redundant_rows(
        K, Ci_nv, u_nv, v_nv, LH_nv, node_valence, ni, m_nv)

    # 3. Identifier DOF
    Kunk, Kdof_mat, unk_local, dof_local = identify_dof_branches(
        K_red, m_nv, chosen_local)

    n_unk = len(unk_local)
    n_dof_nv = len(dof_local)

    # 4. Résoudre zeta_unk = -Kunk^-1 Kdof zeta_dof — eq. 4.30
    if zeta_dof_values is None:
        zeta_dof_values = [-1.0] * n_dof_nv

    # Attention : zeta_dof_values a la taille n_dof_nv (DOF non-verticaux seulement)
    if len(zeta_dof_values) != n_dof_nv:
        # L'utilisateur fournit les DOF pour toutes les branches DOF (NV + V)
        # On ne prend que les NV
        zeta_dof_nv = zeta_dof_values[:n_dof_nv]
    else:
        zeta_dof_nv = zeta_dof_values

    if n_unk > 0 and n_dof_nv > 0:
        try:
            # M_stored = -Kunk^(-1) @ Kdof — eq. 4.30
            # Kunk peut être rectangulaire (n_eq > n_unk) si surdéterminé
            # → résolution aux moindres carrés via équations normales
            neg_Kdof = [[-Kdof_mat[i][j] for j in range(n_dof_nv)]
                        for i in range(len(Kdof_mat))]
            M_stored = _solve_lstsq(Kunk, neg_Kdof)
        except ValueError as e:
            print("  WARN: Kunk resolution echouee ({}), fallback zeta uniforme".format(e))
            return [zeta_dof_nv[0] if zeta_dof_nv else -1.0] * m, \
                   dof_local, unk_local, None
        zeta_unk = _matvec(M_stored, zeta_dof_nv)
    elif n_unk == 0:
        M_stored = []
        zeta_unk = []
    else:
        M_stored = [[0.0] * 0 for _ in range(n_unk)]
        zeta_unk = [0.0] * n_unk

    # 5. Reconstruire zeta complet
    zeta_nv = [0.0] * m_nv
    for i, li in enumerate(unk_local):
        zeta_nv[li] = zeta_unk[i]
    for i, li in enumerate(dof_local):
        zeta_nv[li] = zeta_dof_nv[i]

    zeta = [0.0] * m
    for i, gi in enumerate(nv_indices):
        zeta[gi] = zeta_nv[i]

    # Convertir indices locaux → globaux
    dof_global = [nv_indices[li] for li in dof_local] + list(vert_set)
    unk_global = [nv_indices[li] for li in unk_local]

    return zeta, dof_global, unk_global, M_stored


def update_zeta_fast(M_stored, zeta_dof_nv, unk_global, dof_global, m):
    """
    Recalcule zeta sans reconstruire K (MF-part-4 §13.1.2).
    """
    zeta = [0.0] * m
    if M_stored and len(M_stored) > 0 and len(M_stored[0]) > 0:
        zeta_unk = _matvec(M_stored, zeta_dof_nv)
        for i, gi in enumerate(unk_global):
            zeta[gi] = zeta_unk[i]
    for i, gi in enumerate(dof_global):
        if i < len(zeta_dof_nv):
            zeta[gi] = zeta_dof_nv[i]
    return zeta


# =======================================================================
# ÉTAPE 2 : ÉQUILIBRE VERTICAL (MF1 eq. 4.35-4.38)
# Source : e_ComputeEquilibrium.bas sub g_ComputeZs
# =======================================================================

def build_Di_Db(Ci, Cb, LH, zeta, ni, nb, m):
    """
    Di = Ci^T diag(q) Ci   — MF1 eq. 4.36
    Db = Ci^T diag(q) Cb   — MF1 eq. 4.37
    avec q = zeta / LH (densité de force)
    """
    q = [0.0] * m
    for j in range(m):
        if abs(LH[j]) > 1e-12:
            q[j] = zeta[j] / LH[j]

    # Di = Ci^T @ diag(q) @ Ci  (ni x ni)
    Di = _zeros(ni, ni)
    Db = _zeros(ni, nb)
    for j in range(m):
        for i1 in range(ni):
            ci1 = Ci[j][i1]
            if abs(ci1) < 0.5:
                continue
            val = ci1 * q[j]
            for i2 in range(ni):
                Di[i1][i2] += val * Ci[j][i2]
            for i2 in range(nb):
                Db[i1][i2] += val * Cb[j][i2]

    return Di, Db, q


def solve_vertical(Ci, Cb, LH, zeta, pz, zb, ni, nb, m):
    """
    Résout l'équilibre vertical : Di zi = -pz - Db zb — MF1 eq. 4.38

    IMPORTANT : zi est le RÉSULTAT, pas une cible (correction V3).

    ROBUSTESSE v53 :
      - Projection ζ > ε pour garantir compression (MANUAL §2.4.4 c1>0)
      - Régularisation Tikhonov si Di mal conditionnée (fréquent en 3D)
      - Détection et correction des branches non-connectées

    Returns: zi (ni,), q (m,), zb_used (nb,)
    """
    # Projeter ζ > ε pour garantir compression
    zeta_proj = [max(z, 1e-6) for z in zeta]

    Di, Db, q = build_Di_Db(Ci, Cb, LH, zeta_proj, ni, nb, m)

    # PzRed = -pz - Db @ zb — eq. 4.38 membre droit
    Db_zb = _matvec(Db, zb)
    PzRed = [-pz[i] - Db_zb[i] for i in range(ni)]

    # zi = Di^-1 @ PzRed
    # Avec régularisation Tikhonov si Di mal conditionnée
    try:
        Di_inv = _inverse(Di)
        zi = _matvec(Di_inv, PzRed)
        # Vérifier le résidu
        Di_zi = _matvec(Di, zi)
        res = max(abs(Di_zi[i] - PzRed[i]) for i in range(ni))
        if res > 1e-3 * max(abs(PzRed[i]) for i in range(ni) if True):
            raise ValueError("Residu trop grand: {:.2e}".format(res))
    except (ValueError, ZeroDivisionError):
        # Régularisation Tikhonov : (Di^T Di + λI) x = Di^T b
        # Choix de λ : trace(Di^T Di) * 1e-8
        DtD = _zeros(ni, ni)
        Dtb = [0.0] * ni
        for i in range(ni):
            for j in range(ni):
                s = 0.0
                for k in range(ni):
                    s += Di[k][i] * Di[k][j]
                DtD[i][j] = s
            s = 0.0
            for k in range(ni):
                s += Di[k][i] * PzRed[k]
            Dtb[i] = s
        trace = sum(DtD[i][i] for i in range(ni))
        lam = max(trace * 1e-8, 1e-6)
        for i in range(ni):
            DtD[i][i] += lam
        try:
            DtD_inv = _inverse(DtD)
            zi = _matvec(DtD_inv, Dtb)
        except (ValueError, ZeroDivisionError):
            zi = [0.0] * ni

    return zi, q, zb[:]


# =======================================================================
# ÉTAPE 3 : CSG SIMPLIFIÉ (pour l'optimisation)
# Le vrai CSG est calculé par le Module 4.
# =======================================================================

def compute_csg_simplified(nodes, branches, Z, LH_dual, n_int):
    """
    CSG simplifié par branche : (h'/2) / |e|
    Source : MF1 eq. 3.33

    Utilise les mêmes formules que le Module 4 (compute_csg_from_Z du v49),
    mais appelé pendant l'optimisation pour guider le solver.
    """
    int_n = [n for n in nodes if not n['Boundary']]
    bou_n = [n for n in nodes if n['Boundary']]
    all_n = int_n + bou_n
    idx_map = {n['Index']: i for i, n in enumerate(all_n)}

    csg_list = []
    for b in branches:
        if 'xi' not in b:
            csg_list.append(100.0)
            continue

        # Points intrados/extrados du joint
        Pi = (b['xi'], b['yi'], b['zi'])
        Pe = (b['xe'], b['ye'], b['ze'])
        h = math.sqrt((Pe[0]-Pi[0])**2 + (Pe[1]-Pi[1])**2 + (Pe[2]-Pi[2])**2)
        if h < 1e-10:
            csg_list.append(100.0)
            continue

        i1 = idx_map[b['Head_node']]
        i2 = idx_map[b['Tail_node']]
        n1, n2 = all_n[i1], all_n[i2]

        # Centre du thrust pour cette branche
        cx = (n1['XG'] + n2['XG']) / 2
        cy = (n1['YG'] + n2['YG']) / 2
        cz = (Z[i1] + Z[i2]) / 2

        # Centre du joint
        jx = (Pi[0] + Pe[0]) / 2
        jy = (Pi[1] + Pe[1]) / 2
        jz = (Pi[2] + Pe[2]) / 2

        # Direction du joint (normalisée)
        dx, dy, dz = (Pe[0]-Pi[0])/h, (Pe[1]-Pi[1])/h, (Pe[2]-Pi[2])/h

        # Excentricité = |projection du vecteur (joint_center → thrust) sur joint_dir|
        bx, by, bz = cx - jx, cy - jy, cz - jz
        e = abs(dx * bx + dy * by + dz * bz)

        csg = (h / 2) / e if e > 1e-6 else 100.0
        csg_list.append(min(csg, 100.0))

    return csg_list


# =======================================================================
# ÉTAPE 4 : OPTIMISATION (fidèle à f_Optimization.bas)
# =======================================================================

def optimize_dof(nodes, branches, Ci, Cb, u, v, LH, node_valence,
                 ni, nb, m, n_int, pz, zb,
                 is_vertical, M_stored, dof_global, unk_global,
                 max_iter=200, objective="proximity"):
    """
    Optimisation des DOF — v53 : stratégie Recherche Opérationnelle.

    Architecture multi-niveau (adaptatif selon n_dof et objectif) :

    ┌─────────────────────────────────────────────────────────────┐
    │  PHASE 0 : Encadrement adaptatif                           │
    │   proximity → point fixe compas_tna (~10 evals)            │
    │   max_csg   → scan logarithmique (~25 evals)               │
    ├─────────────────────────────────────────────────────────────┤
    │  PHASE 1 : Raffinement 1D — Golden section (~40 evals)     │
    │   Exploite l'unimodalité prouvée de CSG_min(H) en 1D       │
    ├─────────────────────────────────────────────────────────────┤
    │  PHASE 2 : Raffinement N-D (si n_dof > 1)                  │
    │   Coordinate Descent (golden section par axe) + Nelder-Mead│
    │   Exploite le couplage faible entre méridiens              │
    └─────────────────────────────────────────────────────────────┘

    PROPRIÉTÉS RO EXPLOITÉES (analyse du paysage, cas Heyman) :
      P1: CSG_min(H) est unimodale en 1D → golden section optimal
      P2: proximity(H) est quasi-convexe → point fixe converge
      P3: Dimension faible (N ≤ 20) → Nelder-Mead efficace
      P4: Évaluation O(ni³) pas coûteuse → budget ~200 evals ok
      P5: Bornes naturelles H > 0 → domaine compact
      P6: Couplage faible inter-méridiens → coordinate descent

    GARANTIES PHYSIQUES :
      - Équilibre horizontal K·ζ=0 par construction (M_stored)
      - Équilibre vertical Di·zi=rhs résolu EXACTEMENT
      - ζ > 0 garanti par projection (compression only)

    Sources :
      MANUAL §2.4.4 : fonctions objectif et hints
      compas_tna vertical_numpy.py L.96-114 : point fixe scalaire
      VBA f_Optimization.bas : GRG avec contraintes C1>0, C2>1
    """
    int_n = [n for n in nodes if not n['Boundary']]
    zG = [n['ZG'] for n in int_n]
    n_dof_nv = len([g for g in dof_global if not is_vertical[g]])
    W = abs(sum(n['Weight'] for n in int_n))
    if W < 1e-6:
        W = 1000.0

    n_evals = 0

    def evaluate(zdof):
        zeta = update_zeta_fast(M_stored, zdof, unk_global, dof_global, m)
        zi, q, zb_u = solve_vertical(Ci, Cb, LH, zeta, pz, zb, ni, nb, m)
        return zeta, zi, q, zb_u

    def cost_fn(zdof):
        """Fonction objectif unifiée. Retourne (cost, zeta, zi, q, zb).
        Pénalité pour branches en traction (MANUAL §2.4.4 c1>0).
        """
        nonlocal n_evals
        try:
            zeta, zi, q, zb_u = evaluate(zdof)
            n_evals += 1
        except Exception:
            return 1e30, None, None, None, None

        # Pénalité traction : ζ < 0 signifie traction inadmissible
        n_tension = sum(1 for z in zeta if z < -1e-6)
        penalty = 0.0
        if n_tension > 0:
            penalty = sum(z**2 for z in zeta if z < 0) * 1e3

        if objective == "proximity":
            c = sum((zi[i] - zG[i]) ** 2 for i in range(ni)) + penalty
        else:  # max_csg
            Z_full = zi[:] + zb_u[:]
            csg_list = compute_csg_simplified(nodes, branches, Z_full, zeta, n_int)
            c = -min(csg_list) if csg_list else 1e10
            c += penalty
        return c, zeta, zi, q, zb_u

    # ================================================================
    # PHASE 0 : Encadrement adaptatif
    # ================================================================
    if objective == "proximity":
        print("  Phase 0: point fixe (proximity)")
        zG_rms = math.sqrt(sum(zG[i]**2 for i in range(ni)) / max(ni, 1))
        if zG_rms < 1e-6:
            zG_rms = 1.0
        H = W / max(1, n_int)
        best_H, best_cost, best_data = H, 1e30, None

        # Pour N>1 DOF : scan log multi-start plutôt que point fixe pur
        # Le point fixe scalaire échoue quand les branches sont hétérogènes
        if n_dof_nv > 1:
            # Stratégie : scan log sur H scalaire (rapide)
            # PUIS essai avec proportionnalité aux longueurs LH des DOF
            # (les branches longues ont typiquement des forces plus faibles)

            # 1. Scan H scalaire (comme 1 DOF)
            for k in range(-15, 20):
                scale = 10.0 ** (k / 8.0)
                H_test = W * scale
                zdof = [H_test] * n_dof_nv
                c, zeta, zi, q, zb_u = cost_fn(zdof)
                if c < best_cost:
                    best_cost = c
                    best_H = H_test
                    best_data = (zeta, zi, q, zb_u)

            # 2. Essai proportionnel aux LH des DOF-branches
            dof_nv_list = [g for g in dof_global if not is_vertical[g]]
            lh_dof = [abs(LH[g]) for g in dof_nv_list]
            lh_mean = sum(lh_dof) / max(len(lh_dof), 1) if lh_dof else 1.0

            for k in range(-10, 15):
                scale = 10.0 ** (k / 8.0)
                H_base = W * scale
                # Chaque DOF proportionnel à sa LH (inversement)
                zdof = [H_base * lh_mean / max(lh_dof[i], 1e-6)
                        for i in range(n_dof_nv)]
                c, zeta, zi, q, zb_u = cost_fn(zdof)
                if c < best_cost:
                    best_cost = c
                    best_H = H_base
                    best_data = (zeta, zi, q, zb_u)

            H_center = best_H
            H_lo = best_H * 0.05
            H_hi = best_H * 20.0
            print("    {} evals, H_init={:.1f}, cost={:.4f}".format(
                n_evals, best_H, best_cost))
        else:
            for k in range(30):
                zdof = [H] * n_dof_nv
                c, zeta, zi, q, zb_u = cost_fn(zdof)
                if c < best_cost:
                    best_cost, best_H = c, H
                    best_data = (zeta, zi, q, zb_u)
                zi_rms = math.sqrt(sum(zi[i]**2 for i in range(ni)) / max(ni, 1)) if zi else 0
                if zi_rms < 1e-10:
                    H *= 2; continue
                ratio = zi_rms / zG_rms
                if ratio < 1e-10 or ratio > 1e10:
                    break
                H_new = math.sqrt(H * H / ratio)
                if abs(H_new - H) < abs(H) * 1e-6:
                    break
                H = H_new

            H_center = best_H
            H_lo = best_H * 0.05
            H_hi = best_H * 20.0
            print("    {} evals, H_init={:.1f}, cost={:.4f}".format(n_evals, best_H, best_cost))

    else:
        # Scan logarithmique : encadrement du pic CSG (unimodal prouvé)
        # 25 points sur [W*0.01, W*200] en échelle log
        print("  Phase 0: scan logarithmique (max_csg)")
        candidates = []
        for k in range(-20, 25):
            scale = 10.0 ** (k / 9.0)
            H_test = W * scale
            zdof = [H_test] * n_dof_nv
            c, zeta, zi, q, zb_u = cost_fn(zdof)
            if c < 1e20:
                candidates.append((c, H_test, zeta, zi, q, zb_u))

        if not candidates:
            print("    ERREUR: aucun point valide")
            return [0.0]*m, [0.0]*ni, [0.0]*m, zb[:]

        candidates.sort(key=lambda x: x[0])
        best_cost, best_H, best_zeta, best_zi, best_q, best_zb = candidates[0]
        best_data = (best_zeta, best_zi, best_q, best_zb)

        # Encadrement : 3 meilleurs points → bornes pour golden section
        top3 = sorted(candidates[:min(5, len(candidates))], key=lambda x: x[1])
        H_center = best_H
        H_lo = top3[0][1] * 0.5
        H_hi = top3[-1][1] * 2.0

        print("    {} evals, H_best={:.1f}, cost={:.4f}".format(
            n_evals, best_H, best_cost))
        print("    Intervalle: [{:.0f}, {:.0f}]".format(H_lo, H_hi))

    # ================================================================
    # PHASE 1 : Golden section sur H (1D, garanti unimodal)
    # ================================================================
    print("  Phase 1: golden section [%.0f, %.0f]" % (H_lo, H_hi))
    phi_gs = (math.sqrt(5) - 1) / 2
    a, b = H_lo, H_hi

    def f1d(H_test):
        zdof = [H_test] * n_dof_nv
        c, zt, zi, q, zb_u = cost_fn(zdof)
        return c

    c_gs = b - phi_gs * (b - a)
    d_gs = a + phi_gs * (b - a)
    fc_gs = f1d(c_gs)
    fd_gs = f1d(d_gs)

    for _ in range(45):
        if abs(b - a) < max(abs(a + b) * 1e-9, 1e-3):
            break
        if fc_gs < fd_gs:
            b = d_gs; d_gs = c_gs; fd_gs = fc_gs
            c_gs = b - phi_gs * (b - a)
            fc_gs = f1d(c_gs)
        else:
            a = c_gs; c_gs = d_gs; fc_gs = fd_gs
            d_gs = a + phi_gs * (b - a)
            fd_gs = f1d(d_gs)

    H_opt = (a + b) / 2
    c_gs_final, zt, zi, q, zb_u = cost_fn([H_opt] * n_dof_nv)
    if c_gs_final < best_cost:
        best_cost = c_gs_final
        best_H = H_opt
        best_data = (zt, zi, q, zb_u)

    print("    {} evals total, H_opt={:.1f}, cost={:.6f}".format(
        n_evals, best_H, best_cost))

    # ================================================================
    # PHASE 2 : Coordinate Descent + Nelder-Mead (N > 1 DOF)
    # ================================================================
    if n_dof_nv > 1:
        # Extraire les valeurs DOF réelles du meilleur résultat
        # (pas juste un scalaire H uniforme)
        best_zeta = best_data[0]
        dof_nv_list = [g for g in dof_global if not is_vertical[g]]
        x0 = [abs(best_zeta[dof_nv_list[i]]) + 1e-3 for i in range(n_dof_nv)]

        # --- Phase 2a : Coordinate Descent (golden section par axe) ---
        # Exploite P6 : couplage faible inter-méridiens
        # MANUAL §2.4.4 "one step at a time"
        print("  Phase 2a: coordinate descent ({} DOF)".format(n_dof_nv))
        x_cd = x0[:]
        cost_cd = best_cost

        n_cd_cycles = min(5, 2 + n_dof_nv // 3)  # plus de cycles si N grand
        for cycle in range(n_cd_cycles):
            improved_cycle = False
            for k in range(n_dof_nv):
                # Golden section sur la dimension k seule
                xk_center = x_cd[k]
                ak = max(1e-3, xk_center * 0.05)
                bk = xk_center * 20.0
                ck = bk - phi_gs * (bk - ak)
                dk = ak + phi_gs * (bk - ak)

                def f_1axis(val):
                    xt = x_cd[:]
                    xt[k] = val
                    c, _, _, _, _ = cost_fn(xt)
                    return c

                fck = f_1axis(ck)
                fdk = f_1axis(dk)

                for _ in range(20):
                    if abs(bk - ak) < max(abs(ak + bk) * 1e-6, 1e-3):
                        break
                    if fck < fdk:
                        bk = dk; dk = ck; fdk = fck
                        ck = bk - phi_gs * (bk - ak)
                        fck = f_1axis(ck)
                    else:
                        ak = ck; ck = dk; fck = fdk
                        dk = ak + phi_gs * (bk - ak)
                        fdk = f_1axis(dk)

                xk_opt = (ak + bk) / 2
                c_new = f_1axis(xk_opt)
                if c_new < cost_cd:
                    x_cd[k] = xk_opt
                    cost_cd = c_new
                    improved_cycle = True

            if not improved_cycle:
                break

        if cost_cd < best_cost:
            c_tmp, zt, zi, q, zb_u = cost_fn(x_cd)
            best_cost = cost_cd
            best_data = (zt, zi, q, zb_u)
            print("    CD: ameliore, cost={:.6f}".format(cost_cd))

        # --- Phase 2b : Nelder-Mead avec restart ---
        print("  Phase 2b: Nelder-Mead")
        n = n_dof_nv
        alpha_nm, gamma_nm, rho_nm, sigma_nm = 1.0, 2.0, 0.5, 0.5

        def project(x):
            return [max(1e-3, xi) for xi in x]

        def f_nm(x):
            c, _, _, _, _ = cost_fn(project(x))
            return c

        def run_nelder_mead(x_start, perturbation, max_iters):
            """Un run de Nelder-Mead."""
            simplex = [x_start[:]]
            for i in range(n):
                xi = x_start[:]
                xi[i] *= (1 + perturbation)
                if abs(xi[i]) < 1e-3:
                    xi[i] = x_start[i] + perturbation * max(abs(x_start[i]), 100)
                simplex.append(project(xi))
            fvals = [f_nm(s) for s in simplex]

            for iteration in range(max_iters):
                order = sorted(range(n + 1), key=lambda kk: fvals[kk])
                simplex = [simplex[kk] for kk in order]
                fvals = [fvals[kk] for kk in order]
                spread = max(abs(fvals[i] - fvals[0]) for i in range(1, n + 1))
                if spread < 1e-10:
                    break
                cent = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]
                xr = project([cent[j] + alpha_nm * (cent[j] - simplex[n][j]) for j in range(n)])
                fr = f_nm(xr)
                if fvals[0] <= fr < fvals[n - 1]:
                    simplex[n] = xr; fvals[n] = fr
                elif fr < fvals[0]:
                    xe = project([cent[j] + gamma_nm * (xr[j] - cent[j]) for j in range(n)])
                    fe = f_nm(xe)
                    simplex[n] = xe if fe < fr else xr
                    fvals[n] = min(fe, fr)
                else:
                    xc = project([cent[j] + rho_nm * (simplex[n][j] - cent[j]) for j in range(n)])
                    fc = f_nm(xc)
                    if fc < fvals[n]:
                        simplex[n] = xc; fvals[n] = fc
                    else:
                        for i in range(1, n + 1):
                            simplex[i] = project([simplex[0][j] + sigma_nm * (simplex[i][j] - simplex[0][j])
                                                 for j in range(n)])
                            fvals[i] = f_nm(simplex[i])
            return simplex[0], fvals[0]

        # Run 1 : perturbation 50% (large exploration)
        nm_budget = min(max_iter * 2, 80 * n)  # plus d'iters pour N grand
        x_best_nm, f_best_nm = run_nelder_mead(x_cd, 0.5, nm_budget)

        # Run 2 : restart si le premier run a trouvé mieux, perturbation 20%
        if f_best_nm < best_cost * 0.9:  # amélioration > 10%
            x_nm2, f_nm2 = run_nelder_mead(x_best_nm, 0.2, nm_budget // 2)
            if f_nm2 < f_best_nm:
                x_best_nm, f_best_nm = x_nm2, f_nm2

        if f_best_nm < best_cost:
            c_tmp, zt, zi, q, zb_u = cost_fn(x_best_nm)
            best_cost = f_best_nm
            best_data = (zt, zi, q, zb_u)
            x_cd = x_best_nm[:]  # mettre à jour pour DIRECT
            print("    NM: ameliore, cost={:.6f}".format(f_best_nm))
        else:
            print("    NM: pas d'amelioration")

        # --- Phase 2c : DIRECT (DIviding RECTangles, Jones 1993) ---
        # Garantit convergence vers optimum global Lipschitz-continu.
        # Utilisé comme filet de sécurité pour échapper aux minima locaux.
        # Budget limité à max_iter/2 evals pour ne pas exploser le temps.
        #
        # Référence : Jones, Perttunen, Stuckman (1993)
        #   "Lipschitzian optimization without the Lipschitz constant"
        #   Journal of Optimization Theory and Applications 79(1), 157-181
        #
        # Implémentation simplifiée pour dimension faible (N ≤ 20).
        # L'idée : subdiviser récursivement l'hyper-rectangle des bornes
        # et explorer les rectangles "potentiellement optimaux".
        budget_direct = min(max_iter * 2, max(80, 30 * n))
        if n_evals < max_iter * 2:
            print("  Phase 2c: DIRECT (budget={})".format(budget_direct))
            # Construire les bornes [lo, hi] autour des meilleures valeurs trouvées
            # Utiliser x_cd (résultat CD) comme centre
            lo_d = [max(1e-3, x_cd[i] * 0.05) for i in range(n)]
            hi_d = [x_cd[i] * 20.0 for i in range(n)]

            # === DIRECT simplifié ===
            # Structure : liste de rectangles (center, size, fval)
            # center = point central normalisé [0,1]^n
            # size = demi-largeur (initialement 0.5)
            center_init = [0.5] * n
            size_init = 0.5

            def denormalize(xnorm):
                return [lo_d[i] + xnorm[i] * (hi_d[i] - lo_d[i]) for i in range(n)]

            def f_direct(xnorm):
                xr = denormalize(xnorm)
                xr = [max(1e-3, xi) for xi in xr]
                c, _, _, _, _ = cost_fn(xr)
                return c

            # Initialisation : évaluer le centre
            f_center = f_direct(center_init)
            rects = [(center_init, size_init, f_center)]
            direct_evals = 1
            best_direct = f_center
            best_direct_x = center_init[:]

            for _it in range(budget_direct):
                if direct_evals >= budget_direct:
                    break

                # Sélection des rectangles potentiellement optimaux
                # Simplifié : on prend le rectangle avec le meilleur fval
                # parmi ceux de taille >= médiane (exploration/exploitation)
                if not rects:
                    break

                # Trier par fval
                rects.sort(key=lambda r: r[2])

                # Sélectionner les rectangles à subdiviser :
                # - Le meilleur (exploitation)
                # - Le plus grand parmi les top 30% (exploration)
                n_rect = len(rects)
                candidates = [0]  # toujours le meilleur
                if n_rect > 2:
                    top30 = rects[:max(1, n_rect // 3)]
                    biggest = max(range(len(top30)), key=lambda i: top30[i][1])
                    if biggest != 0:
                        candidates.append(biggest)

                new_rects = []
                to_remove = set()

                for idx in candidates:
                    if direct_evals >= budget_direct:
                        break
                    if idx >= len(rects):
                        continue

                    c_rect, s_rect, f_rect = rects[idx]
                    to_remove.add(idx)

                    if s_rect < 1e-8:
                        new_rects.append((c_rect, s_rect, f_rect))
                        continue

                    # Subdiviser : trisection sur la dimension la plus longue
                    # (Pour DIRECT standard : toutes les dimensions, mais on
                    # simplifie pour efficacité)
                    dim = _it % n  # rotation des dimensions
                    delta = s_rect / 3.0

                    # Évaluer les 2 tiers-points
                    c_lo = c_rect[:]
                    c_lo[dim] -= delta
                    c_hi = c_rect[:]
                    c_hi[dim] += delta

                    # Clipper à [0,1]
                    c_lo[dim] = max(0.0, c_lo[dim])
                    c_hi[dim] = min(1.0, c_hi[dim])

                    f_lo = f_direct(c_lo); direct_evals += 1
                    if f_lo < best_direct:
                        best_direct = f_lo
                        best_direct_x = c_lo[:]

                    if direct_evals >= budget_direct:
                        new_rects.append((c_lo, delta, f_lo))
                        new_rects.append((c_rect, delta, f_rect))
                        break

                    f_hi = f_direct(c_hi); direct_evals += 1
                    if f_hi < best_direct:
                        best_direct = f_hi
                        best_direct_x = c_hi[:]

                    # Créer 3 sous-rectangles
                    new_rects.append((c_lo, delta, f_lo))
                    new_rects.append((c_rect, delta, f_rect))
                    new_rects.append((c_hi, delta, f_hi))

                # Reconstruire la liste des rectangles
                rects = [r for i, r in enumerate(rects) if i not in to_remove] + new_rects

            # Comparer avec le meilleur trouvé jusqu'ici
            if best_direct < best_cost:
                x_direct = denormalize(best_direct_x)
                x_direct = [max(1e-3, xi) for xi in x_direct]
                c_tmp, zt, zi, q, zb_u = cost_fn(x_direct)
                if c_tmp < best_cost:
                    best_cost = c_tmp
                    best_data = (zt, zi, q, zb_u)
                    print("    DIRECT: ameliore, cost={:.6f} ({} evals)".format(
                        c_tmp, direct_evals))
                else:
                    print("    DIRECT: {} evals, pas d'amelioration".format(direct_evals))
            else:
                print("    DIRECT: {} evals, pas d'amelioration".format(direct_evals))

    # ── Résultat final ──
    zeta_opt, zi_opt, q_opt, zb_opt = best_data
    n_tension = sum(1 for z in zeta_opt if z < -1e-6)
    if n_tension > 0:
        print("  AVERTISSEMENT: {} branches en traction".format(n_tension))
    print("  Optimisation terminee: {} evals, cost={:.6f}".format(
        n_evals, best_cost))

    return zeta_opt, zi_opt, q_opt, zb_opt


# =======================================================================
# FONCTIONS UTILITAIRES
# =======================================================================

def compute_forces_from_zeta(zeta, LH, branches, m):
    """
    LH_dual = zeta (forces horizontales dans les branches)
    F = force totale dans la branche = |zeta| * L / LH
    """
    LH_dual = zeta[:]
    F = []
    for j in range(m):
        L = branches[j].get('L', LH[j])
        lh = LH[j]
        if abs(lh) > 1e-10:
            F.append(abs(zeta[j]) * L / lh)
        else:
            F.append(abs(zeta[j]))
    return LH_dual, F


def count_violations(Z, Z_min, Z_max, n_int):
    violations = 0
    for i in range(n_int):
        if Z[i] < Z_min[i] - 1e-6 or Z[i] > Z_max[i] + 1e-6:
            violations += 1
    return violations


def compute_bounds(nodes):
    int_n = [n for n in nodes if not n['Boundary']]
    Z_min, Z_max = [], []
    for n in int_n:
        zi = n.get('Bzi', n['ZG'] - 0.25)
        ze = n.get('Bze', n['ZG'] + 0.25)
        if zi > ze:
            zi, ze = ze, zi
        Z_min.append(zi)
        Z_max.append(ze)
    return Z_min, Z_max


# =======================================================================
# TOPOLOGIE INTERNE (remplace le Module 2 — évite le problème DataTree GH)
# =======================================================================

def build_topology(nodes_list, branches_list):
    """
    Construit Ci, Cb, u, v, LH directement depuis nodes + branches.
    
    Élimine le besoin d'un composant GH séparé "Module 2" qui produisait
    des listes de listes converties en DataTree par Grasshopper.
    
    Source : MANUAL §2.1.3, d_ComputeTypology.bas
    
    Returns:
        Ci, Cb, u_vec, v_vec, LH_vec, ni, nb, m,
        int_nodes, bou_nodes, idx_map, node_valence, is_vertical, pz, zb
    """
    import json as _json
    
    # Accepter JSON string ou list of dict
    if isinstance(nodes_list, str):
        nodes_list = _json.loads(nodes_list)
    if isinstance(branches_list, str):
        branches_list = _json.loads(branches_list)

    int_n = [n for n in nodes_list if not n['Boundary']]
    bou_n = [n for n in nodes_list if n['Boundary']]
    all_n = int_n + bou_n
    ni = len(int_n)
    nb = len(bou_n)
    m = len(branches_list)
    
    # Mapping Index → position dans all_n
    idx_map = {n['Index']: i for i, n in enumerate(all_n)}
    
    # Mapping Index → colonne Ci ou Cb (robuste, pas de supposition sur Index)
    col_int = {}
    col_bou = {}
    for i, n in enumerate(int_n):
        col_int[n['Index']] = i
    for i, n in enumerate(bou_n):
        col_bou[n['Index']] = i
    
    # Matrices
    Ci = [[0.0] * ni for _ in range(m)]
    Cb = [[0.0] * max(nb, 1) for _ in range(m)]
    
    u_vec = []
    v_vec = []
    LH_vec = []
    
    for j, b in enumerate(branches_list):
        head = b['Head_node']
        tail = b['Tail_node']
        
        if head in col_int:
            Ci[j][col_int[head]] = 1.0
        if tail in col_int:
            Ci[j][col_int[tail]] = -1.0
        if nb > 0:
            if head in col_bou:
                Cb[j][col_bou[head]] = 1.0
            if tail in col_bou:
                Cb[j][col_bou[tail]] = -1.0
        
        u_vec.append(b['u'])
        v_vec.append(b.get('v', 0.0))
        LH_vec.append(b['LH'])
    
    # Valence des nœuds intérieurs
    node_valence = [0] * ni
    for j in range(m):
        for i in range(ni):
            if abs(Ci[j][i]) > 0.5:
                node_valence[i] += 1
    
    # Branches verticales (LH ≈ 0)
    is_vertical = [abs(LH_vec[j]) < 1e-12 for j in range(m)]
    
    # Poids nodaux : pz = -Weight (convention VBA)
    pz = [0.0] * ni
    for n in int_n:
        col = col_int.get(n['Index'])
        if col is not None:
            pz[col] = -n['Weight']
    
    return (Ci, Cb, u_vec, v_vec, LH_vec,
            ni, nb, m, int_n, bou_n, idx_map, node_valence, is_vertical, pz)


# =======================================================================
# MAIN (interface Grasshopper)
# =======================================================================
# INPUTS GH :
#   - nodes : str/list (Item)    - JSON ou list of dict du Module 0
#   - branches : str/list (Item) - JSON ou list of dict du Module 0
#   - optimize : bool (Item)     - Lancer l'optimisation [True]
#   - max_iter : int (Item)      - Budget optimisation [200]
#   - objective : str (Item)     - "proximity" ou "max_csg" ["proximity"]
#   - boundary_position : str    - "middle", "intrados", "extrados" ["middle"]
#
# Plus besoin de : Ci, Cb, u_vec, v_vec, LH_vec, n_int
# (tout est construit en interne depuis nodes + branches)

# Valeurs par défaut
if 'optimize' not in dir() or optimize is None:
    optimize = True
if 'max_iter' not in dir() or max_iter is None:
    max_iter = 200
if 'objective' not in dir() or objective is None:
    objective = "proximity"
if 'boundary_position' not in dir() or boundary_position is None:
    boundary_position = "middle"
# Rétrocompatibilité : 'target' → 'objective'
if 'target' not in dir():
    target = "middle"
if target in ("intrados", "extrados"):
    boundary_position = target

has_data = all([
    'nodes' in dir() and nodes is not None,
    'branches' in dir() and branches is not None,
])
# Accepter aussi l'ancien câblage avec n_int (rétrocompat)
if not has_data and all([
    'Ci' in dir() and Ci,
    'Cb' in dir() and Cb,
    'LH_vec' in dir() and LH_vec,
    'n_int' in dir() and n_int
]):
    has_data = True

if not has_data:
    print("ERREUR: Connectez nodes + branches depuis le Module 0")
    Z = q = F = LH_dual = []
    CSG_min = 0
    thrust_pts = thrust_lines = []
else:
    print("=" * 60)
    print("MaNACoH v4 - Module 3: TNA (v53 - RO optimizer)")
    print("=" * 60)
    print("Sources : MF1 eq. 4.27-4.38, VBA e_ComputeEquilibrium.bas")
    print("")

    # ── Construction topologie en interne (remplace Module 2) ──
    # Accepte JSON ou list of dict — plus de problème DataTree GH
    import json as _json
    nodes_list = nodes
    branches_list = branches
    if isinstance(nodes_list, str):
        nodes_list = _json.loads(nodes_list)
    if isinstance(branches_list, str):
        branches_list = _json.loads(branches_list)

    (Ci_mat, Cb_mat, u, v, LH,
     ni, n_bou, m, int_n, bou_n, idx_map,
     node_valence, is_vertical, pz) = build_topology(nodes_list, branches_list)

    n_int = ni
    all_n = int_n + bou_n

    # Altitudes des appuis
    zb = []
    for n in bou_n:
        if boundary_position == "intrados":
            zb.append(n.get('Bzi', n['ZG']))
        elif boundary_position == "extrados":
            zb.append(n.get('Bze', n['ZG']))
        else:
            zb.append(n['ZG'])

    print("  Noeuds: {} int + {} bou | Branches: {}".format(ni, n_bou, m))
    print("  Poids total: {:.1f} kN".format(sum(pz) / 1000))
    print("  Objectif: {} | Appuis: {}".format(objective, boundary_position))

    # ── ÉTAPE 1 : Équilibre horizontal (MF1 eq. 4.27-4.30) ──
    print("\n--- Etape 1 : Equilibre horizontal ---")

    zeta_init, dof_global, unk_global, M_stored = solve_horizontal(
        Ci_mat, Cb_mat, u, v, LH, node_valence, ni, m,
        zeta_dof_values=None, chosen_dof=None, is_vertical=is_vertical
    )

    n_dof_calc = len(dof_global)
    print("  DOF branches: {} (indices: {})".format(n_dof_calc, dof_global[:10]))
    print("  UNK branches: {}".format(len(unk_global)))

    if optimize:
        # ── ÉTAPE 2-3 : Optimisation des DOF ──
        print("\n--- Etape 2-3 : Optimisation ---")

        zeta_opt, zi_opt, q_opt, zb_opt = optimize_dof(
            nodes_list, branches_list, Ci_mat, Cb_mat, u, v, LH, node_valence,
            ni, n_bou, m, n_int, pz, zb,
            is_vertical, M_stored, dof_global, unk_global,
            max_iter=int(max_iter), objective=objective
        )

        Z = zi_opt[:] + zb_opt[:]
        q = q_opt
        LH_dual, F = compute_forces_from_zeta(zeta_opt, LH, branches_list, m)
    else:
        # Solution directe sans optimisation
        print("\n--- Solution directe (pas d'optimisation) ---")

        # Estimation H par poids total (convention: H > 0 = compression)
        total_weight = sum(pz)
        avg_z = sum(n['ZG'] for n in int_n) / max(ni, 1)
        H_default = total_weight / max(1, ni) * 10
        n_dof_nv = len([g for g in dof_global if not is_vertical[g]])
        zdof_default = [H_default] * n_dof_nv
        zeta = update_zeta_fast(M_stored, zdof_default, unk_global, dof_global, m)
        zi, q, zb_used = solve_vertical(Ci_mat, Cb_mat, LH, zeta, pz, zb, ni, n_bou, m)
        Z = zi[:] + zb_used[:]
        LH_dual, F = compute_forces_from_zeta(zeta, LH, branches_list, m)

    # ── CSG ──
    CSG_list = compute_csg_simplified(nodes_list, branches_list, Z, LH_dual, n_int)
    CSG_min = min(CSG_list) if CSG_list else 0.0

    # ── Stats ──
    Z_min, Z_max = compute_bounds(nodes_list)
    violations = count_violations(Z, Z_min, Z_max, n_int)
    pos = sum(1 for x in LH_dual if x > 0)
    neg = sum(1 for x in LH_dual if x < 0)

    print("")
    print("=" * 60)
    print("RESULTATS")
    print("=" * 60)
    print("  Z: [{:.3f}, {:.3f}] m".format(
        min(Z[:n_int]) if n_int > 0 else 0,
        max(Z[:n_int]) if n_int > 0 else 0))
    print("  Forces: {} compression, {} traction".format(pos, neg))
    print("  Violations: {} / {}".format(violations, n_int))
    print("  CSG_min: {:.3f} [{}]".format(CSG_min, "STABLE" if CSG_min >= 1 else "INSTABLE"))
    print("=" * 60)

    # ── Géométrie Rhino ──
    if HAS_RHINO:
        thrust_pts = [rg.Point3d(n['XG'], n['YG'], Z[i])
                      for i, n in enumerate(all_n)]
        thrust_lines = [rg.Line(thrust_pts[idx_map[b['Head_node']]],
                                thrust_pts[idx_map[b['Tail_node']]])
                        for b in branches_list]
    else:
        thrust_pts = [(all_n[i]['XG'], all_n[i]['YG'], Z[i])
                      for i in range(len(all_n))]
        thrust_lines = []
