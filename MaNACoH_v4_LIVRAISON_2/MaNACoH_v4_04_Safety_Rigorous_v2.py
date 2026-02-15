#! python3
"""
MaNACoH v4 - Module 4 : SAFETY (Version Rigoureuse)
=====================================================
Calcul rigoureux des coefficients de sécurité selon Fantin (2017) et Delbecq (1983)

RÉFÉRENCES:
  - Fantin (2017) "Méthode des réseaux de forces" - Thèse ENPC
  - Delbecq (1983) "Analyse de la stabilité des ponts en maçonnerie" - Thèse
  - SETRA (1982) "Ponts en maçonnerie - Guide VOUTE"

=============================================================================
DÉFINITIONS MATHÉMATIQUES (repère local du joint x'y'z')
=============================================================================

1. EFFORTS SUR LE JOINT (repère local)
   - N' : Effort normal (compression positive)
   - M' : Moment fléchissant
   - V' : Effort tranchant

2. EXCENTRICITÉ
   e = M' / N'
   
   Interprétation: distance entre le centre de pression et le centre du joint

3. CSG - COEFFICIENT DE SÉCURITÉ GÉOMÉTRIQUE (Fantin eq. 3.33)
   
   CSG_j = (h'/2) / |e| = h' / |2M'/N'| = N' · h' / (2|M'|)
   
   où h' = hauteur du joint (longueur dans la direction de la normale)
   
   - CSG < 1 : Centre de pression HORS du joint (instable)
   - CSG = 1 : Centre de pression AU BORD (limite)
   - CSG > 1 : Centre de pression DANS le joint (stable)
   
   CSG de la structure: CSG_J = min_j(CSG_j)

4. CSS - COEFFICIENT DE SÉCURITÉ STATIQUE (Fantin eq. 3.22, Delbecq)
   
   Prend en compte la résistance finie à la compression σ_c
   
   Valeurs adimensionnées:
     n' = N' / N_c  où N_c = b' · h' · σ_c
     m' = (2/h') · M' / N_c
   
   Critère de résistance: |m'| ≤ n'(1 - n')
   
   CSS_j = (1 - |m'/n'|) / n' = (1 - 2|e|/h') / n'
   
   - CSS > 1 : Contrainte max < σ_c (admissible)
   - CSS = 1 : Contrainte max = σ_c (limite)
   - CSS < 1 : Contrainte max > σ_c (rupture compression)

5. CONTRAINTE MAXIMALE SUR LE JOINT
   
   Pour une distribution linéaire (triangulaire si e > h'/6):
   
   σ_max = 2N' / (b' · h'_comprimé)
   
   où h'_comprimé = h' - 2|e| (hauteur comprimée si pas de traction)
   
   Si e > h'/2 : le joint s'ouvre, pas d'équilibre possible

6. CRITÈRE DE GLISSEMENT (Coulomb)
   
   |V'| ≤ μ · N'  où μ = tan(φ) ≈ 0.5 à 0.7
   
   CS_glissement = μ · N' / |V'|

=============================================================================
INPUTS (Type Hints Grasshopper):
  - nodes_data : ghdoc (List) - Du Module 2
  - branches_data : ghdoc (List) - Du Module 2
  - Z : ghdoc (List) - Du Module 3
  - LH_dual : ghdoc (List) - Du Module 3 (forces horizontales)
  - sigma_c : float (Item) - Résistance compression [Pa] [10e6 = 10 MPa]
  - mu : float (Item) - Coefficient frottement [0.5]

OUTPUTS DONNÉES:
  - CSG : list - CSG par branche (coefficient géométrique)
  - CSS : list - CSS par branche (coefficient statique)
  - CS_glissement : list - Coefficient sécurité glissement par branche
  - CSG_min : float - CSG minimum (structure)
  - CSS_min : float - CSS minimum (structure)
  
OUTPUTS EFFORTS:
  - N_list : list - Effort normal par branche [N]
  - M_list : list - Moment par branche [N·m]
  - V_list : list - Effort tranchant par branche [N]
  - e_list : list - Excentricité par branche [m]
  - sigma_max_list : list - Contrainte max par branche [Pa]

OUTPUTS VISUALISATION:
  - branch_lines : Line list - Lignes thrust network
  - branch_colors_CSG : Color list - Couleurs selon CSG
  - branch_colors_CSS : Color list - Couleurs selon CSS
  - joint_lines : Line list - Joints [Pi, Pe]
  - joint_colors : Color list - Couleurs joints selon CSG
  - thrust_pts : Point3d list - Points de passage thrust dans joint
"""

import Rhino.Geometry as rg
import System.Drawing as sd
import math

# === VALEURS PAR DÉFAUT ===
if 'sigma_c' not in dir() or sigma_c is None:
    sigma_c = 10e6  # 10 MPa - résistance compression maçonnerie
if 'mu' not in dir() or mu is None:
    mu = 0.5  # tan(27°) - coefficient frottement Coulomb

# === FONCTIONS COULEUR ===

def value_to_color(value, v_min, v_critical, v_good, v_max):
    """
    Gradient de couleur:
    - Rouge foncé : value < v_min
    - Rouge : v_min <= value < v_critical  
    - Orange : value ≈ v_critical
    - Jaune : v_critical < value < v_good
    - Vert : value >= v_good
    """
    if value < v_min:
        return sd.Color.FromArgb(255, 139, 0, 0)      # Rouge foncé (très mauvais)
    elif value < v_critical:
        # Rouge -> Orange
        t = (value - v_min) / (v_critical - v_min) if v_critical > v_min else 0
        r = 255
        g = int(128 * t)
        return sd.Color.FromArgb(255, r, g, 0)
    elif value < v_good:
        # Orange -> Jaune -> Vert clair
        t = (value - v_critical) / (v_good - v_critical) if v_good > v_critical else 0
        r = 255 - int(127 * t)
        g = 128 + int(127 * t)
        return sd.Color.FromArgb(255, r, g, 0)
    elif value < v_max:
        # Vert clair -> Vert
        t = (value - v_good) / (v_max - v_good) if v_max > v_good else 1
        r = int(128 * (1 - t))
        g = 255
        return sd.Color.FromArgb(255, r, g, 0)
    else:
        return sd.Color.FromArgb(255, 0, 255, 0)      # Vert (très bon)


def csg_to_color(csg):
    """Couleur selon CSG: critique à 1.0, bon à 2.0"""
    return value_to_color(csg, 0.5, 1.0, 2.0, 3.0)


def css_to_color(css):
    """Couleur selon CSS: critique à 1.0, bon à 2.0"""
    return value_to_color(css, 0.0, 1.0, 2.0, 4.0)


# === FONCTIONS DE CALCUL DES EFFORTS ===

def compute_joint_geometry(b, nodes, idx_map, all_n):
    """
    Calcule la géométrie du joint pour une branche.
    
    Retourne:
      - h : hauteur du joint (longueur du segment [Pi, Pe])
      - b_width : largeur du joint (estimée)
      - joint_center : centre du joint
      - joint_dir : direction du joint (normalisée)
    """
    if 'xi' not in b:
        return None, None, None, None
    
    Pi = (b['xi'], b['yi'], b['zi'])
    Pe = (b['xe'], b['ye'], b['ze'])
    
    # Hauteur du joint
    h = math.sqrt((Pe[0]-Pi[0])**2 + (Pe[1]-Pi[1])**2 + (Pe[2]-Pi[2])**2)
    
    if h < 1e-10:
        return None, None, None, None
    
    # Centre du joint
    joint_center = ((Pi[0]+Pe[0])/2, (Pi[1]+Pe[1])/2, (Pi[2]+Pe[2])/2)
    
    # Direction du joint (normalisée)
    joint_dir = ((Pe[0]-Pi[0])/h, (Pe[1]-Pi[1])/h, (Pe[2]-Pi[2])/h)
    
    # Largeur du joint (estimée à partir de Joint_area si disponible)
    if 'Joint_area' in b and h > 0:
        b_width = b['Joint_area'] / h
    else:
        b_width = b.get('L', 1.0)  # Fallback: longueur de la branche
    
    return h, b_width, joint_center, joint_dir


def compute_thrust_passage(b, nodes, Z, idx_map, all_n):
    """
    Calcule le point de passage du thrust network dans le joint.
    
    Le thrust passe au milieu de la branche, à la hauteur Z moyenne.
    """
    i1 = idx_map[b['Head_node']]
    i2 = idx_map[b['Tail_node']]
    n1, n2 = all_n[i1], all_n[i2]
    
    # Centre de la branche du thrust network
    cx = (n1['XG'] + n2['XG']) / 2
    cy = (n1['YG'] + n2['YG']) / 2
    cz = (Z[i1] + Z[i2]) / 2
    
    return (cx, cy, cz)


def compute_joint_efforts(b, nodes, Z, LH_dual, idx_map, all_n, nb, branch_pos=0):
    """
    Calcule les efforts (N', M', V') sur le joint dans le repère local.
    
    MÉTHODE:
    1. La force dans la branche est dirigée selon la branche du thrust network
    2. N' = composante normale au joint (perpendiculaire à [Pi,Pe])
    3. V' = composante tangentielle au joint (parallèle à [Pi,Pe])
    4. M' = N' × excentricité (moment par rapport au centre du joint)
    
    Sources: Fantin (2017) eq. 3.33, SETRA (1982) p.45
    """
    # Géométrie du joint
    h, b_width, joint_center, joint_dir = compute_joint_geometry(b, nodes, idx_map, all_n)
    
    if h is None:
        return None
    
    # Point de passage du thrust
    thrust_pt = compute_thrust_passage(b, nodes, Z, idx_map, all_n)
    
    # Indices des noeuds
    i1 = idx_map[b['Head_node']]
    i2 = idx_map[b['Tail_node']]
    n1, n2 = all_n[i1], all_n[i2]
    
    # Direction de la force (le long de la branche du thrust network)
    dx = n2['XG'] - n1['XG']
    dy = n2['YG'] - n1['YG']
    dz = Z[i2] - Z[i1]
    L_thrust = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    if L_thrust < 1e-10:
        return None
    
    force_dir = (dx/L_thrust, dy/L_thrust, dz/L_thrust)
    
    # Force horizontale dans la branche (du Module 3)
    # E5 fix : utiliser branch_pos (0-based), pas b['Index'] (1-based)
    if branch_pos < len(LH_dual):
        LH = abs(LH_dual[branch_pos])
    else:
        LH = 1.0
    
    # Force totale (F = LH / cos(angle_horizontal))
    # F · cos(θ) = LH où θ = angle avec l'horizontale
    LH_branch = b.get('LH', 1.0)
    L_branch = b.get('L', 1.0)
    
    if LH_branch > 1e-10:
        F = LH * L_branch / LH_branch  # F = LH * L / LH_geom
    else:
        F = LH
    
    # Vecteur force
    Fx = F * force_dir[0]
    Fy = F * force_dir[1]
    Fz = F * force_dir[2]
    
    # Normale au joint (perpendiculaire au segment [Pi,Pe])
    # On cherche la direction qui va du centre du joint vers le thrust
    bx = thrust_pt[0] - joint_center[0]
    by = thrust_pt[1] - joint_center[1]
    bz = thrust_pt[2] - joint_center[2]
    
    # Projection sur la direction du joint = excentricité signée
    e_signed = joint_dir[0]*bx + joint_dir[1]*by + joint_dir[2]*bz
    e = abs(e_signed)
    
    # Composante normale de la force (perpendiculaire au joint)
    # Le joint est orienté selon joint_dir, donc la normale est perpendiculaire
    # N' = projection de F sur la direction thrust_pt - joint_center
    dist_to_joint = math.sqrt(bx*bx + by*by + bz*bz)
    
    if dist_to_joint > 1e-10:
        normal_dir = (bx/dist_to_joint, by/dist_to_joint, bz/dist_to_joint)
    else:
        # Si le thrust passe exactement par le centre, utiliser une normale par défaut
        # (perpendiculaire à la direction du joint et horizontale si possible)
        normal_dir = (0, 0, 1)  # Approximation
    
    # Effort normal N' = F · normal_dir (compression positive)
    N_prime = abs(Fx*normal_dir[0] + Fy*normal_dir[1] + Fz*normal_dir[2])
    
    # Effort tranchant V' = F · joint_dir
    V_prime = Fx*joint_dir[0] + Fy*joint_dir[1] + Fz*joint_dir[2]
    
    # Moment M' = N' × e (par rapport au centre du joint)
    M_prime = N_prime * e_signed
    
    return {
        'N': N_prime,          # Effort normal [N]
        'M': M_prime,          # Moment [N·m]
        'V': V_prime,          # Effort tranchant [N]
        'e': e,                # Excentricité [m]
        'e_signed': e_signed,  # Excentricité signée [m]
        'h': h,                # Hauteur du joint [m]
        'b': b_width,          # Largeur du joint [m]
        'F': F,                # Force totale [N]
        'thrust_pt': thrust_pt,
        'joint_center': joint_center
    }


# === FONCTIONS DE CALCUL DES COEFFICIENTS DE SÉCURITÉ ===

def compute_CSG(N, M, h):
    """
    Calcule le CSG - Coefficient de Sécurité Géométrique
    
    CSG = (h/2) / |e| = (h/2) / |M/N| = N·h / (2|M|)
    
    Source: Fantin (2017) eq. 3.33
    """
    if N < 1e-10:
        return 0.0  # Pas d'effort normal = pas de stabilité
    
    e = abs(M / N) if N > 1e-10 else float('inf')
    
    if e < 1e-10:
        return 100.0  # Centre de pression parfaitement centré
    
    csg = (h / 2) / e
    
    return min(csg, 100.0)


def compute_CSS(N, M, h, b, sigma_c):
    """
    Calcule le CSS - Coefficient de Sécurité Statique
    
    CSS = (1 - |m'/n'|) / n' = (1 - 2|e|/h) / n'
    
    où n' = N / (b·h·σ_c)
    
    Source: Fantin (2017) eq. 3.22, Delbecq (1983)
    """
    if N < 1e-10 or h < 1e-10 or b < 1e-10:
        return 0.0
    
    # Effort normal max admissible
    N_c = b * h * sigma_c
    
    # Valeur adimensionnée
    n_prime = N / N_c
    
    if n_prime < 1e-10:
        return 100.0  # Très faible chargement
    
    # Excentricité relative
    e = abs(M / N) if N > 1e-10 else 0
    e_rel = 2 * e / h  # = |m'/n'|
    
    if e_rel >= 1.0:
        # Centre de pression hors du joint
        return 0.0
    
    css = (1 - e_rel) / n_prime
    
    return min(css, 100.0)


def compute_sigma_max(N, M, h, b):
    """
    Calcule la contrainte maximale sur le joint.
    
    Pour une distribution triangulaire (sans traction):
    σ_max = 2N / (b · h_comprimé)
    
    où h_comprimé = h - 2|e| = h(1 - 2|e|/h)
    
    Source: SETRA (1982) p.45
    """
    if N < 1e-10 or h < 1e-10 or b < 1e-10:
        return 0.0
    
    e = abs(M / N) if N > 1e-10 else 0
    
    # Hauteur comprimée
    h_comp = h - 2 * e
    
    if h_comp <= 0:
        # Joint complètement ouvert
        return float('inf')
    
    sigma_max = 2 * N / (b * h_comp)
    
    return sigma_max


def compute_CS_glissement(N, V, mu):
    """
    Calcule le coefficient de sécurité au glissement (Coulomb).
    
    CS = μ·N / |V|
    
    Source: SETRA (1982) p.35
    """
    if abs(V) < 1e-10:
        return 100.0  # Pas d'effort tranchant
    
    if N < 1e-10:
        return 0.0  # Pas d'effort normal = glissement assuré
    
    cs = mu * N / abs(V)
    
    return min(cs, 100.0)


# === MAIN ===

# Accepter soit nodes/branches soit nodes_data/branches_data
if 'nodes_data' in dir() and nodes_data and len(nodes_data) > 0:
    nodes = nodes_data
if 'branches_data' in dir() and branches_data and len(branches_data) > 0:
    branches = branches_data

# Vérifier les données
has_data = all([
    'nodes' in dir() and nodes and len(nodes) > 0,
    'branches' in dir() and branches and len(branches) > 0,
    'Z' in dir() and Z and len(Z) > 0
])

# LH_dual optionnel
if 'LH_dual' not in dir() or LH_dual is None:
    LH_dual = [1.0] * len(branches) if 'branches' in dir() else []

if not has_data:
    print("ERREUR: Connectez nodes_data, branches_data (Module 2), Z (Module 3)")
    CSG = CSS = CS_glissement = []
    CSG_min = CSS_min = 0
    N_list = M_list = V_list = e_list = sigma_max_list = []
    branch_lines = branch_colors_CSG = branch_colors_CSS = []
    joint_lines = joint_colors = thrust_pts = []
else:
    print("="*60)
    print("MaNACoH v4 - Module 4: SAFETY (Rigoureux)")
    print("="*60)
    print("Références: Fantin (2017), Delbecq (1983), SETRA (1982)")
    print("")
    
    # Préparer données
    int_n = [n for n in nodes if not n['Boundary']]
    bou_n = [n for n in nodes if n['Boundary']]
    all_n = int_n + bou_n
    idx_map = {n['Index']: i for i, n in enumerate(all_n)}
    nb = len(branches)
    
    print("Parametres materiaux:")
    print("  sigma_c = {:.1f} MPa (resistance compression)".format(sigma_c / 1e6))
    print("  mu = {:.2f} (coefficient frottement)".format(mu))
    print("")
    
    # === CALCUL DES EFFORTS ET COEFFICIENTS ===
    
    CSG = []
    CSS = []
    CS_glissement = []
    N_list = []
    M_list = []
    V_list = []
    e_list = []
    sigma_max_list = []
    thrust_pts = []
    
    for b_pos, b in enumerate(branches):
        efforts = compute_joint_efforts(b, nodes, Z, LH_dual, idx_map, all_n, nb, branch_pos=b_pos)
        
        if efforts is None:
            CSG.append(1.0)
            CSS.append(1.0)
            CS_glissement.append(1.0)
            N_list.append(0.0)
            M_list.append(0.0)
            V_list.append(0.0)
            e_list.append(0.0)
            sigma_max_list.append(0.0)
            thrust_pts.append(rg.Point3d(0, 0, 0))
            continue
        
        N = efforts['N']
        M = efforts['M']
        V = efforts['V']
        e = efforts['e']
        h = efforts['h']
        b_width = efforts['b']
        
        # Calcul des coefficients
        csg = compute_CSG(N, M, h)
        css = compute_CSS(N, M, h, b_width, sigma_c)
        cs_g = compute_CS_glissement(N, V, mu)
        sigma_max = compute_sigma_max(N, M, h, b_width)
        
        CSG.append(csg)
        CSS.append(css)
        CS_glissement.append(cs_g)
        N_list.append(N)
        M_list.append(M)
        V_list.append(V)
        e_list.append(e)
        sigma_max_list.append(sigma_max)
        
        # Point de passage du thrust
        tp = efforts['thrust_pt']
        thrust_pts.append(rg.Point3d(tp[0], tp[1], tp[2]))
    
    # Valeurs globales
    CSG_min = min(CSG) if CSG else 0
    CSS_min = min(CSS) if CSS else 0
    CS_glissement_min = min(CS_glissement) if CS_glissement else 0
    
    # === VISUALISATION ===
    
    branch_lines = []
    branch_colors_CSG = []
    branch_colors_CSS = []
    
    for i, b in enumerate(branches):
        i1, i2 = idx_map[b['Head_node']], idx_map[b['Tail_node']]
        pt1 = rg.Point3d(all_n[i1]['XG'], all_n[i1]['YG'], Z[i1])
        pt2 = rg.Point3d(all_n[i2]['XG'], all_n[i2]['YG'], Z[i2])
        branch_lines.append(rg.Line(pt1, pt2))
        branch_colors_CSG.append(csg_to_color(CSG[i]))
        branch_colors_CSS.append(css_to_color(CSS[i]))
    
    joint_lines = []
    joint_colors = []
    
    for i, b in enumerate(branches):
        if 'xi' in b:
            Pi = rg.Point3d(b['xi'], b['yi'], b['zi'])
            Pe = rg.Point3d(b['xe'], b['ye'], b['ze'])
            joint_lines.append(rg.Line(Pi, Pe))
            joint_colors.append(csg_to_color(CSG[i]))
    
    # === STATISTIQUES ===
    
    print("="*60)
    print("RESULTATS")
    print("="*60)
    
    # CSG
    csg_stable = sum(1 for c in CSG if c >= 1.0)
    print("\nCSG (Coefficient Securite Geometrique):")
    print("  CSG_min = {:.3f}".format(CSG_min))
    print("  CSG_max = {:.3f}".format(max(CSG) if CSG else 0))
    print("  Branches stables (CSG>=1): {}/{}".format(csg_stable, nb))
    print("  Status: {}".format("STABLE" if CSG_min >= 1.0 else "INSTABLE"))
    
    # CSS
    css_stable = sum(1 for c in CSS if c >= 1.0)
    print("\nCSS (Coefficient Securite Statique):")
    print("  CSS_min = {:.3f}".format(CSS_min))
    print("  CSS_max = {:.3f}".format(max(CSS) if CSS else 0))
    print("  Branches admissibles (CSS>=1): {}/{}".format(css_stable, nb))
    print("  Status: {}".format("ADMISSIBLE" if CSS_min >= 1.0 else "RUPTURE COMPRESSION"))
    
    # Glissement
    cs_g_stable = sum(1 for c in CS_glissement if c >= 1.0)
    print("\nGlissement (Coulomb):")
    print("  CS_min = {:.3f}".format(CS_glissement_min))
    print("  Branches stables: {}/{}".format(cs_g_stable, nb))
    print("  Status: {}".format("PAS DE GLISSEMENT" if CS_glissement_min >= 1.0 else "RISQUE GLISSEMENT"))
    
    # Contraintes
    sigma_max_global = max(sigma_max_list) if sigma_max_list else 0
    print("\nContraintes:")
    print("  sigma_max = {:.2f} MPa".format(sigma_max_global / 1e6))
    print("  sigma_c   = {:.2f} MPa".format(sigma_c / 1e6))
    print("  Ratio     = {:.2f}".format(sigma_max_global / sigma_c if sigma_c > 0 else 0))
    
    # Excentricités
    e_max = max(e_list) if e_list else 0
    print("\nExcentricites:")
    print("  e_max = {:.4f} m".format(e_max))
    
    print("\n" + "="*60)
    print("LEGENDE COULEURS:")
    print("  Rouge fonce : CSG/CSS < 0.5 (tres critique)")
    print("  Rouge       : 0.5 <= CSG/CSS < 1.0 (instable/inadmissible)")
    print("  Orange      : CSG/CSS ≈ 1.0 (limite)")
    print("  Jaune       : 1.0 < CSG/CSS < 2.0 (acceptable)")
    print("  Vert        : CSG/CSS >= 2.0 (bon)")
    print("="*60)
