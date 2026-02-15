#!/usr/bin/env python3
"""
MaNACoH v4 — Module 1 : GEOMETRIE PARAMETRIQUE
================================================
Génère les géométries paramétriques pour la Thrust Network Analysis.

FIDÈLE AU VBA :
  - c1_CreateModelParameters.bas : paramètres pré-définis
  - c1_CreateModelVinAbe.bas : a3_CreateRegularGeo (création noeuds/branches)
  - c1_CreateModel.bas : a3_CreateRegularGeo wrapper

GÉOMÉTRIES SUPPORTÉES :
  - Arc cylindrique (Cyl) : phi ∈ [phi0, phi1], y ∈ [y0, y1]
  - Dôme sphérique (Sphe) : phi ∈ [phi0, phi1], theta ∈ [theta0, theta1]

CAS PRÉDÉFINIS (MANUAL Chapter 4) :
  "heyman_15"     : Fig 4.27, 15 voussoirs, rm=10, t/R=0.15
  "heyman_tmin"   : Fig 4.28, 15 voussoirs, rm=10, t/R=0.1075
  "oikono_min"    : Fig 4.31, 15 voussoirs, rm=10, t/R=0.25, H_min
  "smars"         : Fig 4.32, 16 voussoirs, rm=10, friction 20°

CONVENTION (MF-part-4.pdf §13.1.1) :
  - Noeud i sur le rayon moyen rm à l'angle phi_i
  - Intrados à r_int = rm - t/2, extrados à r_ext = rm + t/2
  - Poids exact du voussoir (formules intégrales, Appendix)
  - Centre de gravité exact (formules intégrales, Appendix)

OUTPUTS (format dict, compatible Module 2 + Module 3) :
  nodes : list of dict {Index, XG, YG, ZG, Boundary, Weight, ...}
  branches : list of dict {Index, Head_node, Tail_node, u, v, LH, ...}

Auteur  : Claude (Pipeline MaNACoH v4)
Date    : 2026-02-14
Licence : MIT
"""

import math


# ==================================================================
# FORMULES EXACTES DE GÉOMÉTRIE (Appendix du MANUAL + MF-part-4.pdf)
# ==================================================================

def _cg_cylindrical_voussoir(phi_a, phi_b, ri, re, depth=1.0,
                              simplify=False):
    """
    Centre de gravité et poids d'un voussoir cylindrique.

    Source : c1_CreateModelVinAbe.bas L.650-700 + MF-part-4.pdf Appendix A
    
    Voussoir entre phi_a et phi_b, entre ri et re, profondeur depth.
    Convention : z vers le haut, x horizontal.
    
    Volume = depth × ∫(phi_a→phi_b) ∫(ri→re) r dr dphi
           = depth × (re² - ri²)/2 × (phi_b - phi_a)
    
    xG = ∫ x dV / V = ∫ r·sin(phi) · r dr dphi · depth / V
    zG = ∫ z dV / V = ∫ r·cos(phi) · r dr dphi · depth / V
    """
    dphi = phi_b - phi_a
    rm = (ri + re) / 2
    t = re - ri

    # Volume exact
    vol = depth * (re**2 - ri**2) / 2 * dphi

    if simplify:
        # Formule simplifiée (SimplifyCentreGravityFormulas = True)
        phi_m = (phi_a + phi_b) / 2
        xG = rm * math.sin(phi_m)
        zG = rm * math.cos(phi_m)
    else:
        # Formule exacte
        # ∫(phi_a→phi_b) sin(phi) dphi = cos(phi_a) - cos(phi_b)
        # ∫(phi_a→phi_b) cos(phi) dphi = sin(phi_b) - sin(phi_a)
        # ∫(ri→re) r² dr = (re³ - ri³)/3
        r3 = (re**3 - ri**3) / 3
        r2 = (re**2 - ri**2) / 2
        
        sin_int = math.cos(phi_a) - math.cos(phi_b)
        cos_int = math.sin(phi_b) - math.sin(phi_a)
        
        xG = depth * r3 * sin_int / vol
        zG = depth * r3 * cos_int / vol

    return xG, 0.0, zG, vol


def _cg_spherical_voussoir(phi_a, phi_b, theta_a, theta_b, ri, re,
                            simplify=False):
    """
    Centre de gravité et volume d'un voussoir sphérique.

    Source : c1_CreateModelVinAbe.bas L.700+ + MF-part-4.pdf Appendix A
    
    Volume = ∫(phi_a→phi_b) ∫(theta_a→theta_b) ∫(ri→re) r² sin(phi) dr dtheta dphi
           = (re³ - ri³)/3 × (theta_b - theta_a) × (cos(phi_a) - cos(phi_b))
    """
    dtheta = theta_b - theta_a
    r3 = (re**3 - ri**3) / 3
    cos_int = math.cos(phi_a) - math.cos(phi_b)

    vol = r3 * dtheta * cos_int

    if simplify:
        phi_m = (phi_a + phi_b) / 2
        theta_m = (theta_a + theta_b) / 2
        rm = (ri + re) / 2
        xG = rm * math.sin(phi_m) * math.cos(theta_m)
        yG = rm * math.sin(phi_m) * math.sin(theta_m)
        zG = rm * math.cos(phi_m)
    else:
        # Formules intégrales exactes
        r4 = (re**4 - ri**4) / 4
        
        # ∫ sin²(phi) dphi = (phi - sin(phi)cos(phi))/2
        def int_sin2(a, b):
            return (b - a) / 2 - (math.sin(2*b) - math.sin(2*a)) / 4
        
        # ∫ sin(phi)cos(phi) dphi = sin²(phi)/2
        def int_sincos(a, b):
            return (math.sin(b)**2 - math.sin(a)**2) / 2

        sin2_int = int_sin2(phi_a, phi_b)
        sincos_int = int_sincos(phi_a, phi_b)

        sin_theta_int = math.sin(theta_b) - math.sin(theta_a)
        cos_theta_int = -(math.cos(theta_b) - math.cos(theta_a))
        # Note: ∫cos(theta)dtheta = sin(theta)

        xG = r4 * sin2_int * (math.sin(theta_b) - math.sin(theta_a)) / vol
        yG = r4 * sin2_int * (-(math.cos(theta_b) - math.cos(theta_a))) / vol
        zG = r4 * sincos_int * dtheta / vol

    return xG, yG, zG, vol


# ==================================================================
# GENERATEUR ARC CYLINDRIQUE
# ==================================================================

def create_cylindrical_arch(n=15, m=1, phi0=0.0, phi1=None,
                             y0=0.0, y1=1.0,
                             r=10.0, thickness=None, t_over_R=0.15,
                             rho=1000.0, g=9.81,
                             simplify_cg=False,
                             boundary_on_extrados=True):
    """
    Arc cylindrique (ou voûte en berceau si m > 1).
    
    Source : c1_CreateModelParameters.bas L.38-60 (default values)
             c1_CreateModelVinAbe.bas L.60-400 (a3_CreateRegularGeo)
    
    Params :
      n : nombre de voussoirs le long de phi (default 15)
      m : nombre de voussoirs le long de y (default 1 = arc 2D)
      phi0, phi1 : bornes angulaires en radians [0, π/2]
      r : rayon moyen
      thickness : épaisseur (si None, calculé depuis t_over_R)
      t_over_R : ratio épaisseur/rayon
      rho : masse volumique [kg/m³]
      boundary_on_extrados : MeridianBoundaryNodesOnExtrados (MANUAL §2.1.2)
    
    Returns :
      nodes, branches : listes de dicts compatibles Module 2/3
    """
    if phi1 is None:
        phi1 = math.pi / 2

    if thickness is None:
        thickness = t_over_R * r

    ri = r - thickness / 2
    re = r + thickness / 2

    nodes = []
    branches = []
    node_id = 1
    
    # Nombre de noeuds : (n+1) le long de phi × (m+1) le long de y
    # Pour m=1 (arc 2D) : (n+1) × 2 = n+1 noeuds sur 1 tranche
    # Mais en 2D, on n'a que n+1 noeuds au total (1 tranche y)
    # → Pour reproduire le VBA exactement : m=1, y dimension ignorée

    dphi = (phi1 - phi0) / n
    dy = (y1 - y0) / max(m, 1)

    # Grille de noeuds [i_phi, j_y]
    # i_phi : 0 → n (n+1 noeuds le long du méridien)
    # j_y : 0 → m (m+1 noeuds le long de y, mais pour m=1, on collapse)
    
    # Pour un arc 2D (m=1), on a n+1 noeuds et n branches
    # Boundary : phi=phi0 (clé) et phi=phi1 (naissance)
    
    # Noeud boundary position (MANUAL Fig 2.3)
    if boundary_on_extrados:
        r_bou = re
    else:
        r_bou = r  # mid-plane
    
    node_grid = {}  # (i_phi, j_y) → node_id
    
    for j_y in range(m + 1):
        y = y0 + j_y * dy
        for i_phi in range(n + 1):
            phi = phi0 + i_phi * dphi

            # Déterminer si boundary
            is_bou_phi = (i_phi == 0 or i_phi == n)
            is_bou_y = (m > 1 and (j_y == 0 or j_y == m))
            is_bou = is_bou_phi  # Pour arc 2D
            
            if is_bou:
                r_node = r_bou
            else:
                r_node = r
            
            x = r_node * math.sin(phi)
            z = r_node * math.cos(phi)
            
            # Intrados/extrados radial du noeud (pour CSG)
            xi = ri * math.sin(phi)
            zi_int = ri * math.cos(phi)
            xe = re * math.sin(phi)
            ze_ext = re * math.cos(phi)

            # Poids : voussoir centré sur ce noeud (sauf boundary)
            if is_bou:
                weight = 0.0
                xG, yG, zG = x, y, z
            else:
                # Voussoir entre phi - dphi/2 et phi + dphi/2
                pa = phi - dphi / 2
                pb = phi + dphi / 2
                xG, yG, zG, vol = _cg_cylindrical_voussoir(
                    pa, pb, ri, re, depth=dy if m > 1 else (y1-y0),
                    simplify=simplify_cg)
                weight = rho * g * vol
            
            node = {
                'Index': node_id,
                'XG': xG if not is_bou else x,
                'YG': yG if not is_bou else y,
                'ZG': zG if not is_bou else z,
                'Boundary': is_bou,
                'Weight': weight,
                # Intrados/extrados verticaux (pour CSG sur noeuds)
                'BXi': xi, 'BYi': y, 'BZi': zi_int,
                'BXe': xe, 'BYe': y, 'BZe': ze_ext,
                # Épaisseur locale
                'Thickness': thickness,
                # Angles
                'phi': phi, 'y_coord': y,
            }
            nodes.append(node)
            node_grid[(i_phi, j_y)] = node_id
            node_id += 1

    # Branches le long du méridien (phi direction)
    branch_id = 1
    for j_y in range(m + 1):
        for i_phi in range(n):
            h = node_grid[(i_phi, j_y)]
            t_node = node_grid[(i_phi + 1, j_y)]
            nh = nodes[h - 1]
            nt = nodes[t_node - 1]

            # Joint au milieu entre les 2 noeuds, sur intrados/extrados
            phi_mid = (nh['phi'] + nt['phi']) / 2
            
            b = _make_branch(branch_id, nh, nt, phi_mid, ri, re)
            branches.append(b)
            branch_id += 1
    
    # Branches le long de y (hoop, si m > 1)
    if m > 1:
        for j_y in range(m):
            for i_phi in range(n + 1):
                h = node_grid[(i_phi, j_y)]
                t_node = node_grid[(i_phi, j_y + 1)]
                nh = nodes[h - 1]
                nt = nodes[t_node - 1]
                
                phi_mid = nh['phi']
                b = _make_branch(branch_id, nh, nt, phi_mid, ri, re)
                branches.append(b)
                branch_id += 1

    # Filtrer pour m=1 : ne garder que j_y=0
    if m <= 1:
        # On a créé 2 tranches y (j_y=0 et j_y=1) mais n'en veut qu'une
        # En fait pour m=1 on a m+1=2 tranches, on ne garde que j_y=0
        kept_ids = set(node_grid[(i, 0)] for i in range(n + 1))
        nodes = [n for n in nodes if n['Index'] in kept_ids]
        branches = [b for b in branches 
                   if b['Head_node'] in kept_ids and b['Tail_node'] in kept_ids]

    return nodes, branches


def _make_branch(branch_id, nh, nt, phi_mid, ri, re):
    """Crée une branche avec ses joints intrados/extrados."""
    # Convention VBA (MANUAL Table 2.2) : u = x_head - x_tail
    dx = nh['XG'] - nt['XG']
    dy = nh['YG'] - nt['YG']
    dz = nh['ZG'] - nt['ZG']
    lh = math.sqrt(dx**2 + dy**2)
    L = math.sqrt(dx**2 + dy**2 + dz**2)

    # Joint : sur intrados et extrados à phi_mid
    ji_x = ri * math.sin(phi_mid)
    ji_z = ri * math.cos(phi_mid)
    je_x = re * math.sin(phi_mid)
    je_z = re * math.cos(phi_mid)

    return {
        'Index': branch_id,
        'Head_node': nh['Index'],
        'Tail_node': nt['Index'],
        'u': dx, 'v': dy, 'LH': max(lh, 1e-12), 'L': L,
        'xi': ji_x, 'yi': (nh['YG'] + nt['YG']) / 2, 'zi': ji_z,
        'xe': je_x, 'ye': (nh['YG'] + nt['YG']) / 2, 'ze': je_z,
    }


# ==================================================================
# CAS PRÉDÉFINIS (MANUAL Chapter 4)
# ==================================================================

PRESETS = {
    # MANUAL Fig 4.27 : épaisseur minimale d'un arc
    "heyman_15": dict(
        n=15, m=1, phi0=0, phi1=math.pi/2,
        r=10, t_over_R=0.15, rho=1000,
        simplify_cg=False, boundary_on_extrados=True,
        description="Heyman Fig 4.27: t/R=0.15, CSG_J=1.3952",
        expected_csg=1.3952,
    ),
    # MANUAL Fig 4.28 : épaisseur minimale
    "heyman_tmin": dict(
        n=15, m=1, phi0=0, phi1=math.pi/2,
        r=10, t_over_R=0.1075, rho=1000,
        simplify_cg=False, boundary_on_extrados=True,
        description="Heyman Fig 4.28: t/R=0.1075, CSG_J≈1.0004",
        expected_csg=1.0004,
    ),
    # MANUAL Fig 4.31 : poussée minimale
    "oikono_min": dict(
        n=15, m=1, phi0=0, phi1=math.pi/2,
        r=10, t_over_R=0.25, rho=1000,
        simplify_cg=False, boundary_on_extrados=True,
        description="Oikonomopoulou Fig 4.31: min thrust, t/R=0.25",
        expected_csg=None,
    ),
    # MANUAL Fig 4.32 : Smars friction
    "smars": dict(
        n=16, m=1, phi0=-math.pi/2, phi1=math.pi/2,
        r=10, thickness=2*10*(1.5-1)/(1.5+1), rho=1000,
        simplify_cg=False, boundary_on_extrados=True,
        description="Smars 2000: full arch, friction 20°",
        expected_csg=None,
    ),
}


def create_preset(name):
    """Crée une géométrie prédéfinie. Retourne (nodes, branches, info)."""
    if name not in PRESETS:
        raise ValueError("Preset '{}' not found. Available: {}".format(
            name, list(PRESETS.keys())))
    p = PRESETS[name]
    desc = p.get('description', '')
    expected = p.get('expected_csg', None)
    build_args = {k: v for k, v in p.items()
                  if k not in ('description', 'expected_csg')}
    nodes, branches = create_cylindrical_arch(**build_args)
    return nodes, branches, {'description': desc, 'expected_csg': expected}


# ==================================================================
# TEST AUTONOME
# ==================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Module 1 — Test géométrie paramétrique")
    print("=" * 60)

    for name in ["heyman_15", "heyman_tmin", "smars"]:
        nodes, branches, info = create_preset(name)
        int_n = [n for n in nodes if not n['Boundary']]
        bou_n = [n for n in nodes if n['Boundary']]
        W = sum(n['Weight'] for n in int_n)
        print("\n  {}: {}".format(name, info['description']))
        print("    {} noeuds ({} int + {} bou), {} branches".format(
            len(nodes), len(int_n), len(bou_n), len(branches)))
        print("    Poids total: {:.1f} N".format(W))
        if bou_n:
            print("    zb = {}".format(["{:.4f}".format(n['ZG']) for n in bou_n]))
        if info['expected_csg']:
            print("    CSG attendu: {}".format(info['expected_csg']))
