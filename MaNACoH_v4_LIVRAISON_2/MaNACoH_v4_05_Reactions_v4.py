#! python3
"""
MaNACoH v4 - Module 5 : REACTIONS D'APPUI (Version 3 - Équilibre Direct)
=========================================================================
Calcul des réactions d'appui par ÉQUILIBRE NODAL DIRECT

RÉFÉRENCES:
  - Fantin (2017) "Méthode des réseaux de forces" §4.3.4-4.3.5
  - SETRA (1982) "Ponts en maçonnerie" p.92

=============================================================================
THÉORIE - APPROCHE PAR ÉQUILIBRE NODAL
=============================================================================

Pour chaque noeud INTÉRIEUR i, l'équilibre des forces est:

  Σ F_branches + P_i = 0

où:
  - F_branches = somme des forces dans les branches connectées
  - P_i = (0, 0, -Weight_i) = poids du bloc (vers le bas)

Pour les noeuds BOUNDARY, il n'y a pas de poids propre associé,
donc la réaction R_b équilibre exactement les forces des branches:

  R_b = - Σ F_branches_connectées

CALCUL DES FORCES DANS LES BRANCHES:
La force dans une branche est PARALLÈLE à la branche (réseau de forces).

Si on connaît q (densité de force) et Z, alors:
  - Force horizontale: l*_H = q · l_H
  - Force totale: l* = l*_H · l / l_H = q · l
  - Composantes: (u*, v*, w*) = l* · (u, v, w) / l

VÉRIFICATION GLOBALE:
  Σ R_V (toutes réactions) = Σ Poids (tous blocs intérieurs)

=============================================================================
INPUTS:
  - nodes : list/str - Noeuds du Module 1
  - branches : list/str - Branches du Module 1  
  - Z : list - Altitudes du thrust network (Module 3)
  - q : list - Densités de force (Module 3) [REQUIS]
  - vector_scale : float - Échelle des vecteurs [optionnel, défaut=auto]

OUTPUTS DONNÉES:
  - reactions : list - Dict par noeud d'appui
  - total_RV : float - Somme réactions verticales
  - total_weight : float - Poids total
  - equilibrium_ok : bool - Équilibre vérifié?

OUTPUTS VISUALISATION:
  - reaction_pts : Point3d[] - Points d'application
  - reaction_lines : Line[] - Lignes réaction totale R
  - reaction_H_lines : Line[] - Lignes composante horizontale RH
  - reaction_V_lines : Line[] - Lignes composante verticale RV
  - text_dots_RH : TextDot[] - Labels "RH=xxx kN"
  - text_dots_RV : TextDot[] - Labels "RV=xxx kN"
  - text_dots_R : TextDot[] - Labels "R=xxx kN"
"""

import Rhino.Geometry as rg
import math
import json

# === PARAMÈTRE ÉCHELLE ===
if 'vector_scale' not in dir() or vector_scale is None:
    vector_scale = 0  # 0 = auto

# === FONCTIONS ===

def compute_branch_force_from_q(b, n1, n2, Z, q_vec, idx_map, branch_idx):
    """
    Calcule la force dans une branche à partir de la densité de force q.
    
    F = q · L (dans la direction de la branche)
    
    Composantes:
      u* = F · (u/L) = q · u
      v* = F · (v/L) = q · v
      w* = F · (w/L) = q · w
    
    où (u, v, w) sont les composantes de la branche dans le réseau de FORCES.
    """
    # Indices des noeuds
    i1 = idx_map[b['Head_node']]
    i2 = idx_map[b['Tail_node']]
    
    # Composantes de la branche du thrust (réseau de forces)
    # Convention MANUAL Table 2.2 : u = x_head - x_tail
    # E7 fix : cohérent avec E1
    u = n1['XG'] - n2['XG']  # x_head - x_tail
    v = n1['YG'] - n2['YG']  # y_head - y_tail
    w = Z[i1] - Z[i2]        # z_head - z_tail
    
    # Densité de force - utiliser branch_idx (position dans la liste)
    if q_vec and branch_idx < len(q_vec):
        q = q_vec[branch_idx]
    else:
        q = 1.0  # Fallback
    
    # Composantes de la force: F_comp = q * branch_comp
    # (car F = q*L et F_comp/F = branch_comp/L)
    u_star = q * u
    v_star = q * v
    w_star = q * w
    
    return u_star, v_star, w_star


def compute_reactions_from_equilibrium(nodes_list, branches_list, Z, q_vec):
    """
    Calcule les réactions d'appui par équilibre nodal.
    
    Pour chaque noeud boundary, la réaction équilibre
    les forces des branches connectées.
    """
    # Séparer noeuds
    int_n = [n for n in nodes_list if not n['Boundary']]
    bou_n = [n for n in nodes_list if n['Boundary']]
    all_n = int_n + bou_n
    
    idx_map = {n['Index']: i for i, n in enumerate(all_n)}
    nmap = {n['Index']: n for n in nodes_list}
    
    # DEBUG: Afficher infos sur q
    print("\nDEBUG - Vecteur q:")
    print("  Longueur q: {}".format(len(q_vec) if q_vec else 0))
    print("  Nb branches: {}".format(len(branches_list)))
    if q_vec and len(q_vec) > 0:
        print("  q[0:5]: {}".format(q_vec[:min(5, len(q_vec))]))
        print("  q min/max: {:.2f} / {:.2f}".format(min(q_vec), max(q_vec)))
    
    # Initialiser les sommes de forces par noeud
    node_forces = {n['Index']: [0.0, 0.0, 0.0] for n in nodes_list}
    
    # DEBUG: compter les branches connectées aux appuis
    boundary_branches = 0
    
    # Parcourir les branches et accumuler les forces
    for branch_idx, b in enumerate(branches_list):
        head_idx = b['Head_node']
        tail_idx = b['Tail_node']
        n1 = nmap[head_idx]
        n2 = nmap[tail_idx]
        
        # Force dans la branche (Head -> Tail)
        u_star, v_star, w_star = compute_branch_force_from_q(
            b, n1, n2, Z, q_vec, idx_map, branch_idx
        )
        
        # DEBUG: afficher pour les premières branches aux appuis
        if (n1['Boundary'] or n2['Boundary']) and boundary_branches < 3:
            print("  Branch {}: Head={}, Tail={}, Force=({:.1f}, {:.1f}, {:.1f})".format(
                branch_idx, head_idx, tail_idx, u_star, v_star, w_star
            ))
            boundary_branches += 1
        
        # La force sur le noeud Head = +C[head] * F = +1 * (u*, v*, w*)
        # La force sur le noeud Tail = +C[tail] * F = -1 * (u*, v*, w*)
        # E7 fix : cohérent avec C[j,head]=+1, C[j,tail]=-1
        # et u = x_head - x_tail (MANUAL Table 2.2)
        node_forces[head_idx][0] += u_star
        node_forces[head_idx][1] += v_star
        node_forces[head_idx][2] += w_star
        
        node_forces[tail_idx][0] -= u_star
        node_forces[tail_idx][1] -= v_star
        node_forces[tail_idx][2] -= w_star
    
    # Vérifier l'équilibre des noeuds intérieurs
    # Équilibre : Σ C[j,i] * F_j + pz_i = 0
    # Avec pz = -Weight (VBA d_ComputeTypology L.244) :
    #   Σ fz_branches - Weight = 0  →  fz = Weight
    # E8 fix : résidu = |fz - Weight| (correct, commentaire clarifié)
    equilibrium_errors = []
    for n in int_n:
        fx, fy, fz = node_forces[n['Index']]
        weight = n['Weight']
        # Résidu : fx=0, fy=0, fz=Weight pour un noeud en équilibre
        residual = math.sqrt(fx**2 + fy**2 + (fz - weight)**2)
        equilibrium_errors.append(residual)
    
    avg_error = sum(equilibrium_errors) / len(equilibrium_errors) if equilibrium_errors else 0
    max_error = max(equilibrium_errors) if equilibrium_errors else 0
    
    # Pour les noeuds boundary, la réaction = -Σ F_branches
    reactions = []
    for n in bou_n:
        fx, fy, fz = node_forces[n['Index']]
        
        # Réaction = opposé des forces des branches
        # (les branches tirent sur l'appui, l'appui pousse en retour)
        Rx = -fx
        Ry = -fy
        Rz = -fz
        
        RH = math.sqrt(Rx**2 + Ry**2)
        R = math.sqrt(Rx**2 + Ry**2 + Rz**2)
        angle = math.degrees(math.atan2(abs(Rz), RH)) if RH > 1e-10 else 90
        
        reactions.append({
            'Index': n['Index'],
            'X': n['XG'], 'Y': n['YG'], 'Z': n['ZG'],
            'Rx': Rx, 'Ry': Ry, 'Rz': Rz,
            'RH': RH, 'RV': Rz, 'R': R,
            'angle_deg': angle
        })
    
    return reactions, avg_error, max_error


# === MAIN ===

# Parser entrées
nodes_list = None
branches_list = None

if 'nodes' in dir() and nodes:
    if isinstance(nodes, str):
        try:
            nodes_list = json.loads(nodes)
        except:
            pass
    elif isinstance(nodes, list):
        nodes_list = nodes

if 'branches' in dir() and branches:
    if isinstance(branches, str):
        try:
            branches_list = json.loads(branches)
        except:
            pass
    elif isinstance(branches, list):
        branches_list = branches

# q est REQUIS - mais on peut le calculer depuis LH_dual si besoin
q_vec = None
if 'q' in dir() and q:
    if isinstance(q, list) and len(q) > 0 and q[0] is not None:
        q_vec = q

# Fallback: calculer q depuis LH_dual et LH_vec
if not q_vec and 'LH_dual' in dir() and LH_dual:
    print("INFO: q non fourni, calcul depuis LH_dual...")
    # q = LH_dual / LH
    if 'branches' in dir() and branches:
        branches_list_temp = branches if isinstance(branches, list) else json.loads(branches)
        q_vec = []
        for i, b in enumerate(branches_list_temp):
            LH = b.get('LH', 1.0)
            if LH > 1e-10 and i < len(LH_dual):
                q_vec.append(LH_dual[i] / LH)
            else:
                q_vec.append(1.0)
        print("  q calculé: {} valeurs".format(len(q_vec)))

has_data = all([
    nodes_list and len(nodes_list) > 0,
    branches_list and len(branches_list) > 0,
    'Z' in dir() and Z and len(Z) > 0,
    q_vec and len(q_vec) > 0  # q_vec peut venir de q ou de LH_dual
])

if not has_data:
    print("ERREUR: Connectez nodes, branches, Z, et (q OU LH_dual)")
    print("  Option 1: Connectez 'q' du Module 3 (recommandé)")
    print("  Option 2: Connectez 'LH_dual' du Module 3 (fallback)")
    reactions = []
    total_RH = total_RV = total_weight = 0
    equilibrium_ok = False
    equilibrium_error = 100
    reaction_pts = reaction_vectors = reaction_lines = []
    reaction_H_vectors = reaction_V_vectors = []
    reaction_H_lines = reaction_V_lines = []
    text_dots_RH = text_dots_RV = text_dots_R = []
    arrow_H_pts = arrow_V_pts = arrow_R_pts = []
else:
    print("="*60)
    print("MaNACoH v4 - Module 5: REACTIONS D'APPUI (Equilibre Direct)")
    print("="*60)
    print("Méthode: Calcul par équilibre nodal avec q")
    print("")
    
    # Calculer réactions
    reactions, avg_error, max_error = compute_reactions_from_equilibrium(
        nodes_list, branches_list, Z, q_vec
    )
    
    # Poids total des noeuds intérieurs
    int_nodes = [n for n in nodes_list if not n['Boundary']]
    total_weight = sum(n['Weight'] for n in int_nodes)
    
    # Totaux réactions
    total_Rx = sum(r['Rx'] for r in reactions)
    total_Ry = sum(r['Ry'] for r in reactions)
    total_Rz = sum(r['Rz'] for r in reactions)
    total_RH = math.sqrt(total_Rx**2 + total_Ry**2)
    total_RV = total_Rz
    
    # Erreur d'équilibre global
    equilibrium_error = 100 * abs(total_RV - total_weight) / total_weight if total_weight > 0 else 100
    equilibrium_ok = equilibrium_error < 5
    
    # === VISUALISATION ===
    # Échelle des vecteurs
    max_R = max(r['R'] for r in reactions) if reactions else 1
    
    if vector_scale and vector_scale > 0:
        # Échelle fournie par l'utilisateur
        scale = vector_scale
        print("Échelle vecteurs: {} (utilisateur)".format(scale))
    else:
        # Échelle automatique (vecteur max ≈ 2-3m)
        scale = 2.5 / max_R if max_R > 0 else 1
        print("Échelle vecteurs: {:.2e} (auto)".format(scale))
    
    # Listes de sortie
    reaction_pts = []       # Points d'application
    reaction_vectors = []   # Vecteurs réaction totale
    reaction_H_vectors = [] # Vecteurs composante horizontale
    reaction_V_vectors = [] # Vecteurs composante verticale
    reaction_lines = []     # Lignes réaction totale
    reaction_H_lines = []   # Lignes composante horizontale
    reaction_V_lines = []   # Lignes composante verticale
    
    # Labels pour annotation
    label_pts = []          # Points pour les labels
    label_RH_texts = []     # Textes RH
    label_RV_texts = []     # Textes RV
    label_R_texts = []      # Textes R total
    
    # Points pour les flèches (extrémités)
    arrow_H_pts = []        # Extrémités flèches H
    arrow_V_pts = []        # Extrémités flèches V
    arrow_R_pts = []        # Extrémités flèches R
    
    for r in reactions:
        # Point d'application (noeud d'appui)
        pt = rg.Point3d(r['X'], r['Y'], r['Z'])
        reaction_pts.append(pt)
        
        # === RÉACTION TOTALE ===
        vec_R = rg.Vector3d(r['Rx']*scale, r['Ry']*scale, r['Rz']*scale)
        reaction_vectors.append(vec_R)
        
        pt_R_end = rg.Point3d(
            r['X'] + r['Rx']*scale,
            r['Y'] + r['Ry']*scale,
            r['Z'] + r['Rz']*scale
        )
        reaction_lines.append(rg.Line(pt, pt_R_end))
        arrow_R_pts.append(pt_R_end)
        
        # === COMPOSANTE HORIZONTALE (RH) ===
        vec_H = rg.Vector3d(r['Rx']*scale, r['Ry']*scale, 0)
        reaction_H_vectors.append(vec_H)
        
        pt_H_end = rg.Point3d(
            r['X'] + r['Rx']*scale,
            r['Y'] + r['Ry']*scale,
            r['Z']  # Reste à Z du point d'appui
        )
        reaction_H_lines.append(rg.Line(pt, pt_H_end))
        arrow_H_pts.append(pt_H_end)
        
        # === COMPOSANTE VERTICALE (RV) ===
        vec_V = rg.Vector3d(0, 0, r['Rz']*scale)
        reaction_V_vectors.append(vec_V)
        
        # RV part de l'extrémité de RH pour former le triangle
        pt_V_start = pt_H_end
        pt_V_end = pt_R_end  # Rejoint l'extrémité de R
        reaction_V_lines.append(rg.Line(pt_V_start, pt_V_end))
        arrow_V_pts.append(pt_V_end)
        
        # === LABELS ===
        # Position du label RH (milieu de la flèche H)
        pt_label_H = rg.Point3d(
            (r['X'] + pt_H_end.X) / 2,
            (r['Y'] + pt_H_end.Y) / 2,
            r['Z'] - 0.1  # Légèrement en dessous
        )
        label_pts.append(pt_label_H)
        label_RH_texts.append("RH={:.0f}kN".format(r['RH']/1000))
        
        # Position du label RV (milieu de la flèche V)
        pt_label_V = rg.Point3d(
            pt_H_end.X + 0.1,  # Légèrement décalé
            pt_H_end.Y,
            (pt_V_start.Z + pt_V_end.Z) / 2
        )
        label_RV_texts.append("RV={:.0f}kN".format(abs(r['RV'])/1000))
        
        # Label R total (à l'extrémité)
        label_R_texts.append("R={:.0f}kN ({:.0f}°)".format(r['R']/1000, r['angle_deg']))
    
    # === Créer les TextDots pour Grasshopper ===
    # (Les TextDots s'affichent toujours face à la caméra)
    text_dots_RH = []
    text_dots_RV = []
    text_dots_R = []
    
    for i, r in enumerate(reactions):
        # Position pour RH (au milieu de la flèche horizontale)
        pt_h = rg.Point3d(
            r['X'] + r['Rx']*scale*0.5,
            r['Y'] + r['Ry']*scale*0.5,
            r['Z']
        )
        text_dots_RH.append(rg.TextDot("RH={:.0f}".format(r['RH']/1000), pt_h))
        
        # Position pour RV (au milieu de la flèche verticale)
        pt_v = rg.Point3d(
            r['X'] + r['Rx']*scale,
            r['Y'] + r['Ry']*scale,
            r['Z'] + r['Rz']*scale*0.5
        )
        text_dots_RV.append(rg.TextDot("RV={:.0f}".format(abs(r['RV'])/1000), pt_v))
        
        # Position pour R (à l'extrémité)
        pt_r = rg.Point3d(
            r['X'] + r['Rx']*scale*1.1,
            r['Y'] + r['Ry']*scale*1.1,
            r['Z'] + r['Rz']*scale*1.1
        )
        text_dots_R.append(rg.TextDot("{:.0f}kN".format(r['R']/1000), pt_r))
    
    # === AFFICHAGE ===
    print("{:>5} {:>10} {:>10} {:>10} {:>10} {:>8}".format(
        "Node", "RH [kN]", "RV [kN]", "R [kN]", "Angle", "Dir"
    ))
    print("-"*60)
    
    for r in reactions:
        dir_angle = math.degrees(math.atan2(r['Ry'], r['Rx'])) if r['RH'] > 0 else 0
        print("{:>5} {:>10.1f} {:>10.1f} {:>10.1f} {:>7.1f}° {:>7.1f}°".format(
            r['Index'],
            r['RH']/1000, r['RV']/1000, r['R']/1000,
            r['angle_deg'], dir_angle
        ))
    
    print("-"*60)
    print("{:>5} {:>10.1f} {:>10.1f}".format("TOTAL", total_RH/1000, total_RV/1000))
    
    # Vérification
    print("")
    print("="*60)
    print("VERIFICATION EQUILIBRE")
    print("="*60)
    print("Poids total noeuds int.: {:>10.1f} kN".format(total_weight/1000))
    print("Σ RV (réactions)       : {:>10.1f} kN".format(total_RV/1000))
    print("Écart                  : {:>10.1f} kN ({:.1f}%)".format(
        abs(total_RV - total_weight)/1000, equilibrium_error
    ))
    print("")
    print("Erreur équilibre nodal:")
    print("  Moyenne : {:.1f} N".format(avg_error))
    print("  Maximum : {:.1f} N".format(max_error))
    print("")
    
    if equilibrium_ok:
        print("✓ EQUILIBRE GLOBAL VÉRIFIÉ")
    else:
        print("⚠️ EQUILIBRE NON SATISFAIT")
        print("")
        print("CAUSE PROBABLE:")
        print("  Le thrust network Z a été modifié par l'optimisation")
        print("  et ne satisfait plus exactement l'équilibre vertical.")
        print("  Les réactions affichées sont celles de l'équilibre nodal,")
        print("  pas les vraies réactions sur les appuis.")
    
    print("="*60)
