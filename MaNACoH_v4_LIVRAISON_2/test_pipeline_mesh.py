#!/usr/bin/env python3
"""
TEST PIPELINE MESH → TNA (hors Rhino)
=======================================
Simule le workflow complet :
  1. Un mesh (simulé comme liste de vertices + faces)
  2. Paramètres physiques : t, rho, g
  3. Appuis : indices de vertex
  4. → Module 0 v4 logic → Module 2 → Module 3 → CSG

Cas test : arc semi-circulaire 15 voussoirs (Heyman)
  reconstruit "comme si" l'utilisateur l'avait modélisé en Rhino
"""
import math, time

# Charger Module 3
with open('/home/claude/MaNACoH_v4_03_TNA_v51.py', 'r') as f:
    src = f.read()
cut = src.find("# Valeurs par défaut")
if cut > 0:
    src = src[:cut]
exec(src)


# ======================================================================
# SIMULER UN MESH RHINO (arc semi-circulaire sur la surface moyenne)
# ======================================================================

def simulate_mesh_arch(n=15, R=10.0, t_val=1.5, phi0_deg=0, phi1_deg=90):
    """
    Simule un mesh 1D (chaîne de segments) d'arc semi-circulaire.
    Retourne :
      vertices : [(x,y,z), ...] — positions sur la surface moyenne
      edges : [(i,j), ...] — connectivité
      normals : [(nx,ny,nz), ...] — normales (radiales sortantes)
      faces : [(i,j,k), ...] — faces triangulaires (pour aires)
    """
    phi0 = math.radians(phi0_deg)
    phi1 = math.radians(phi1_deg)
    dphi = (phi1 - phi0) / n

    vertices = []
    normals = []
    for i in range(n + 1):
        phi = phi0 + i * dphi
        x = R * math.sin(phi)
        z = R * math.cos(phi)
        vertices.append((x, 0.0, z))
        # Normale radiale (pointe vers l'extérieur)
        normals.append((math.sin(phi), 0.0, math.cos(phi)))

    edges = [(i, i + 1) for i in range(n)]

    # Faces fictives pour le calcul d'aire tributaire
    # Pour un arc 1D, l'aire tributaire = dphi × R × depth
    # On simule des triangles dégénérés
    faces = []  # pas de faces réelles en 1D

    return vertices, edges, normals, faces


def mesh_to_model(vertices, edges, normals, supports, t_val, rho_val,
                  g_val=9.81, depth=1.0):
    """
    Reproduit la logique du Module 0 v4 build_model() SANS Rhino.

    Pour chaque nœud intérieur :
      Weight = rho × g × aire_tributaire × t
      Bzi = ZG - t/2, Bze = ZG + t/2

    Pour chaque branche :
      Joint = milieu ± t/2 × n̂
    """
    nv = len(vertices)
    ri = (t_val)  # épaisseur, pas rayon

    # Séparer intérieurs / boundary
    interior = [i for i in range(nv) if i not in supports]
    boundary = [i for i in range(nv) if i in supports]

    # Index nœud (int d'abord, 1-based)
    v2n = {}
    idx = 1
    for vi in interior:
        v2n[vi] = idx; idx += 1
    ni = len(interior)
    for vi in boundary:
        v2n[vi] = idx; idx += 1
    nb = len(boundary)

    # Aire tributaire 1D : longueur des segments adjacents / 2
    seg_lengths = {}
    for (vi1, vi2) in edges:
        v1, v2 = vertices[vi1], vertices[vi2]
        L = math.sqrt(sum((a-b)**2 for a, b in zip(v1, v2)))
        seg_lengths.setdefault(vi1, []).append(L)
        seg_lengths.setdefault(vi2, []).append(L)

    def trib_area(vi):
        """Aire tributaire = somme(L_adj)/2 × depth."""
        segs = seg_lengths.get(vi, [])
        return sum(segs) / 2.0 * depth if segs else 0.0

    # Nœuds
    nodes = []
    for vi in interior:
        x, y, z = vertices[vi]
        nx, ny, nz = normals[vi]
        area = trib_area(vi)
        w = rho_val * g_val * area * t_val

        nodes.append({
            'Index': v2n[vi], 'VertexIndex': vi,
            'XG': x, 'YG': y, 'ZG': z,
            'Boundary': False, 'Weight': w, 'Thickness': t_val,
            'Bzi': z - t_val / 2, 'Bze': z + t_val / 2,
            'BXi': x - t_val/2*nx, 'BYi': y - t_val/2*ny, 'BZi': z - t_val/2*nz,
            'BXe': x + t_val/2*nx, 'BYe': y + t_val/2*ny, 'BZe': z + t_val/2*nz,
        })
    for vi in boundary:
        x, y, z = vertices[vi]
        nx, ny, nz = normals[vi]
        nodes.append({
            'Index': v2n[vi], 'VertexIndex': vi,
            'XG': x, 'YG': y, 'ZG': z,
            'Boundary': True, 'Weight': 0.0, 'Thickness': t_val,
            'Bzi': z - t_val / 2, 'Bze': z + t_val / 2,
            'BXi': x - t_val/2*nx, 'BYi': y - t_val/2*ny, 'BZi': z - t_val/2*nz,
            'BXe': x + t_val/2*nx, 'BYe': y + t_val/2*ny, 'BZe': z + t_val/2*nz,
        })

    # Branches
    branches = []
    bidx = 1
    for (vi1, vi2) in edges:
        v1 = vertices[vi1]; v2 = vertices[vi2]
        dx, dy, dz = v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]
        L = math.sqrt(dx**2 + dy**2 + dz**2)
        LH = math.sqrt(dx**2 + dy**2)

        mx = (v1[0]+v2[0])/2; my = (v1[1]+v2[1])/2; mz = (v1[2]+v2[2])/2
        n1 = normals[vi1]; n2 = normals[vi2]
        nmx = (n1[0]+n2[0])/2; nmy = (n1[1]+n2[1])/2; nmz = (n1[2]+n2[2])/2
        nml = math.sqrt(nmx**2+nmy**2+nmz**2)
        if nml > 1e-12:
            nmx/=nml; nmy/=nml; nmz/=nml

        branches.append({
            'Index': bidx,
            'Head_node': v2n[vi1], 'Tail_node': v2n[vi2],
            'u': dx, 'v': dy, 'w': dz,
            'L': L, 'LH': max(LH, 1e-12),
            'xi': mx - t_val/2*nmx, 'yi': my - t_val/2*nmy, 'zi': mz - t_val/2*nmz,
            'xe': mx + t_val/2*nmx, 'ye': my + t_val/2*nmy, 'ze': mz + t_val/2*nmz,
            'Surface': L * t_val,
        })
        bidx += 1

    return nodes, branches, ni, nb


def run_tna(nodes, branches, ni, nb, objective="max_csg", max_iter=200):
    """Module 2 → Module 3 complet."""
    int_n = [n for n in nodes if not n['Boundary']]
    bou_n = [n for n in nodes if n['Boundary']]
    all_n = int_n + bou_n
    idx_map = {n['Index']: i for i, n in enumerate(all_n)}
    m = len(branches)

    Ci = [[0.0]*ni for _ in range(m)]
    Cb = [[0.0]*nb for _ in range(m)]
    u_vec, v_vec, LH_vec = [], [], []

    for j, b in enumerate(branches):
        h, t_n = idx_map[b['Head_node']], idx_map[b['Tail_node']]
        if h < ni: Ci[j][h] = +1.0
        else: Cb[j][h - ni] = +1.0
        if t_n < ni: Ci[j][t_n] = -1.0
        else: Cb[j][t_n - ni] = -1.0
        u_vec.append(b['u']); v_vec.append(b['v']); LH_vec.append(b['LH'])

    pz = [0.0]*ni
    for nd in int_n:
        idx = idx_map[nd['Index']]
        if idx < ni:
            pz[idx] = -nd['Weight']
    zb = [nd['ZG'] for nd in bou_n]

    node_valence = [0]*ni
    for j in range(m):
        for i in range(ni):
            if abs(Ci[j][i]) > 0.5:
                node_valence[i] += 1
    is_vertical = [abs(LH_vec[j]) < 1e-12 for j in range(m)]

    zeta_init, dof_global, unk_global, M_stored = solve_horizontal(
        Ci, Cb, u_vec, v_vec, LH_vec, node_valence, ni, m)

    t0 = time.time()
    zeta_opt, zi_opt, q_opt, zb_opt = optimize_dof(
        nodes, branches, Ci, Cb, u_vec, v_vec, LH_vec, node_valence,
        ni, nb, m, ni, pz, zb,
        is_vertical, M_stored, dof_global, unk_global,
        max_iter=max_iter, objective=objective)
    dt = time.time() - t0

    Z = zi_opt[:] + zb_opt[:]
    zG = [nd['ZG'] for nd in int_n]
    csg_list = compute_csg_simplified(nodes, branches, Z, zeta_opt, ni)
    csg_min = min(csg_list) if csg_list else 0
    prox = sum((zi_opt[i]-zG[i])**2 for i in range(ni))

    return {
        'zeta': zeta_opt, 'zi': zi_opt, 'zb': zb_opt,
        'csg_list': csg_list, 'csg_min': csg_min,
        'prox': prox, 'dt': dt,
        'n_dof': len([g for g in dof_global if not is_vertical[g]]),
    }


# ======================================================================
# TEST 1 : Arc Heyman depuis "mesh"
# ======================================================================

print("=" * 70)
print("TEST PIPELINE : Mesh → Modèle → TNA → CSG")
print("=" * 70)

R = 10.0
t_val = 0.15 * R  # t/R = 0.15
rho_val = 1000.0

print("\n--- Cas 1 : Arc Heyman (mesh simulé) ---")
print("  R={}, t={}, t/R={}, rho={}".format(R, t_val, t_val/R, rho_val))

# Simuler mesh
vertices, edges, normals, faces = simulate_mesh_arch(
    n=15, R=R, t_val=t_val, phi0_deg=0, phi1_deg=90)
supports = {0, 15}  # clé et naissance

print("  Mesh: {} vertices, {} edges, supports={}".format(
    len(vertices), len(edges), supports))

# Module 0 v4 : mesh → modèle
nodes, branches, ni, nb = mesh_to_model(
    vertices, edges, normals, supports, t_val, rho_val)
W = sum(nd['Weight'] for nd in nodes if not nd['Boundary'])
print("  Modèle: {} int + {} bou, {} branches, W={:.1f} N".format(
    ni, nb, len(branches), W))

# Module 3 : proximity
r_prox = run_tna(nodes, branches, ni, nb, "proximity")
print("  Proximity: cost={:.4f}, CSG={:.4f}, dt={:.3f}s".format(
    r_prox['prox'], r_prox['csg_min'], r_prox['dt']))

# Module 3 : max_csg
r_csg = run_tna(nodes, branches, ni, nb, "max_csg")
print("  Max CSG:   CSG={:.4f}, dt={:.3f}s".format(
    r_csg['csg_min'], r_csg['dt']))
print("  CSG attendu: 1.3952 (MANUAL)")
print("  t/R/CSG = {:.5f} (attendu: 0.10751)".format(
    t_val/R / r_csg['csg_min'] if r_csg['csg_min'] > 0 else 0))

# CSG par joint
print("\n  CSG par joint:")
for j, c in enumerate(r_csg['csg_list']):
    b = branches[j]
    mx = (b['xi']+b['xe'])/2
    mz = (b['zi']+b['ze'])/2
    phi = math.degrees(math.atan2(mx, mz)) if (mx != 0 or mz != 0) else 0
    marker = " ◄" if abs(c - r_csg['csg_min']) < 0.02 else ""
    print("    J{:2d} (φ={:5.1f}°): CSG={:.4f}{}".format(j, phi, c, marker))


# ======================================================================
# TEST 2 : Même arc mais épaisseur variable (plus épais aux naissances)
# ======================================================================

print("\n--- Cas 2 : Épaisseur variable ---")

t_list = []
for i in range(16):
    phi = i * 90.0 / 15
    # Plus épais aux naissances (phi=90°)
    t_var = t_val * (1.0 + 0.5 * math.sin(math.radians(phi)))
    t_list.append(t_var)
print("  t varie de {:.3f} à {:.3f} m".format(min(t_list), max(t_list)))

# Reconstruire avec épaisseur variable
nodes2 = []
for nd in nodes:
    nd2 = dict(nd)
    vi = nd['VertexIndex']
    ti = t_list[vi]
    if not nd2['Boundary']:
        # Recalculer le poids avec la nouvelle épaisseur
        # On cherche l'aire tributaire depuis le poids original
        area = nd['Weight'] / (rho_val * 9.81 * t_val)
        nd2['Weight'] = rho_val * 9.81 * area * ti
    nd2['Thickness'] = ti
    nd2['Bzi'] = nd2['ZG'] - ti/2
    nd2['Bze'] = nd2['ZG'] + ti/2
    nodes2.append(nd2)

branches2 = []
for b in branches:
    b2 = dict(b)
    # Recalculer joints avec épaisseur moyenne des 2 extrémités
    vi1 = [nd for nd in nodes if nd['Index'] == b['Head_node']][0]['VertexIndex']
    vi2 = [nd for nd in nodes if nd['Index'] == b['Tail_node']][0]['VertexIndex']
    tm = (t_list[vi1] + t_list[vi2]) / 2
    # Recalculer les joints (même logique que build_model)
    v1 = vertices[vi1]; v2 = vertices[vi2]
    mx = (v1[0]+v2[0])/2; mz = (v1[2]+v2[2])/2
    n1 = normals[vi1]; n2 = normals[vi2]
    nmx = (n1[0]+n2[0])/2; nmz = (n1[2]+n2[2])/2
    nml = math.sqrt(nmx**2+nmz**2)
    if nml > 1e-12:
        nmx/=nml; nmz/=nml
    b2['xi'] = mx - tm/2*nmx; b2['zi'] = mz - tm/2*nmz
    b2['xe'] = mx + tm/2*nmx; b2['ze'] = mz + tm/2*nmz
    branches2.append(b2)

r2 = run_tna(nodes2, branches2, ni, nb, "max_csg")
print("  Max CSG (t variable): {:.4f}, dt={:.3f}s".format(
    r2['csg_min'], r2['dt']))
print("  (épaisseur variable → CSG change car joints plus larges)")


# ======================================================================
# TEST 3 : Arc plein (phi: -90° → +90°) — Smars-like
# ======================================================================

print("\n--- Cas 3 : Arc plein -90°/+90° (type Smars) ---")

vertices3, edges3, normals3, faces3 = simulate_mesh_arch(
    n=16, R=10.0, t_val=2*10*(1.5-1)/(1.5+1), phi0_deg=-90, phi1_deg=90)
t3 = 2 * 10.0 * (1.5 - 1) / (1.5 + 1)
supports3 = {0, 16}

nodes3, branches3, ni3, nb3 = mesh_to_model(
    vertices3, edges3, normals3, supports3, t3, 1000.0)
W3 = sum(nd['Weight'] for nd in nodes3 if not nd['Boundary'])
print("  {} int + {} bou, {} branches, W={:.1f} N, t={:.3f}".format(
    ni3, nb3, len(branches3), W3, t3))

r3 = run_tna(nodes3, branches3, ni3, nb3, "max_csg")
print("  Max CSG: {:.4f}, dt={:.3f}s".format(r3['csg_min'], r3['dt']))


# ======================================================================
# RÉSUMÉ
# ======================================================================

print("\n" + "=" * 70)
print("RÉSUMÉ")
print("=" * 70)
print("{:>30s} {:>6s} {:>8s} {:>8s}".format("Test", "DOF", "CSG", "Temps"))
print("{:>30s} {:6d} {:8.4f} {:7.3f}s".format(
    "Heyman mesh proximity", r_prox['n_dof'], r_prox['csg_min'], r_prox['dt']))
print("{:>30s} {:6d} {:8.4f} {:7.3f}s".format(
    "Heyman mesh max_csg", r_csg['n_dof'], r_csg['csg_min'], r_csg['dt']))
print("{:>30s} {:6d} {:8.4f} {:7.3f}s".format(
    "Heyman t_variable max_csg", r2['n_dof'], r2['csg_min'], r2['dt']))
print("{:>30s} {:6d} {:8.4f} {:7.3f}s".format(
    "Smars plein max_csg", r3['n_dof'], r3['csg_min'], r3['dt']))
