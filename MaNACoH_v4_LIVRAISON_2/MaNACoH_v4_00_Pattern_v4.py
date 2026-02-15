#! python3
"""
MaNACoH v4 - Module 0 : MESH → MODELE COMPLET (v4)
====================================================
Lit un mesh Rhino (surface moyenne de la voûte), ajoute l'épaisseur,
la densité et les conditions d'appui, et produit un modèle complet
prêt pour l'analyse TNA (Modules 2→3→4→5).

WORKFLOW :
  1. Modéliser la surface moyenne de la voûte comme un Mesh dans Rhino
  2. Brancher le mesh + déclarer les appuis (points, courbe, ou Z auto)
  3. Entrer l'épaisseur t [m] et la densité rho [kg/m³]
  4. → Le module produit nodes, branches avec poids, joints, bornes CSG

CALCULS PHYSIQUES (absents du v3) :
  - Weight = rho × g × Aire_tributaire × t       [MANUAL Table 2.1]
  - Bzi = ZG - t/2 (projection verticale intrados) [MANUAL Table 2.1]
  - Bze = ZG + t/2 (projection verticale extrados) [MANUAL Table 2.1]
  - xi,yi,zi = milieu_arête - t/2 × n̂            [MANUAL Table 2.2]
  - xe,ye,ze = milieu_arête + t/2 × n̂            [MANUAL Table 2.2]
  - Surface_joint ≈ longueur_arête × t            [MANUAL Table 2.2]

MODES D'APPUI (par priorité) :
  1. support_idx : indices vertex directs
  2. support_pts : Point3d → nœuds les plus proches
  3. support_crv : Curve → nœuds à distance < tol
  4. Z auto : Z ≤ Zmin + tol

ÉPAISSEUR :
  - t = scalaire → uniforme
  - t = liste (1/vertex) → variable

INPUTS Grasshopper :
  - mesh : Mesh (Item)       - Surface moyenne voûte
  - t : float/list (Item)    - Épaisseur [m]
  - rho : float (Item)       - Masse volumique [kg/m³]
  - g : float (Item)         - Gravité [9.81]
  - support_pts : Point3d (List)  - Points appui [opt.]
  - support_crv : Curve (Item)    - Courbe appui [opt.]
  - support_idx : int (List)      - Indices vertex [opt.]
  - support_tol : float (Item)    - Tolérance [0.1]
  - flip : bool (Item)            - Inverser normales [False]

OUTPUTS Grasshopper :
  - nodes : str - JSON nœuds (compatible Module 2)
  - branches : str - JSON branches (compatible Module 2)
  - n_int : int - Nombre nœuds intérieurs
  - nodes_int : Point3d list - Nœuds intérieurs (visu)
  - nodes_bou : Point3d list - Nœuds boundary (visu)
  - form_diagram : Line list - Projection XY
  - joint_lines : Line list - Joints [Pi, Pe] (visu)
  - info : str - Résumé modèle

Sources :
  MANUAL §2.1.2-2.1.3 : création du modèle
  MANUAL Table 2.1, 2.2 : colonnes N et B
  MF1 §3.3 eq.3.33 : CSG = h'/(2|e|)
  c2_Blocks.bas : centres de gravité et poids
"""

import Rhino.Geometry as rg
import Rhino
import math
import json

# ====================== DÉFAUTS ======================
if 't' not in dir() or t is None:
    t = 0.3
if 'rho' not in dir() or rho is None:
    rho = 2000.0
if 'g' not in dir() or g is None:
    g = 9.81
if 'support_tol' not in dir() or support_tol is None:
    support_tol = 0.1
if 'flip' not in dir() or flip is None:
    flip = False


# ====================== UTILITAIRES ======================

def to_mesh(obj):
    """Convertit input en Mesh Rhino."""
    if obj is None:
        return None
    if isinstance(obj, rg.Mesh):
        return obj
    try:
        doc_obj = Rhino.RhinoDoc.ActiveDoc.Objects.Find(obj)
        if doc_obj and isinstance(doc_obj.Geometry, rg.Mesh):
            return doc_obj.Geometry.DuplicateMesh()
    except:
        pass
    return None


def mesh_edges(m):
    """Arêtes du mesh → {(v1,v2): [face_indices]}."""
    edges = {}
    for fi in range(m.Faces.Count):
        f = m.Faces[fi]
        vs = [f.A, f.B, f.C] + ([f.D] if f.IsQuad else [])
        for i in range(len(vs)):
            e = tuple(sorted([vs[i], vs[(i+1) % len(vs)]]))
            edges.setdefault(e, []).append(fi)
    return edges


def vertex_area(m, vi):
    """Aire tributaire = Σ (aire_face / nb_vertex_face)."""
    area = 0.0
    for fi in range(m.Faces.Count):
        f = m.Faces[fi]
        vs = [f.A, f.B, f.C] + ([f.D] if f.IsQuad else [])
        if vi not in vs:
            continue
        pts = [m.Vertices[v] for v in vs]
        fa = 0.0
        for i in range(1, len(pts) - 1):
            v1 = rg.Vector3d(pts[i].X - pts[0].X, pts[i].Y - pts[0].Y, pts[i].Z - pts[0].Z)
            v2 = rg.Vector3d(pts[i+1].X - pts[0].X, pts[i+1].Y - pts[0].Y, pts[i+1].Z - pts[0].Z)
            fa += rg.Vector3d.CrossProduct(v1, v2).Length / 2
        area += fa / len(vs)
    return area


def is_naked(vi, edges):
    """Vertex sur bord libre (arête adjacente à 1 seule face)."""
    for (a, b), faces in edges.items():
        if vi in (a, b) and len(faces) == 1:
            return True
    return False


# ====================== APPUIS ======================

def find_supports(m, edges, pts_in=None, crv_in=None, idx_in=None, tol=0.1):
    """
    Priorité : idx > pts > crv > Z auto.
    Retourne (list_vertex_indices, méthode_str).
    """
    nv = m.Vertices.Count

    # 1. Indices directs
    if idx_in is not None:
        try:
            il = [int(i) for i in list(idx_in)]
            if il:
                return il, "indices directs ({})".format(len(il))
        except:
            pass

    # 2. Points
    if pts_in is not None:
        try:
            pl = list(pts_in)
            if pl:
                sup = []
                for vi in range(nv):
                    v = m.Vertices[vi]
                    for p in pl:
                        if math.sqrt((v.X-p.X)**2+(v.Y-p.Y)**2+(v.Z-p.Z)**2) < tol:
                            sup.append(vi)
                            break
                return sup, "points ({} pts)".format(len(pl))
        except:
            pass

    # 3. Courbe
    if crv_in is not None:
        try:
            sup = []
            for vi in range(nv):
                pt = rg.Point3d(m.Vertices[vi].X, m.Vertices[vi].Y, m.Vertices[vi].Z)
                ok, param = crv_in.ClosestPoint(pt)
                if ok and pt.DistanceTo(crv_in.PointAt(param)) < tol:
                    sup.append(vi)
            return sup, "courbe"
        except:
            pass

    # 4. Z auto
    zmin = min(m.Vertices[i].Z for i in range(nv))
    return [i for i in range(nv) if m.Vertices[i].Z <= zmin + tol], \
           "Z auto (zmin={:.2f})".format(zmin)


# ====================== ÉPAISSEUR ======================

def thickness_at(vi, t_in, nv):
    """Épaisseur au vertex vi (uniforme ou variable)."""
    if isinstance(t_in, (int, float)):
        return float(t_in)
    try:
        tl = list(t_in)
        if len(tl) >= nv:
            return float(tl[vi])
        return float(tl[0]) if tl else 0.3
    except:
        return float(t_in)


# ====================== CONSTRUCTION ======================

def build_model(m, supports, t_in, rho_v, g_v, do_flip):
    """
    Mesh + paramètres physiques → modèle TNA complet.

    Nœuds intérieurs :
      Weight = rho × g × aire_tributaire × t
      Bzi = ZG - t/2,  Bze = ZG + t/2  (projection verticale)

    Branches :
      Joint : milieu_arête ± t/2 × normale_moyenne
      Surface ≈ L_arête × t

    Sources : MANUAL §2.1.3, Table 2.1, Table 2.2
    """
    if do_flip:
        m.Flip(True, True, True)
    m.Normals.ComputeNormals()

    edges = mesh_edges(m)
    nv = m.Vertices.Count

    interior, boundary = [], []
    for vi in range(nv):
        if vi in supports or is_naked(vi, edges):
            boundary.append(vi)
        else:
            interior.append(vi)

    # Mapping vertex → index nœud (int d'abord, 1-based)
    v2n = {}
    idx = 1
    for vi in interior:
        v2n[vi] = idx; idx += 1
    ni = len(interior)
    for vi in boundary:
        v2n[vi] = idx; idx += 1
    nb = len(boundary)

    def unit_normal(vi):
        """Normale unitaire au vertex, orientée z>0."""
        n = m.Normals[vi]
        nl = math.sqrt(n.X**2 + n.Y**2 + n.Z**2)
        if nl < 1e-12:
            return 0.0, 0.0, 1.0
        nx, ny, nz = n.X/nl, n.Y/nl, n.Z/nl
        if nz < 0:
            nx, ny, nz = -nx, -ny, -nz
        return nx, ny, nz

    # ---- NOEUDS ----
    nodes = []
    for vi in interior:
        v = m.Vertices[vi]
        nx, ny, nz = unit_normal(vi)
        ti = thickness_at(vi, t_in, nv)
        area = vertex_area(m, vi)
        w = rho_v * g_v * area * ti

        nodes.append({
            'Index': v2n[vi], 'VertexIndex': vi,
            'XG': float(v.X), 'YG': float(v.Y), 'ZG': float(v.Z),
            'nx': nx, 'ny': ny, 'nz': nz,
            'Area': area, 'Boundary': False,
            'Weight': w, 'Thickness': ti,
            'Bzi': v.Z - ti/2.0,  'Bze': v.Z + ti/2.0,
            'BXi': v.X - ti/2*nx, 'BYi': v.Y - ti/2*ny, 'BZi': v.Z - ti/2*nz,
            'BXe': v.X + ti/2*nx, 'BYe': v.Y + ti/2*ny, 'BZe': v.Z + ti/2*nz,
        })

    for vi in boundary:
        v = m.Vertices[vi]
        nx, ny, nz = unit_normal(vi)
        ti = thickness_at(vi, t_in, nv)

        nodes.append({
            'Index': v2n[vi], 'VertexIndex': vi,
            'XG': float(v.X), 'YG': float(v.Y), 'ZG': float(v.Z),
            'nx': nx, 'ny': ny, 'nz': nz,
            'Area': vertex_area(m, vi), 'Boundary': True,
            'Weight': 0.0, 'Thickness': ti,
            'Bzi': v.Z - ti/2.0, 'Bze': v.Z + ti/2.0,
            'BXi': v.X - ti/2*nx, 'BYi': v.Y - ti/2*ny, 'BZi': v.Z - ti/2*nz,
            'BXe': v.X + ti/2*nx, 'BYe': v.Y + ti/2*ny, 'BZe': v.Z + ti/2*nz,
        })

    # ---- BRANCHES ----
    branches = []
    bidx = 1
    for (vi1, vi2) in edges.keys():
        vh, vt = m.Vertices[vi1], m.Vertices[vi2]
        dx, dy, dz = vh.X - vt.X, vh.Y - vt.Y, vh.Z - vt.Z
        L = math.sqrt(dx**2 + dy**2 + dz**2)
        LH = math.sqrt(dx**2 + dy**2)

        # Milieu arête
        mx = (vh.X + vt.X) / 2.0
        my = (vh.Y + vt.Y) / 2.0
        mz = (vh.Z + vt.Z) / 2.0

        # Normale moyenne au joint
        n1x, n1y, n1z = unit_normal(vi1)
        n2x, n2y, n2z = unit_normal(vi2)
        nmx, nmy, nmz = (n1x+n2x)/2, (n1y+n2y)/2, (n1z+n2z)/2
        nml = math.sqrt(nmx**2 + nmy**2 + nmz**2)
        if nml < 1e-12:
            nml = 1.0
        nmx /= nml; nmy /= nml; nmz /= nml

        # Épaisseur au joint
        t1 = thickness_at(vi1, t_in, nv)
        t2 = thickness_at(vi2, t_in, nv)
        tm = (t1 + t2) / 2.0

        branches.append({
            'Index': bidx,
            'Head_node': v2n[vi1], 'Tail_node': v2n[vi2],
            'u': dx, 'v': dy, 'w': dz,
            'L': L, 'LH': max(LH, 1e-12),
            'xi': mx - tm/2*nmx, 'yi': my - tm/2*nmy, 'zi': mz - tm/2*nmz,
            'xe': mx + tm/2*nmx, 'ye': my + tm/2*nmy, 'ze': mz + tm/2*nmz,
            'Surface': L * tm,
        })
        bidx += 1

    return {'nodes': nodes, 'branches': branches, 'n_int': ni, 'n_bou': nb, 'mesh': m}


# ====================== MAIN ======================

m_in = to_mesh(mesh) if 'mesh' in dir() else None

if m_in is None:
    print("ERREUR: Mesh invalide ou non connecte.")
    print("  Branchez un Mesh Rhino sur 'mesh' (Type Hint = Mesh)")
    nodes = branches = "[]"
    n_int = 0
    nodes_int = nodes_bou = form_diagram = joint_lines = []
    info = "ERREUR: pas de mesh"
else:
    print("=" * 50)
    print("MaNACoH v4 - Module 0: MESH -> MODELE (v4)")
    print("=" * 50)

    edg = mesh_edges(m_in)
    has_idx = 'support_idx' in dir() and support_idx is not None
    has_pts = 'support_pts' in dir() and support_pts is not None
    has_crv = 'support_crv' in dir() and support_crv is not None

    supports, method = find_supports(
        m_in, edg,
        pts_in=support_pts if has_pts else None,
        crv_in=support_crv if has_crv else None,
        idx_in=support_idx if has_idx else None,
        tol=support_tol)

    model = build_model(m_in, supports, t, rho, g, flip)
    n_int = model['n_int']
    n_bou = model['n_bou']

    # Sérialiser
    nodes = json.dumps(model['nodes'])
    branches = json.dumps(model['branches'])

    # Visuels
    nmap = {nd['Index']: nd for nd in model['nodes']}
    nodes_int = [rg.Point3d(nd['XG'], nd['YG'], nd['ZG'])
                 for nd in model['nodes'] if not nd['Boundary']]
    nodes_bou = [rg.Point3d(nd['XG'], nd['YG'], nd['ZG'])
                 for nd in model['nodes'] if nd['Boundary']]
    form_diagram = [rg.Line(
        rg.Point3d(nmap[b['Head_node']]['XG'], nmap[b['Head_node']]['YG'], 0),
        rg.Point3d(nmap[b['Tail_node']]['XG'], nmap[b['Tail_node']]['YG'], 0)
    ) for b in model['branches']]
    joint_lines = [rg.Line(
        rg.Point3d(b['xi'], b['yi'], b['zi']),
        rg.Point3d(b['xe'], b['ye'], b['ze'])
    ) for b in model['branches']]

    W_total = sum(nd['Weight'] for nd in model['nodes'] if not nd['Boundary'])
    t_str = "{:.3f}m".format(t) if isinstance(t, (int, float)) else "variable"
    info_lines = [
        "Vertices: {}".format(m_in.Vertices.Count),
        "Noeuds: {} int + {} bou".format(n_int, n_bou),
        "Branches: {}".format(len(model['branches'])),
        "Appuis ({}): {}".format(method, len(supports)),
        "Epaisseur: {}".format(t_str),
        "Densite: {} kg/m3".format(rho),
        "Poids total: {:.1f} N ({:.1f} kN)".format(W_total, W_total/1000),
    ]
    info = "\n".join(info_lines)
    for l in info_lines:
        print("  " + l)
