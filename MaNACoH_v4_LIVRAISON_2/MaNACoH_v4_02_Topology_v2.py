#! python3
"""
MaNACoH v4 - Module 2 : TOPOLOGY
================================
INPUTS (Type Hints Grasshopper):
  - nodes : str (Item) - JSON du Module 1
  - branches : str (Item) - JSON du Module 1
  - n_int : int (Item) - Du Module 1

OUTPUTS:
  - Ci : list - Matrice branch-node (intérieurs) [n_branches x n_int]
  - Cb : list - Matrice branch-node (boundary) [n_branches x n_bou]
  - u_vec : list - Composantes u des branches
  - v_vec : list - Composantes v des branches  
  - w_vec : list - Composantes w des branches
  - L_vec : list - Longueurs des branches
  - LH_vec : list - Longueurs horizontales
  - nodes_data : list - Noeuds désérialisés (pour modules suivants)
  - branches_data : list - Branches désérialisées (pour modules suivants)
"""

import json

# === FONCTIONS ===

def parse_json_input(data):
    """Parse un input JSON string ou retourne la liste si déjà parsée"""
    if data is None:
        return None
    
    # Déjà une liste de dicts
    if isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
        return data
    
    # String JSON
    if isinstance(data, str):
        try:
            parsed = json.loads(data)
            if isinstance(parsed, list):
                return parsed
        except Exception as e:
            print("Erreur JSON: {}".format(e))
            return None
    
    # Autre cas: essayer de convertir
    try:
        # Peut-être une liste Grasshopper
        result = list(data)
        if len(result) > 0 and isinstance(result[0], dict):
            return result
    except:
        pass
    
    return None

# === MAIN ===

print("="*50)
print("MaNACoH v4 - Module 2: TOPOLOGY")
print("="*50)

# Debug: voir ce qu'on reçoit
if 'nodes' in dir():
    print("Type nodes: {}".format(type(nodes)))
    if isinstance(nodes, str):
        print("nodes (str): {}...".format(nodes[:100] if len(nodes) > 100 else nodes))
if 'branches' in dir():
    print("Type branches: {}".format(type(branches)))

# Parser les inputs
nodes_data = parse_json_input(nodes) if 'nodes' in dir() else None
branches_data = parse_json_input(branches) if 'branches' in dir() else None

has_n_int = 'n_int' in dir() and n_int is not None and n_int > 0

if nodes_data is None:
    print("ERREUR: nodes non valide ou non connecte")
if branches_data is None:
    print("ERREUR: branches non valide ou non connecte")
if not has_n_int:
    print("ERREUR: n_int non valide ou non connecte")

has_data = nodes_data is not None and branches_data is not None and has_n_int

if not has_data:
    Ci = Cb = []
    u_vec = v_vec = w_vec = L_vec = LH_vec = []
    nodes_data = branches_data = []
else:
    print("Noeuds: {} | Branches: {}".format(len(nodes_data), len(branches_data)))
    
    nb = len(branches_data)
    nn = len(nodes_data)
    n_bou = nn - n_int
    
    # Mapping Index -> colonne (E3 fix : robuste, pas de supposition sur Index)
    # Source : MANUAL §2.1 "interior nodes have lower indices than boundary"
    col_int = {}
    col_bou = {}
    
    int_count = 0
    bou_count = 0
    for n in nodes_data:
        if not n['Boundary']:
            col_int[n['Index']] = int_count
            int_count += 1
        else:
            col_bou[n['Index']] = bou_count
            bou_count += 1
    
    # Vérification de cohérence
    if int_count != n_int:
        print("ATTENTION: n_int={} mais {} noeuds interieurs trouves".format(n_int, int_count))
        n_int = int_count
    
    # Matrices
    Ci = [[0.0] * n_int for _ in range(nb)]
    Cb = [[0.0] * max(n_bou, 1) for _ in range(nb)]
    
    # Vecteurs
    u_vec = []
    v_vec = []
    w_vec = []
    L_vec = []
    LH_vec = []
    
    for j, b in enumerate(branches_data):
        head = b['Head_node']
        tail = b['Tail_node']
        
        if head in col_int:
            Ci[j][col_int[head]] = 1.0
        if tail in col_int:
            Ci[j][col_int[tail]] = -1.0
        
        if n_bou > 0:
            if head in col_bou:
                Cb[j][col_bou[head]] = 1.0
            if tail in col_bou:
                Cb[j][col_bou[tail]] = -1.0
        
        u_vec.append(b['u'])
        v_vec.append(b['v'])
        w_vec.append(b.get('w', 0.0))
        L_vec.append(b['L'])
        LH_vec.append(b['LH'])
    
    print("Ci: {} x {} | Cb: {} x {}".format(nb, n_int, nb, n_bou))
