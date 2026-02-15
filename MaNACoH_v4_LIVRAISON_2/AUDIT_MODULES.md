# MaNACoH v4 — AUDIT FINAL DES MODULES
# =======================================
# Date : 2026-02-14
# Auditeur : Claude

## SYNTHÈSE

| Module | Fichier | Lignes | Bugs | Gravité max |
|--------|---------|--------|------|-------------|
| 0 v4 — Pattern+Physique | `_00_Pattern_v4.py` | 396 | 0 | — |
| 1 v1 — Géométrie param. | `_01_Geometry_v1.py` | 405 | 3 | **HAUTE** |
| 2 v2 — Topologie | `_02_Topology_v2.py` | 149 | 1 | MOYENNE |
| 3 v53 — TNA Optimiseur | `_03_TNA_v53.py` | 1342 | 1 | BASSE |
| 4 v2 — Safety | `_04_Safety_Rigorous_v2.py` | 599 | 0 | — |
| 5 v4 — Réactions | `_05_Reactions_v4.py` | 475 | 0 | — |


## CHECKLIST PAR MODULE

### ═══ MODULE 0 v4 : PATTERN + MODÈLE PHYSIQUE ═══

**Rôle** : Mesh Rhino → nœuds/branches avec poids, joints, bornes CSG
**Environnement** : Grasshopper (IronPython, Rhino.Geometry)

**INPUTS (GH)**
| Input | Type | Défaut | Vérifié |
|-------|------|--------|---------|
| mesh | Mesh (Item) | — | ✅ Validation to_mesh() |
| t | float/list (Item) | 0.3 | ✅ thickness_at() gère scalaire/liste |
| rho | float (Item) | 2000 | ✅ |
| g | float (Item) | 9.81 | ✅ |
| support_pts | Point3d (List) | None | ✅ 4 modes de détection |
| support_crv | Curve (Item) | None | ✅ |
| support_idx | int (List) | None | ✅ priorité max |
| support_tol | float (Item) | 0.1 | ✅ |
| flip | bool (Item) | False | ✅ |

**OUTPUTS (GH)**
| Output | Type | Vérifié |
|--------|------|---------|
| nodes | str (JSON) | ✅ Toutes les clés pour Module 2/3 |
| branches | str (JSON) | ✅ Toutes les clés pour Module 2/3 |
| n_int | int | ✅ |
| nodes_int | Point3d list | ✅ |
| nodes_bou | Point3d list | ✅ |
| form_diagram | Line list | ✅ Projection XY |
| joint_lines | Line list | ✅ **NOUVEAU** (intrados→extrados) |
| info | str | ✅ Résumé lisible |

**CONVENTION u** : u = x_head - x_tail ✅ (L.302)
**POIDS** : Weight = rho × g × aire_trib × t ✅ (L.261)
**BORNES** : Bzi = ZG - t/2, Bze = ZG + t/2 ✅ (L.264)
**JOINTS** : milieu_arête ± t/2 × n̂ ✅ (L.311-312)
**NORMALE** : orientée nz > 0 ✅ (L.241)
**ERREUR** : message clair si mesh absent ✅

**Verdict : ✅ OK — Aucun bug détecté**

---

### ═══ MODULE 1 v1 : GÉOMÉTRIE PARAMÉTRIQUE ═══

**Rôle** : Générer géométries analytiques (arcs, dômes) pour benchmarks
**Environnement** : Python pur (pas de Rhino)

#### BUG B1.1 — HAUTE : Convention u INVERSÉE
**Fichier** : L.306-308 dans `_make_branch()`
```python
# ACTUEL (FAUX) :
dx = nt['XG'] - nh['XG']  # tail - head = INVERSÉ
dy = nt['YG'] - nh['YG']

# CORRECT :
dx = nh['XG'] - nt['XG']  # head - tail (convention VBA MANUAL)
dy = nh['YG'] - nt['YG']
```
**Impact** : u est de signe opposé → dans Module 2, Ci·u inverse → K inversé
→ les ζ changent de signe → impact sur le CSG et les forces.
**Pourquoi ça marchait** : en 1D (arc symétrique), le changement de signe
de u est compensé par le changement de signe de ζ → CSG identique.
Mais pour un réseau 2D/3D asymétrique, ça casserait.

#### BUG B1.2 — BASSE : `create_preset()` mute PRESETS via `pop()`
**Fichier** : L.374-375
```python
desc = p.pop('description', '')     # MUTE le dict global !
expected = p.pop('expected_csg', None)
```
**Impact** : Fonctionne car `p.pop()` + `p['description'] = desc` L.378-379
remet la valeur. Mais thread-unsafe et malpropre.
**Fix** : Utiliser `p.get()` + copie locale.

#### BUG B1.3 — BASSE : Docstring annonce wolfe/reese non implémentés
**Fichier** : L.22-23
```
"wolfe" : Fig 4.22, 18 voussoirs, dome R=10
"reese" : Fig 4.24, dome R=10, 10 meridiens
```
**Impact** : KeyError si l'utilisateur tente `create_preset('wolfe')`.
**Fix** : Retirer de la docstring ou implémenter.

#### INCOMPLET : Pas de géométrie sphérique
Le Module 1 n'implémente que l'arc cylindrique (pas le dôme sphérique).
La fonction `_cg_spherical_voussoir()` existe mais n'est appelée nulle part.

**Verdict : ⚠️ 3 bugs (1 HAUTE, 2 BASSES) + incomplet dôme**

---

### ═══ MODULE 2 v2 : TOPOLOGIE ═══

**Rôle** : Parser JSON → matrices Ci, Cb, vecteurs u, v, LH
**Environnement** : Grasshopper

#### BUG B2.1 — MOYENNE : Nom d'output `nodes_data` ≠ attendu `nodes`
**Fichier** : L.71, 88
Le Module 3 MAIN attend `nodes` et `branches` (pas `nodes_data`/`branches_data`).
En Grasshopper c'est géré par les noms d'output du composant, mais c'est
source de confusion si on câble mal.

**Recommandation** : Ajouter des alias ou renommer les outputs.

**Autres vérifications** :
- Convention Ci[j, head]=+1, Ci[j, tail]=-1 : ✅ (L.132-134)
- u provient de b['u'] directement : ✅ (L.142)
- Gestion n_int incohérent : ✅ Warning + autocorrection (L.112-114)
- Erreur messages : ✅ Clairs

**Verdict : ⚠️ 1 bug mineur (nommage outputs)**

---

### ═══ MODULE 3 v53 : TNA + OPTIMISEUR ═══

**Rôle** : Équilibre horizontal + vertical + optimisation DOF
**Environnement** : Python pur + optionnel Rhino

#### BUG B3.1 — BASSE : Docstring dit "v50" mais c'est v53
**Fichier** : L.3
```
Module 3 : EQUILIBRIUM TNA (v50 - Refonte théorique)
```
Devrait être v53.

**FONCTIONS CRITIQUES VÉRIFIÉES** :
| Fonction | Rôle | Vérifié |
|----------|------|---------|
| build_K_matrix | K = [CiᵀUL⁻¹; CiᵀVL⁻¹] | ✅ MF1 eq.4.27 |
| remove_redundant_rows | Nettoyage lignes dépendantes | ✅ |
| identify_dof_branches | Identification DOF | ✅ Gram-det |
| solve_horizontal | ζ_unk = -K⁻¹·K_dof·ζ_dof | ✅ MF1 eq.4.30 |
| update_zeta_fast | Reconstruction rapide ζ | ✅ |
| build_Di_Db | Di, Db depuis q | ✅ MF1 eq.4.36-37 |
| solve_vertical | zi = Di⁻¹(-pz-Db·zb) | ✅ MF1 eq.4.38 |
| compute_csg_simplified | CSG = h'/(2|e|) | ✅ MF1 eq.3.33 |
| optimize_dof | Phase 0/1/2 | ✅ Testé multi-DOF |

**CONVENTION pz** : pz = -Weight (négatif) ✅ (L.1252)
**CONVENTION ζ** : ζ > 0 = compression ✅
**FALLBACK lstsq** : ✅ Pour Kunk rectangulaire (B.E7)
**Phase 0 multi-DOF** : scan log + proportionnel LH ✅
**Phase 2 CD** : bornes élargies [0.05×, 20×] ✅
**DIRECT** : Jones 1993 simplifié ✅

**MAIN GH** :
| Variable attendue | Source | Vérifié |
|-------------------|--------|---------|
| nodes | Module 0/1 | ✅ (list of dict) |
| branches | Module 0/1 | ✅ |
| Ci, Cb | Module 2 | ✅ |
| u_vec, v_vec | Module 2 | ✅ |
| LH_vec | Module 2 | ✅ |
| n_int | Module 0/2 | ✅ |
| optimize | bool | ✅ défaut True |
| max_iter | int | ✅ défaut 200 |
| objective | str | ✅ "proximity"/"max_csg" |
| boundary_position | str | ✅ "middle"/"intrados"/"extrados" |

**OUTPUTS GH** :
| Output | Type | Vérifié |
|--------|------|---------|
| Z | list float | ✅ ni + nb altitudes |
| q | list float | ✅ m densités de force |
| F | list float | ✅ m forces |
| LH_dual | list float | ✅ m forces horizontales |
| CSG_min | float | ✅ |
| thrust_pts | Point3d/tuple | ✅ (conditionnel Rhino) |
| thrust_lines | Line list | ✅ (conditionnel Rhino) |

**Verdict : ✅ OK — 1 bug cosmétique (version docstring)**

---

### ═══ MODULE 4 v2 : SAFETY ═══

**Rôle** : CSG, CSS, glissement, contraintes par joint
**Environnement** : Grasshopper (Rhino + System.Drawing)

**FORMULES VÉRIFIÉES** :
- CSG = h'/(2|e|) : ✅ MF1 eq.3.33
- CSS = (1-2|e|/h')/n' : ✅ Fantin eq.3.22
- σ_max = 2N'/(b'·h'_comprimé) : ✅
- CS_glissement = μ·N'/|V'| : ✅ Coulomb

**INPUTS** : nodes_data, branches_data, Z, LH_dual, sigma_c, mu ✅
**OUTPUTS** : CSG, CSS, CS_glissement, N/M/V/e/σ_max, couleurs ✅

**Note** : Requiert Rhino (System.Drawing pour les couleurs).
Pas de fallback Python pur.

**Verdict : ✅ OK**

---

### ═══ MODULE 5 v4 : RÉACTIONS ═══

**Rôle** : Réactions d'appui par équilibre nodal
**Environnement** : Grasshopper (Rhino)

**INPUTS** : nodes, branches, Z, q, vector_scale ✅
**OUTPUTS** : reactions, total_RV, total_weight, equilibrium_ok,
              vecteurs visuels (R, RH, RV), TextDots ✅

**VÉRIFICATION GLOBALE** : Σ RV = Σ Poids ✅

**Verdict : ✅ OK**

---

## MATRICE DE COMPATIBILITÉ ENTRÉES/SORTIES

```
Module 0 v4 ──→ Module 2 v2 ──→ Module 3 v53 ──→ Module 4 v2
  nodes(JSON)      Ci,Cb           Z,q,F            CSG,CSS
  branches(JSON)   u,v,LH_vec     LH_dual          N,M,V,e
  n_int            nodes_data*     CSG_min          couleurs
                   branches_data*  thrust_pts       ──→ Module 5 v4
                                   thrust_lines         reactions
                                                        équilibre

Module 1 v1 ──→ (même interface que Module 0, sans JSON)
  nodes(list)
  branches(list)
```

*Nommage : Module 2 produit `nodes_data` mais Module 3 attend `nodes`.
En Grasshopper, c'est le nom de l'output du composant qui compte,
donc il suffit de nommer correctement les pins.

## BUGS À CORRIGER

| ID | Module | Gravité | Description | Fix |
|----|--------|---------|-------------|-----|
| B1.1 | Mod 1 | HAUTE | u = tail-head au lieu de head-tail | Inverser dx, dy dans _make_branch |
| B1.2 | Mod 1 | BASSE | create_preset mute PRESETS via pop | Utiliser get + copie |
| B1.3 | Mod 1 | BASSE | Docstring annonce wolfe/reese absents | Retirer ou implémenter |
| B2.1 | Mod 2 | MOYENNE | Output nodes_data ≠ nodes attendu | Renommer ou documenter |
| B3.1 | Mod 3 | COSM. | Docstring dit v50 au lieu de v53 | Mettre à jour |

## RECOMMANDATIONS

1. **Corriger B1.1 en priorité** (inversion u dans Module 1)
2. Le Module 1 devrait être utilisé UNIQUEMENT pour les benchmarks
   analytiques. Le workflow production passe par Module 0 v4 (mesh Rhino).
3. Nettoyer les anciennes versions (v51, v52 du Module 3) de /outputs.
4. Ajouter un test d'intégration automatique Module 0→2→3→4→5.
