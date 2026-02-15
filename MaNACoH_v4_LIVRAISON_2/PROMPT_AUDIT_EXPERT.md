# PROMPT AUDIT + REFONTE — MaNACoH v4

## Contexte

Tu audites et complètes le pipeline Python **MaNACoH v4** — un portage en Python/Grasshopper du logiciel Excel/VBA MaNACoH v2.0.5 d'analyse de stabilité des voûtes en maçonnerie par la méthode des réseaux de forces (Thrust Network Analysis).

Le code VBA original est dans les fichiers projet (.bas). La documentation théorique est dans les PDF projet (MF1.pdf à MF4.pdf = thèse Fantin 2017, MANUAL = MaNACoH_v2_0_5_MANUAL.pdf, SETRA 1982 = ponts_en_maconnerie...).

## OBJECTIF PRIORITAIRE : Module 1 universel

Le Module 1 actuel ne génère que des **arcs cylindriques 1D**. Il doit devenir un **générateur paramétrique universel** comme la fonction VBA `a3_CreateRegularGeo` dans `c1_CreateModel.bas`.

### Architecture VBA à reproduire

La fonction `a3_CreateRegularGeo(TypeOfGeo, n, m, u0, u1, v0, v1, ...)` :
- Crée une grille **(n+2) × (m+2)** de nœuds : i=1..n+2 (le long de u), j=1..m+2 (le long de v)
- i=1 et i=n+2 sont boundary (u0, u1), j=1 et j=m+2 sont boundary (v0, v1)
- Nœuds intérieurs : i=2..n+1, j=2..m+1 → total n×m nœuds intérieurs
- Index nœud : p = (j-1)×(n+2) + i
- **Branches méridien** (le long de u/phi) : head=p, tail=p+1
- **Branches parallèle/hoop** (le long de v/theta) : head=p, tail=p+(n+2)
- Chaque type a sa propre fonction `BlockCarac[Type]` pour le CG et le poids

### Types de géométrie (c1_CreateModelParameters.bas)

| Type | Paramètres u,v | Formules CG | Fichier VBA |
|------|----------------|-------------|-------------|
| **Cyl** | u=phi, v=y | `BlockCaracCyl` c2_Blocks L.1038 | c1_CreateModelParameters L.5 |
| **Sphe** | u=phi, v=theta | `BlockCaracSphe` c2_Blocks L.253 | c1_CreateModelParameters L.190 |
| **Oval** | u=phi, v=theta | `BlockCaracOval` c2_Blocks | c1_CreateModelParameters L.258 |
| **Pav** | u=phi, v=theta | `BlockCaracPav` c2_Blocks | c1_CreateModelParameters L.137 |
| **Fla** | u=x, v=y | `BlockCaracFla` c2_Blocks | c1_CreateModelParameters L.338 |
| **Hel** | u=r, v=theta | `BlockCaracHel` c2_Blocks | c1_CreateModelParameters L.81 |

### Formules clés BlockCaracSphe (c2_Blocks.bas L.253-360)

```
Volume = rho/3 × (re³ - ri³) × ΔΘ × |cos(φ₀) - cos(φ₁)|

rG·sin(φG) = 3/8 × (re⁴-ri⁴)/(re³-ri³) × sin(ΔΘ/2)/ΔΘ ×
             [(φ₁-φ₀)×2 - sin(2φ₁) + sin(2φ₀)] / [-cos(φ₁) + cos(φ₀)]

rG·cos(φG) = 3/16 × (re⁴-ri⁴)/(re³-ri³) ×
             [-cos(2φ₁) + cos(2φ₀)] / [-cos(φ₁) + cos(φ₀)]

XG = rG·sin(φG)·cos(θG),  YG = rG·sin(φG)·sin(θG),  ZG = rG·cos(φG)

Bzi = ri·cos(φG),  Bze = re·cos(φG)  (projection radiale)
zi_vertical = ri·cos(asin(rG·sin(φG)/ri))  (projection verticale sur intrados)
```

### Conditions aux limites (VBA c1_CreateModel.bas L.640-700)

Les boundary sont contrôlés par des flags :
- `u0Boundary_BranchesExist` : branches le long de u en i=1 (phi_min boundary)
- `u1Boundary_BranchesExist` : branches le long de u en i=n+2 (phi_max boundary)  
- `v0Boundary_BranchesExist` : branches le long de v en j=1 (theta_min boundary) = **hoop branches**
- `v1Boundary_BranchesExist` : branches le long de v en j=m+2 (theta_max boundary)
- `MeridianBoundaryNodesOnExtrados` : nœuds boundary sur l'extrados

Pour un dôme sphérique typique (c1_CreateModelParameters L.230) :
- HoopBranchesExist = True
- n=30 voussoirs méridien, m=1 voussoir parallèle (lune complète)
- phi: 0→90°, theta: -0.5°→+0.5° (lune mince)

### Ce qui manque dans le Module 1 actuel

1. ❌ Pas de grille 2D (n×m) — seulement une chaîne 1D
2. ❌ Pas de branches parallèles/hoop
3. ❌ Pas de `BlockCaracSphe` — seulement cylindrique
4. ❌ Pas de gestion des 4 bords boundary (u0, u1, v0, v1)
5. ❌ Pas de joints méridien et parallèle séparés

### Outputs requis (compatibles Module 3)

Le Module 3 attend `nodes` et `branches` (JSON ou list of dict) avec :

**Nœuds** : Index, XG, YG, ZG, Boundary, Weight, Thickness, Bzi, Bze, BXi, BYi, BZi, BXe, BYe, BZe

**Branches** : Index, Head_node, Tail_node, u, v, LH, L, xi, yi, zi, xe, ye, ze

### Presets de validation (MANUAL Chapter 4)

| Preset | Type | n | m | Paramètres | CSG attendu |
|--------|------|---|---|------------|-------------|
| heyman_15 | Cyl | 15 | 1 | R=10, t/R=0.15 | 1.3952 |
| heyman_tmin | Cyl | 15 | 1 | R=10, t/R=0.1075 | ~1.0 |
| wolfe | Sphe | 18 | ? | R=10 | MANUAL Fig 4.22 |
| levy | Sphe | 30 | 1 | R=10, t/R=0.042 | MANUAL §4.3 |

## Les 6 modules à auditer

(suite identique au prompt précédent...)

| Module | Fichier | Rôle |
|--------|---------|------|
| 0 v4 | `MaNACoH_v4_00_Pattern_v4.py` | Mesh Rhino → nœuds/branches + poids, joints, bornes CSG |
| 1 v1 | `MaNACoH_v4_01_Geometry_v1.py` | Géométries paramétriques analytiques (arcs, presets Heyman) |
| 2 v2 | `MaNACoH_v4_02_Topology_v2.py` | JSON → matrices Ci, Cb, vecteurs u, v, LH |
| 3 v53 | `MaNACoH_v4_03_TNA_v53.py` | Équilibre horizontal + vertical + optimisation DOF |
| 4 v2 | `MaNACoH_v4_04_Safety_Rigorous_v2.py` | CSG, CSS, glissement par joint |
| 5 v4 | `MaNACoH_v4_05_Reactions_v4.py` | Réactions d'appui par équilibre nodal |

## Chaîne de calcul théorique (MF1.pdf §4.3.5)

```
1. K = [Ciᵀ·diag(u/LH) ; Ciᵀ·diag(v/LH)]     → eq. 4.27
2. Identifier DOF/UNK dans K                      → eq. 4.28-4.29
3. ζ_unk = -Kunk⁻¹ · Kdof · ζ_dof                → eq. 4.30
4. q = ζ / LH                                     → densité de force
5. Di = Ciᵀ·diag(q)·Ci,  Db = Ciᵀ·diag(q)·Cb    → eq. 4.36-4.37
6. zi = Di⁻¹·(-pz - Db·zb)                        → eq. 4.38
7. CSG_j = h'_j / (2·|e_j|)                       → MF1 eq. 3.33
8. Optimiser ζ_dof pour max(min(CSG_j))            → f_Optimization.bas
```

## Conventions critiques (sources de bugs passés)

- **u = x_head − x_tail** (MANUAL Table 2.2, VBA convention)
- **Ci[j, head] = +1, Ci[j, tail] = −1** (matrice branche-nœud)
- **pz = −Weight** (négatif = poids vers le bas, convention VBA d_ComputeTypology.bas L.244)
- **ζ > 0 = compression** (densité de force positive = branche comprimée)
- **LH = √(u² + v²)** (longueur horizontale, jamais négative)
- **Orientation arêtes arbitraire** (pas de contrainte Head < Tail, invariance prouvée)

## Checklist d'audit demandée

### A. Fidélité théorique
Pour chaque formule dans le code, vérifier :
1. Correspondance avec l'équation source (MF1, MANUAL, VBA)
2. Conventions de signe respectées
3. Dimensions/unités cohérentes
4. Domaine de validité respecté

### B. Cohérence inter-modules
Pour chaque interface Module N → Module N+1 :
1. Les clés de dict produites par N sont-elles toutes consommées par N+1 ?
2. Les noms d'inputs/outputs GH correspondent-ils ?
3. Les conventions (signe u, signe pz, orientation) sont-elles identiques ?

### C. Robustesse numérique
1. Division par zéro protégée (LH=0 pour branches verticales)
2. Inversion de matrice : cas singulier géré (lstsq fallback)
3. Valeurs extrêmes : ζ très grand/petit, pz=0, t=0
4. Branches verticales : exclues des DOF ?

### D. Conformité VBA
Pour les fonctions critiques, comparer ligne à ligne avec le VBA :
- `e_ComputeEquilibrium.bas` ↔ Module 3 solve_horizontal/solve_vertical
- `f_Optimization.bas` ↔ Module 3 optimize_dof
- `c2_Blocks.bas` ↔ Module 1 _cg_cylindrical_voussoir
- `d_ComputeTypology.bas` ↔ Module 2 construction Ci/Cb
- `g_GraphicalOutput.bas` ↔ Module 4 compute_csg

### E. Points sensibles identifiés (bugs corrigés récemment)
- E1 : Orientation arêtes — était forcée Head<Tail, maintenant arbitraire
- E3 : Mapping Index→colonne — ne plus supposer Index = position
- E5 : u_vec obligatoire — plus de fallback u=LH (perd le signe)
- E7 : Kunk rectangulaire — résolu par lstsq, vérifier que le fallback fonctionne
- E8 : Boundary position — intrados/extrados/middle pour zb
- B1.1 : Convention u dans Module 1 — était inversée (tail−head), corrigée en head−tail

### F. Qualité logicielle
1. Chaque fonction a-t-elle une docstring avec source ?
2. Les inputs GH ont-ils des défauts raisonnables ?
3. Les messages d'erreur sont-ils clairs ?
4. Le code est-il testable hors Rhino ? (Modules 0/4/5 = non, Modules 1/2/3 = oui)

## Format de réponse attendu

Pour chaque module, produire :

```
MODULE X : [nom]
  THÉORIE : ✅/⚠️/❌ [détail par formule]
  INTERFACES : ✅/⚠️/❌ [compatibilité entrées/sorties]
  NUMÉRIQUE : ✅/⚠️/❌ [protections, cas limites]
  VBA : ✅/⚠️/❌ [conformité avec le code original]
  BUGS TROUVÉS : [liste avec gravité HAUTE/MOYENNE/BASSE]
```

Terminer par un tableau de synthèse et les corrections proposées.

## Important

- **NE PAS halluciner** de formules. Si tu ne trouves pas la source dans le corpus, dis-le.
- **Citer les sources** : MF1 eq.X.XX, MANUAL §X.X, VBA fichier L.XXX
- **Vérifier en exécutant** : lancer les tests existants, créer des cas limites.
- Le Module 3 fait 1342 lignes — c'est le cœur, il mérite 60% de l'effort d'audit.
