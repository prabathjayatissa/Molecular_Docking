
# Molecular Docking Practical Guide

## Beginner-Level Practical for Bachelor Students

---

# Aim of the Practical

The aim of this practical is to introduce students to the basics of molecular docking by studying how small molecules (ligands) bind to proteins.

Students will learn how to:

- Obtain protein structures from the Protein Data Bank (PDB)
- Prepare proteins and ligands for docking
- Perform docking using AutoDock Vina/PyRx
- Visualize protein–ligand interactions
- Interpret docking scores and binding modes

---

# Software Required

| Software | Purpose |
|---|---|
| PyRx | Docking platform |
| AutoDock Vina | Docking calculations |
| PyMOL | Structure visualization |
| Discovery Studio Visualizer | Interaction analysis |

---

# Protein–Ligand Systems Used

| Experiment | Protein | Ligand | PDB ID |
|---|---|---|---|
| 1 | Carbonic Anhydrase II | Acetazolamide | 3HS4 |
| 2 | HIV-1 Protease | Amprenavir | 1HPV |
| 3 | COX-2 | Celecoxib | 3LN1 |
| 4 | Estrogen Receptor Alpha | 4-Hydroxytamoxifen | 3ERT |

---

# Theory

## What is Molecular Docking?

Molecular docking is a computational method used to predict how a ligand binds to a protein active site.

---

## Key Terms

| Term | Meaning |
|---|---|
| Protein | Biological target molecule |
| Ligand | Small molecule/drug |
| Active Site | Region where ligand binds |
| Binding Affinity | Strength of interaction |
| Docking Score | Estimated binding energy |

> Lower docking scores generally indicate stronger binding.

---

# General Practical Workflow

1. Download protein structure from PDB
2. Prepare protein structure
3. Extract ligand
4. Load files into PyRx
5. Define docking grid
6. Run docking
7. Visualize interactions
8. Interpret docking score

---

# Experiment 1: Carbonic Anhydrase II – Acetazolamide

**Protein:** Human Carbonic Anhydrase II  
**Ligand:** Acetazolamide  
**PDB ID:** `3HS4`

---

## Step 1: Download Protein Structure

1. Visit the RCSB Protein Data Bank:
   https://www.rcsb.org

2. Search for:
   `3HS4`

3. Download the structure in `.pdb` format.

---

## Step 2: Open Structure in PyMOL

1. Launch PyMOL
2. Open the downloaded `.pdb` file
3. Observe:
   - Protein chains
   - Bound ligand
   - Zinc ion in the active site

---

## Step 3: Protein Preparation

### Remove Water Molecules

Use the PyMOL command:

```python
remove solvent
```

### Remove Existing Ligand

```python
remove resn AZM
```

### Save Prepared Protein

- File → Save Molecule
- Save as:

```text
protein_clean.pdb
```

---

## Step 4: Ligand Preparation

1. Open the original structure again
2. Select ligand:

```python
select ligand, resn AZM
```

3. Export ligand as:

```text
ligand.pdb
```

---

## Step 5: Load Files into PyRx

1. Open PyRx
2. Import:
   - Prepared protein
   - Ligand file

3. Convert both files to `.pdbqt`

---

## Step 6: Set Docking Grid

1. Open Vina Wizard
2. Center grid around original ligand binding site
3. Adjust the grid box to fully cover the active site

---

## Step 7: Run Docking

1. Click **Run Vina**
2. Wait for docking completion
3. Record:
   - Best docking score
   - Number of poses

---

## Step 8: Visualization of Results

Open docked complex in:
- PyMOL
- Discovery Studio Visualizer

Observe:
- Hydrogen bonds
- Ligand orientation
- Zinc interaction

---

# Student Observation Table

| Parameter | Observation |
|---|---|
| Protein Used | |
| Ligand Used | |
| Docking Score (kcal/mol) | |
| Number of Hydrogen Bonds | |
| Important Residues | |

---

# Experiment 2: HIV-1 Protease – Amprenavir

**Protein:** HIV-1 Protease  
**Ligand:** Amprenavir  
**PDB ID:** `1HPV`

---

## Important Notes

- Active site contains catalytic residues:
  - Asp25
  - Asp25'

- Observe symmetric binding pocket

---

## Expected Learning Outcomes

- Protease inhibition
- Drug binding specificity
- Protein–ligand interactions

---

# Experiment 3: COX-2 – Celecoxib

**Protein:** Cyclooxygenase-2 (COX-2)  
**Ligand:** Celecoxib  
**PDB ID:** `3LN1`

---

## Important Notes

- Observe hydrophobic interactions
- Compare ligand orientation with original pose

---

## Expected Learning Outcomes

- Selective inhibition
- Drug–target interactions
- Binding pocket analysis

---

# Experiment 4: Estrogen Receptor Alpha – Tamoxifen

**Protein:** Estrogen Receptor Alpha  
**Ligand:** 4-Hydroxytamoxifen  
**PDB ID:** `3ERT`

---

## Important Notes

- Large hydrophobic cavity
- Easy visualization of ligand pose

---

## Expected Learning Outcomes

- Hormone receptor binding
- Shape complementarity
- Ligand orientation analysis

---

# Result

Students should obtain:

- Docking scores
- Predicted ligand poses
- Interaction profiles
- Hydrogen bonding patterns

---

# Discussion Questions

1. What is the purpose of molecular docking?
2. Why are water molecules removed before docking?
3. What does a negative docking score indicate?
4. Which protein showed the strongest ligand binding?
5. Why is ligand orientation important?

---

# Precautions

- Always remove water molecules before docking
- Ensure correct ligand extraction
- Define grid box properly
- Use identical docking parameters for comparison

---

# Conclusion

This practical introduced the basic workflow of molecular docking, including:

- Protein preparation
- Ligand preparation
- Docking execution
- Interaction analysis

Students learned how computational methods can predict protein–ligand binding and assist in modern drug discovery research.

---

# References

1. Trott O, Olson AJ. AutoDock Vina: Improving the speed and accuracy of docking.
2. Protein Data Bank (PDB): https://www.rcsb.org
3. PyRx Virtual Screening Tool
4. PyMOL Molecular Graphics System

---
