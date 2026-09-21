#!/usr/bin/env python3

import json
import os
from copy import deepcopy

import pandas as pd
from rdkit import Chem


# ============================================================
# User Settings
# ============================================================

CSV_FILE = "ligands.csv"
TEMPLATE_JSON = "PB1_template.json"
OUTPUT_DIR = "AF3_inputs"

# ============================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# Function: Remove explicit hydrogens and canonicalize SMILES
# ============================================================

def clean_smiles(smiles):
    """
    Convert hydrogen-explicit SMILES into standard SMILES.

    Steps:
      1. Parse original SMILES with RDKit
      2. Remove explicit H atoms
      3. Sanitize molecule
      4. Generate canonical isomeric SMILES

    Returns:
      cleaned_smiles, error_message
    """

    if pd.isna(smiles):
        return None, "SMILES is missing"

    smiles = str(smiles).strip()

    if not smiles:
        return None, "SMILES is empty"

    try:
        # Parse the original SMILES
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            return None, "RDKit could not parse SMILES"

        # Remove explicit hydrogens
        mol = Chem.RemoveHs(mol)

        # Re-sanitize after hydrogen removal
        Chem.SanitizeMol(mol)

        # Generate canonical SMILES while retaining stereochemistry
        cleaned = Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=True
        )

        if not cleaned:
            return None, "RDKit generated empty SMILES"

        return cleaned, None

    except Exception as e:
        return None, str(e)


# ============================================================
# Read template JSON
# ============================================================

with open(TEMPLATE_JSON, "r") as f:
    template = json.load(f)


# ============================================================
# Read ligand table
# ============================================================

ligands = pd.read_csv(CSV_FILE)

required_columns = ["Ligand", "SMILES"]

for col in required_columns:
    if col not in ligands.columns:
        raise ValueError(
            f"Missing required column: {col}"
        )

print(f"Found {len(ligands)} ligands.")
print()


# ============================================================
# Prepare failure log
# ============================================================

failed = []


# ============================================================
# Generate one JSON per ligand
# ============================================================

generated = 0

for index, row in ligands.iterrows():

    ligand_name = str(row["Ligand"]).strip()
    original_smiles = row["SMILES"]

    print(
        f"[{index + 1}/{len(ligands)}] "
        f"Processing {ligand_name}"
    )

    # --------------------------------------------------------
    # Clean / hydrogen-depleted SMILES
    # --------------------------------------------------------

    cleaned_smiles, error = clean_smiles(original_smiles)

    if error:

        print(f"    ERROR: {error}")

        failed.append({
            "Ligand": ligand_name,
            "Original_SMILES": original_smiles,
            "Error": error
        })

        continue

    print(f"    Original : {original_smiles}")
    print(f"    Cleaned  : {cleaned_smiles}")

    # --------------------------------------------------------
    # Copy template
    # --------------------------------------------------------

    data = deepcopy(template)

    # Rename prediction
    data["name"] = ligand_name

    # --------------------------------------------------------
    # Find ligand block and replace SMILES
    # --------------------------------------------------------

    ligand_found = False

    for item in data["sequences"]:

        if "ligand" in item:

            item["ligand"]["smiles"] = cleaned_smiles
            ligand_found = True
            break

    if not ligand_found:

        raise RuntimeError(
            "No ligand entry found in template JSON."
        )

    # --------------------------------------------------------
    # Output filename
    # --------------------------------------------------------

    outfile = os.path.join(
        OUTPUT_DIR,
        f"{ligand_name}.json"
    )

    # --------------------------------------------------------
    # Write JSON
    # --------------------------------------------------------

    with open(outfile, "w") as f:
        json.dump(
            data,
            f,
            indent=2
        )

    generated += 1

    print(f"    JSON     : {outfile}")
    print()


# ============================================================
# Write failure report
# ============================================================

if failed:

    failed_df = pd.DataFrame(failed)

    failed_file = os.path.join(
        OUTPUT_DIR,
        "failed_smiles.csv"
    )

    failed_df.to_csv(
        failed_file,
        index=False
    )

else:
    failed_file = None


# ============================================================
# Summary
# ============================================================

print("=" * 60)
print("DONE")
print("=" * 60)

print(f"Input ligands       : {len(ligands)}")
print(f"AF3 JSON files      : {generated}")
print(f"Failed SMILES       : {len(failed)}")
print(f"Output directory    : {OUTPUT_DIR}")

if failed_file:
    print(f"Failure report      : {failed_file}")
