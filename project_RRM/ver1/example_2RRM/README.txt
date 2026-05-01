all_atom MD Setup — 2026-03-28 12:05
Engine    : OpenMM
Structure : fold_2026_03_23_21_01_model_4_from_cif.pdb
FF        : charmm36.xml, charmm36/water.xml

Run on server:
  conda activate allatom
  python run_openmm.py

NOTE: Small molecule ligands are NOT supported by this script.
      They require separate parameterization (GAFF2 or CGenFF).
