# L-ARRAP
Long-read based Antibiotic Resistome Risk Assessment Pipeline (L-ARRAP) calculates the Long-read based Antibiotic Resistome Risk Index (L-ARRI) to 
quantify antibiotic resistome risks.

Building upon our previous ARRI framework, L-ARRAP leverages long-read sequencing advantages to concurrently identify ARGs, 
mobile genetic elements, and human bacterial pathogens, integrating their interactions for risk scoring.

# Instructions

Install dependencies:

conda install -c bioconda chopper minimap2 last centrifuge seqtk seqkit
pip install argparse

Prepare the database:

database/
├── centrifuge/
│   └── p_compressed+h+v
├── minimap2/
│   └── SARG_20211207_14210_filter.ffn
├── mobileOG/
│   └── mobileOG-db_beatrix-1.6.All.faa
├── HBP/
│   └── new-species-all.txt
└── structure/
    └── structure_20181107.LIST

Using L-ARRAp:

# Custom parameter example:
python L-ARRI.py \
  -i my_raw_data \
  -t 32 \
  -a_id 0.8 \
  -a_cov 0.85 \
  -m_id 0.8 \
  -m_cov 0.85 \
  -p map-pb/map-ont


Parameter default value description:

-i, --input rawdata The input directory containing the FASTQ file
-t, -- number of CPU threads used by threads 56
-a_id, --arg_identity 0.75 The minimum identity threshold extracted by ARG
-a_cov, --arg_coverage 0.9 The minimum coverage threshold extracted by ARG
-m_id, --mge_identity 0.75 The minimum identity threshold extracted by MGE
-m_cov, --mge_coverage 0.9 The minimum coverage threshold extracted by MGE

Output files:

├── 3_hbp_abundance/ # HBP-related abundance results
├── ARG_abundance/ # ARG abundance results
├── ARG_raw_result/ # Raw ARG alignment result
├── ARRI/ # ARRI index results
├── Nanofilt_result/ # sequence after quality control
├── all_centrifuge/ # Centrifuge classification results
├── arg-mge/ # ARG-MGE associative results
├── arg_hbp/ # ARG-HBP results (reserved)
├── extract_ARG/ # The extracted ARG result
├── extract_mges/ # Extracted MGE results
├── mges_abundance/ # MGE abundance results
├── raw_mges/ # Raw MGE comparison results
└─ sum_length.txt # Sample sequence length summary
