
# SLiMFold Pipeline

The **SLiMFold** pipeline integrates multiple bioinformatics tools to identify, filter, and predict short linear motifs (SLiMs), based on an initial sequence alignment. The pipeline was run in three main steps: **Prerun**, **ColabFold_looped**, and **Postanalysis**.

![Alt text](images/Pipeline.png)

---

## Installation & Environment

<details>
  
1. **Clone This Repo & Create the Conda Environment**  
   ```bash
   git clone https://github.com/YourUserName/SLiM_AF2_screen.git
   cd SLiM_AF2_screen

   conda env create -f slim_env.yml
   conda activate SLiM_AF2_screen
   ```
2. **Register as Jupyter Kernel** (optional, but recommended)
   ```bash
   python -m ipykernel install --user --name SLiM_AF2_screen --display-name "SLiM_AF2_screen"
   ```
3. **Install External Tools**  
   - **PsiPred 4.0**: [psipred GitHub](https://github.com/psipred/psipred)  
   - **IUPred3**: [iupred3.elte.hu](https://iupred3.elte.hu/download_new)  
   - **Databases**: [UniRef90 in .fasta.gz](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/) & [NCBI protein dataset in .fasta](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/)

</details>


---


## 1. Prerun.ipynb

> This notebook is where the pipeline begins. It creates a position-specific scoring matrix (PSSM) from user-provided input SLiMs and searches the proteome to find candidate sequences. Key filters (IUPRED, ANCHOR, PSIPRED, etc.) are applied to reduce false positives, and each hit is paired with the user’s chosen “bait” (e.g., human actin). Multiple sequence alignments (MSAs) are then generated in parallel and reformatted for downstream structure prediction.

<details>
  <summary>Details in <code>1.Prerun.ipynb</code></summary>

1. **Folder and pathway setup**  
   - Please define the paths iupred_path, psipred_path, NCBI_protein_database, uniref90_path, reformat_path and your bait_sequence. 
   - Automatically creates a consistent project folder structure.  
   - Requires the user’s environment to be active (e.g. `conda activate SLiM_AF2_screen`).
   - Move your inital aligment to the **Input Folder** and rename it to **input.fasta** (see example folder). 

2. **PSSM Generation with BLOSUM62**  
   - Uses user-provided SLiMs (aligned FASTA in input.fasta) and the BLOSUM62 substitution matrix.  
   - Produces a PSSM cutoff (default for BLOSUM62 is set to 10).

3. **Proteome Search**  
   - (A) Prompts the user first to define the probable secondary structure involved in the interaction. User can choose bewtween 'helix', 'strand', 'coil' or 'unknown'. 
   - (B) Then scores the human proteome (or your organism of choice) using the PSSM, as well as IUPRED, ANCHOR, PSIPRED. Retains only hits meeting specified cutoffs for PSSM, IUPRED, ANCHOR, PSIPRED, etc. Extends each hit by ±20 residues to capture potential context.
   - (C) Removes identical sequences found to avoid running them through jackhmmer and colabfold multiple times.
   - (Optional, if not first iteration) Compare the PSSM-hits of two iterations and write the unique hits to a new FASTA file. Please ignore this cell in case you are running the first iteration.

4. **Bait Fusion and Prey-Bait Preparation**
   - The predefined bait sequence is appended to each unique hit, separated by a colon (> header as peptide:bait).

5. **Split Prey-Bait pairs into individual FASTA files for ColabFold input** 

6. **Multiple Sequence Alignment for Bait** 
   - Runs jackhmmer for the bait sequence, with modified filters, against the UniRef90 database to identify homologs and generate a .sto alignment file.

7. **Multiple Sequence Alignment for Peptides** 
   - Runs jackhmmer for each peptide, with modified filters, against the UniRef90 database to identify homologs and generate a .sto alignment file.
   - Uses parallel processing to speed up computation — both the number of CPU cores per search and the number of parallel processes can be adjusted by the user.
   - Automatically tracks remaining peptides, so the run can resume from where it left off using the *input_remaining.fasta* file in case of interruption.

8. **Converts the .sto to .a3m** 

9. **Sort and Deduplicate .a3m Files Based on Sequence Identity**
   - Sorts all .a3m files by global sequence identity to the reference (first) sequence, placing the most similar sequences at the top to improve MSA quality for structure prediction.
   - For the bait MSA (bait_sequence.a3m), the user is prompted whether they want to sort it.
   - Deduplicates the bait .a3m file (based on exact sequence match) to remove redundant homologs, ensuring higher sequence diversity and enhancing co-evolutionary signal strength for better complex prediction accuracy.

10. **Trims the MSA**
    - Reduces the size of each .a3m file by keeping only the first N sequences (default: 2048).

11. **Combines Bait and Prey MSAs for ColabFoldA**

</details>

---

## 2. ColabFold_looped.ipynb

> **Description**  
> This notebook is a **modified version** of the [ColabFold batch pipeline](https://github.com/sokrypton/ColabFold) originally developed by the Steinegger lab. In our pipeline, **ColabFold_looped.ipynb** automates structure predictions of candidate motif–bait pairs (e.g., SLiM–actin) by:
>
> 1. Accepting custom MSAs from the `Prerun` stage.  
> 2. Looping through all FASTA files in a specified folder.  
> 3. Integrating multiple sequences (SLiM + bait) for AlphaFold2 Multimer predictions.

<details>
  <summary>Key Steps in <code>ColabFold_looped.ipynb</code></summary>

1. **Batch Processing & Automation**  
   - Scans a given folder (e.g., in Google Drive) for all `.fasta` files.  
   - Automatically runs AlphaFold2 Multimer for each file–MSA pair.  
   - No manual input required for each sequence, speeding up large-scale screening.

2. **Integration of Custom MSAs**  
   - If `.a3m` alignment files exist from **Prerun** (stored in a designated “MSA” folder), the notebook **matches each `.a3m`** to its corresponding `.fasta` by name.  
   - These custom MSAs help improve structure predictions by providing more specific alignments.

3. **Configurable Output**  
   - Allows custom output folders to keep results organized.  
   - Lets you adjust seeds and other AlphaFold2 parameters, enabling fine-tuning of your predictions.

4. **AlphaFold2 Multimer Predictions**  
   - Uses the combined sequences (motif + bait) to generate 3D models.  
   - Logs pLDDT, pTM, and ipTM scores to measure model confidence and interface quality.

</details>

---

## 3. Postanalysis.ipynb

> **Description**  
> Once you have raw predictions from **ColabFold_looped**, the **Postanalysis** stage aggregates, filters, and clusters candidate structures to identify meaningful F-actin-binding SLiMs. This notebook helps evaluate each predicted motif’s structural reliability and organizes results for detailed inspection and downstream analyses (e.g., functional enrichment).

<details>
  <summary>Key Steps in <code>Postanalysis.ipynb</code></summary>

1. **Loading & Filtering Models**  
   - Reads all predicted structures and extracts:  
     - pLDDT: Per-residue confidence.  
     - pTM & ipTM: Global and interface metrics indicating interchain confidence.  
   - Excludes structures with ipTM < 0.6 (default), which generally indicates poor interface reliability.

2. **Structural Alignment & RMSD**  
   - Aligns each candidate structure to a reference PDB (e.g., ITPKA–actin complex).  
   - Calculates RMSD over the alpha-carbon atoms in the motif region (P1–P9).  
     - RMSD quantifies how closely the predicted motif aligns to the known reference.

3. **Angular Measurements**  
   - Computes φ (azimuth) and θ (polar) angles to represent the motif’s orientation.  
   - Calculates helix polarity to capture directionality.  
   - Generates Δφ and Δθ values by comparing each predicted motif’s orientation to the reference.

4. **Clustering**  
   - Performs HDBSCAN clustering using RMSD, Δφ, Δθ, and polarity as features.  
   - Chooses optimal clustering parameters (e.g., minimum cluster size, minimum samples) based on silhouette score, Davies-Bouldin index, and Calinski-Harabasz index.

5. **Cluster Examination & Data Export**  
   - Structures in each cluster are exported to `.pml` files for inspection in PyMOL.  
   - Corresponding sequences are compiled into FASTA files, enabling:  
     - Sequence logo generation to identify conserved positions.  
     - Gene ontology (GO) enrichment analysis (by mapping each sequence to its gene via NCBI).

</details>

---


## Citation

If you use this pipeline in published research, please cite:
- Your own manuscript
- Tools like AlphaFold2, ColabFold, IUPRED, ANCHOR, PSIPRED, HDBSCAN

---

