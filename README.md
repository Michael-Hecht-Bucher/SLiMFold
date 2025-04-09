
# SLiMFold Pipeline

The **SLiMFold** pipeline integrates multiple bioinformatics tools to identify, filter, and predict short linear motifs (SLiMs), based on an initial sequence alignment. The pipeline is run in three main steps: **Prerun**, **ColabFold_looped**, and **Postanalysis**. Note that for very large bait proteins (> 2,950 residues), a H100 GPU might be required, which is not directly accessible within ColabFold. 

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

> This notebook is where the pipeline begins. It creates a position-specific scoring matrix (PSSM) from user-provided input SLiMs and searches the proteome to find candidate sequences. Key filters (PSSM score, IUPRED, ANCHOR and PSIPRED) are applied to reduce false positives, and each hit is paired with the user’s chosen “bait” (e.g., human actin). Multiple sequence alignments (MSAs) are then generated in parallel and reformatted for downstream structure prediction.

<details>
  <summary>Details</summary>
  
0. **Open Prerun.ipynb**
  
1. **Folder and pathway setup**
   - Select the kernel ```SLiM_AF2_screen```
   - Define the paths ```iupred_path```, ```psipred_path```, ```NCBI_protein_database```, ```uniref90_path```, ```reformat_path``` and your ```bait_sequence```. 
   - Execute the cell, enter a project name in the prompt. A consistent project folder structure will be automatically created.  
   - Move your initial FASTA-file to the **Input Folder** and rename it to **input.fasta**. Please make sure that input sequences contain only **the motif without flanking residues** (see example folder). Input motifs should have the same sequence length! 

2. **PSSM Generation with BLOSUM62**  
   - Uses input.fasta and the BLOSUM62 substitution matrix to generate an initial position-specific scoring matrix (PSSM) as CSV-file-output (stored in ```{project_name}/Output/pssm_BLOSUM62.csv```)

3. **Proteome Search**
   - (A) Defines several thresholds for subsequent motif identification (```pssm_cutoff```, ```iupred_cutoff```, ```anchor_cutoff```, secondary structure cutoffs for helix, strand, coil or unknown). Prompts to define the probable secondary structure (of the motif) involved in the interaction. Choose bewtween 'helix', 'strand', 'coil' or 'unknown'.
     
   - (B) Scores the human proteome (or your proteome of choice) using the PSSM, as well as IUPRED, ANCHOR, PSIPRED. Retains only hits meeting specified cutoffs. Extends each hit by ±20 residues to capture potential context (can be modified by changing ```flanking_aa_size```). This will produce an output FASTA-file containing identified hits (stored in ```{project_name}/Output/PSSM_Hits/Hits.fasta```)
     
   - (C) Removes identical sequences to avoid running them through jackhmmer and colabfold multiple times. This will produce another FASTA-file containing only non-redundant hits (stored in ```{project_name}/Output/PSSM_Hits/Hits_nonred.fasta```)
     
   - (Optional, if not first iteration): Compare the PSSM-hits of two iterations and write the unique hits to a new FASTA file. Please ignore this cell in case you are running the first iteration.

4. **Bait Fusion and Prey-Bait Preparation**
   - Input: ```Hits_nonred.fasta``` generated in the previous step.
   - The predefined bait sequence is appended to each unique hit, separated by a colon (> header as peptide:bait).
   - Outputs a formatted FASTA-file ```{project_name}/Output/PSSM_Hits/PreyBait.fasta```

5. **Split PreyBait.fasta into individual FASTA files for ColabFold input** 
   - Input: ```PreyBait.fasta``` generated in the previous step.
   - Creates for each PreyBait Sequence pair an individual FASTA file (stored in ```{project_name}/Output/Fasta/```)

6. **Multiple Sequence Alignment for Bait** 
   - The predefined bait sequence is run with jackhmmer (with modified filters) against the UniRef90 database to identify homologs and generates a .sto alignment file (stored in ```{project_name}/Output/MSA/sto```).
   - The filters can be modified by changing ```-E```, ```-N```, ```-F1```, ```-F2``` or ```-F3``` 

7. **Multiple Sequence Alignment for Peptides** 
   - Each peptide is run with jackhmmer (with modified filters) against the UniRef90 database to identify homologs and generate a .sto alignment file (stored in ```{project_name}/Output/MSA/sto```).
   - The filters can be modified by changing ```-E```, ```-N```, ```-F1```, ```-F2``` or ```-F3``` 
   - To speed up computation parallel processing is used. Both, the number of CPU cores per search (```num_cpus_per_process```)  and the number of parallel processes (```num_processes```) can be adjusted.
   - Automatically tracks remaining peptides, so the run can resume from where it left off using the ```input_remaining.fasta``` file in case of interruption.

8. **Converts the .sto to .a3m** 

9. **Sort and Deduplicate .a3m Files Based on Sequence Identity**
   - Sorts all .a3m files by global sequence identity to the reference (first) sequence, placing the most similar sequences at the top to improve MSA quality for structure prediction.
   - For the bait MSA (bait_sequence.a3m), the user is prompted whether they want to sort it.
   - Deduplicates the bait .a3m file (based on exact sequence match) to remove redundant homologs, ensuring higher sequence diversity and enhancing co-evolutionary signal strength for better complex prediction accuracy.

10. **Trims the MSA**
    - Reduces the size of each .a3m file by keeping only the first N sequences (default: 2048).

11. **Combines Bait and Prey MSAs for ColabFold**

</details>

---

## 2. ColabFold_looped.ipynb

> This notebook is a modified version of the [ColabFold batch pipeline](https://github.com/sokrypton/ColabFold) originally developed by the Steinegger lab. In our pipeline, ColabFold_looped.ipynb automates structure predictions of candidate motif–bait pairs (e.g., SLiM–actin). The key difference from the original ColabFold batch notebook is the ability to loop through multiple FASTA files and their associated custom A3M files, while also allowing users to specify the number of seeds for increased model diversity.

<details>
  <summary>Details</summary>

1. **Preparation**  
   - Upload the FASTA files (Output/Fasta/) and the custom MSAs (Output/combined_a3m/) you generated in Prerun.ipynb to your Google Drive.
   - Open ColabFold_looped.ipynb in Google Colab, connect to a runtime, and select a GPU (we recommend using an A100 for faster inference).
   - Set the paths to your uploaded:
     -   fasta_directory (FASTA files)
     -   msa_directory (custom A3Ms)
     -   result_directory (where predictions will be saved)
   - Under **msa_mode**, choose whether to use
     -   custom – uses your **precomputed MSAs** (recommended for peptides)
     -   MMseqs2_uniref_env – generates MSAs on the fly (often insufficient for short sequences like peptides)

2. **Running the Prediction**  
   - Run the main prediction cell: the script will automatically loop through all FASTA files, pair them with the corresponding .a3m files, and run structure prediction for each pair.  
   - Prediction results are saved in your defined result_directory.
   - If the Colab runtime disconnects (e.g., after 24 hours), don't worry:
     -   Already processed FASTA files are moved into a done folder.
     -   Simply reconnect to the notebook and rerun the prediction cell to continue from where it left off.


</details>

---

## 3. Postanalysis.ipynb

> Once you have raw predictions from **ColabFold_looped**, the **Postanalysis** stage aggregates, filters, and clusters candidate structures to identify meaningful F-actin-binding SLiMs. This notebook helps evaluate each predicted motif’s structural reliability and organizes results for detailed inspection and downstream analyses (e.g., functional enrichment).

<details>
  <summary>Details</summary>

1. **Folder and pathway setup**
   - Please download the results (zip files) and place them in a folder 
   - Please define the paths where
     - the zip files are stored (zip_files_folder)
     - the fastas are stored (Output/Fastas)
     - the reference_pdb_path (this is important, for RMSD and angle calculation. The pdb should be in the same format (meaning number of residues and chains, as well as reihenfolge) as your predicted structures
     - the outputs are stored (output_directory)
   - Automatically creates a consistent project folder structure.

2. **Unpacking**
   - Unpacks all the zip files

3. **Analysis of Model Metrics and Structural Comparisons**
   - Unpacks all the zip files
   - Reads all predicted structures and extracts:  
     - pLDDT: Per-residue confidence.  
     - pTM & ipTM: Global and interface metrics indicating interchain confidence. For this value the mean of the Top 3 models are calculated. 
   - Also calculates RMSD (Root Mean Square Deviation)(Calculates RMSD over the alpha-carbon atoms in the motif region (P1–P9)), spherical angles (φ (azimuth) and θ (polar))(Generates Δφ and Δθ values by comparing each predicted motif’s orientation to the reference.), and helix polarity by comparing the predicted models against a reference structure.

4. **Filter Combined Results by ipTM Cutoff**
   - Excludes structures with ipTM < 0.6 (default), which generally indicates poor interface reliability.
   - Creates a scatter plot, showing Mean ipTM vs. RMSD. 

5. **Visualization of 2D and 3D Scatter Plots for Protein Metrics**
   - Visualizes the relationships between three angular dimensions (Delta Angles Theta and Phi and Helix Polarity) and Mean RMSD values of the predicted Hits.
   - The data is displayed using 2D and 3D scatter plots, with a blue-to-red colormap for the RMSD.

6. **Optimizing Clustering Parameters with differenet Algorithms and Evaluating Cluster Quality**
   - This script searches for the optimum clustering parameters of KMeans, Agglomerative and HDBScan  and evaluates cluster quality using various metrics such as silhouette score, Calinski-Harabasz score, and Davies-Bouldin score. The results are visualized to help determine the best clustering configuration. Insights from these scores can guide the selection of n cluster size, min_cluster_size and min_samples.

7. **Clustering**
   - You can choose between three clustering methods, and choose the size based on the above calcualted silhouette score, Davies-Bouldin index, and Calinski-Harabasz index:
     - (A) Kmeans: Needs cluster size as input. Maybe you can add like one sentence to this methdod.
     - (B) Agglomerative: Needs cluster size as input. Maybe you can add like one sentence to this methdod.
     - (C) HDBScan: Needs minimum cluster size and minimum samples as input. Maybe you can add like one sentence to this methdod.
   - All perform clustering using RMSD, Δφ, Δθ, and polarity as features.
   - 2D and 3D scatter plots are created for cluster visualization and Cluster-wise metrics (PSSM Score, ipTM Score, RMSD, IUPRED and ANCHOR Score)

HIER STEHEN GEBLIEBEN

3. **Structural Alignment & RMSD**  
   - Aligns each candidate structure to a reference PDB (e.g., ITPKA–actin complex).  
   - Calculates RMSD over the alpha-carbon atoms in the motif region (P1–P9).  
     - RMSD quantifies how closely the predicted motif aligns to the known reference.

4. **Angular Measurements**  
   - Computes φ (azimuth) and θ (polar) angles to represent the motif’s orientation.  
   - Calculates helix polarity to capture directionality.  
   - Generates Δφ and Δθ values by comparing each predicted motif’s orientation to the reference.

5. **Clustering**  
   - Performs HDBSCAN clustering using RMSD, Δφ, Δθ, and polarity as features.  
   - Chooses optimal clustering parameters (e.g., minimum cluster size, minimum samples) based on silhouette score, Davies-Bouldin index, and Calinski-Harabasz index.

6. **Cluster Examination & Data Export**  
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

