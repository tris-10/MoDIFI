# MoDIFI: Multi-omics Differential Inference for Functional Interpretation

MoDIFI is a Nextflow DSL2 pipeline that integrates ATAC-seq, RNA-seq, Hi-C, and promoter annotations to calculate Multi-omics Differential Inference for Functional Interpretation (MoDIFI) scores.
It supports modular workflows, containerized execution (Docker/Singularity), and flexible re-analysis with customizable sample pairings.

### Features

- End-to-end integration of ATAC-seq, RNA-seq, Hi-C, and promoter data
- Modular workflow in Nextflow DSL2
- Portable with Docker or Singularity containers
- Automatic generation of SampleInfo.tsv and SamplePair.tsv
- Re-analysis workflow recalMoDIFI to test new combinations without rerunning the full pipeline


### Requirements

- Nextflow ≥ 22.10
- Java ≥ 11
- Docker or Singularity/Apptainer (recommended for HPC)


### Installation
Clone the repository:
<pre>
git clone https://github.com/your-org/MoDIFI.git
cd MoDIFI
</pre>

### Directory Layout
<pre>
MoDIFI/
├─ modifi.nf            # Main Nextflow pipeline
├─ modifi_example.config        # Configuration file
├─ modifi.sif         # Singularity container
├─ scripts/           # Pipeline scripts
│   ├─ *.py           # Python scripts
│   └─ *.R            # R scripts
├─ resources/
│   ├─ SampleInfo.tsv   # Auto-generated
│   ├─ SamplePair.tsv   # Auto-generated (editable for recalMoDIFI)
│   ├─ hg38_annotation.txt   # Gene annotation file required for running RNA-seq DESeq2 analysis.
│   ├─ Gnocchi.tsv      # prior information
│   └─ Promoter/Promoter.tsv     # External promoter annotations
│   └─ imr90/
│       ├─ atac_seq/
│       │   ├─ *.bed.gz or *.narrowPeak.gz
│       │   └─ bam/*.bam
│       ├─ hic/*.bedpe.gz      # Hi-C data
│       └─ rna_seq/*.genes.results or *.tsv files
└─ output/
</pre>

### Inputs includes
Update the input paths and parameters in <b> dact.config </b>.
#### ATAC-seq:
  - Peak files:
    <pre> ATACBEDFile= absolute path </pre>
  - BAM files:
    <pre> ATACBAMFiles= absolute path </pre>   
  - Label for DESeq2 outputs:
    <pre> ATACSeq='ATACseq' </pre>
  - Filter out low quality variants:
    <pre> atac_minQ=5 </pre> 
  - The column used for merging:
    <pre> ATAC_Key_col='Region' </pre>      
#### RNA-seq:
  - Gene expression files:
    <pre> RNACountFile= absolute path  </pre>
  - Label for DESeq2 outputs:
    <pre> RRNASeq='RNAseq'   </pre>
  - The column of gene IDs for DEseq anlysis:
    <pre> RNA_Key_col='GeneID'   </pre>
  - The column of RNA counts for DEseq anlysis:
    <pre> RNA_quantification='expected_count'   </pre>  
#### Hi-C:
  - HiC files:
    <pre> HiCLoopsFile= absolute path  </pre>

#### Others:
  - Gene annotation:
    <pre> RNA_ANN_File="${resources_dir}/hg38_annotation.txt" </pre>
  - Prior information:
    <pre> prior_file="${resources_dir}/Gnocchi.tsv" </pre>
  - Prior information:
  <pre> SamplePair = "${resources_dir}/SamplePair.tsv" </pre>

### Workflow Overview
1. ATAC-seq peaks linked with promoters
2. ATAC-promoter links intersected with Hi-C loops
3. RNA-seq fold changes mapped to gene IDs in loops
4. Normalization restricted to Hi-C regions
5. MoDIFI calculation

### Running the Pipeline

Run the full workflow:

<pre> bash nextflow run modifi.nf -c modifi_example.config </pre>

### Outputs include:

- ATAC_counts.txt, ATAC_ann.txt, ATAC_conds.txt
- RNA_counts.txt, RNA_ann.txt, RNA_conds.txt
- Comparison results: [Target]_vs_[Reference]_ATACseq.txt, [Target]_vs_[Reference]_RNAseq.txt, etc.
- MoDIFI results: MoDIFI_all_[Target]_vs_[Reference].tsv, MoDIFI_loop_[Target]_vs_[Reference].tsv 

### Re-running with recalMoDIFI

After the full run, edit resources/SamplePair.tsv to define new Target–Reference comparisons.

Example:
| Target  | Reference         | Check |
|---------|-------------------|-------|
| GM12878 | IMR90             | PASS  |

Then run:

<pre> nextflow run modifi.nf -c modifi_example.config -entry recalMoDIFI </pre>


