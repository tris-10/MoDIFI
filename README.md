# MoDIFI: Multi-omics Differential Inference for Functional Interpretation

MoDIFI is a Nextflow DSL2 pipeline that integrates ATAC-seq, RNA-seq, Hi-C, and promoter annotations to calculate Differential ATAC–Chromatin–Transcriptome (DACT) scores.
It supports modular workflows, containerized execution (Docker/Singularity), and flexible re-analysis with customizable sample pairings.

### Features

- End-to-end integration of ATAC-seq, RNA-seq, Hi-C, and promoter data
- Modular workflow in Nextflow DSL2
- Portable with Docker or Singularity containers
- Automatic generation of SampleInfo.tsv and SamplePair.tsv
- Re-analysis workflow recalDACT to test new combinations without rerunning the full pipeline


### Requirements

- Nextflow ≥ 23.10
- Java ≥ 11
- Docker or Singularity/Apptainer (recommended for HPC)


### Installation
Clone the repository:
<pre>
git clone https://github.com/your-org/MoDIFI.git
cd DACT
</pre>

### Directory Layout
<pre>
DACT/
├─ modifi_example.nf            # Main Nextflow pipeline
├─ modifi_example.config        # Configuration file
├─ images/
│   └─ modifi.sif       # Singularity container
├─ scripts/           # Pipeline scripts
│   ├─ *.py           # Python scripts
│   └─ *.R            # R scripts
├─ resources/
│   ├─ SampleInfo.tsv   # Auto-generated
│   ├─ SamplePair.tsv   # Auto-generated (editable for recalDACT)
│   └─ Promoter/Promoter.tsv     # External promoter annotations
│   └─ imr90/
│       ├─ atac_seq/
│       │   ├─ *.bed.gz
│       │   ├─ bam/rep1/*.bam
│       │   └─ bam/rep2/*.bam
│       ├─ hic/*.bedpe.gz      # Hi-C data
│       └─ rna_seq/
│           ├─ rep1/*.tsv
│           └─ rep2/*.tsv
└─ output/
</pre>

### Inputs includes
Update the input paths and parameters in <b> dact.config </b>.
#### ATAC-seq:
  - Peak files:
    <pre> ATACpeakFile='/atac_seq/*.bed.gz' </pre>
  - Label for DESeq2 outputs:
    <pre> ATACSeq='ATACseq' </pre>
  - Quality score filter (remove peaks with Q < 5):
    <pre> atac_minQ=5 </pre>        
#### RNA-seq:
  - Gene expression files:
    <pre> RNAFilePattern='/*tsv' </pre>
  - Label for DESeq2 outputs:
    <pre> RNASeq='RNAseq'   </pre>
  - Gene annotation:
    <pre> RNA_ANN_File="${resources_dir}/hg38_annotation.txt" </pre>
  - Quantification type:
    <pre> RNA_quantification='expected_count' </pre>

### Workflow Overview
1. ATAC-seq peaks linked with promoters
2. ATAC-promoter links intersected with Hi-C loops
3. RNA-seq fold changes mapped to gene IDs in loops
4. Normalization restricted to Hi-C regions
5. MoDIFI calculation

### Running the Pipeline

Run the full workflow:

<pre> bash nextflow run modifi_example.nf -c modifi_example.config </pre>

### Outputs include:

- ATAC_counts.txt, ATAC_ann.txt, ATAC_conds.txt
- RNA_counts.txt, RNA_ann.txt, RNA_conds.txt
- Comparison results: GM12878_vs_IMR90_ATACseq.txt, GM12878_vs_IMR90_RNAseq.txt, etc.
