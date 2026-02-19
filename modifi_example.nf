nextflow.enable.dsl=2

process extractPeaks {
    tag "$cell"

    container 'modifi.sif'

    publishDir "${params.output_dir}", mode: 'copy', saveAs: { "${cell}.narrowPeak" }

    input:
    val cell

    output:
    path "${cell}.narrowPeak", emit: peaks, optional:true

    script:
    """
    zcat ${params.resources_dir}/${cell.toLowerCase()}${params.ATACpeakFile} | sed "s/Peak/${cell}_Peak/g" > ${cell}.narrowPeak
    """
}


process ConcatenateAndSortPeaks {
    container 'modifi.sif'
    
    publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename -> filename }

    input:
    path peaks

    output:
    path "merge.bed", emit: merge_bed

    script:
    """
    cat ${peaks.join(' ')} > all.narrowPeak
    bedtools sort -i all.narrowPeak > sort.narrowPeak
    bedtools merge -i sort.narrowPeak -c 4,9 -o collapse,min | \\
    awk '{
        split(\$4, a, ",");
        delete seen;
        out = "";
        for (i in a) {
            if (!(a[i] in seen)) {
                seen[a[i]] = 1;
                out = out ? out "," a[i] : a[i];
            }
        }
        print \$1, \$2, \$3, out, \$5;
    }' OFS="\\t" > tmp_merge.bed
    # Filter by column 5 and then remove the quality column
    awk '\$5 > ${params.atac_minQ}' tmp_merge.bed | cut -f1-4 > merge.bed    
    """
}


process GenerateDESeq2InputsCounts {
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename -> "${cell}_R${RepNum}_counts.txt" }

    input:
    tuple path(merge_bed), val(cell), val(RepNum)

    output:
    path "${cell}_R${RepNum}_counts.txt", emit: ATAC_counts, optional:true
    
    script:
    """
    set -euxo pipefail
    BAM_DIR="${params.resources_dir}/${cell.toLowerCase()}/atac_seq/bam/rep${RepNum}"
    echo "Checking BAM files in \${BAM_DIR}" | tee -a ${cell}_R${RepNum}_process.log
    OUTPUT_FILE="${cell}_R${RepNum}_counts.txt"

    if [ -d "\${BAM_DIR}" ]; then
        BAM_FILES=(\${BAM_DIR}/*.bam)
        if [ \${#BAM_FILES[@]} -eq 0 ]; then
            echo "No BAM files found in \${BAM_DIR}" | tee -a ${cell}_R${RepNum}_process.log
        else
            for bamFile in "\${BAM_FILES[@]}"; do
                if [ -f "\${bamFile}" ]; then
                    if [ ! -f "\${bamFile}.bai" ]; then
                        echo "Index file missing for \${bamFile}, generating now..." | tee -a ${cell}_R${RepNum}_process.log
                        samtools index \${bamFile}
                    fi
		    if [ "$RepNum" -le ${params.RepNUM_atac["${cell}"]} ]; then
                    	python ${params.script_dir}/atac_seq_counts.py ${merge_bed} \${bamFile} \${OUTPUT_FILE} ${cell}_R${RepNum} --minQ ${params.atac_minQ}
                    fi
		else
                    echo "BAM file \${bamFile} not found" | tee -a ${cell}_R${RepNum}_process.log
                fi
            done
        fi
    else
        echo "Directory \${BAM_DIR} does not exist" | tee -a ${cell}_R${RepNum}_process.log
    fi
    """
}



process GenerateDESeq2InputsAnnCOUNTs{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename -> filename }

    input:
    path merge_bed
    path ATAC_counts 

    output:
    path "ATAC_counts.txt", emit: ATAC_Counts
    path "ATAC_ann.txt", emit: ATAC_ann

    script:
    """
    echo 'Regions' > region_names.txt
    cut -f 4 ${merge_bed} >> region_names.txt
    paste region_names.txt *R*_counts.txt > ATAC_counts.txt

    echo -e 'Regions\\tCoords' > ATAC_ann.txt
    awk '{ printf "%s\\t%s:%s-%s\\n", \$4, \$1, \$2, \$3 }' ${merge_bed} >> ATAC_ann.txt
    
    """    
}

process GenerateDESeq2InputsConds{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename -> filename }

    input:
    path ATAC_Counts

    output:
    path "ATAC_conds.txt", emit: ATAC_conds

    script:
    """
    python ${params.script_dir}/getATAC_conds.py -i "${ATAC_Counts}" 
    """

}

process runDEseq{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', pattern:"*_vs_*_${data_type}.txt", saveAs: { filename -> filename }

    input:
    path ATAC_Counts_file
    path ATAC_conds_file
    path ATAC_ann_file
    val ANN_COL
    val data_type

    output:
    path "*_vs_*_${data_type}.txt", emit: DEseq_atac, optional:true

    script:
    """
    Rscript ${params.script_dir}/atac_deseq.R ${ATAC_Counts_file} ${ATAC_conds_file} ${ATAC_ann_file} "\$PWD" ${ANN_COL} ${data_type}
    """    
}

process runDEseqRSEM{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', pattern:"*_vs_*_${data_type}.txt", saveAs: { filename -> filename }

    input:
    path RNA_conds_file
    path RNA_ann_file
    val ANN_COL
    val data_type

    output:
    path "*_vs_*_${data_type}.txt", emit: DEseq_rna, optional:true
          
    script:
    """   
    Rscript ${params.script_dir}/rna_deseq.R ${RNA_conds_file} ${RNA_ann_file} "\$PWD" ${ANN_COL} ${data_type}
    """    
}

process getDEseqInputForRNAseq{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', pattern:"RNA_*.txt", saveAs: { filename -> filename }

    output:
    path "RNA_counts.txt", emit: RNA_counts, optional:true
    path "RNA_conds.txt", emit: RNA_conds, optional:true
    path "RNA_ann.txt", emit: RNA_ann, optional:true

    script:
    """
    python ${params.script_dir}/GenerateDESeq2InputsFromRNAseq.py \\
	-i "${params.resources_dir}/" \\
	-l "${params.RepNUM_rna.keySet()}" \\
        -r "${params.RepNUM_rna.values()}" \\
        -p "${params.RNAFilePattern}" \\
        --RNA_COL "${params.RNA_quantification}" 
    """
}

process getDEseqInputForRNAseqRSEM{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', pattern:"RNA_*.txt", saveAs: { filename -> filename }

    output:
    path "RNA_conds.txt", emit: RNA_conds, optional:true

    script:
    """
    python ${params.script_dir}/GenerateDESeq2InputsForRNAseqRSEM.py \\
        -i "${params.resources_dir}/" \\
        -l "${params.RepNUM_rna.keySet()}" \\
        -r "${params.RepNUM_rna.values()}" \\
        -p "${params.RNAFilePattern}" 
    """
}


process generateSampleInfo{
    container 'modifi.sif'
    publishDir "${params.resources_dir}", mode: 'copy', pattern:"*Sample*.tsv", saveAs: { filename -> filename }
    publishDir "${params.output_dir}", mode: 'copy', pattern:"*piece*.tsv", saveAs: { filename -> filename }


    input:
    path DEseq_atac
    path DEseq_rna

    output:
    path "SamplePair.tsv", emit: SamplePair, optional:true
    path "SampleInfo.tsv", emit: SampleInfo, optional: true
    path "piece*.tsv", emit: Piece, optional: true    

    script:
    """
    python ${params.script_dir}/CheckSamplesForMoDIFI.py \\
	-i "${params.resources_dir}/" \\
	-o "${params.output_dir}/" \\
	--LABEL_RNA "${params.RNASeq}" \\
	--LABEL_ATAC "${params.ATACSeq}" \\
        -sp ${params.SamplePair}
    """
}

process generateSampleInfo_noInput{
    container 'modifi.sif'
    publishDir "${params.resources_dir}", mode: 'copy', pattern:"*Sample*.tsv", saveAs: { filename -> filename }
    publishDir "${params.output_dir}", mode: 'copy', pattern:"*piece*.tsv", saveAs: { filename -> filename }


   // input:
   // path DEseq_atac
   // path DEseq_rna

    output:
    path "SamplePair.tsv", emit: SamplePair, optional:true
    path "SampleInfo.tsv", emit: SampleInfo, optional: true
    path "piece*.tsv", emit: Piece, optional: true

    script:
    """
    python ${params.script_dir}/CheckSamplesForMoDIFI.py \\
        -i "${params.resources_dir}/" \\
        -o "${params.output_dir}/" \\
        --LABEL_RNA "${params.RNASeq}" \\
        --LABEL_ATAC "${params.ATACSeq}" \\
	-sp ${params.SamplePair} 
    """
}

process MapATACWithPRO{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', pattern:"PRO_ATAC*.tsv", saveAs: { filename -> filename }

    input:
    path SampleInfo
    path SamplePair

    output:
    path "PRO_ATAC*tsv", emit: PRO_ATAC, optional:true

    script:
    """
    python ${params.script_dir}/MapATACWithPRO.py \\
        -i "${params.resources_dir}/" \\
        -o "${params.output_dir}/" \\
	-r 5000 \\
	-m 0.5 \\
        --LABEL_RNA "${params.RNASeq}" \\
        --LABEL_ATAC "${params.ATACSeq}" \\
	-si "${params.SampleInfo}" \\
	-sp "${params.SamplePair}"
    """
}

process MapATACWithHiCWithPRO{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', pattern:"ToBacon_*.tsv", saveAs: { filename -> filename }
    publishDir "${params.output_dir}", mode: 'copy', pattern:"atacWithHiC*.tsv", saveAs: { filename -> filename }
    publishDir "${params.output_dir}", mode: 'copy', pattern:"PRO_ATAC*.tsv", saveAs: { filename -> filename }

    input:
    path piece_file

    output:
    path "ToBacon_*tsv", emit: ToBacon, optional:true
    path "atacWithHiC*tsv", emit: HiC, optional:true
    path "PRO_ATAC*tsv", emit: PRO_ATAC, optional:true
    
    script:
    """
    echo "Using piece file: ${piece_file}"
    python ${params.script_dir}/MapATACWithHiCWithPRO.py \\
        -i "${params.resources_dir}/" \\
        -o "${params.output_dir}/" \\
        --LABEL_RNA "${params.RNASeq}" \\
        --LABEL_ATAC "${params.ATACSeq}" \\
        -si "${params.SampleInfo}" \\
        -sp "${params.output_dir}/${piece_file}" \\
        --PRO_Region 5000 \\
        --PRO_minOL 0.5
    """
}

process AdjustZInflation{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', pattern:"ToBacon*bacon*.tsv", saveAs: { filename -> filename }

    input:
    path ToBacon

    output:
    path "ToBacon*bacon*tsv", emit: adjZ, optional:true

    script:
    """
    Rscript ${params.script_dir}/baconForDEseq.R ${ToBacon}
    """
}

process calMoDIFY{
    container 'modifi.sif'
    publishDir "${params.output_dir}", mode: 'copy', pattern:"MoDIFI*.tsv", saveAs: { filename -> filename }

    input:
    path adjZ

    output:
    path "MoDIFI*tsv", emit: MoDIFI, optional:true

    script:
    """
    echo "Z scores were adjusted: ${adjZ}"
    python ${params.script_dir}/calMoDIFY.py \\
        -i "${params.output_dir}/${adjZ}" \\
        -o "${params.output_dir}/" \\
        -p "${params.prior_file}" 
    """
}


process generateIntermediateFile {
    //publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename -> filename }

    input:
    val atacBamfile 

    output:
    path "intermediate.txt"

    script:
    """
    touch intermediate.txt
    for cell in \$(echo \${atacBamfile.keySet().join(' ')})
    do
        idx=1
        for path in \$(echo \${atacBamfile[cell].join(' ')})
        do
            echo "Processing cell: \${cell}, path: \${path}" >> intermediate.txt
            find "${path}" -name '*.bam' | while read bamfile
            do
                echo \${cell} \${idx} \"\${bamfile}\" >> intermediate.txt
            done
            ((idx++))
        done
    done
    """
}


workflow ATACseq_workflow{
    // Processing ATAC-seq signals and get three inputfiles for DEseq
    cellChannel = Channel.from(params.RepNUM_atac.keySet())
    Peaks = extractPeaks(cellChannel).peaks.collect()
    Merge_bed = ConcatenateAndSortPeaks(Peaks).collect()
    crossProductChannel = Merge_bed.combine(cellChannel).combine(Channel.from(1..params.RepNUM_atac.values().max()))
    crossProductChannel.view { item ->
        println "Bed File: ${item[0]}, Cell: ${item[1]}, RepNum: ${item[2]}"
    }
    ATAC_COUNT = GenerateDESeq2InputsCounts(crossProductChannel).ATAC_counts.collect()
    DEseqInputFile = GenerateDESeq2InputsAnnCOUNTs(Merge_bed, ATAC_COUNT)
    ATAC_CONDS = GenerateDESeq2InputsConds(DEseqInputFile.ATAC_Counts)
    DEseq_atac = runDEseq(DEseqInputFile.ATAC_Counts, ATAC_CONDS, DEseqInputFile.ATAC_ann, params.ATAC_Key_col,params.ATACSeq)

    emit:
    DEseq_atac
}

workflow RNAseq_workflow {
    // DEseq for RNA-seq
    DEseqInputRNA = getDEseqInputForRNAseqRSEM()
    DEseq_rna = runDEseqRSEM(DEseqInputRNA.RNA_conds, params.RNA_ANN_File, params.RNA_Key_col, params.RNASeq) 
    DEseq_rna.view()
    
    emit:
    DEseq_rna
}

workflow recalMoDIFI{
    InfoPair = generateSampleInfo_noInput()
    InfoPair.Piece.flatten()
        .view { item ->
            println "File: ${item}"
        }
        .set { pieceInput }
    Combine = MapATACWithHiCWithPRO(pieceInput)
    adjZ = AdjustZInflation(Combine.ToBacon)
    calMoDIFY(adjZ)
}


workflow {
    DEseq_atac = ATACseq_workflow()
    DEseq_rna = RNAseq_workflow()
    DEseq_rna.view()
    
    //DEseq_rna.view { item ->
    //        println "File: ${item}"
    //    }
    //    .set { DEseq_rnaInput }
    InfoPair = generateSampleInfo(DEseq_atac, DEseq_rna)
    //InfoPair = generateSampleInfo_noInput()
    InfoPair.Piece.flatten()
        .view { item -> 
            println "File: ${item}" 
        }
        .set { pieceInput }
    Combine = MapATACWithHiCWithPRO(pieceInput)
    adjZ = AdjustZInflation(Combine.ToBacon)
    calMoDIFY(adjZ)
}

