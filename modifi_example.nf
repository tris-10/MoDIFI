nextflow.enable.dsl=2

process extractPeaks {
    tag "$cell"

    container 'modifi.sif'

    publishDir "${params.temp_dir}", mode: 'copy', saveAs: { "${cell}.narrowPeak" }

    input:
    tuple val(cell), path(BedPath)

    output:
    path "${cell}.narrowPeak", emit: peaks, optional:true

    script:
    """
    zcat ${BedPath} | sed "s/Peak/${cell}_Peak/g" > ${cell}.narrowPeak
    """
}


process ConcatenateAndSortPeaks {
    container 'modifi.sif'
    
    publishDir "${params.temp_dir}", mode: 'copy', saveAs: { filename -> filename }

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
    publishDir "${params.temp_dir}", mode: 'copy', saveAs: { filename -> "${cell}_R${RepNum}_counts.txt" }

    input:
    tuple path(merge_bed), val(cell), path(BAM_FILE), val(RepNum)

    output:
    path "${cell}_R${RepNum}_counts.txt", emit: ATAC_counts, optional:true

    script:
    """
    set -euxo pipefail
    echo "Checking BAM files in ${BAM_FILE}" | tee -a ${cell}_R${RepNum}_process.log
    OUTPUT_FILE="${cell}_R${RepNum}_counts.txt"

    if [ -f "${BAM_FILE}" ]; then
    	if [ ! -f "${BAM_FILE}.bai" ]; then
        	echo "Index file missing for ${BAM_FILE}, generating now..." | tee -a ${cell}_R${RepNum}_process.log
                samtools index ${BAM_FILE}
        fi
        python ${params.script_dir}/atac_seq_counts.py ${merge_bed} ${BAM_FILE} \${OUTPUT_FILE} ${cell}_R${RepNum} --minQ ${params.atac_minQ}
    else
    	echo "BAM file ${BAM_FILE} not found" | tee -a ${cell}_R${RepNum}_process.log
    fi
    """
}


process GenerateDESeq2InputsAnnCOUNTs{
    container 'modifi.sif'
    publishDir "${params.temp_dir}", mode: 'copy', saveAs: { filename -> filename }

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
    publishDir "${params.temp_dir}", mode: 'copy', saveAs: { filename -> filename }

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
    //publishDir "${params.output_dir}", mode: 'copy', pattern:"*_vs_*_${data_type}.txt", saveAs: { filename -> filename }
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"*_vs_*_${data_type}.txt", saveAs: { filename -> filename }

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
    //publishDir "${params.output_dir}", mode: 'copy', pattern:"*_vs_*_${data_type}.txt", saveAs: { filename -> filename }
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"*_vs_*_${data_type}.txt", saveAs: { filename -> filename }


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
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"RNA_*.txt", saveAs: { filename -> filename }

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

process getDEseqInputForRNAseqCONDS{
    container 'modifi.sif'
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"RNA_*.txt", saveAs: { filename -> filename }

    output:
    path "RNA_conds.txt", emit: RNA_conds, optional:true

    script:
    """
    python ${params.script_dir}/GenerateDESeq2InputsForRNAseqConds.py \\
        -i "${params.resources_dir}/" \\
        -l "${params.RNACountFile.keySet()}" \\
        -f "${params.RNACountFile.values()}" 
    """
}


process generateSampleInfo{
    container 'modifi.sif'
    publishDir "${params.resources_dir}", mode: 'copy', pattern:"Sample*.tsv", saveAs: { filename -> filename }
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"*piece*.tsv", saveAs: { filename -> filename }
    //publishDir "${params.temp_dir}", mode: 'copy', pattern:"*_MergedLoop*.tsv", saveAs: { filename -> filename }

    input:
    path DEseq_atac
    path DEseq_rna

    output:
    path "SamplePair*tsv", emit: SamplePair, optional:true
    path "SampleInfo*tsv", emit: SampleInfo, optional: true
    path "piece*.tsv", emit: Piece, optional: true    
    path "*_MergedLoop*.tsv", emit: HiCPath, optional: true

    script:
    def samplePairArg = params.SamplePair ? "\"${params.SamplePair}\"" : "None"
    def sampleInfoArg = params.SampleInfo ? "\"${params.SampleInfo}\"" : "None"
    """
    python ${params.script_dir}/CheckSamplesForMoDIFI.py \\
	-i "${params.resources_dir}/" \\
	-o "${params.temp_dir}/" \\
	-l "${params.HiCLoopsFile.keySet()}" \\
	-c "${params.HiCLoopsFile.values()}" \\
	--LABEL_RNA "${params.RNASeq}" \\
	--LABEL_ATAC "${params.ATACSeq}" \\
        -sp ${samplePairArg} \\
        -si ${sampleInfoArg}
    """
}

process generateSampleInfo_noInput{
    container 'modifi.sif'
    publishDir "${params.resources_dir}", mode: 'copy', pattern:"Sample*.tsv", saveAs: { filename -> filename }
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"*piece*.tsv", saveAs: { filename -> filename }
    //publishDir "${params.temp_dir}", mode: 'copy', pattern:"*_MergedLoop*.tsv", saveAs: { filename -> filename }

   // input:
   // path DEseq_atac
   // path DEseq_rna

    output:
    path "SamplePair*tsv", emit: SamplePair, optional:true
    path "SampleInfo*tsv", emit: SampleInfo, optional: true
    path "piece*.tsv", emit: Piece, optional: true
    path "*_MergedLoop*.tsv", emit: HiCPath, optional: true

    script:
    def samplePairArg = params.SamplePair ? "\"${params.SamplePair}\"" : "None"
    def sampleInfoArg = params.SampleInfo ? "\"${params.SampleInfo}\"" : "None"
    """
    python ${params.script_dir}/CheckSamplesForMoDIFI.py \\
        -i "${params.resources_dir}/" \\
        -o "${params.temp_dir}/" \\
        -l "${params.HiCLoopsFile.keySet()}" \\
        -c "${params.HiCLoopsFile.values()}" \\
        --LABEL_RNA "${params.RNASeq}" \\
        --LABEL_ATAC "${params.ATACSeq}" \\
	-sp ${samplePairArg} \\
	-si ${sampleInfoArg}
    """
}


process MapATACWithPRO{
    container 'modifi.sif'
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"PRO_ATAC*.tsv", saveAs: { filename -> filename }

    input:
    tuple path(piece_file), val(chrom)

    output:
    path "PRO_ATAC*tsv", emit: PRO_ATAC, optional:true

    script:
    def sampleInfoArg = params.SampleInfo ? "\"${params.SampleInfo}\"" : "\"${params.resources_dir}/SampleInfo.tsv\""
    """
    echo "Using piece file: ${piece_file}"
    python ${params.script_dir}/MapATACWithPROPerChrom.py \\
        -i "${params.resources_dir}/" \\
        -o "${params.temp_dir}/" \\
        --LABEL_RNA "${params.RNASeq}" \\
        --LABEL_ATAC "${params.ATACSeq}" \\
        -si "${sampleInfoArg}" \\
        -sp "${piece_file}" \\
	-c ${chrom} \\
        --PRO_Region 5000 \\
        --PRO_minOL 0.5
    """
}

process combinPRO_ATAC{
    container 'modifi.sif'
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"PRO_ATAC*.tsv", saveAs: { filename -> filename }

    input:
    path PRO_ATAC

    output:
    path "PRO_ATAC*tsv", emit: PRO_ATAC_COM, optional:true

    script:
    """
    echo "Files are combining: ${PRO_ATAC}"
    python ${params.script_dir}/CombineFiles.py \\
        -i "${PRO_ATAC}" \\
        -o "${params.temp_dir}/" 
    """
}

process MapATACWithHiCWithPRO{
    container 'modifi.sif'
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"ToBacon_*.tsv", saveAs: { filename -> filename }
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"atacWithHiC*.tsv", saveAs: { filename -> filename }
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"PRO_ATAC*.tsv", saveAs: { filename -> filename }

    input:
    tuple path(piece_file), path(ALL_PRO_ATAC), path(hicpathInput)

    output:
    path "ToBacon_*tsv", emit: ToBacon, optional:true
    path "atacWithHiC*tsv", emit: HiC, optional:true
    //path "PRO_ATAC*tsv", emit: PRO_ATAC, optional:true 

    script:
    def sampleInfoArg = params.SampleInfo ? "\"${params.SampleInfo}\"" : "\"${params.resources_dir}/SampleInfo.tsv\""
    """
    echo "Using piece file: ${piece_file}"
    python ${params.script_dir}/MapATACWithHiCWithPROperLoopFile.py \\
        -i "${params.resources_dir}/" \\
        -o "${params.temp_dir}/" \\
        --LABEL_RNA "${params.RNASeq}" \\
        --LABEL_ATAC "${params.ATACSeq}" \\
        -si "${sampleInfoArg}" \\
        -sp "${piece_file}" \\
        --PRO_Region 5000 \\
        --PRO_minOL 0.5 \\
        --PRO_ATAC "${ALL_PRO_ATAC}" \\
        --hicloop ${hicpathInput}
    """
}

process combinAdjZFile{
    container 'modifi.sif'
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"ToBacon_*.tsv", saveAs: { filename -> filename }

    input:
    tuple path(piece_file), path(ToBacon)

    output:
    path "ToBacon_*tsv", emit: ToBaconFinal, optional:true

    script:
    """
    echo "ToBacon files are combing based on ${piece_file}"
    python ${params.script_dir}/CombineFiles.py \\
        -i "${params.temp_dir}/" \\
        -o "${params.temp_dir}/" \\
        -sp "${piece_file}" \\
        -p "bacon"
    """
}

process AdjustZInflation{
    container 'modifi.sif'
    publishDir "${params.temp_dir}", mode: 'copy', pattern:"ToBacon*bacon*.tsv", saveAs: { filename -> filename }

    input:
    path ToBaconFinal 

    output:
    path "ToBacon*bacon*tsv", emit: adjZ, optional:true

    script:
    """
    Rscript ${params.script_dir}/baconForDEseq.R ${ToBaconFinal}
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
        -i "${params.temp_dir}/${adjZ}" \\
        -o "${params.temp_dir}/" \\
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
workflow test{
	//cellChannel = Channel.fromList(params.ATACBAMFiles_B.values())
	CellBAMM = Channel.fromList(
    	params.ATACBAMFiles_B
        .collectMany { name, bams ->
            def sorted = bams.toList().sort { it.toString() }
            sorted.withIndex().collect { bam, idx ->
                tuple(name, file(bam, checkIfExists: true), idx + 1)
            }
        }
)
	CellBAMM.view { item ->
                println "Cell: ${item[0]}, Path: ${item[1]}, N: ${item[2]}"
        }
}



workflow tt {
	
	CellBAMM = Channel.fromList(
  	params.ATACBAMFiles.collectMany { cell, bams ->
    	def list = (bams instanceof Collection) ? bams : [ bams ]
    	list.withIndex().collect { b, i -> tuple(cell, file(b, checkIfExists: true), i + 1) }
  	}
	)
	CellBAMM.view { item ->
                println "Cell: ${item[0]}, Path: ${item[1]}, N: ${item[2]}"
        }

}

workflow ATACseq_workflow{
        // Processing ATAC-seq signals and get three inputfiles for DEseq
	cellChannel = Channel.from(params.ATACBEDFile.keySet())
	rep = Channel.fromList(params.ATACBEDFile.values())
        CellRep = Channel
        .fromList(params.ATACBEDFile.collect { cell, bed ->
            tuple(cell, file(bed, checkIfExists: true))
        })
	Cell = Channel.from(params.ATACBAMFiles.keySet())
	BAM = Channel.from(params.ATACBAMFiles.values())
	CellBAMM = Channel.fromList(params.ATACBAMFiles.collect { name, bams ->
        tuple(name, file(bams, checkIfExists: true))})
    	.flatMap { name, bams ->
        def sorted = bams.sort { it.toString() }
	sorted.withIndex().collect { b, i -> tuple(name, file(b, checkIfExists: true), i + 1) }
    	}
	Peaks = extractPeaks(CellRep).peaks.collect()
	Merge_bed = ConcatenateAndSortPeaks(Peaks).collect()
	crossProductChannel = Merge_bed.combine(CellBAMM)
        crossProductChannel.view { item ->
                println "Bed File: ${item[0]}, Cell: ${item[1]}, Path: ${item[2]},  RepNum: ${item[3]}"
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
    DEseqInputRNA = getDEseqInputForRNAseqCONDS()
    DEseq_rna = runDEseqRSEM(DEseqInputRNA.RNA_conds, params.RNA_ANN_File, params.RNA_Key_col, params.RNASeq) 
    DEseq_rna.view()
    
    emit:
    DEseq_rna
}


workflow ttR{
    InfoPair = generateSampleInfo_noInput()
    InfoPair.Piece.flatten()
        .view { item ->
            println "File: ${item}"
        }
        .set { pieceInput }
    InfoPair.HiCPath.flatten()
        .view { item ->
           // println "HiC bedFile: ${item}"
        }
        .set { hicpathInput }

    indexChannel = Channel.of(1..24)
    PRO_INPUT = pieceInput.combine(indexChannel)
    .map{x->
	tuple(x[0], x[1])
	}
    PRO_ATAC = MapATACWithPRO(PRO_INPUT)
    PRO_ATAC_COM = combinPRO_ATAC(PRO_ATAC.collect())
}

workflow recalMoDIFI{
    InfoPair = generateSampleInfo_noInput()
    InfoPair.Piece.flatten()
        .view { item ->
            println "File: ${item}"
        }
        .set { pieceInput }
    InfoPair.HiCPath.flatten()
        .view { item ->
           // println "HiC bedFile: ${item}"
        }
        .set { hicpathInput }
    indexChannel = Channel.of(1..24)
    PRO_INPUT = pieceInput.combine(indexChannel)
    .map{x->
        tuple(x[0], x[1])
        }
    PRO_ATAC = MapATACWithPRO(PRO_INPUT)
    PRO_ATAC_COM = combinPRO_ATAC(PRO_ATAC.collect())

    ALL_PRO_ATAC = PRO_ATAC_COM.collect()
    ALL_PRO_ATAC.view()
    COMBINE_INPUT = pieceInput
    .combine(ALL_PRO_ATAC)
    .combine(hicpathInput)
    .map { x ->
        //println "DEBUG: ${x}"
        //println "DEBUG CLASS: ${x.getClass()}"
        def piece = x[0]
        def hicpath = x[-1]
        def all_pro_atac = x[1..-2]
        tuple(piece, all_pro_atac, hicpath)
    }
    CombineToBacon = MapATACWithHiCWithPRO(COMBINE_INPUT)
    AllCombineToBacon = CombineToBacon.ToBacon.collect()
    ALLHIC_INPUT = pieceInput
    .combine(AllCombineToBacon)
    .map {x->
        tuple(x[0], x[1..-1])
    }
    ToBaconFinal = combinAdjZFile(ALLHIC_INPUT)
    adjZ = AdjustZInflation(ToBaconFinal)
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
    InfoPair.HiCPath.flatten()
        .view { item ->
           // println "HiC bedFile: ${item}"
        }
        .set { hicpathInput }
    indexChannel = Channel.of(1..24)
    PRO_INPUT = pieceInput.combine(indexChannel)
    .map{x->
        tuple(x[0], x[1])
        }
    PRO_ATAC = MapATACWithPRO(PRO_INPUT)
    PRO_ATAC_COM = combinPRO_ATAC(PRO_ATAC.collect())

    ALL_PRO_ATAC = PRO_ATAC_COM.collect()
    ALL_PRO_ATAC.view()
    COMBINE_INPUT = pieceInput
    .combine(ALL_PRO_ATAC)
    .combine(hicpathInput)
    .map { x ->
       // println "DEBUG: ${x}"
       // println "DEBUG CLASS: ${x.getClass()}"
        def piece = x[0]
        def hicpath = x[-1]
        def all_pro_atac = x[1..-2]
        tuple(piece, all_pro_atac, hicpath)
    }
    CombineToBacon = MapATACWithHiCWithPRO(COMBINE_INPUT) 

    AllCombineToBacon = CombineToBacon.ToBacon.collect()
    ALLHIC_INPUT = pieceInput
    .combine(AllCombineToBacon)
    .map {x->
        tuple(x[0], x[1..-1])
    }
    ToBaconFinal = combinAdjZFile(ALLHIC_INPUT)
    adjZ = AdjustZInflation(ToBaconFinal)
    calMoDIFY(adjZ)
}

