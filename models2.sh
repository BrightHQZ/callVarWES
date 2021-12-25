#!/bin/bash

function error_exit() {
  local dateStr=$(date +"%Y-%m-%d %H:%M:%S");
  echo -e $dateStr":\t$2" >> $1_error.log;
  exit 1
}

function write_log() {
    local dateStr=$(date +"%Y-%m-%d %H:%M:%S");
    local wTime=$(echo $(($(date +%s -d "$dateStr") - $(date +%s -d "$3 $4"))) | awk '{t=split("60 s 60 m 24 h 999 d",a);for(n=1;n<t;n+=2){if($1==0)break;s=$1%a[n]a[n+1]s;$1=int($1/a[n])}printf s}')
    #$3 is date $4 is time
    echo -e $dateStr":\t$wTime\t$2" >> $1.log;
}

function checkFamily() {
    local inFile1=$1;
    local inFile2=$2;
    local inSample=$3;
    local infileType=$4;
    local res="";
    if [ -z $inFile1 ]; then
        res="Error: -1 is necessary!";
    else
        inFileA=($(echo $inFile1 | tr "," "\n"));
        inSamples=($(echo $inSample | tr "," "\n"));
        if [[ "$infileType" == "bam" || "$infileType" == "vcf" ]]; then
            if [ ${#inFileA[@]} == ${#inSamples[@]} ]; then
                res=${infileType^^}
            else
                res="Error: the file and sample count should equaly in -1 and sampleID!";
            fi;
        elif [ "$infileType" == "fq" ]; then
            if [ "$inFile2" == "|" ] && [ ${#inFileA[@]} == ${#inSamples[@]} ]; then
                res="SE";
            else
                inFileB=($(echo $inFile2 | tr "," "\n"));
                if [ ${#inFileA[@]} == ${#inFileB[@]} ] && [ ${#inFileA[@]} == ${#inSamples[@]} ]; then
                    res="PE"
                else
                    res="Error: the file and sample count should equaly in -1, -2 and sampleID!";
                fi;
            fi;
        fi;
    fi;
    echo $res;
}

function qualityFASTQ_SE() {
    local SOAPnuke=$1;
    local inFile=$2;
    local outDir=$3;
    local sample=$4;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$SOAPnuke filter -n 0.1 -q 0.1 -l 5 -G 2 --seqType 1 -Q 2 -1 $inFile -o $outDir -C $sample.clean.fq.gz -T 10";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$sample.clean.fq.gz";
}

function qualityFASTQ_PE() {
    local SOAPnuke=$1;
    local inFile1=$2;
    local inFile2=$3;
    local outDir=$4;
    local sample=$5;
    local outLog=$6;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$SOAPnuke filter -n 0.1 -q 0.1 -l 5 -G 2 --seqType 1 -Q 2 -1 $inFile1 -2 $inFile2 -o $outDir -C $sample.1.clean.fq.gz -D $sample.2.clean.fq.gz -T 10" ;
    printf "$command\n" >&2;
    $command || error_exit "$outLog/$sample" "$sample\t$command failed";
    write_log "$outLog/$sample" "$sample\t$command" $pStart;
    echo "OK";
}


function aligmentFQ_SE() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local ref=$4;
    local bwa=$5;
    local mapScore=$6
    local samtools=$7;
    local outLog=$8;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";

    command="$bwa mem -a -T $mapScore -t 5 -M -R '@RG\tID:$sample\tPL:illumina\tPU:$sample\tLB:$sample\tSM:$sample\tCN:SZSL' $ref $outDir/$inFile | $samtools view -F 4 -F 256 -q 10 -S -b -o $outLog/$sample.bam";
    printf "$command\n" >&2;
    #mapping reads to reference genetic and then removing unmapped reads and low quality reads 
    #(-F 4 (read1 unmap): removeing unmapped reads; -q 10: the mapped quality lower than 10) 
    #-F 256 (Filter mul-aligment reads)
    bwa mem -a -T $mapScore -t 5 -M -R "@RG\tID:$sample\tPL:illumina\tPU:$sample\tLB:$sample\tSM:$sample\tCN:SZSL" $ref $outDir/$inFile | $samtools view -F 4 -F 256 -q 10 -S -b -o $outDir/$sample.bam || error_exit "$outLog/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    pStart=$(date +"%Y-%m-%d %H:%M:%S");
    #command="$samtools sort $outDir/$sample.bam -T $outDir/tmp -@ 5 -o $outDir/$sample.sort.bam";
    #$command || error_exit "$outDir/$sample" "$sample\t$command failed";
    #printf "$command\n" >&2;
    #write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$sample.bam";
}

function aligmentFQ_PE() {
    local inFile1=$1;
    local inFile2=$2;
    local outDir=$3;
    local sample=$4;
    local ref=$5;
    local bwa=$6;
    local mapScore=$7
    local samtools=$8;
    local outLog=$9;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";

    #Aligment command.
    #command="$bwa mem -a -T $mapScore -t 5 -M -R '@RG\tID:$sample\tPL:illumina\tPU:$sample\tLB:$sample\tSM:$sample\tCN:SZSL' $ref $inFile1 $inFile2 | samtools view -S -b -o $outDir/$sample.bam";
    command="$bwa mem -a -T $mapScore -t 10 -M -R '@RG\tID:$sample\tPL:illumina\tPU:$sample\tLB:$sample\tSM:$sample\tCN:SZSL' $ref $inFile1 $inFile2 | $samtools view -F 4 -F 8 -F 12 -F 256 -q 10 -S -b -o $outDir/$sample.bam";
    printf "$command\n" >&2;
    #mapping reads to reference genetic and then removing unmapped reads and low quality reads 
    #(-F 4 (read1 unmap) -F 8 (read2 unmap) -F 12 (both unmap): removeing unmapped reads; -q 10: the mapped quality lower than 10) 
    #-F 256 (Filter mul-aligment reads)
    bwa mem -a -T $mapScore -t 10 -M -R "@RG\tID:$sample\tPL:illumina\tPU:$sample\tLB:$sample\tSM:$sample\tCN:SZSL" $ref $inFile1 $inFile2 | $samtools view -F 4 -F 8 -F 12 -F 256 -q 10 -S -b -o $outDir/$sample.bam || error_exit "$outLog/$sample" "$sample\t$command failed";
    write_log "$outLog/$sample" "$sample\t$command" $pStart;
    #pStart=$(date +"%Y-%m-%d %H:%M:%S");
    #command="$samtools sort $outDir/$sample.bam -T $outDir/tmp -@ 5 -o $outDir/$sample.sort.bam";
    #printf "$command\n" >&2;
    #$command || error_exit "$outLog/$sample" "$sample\t$command failed";
    #write_log "$outLog/$sample" "$sample\t$command" $pStart;
    echo "$sample.bam";
}

function sortBam() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local samtools=$4;
    local outLog=$5;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$samtools sort $inFile -T $outDir/tmp -@ 5 -o $outDir/$sample.sort.bam";
    printf "$command\n" >&2;
    $command || error_exit "$outLog/$sample" "$sample\t$command failed";
    write_log "$outLog/$sample" "$sample\t$command" $pStart;
    echo "$sample.sort.bam";
}

function markDup() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local sambamba=$4;
    local outLog=$5;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$sambamba markdup $inFile $outDir/$sample.sort.markdup.bam -t 5 --tmpdir=$outDir/tmp";
    printf "$command\n" >&2;
    $command || error_exit "$outLog/$sample" "$sample\t$command failed";
    write_log "$outLog/$sample" "$sample\t$command" $pStart;
    echo "$sample.sort.markdup.bam";
}

#Quality bam
function qualityBAM() {
    local bamdst=$1;
    local inFile=$2;
    local outDir=$3
    local bedFile=$4;
    local command="";
    command="$bamdst -p $bedFile -o $outDir $inFile";
    printf "$command\n" >&2;
    $command || error_exit "$outLog/$sample" "$sample\t$command failed";
} 

function BaseRecalibrator() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local GAKT=$4;
    local ref=$5;
    local vcf1=$6;
    local vcf2=$7;
    local vcf3=$8;
    local outLog=$9;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx1g -jar $GAKT BaseRecalibrator -R $ref -I $inFile -O $outDir/$sample.table --known-sites $vcf1 --known-sites $vcf2 --known-sites $vcf3";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    pStart=$(date +"%Y-%m-%d %H:%M:%S");
    command="java -Xmx1g -jar $GAKT ApplyBQSR -R $ref -I $inFile -O $outDir/$sample.sort.markdup.recal.bam --bqsr-recal-file $outDir/$sample.table";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$sample.sort.markdup.recal.bam";
}


function CalibrateDragstrModel() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local GATK=$4;
    local inDrageStr=$5;
    local outLog=$6;
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar $GATK CalibrateDragstrModel";
    command="$command -I $inFile ";
    command="$command -O $outDir/$sample.drags"
    command="$command -R $ref -str $inDrageStr" ;
    printf "$command\n" >&2;
    $command || error_exit $outLog "$sample\t$command failed";
    write_log $outLog "$sample\t$command" $pStart;
    echo "$outDir/$sample.drags";
}

function callSNV_chr_germline_gvcf() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local chr=$4;
    local intervals=$5;
    local GATK=$6;
    local ref=$7;
    local dbSNP=$8;
    local inDrageFile=$9;
    local outLog=${10};
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar $GATK HaplotypeCaller";
    command="$command -I $inFile ";
    command="$command -O $outDir/$sample.$chr.g.vcf.gz"
    command="$command -R $ref -D $dbSNP " ;
    command="$command -ERC GVCF -L $intervals --dragstr-params-path $inDrageFile" ;
    printf "$command\n" >&2;
    $command || error_exit $outLog "$sample\t$command failed";
    write_log $outLog "$sample\t$command" $pStart;
    echo "$outDir/$sample.$chr.g.vcf.gz";
}

function GatherVcfs () {
    local outDir=$1;
    local sample=$2;
    local GATK=$3;
    local outLog=$4;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar $GATK GatherVcfs ";
    command="$command -I $outDir/$sample.chr1.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr2.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr3.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr4.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr5.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr6.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr7.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr8.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr9.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr10.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr11.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr12.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr13.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr14.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr15.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr16.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr17.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr18.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr19.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr20.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr21.g.vcf.gz ";
    command="$command -I $outDir/$sample.chr22.g.vcf.gz ";
    command="$command -I $outDir/$sample.chrX.g.vcf.gz ";
    command="$command -I $outDir/$sample.chrY.g.vcf.gz ";
    command="$command -O $outDir/$sample.g.vcf.gz ";
    printf "$command\n" >&2;
    $command || error_exit $outLog "$sample\t$command failed";
    write_log $outLog "$sample\t$command" $pStart;
    #command="tabix -p vcf $outDir/$sample.g.vcf.gz";
    #printf "$command\n" >&2;
    #$command || error_exit "$outLog/$sample" "$sample\t$command failed";
    #write_log "$outLog/$sample" "$sample\t$command" $pStart;
    #echo "\n$sample.g.vcf.gz";
}

function IndexFeatureFile () {
    local inFile=$1;
    local GATK=$2;
    local outLog=$3;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar $GATK IndexFeatureFile ";
    command="$command -I $inFile ";
    printf "$command\n" >&2;
    $command || error_exit $outLog "$sample\t$command failed";
    write_log $outLog "$sample\t$command" $pStart;
}

function CombineGVCFs() {
    local outDir=$1;
    local inFile=$2;
    local sample=$3;
    local GATK=$4;
    local ref=$5;
    local outLog=$6;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar $GATK CombineGVCFs ";
    inFiles=($(echo $inFile | tr "," "\n"));
    for index in "${!inFiles[@]}"; do
        command="$command -V ${inFiles[index]} ";
    done;
    command="$command -R $ref -O $outDir/$sample.g.vcf.gz ";
    printf "$command\n" >&2;
    $command || error_exit error_exit $outLog "$sample\t$command failed";
    write_log $outLog "$sample\t$command" $pStart;
    #command="tabix -p vcf $outDir/$sample.g.vcf.gz";
    #printf "$command\n" >&2;
    #$command || error_exit "$outLog/$sample" "$sample\t$command failed";
    #write_log "$outLog/$sample" "$sample\t$command" $pStart;
    echo "$sample.g.vcf.gz";
}

function gvcfTOvcf {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local GATK=$4;
    local ref=$5;
    local outLog=$6;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar $GATK GenotypeGVCFs";
    command="$command -V $inFile"
    command="$command -O $outDir/$sample.vcf.gz"
    command="$command -R $ref";
    printf "$command\n" >&2;
    $command || error_exit $outLog "$sample\t$command failed";
    write_log $outLog "$sample\t$command" $pStart;
    echo "$sample.vcf.gz";
}


function anntationSNP() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local ver=$4;
    local annovar=$5;
    local outLog=$6;
    local refGene=$7;
    local cytoBand=$8;
    local gnomad211=$9;
    local clinvar=${10};
    local avsnp=${11};
    local dbscsnv=${12};
    command="";
    if [ $ver == "hg38" ]; then
        command="perl $annovar/table_annovar.pl $inFile $annovar/humandb/ -buildver $ver -out $outDir/$sample.annovar ";
        command="$command -protocol $refGene,$cytoBand,$gnomad211,$clinvar,$avsnp,$dbscsnv -operation g,r,f,f,f,f -nastring . -vcfinput";
    elif [ $ver == "hg19" ]; then 
        command="perl $annovar/table_annovar.pl $inFile $annovar/humandb/ -buildver $ver -out $outDir/$sample.annovar ";
        command="$command -protocol $refGene,$cytoBand,$gnomad211,$clinvar,$avsnp,$dbscsnv -operation g,r,f,f,f,f -nastring . -vcfinput";
    fi;
    printf "$command\n" >&2;
    $command || error_exit "$outLog/$sample" "$sample\t$command failed";
    write_log "$outLog/$sample" "$sample\t$command" $pStart;
    echo "$sample.annovar";
}

function clearTempFile {
    local outDir=$1;
    local sample=$2;
    local outLog=$3;
    local inVer=$4;
    local clearRoot="";
    local command="";
    
    #Clean the direct of bam
    clearRoot="$outDir/bam/$sample";
    command="rm -rf $clearRoot/$sample.bam $clearRoot/$sample.sort.bam $clearRoot/tmp"
    printf "$command\n" >&2;
    $command || error_exit $outLog "$sample\t$command failed";
    #Clean the direct of gvcf
    clearRoot="$outDir/gvcf/$sample";
    command="rm $clearRoot/$sample.hg19.chr*"
    printf "$command\n" >&2;
    $command || error_exit $outLog "$sample\t$command failed";
    #Clean the direct of vcf
    clearRoot="$outDir/vcf/$sample";
    command="rm $clearRoot/$sample.annovar.avinput $clearRoot/$sample.$inVer""_avsnp150* $clearRoot/$sample.$inVer""_gnom*"
    command="$command $clearRoot/$sample.$inVer""_clinvar* $clearRoot/$sample.$inVer""_dbscsnv* $clearRoot/$sample.$inVer""_cyto*"
    command="$command $clearRoot/$sample.annovar.refGen*"
    printf "$command\n" >&2;
    $command || error_exit $outLog "$sample\t$command failed";
    echo "Clearnning temp file finished!";
}

###################################################################################################################



function blastn() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local blastn=$4;
    local blastdb=$5;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$blastn -query $outDir/$sample.fa -db $blastdb -perc_identity 97  -evalue 1e-5 -outfmt \"6 qseqid sseqid pident nident length mismatch positive gapopen gaps ppos qframe sframe sstrand qcovs qstart qend qseq sstart send sseq evalue bitscore score\" -max_hsps 1 |";
    command="$command awk '{if(\$5>32){print \$0}}' > $outDir/$sample.blast_tab"
    printf "$command\n" >&2;       
    $blastn -query $outDir/$sample.fa -db $blastdb -perc_identity 97  -evalue 1e-5 -outfmt "6 qseqid sseqid pident nident length mismatch positive gapopen gaps ppos qframe sframe sstrand qcovs qstart qend qseq sstart send sseq evalue bitscore score" -max_hsps 1 | awk '{if($5>32){print $0}}' > $outDir/$sample.blast_tab || error_exit "$outDir/$sample" "$sample\t$command failed";
    #$command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$outDir/$sample.blast_tab";
}



function mappedReads() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local samtools=$4;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$samtools view -Sb -F 1284 -@ 5 -q 1 $outDir/$inFile -o $outDir/$sample.mapped.bam";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    pStart=$(date +"%Y-%m-%d %H:%M:%S");
    command="$samtools index $outDir/$sample.mapped.bam";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$sample.mapped.bam";
}

function un-mappedReads() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local samtools=$4;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$samtools view -Sb -f 4 -@ 5 $outDir/$inFile -o $outDir/unmapped_$sample.bam";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    #samtools view -Sb -f 4 -@ 5 -q 1 $outDir/$inFile -o $outDir/$sample.unmapped.bam || error_exit $sample "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "unmapped_$sample.bam";
}

function bamTofq() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local samtools=$4;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$samtools fastq -0 $outDir/$sample.fq.gz -c 9 $outDir/$inFile";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$sample.fq.gz";
}

function bamTofa() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local samtools=$4;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="";
    command="$samtools fasta -0 $outDir/$sample.fa $outDir/$inFile";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$sample.fa";
}



# function BaseRecalibrator() {
#     local inFile=$1;
#     local outDir=$2;
#     local sample=$3;
#     local GAKT=$4;
#     local ref=$5;
#     local vcf1=$6;
#     local vcf2=$7;
#     local vcf3=$8;
#     local sample2=$9;
#     local pStart=$(date +"%Y-%m-%d %H:%M:%S");
#     local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx1g -jar $GAKT BaseRecalibrator -R $ref -I $outDir/$inFile -O $outDir/$sample.table --known-sites $vcf1 --known-sites $vcf2 --known-sites $vcf3";
#     printf "$command\n" >&2;
#     $command || error_exit "$outDir/$sample2" "$sample2\t$command failed";
#     write_log "$outDir/$sample2" "$sample2\t$command" $pStart;
#     pStart=$(date +"%Y-%m-%d %H:%M:%S");
#     command="java -Xmx1g -jar $GAKT ApplyBQSR -R $ref -I $outDir/$inFile -O $outDir/$sample.mapped.recal.bam --bqsr-recal-file $outDir/$sample.table";
#     printf "$command\n" >&2;
#     $command || error_exit "$outDir/$sample2" "$sample2\t$command failed";
#     write_log "$outDir/$sample2" "$sample2\t$command" $pStart;
#     echo "$sample.mapped.recal.bam";
# }

function callSNV_callSNV_BCFTOOLS() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local ref=$4;
    local bcftools=$5;
    local outLog=$6;
    local command="$bcftools mpileup -f $ref -a FORMAT/DP4 $inFile | $bcftools call -mv | $bcftools filter -i 'QUAL>20 && DP>10 && DP4[2]/(DP4[2]+DP4[3])<0.8 && DP4[2]/(DP4[2]+DP4[3])>0.2' -O v -o $outDir/$sample.bcftools.vcf";
    printf "$command\n" >&2;
    $bcftools mpileup -f $ref -a FORMAT/DP4 $inFile | $bcftools call -mv | $bcftools filter -i 'QUAL>20 && DP>10 && DP4[2]/(DP4[2]+DP4[3])<0.8 && DP4[2]/(DP4[2]+DP4[3])>0.2' -O v -o $outDir/$sample.bcftools.vcf || error_exit "$outLog/$sample" "$sample\t$command failed";
    write_log "$outLog/$sample" "$sample\t$command" $pStart;
    # command="$bcftools sort -O v -o $outDir/$sample.bcftools.sorted.vcf" 
    # printf "$command\n" >&2;
    # $command || error_exit "$outLog/$sample" "$sample\t$command failed";
    # write_log "$outLog/$sample" "$sample\t$command" $pStart;
    echo "$sample.bcftools.vcf";
}



function callSNV_chr_germline_vcf() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local chr=$4;
    local GATK=$5;
    local ref=$6;
    local dbSNP=$7;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    if [ $chr -eq 24 ]; then
        chr="chrY";
    elif [ $chr -eq 23 ]; then
        chr="chrX";
    else
        chr="chr$i";
    fi;
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx4g -jar $GATK HaplotypeCaller";
    command="$command -I $inFile"
    command="$command -O $outDir/$chr.vcf.gz -stand-call-conf 20"
    command="$command -R $ref -D $dbSNP";
    command="$command -L $chr";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$outDir/$chr.vcf.gz";
}


function callSNV_germline() {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local GATK=$4;
    local ref=$5;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx4g -jar $GATK HaplotypeCaller";
    command="$command -I $outDir/$sample.mapped.recal.bam"
    command="$command -O $outDir/$sample.g.vcf.gz"
    command="$command -R $ref";
    command="$command -ERC GVCF";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$outDir/$sample.g.vcf.gz";
}

function gvcfTOdb {
    local inFile=$1;
    local outDir=$2;
    local sample=$3;
    local GATK=$4;
    local dbPath=$5;
    local intervals=$6
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar -Xmx128m -Xms64m $GATK GenomicsDBImport";
    command="$command -V $inFile -L $intervals"
    if [ ! -d $dbPath ]; then
        command="$command --genomicsdb-workspace-path $dbPath"
    else 
        command="$command --genomicsdb-update-workspace-path $dbPath"
    fi;
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$outDir/$sample.g.vcf.gz";
}


function gvcfTOdb_M {
    local inFile=$1;
    local outDir=$2;
    local GATK=$3;
    local dbPath=$4;
    local intervals=$5
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar -Xmx128m -Xms64m $GATK GenomicsDBImport";
    command="$command $inFile -L $intervals --tmp-dir=$outDir"
    if [ ! -d $dbPath ]; then
        command="$command --genomicsdb-workspace-path $dbPath"
    else 
        command="$command --genomicsdb-update-workspace-path $dbPath"
    fi;
    printf "$command\n" >&2;
    $command #|| error_exit "$outDir/$sample" "$sample\t$command failed";
    #write_log "$outDir/$sample" "$sample\t$command" $pStart;
}



function VariantRecalibrator() {
    local inF=$1;
    local outDir=$2;
    local sample=$3;
    local ref=$4;
    local GATK=$5;
    local vcf1=$6;
    local vcf2=$7;
    local vcf3=$8; 
    local vcf4=$9;
    local mode=$10;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -Xmx4g -jar $GATK VariantRecalibrator";
    command="$command -R $refs";
    command="$command -V $inF";
    command="$command --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $vcf3";
    command="$command --resource:omni,known=false,training=true,truth=false,prior=12.0 $vcf2";
    command="$command --resource:1000G,known=false,training=true,truth=false,prior=10.0 $vcf4";
    command="$command --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $vcf1";
    command="$command -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR";
    command="$command -mode $mode";
    if [ $mode == "SNP" ]; then
        command="$command -O $outDir/$sample.$mode.recal";
        command="$command --tranches-file $outDir/$sample.$mode.tranches";
        command="$command --rscript-file $outDir/$sample.$mode.plots.R";
    else
        command="$command -O $outDir/$sample.SNP_$mode.recal";
        command="$command --tranches-file $outDir/$sample.SNP_$mode.tranches";
        command="$command --rscript-file $outDir/$sample.SNP_$mode.plots.R";
    fi
    printf "$command\n";
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    pStart=$(date +"%Y-%m-%d %H:%M:%S");
    command="java -Xmx4g -jar $GATK ApplyVQSR";
    command="$command -R $ref";
    command="$command -V $inF";
    if  [ $mode == "SNP" ]; then
        command="$command -O $outDir/$sample.$mode.VSQR.vcf.gz";
        command="$command --tranches-file $outDir/$sample.$mode.tranches";
        command="$command --recal-file $outDir/$sample.$mode.recal";
    else
        command="$command -O $outDir/$sample.SNP_$mode.VSQR.vcf.gz";
        command="$command --tranches-file $outDir/$sample.SNP_$mode.tranches";
        command="$command --recal-file $outDir/$sample.SNP_$mode.recal";
    fi
    command="$command --truth-sensitivity-filter-level 99.0";
    command="$command -mode $mode";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    if [ $mode == "SNP" ]; then
        echo "$outDir/$sample.$mode.VSQR.vcf.gz";
    else
        echo "$outDir/$sample.SNP_$mode.VSQR.vcf.gz";
    fi
}



function CombineVCFs_chr() {
    local outDir=$1;
    local sample=$2;
    local bcftools=$3;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="$bcftools concat ";
    command="$command $outDir/chr1.vcf.gz ";
    command="$command $outDir/chr2.vcf.gz ";
    command="$command $outDir/chr3.vcf.gz ";
    command="$command $outDir/chr4.vcf.gz ";
    command="$command $outDir/chr5.vcf.gz ";
    command="$command $outDir/chr6.vcf.gz ";
    command="$command $outDir/chr7.vcf.gz ";
    command="$command $outDir/chr8.vcf.gz ";
    command="$command $outDir/chr9.vcf.gz ";
    command="$command $outDir/chr10.vcf.gz ";
    command="$command $outDir/chr11.vcf.gz ";
    command="$command $outDir/chr12.vcf.gz ";
    command="$command $outDir/chr13.vcf.gz ";
    command="$command $outDir/chr14.vcf.gz ";
    command="$command $outDir/chr15.vcf.gz ";
    command="$command $outDir/chr16.vcf.gz ";
    command="$command $outDir/chr17.vcf.gz ";
    command="$command $outDir/chr18.vcf.gz ";
    command="$command $outDir/chr19.vcf.gz ";
    command="$command $outDir/chr20.vcf.gz ";
    command="$command $outDir/chr21.vcf.gz ";
    command="$command $outDir/chr22.vcf.gz ";
    command="$command $outDir/chrX.vcf.gz ";
    command="$command $outDir/chrY.vcf.gz ";
    command="$command -O z -o $outDir/$sample.vcf.gz ";
    # local command="java -Xmx4g -jar $GATK CombineGVCFs";
    # command="$command -R $ref";
    # command="$command -O $outDir/$sample.g.vcf.gz ";
    # command="$command -V $outDir/chr1.g.vcf.gz ";
    # command="$command -V $outDir/chr2.g.vcf.gz ";
    # command="$command -V $outDir/chr3.g.vcf.gz ";
    # command="$command -V $outDir/chr4.g.vcf.gz ";
    # command="$command -V $outDir/chr5.g.vcf.gz ";
    # command="$command -V $outDir/chr6.g.vcf.gz ";
    # command="$command -V $outDir/chr7.g.vcf.gz ";
    # command="$command -V $outDir/chr8.g.vcf.gz ";
    # command="$command -V $outDir/chr9.g.vcf.gz ";
    # command="$command -V $outDir/chr10.g.vcf.gz ";
    # command="$command -V $outDir/chr11.g.vcf.gz ";
    # command="$command -V $outDir/chr12.g.vcf.gz ";
    # command="$command -V $outDir/chr13.g.vcf.gz ";
    # command="$command -V $outDir/chr14.g.vcf.gz ";
    # command="$command -V $outDir/chr15.g.vcf.gz ";
    # command="$command -V $outDir/chr16.g.vcf.gz ";
    # command="$command -V $outDir/chr17.g.vcf.gz ";
    # command="$command -V $outDir/chr18.g.vcf.gz ";
    # command="$command -V $outDir/chr19.g.vcf.gz ";
    # command="$command -V $outDir/chr20.g.vcf.gz ";
    # command="$command -V $outDir/chr21.g.vcf.gz ";
    # command="$command -V $outDir/chr22.g.vcf.gz ";
    # command="$command -V $outDir/chrX.g.vcf.gz ";
    # command="$command -V $outDir/chrY.g.vcf.gz ";
    printf "$command\n" >&2;
    $command || error_exit error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    rm $outDir/chr*.vcf*;
    command="tabix -p vcf $outDir/$sample.vcf.gz";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$outDir/$sample.vcf.gz";
}

function VCF_FormDB_chr() {
    local outDir=$1;
    local snpDB=$2;
    local chr=$3;
    local GATK=$4;
    local ref=$5;
    local dbsnp=$6;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx4g -jar $GATK GenotypeGVCFs -R $ref -D $dbsnp -V gendb://$snpDB -L $chr";
    command="$command --tmp-dir $outDir/tmp -O $outDir/$chr.vcf.gz ";
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$outDir/$chr.vcf.gz";
}

function VCF_FormFile {
    local inF=$1;
    local outDir=$2;
    local sample=$3;
    local GATK=$4;
    local ref=$5;
    local dbsnp=$6;
    local outLog=$7;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    if [ ! -d "$outDir/tmp" ]; then
        mkdir -p "$outDir/tmp" || error_exit "$outLog/tmp" "$sample\tMake temp dir fileed";
    fi; 
    local command="java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx4g -jar $GATK GenotypeGVCFs -R $ref -D $dbsnp -V $inF";
    command="$command --tmp-dir $outDir/tmp -O $outDir/$sample.gatk.vcf -stand-call-conf 20";
    printf "$command\n" >&2;
    $command || error_exit "$outLog/$sample" "$sample\t$command failed";
    write_log "$outLog/$sample" "$sample\t$command" $pStart;
    echo "$sample.gatk.vcf";
}

function MergeSamFiles_Chr() {
    local outDir=$1;
    local sample=$2; 
    local GATK=$3;
    local pStart=$(date +"%Y-%m-%d %H:%M:%S");
    local command="java -Xmx4g -jar $GATK MergeSamFiles";
    command="$command -I $outDir/chr1.mapped.recal.bam \\";
    command="$command -I $outDir/chr2.mapped.recal.bam \\";
    command="$command -I $outDir/chr3.mapped.recal.bam \\";
    command="$command -I $outDir/chr4.mapped.recal.bam \\";
    command="$command -I $outDir/chr5.mapped.recal.bam \\";
    command="$command -I $outDir/chr6.mapped.recal.bam \\";
    command="$command -I $outDir/chr7.mapped.recal.bam \\";
    command="$command -I $outDir/chr8.mapped.recal.bam \\";
    command="$command -I $outDir/chr9.mapped.recal.bam \\";
    command="$command -I $outDir/chr10.mapped.recal.bam \\";
    command="$command -I $outDir/chr11.mapped.recal.bam \\";
    command="$command -I $outDir/chr12.mapped.recal.bam \\";
    command="$command -I $outDir/chr13.mapped.recal.bam \\";
    command="$command -I $outDir/chr14.mapped.recal.bam \\";
    command="$command -I $outDir/chr15.mapped.recal.bam \\";
    command="$command -I $outDir/chr16.mapped.recal.bam \\";
    command="$command -I $outDir/chr17.mapped.recal.bam \\";
    command="$command -I $outDir/chr18.mapped.recal.bam \\";
    command="$command -I $outDir/chr19.mapped.recal.bam \\";
    command="$command -I $outDir/chr20.mapped.recal.bam \\";
    command="$command -I $outDir/chr21.mapped.recal.bam \\";
    command="$command -I $outDir/chr22.mapped.recal.bam \\";
    command="$command -I $outDir/chrX.mapped.recal.bam \\";
    command="$command -I $outDir/chrY.mapped.recal.bam \\";
    command="$command -O $outDir/$sample.mapped.recal.bam --CREATE_INDEX -SO"
    printf "$command\n" >&2;
    $command || error_exit "$outDir/$sample" "$sample\t$command failed";
    write_log "$outDir/$sample" "$sample\t$command" $pStart;
    echo "$outDir/$sample.mapped.recal.bam";
}

