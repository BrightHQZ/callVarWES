#!/bin/bash

source ./models2.sh; 

function call_SNP() {
    local inF1=$1;
    local outDir1=$2;
    local sample1=$3;
    local ref1=$4;
    local vcf11=$5;
    local GATK1=$6
    local intervals_exons1=$7;
    local inDrageStr=$8;
    local inLog=$9;

    #Prepared threads for call SNP.
    starttime=$(date +%s)
    export starttime
    printf "$outDir1\n" >&2;
    tmp_fifofile="$outDir1/$.fifo"
    printf "$tmp_fifofile\n" >&2;
    mkfifo $tmp_fifofile;
    exec 6<>$tmp_fifofile;
    rm $tmp_fifofile;

    for i in {1..6}
    do
        echo
    done >&6

    #call snp for chromosomes
    for i in {1..24}
    do
        read -u6
        {
            if [ $i -eq 24 ]; then
                printf $(callSNV_chr_germline_gvcf $inF1 $outDir1 $sample1 "chrY" "$intervals_exons1/chrY.bed" $GATK1 $ref1 $vcf11 $inDrageStr $inLog) >&2;
            elif [ $i -eq 23 ]; then
                printf $(callSNV_chr_germline_gvcf $inF1 $outDir1 $sample1 "chrX" "$intervals_exons1/chrX.bed" $GATK1 $ref1 $vcf1 $inDrageStr $inLog) >&2;
            else
                printf $(callSNV_chr_germline_gvcf $inF1 $outDir1 $sample1 "chr$i" "$intervals_exons1/chr$i.bed" $GATK1 $ref1 $vcf1 $inDrageStr $inLog) >&2;
            fi;
            echo >&6;   
        } &
    done
    wait
    exec 6>&-
    echo "OK";
}

#global var
outDir=sample=remote=proc=treeRoot=refV="";
#mapModel marked the beginning step: "BAM" names the input file is bam file which should be used in call snp and annotation.
mapModel="": 
#
fileListA=();
fileListB=();
samples=();
#tools
SOAPnuke="/data/bioTools/bin/SOAPnuke/SOAPnuke";
GATK="/data/bioTools/bin/gatk4/gatk-package-4.2.2.0-local.jar";
BWA="bwa";
SAMTOOLS="samtools";
sambamba="/data/bioTools/bin/sambamba";
bamdst="/data/bioTools/bin/bamdst/bamdst";
bcftools="bcftools"
annovar="/data/bioTools/bin/annovar"
refGene="refGene";
cytoBand="cytoBand";
gnomad211="gnomad211_genome"
clinvar="clinvar_20210501";
avsnp="avsnp150";
dbscsnv="dbscsnv11";
#reference
ref_hg38="/data/bioTools/resource/callSNP/hg38/index_BWA/Homo_sapiens_assembly38.fasta";
ref_hg19="/data/bioTools/resource/callSNP/hg19/index_BWA/ucsc.hg19.fasta";
#vcf
vcf38A="/data/bioTools/resource/callSNP/hg38/gatk/Homo_sapiens_assembly38.dbsnp138.vcf";
vcf38B="/data/bioTools/resource/callSNP/hg38/gatk/1000G_omni2.5.hg38.vcf";
vcf38C="/data/bioTools/resource/callSNP/hg38/gatk/hapmap_3.3.hg38.vcf";

vcf19A="/data/bioTools/resource/callSNP/hg19/gatk/gz/dbsnp_138.hg19.vcf";
vcf19B="/data/bioTools/resource/callSNP/hg19/gatk/gz/1000G_phase1.indels.hg19.vcf";
vcf19C="/data/bioTools/resource/callSNP/hg19/gatk/gz/Mills_and_1000G_gold_standard.indels.hg19.vcf";
#perl path
perlPath="/data/GWASCall/SH";
#blast path
blastn="/data/bioTools/bin/blast/bin/blastn";
blastdb="/data/bioTools/resource/blast/virus/virusRef";
#exons_region
exons19="/data/bioTools/resource/callSNP/hg19/gatk/regions/exons"
exons38="/data/bioTools/resource/callSNP/hg38/gatk/regions/exons"
#varDB
varDB="/data/varDB_wes";
intervals_genome="/data/bioTools/resource/index/grch38/GRCh38.bed";
intervals_hg19_genome="/data/bioTools/resource/index/hg19/hg19.bed";
#call snp in drageModel
DrageStr="/data/bioTools/resource/callSNP/hg19/gatk/regions/hg19_drags.str.zip"

while getopts ":f:o:s:r:v:t:1:2:S:V:" opt
do
    case $opt in
        f)
        fileType=$OPTARG;
        ;;
        t)
        treeRoot=$OPTARG; #data dir root
        ;;
        o)
        outDir=$OPTARG;
        ;;
        s)
        sample=$OPTARG
        ;;
        r)
        remote=$OPTARG
        ;;
        1)
        fastq1=$OPTARG
        ;;
        2)
        fastq2=$OPTARG
        ;;
        v)
        refV=$OPTARG
        ;;
        S)
        save=$OPTARG
        ;;
        V)
        probeV=$OPTARG
        ;;
        ?)
        echo "未知参数$OPTARG"
        exit 1;;
    esac
done

if ([ -z $fastq1 ] && [ -z $fastq2 ]) || [ -z $outDir ] || [ -z $sample ] || [ -z $treeRoot ] || [ -z $fileType ] || [ -z $refV ]; then
    echo "The base command line as following :";
    echo "bash callVariationFor_WES.sh -t xxxx -1 file1 -2 file2 -o dir -s sampleID -f fq -v hg19 -V A;"
    echo "bash callVariationFor_WES.sh -t xxxx -1 file1,file1.1,file1.2 -2 file2,file2.1,file2.2 -o dir -s sampleID,sampleID1,sampleID2 -f fq -v hg19 -V A -F Y;"
    echo "bash callVariationFor_WES.sh -t xxxx -1 file1 -o dir -s sampleID -f fq -v hg38;"
    echo "bash callVariationFor_WES.sh -t xxxx -1 file1.bam -o dir -s sampleID -f bam -v hg38;"
    echo "bash callVariationFor_WES.sh -t xxxx -1 file1.vcf -o dir -s sampleID -f vcf -v hg38;"
    echo "-r usename@1.1.1.1:/xxx/xxxx/ddddd.fq/bam (The file exists in the remote computer)";
    echo "-V is the type of catch chip of exon: A is Agilent, IDTV1 is IDT-V1"
    echo "Mul-input file names the cases will be analysised together as family. The input files should liks this: file1,file1.1,file1.2 and sampleID like this: sample1,sample2,sample3."
    echo "In the model of mul-input, the count of files and samples should equal."
    echo "Input values: -1 $fastq1  -2 $fastq2  -o $outDir  -s $sample  -t $treeRoot  -f $fileType  -v $refV  -V $probeV"
    exit;
fi;

echo "Input values: -1 $fastq1  -2 $fastq2  -o $outDir  -s $sample  -t $treeRoot  -f $fileType  -v $refV  -V $probeV"

if [ "$refV" != "hg19" ] && [ "$refV" != "hg38" ]; then
    echo "-v hg19 or -v hg38 were ok but others were fail.";
    exit;
fi;

if [ "$fileType" != "fq" ] && [ "$fileType" != "bam" ] && [ "$fileType" != "vcf" ]; then 
    echo "-f fq, -f bam, -f vcf were ok but others were fail.";
    exit;
fi;

if [ -z $fastq2 ]; then
    fastq2="|";
fi;

mapModel=$(checkFamily $fastq1 $fastq2 $sample $fileType)
if [ "$mapModel" != "SE" ] && [ "$mapModel" != "PE" ] && [ "$mapModel" != "BAM" ] && [ "$mapModel" != "VCF" ]; then
    echo $mapModel;
    exit;
elif [ "$mapModel" == "BAM" ] || [ "$mapModel" == "VCF" ] || [ "$mapModel" == "SE" ]; then
    fileListA=($(echo $fastq1 | tr "," "\n"));
    samples=($(echo $sample | tr "," "\n"));
    for index in "${!fileListA[@]}"; do
        fileTemp=${fileListA[index]};
        fileListA[index]=$treeRoot/${fileListA[index]};
        if [ ! -f ${fileListA[index]} ]; then
            echo "The $mapModel file of -1 ${fileListA[index]} is not exist!";
            exit;
        fi;
        mkdir -p "$outDir/bam/${samples[index]}" || error_exit $outDir "${samples[index]}\t$command fileed";
        mkdir -p "$outDir/gvcf/${samples[index]}" || error_exit $outDir "${samples[index]}\t$command fileed";
        mkdir -p "$outDir/vcf/${samples[index]}" || error_exit $outDir "${samples[index]}\t$command fileed";
    done;
elif [ "$mapModel" == "PE" ]; then
    fileListA=($(echo $fastq1 | tr "," "\n"));
    fileListB=($(echo $fastq2 | tr "," "\n"));
    samples=($(echo $sample | tr "," "\n"));
    #echo ${#samples[@]};
    for index in "${!fileListA[@]}"; do
       #echo "$treeRoot/${fileListA[$index]}";
       fileListA[index]=$treeRoot/${fileListA[index]};
       fileListB[index]=$treeRoot/${fileListB[index]};
       #echo ${fileListA[index]};
        if [ ! -f ${fileListA[index]} ] || [ ! -f ${fileListB[index]} ]; then
            echo "The files: ${fileListA[index]} or ${fileListB[index]} is not exists!";
            exit;
        fi;
        mkdir -p "$outDir/bam/${samples[index]}" || error_exit $outDir "${samples[index]}\t$command fileed";
        mkdir -p "$outDir/gvcf/${samples[index]}" || error_exit $outDir "${samples[index]}\t$command fileed";
        mkdir -p "$outDir/vcf/${samples[index]}" || error_exit $outDir "${samples[index]}\t$command fileed";
    done
fi;

#outDir="$treeRoot/$outDir";
#echo $outDir;

#Check out put dir.

if [ ! -d "$outDir/log" ]; then
    mkdir -p "$outDir/log" || error_exit $outDir "$sample\t$command fileed";
fi;

if [ ! -d "$outDir" ]; then
    mkdir -p "$outDir" || error_exit "$outDir/log/$sample" "$sample\t$command failed";
fi

if [ ! -d "$outDir/bam" ]; then
    mkdir -p "$outDir/bam" || error_exit "$outDir/log/$sample" "$sample\t$command fileed";
fi;

if [ ! -d "$outDir/gvcf" ]; then
    mkdir -p "$outDir/gvcf" || error_exit "$outDir/log/$sample" "$sample\t$command fileed";
fi;

if [ ! -d "$outDir/vcf" ]; then
    mkdir -p "$outDir/vcf" || error_exit "$outDir/log/$sample" "$sample\t$command fileed";
fi;

#Check file exists and aligment model.
if [ ! -z $remote ]; then
    if [ !-n $fastq1 ]; then
        fastq1=$(transFQ $fastq1 $outDir $sample"1" $remote);
    fi;
    if [ !-n $fastq2 ]; then
        fastq2=$(transFQ $fastq2 $outDir $sample"2" $remote);
    fi;
fi;

if [ -z "$probeV" ]; then
    probeV="common";
elif [ "$probeV" == "A" ]; then
    probeV="Agilent";
elif [ "$probeV" == "IDTV1" ]; then
    probeV="IDT-V1"
fi;

#Check refgenes version
caputerBed="";
if [ "$refV" == "hg19" ]; then
    ref=$ref_hg19;
    vcf1=$vcf19A;
    intervals_exons="$exons19/$probeV";
    caputerBed="$exons19/$probeV/$probeV.bed"
    intervals_genome="/data/bioTools/resource/index/hg19/hg19.bed";
else
    ref=$ref_hg38;
    refV="hg38";
    intervals_exons=$exons38;
    caputerBed="$exons38/$probeV/$probeV.bed"
    intervals_genome="";
fi;

echo "Mapping model is "$mapModel;

starttime=$(date +%s)
export starttime
mkdir -p "$outDir/bam/$starttime";
printf "$outDir/bam/$starttime\n" >&2;
tmp_fifofile="$outDir/bam/$starttime/$.fifo"
printf "$tmp_fifofile\n" >&2;
mkfifo $tmp_fifofile;
exec 3<>$tmp_fifofile;
rm -r "$outDir/bam/$starttime";

for i in {1..3}
do
    echo
done >&3

inF="";
for index in "${!fileListA[@]}"; 
do
    read -u3
    {
        sample=${samples[index]};
        if [ "$mapModel" == "SE" ]; then
            fastq1=${fileListA[index]};
            inF=$(qualityFASTQ_SE $SOAPnuke $fastq1 $outDir"/bam/"$sample $sample $outDir"/log");
            inF=$(aligmentFQ_SE "$outDir/bam/$sample/$sample.1.clean.fq.gz" $outDir"/bam/"$sample $sample $ref $BWA 30 $SAMTOOLS $outDir"/log");
            inF=$outDir"/bam/"$inF;
        elif [ "$mapModel" == "PE" ]; then
            fastq1=${fileListA[index]};
            fastq2=${fileListB[index]};
            inF=$(qualityFASTQ_PE $SOAPnuke $fastq1 $fastq2 $outDir"/bam/"$sample $sample $outDir"/log");
            #In this step the unmapped and mul-aligment reads were filtered by samtools. The detail shown the function of "aligmentFQ_P" in models2.sh 
            inF=$(aligmentFQ_PE "$outDir/bam/$sample/$sample.1.clean.fq.gz" "$outDir/bam/$sample/$sample.2.clean.fq.gz" $outDir"/bam/"$sample $sample $ref $BWA 30 $SAMTOOLS $outDir"/log");
            inF=$outDir"/bam/"$sample/$inF;
        elif [ "$mapModel" == "BAM" ]; then
            inF=${fileListA[index]}
        fi;
        if [ "$mapModel" != "VCF" ]; then 
            #inF=$outDir"/bam/sample202112.bam";
            inF=$(sortBam $inF $outDir"/bam/"$sample $sample $SAMTOOLS $outDir"/log")
            inF=$(markDup $outDir"/bam/"$sample/$inF $outDir"/bam/"$sample $sample $sambamba $outDir"/log");
            inF="$outDir/bam/$sample/$sample.sort.markdup.bam";
            $(qualityBAM $bamdst $inF "$outDir/bam/"$sample $caputerBed);
            if [ "$refV" == "hg19" ]; then
                inF=$(BaseRecalibrator $outDir/bam/$sample/$sample.sort.markdup.bam $outDir/bam/$sample $sample $GATK $ref $vcf19A $vcf19B $vcf19C $outDir/log/$sample)
            else
                inF=$(BaseRecalibrator $outDir/bam/$sample/$sample.sort.markdup.bam $outDir/bam/$sample $sample $GATK $ref $vcf38A $vcf38B $vcf38C $outDir/log/$sample)
            fi;
        fi;
        echo >&3;
    } &
done;
wait
exec 3>&-
echo "OK";

#Call snp using GATK;
if [ "$mapModel" != "VCF" ]; then
    #Call snp using bcftools
    #echo "Calling snp using bcftools...";
    #snpForbcftools=$(callSNV_callSNV_BCFTOOLS $inF $outDir"/vcf" $sample.$refV $ref $bcftools $outDir"/log");
    #echo "Call snp using bcftools finished!"
    #Collection the bam file for call snp;
    inF="";
    for bamFile in "${samples[@]}"; do
        inF="$outDir/bam/$bamFile/$bamFile.sort.markdup.recal.bam"
        sample=$bamFile;
        DrageFile=$(CalibrateDragstrModel $inF $outDir/bam/$bamFile $sample $GATK $DrageStr $outDir/log/$sample)
        #DrageFile="/data/snpForWes4/bam/MQX_YKD2821/MQX_YKD2821.drags";
        #Call snp using GATK;
        echo "Calling snp using GATK...$ref";
        callRes=$(call_SNP $inF $outDir/gvcf/$sample $sample.$refV $ref $vcf1 $GATK $intervals_exons $DrageFile $outDir"/log/"$sample);
        if [ "$callRes" == "OK" ]; then
            #process g.vcf.gz files
            $(GatherVcfs $outDir/gvcf/$sample $sample.$refV $GATK $outDir/log/$sample);
            inF="$sample.$refV.g.vcf.gz";
            $(IndexFeatureFile "$outDir/gvcf/$sample/$inF" $GATK $outDir/log/$sample);
        fi;
    done;

    if [ ${#samples[@]} == 1 ]; then
        inF="$outDir/gvcf/$sample/$inF";
        inF=$(gvcfTOvcf $inF $outDir/vcf/$sample $sample.$refV $GATK $ref $outDir/log);
        echo "Call snp using GATK finished!"
    else
        inF="";
        for gvcfFile in "${samples[@]}"; do
            inF="$inF$outDir/gvcf/$gvcfFile/$gvcfFile.$refV.g.vcf.gz,"
        done;
        sample=${samples[0]}
        inF=$(CombineGVCFs "$outDir/gvcf/$sample" $inF $sample.$refV.comb $GATK $ref $outDir"/log/"$sample);
        inF="$outDir/gvcf/$sample/$inF";
        inF=$(gvcfTOvcf $inF "$outDir/vcf/$sample" $sample.$refV.comb $GATK $ref $outDir"/log/"$sample);
        echo "Call snp using GATK finished!"
    fi;
    
    inF="$outDir/vcf/$sample/$inF";
    annoGATK=$(anntationSNP $inF "$outDir/vcf/$sample" $sample $refV $annovar $outDir"/log" $refGene $cytoBand $gnomad211 $clinvar $avsnp $dbscsnv);
    #echo "All analysis were finished and the annotated file in the directer of $outDir/vcf (names: $annoBcftools, $annoGATK)";
    echo "All analysis were finished and the annotated file in the directer of $outDir/vcf (names: $annoGATK)";

    for sampleT in "${samples[@]}"; do
        echo $(clearTempFile $outDir $sampleT $outDir"/log/"$sampleT $refV);
    done;
else
    for index in "${!fileListA[@]}"; do
        fastq1=${fileListA[index]};
        sample=${samples[index]};
        annoBcftools=$(anntationSNP $fastq1 "$outDir/vcf" $sample.$refV $refV $annovar $outDir"/log" $refGene $cytoBand $gnomad211 $clinvar $avsnp $dbscsnv);
        echo "All analysis were finished and the annotated file in the directer of $outDir/vcf (names: $annoBcftools)";
    done;
fi;






    #Prepared threads for call SNP.
    # starttime=$(date +%s)
    # export starttime
    # echo $outDir;
    # tmp_fifofile="$outDir/$.fifo"
    # echo $tmp_fifofile
    # mkfifo $tmp_fifofile
    # exec 6<>$tmp_fifofile
    # rm $tmp_fifofile

    # for i in {1..6}
    # do
    #     echo
    # done >&6

    # #call snp for chromosomes
    # for i in {1..24}
    # do
    #     read -u6
    #     {
    #         if [ $i -eq 24 ]; then
    #             callSNV_chr_germline_gvcf $inF "$outDir/gvcf" $sample "chrY" "$intervals_exons/chrY.bed" $GATK $ref $vcf1 $outDir"/log"
    #         elif [ $i -eq 23 ]; then
    #             callSNV_chr_germline_gvcf $inF "$outDir/gvcf" $sample "chrX" "$intervals_exons/chrX.bed" $GATK $ref $vcf1 $outDir"/log"
    #         else
    #             callSNV_chr_germline_gvcf $inF "$outDir/gvcf" $sample "chr$i" "$intervals_exons/chr$i.bed" $GATK $ref $vcf1 $outDir"/log"
    #         fi;
    #         echo >&6;   
    #     } &
    # done
    # wait
    # exec 6>&-

    #process g.vcf.gz files
    #inF=$(CombineGVCFs_chr "$outDir/gvcf" $sample $bcftools $outDir"/log")
    #inF="$outDir/gvcf/$inF";
    #inF=$(VCF_FormFile $inF "$outDir/vcf" $sample $GATK $ref $vcf1 $outDir"/log")
#fi;

#echo $inF1;

#if [ ! -n $outDir ] || [ -z $outDir ] || [ ! -n $sample ] || [ -z $sample ]; then
#    echo "Example: bash NIPT.sh -i file -o dir -s sampleID -v Y -p bam -r usename@1.1.1.1: -j bam"
#    echo "-i: input file name and abosolute path"
#    echo "-o: output dir"
#    echo "-s: sample id"
#    echo "-v: either \"Y\" or \"N\" (using mapped reads to math virus database)";
#    echo "-p: either \"bam\" (to finish pre-process for call variation) or \"vcf\"(pre-process and call variation)";
#    echo "    the default value of -p is bam"
#    echo "-r: the use name and ip address of remote computer";
#    echo "-j: name is 'jump', the value of it is only 'bam' but the default value is null"
#    echo "Note:"
#    echo "1: the value of -i, -o and -s were necessary but others ignore."
#    echo "2: if you want to get vcf file then -p is necessary."
#    echo "3: if the value of -r is not null then the value of -i should the file's abosolut path in remote computer!"
#    echo "4: the program doesn't support password type of remote computer. So, please sure password is not necessary when "
#    echo "   transform file form remote computer to your compute is not necessary. Testing using by 'scp' command."
#    exit;
#fi;



#if [ $(checkStr $remote) == "1" ]; then
#    if [ !-n $fastq1 ]; then
#        fastq1=$(transFQ $fastq1 $outDir $sample"1" $remote);
#    fi;
#    if [ !-n $fastq2 ]; then
#        fastq2=$(transFQ $fastq2 $outDir $sample"2" $remote);
#    fi;
#fi;

#if [ !-n $fastq1 ]; then
#    fastq1=$(qualityFASTQ $SOAPnuke $fastq1 $outDir $sample"1");
#    echo $fastq1;
#fi;
#if [ !-n $fastq2 ]; then
#    fastq2=$(qualityFASTQ $SOAPnuke $fastq2 $outDir $sample"2");
#    echo $fastq2;
#fi;



# inF=$(aligmentFQ $fastq1 $fastq2 $outDir $sample $ref $BWA 30 $SAMTOOLS);
# echo $inF;
# inF=$(markDup $inF $outDir $sample $sambamba);
# echo $inF;

#process bam for mapped reads in human genome
#inF=$(mappedReads $inF $outDir $sample $SAMTOOLS);
#echo $inF;
#bamQC report
#pStart=$(date +"%Y-%m-%d %H:%M:%S");
#command="perl $perlPath/bamStat.pl -i $outDir/$sample.sort.markdup.bam -o $outDir";
#printf "$command\n"
#$command;
#perl "$perlPath/bamStat.pl" -i $outDir/$sample.sort.markdup.bam -o $outDir
#write_log "$outDir/$sample" "$sample\t$command" $pStart;
#pStart=$(date +"%Y-%m-%d %H:%M:%S");
#command="perl $perlPath/plot.in.pl $outDir $outDir"
#printf "$command\n"
#$command;
#write_log "$outDir/$sample" "$sample\t$command" $pStart;
#clear temp file
#command="rm -fr $outDir/tmp $outDir/$sample.table $outDir/$sample.bam $outDir/$sample.sort* $outDir/*.txt $outDir/$sample.fq.gz $outDir/$sample.mapped.bam $outDir/$sample.mapped.bam.bai"
#printf "$command\n"
#$command;


# starttime=$(date +%s)
# export starttime
# echo $outDir;
# tmp_fifofile="$outDir/$.fifo"
# echo $tmp_fifofile
# mkfifo $tmp_fifofile
# exec 6<>$tmp_fifofile
# rm $tmp_fifofile

# for i in {1..6}
# do
#     echo
# done >&6

# #call snp for chromosomes
# for i in {1..24}
# do
#     read -u6
#     {
#         callSNV_chr_germline_gvcf $inF $outDir $sample $i $GATK $ref $vcf1
#         echo >&6;   
#     } &
# done
# wait
# exec 6>&-
# #process g.vcf.gz files
# inF=$(CombineGVCFs_chr $outDir $sample $bcftools)
# rm $outDir/chr*.g.vcf*;

# if [ "$save" == "db" ]; then
#     gvcfTOdb $inF $outDir $sample $GATK $varDB $intervals_genome
# elif [ "$save" == "vcf" ]; then
#     VCF_FormFile $inF $outDir $sample $GATK $ref $dbsnp;
# elif [ "$save" == "both" ]; then
#     gvcfTOdb $inF $outDir $sample $GATK $varDB $intervals_genome
#     VCF_FormFile $inF $outDir $sample $GATK $ref $dbsnp;
# fi;

    #inF=$(VariantRecalibrator $inF $outDir $sample $ref $GATK $vcf1 $vcf2 $vcf3 $vcf4 "SNP");
    #inF=$(VariantRecalibrator $inF $outDir $sample $ref $GATK $vcf1 $vcf2 $vcf3 $vcf4 "INDEL"); 
