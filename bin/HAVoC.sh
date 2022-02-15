#!/bin/bash

print_fancy_box() {
    label="$(printf $fastq | sed "s/.*\///g;s/_.*//g")"
    len="$(printf $label | wc -c)"
    box="$(expr $len + 16)"
    row="$(expr $box + 1)"
    space="$(expr $row - $box)"
    printf "┏"
    printf "━%.0s" $(eval "echo {1.."$(($box))"}")
    printf "┓\n"
    printf "┃ Processing $label..."
    printf " %.0s" $(eval "echo {1.."$(($space))"}")
    printf "┃\n"
    printf "┗"
    printf "━%.0s" $(eval "echo {1.."$(($box))"}")
    printf "┛\n"
}

print_instructions() {
    printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━INSTRUCTIONS━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    printf "\e[4mUsage:\n\e[0m"
    printf "    HAVoC.sh -n adapter -r reference -p tools_prepro -a tools_aligner -s tools_sam -m min_coverage -o run_pangolin read1 read2  \n"
    printf "\e[4m\nFollowing options can be changed in at commandline (only choose one):\n\e[0m"
    printf "    FASTQ preprocessing..........tools_prepro=\"fastp* or trimmomatic\"\n"
    printf "    Read aligner.................tools_aligner=\"bowtie or bwa*\"\n"
    printf "    SAM/BAM processing...........tools_sam=\"sambamba* or samtools\"\n"
    printf "    Masking threshold............min_coverage=30*\n"
    printf "    Run pangolin.................run_pangolin=\"yes* or no\"\n"
    printf "\n* By default.\n"
    printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n" && exit
}

####################################Options####################################
# Parsing options

OPTIND=1 # Reset OPTIND
while getopts :n:r:p:a:s:m:o:h opt
    do
        case $opt in
            n) adapter=$OPTARG;;
            r) reference=$OPTARG;;
            p) tools_prepro=$OPTARG;;
            a) tools_aligner=$OPTARG;;
            s) tools_sam=$OPTARG;;
            m) min_coverage=$OPTARG;;
            o) run_pangolin=$OPTARG;;
            h) print_instructions;;
        esac
    done

shift $(($OPTIND -1))

thread_num=8

#########################Unzip files and shorten names#########################
for old_name in $(ls *.fastq*); do
    if ! [ -f "$old_name" ]; then
        printf "No FASTQ files. Stopping assembly.\n" && exit
    elif [[ "$old_name" =~ \.gz$ ]]; then
        new_name="$(printf $old_name | sed "s/.*\///g" | sed "s/_.*R1.*\./_R1.fastq./" | sed "s/_.*R2.*\./_R2.fastq./")"
        cp $old_name $new_name
    else
        new_name="$(printf $old_name | sed "s/.*\///g" | sed "s/_.*R1.*\./_R1./" | sed "s/_.*R2.*\./_R2./")"
        cp $old_name $new_name
    fi
done

#############################Assembly and analyses#############################
fastq=$(ls *_R1.fastq*)
print_fancy_box
label="$(printf $fastq | sed "s/_R1.*//g" | sed "s/.*\///g")"

if ! [ -f "$label"_R2.fastq* ]; then
    printf "Missing R2 FASTQ file for $label. Skipping assembly.\n\n" && exit
fi

# Localize name of file
grep ">" $reference | sed "s/>.*/>$label/g" > "$label"_reference.fa
awk 'NR>1{printf "%s",$0}' $reference >> "$label"_reference.fa


if [ $tools_prepro = "fastp" ]; then
    printf "\e[4mRunning fastp...\n\e[0m"
    fastp --thread $thread_num -i "$label"_R1.fastq* -I "$label"_R2.fastq* -o "$label"_trimmed_1P -O "$label"_trimmed_2P --unpaired1 "$label"_trimmed_1U --unpaired2 "$label"_trimmed_2U -q 15 -u 40 -l 25 --cut_right --cut_window_size 20 --cut_mean_quality 30 --correction
    mv fastp.html fastp_"$label".html
    mv fastp.json fastp_"$label".json
else
    printf "\e[4mRunning Trimmomatic...\n\e[0m"
    trimmomatic PE -threads $thread_num -phred33 $directory/"$label"_R1.fastq* $directory/"$label"_R2.fastq* "$label"_trimmed_1P "$label"_trimmed_1U "$label"_trimmed_2P "$label"_trimmed_2U ILLUMINACLIP:$adapter:2:30:10 MINLEN:25 SLIDINGWINDOW:20:30 2> "$label"_trim_out.log
fi

if [ $tools_aligner = "bwa" ]; then
    printf "\e[4m\nRunnning BWA-MEM...\n\e[0m"
    bwa index "$label"_reference.fa
    bwa mem -t $thread_num -v 1 "$label"_reference.fa "$label"_trimmed_1P "$label"_trimmed_2P > "$label".sam
    rm "$label"_trimmed_*
else
    printf "\e[4m\nRunning Bowtie 2...\n\e[0m"
    bowtie2-build -q --threads $thread_num "$label"_reference.fa "$label"
    bowtie2 -p $thread_num -x "$label" -1 "$label"_trimmed_1P -2"$label"_trimmed_2P > "$label".sam
    rm "$label"*.bt2 "$label"_trimmed_*
fi

printf "\e[4m\nSorting, filling in mate coordinates, marking duplicate alignments, and indexing BAM file...\n\e[0m"
if [ $tools_sam = "sambamba" ]; then
    sambamba view --sam-input -o "$label".bam -f bam -t $thread_num "$label".sam
    sambamba sort -n -o "$label"_namesort.bam -t $thread_num "$label".bam
    rm "$label".sam

    samtools fixmate -m "$label"_namesort.bam "$label"_fixmate.bam -@ $thread_num

    sambamba sort -o "$label"_sorted.bam --tmpdir "$label"_tmpdir -t $thread_num "$label"_fixmate.bam

    sambamba markdup -r -t $thread_num --overflow-list-size=500000 "$label"_sorted.bam "$label"_markdup.bam

    sambamba index -t $thread_num "$label"_markdup.bam
else
    samtools sort -n -O bam -o "$label"_namesort.bam "$label".sam -@ $thread_num
    rm "$label".sam

    samtools fixmate -m "$label"_namesort.bam "$label"_fixmate.bam -@ $thread_num

    samtools sort -O bam -o "$label"_sorted.bam -T "$label"_temp.txt "$label"_fixmate.bam -@ $thread_num
    samtools markdup -r -s "$label"_sorted.bam "$label"_markdup.bam -@ $thread_num

    samtools index "$label"_markdup.bam -@ $thread_num
fi

printf "\e[4m\nMasking low coverage regions...\n\e[0m"
bedtools genomecov -ibam "$label"_markdup.bam -bga | awk -v cov="$min_coverage" '$4<cov' | bedtools merge -i - 1>"$label"_lowcovmask.bed
bedtools maskfasta -fi "$label"_reference.fa -bed "$label"_lowcovmask.bed -fo "$label"_reference_masked.fa
mv -f "$label"_reference_masked.fa "$label"_reference.fa
printf "Masked cov<$min_coverage regions in "$label"_reference.fa\n"

printf "\e[4m\nVariant calling with lofreq...\n\e[0m"
lofreq indelqual "$label"_markdup.bam --dindel -f "$label"_reference.fa -o "$label"_indel.bam
samtools index "$label"_indel.bam -@ $thread_num
lofreq call --call-indels -f "$label"_reference.fa -o "$label"_indel.vcf "$label"_indel.bam
lofreq filter -i "$label"_indel.vcf -o "$label"_indel_flt.vcf --af-min 0.5 --cov-min 10

bgzip -c "$label"_indel_flt.vcf 1> "$label"_indel_flt.vcf.gz
tabix "$label"_indel_flt.vcf.gz
bcftools consensus -f "$label"_reference.fa "$label"_indel_flt.vcf.gz -o "$label"_consensus.fa


if [ $run_pangolin = "yes" ]; then
    printf "\e[4m\nLineage mapping with Pangolin...\n\e[0m"
    pangolin --no-temp "$label"_consensus.fa --outfile "$label"_pangolin_lineage.csv -t $thread_num
    mv pangolearn_assignments.csv "$label"_pangolearn_assignments.csv
fi




printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
printf "All done! "
printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"

#################################Stop and exit#################################
exit
