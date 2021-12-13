#!/bin/bash -l

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
    printf "    sh HAVoC.sh [FASTQ directory]\n"
    printf "\nNexteraPE-PE.fa and ref.fa need to be in the same directory as HAVoC.\n"
    printf "\e[4m\nFollowing options can be changed in script (only choose one):\n\e[0m"
    printf "    Number of threads used.......thread_num=8*\n"
    printf "    FASTQ preprocessing..........tools_prepro=\"fastp* or trimmomatic\"\n"
    printf "    Read aligner.................tools_aligner=\"bowtie or bwa*\"\n"
    printf "    SAM/BAM processing...........tools_sam=\"sambamba* or samtools\"\n"
    printf "    Masking threshold............min_coverage=30*\n"
    printf "    Run pangolin.................run_pangolin=\"yes* or no\"\n"
    printf "    Run in CSC...................run_in_csc=\"yes or no*\"\n"
    printf "\n* By default.\n"
    printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n" && exit
}

###################################Set time###################################
total_time="$(date +%s)"
start_time="$(date +%s)"
####################################Options####################################
thread_num=8
tools_prepro="fastp"
tools_aligner="bwa"
tools_sam="sambamba"
min_coverage=30
run_pangolin="yes"
run_in_csc="no"
######################################Help######################################
if [ -z $1 ] || [ $1 = "-help" ] || [ $1 = "-h" ]; then
    print_instructions
fi
#############################Check necessary files#############################
path_to_files="$(printf $0)"
if ! [ -d "$1" ]; then
    printf "\e[1mError!\e[0m Given FASTQ directory doesn't exist.\n"
    print_instructions
fi
adapter="NexteraPE-PE.fa"
if ! [ -f $adapter ]; then
    printf "\e[1mError!\e[0m Missing adapter file. Stopping assembly.\n"
    print_instructions
fi
reference="ref.fa"
if ! [ -f $reference ]; then
    printf "\e[1mError!\e[0m Missing reference file. Stopping assembly.\n"
    print_instructions
fi
################################Choose directory################################
directory="$(printf $1)"
#########################Unzip files and shorten names#########################
for old_name in $directory/*.fastq*; do
    if ! [ -f "$old_name" ]; then
        printf "No FASTQ files. Stopping assembly.\n" && exit
    elif [[ "$old_name" =~ \.gz$ ]]; then
        new_name="$(printf $old_name | sed "s/.*\///g" | sed "s/_.*R1.*\./_R1.fastq./" | sed "s/_.*R2.*\./_R2.fastq./")"
        mv -f $old_name $directory/$new_name
    else
        new_name="$(printf $old_name | sed "s/.*\///g" | sed "s/_.*R1.*\./_R1./" | sed "s/_.*R2.*\./_R2./")"
        mv -f $old_name $directory/$new_name
    fi
done >/dev/null 2>&1
###############################Load modules/packages###########################
if [ $run_in_csc = "yes" ]; then
    printf "\e[4mLoading modules...\n\e[0m"
    module purge >/dev/null 2>&1
    module load biokit lofreq >/dev/null 2>&1
    if [ $run_pangolin = "yes" ]; then
        module load bioconda >/dev/null 2>&1
        conda activate pangolin2.2.2
        pangolin --update
    fi
    printf "Done!\n\n"
else
    printf "\e[4mLoading packages...\n\e[0m"
    #source ~/miniconda3/etc/profile.d/conda.sh
    #conda activate havoc
    pangolin --update
    printf "Done!\n\n"
fi
#############################Assembly and analyses#############################
(
    for fastq in *_R1.fastq*; do
        start_time2="$(date +%s)"

        label="$(printf $fastq | sed "s/_R1.*//g" | sed "s/.*\///g")"

        if ! [ -f $directory/"$label"_R2.fastq* ]; then
            printf "Missing R2 FASTQ file for $label. Skipping assembly.\n\n" && continue
        fi

        echo grep ">" $reference | sed "s/>.*/>$label/g" > "$label"_reference.fa
        awk 'NR>1{printf "%s",$0}' $reference >> "$label"_reference.fa

        if [ $tools_prepro = "fastp" ]; then
            printf "\e[4mRunning fastp...\n\e[0m"
            start_time3="$(date +%s)"
            fastp --thread $thread_num -i "$label"_R1.fastq* -I "$label"_R2.fastq* -o "$label"_trimmed_1P -O "$label"_trimmed_2P --unpaired1 "$label"_trimmed_1U --unpaired2 "$label"_trimmed_2U -q 15 -u 40 -l 25 --cut_right --cut_window_size 20 --cut_mean_quality 30 --correction
            #stop_time3="$(expr "$(date +%s)" - $start_time3)" && printf "fastp time : %.0f sec\n" $stop_time3
        else
            printf "\e[4mRunning Trimmomatic...\n\e[0m"
            start_time3="$(date +%s)"
            trimmomatic PE -threads $thread_num -phred33 $directory/"$label"_R1.fastq* $directory/"$label"_R2.fastq* "$label"_trimmed_1P "$label"_trimmed_1U "$label"_trimmed_2P "$label"_trimmed_2U ILLUMINACLIP:$adapter:2:30:10 MINLEN:25 SLIDINGWINDOW:20:30
            stop_time3="$(expr "$(date +%s)" - $start_time3)" && printf "Trimmomatic time : %.0f sec\n" $stop_time3
        fi

        if [ $tools_aligner = "bwa" ]; then
            printf "\e[4m\nRunnning BWA-MEM...\n\e[0m"
            start_time3="$(date +%s)"
            bwa index "$label"_reference.fa
            bwa mem -t $thread_num -v 1 "$label"_reference.fa "$label"_trimmed_1P "$label"_trimmed_2P > "$label".sam
            rm "$label"_trimmed_*
            # stop_time3="$(expr "$(date +%s)" - $start_time3)" && printf "BWA-MEM time : %.0f sec\n" $stop_time3
        else
            printf "\e[4m\nRunning Bowtie 2...\n\e[0m"
            start_time3="$(date +%s)"
            bowtie2-build -q --threads $thread_num "$label"_reference.fa "$label"
            bowtie2 -p $thread_num -x "$label" -1 "$label"_trimmed_1P -2"$label"_trimmed_2P > "$label".sam
            rm "$label"*.bt2 "$label"_trimmed_*
            stop_time3="$(expr "$(date +%s)" - $start_time3)"
            printf "Bowtie 2 time : %.0f min " "$(printf "$(expr $stop_time3 / 60)")"
            printf "%.0f sec\n\n" "$(printf "$(expr $stop_time3 % 60)")"
        fi

        printf "\e[4m\nSorting, filling in mate coordinates, marking duplicate alignments, and indexing BAM file...\n\e[0m"
        start_time3="$(date +%s)"
        if [ $tools_sam = "sambamba" ]; then
            start_time4="$(date +%s)"
            sambamba view --sam-input -o "$label".bam -f bam -t $thread_num "$label".sam
            sambamba sort -n -o "$label"_namesort.bam -t $thread_num "$label".bam
            rm "$label".sam
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Sambamba sort time : %.0f sec\n" $stop_time4

            start_time4="$(date +%s)"
            samtools fixmate -m "$label"_namesort.bam "$label"_fixmate.bam -@ $thread_num
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Samtools fixmate time : %.0f sec\n" $stop_time4

            start_time4="$(date +%s)"
            sambamba sort -o "$label"_sorted.bam --tmpdir "$label"_tmpdir -t $thread_num "$label"_fixmate.bam
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Sambamba sort time : %.0f sec\n" $stop_time4

            start_time4="$(date +%s)"
            sambamba markdup -r -t $thread_num --overflow-list-size=500000 "$label"_sorted.bam "$label"_markdup.bam
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Sambamba markdup time : %.0f sec\n" $stop_time4

            start_time4="$(date +%s)"
            sambamba index -t $thread_num "$label"_markdup.bam
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Samtools index time : %.0f sec\n\n" $stop_time4
            stop_time3="$(expr "$(date +%s)" - $start_time3)"
            printf "Sambamba/Samtools time : %.0f min " "$(printf "$(expr $stop_time3 / 60)")"
            printf "%.0f sec\n\n" "$(printf "$(expr $stop_time3 % 60)")"
        else
            start_time4="$(date +%s)"
            samtools sort -n -O bam -o "$label"_namesort.bam "$label".sam -@ $thread_num
            rm "$label".sam
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Samtools sort time : %.0f sec\n" $stop_time4

            start_time4="$(date +%s)"
            samtools fixmate -m "$label"_namesort.bam "$label"_fixmate.bam -@ $thread_num
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Samtools fixmate time : %.0f sec\n" $stop_time4

            start_time4="$(date +%s)"
            samtools sort -O bam -o "$label"_sorted.bam -T "$label"_temp.txt "$label"_fixmate.bam -@ $thread_num
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Samtools sort time : %.0f sec\n" $stop_time4

            start_time4="$(date +%s)"
            printf "\n"
            samtools markdup -r -s "$label"_sorted.bam "$label"_markdup.bam -@ $thread_num
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Samtools markdup time : %.0f sec\n" $stop_time4

            start_time4="$(date +%s)"
            samtools index "$label"_markdup.bam -@ $thread_num
            stop_time4="$(expr "$(date +%s)" - $start_time4)" && printf "   Samtools index time : %.0f sec\n" $stop_time4
            stop_time3="$(expr "$(date +%s)" - $start_time3)" && printf "Samtools time : %.0f sec\n" $stop_time3
        fi

        printf "\e[4m\nMasking low coverage regions...\n\e[0m"
        bedtools genomecov -ibam "$label"_markdup.bam -bga | awk -v cov="$min_coverage" '$4<cov' | bedtools merge -i - 1>"$label"_lowcovmask.bed
        bedtools maskfasta -fi "$label"_reference.fa -bed "$label"_lowcovmask.bed -fo "$label"_reference_masked.fa
        mv -f "$label"_reference_masked.fa "$label"_reference.fa
        printf "Masked cov<$min_coverage regions in "$label"_reference.fa\n"

        printf "\e[4m\nVariant calling with lofreq...\n\e[0m"
        start_time3="$(date +%s)"
        lofreq indelqual"$label"_markdup.bam --dindel -f "$label"_reference.fa -o "$label"_indel.bam
        samtools index "$label"_indel.bam -@ $thread_num
        if [ $run_in_csc = "yes" ]; then
            lofreq call-parallel --pp-threads $thread_num --call-indels -f "$label"_reference.fa -o "$label"_indel.vcf "$label"_indel.bam
        else
            lofreq call --call-indels -f "$label"_reference.fa -o "$label"_indel.vcf "$label"_indel.bam
        fi
        lofreq filter -i "$label"_indel.vcf -o "$label"_indel_flt.vcf --af-min 0.5 --cov-min 10

        bgzip -c "$label"_indel_flt.vcf 1> "$label"_indel_flt.vcf.gz
        tabix "$label"_indel_flt.vcf.gz
        bcftools consensus -f "$label"_reference.fa "$label"_indel_flt.vcf.gz -o "$label"_consensus.fa
        stop_time3="$(expr "$(date +%s)" - $start_time3)" && printf "Variant calling time : %.0f sec\n" $stop_time3


        if [ $run_pangolin = "yes" ]; then
            if [ $run_in_csc = "yes" ]; then
                printf "\e[4m\nLineage mapping with Pangolin...\n\e[0m"
                start_time3="$(date +%s)"
                pangolin "$label"_consensus.fa -o "$label" --outfile "$label"_pangolin_lineage.csv -t $thread_num #2>$directory/"$label"_pangolin.log
                stop_time3="$(expr "$(date +%s)" - $start_time3)"
                printf "Pangolin time : %.0f min " "$(printf "$(expr $stop_time3 / 60)")"
                printf "%.0f sec\n\n" "$(printf "$(expr $stop_time3 % 60)")"
            else
                printf "\e[4m\nLineage mapping with Pangolin...\n\e[0m"
                start_time3="$(date +%s)"
                pangolin "$label"_consensus.fa -o "$label" --outfile "$label"_pangolin_lineage.csv -t $thread_num #2>$directory/"$label"_pangolin.log
                stop_time3="$(expr "$(date +%s)" - $start_time3)"
                printf "Pangolin time : %.0f min " "$(printf "$(expr $stop_time3 / 60)")"
                printf "%.0f sec\n\n" "$(printf "$(expr $stop_time3 % 60)")"
            fi
        fi



        stop_time1="$(expr "$(date +%s)" - $start_time2)"
        printf "━━━━━━━━━━━━━━━━━━━━━━━━\n"
        printf "Run time : %.0f min " "$(printf "$(expr $stop_time1 / 60)")"
        printf "%.0f sec\n\n" "$(printf "$(expr $stop_time1 % 60)")"
    done

    stop_time2="$(expr "$(date +%s)" - $total_time)"
    printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    printf "All done! "
    printf "Total run time : %.0f hrs " "$(printf "$(expr $stop_time2 / 3600)")"
    printf "%.0f min " "$(printf "$(expr $stop_time2 / 60 % 60)")"
    printf "%.0f sec\n" "$(printf "$(expr $stop_time2 % 60)")"
    printf "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
) 2>&1 | tee "HAVoC_run.log"
#################################Stop and exit#################################
exit
