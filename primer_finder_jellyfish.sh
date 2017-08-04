#!/bin/bash

:<<'END'

Wiki
The 13 serotypes of L. monocytogenes can cause disease, 
but more than 90% of human isolates belong to only three serotypes: 1/2a, 1/2b, and 4b.

Borucki and Call 2003

Target ForwardPrimer ReversePrimer
Division_I_III CGATATTTTATCTACTTTGTCA TTGCTCCAAAGCAGGGCAT 
Division_II GCGGAGAAAGCTATCGCA TTGTTCAAACATAGGGCTA
1/2a_3a TTACTAGATCAAACTGCTCC AAGAAAAGCCCCTCGTCC 
1/2b_3b AAAGTGAGTTCTTACGAGATTT AATTAGGAAATCGACCTTCT 

Chen and Knabel 2007

ECI AATAGAAATAAGCGGAAGTGT TTATTTCCTGTCGGCTTAG         -> ECI=Epidemic Clone I -> small number of closely-related 4b
ECII ATTATGCCAAGTGGTTACGGA ATCTGTTTGCGAGACCGTGTC    -> 
ECIII TTGCTAATTCTGATGCGTTGG GCGCTAGGGAATAGTAAAGG
4b AGTGGACAATTGATTGGTGAA CATCCATCCCTTACTTTGGAC        -> same primers as 4b_4d_4e below
1/2a GAGTAATTATGGCGCAACATC CCAATCGCGTGAATATCGG
Listeria_monocytogenes TGTCCAGTTCCATTTTTAACT TTGTTGTTCTGCTGTACGA
Listeria_spp ATGAATATGAAAAAAGCAAC TTATACGCGACCGAAGCCAAC

Doumith et al 2004

1/2a_1/2c_3a_3c AGGGCTTCAAGGACTTACCC ACGATTTCTGCTTGCCATTC
1/2c_3c AGGGGTCTTAAATCCTGGAA CGGCTTGTTCGGCATACTTA
1/2b_3b_4b_4d_4e AGCAAAATGCCAAAACTCGT CATCACTAAAGCCTCCCATTG
4b_4d_4e AGTGGACAATTGATTGGTGAA CATCCATCCCTTACTTTGGAC
Listeria_spp GCTGAAGAGATTGCGAAAGAAG CAAAGAAACCTTGGATTTGCGG

Tao et al 2015

Listeria_monocytogenes TGATGAAATAAAGGTCCACG CAAGCCATAATGAACAAACG
END


#file containing the primers to verify
#3 tab-separated columns:
#Serotype        Fwd_primer_sequence        Rev_primer_sequence




############
#          #
#   Help   #
#          #
############


function print_help() {
    echo "\
Usage: bash primer_finder_jellyfish.sh [-h] -f|-q -p primers.txt -o report.txt sample_file_1 [sample_file_n]

Mandatory flags:
    
    -f      Input file is fasta format. Choose \"-f\" or \"-q\".
            Multiple files are accepted and must be separated by a \"space\" character.
            Multiple files will be concatenated.

    -q      Input file is fastq format. Choose \"-f\" or \"-q\".
            Multiple files are accepted and must be separated by a \"space\" character.
            Multiple files will be concatenated.

    -p      Input primer sequence file to test.
            Format: 3 tab-separated columns.
            Colunm 1: sequence name.
            Colunm 2: forward primer sequence.
            Column 3: reverse primer sequece.

    -o      Output folder.

Optional flag:
    
    -h      Print this help.\
"
exit 1
}


#http://wiki.bash-hackers.org/howto/getopts_tutorial
# Default values
fasta=0
fastq=0
primers=""
output=""

options=':fqp:o:h'

while getopts "$options" opt; do
    case "$opt" in
        f) fasta=1;;
        q) fastq=1;;
        p) primers="$OPTARG";;
        o) output="$OPTARG";;
        h) print_help;;
        \?) echo "Invalid option: -"$OPTARG"" >&2; print_help;;
        :) echo "Option -"$OPTARG" requires an argument." >&2; print_help;;
    esac
done

# Get frags and arguments
shift $(($OPTIND - 1))


# Exit if argument is missing
if [ -z "$fasta" ] && [ -z "$fastq" ] || [ -z "$primers" ] || [ -z "$output" ]; then 
    echo "All \"-f or -q\", \"-p\" and \"-o\" options are mandatory"
    print_help
fi

# Check if only one of "-f" or "-q" is being used
if [ "$fasta" -eq 1 ] && [ "$fastq" -eq 1 ]; then
    echo "Only one type of input file can be selected: fasta (\"-f\") or fastq (\"-q\")."
    print_help
fi

# Check if primer file exists and in valid format
if [ ! -s "$primers" ]; then
    echo "The primer file provided does not exist."
    print_help
fi

# Check if primer file in valid format (3 columns)
while read line; do
    nFields=$(awk '{print NF}' <<< "$line")
    # echo "$nFields"
    if [ "$nFields" -ne 3 ]; then
        echo "Invalid primer input file format."
        print_help
    fi
done < "$primers"

# Check if output folder already exists
if [ -d "$output" ]; then
    echo "Output folder already exists. Output files will be overwritten."
else
    mkdir -p "$output"
fi

# Array to store the input sequence types
seq_types=()

function check_format()
{
    if [ "${1##*.}" == "gz" ]; then #is it gzipped?
        line1="$(zcat "$1" | head -n 1)" #assuming the first line is a sequence header
    else
        line1="$(cat "$1" | head -n 1)"
    fi

    char1="${line1:0:1}"

    if [ $(grep -E "fastq|fq" <<< "$1") ]; then #is input fastq file? Check file extension fisrt.
        if [ "$char1" == "@" ]; then
            seq_types+=("fastq")
        else
            echo "Invalid fastq file: "$1""
            exit 1
        fi
    else #it's a fasta
        if [ "$char1" = ">" ]; then
            seq_types+=("fasta")
        else
            echo "Invalid fasta file: "$1""
            exit 1
        fi
    fi
}

#check how many input files
nArgs=$#
input_seq_files=()
for f in "${@}"; do
    input_seq_files+=("$f")
done
# echo "${input_seq_files[@]}"

# Check format for each input sequence file
for seq in "${input_seq_files[@]}"; do
    check_format "$seq"
done

# check if all sequences the same type
all_types=($(echo "${seq_types[@]}" | tr " " "\n" | sort | uniq | tr "\n" " "))
# echo "${#all_types[@]}"

if [ "${#all_types[@]}" -gt 1 ]; then
    echo "Both fasta and fastq files were detected as input."
    print_help
fi


###########
#         #
#   I/O   #
#         #
###########


#fastq files

name="$(basename "$1")" #name without path
sampleName="$(cut -d "_" -f 1 <<< "${name%%.*}")" #nameR1 with the ".fastq.gz"


#################
#               #
#   Resources   #
#               #
#################


cpu=$(nproc)  # All CPUs
mem="$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100))" #85% of total available memory


#############################
#                           #
#   Primer Sequence Files   #
#                           #
#############################


#removes length-specific primer files. Because appending.
[[ -e "${output}"/"${sampleName}"_primers_??-mer.txt ]] && rm "${output}"/"${sampleName}"_primers_??-mer.txt

declare -a sizes=()
#make list of primer sizes
while IFS= read -r line || [[ -n "$line" ]]; do
    line=$(sed -e 's/ /\t/' -e 's/\r//g' -e 's/\n//g' <<< "$line") #transform space to tabs, remove carriage return
    #echo "$line"
    [ -n "$line" ] || continue #skip if line is empty
    
    #Primer sequences
    f_primer=$(cut -f 2 <<< "$line")
    r_primer=$(cut -f 3 <<< "$line")

    #Primer lengths
    f_primerLen=$(expr length "$f_primer")
    r_primerLen=$(expr length "$r_primer")

    sizes+=("$f_primerLen" "$r_primerLen")

    #copy line in a primer seqence in file for each length
    echo -e "$f_primer" >> "${output}"/"${sampleName}"_primers_"${f_primerLen}"-mer.txt
    echo -e "$r_primer" >> "${output}"/"${sampleName}"_primers_"${r_primerLen}"-mer.txt
done < "$primers"

#Remove duplicates and sort array
eval sizes=($(printf "%q\n" "${sizes[@]}" | sort -u))
#echo "${sizes[@]}"


#################
#               #
#   Jellyfish   #
#               #
#################


#Make kmer count database with Jellyfish
# cat $(echo "${input_seq_files[@]}")

#mem2000000
if [ "$fastq" -eq 1 ]; then
    for s in "${sizes[@]}"; do
        # Count k-mer
        zcat $(echo "${input_seq_files[@]}") \
        | jellyfish count \
            -s "$mem" \
            -m "$s" \
            -C \
            -t "$cpu" \
            -o "${output}"/"${sampleName}"_"${s}"-mer.counts \
            -L 3 \
            /dev/stdin

        # # Using Bloom filter with the 2 pass method
        # zcat < "$r1" < "$r2" \
        # | jellyfish bc \
        #     -s "$mem" \
        #     -m "$s" \
        #     -C \
        #     -t "$cpu" \
        #     -o "${output}"/"${sampleName}"_"${s}"-mer.bc \
        #     /dev/stdin
        # # count
        # zcat < "$r1" < "$r2" \
        # | jellyfish count \
        #     -s "$mem" \
        #     -m "$s" \
        #     -C \
        #     -t "$cpu" \
        #     --bc "${output}"/"${sampleName}"_"${s}"-mer.bc \
        #     -o "${output}"/"${sampleName}"_"${s}"-mer.counts \
        #     /dev/stdin
    done

elif [ "$fasta" -eq 1 ]; then
    for s in "${sizes[@]}"; do
        cat $(echo "${input_seq_files[@]}") \
        | jellyfish count \
            -s "$mem" \
            -m "$s" \
            -C \
            -t "$cpu" \
            -o "${output}"/"${sampleName}"_"${s}"-mer.counts \
            /dev/stdin
        done
else  # should not get here
    echo "Please use a fastq or a fasta file as input. Something went wrong."
    exit 1
fi


######################
#                    #
#   Coverage stats   #
#                    #
######################


: <<'END'
#make histogram of kmer abundance distribution
jellyfish histo \
    "${output}"/"${sampleName}"_"${size}"-mer.counts* \
    > "${output}"/"${sampleName}"_"${size}"-mer.histo

#find peak in excel - plot form 3 to 100
#have to find the peak automatically and produce a figure to support it

#estimage coverage and genome size
estimate_genome_size.pl \
    --kmer="$size" \
    --peak=45 \
    --fastq="$r1" "$r2"
END


#######################
#                     #
#   Check for k-mer   #
#                     #
#######################


#check for presence of kmer (count)
for s in "${sizes[@]}";do
    jellyfish query \
        "${output}"/"${sampleName}"_"${s}"-mer.counts \
        $(cat "${output}"/"${sampleName}"_primers_"${s}"-mer.txt | tr "\n" " ") \
        | tee "${output}"/insilicoPCR_"${s}"-mer_output.txt
done


##############
#            #
#   Report   #
#            #
##############


function reverse_complement ()
{
    #make sure all character are to uppercase
    upper=$(awk '{print toupper($0)}' <<< "$1")  # $1 is a string of ATCG
    rc=$(echo "$upper" | rev | tr "ATGC" "TACG")
    echo "$rc"
}


[ -e "${output}"/report.tsv ] && rm "${output}"/report.tsv
while IFS= read -r line || [[ -n "$line" ]]; do
    line=$(sed -e 's/ /\t/' -e 's/\r//g' -e 's/\n//g' <<< "$line") #transform space to tabs, remove carriage return
    
    #Primer sequences
    name=$(cut -f 1 <<< "$line")

    f_name=""${name}"_F"
    r_name=""${name}"_R"

    f_primer=$(cut -f 2 <<< "$line")
    r_primer=$(cut -f 3 <<< "$line")

    #get the counts from the jellyfish query output files
    f_counts=$(cat "${output}"/insilicoPCR_*-mer_output.txt | grep -wF "$f_primer" | cut -d " " -f 2)
    r_counts=$(cat "${output}"/insilicoPCR_*-mer_output.txt | grep -wF "$r_primer" | cut -d " " -f 2)

    # try the reverse complment of the primer
    if [ ! "$f_counts" ]; then
        r_c=$(reverse_complement "$f_primer")
        f_counts=$(cat "${output}"/insilicoPCR_*-mer_output.txt | grep -wF "$r_c" | cut -d " " -f 2)
        [ "$f_counts" ] || f_counts=0  # if still not present: error!?
    fi

    if [ ! "$r_counts" ]; then
        r_c=$(reverse_complement "$r_primer")
        r_counts=$(cat "${output}"/insilicoPCR_*-mer_output.txt | grep -wF "$r_c" | cut -d " " -f 2)
        [ "$r_counts" ] || r_counts=0  # if still not present: error!?
    fi

    printf "%s\t%s\t%s\n" "$f_name" "$f_primer" "$f_counts" | tee -a "${output}"/report.tsv
    printf "%s\t%s\t%s\n" "$r_name" "$r_primer" "$r_counts" | tee -a "${output}"/report.tsv
done < "$primers"
