# Default values
outdir="$(pwd)"
barcode_range="01-96"
maxlen=2000
minq=10
mindepth=5
cpus=8
minlen=300  # Default minlen value
minalen=200
clair3_model_name="r1041_e82_400bps_hac_v430"

# Function to display usage information
usage() {
    echo "Usage: $0 [-i fastqdir] [-o outdir][-r ref_genome] [-p primers_fasta] [-b barcode_range] [-l minlen] [-x maxlen] [-a alen] [-q minq] [-d mindepth] [-c cpus]  [-m clair3_model_name]   [-h]"
    echo "  -i <inputdir>      input directory contains fastq files (default: current directory)"
    echo "  -o <outdir>         Working directory (default: current directory)"
    echo "  -r <ref_genome>      Reference genome file (default: ${ref_genome})"
    echo "  -p <primers_fasta>   Primers FASTA file (default: ${primers_fasta})"
    echo "  -b <barcode_range>   Barcode range (default: 01-96)"
    echo "  -l <minlen>          Min read length after chopping the adapter (default: 300)"
    echo "  -x <maxlen>          Maximum read length after chopping off the adapter (default: 2000)"
    echo "  -a <minalen>            Minimum length after trimming off primer  (default: 200)"               
    echo "  -q <minq>            Minimum quality score (default: 10)"
    echo "  -d <mindepth>        Minimum depth for variant calling (default: 5)"
    echo "  -c <cpus>            Number of CPUs to use (default: 8)"
    echo "  -m <clair3_model_name> clair model (default: the_script_dir/progs/clair3_models/r1041_e82_400bps_hac_v430)"
    echo "  -h                   Display this help message"
    echo
    echo "  Example:"
    echo "      $0 -i ./fastq -o results -b 73-85 -l 300 -x 2000 -a 200 -q 20 -d 5 -c 8 -m r1041_e82_400bps_hac_v400"
    echo
    exit 1
}

# Parse command-line options
while getopts "i:o:r:p:b:l:x:a:q:d:c:m:h" opt; do
    case ${opt} in
        i) inputdir="${OPTARG}" ;;
        o) outdir="${OPTARG}" ;;
        r) ref_genome="${OPTARG}" ;;
        p) primers_fasta="${OPTARG}" ;;
        b) barcode_range="${OPTARG}" ;;
        l) minlen="${OPTARG}" ;;  # This is the key line
        x) maxlen="${OPTARG}" ;;
        a) minalen="${OPTARG}" ;;
        q) minq="${OPTARG}" ;;
        d) mindepth="${OPTARG}" ;;
        c) cpus="${OPTARG}" ;;
        m) clair3_model_name="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check if minlen has been set, and assign default if not
if [ -z "$minlen" ]; then
    minlen=500  # default minlen value
    echo "Warning: minlen not specified. Using default value: $minlen"
fi
# Check if minlen has been set, and assign default if not
if [ -z "$minalen" ]; then
    minalen=300  # default minlen value
    echo "Warning: minlen not specified. Using default value: $minlen"
fi

echo "minlen is set to: $minlen"  # This line will print the value of minlen
echo "minalen is set to: $minalen"  # This line will print the value of minalen
# Validate barcode range format
if ! [[ $barcode_range =~ ^[0-9]+-[0-9]+$ ]]; then
    echo "Error: Invalid barcode range format. Expected format: start-end (e.g., 01-96)."
    usage
fi

# Ensure required options are provided
if [[ -z "$inputdir" ]]; then
    echo "Error: The input directory (-i) containing fastq files is required."
    usage
fi

# Get the directory of the current script
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Output the script directory
echo "Script directory: $SCRIPT_DIR"

progs="${SCRIPT_DIR}/progs"
refs="${SCRIPT_DIR}/refs"
clair3_model_path="${progs}/clair3_models/${clair3_model_name}"
ref_genome="${refs}/human_pev_refgenomes.fasta"
primers_fasta="${refs}/human_pev_primers.fasta"
primers_on_ref_bed="${refs}/locate_primers_on_ref.bed"

# Check if the clair3 model exist
if [[ ! -d "${clair3_model_path}" ]]; then
    echo "Clair3 model does not exist: ${clair3_model_path}, exit"
    exit
fi

# List of Singularity image files and their corresponding URLs
declare -A images
base_url="https://depot.galaxyproject.org/singularity"
images["seqkit"]="${base_url}/seqkit:2.9.0--h9ee0642_0"
images["nanoplot"]="${base_url}/nanoplot:1.43.0--pyhdfd78af_1"
images["porechop"]="${base_url}/porechop:0.2.4--py310h184ae93_9"
images["chopper"]="${base_url}/chopper:0.9.0--hdcf5f25_1"
images["hostile"]="${base_url}/hostile:2.0.0--pyhdfd78af_0"
images["csvtk"]="${base_url}/csvtk:0.31.0--h9ee0642_0"
images["minimap2"]="${base_url}/minimap2:2.28--he4a0461_3"
images["samtools"]="${base_url}/samtools:1.21--h96c455f_1"
images["ivar"]="${base_url}/ivar:1.4.3--h43eeafb_0"
images["clair3"]="${base_url}/clair3:1.0.8--py39h46983ab_2"
images["bcftools"]="${base_url}/bcftools:1.21--h3a4d415_1"
images["bedtools"]="${base_url}/bedtools:2.31.1--h13024bc_3"
images["vcfpy"]="${base_url}/vcfpy:0.13.8--pyhdfd78af_0"

declare -A tools

# Loop through the list of images and check whether images exist
for image in "${!images[@]}"; do
    url="${images[$image]}"
    file_name=$(basename "$url")
    image_path="${progs}/${file_name}"
    tools[$image]=${file_name}

    # Check if the image already exists
    if [ ! -f "$image_path" ]; then
        echo "Image '$image' does not exist. Downloading from $url..."
        wget -P ${progs} "$url"

        # Check if the download was successful
        if [ $? -eq 0 ]; then
        echo "Image '$image' downloaded successfully."
        else
        echo "Error downloading '$image'." >&2
        exit 1
        fi
    else
        echo "Image '$image' already exists in $outdir."
    fi
done

seqkit="singularity run ${progs}/${tools["seqkit"]} seqkit"
nanoplot="singularity run ${progs}/${tools["nanoplot"]} NanoPlot"
porechop="singularity run ${progs}/${tools["porechop"]} porechop"
chopper="singularity run ${progs}/${tools["chopper"]} chopper"
hostile="singularity run ${progs}/${tools["hostile"]} hostile"
csvtk="singularity run ${progs}/${tools["csvtk"]} csvtk"
minimap2="singularity run ${progs}/${tools["minimap2"]} minimap2"
samtools="singularity run ${progs}/${tools["samtools"]} samtools"
ivar="singularity run ${progs}/${tools["ivar"]} ivar"
bedtools="singularity run ${progs}/${tools["bedtools"]} bedtools"
clair3="singularity run ${progs}/${tools["clair3"]}  run_clair3.sh"
bcftools="singularity run ${progs}/${tools["bcftools"]} bcftools"
vcfpy="singularity run ${progs}/${tools["vcfpy"]} python"
reformatFasta="perl ${SCRIPT_DIR=}/scripts/reformatFasta.pl"

# Print working directory and results directory
echo "Current working directory: $outdir"
echo "Starting run at: $(date)"
echo "Using barcode range: $barcode_range"
echo "Using reference genome: $ref_genome"
echo "Using primers FASTA: $primers_fasta"

echo "Results directory: $outdir"
[ ! -d "$outdir" ] && mkdir -p "$outdir"

#produce primer matching reference bed file for trimming off primer from alignment later on
cmd="cat ${ref_genome} | ${seqkit} locate -d -i -f ${primers_fasta} --bed > ${primers_on_ref_bed}"
echo $cmd
echo ""
if [ ! -f "${primers_on_ref_bed}" ]; then
    eval "$cmd"
fi

# Iterate over the specified barcode range
for var in $(seq ${barcode_range//-/ }); do
    prefix="barcode$var"

    if [ $var -le 9 ]; then
        prefix="barcode0$var"
    fi

    
    # raw stats
    [ ! -d "${outdir}/${prefix}/raw" ] && mkdir -p ${outdir}/$prefix/raw
    cmd="${seqkit} stats -a -T -j 8 ${inputdir}/${prefix}.fastq.gz --out-file ${outdir}/${prefix}/raw/${prefix}.seqkit_stats.tsv"
    echo $cmd
    echo ""
    if [ ! -f "${outdir}/${prefix}/raw/${prefix}.seqkit_stats.tsv" ]; then
        eval "$cmd"
    fi
    
    ###################################### qc ##################
    #reads quality assesment

    # cmd="${nanoplot} --threads ${cpus} --fastq ${outdir}/fastq/${prefix}.fastq.gz --outdir ${outdir}/${prefix}/qc/nanoplot_raw -f json"
    # echo $cmd
    # echo ""
    # eval "$cmd"

    # chop adapters
    [ ! -d "${outdir}/${prefix}/qc/porechop" ] && mkdir -p ${outdir}/$prefix/qc/porechop

    cmd="${porechop} --threads ${cpus} --input ${inputdir}/${prefix}.fastq.gz \
    --output ${outdir}/${prefix}/qc/porechop/${prefix}.porechop.fastq.gz"

    echo $cmd
    echo ""
    if [ ! -f "${outdir}/${prefix}/qc/porechop/${prefix}.porechop.fastq.gz" ]; then
        eval "$cmd"
    fi

    cmd="${seqkit} stats -a -T -j 8 \
    ${outdir}/${prefix}/qc/porechop/${prefix}.porechop.fastq.gz \
    --out-file ${outdir}/${prefix}/qc/porechop/${prefix}.porechop.seqkit_stats.tsv"
    echo $cmd
    echo ""
    if [ ! -f "${outdir}/${prefix}/qc/porechop/${prefix}.porechop.seqkit_stats.tsv" ]; then
        eval "$cmd"
    fi

    # quality filtering
    [ ! -d "${outdir}/$prefix/qc/chopper" ] && mkdir -p ${outdir}/$prefix/qc/chopper
   
    cmd="zcat ${outdir}/${prefix}/qc/porechop/${prefix}.porechop.fastq.gz \
    | ${chopper} -q ${minq}  --minlength ${minlen}  --maxlength ${maxlen} \
    > ${outdir}/${prefix}/qc/chopper/${prefix}.chopper.fastq"
    echo $cmd
    echo ""
    if [ ! -f "${outdir}/${prefix}/qc/chopper/${prefix}.chopper.fastq" ]; then
        eval "$cmd"
    fi
    
    cmd="${seqkit} stats -a -T -j ${cpus} \
    --out-file ${outdir}/${prefix}/qc/chopper/${prefix}.chopper.seqkit_stats.tsv \
    ${outdir}/${prefix}/qc/chopper/${prefix}.chopper.fastq"
    echo $cmd
    echo ""
    if [ ! -f "${outdir}/${prefix}/qc/chopper/${prefix}.chopper.seqkit_stats.tsv" ]; then
        eval "$cmd"
    fi

    #dehost
    [ ! -d "${outdir}/$prefix/qc/dehost" ] && mkdir -p ${outdir}/$prefix/qc/dehost
    cmd="${hostile} clean --force --threads ${cpus} \
    --fastq1 ${outdir}/${prefix}/qc/chopper/${prefix}.chopper.fastq \
    --index human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401 \
    -o ${outdir}/${prefix}/qc/dehost"
    echo $cmd
    echo ""
    if [ ! -f "${outdir}/${prefix}/qc/dehost/${prefix}.chopper.clean.fastq.gz" ]; then
        eval "$cmd"
    fi

    #dehosted sequence stats
    cmd="${seqkit} stats -a -T -j ${cpus} \
    --out-file ${outdir}/${prefix}/qc/dehost/${prefix}.chopper.clean.seqkit_stats.tsv \
    ${outdir}/${prefix}/qc/dehost/${prefix}.chopper.clean.fastq.gz"
    echo $cmd
    echo ""
    if [ ! -f "${outdir}/${prefix}/qc/dehost/${prefix}.chopper.clean.seqkit_stats.tsv" ]; then
        eval "$cmd"
    fi

    cmd="${csvtk} concat --out-file ${outdir}/${prefix}/${prefix}_seqkit_stats.tsv \
    ${outdir}/${prefix}/raw/*.seqkit_stats.tsv ${outdir}/${prefix}/qc/porechop/*.seqkit_stats.tsv  \
    ${outdir}/${prefix}/qc/chopper/*.seqkit_stats.tsv ${outdir}/${prefix}/qc/dehost/*.seqkit_stats.tsv"
    echo $cmd
    echo ""
    if [ ! -f "${outdir}/${prefix}/${prefix}_seqkit_stats.tsv" ]; then
        eval "$cmd"
    fi

    # mapping
    [ ! -d "${outdir}/mapping" ] && mkdir -p ${outdir}/$prefix/mapping
    [ ! -d "${outdir}/$prefix/variants" ] && mkdir -p ${outdir}/$prefix/variants
    
    cmd="${minimap2} -ax map-ont -t ${cpus} ${ref_genome} \
    ${outdir}/${prefix}/qc/dehost/${prefix}.chopper.clean.fastq.gz \
    > ${outdir}/${prefix}/mapping/${prefix}.amplicon.sam;"
    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/${prefix}/mapping/${prefix}.amplicon.sam" ]; then
        eval "$cmd"
    fi

    cmd=" ${samtools} sort --reference ${ref_genome} \
    -o ${outdir}/${prefix}/mapping/${prefix}.amplicon.sorted.bam \
    ${outdir}/${prefix}/mapping/${prefix}.amplicon.sam; \
    ${samtools} index ${outdir}/${prefix}/mapping/${prefix}.amplicon.sorted.bam"
    
    
    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/${prefix}/mapping/${prefix}.amplicon.sorted.bam.bai" ]; then
        eval "$cmd"
    fi

    #soft clip the primers regions from the begining of the reads
    cmd="${ivar} trim -b ${primers_on_ref_bed} \
    -p ${outdir}/${prefix}/mapping/${prefix}.ivar_trim \
    -i ${outdir}/${prefix}/mapping/${prefix}.amplicon.sorted.bam \
    -q 1 -k -e -m ${minalen}"

    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/${prefix}/mapping/${prefix}.ivar_trim.bam" ]; then
        eval "$cmd"
    fi

    #filter out unmapped, 
    #not primary alignment, 
    #read fails platform/vendor quality checks, 
    #and read is PCR or optical duplicate supplementary alignment
    # 3844
    cmd="${samtools} view -T ${ref_genome} -h -F 3588 \
    ${outdir}/${prefix}/mapping/${prefix}.ivar_trim.bam \
    -o ${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered.bam"
    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered.bam" ]; then
        eval "$cmd"
    fi


    cmd="${samtools} sort \
    -o ${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered_sorted.bam \
    ${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered.bam; \
    ${samtools} index ${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered_sorted.bam;"
    
    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered_sorted.bam.bai" ]; then
        eval "$cmd"
    fi

    cmd=" ${samtools} coverage -d 0  \
    -o ${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered.coverage.txt \
    ${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered_sorted.bam"
    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered.coverage.txt" ]; then
        eval "$cmd"
    fi

    
    
    # mask reference
    cmd="${bedtools} genomecov -ibam ${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered_sorted.bam \
    -bga | awk '\$4<5' | ${bedtools} maskfasta -fi ${ref_genome} -bed - \
    -fo ${outdir}/${prefix}/mapping/masked_ref.fasta; \
    ${samtools} faidx ${outdir}/${prefix}/mapping/masked_ref.fasta"

    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/${prefix}/mapping/masked_ref.fasta" ]; then
        eval "$cmd"
    fi

    # assembly the whole genome 
    cmd="${clair3} \
    --bam_fn=${outdir}/${prefix}/mapping/${prefix}.amplicon.filtered_sorted.bam\
    --ref_fn=${outdir}/${prefix}/mapping/masked_ref.fasta \
    --model_path=${clair3_model_path} \
    --threads=${cpus} \
    --output=${outdir}/$prefix/variants \
    --platform="ont" \
    --include_all_ctgs \
    --haploid_precise  \
    --no_phasing_for_fa \
    --chunk_size=10000 \
    --var_pct_full=1 \
    --ref_pct_full=1 \
    --var_pct_phasing=1"

    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/$prefix/variants/merge_output.vcf.gz" ]; then
        eval "$cmd"
    fi
    

    
    #norm
    cmd="${bcftools} norm \
    --fasta-ref ${outdir}/${prefix}/mapping/masked_ref.fasta \
    --output ${outdir}/$prefix/variants/norm.vcf.gz \
    --check-ref w --output-type z --write-index=tbi \
    --threads 6 \
    ${outdir}/$prefix/variants/merge_output.vcf.gz"
    echo $cmd
    echo ""
    if [ ! -f "${outdir}/$prefix/variants/norm.vcf.gz" ]; then
        eval "$cmd"
    fi
    
    

    cmd="cat ${outdir}/${prefix}/mapping/masked_ref.fasta | ${bcftools} consensus \
        ${outdir}/$prefix/variants/norm.vcf.gz \
        --iupac-codes --mark-del '-' --mark-ins lc --mark-snv lc \
        > ${outdir}/$prefix/variants/${prefix}.consensus.fa"
    
    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/$prefix/variants/${prefix}.consensus.fa" ]; then
        eval "$cmd"
    fi

    cmd="${seqkit} \
    fx2tab \
    --only-id --name --length -C ATCG -C RYSWKMBDHV -C N -H \
    --threads 6 \
    ${outdir}/$prefix/variants/${prefix}.consensus.fa \
    -o ${outdir}/$prefix/variants/${prefix}.consensus_stats.txt"
    
    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/$prefix/variants/${prefix}.consensus_stats.txt" ]; then
        eval "$cmd"
    fi


    if [ ! -f "${outdir}/$prefix/variants/${prefix}.consensus.fa" ]; then
        eval "$cmd"
    fi

    cmd="${reformatFasta} --prefix ${prefix} \
    -i ${outdir}/$prefix/variants/${prefix}.consensus.fa \
    > ${outdir}/$prefix/variants/${prefix}.consensus_reformat.fa"

    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/$prefix/variants/${prefix}.consensus_reformat.fa" ]; then
        eval "$cmd"
    fi

    cmd="${seqkit} \
    fx2tab \
    --only-id --name --length -C ATCG -C RYSWKMBDHV -C N -H \
    --threads 6 \
    ${outdir}/$prefix/variants/${prefix}.consensus_reformat.fa \
    -o ${outdir}/$prefix/variants/${prefix}.consensus_reformat_stats.txt"

    echo "$cmd"
    echo ""
    if [ ! -f "${outdir}/$prefix/variants/${prefix}.consensus_reformat_stats.txt" ]; then
        eval "$cmd"
    fi

   
done

# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: $(date)"