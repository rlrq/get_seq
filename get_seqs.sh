#!/bin/bash

## TODO: accommodate GFF3 instead of requiring the GFF2BED version
ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

## if no arguments, show manual
if [[ $# -eq 0 ]]; then
    man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1
    exit 1
fi


ACC_DEFAULT='ref'
BED_DEFAULT='/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed'
RPS_DB_DEFAULT='/mnt/chaelab/shared/blastdb/Cdd.v3.18/Cdd'
REFERENCE_DEFAULT='/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta'
DIR_DEFAULT="/mnt/chaelab/$(whoami)/get_seq"
FEATURE_DEFAULT='gene'
ACC_LOOKUP_DEFAULT="${SCRIPT_DIR}/data/1135acc.csv"
ATTR_MOD_DEFAULT='{}'
DOMAIN_DEFAULT=''
DOMAIN_F_DEFAULT=''
GENETIC_CODE_DEFAULT=1
MODE_DEFAULT='col0'

while (( "$#" )); do
    case "$1" in
        -g|--gene) GENE="${2}";;
        -a|--acc) ACC="${2}";;
        -i|--input) ACCS_F="$(realpath ${2})";;
        -f|--feature) FEATURE="${2}";;
        -d|--dir) DIR="${2}";;
        -c|--chrom) CHROM="${2}";;
        -s|--start) START="${2}";; ## 1-indexed
        -e|--end) END="${2}";; ## 1-indexed
        -b|--bed) BED="${2}";;
        -r|--reference) REFERENCE="$(realpath ${2})"; MODE='ref';;
        # -r|--reference) REFERENCE="${2}"; MODE='ref';; ## now acc-lookup and ACC_LOOKUP
        -o|--out) FOUT="$(realpath ${2})";;
        --append) APPEND="${2}";;
        --gff) GFF="$(realpath ${2})";;
        --genetic-code) GENETIC_CODE="${2}";;
        # --fasta) FASTA="${2}";;
        --domain) DOMAIN="${2}";;
        --domain-file) DOMAIN_F="$(realpath ${2})";;
        --domain-db) RPS_DB="$(realpath ${2})";;
        --input-dir) INPUT_DIR="$(realpath ${2})";;
        --bed-out) FOUT_BED="$(realpath ${2})";;
        --qname-dname) QDNAME="${2}";;
        --qstart-qend) QSTARTEND="${2}";;
        --acc-lookup) ACC_LOOKUP="$(realpath ${2})"; MODE='ref';;
        --attr-mod) ATTR_MOD="${2}";; ## e.g. rice's GFF uses Locus_id instead of Parent to indicate gene for mRNA entries, so --attr-mod "{'mRNA': {'Parent': 'Locus_id'}}" would help to map it properly
        --no-bed) NO_BED='True';;
        --merge) MERGE='True';;
        --translate) TRANSLATE='True';;
        --no-header) HEADER='False';;
        --complete) COMPLETE='True';;
        --adjust-dir) ADJ_DIR='True';;
        --by-gene) BY_GENE='True';;
        --filter-fasta)FILTER_FASTA='True';; ## if true, fasta file will be filtered for relevant molecules. Useful if fasta file has many small molecules but only a handful of them contain the desired gene(s)
        -h|--help) man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1; exit 0;;
        --readme) cat ${SCRIPT_DIR}/README; exit 0;;
    esac
    shift
done

ACC="${ACC:-${ACC_DEFAULT}}"
DIR="${DIR:-${DIR_DEFAULT}}"
if [ -z "${GFF}" ]; then
    BED="${BED:-${BED_DEFAULT}}"
fi
FEATURE="${FEATURE:-${FEATURE_DEFAULT}}"
ACC_LOOKUP="${ACC_LOOKUP:-${ACC_LOOKUP_DEFAULT}}"
# REFERENCE="${REFERENCE:-FASTA}"
REFERENCE="${REFERENCE:-${REFERENCE_DEFAULT}}"
RPS_DB="${RPS_DB:-${RPS_DB_DEFAULT}}"
HEADER="${HEADER:-True}"
COMPLETE="${COMPLETE:-False}"
MERGE="${MERGE:-False}"
TRANSLATE="${TRANSLATE:-False}"
SINGLE_FILE="${SINGLE_FILE:-False}"
ADJ_DIR="${ADJ_DIR:-False}"
BY_GENE="${BY_GENE:-False}"
FILTER_FASTA="${FILTER_FASTA:-False}"
DOMAIN="${DOMAIN:-${DOMAIN_DEFAULT}}"
DOMAIN_F="${DOMAIN_F:-${DOMAIN_F_DEFAULT}}"
GENETIC_CODE="${GENETIC_CODE:-${GENETIC_CODE_DEFAULT}}"
MODE="${MODE:-${MODE_DEFAULT}}"
ATTR_MOD="${ATTR_MOD:-${ATTR_MOD_DEFAULT}}"

## if using user-provided reference (i.e. not the built-in A. thaliana accs)
if [ "${MODE}" == "${MODE_DEFAULT}" ] && [ -z "${ACCS}" ]; then
    ACCS='ref'
fi

## throw error if directory not provided
if [ -z "${DIR}" ]; then
    echo "Directory required. Please provide directory using '-d <path to directory>'"
    exit 1
elif [ -z "${ACCS_F}" ] && [ -z "${ACC}" ] && [ -z "${INPUT_DIR}" ]; then
    echo "Accession name(s) or number(s) required. Please provide a file of accession name(s)/number(s) using '-i <path to file>' or provide a single accession name or number using '-a <accession name/number>', or a directory of files of accession numbers (named according to gene ID) using '--input-dir <path to directory>'"
    exit 1
## else if any combination of 2 of ACCS_F, ACC, and INPUT_DIR are provided for some reason
elif (! [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ]) || \
         (! [ -z "${ACCS_F}" ] && ! [ -z "${INPUT_DIR}" ]) || \
         (! [ -z "${ACC}" ] && ! [ -z "${INPUT_DIR}" ]); then
    echo "Please only use either '-a <accession name/number>', '-i <path to file>', or '--input-dir <path to directory>', and not 2 or more at the same time. These parameters are mutually exclusive."
    exit 1
## throw error if GENE or INPUT_DIR or CHROM+START+END are not provided
elif [ -z "${GENE}" ] &&
       ([ "${MODE}" != "${MODE_DEFAULT}" ] || ([ "${MODE}" == "${MODE_DEFAULT}" ] && [ -z "${INPUT_DIR}" ])) && \
       ([ -z "${CHROM}" ] || [ -z "${START}" ] || [ -z "${END}" ]); then
    echo "Gene ID(s) or genome range required. Please provide gene ID(s) using '-g <gene ID(s)>', a directory of files of accession numbers (name according to gene ID) using '--input-dir <path to directory>', or a genome range using '-c <chromosome> -s <start position (inclusive)> -e <end position (inclusive)>'"
    exit 1
## else if any combination of 2 of GENE or INPUT_DIR or CHROM+START+END are provided for some reason
elif (! [ -z "${GENE}" ] && ! [ -z "${INPUT_DIR}" ]) || \
         (! [ -z "${GENE}" ] && (! [ -z "${CHROM}" ] || ! [ -z "${START}" ] || ! [ -z "${END}" ])) || \
         (! [ -z "${INPUT_DIR}" ] && (! [ -z "${CHROM}" ] || ! [ -z "${START}" ] || ! [ -z "${END}" ])); then
    echo "Please only use either '-g <gene ID(s)>', '--input-dir <path to directory>', or '-c <chromosome> -s <start position (inclusive)> -e <end position (inclusive)>', and not 2 or more at the same time. These parameters are mutually exclusive."
    exit 1
## else if genomic range is provided but information is incomplete
elif [ -z "${GENE}" ] && [ -z "${INPUT_DIR}" ] && \
         ([ -z "${CHROM}" ] || [ -z "${START}" ] || [ -z "${END}" ]) && \
         (! [ -z "${CHROM}" ] || ! [ -z "${START}" ] || ! [ -z "${END}" ]); then
    echo "Please provide all of the required parameters to define a genomic range using '-c <chromosome>', '-s <start position(inclusive)> -e <end position (inclusive)>'"
## else if GFF AND BED are provided
elif ! [ -z "${GFF}" ] && ! [ -z "${BED}" ]; then
    echo "Please only use '--gff <GFF3 file>' OR '--bed <GFF3 file converted to BED format by gff2bed>', not both."
fi

## convert values of variables storing file and directory paths to absolute paths
path_vars_to_convert=( "GENE" "ACCS_F" "DIR" "BED" "ACC_LOOKUP" "FOUT" "DOMAIN_F" "RPS_DB" "INPUT_DIR" "FOUT_BED" )
for varname in ${path_vars_to_convert[@]}; do
    if ! [ -z "${!varname}" ] && ([ -f "${!varname}" ] || [ -d "${!varname}" ]); then
        eval ${varname}="$(realpath ${!varname})"
    fi
done

## move to output dir, create if doesn't exist
mkdir -p ${DIR}
cd ${DIR}
DIR=$(pwd)
echo "Output files will be generated in $(pwd)"
## create temporary directory
TMPDIR=$(mktemp -d --tmpdir=${DIR})
mkdir -p ${TMPDIR}
trap 'rm -rf -- "${TMPDIR}"' 0 1 2 15 EXIT ## delete temporary directory on exit or interruption
# trap 'exit 2' 1 2 15 ## exit when interrupted
# mkdir ${TMPDIR}

## filter BED/GFF file
tmp_bed=$(mktemp -p ${TMPDIR})
tmp_gid=$(mktemp -p ${TMPDIR})
# trap 'rm "${tmp_bed}" "${tmp_gid}"' EXIT ## apparently a later trap will override an earlier one
if [ -f "${GENE}" ]; then
    cat ${GENE} > ${tmp_gid}
else
    echo ${GENE} | tr ',' '\n' > ${tmp_gid}
fi
if ! [ -z "${GFF}" ]; then
    BED=${GFF}
fi
/mnt/chaelab/rachelle/src/extract_gff_features.py ${BED} ${tmp_gid} ${tmp_bed} BED
BED=${tmp_bed}

echo "BED: ${BED}  ATTR_MOD: ${ATTR_MOD}"
# echo "$(head ${tmp_bed})"

if [ "${FILTER_FASTA}" == 'True' ]; then
    tmp_ref=$(mktemp -p ${TMPDIR})
    python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); from fasta_manip import fasta_to_dict, dict_to_fasta; from data_manip import splitlines; chroms = set(x.split('\t')[0] for x in splitlines('${tmp_bed}')); seqs = {k: v for k, v in fasta_to_dict('${REFERENCE}').items() if k in chroms}; dict_to_fasta(seqs, '${tmp_ref}')"
    REFERENCE=${tmp_ref}
fi

# start_t=$(date '+%Y-%m-%d %T') ## for combining files
## ## start parsing data
extra=''
if ! [ -z "${QDNAME}" ]; then
    extra+=", qname_dname=${QDNAME}"
fi
if ! [ -z "${QSTARTEND}" ]; then
    extra+=", qstart_qend=${QSTARTEND}"
fi
if ! [ -z "${REFERENCE}" ]; then
    extra+=", ref_fasta_files='${REFERENCE}'"
fi
if ! [ -z "${FOUT}" ] || ! [ -z "${APPEND}" ]; then
    fa_fnames=${DIR}/fa_fnames_$(date '+%s%3N').txt
    extra+=", store_fasta='${fa_fnames}'"
fi
if ! [ -z "${FOUT_BED}" ] || ! [ -z "${NO_BED}" ]; then
    bed_fnames=${DIR}/bed_fnames_$(date '+%s%3N').txt
    extra+=", store_bed='${bed_fnames}'"
fi

## if not using --input-dir
if [ -z "${INPUT_DIR}" ]; then
    ## if accessions are provided (i.e. not requesting ref Col-0), reformat them into new variable $ACCS_INPUT
    if ! [ -z "${ACCS_F}" ] || (! [ -z "${ACC}" ] && ! [ "${ACC}" == 'ref' ]); then
        if [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ]; then
            ACCS_INPUT="('$(echo ${ACC} | sed 's/,/\x27,\x27/g')',)"
        ## else if a file of accessions is provided
        elif ! [ -z "${ACCS_F}" ] && [ -z "${ACC}" ]; then
            sed -i 's/\r//g' ${ACCS_F}
            ACCS_INPUT="('$(awk 'NF > 0' ${ACCS_F} | awk 'NR>1{print PREV} {PREV=$0} END{printf("%s",$0)}' | cat | tr '\n' ',' | sed 's/,/\x27,\x27/g')',)"
        fi
        echo "ACCS_INPUT: ${ACCS_INPUT}"
    fi
    ## if genomic range provided
    if ! [ -z "${CHROM}" ]; then
        ## if 'ref' is provided as $ACC
        if [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ] && [ "${ACC}" == 'ref' ]; then
            python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_ref_by_range('${CHROM}', $(( ${START}-1 )), ${END}, '${DIR}', mode='${MODE}')"
        else
            python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_1001pseudogenome_by_range(${ACCS_INPUT}, '${CHROM}', ${START}, ${END}, '${DIR}')"
        fi
    ## else if genes provided
    else
        if [ -f "${GENE}" ]; then
            readarray -d ',' -t parsed_genes <<< "$(tr '\n' ',' < ${GENE} | sed 's/,\+$//g')"
        else
            readarray -d , -t parsed_genes <<< "$(echo ${GENE} | sed 's/\n//g')"
        fi
        for parsed_gene in "${parsed_genes[@]}"; do
            ## if 'ref' is provided as $ACC
            if [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ] && [ "${ACC}" == 'ref' ]; then
                python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_ref_by_gene(\"\"\"${parsed_gene}\"\"\".replace('\n',''), '${FEATURE}', '${DIR}', bed='${BED}', complete=(${COMPLETE}), domain='${DOMAIN}', domain_f='${DOMAIN_F}', merge=(${MERGE}), translate=(${TRANSLATE}), adj_dir=(${ADJ_DIR}), by_gene=(${BY_GENE}), attribute_mod=${ATTR_MOD}, genetic_code='${GENETIC_CODE}', mode='${MODE}'${extra})"
            else
                python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_1001pseudogenome_by_gene(${ACCS_INPUT}, \"\"\"${parsed_gene}\"\"\".replace('\n',''), '${FEATURE}', '${DIR}', bed='${BED}', acc_ref_file='${ACC_LOOKUP}', header=(${HEADER}), complete=(${COMPLETE}), domain='${DOMAIN}', domain_f='${DOMAIN_F}', adj_dir=(${ADJ_DIR})${extra})"
            fi
        done
    fi
else
    for file in ${INPUT_DIR}/AT?G*; do
        fname=$(basename -- "${file}")
        GENE="${fname%.*}"
        sed -i 's/\r//g' ${file}
        ACCS_INPUT="('$(awk 'NF > 0' ${file} | awk 'NR>1{print PREV} {PREV=$0} END{printf("%s",$0)}' | cat | tr '\n' ',' | sed 's/,/\x27,\x27/g')')"
        python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); import get_seqs_functions; get_seqs_functions.get_1001pseudogenome_by_gene(${ACCS_INPUT}, '${GENE}', '${FEATURE}', '${DIR}', bed='${BED}', acc_ref_file='${ACC_LOOKUP}', header=(${HEADER}), complete=(${COMPLETE}), domain='${DOMAIN}', domain_f='${DOMAIN_F}', adj_dir=(${ADJ_DIR})${extra})"
    done
fi

## cat .fasta files if output file path provided
if ( ! [ -z "${FOUT}" ] || ! [ -z "${APPEND}" ] ) && [ -f ${fa_fnames} ]; then
    ## if using FOUT
    if ! [ -z "${FOUT}" ]; then
       if [ $(dirname "${FOUT}") == '.' ]; then
           fout=${DIR}/${FOUT}
       else
           fout=${FOUT}
       fi
       cat $(sort ${fa_fnames} | uniq) > ${fout}
       echo "${num_files} FASTA files concatenated into ${fout}"
    fi
    ## if using APPEND
    if ! [ -z "${APPEND}" ]; then
        if [ $(dirname "${APPEND}") == '.' ]; then
            append=${DIR}/${APPEND}
        else
            append=${APPEND}
        fi
        ## if file doesn't end with newline, add newline
        if [ -f "${append}" ] && ! [[ $(tail -c1 ${append} | wc -l) -gt 0 ]]; then
            printf '\n' >> ${append}
        fi
        cat $(sort ${fa_fnames} | uniq) >> ${append}
        echo "${num_files} FASTA files concatenated into ${append}"
    fi
    # cat $(sort ${fa_fnames} | uniq) > ${fout}
    num_files=$(wc -l ${fa_fnames})
    rm $(sort ${fa_fnames} | uniq)
    rm ${fa_fnames}
fi

## delete .bed files if no-bed raised
if ! [ -z "${NO_BED}" ]; then
    rm $(sort ${bed_fnames} | uniq)
    rm ${bed_fnames}
    find ${DIR}/bed -type d -empty -delete ## remove ./bed if empty
    # rm $(find ${DIR}/bed/*.bed -maxdepth 1 -type f -newermt "${start_t}")
    # find ${DIR}/bed -type d -empty -delete ## remove ./bed if empty
elif ! [ -z "${FOUT_BED}" ] && [ -f ${bed_fnames} ]; then ## cat .bed files if output .bed path provided
    if [ $(dirname "${FOUT_BED}") == '.' ]; then
        fout=${DIR}/${FOUT_BED}
    else
        fout=${FOUT_BED}
    fi
    cat $(sort ${bed_fnames} | uniq) > ${fout}
    rm $(sort ${bed_fnames} | uniq) ## previously, all $(sort ${xx_fname} | uniq) were $(cat ${xx_fname})
    rm ${bed_fnames}
    find ${DIR}/bed -type d -empty -delete ## remove ./bed if empty
    echo "BED files concatenated into ${fout}"
fi

exit 0
