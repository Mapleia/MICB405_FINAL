# !/bin/bash

# list the files in the directory as an array of string
declare -A fileGroups

# for each file, group file by string before _
# then create an array of them

LEAD_SRA_DIR="./test_sra_script"
PROJECT_DIR="."

for file in $LEAD_SRA_DIR/*; do
    base=$(basename "$file")
    IFS='_' read -ra group <<< "$base"
    fileGroups["$group"]+="$file "
done

star_group_path="$PROJECT_DIR/alignments"
echo "STAR group folder: $star_group_path"

for i in "${!fileGroups[@]}"; do
    group=${fileGroups[$i]}
    echo "Group: $group"

    bioproject="PRJNA576974"

    bioproject_arr=("PRJNA608616" "PRJNA735693" "PRJNA335844" "PRJNA394256" "PRJNA576974")
    srr_group_arr=("SRR112495" "SRR147469" "SRR3987" "SRR59882" "SRR102697")

    for str in ${!bioproject_arr[@]}; do
        if [[ $group == *${srr_group_arr[$str]}* ]]; then
            bioproject=${bioproject_arr[$str]}
        fi
    done

    echo "Bioproject: $bioproject"

    mkdir -v -p $star_group_path/$bioproject # -m 770


    IFS=' ' read -r -a array <<< "${fileGroups[$i]}"
    for item in "${array[@]}"; do
        echo "  File: $item"
    done

    # run STAR
    # STAR \
    # --genomeDir $STAR_INDEX_DIR \
    # --readFilesIn $array \
    # --readFilesCommand zcat \
    # --outFileNamePrefix "$star_group_path/$bioproject/$group/star_${group}_sample_" \
    # --outSAMtype BAM SortedByCoordinate \
    # --runThreadN 16 \
    # --quantMode GeneCounts
done