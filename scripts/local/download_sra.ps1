prefetch --option-file "./references/srr_id_list.txt" --output-directory "sra_files"
[string[]]$arrayFromFile = Get-Content -Path './references/srr_id_list.txt'

for ($srr_id in $arrayFromFile) {
    fasterq-dump --outdir "./sra_files" $srr_id
}