# read file in ./references/..._srr_id_list.txt as an array

$project_id_file='./references/PRJNA335844_srr_id_list.txt'

[string[]]$arrayFromFile = Get-Content -Path $project_id_file

# create a file for each SRR id in ./references/..._srr_id_list.txt
for ($i=0; $i -lt $arrayFromFile.Length; $i++) {
    $srr_ftp_path = $arrayFromFile[$i]
    $srr_id = $srr_ftp_path.Split("/")[-1]
    # remove the .gz extension
    $srr_id = $srr_id.Split(".")[0]
    echo $srr_id
    $file_name = "./test_sra_script/$srr_id.fastq"
    New-Item -Path $file_name -ItemType File
}
echo $arrayFromFile
