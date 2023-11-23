# !/bin/bash
while IFS= read -r line; do
  echo "Downloading from $line..."
  wget -nc $line
done < "ftp_list.txt"
