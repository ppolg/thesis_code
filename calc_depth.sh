#! /bin/bash -p

bam_dir=$1

extension=".tsv"
shopt -s nullglob
for bam_path in "$bam_dir"/*.bam; do
    bam_file=${bam_path##*/}
    bam_name=${bam_file%.bam}
    echo "Calculating depth for $bam_file"
    samtools depth -o "$bam_name$extension" "$bam_file"
    echo "saved as $bam_name$extension"
done
echo "All done! :)"