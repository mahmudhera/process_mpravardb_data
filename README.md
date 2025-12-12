# process_mpravardb_data

```bash
conda create -n basics numpy pandas biopython -c bioconda -c conda-forge -y
conda activate basics
python process_mpravardb_csv.py K562.csv K562_processed.csv /projects/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.fa /projects/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
python test_processed_file.py --input_csv K562_processed.csv
```