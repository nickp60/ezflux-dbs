
### Getting data and building Sylph db


The AGORA2 metadata list ~7k genomes, which are a mixure of genbank, refseq and unassambled reads. the "Genome Link" column has some accessions with full FTP links, wheras others just have the link to the assembly webpage.  The script `./scripts/split_accs_and_sras.R` fetches and cleans up the file, outputting two files: `agora2_accs` which contains the assembly accesions, and `agora2_sras` which contains the SRA accessions.

```
snakemake --snakefile agora2/01-get_agora2_genomes/Snakefile --config agora2_sras=$PWD/agora2/01-get_agora2_genomes/agora2_sras agora2_accs=$PWD/agora2/01-get_agora2_genomes/agora2_accs --directory /data1/collab004/resources/dbs/sylph/agora2_clean_full/ --rerun-triggers mtime  --rerun-incomplete -n
```

## Note on Failures

this successfully downloads 7182 genomes; some may have been removed from NCBI, causing the descrepency between that number and whats listed in the papers supplementary


### Creating 16S DB


```
snakemake --snakefile agora2/02-make_agora2_blastdb/Snakefile --config agora_genomes=/data1/collab004/resources/dbs/sylph/agora2_clean_full/genomes/ dbname="agora2_16S"   --directory /data1/collab004/resources/dbs/sylph/agora2_clean_full/ --rerun-triggers mtime  --rerun-incomplete -n

```
