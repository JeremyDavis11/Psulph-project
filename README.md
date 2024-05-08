# Undergraduate-Thesis-JDavis
Script documentation for a undergraduate bioinformatics project

# Thesis project readme

## 01_mapping

mapping raw Pacbio read files to a *Poecilia mexicana* reference genome.

```
#!/bin/bash
#SBATCH --job-name=map 
#SBATCH --time=2-00:00:00
#SBATCH --partition= 
#SBATCH --mail-user= 
#SBATCH --mail-type=
#SBATCH --output=Map_psulph-to-mex_%j.out 
#SBATCH --error=Map_psulph-to-mex_%j.err 
#SBATCH --ntasks=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=20GB 

module load minimap2/m2-2.17 hb hb-gnu

#Submit as sbatch 01_map.sh genome reads output_sam 
#-ax asm20 is for hifi reads
```
---
## 02_flagstat

obtain percent of genome mapped from a bam file 

```
#!/bin/bash
#SBATCH --job-name=samtools-flagstat # Name the job
#SBATCH --time=0-10:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --partition=Instruction # Run on the instruction partition
#SBATCH --mail-user=jealdavi@ucsc.edu # Send updates to email
#SBATCH --mail-type=ALL # send all types of updates
#SBATCH --output=out/flagstat_Psulph_%j.out # output file
#SBATCH --error=err/flagstat_Psulph_%j.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=8 # 2 CPUs per job
#SBATCH --mem=20GB # Memory limit 20GB

module load samtools

samtools flagstat /hb/home/jealdavi/thesis_project/01_mapping/bam-files/map-psulph-to-mex_sorted.bam
```
---

##03_Convert
Converting the sam moutput file from minimap to a sorted and indexed bam file compatable with Sniffles2.

``` 
#!/bin/bash
#SBATCH --job-name=convert 
#SBATCH --time=2-00:00:00 
#SBATCH --partition=128x24 
#SBATCH --mail-user= 
#SBATCH --mail-type=
#SBATCH --output=convert_sulph-2-mex_%j.out 
#SBATCH --error=convert_sulph-2-mex_%j.err
#SBATCH --ntasks=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=20GB 

module load samtools/samtools-1.17

samtools view -bS ${1} > ${1}.bam

samtools sort ${1}.bam -o ${1}_sorted.bam

samtools index ${1}_sorted.bam`
```

---
##04_variant-calling

Sniffles script that calls variant between bam files. Outputs a vcf file with variant information

```
#!/bin/bash
#SBATCH --partition=128x24
#SBATCH --job-name=sniffles
#SBATCH --output=out/%x_%j.out
#SBATCH --error=err/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --ntasks=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G

module load sniffles/sniffles-2.3.2 hb-gnu

sniffles --input /hb/home/jealdavi/thesis_project/01_mapping/bam-files/map-pmex-ref_sorted.bam --vcf Pmex2Ref_sniffles.vcf.gz --threads 8 --reference /hb/groups/kelley_lab/poeciliids/hifi_genomes/02_nuc_assembly/02_HiFiasm/assemblies/m84066_231208_213947_s3.hifi_reads.bc2031.asm.bp.p_ctg.fa --minsvlen 4
```
---
##05_intersect

screens for and displays the overlap between genomic features between the genome and reference. This script was run with 3 seperate vcf files containing variant information on exon, introns, and genes respectively. This script was also run for each vcf file with the [-v] option that writes the original bed file entry for each overlap to count SVs in intergenic regions.

```
#!/bin/bash
#SBATCH --job-name=map 
#SBATCH --time=2-00:00:00 
#SBATCH --partition=128x24 
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --output=out/intersect_intron_%j.out 
#SBATCH --error=err/intersect_intron_%j.err 
#SBATCH --ntasks=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=20GB 

module load bedtools/bedtools-2.26.0 hb hb-gnu 

bedtools intersect [-v] -a /hb/home/jealdavi/thesis_project/03_intersect/Pmex_NS_braker_inexon.bed -b /hb/home/jealdavi/thesis_project/03_intersect/bedfiles/intron_regions.bed
```
---
## 06_unique-SVs

filter overlap output files for unique lines as to not count SVs that overlap with multiple features


```
#!/bin/bash
#SBATCH --job-name=uniq
#SBATCH --time=2-00:00:00 
#SBATCH --partition=128x24 
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --output=out/uniq-SV_%j.out 
#SBATCH --error=err/uniq-SV_%j.err 
#SBATCH --ntasks=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=8 
#SBATCH --mem=20GB 

cat | grep -oP ‘SVTYPE=[^;]*’ | sort | uniq -c 

```

