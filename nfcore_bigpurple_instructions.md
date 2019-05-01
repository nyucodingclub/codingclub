## Running nf-core pipelines on BigPurple
nf-core provides easy to use bioinformatics pipelines like: 
| pipeline            | description |
|---|:---:|
| nf-core/bcellmagic  | B cell repertoire analysis pipeline with immcantation framework.|
| nf-core/nascent     | Nascent Transcription Processing Pipeline |
| nf-core/atacseq     | ATAC-seq peak-calling and differential analysis pipeline.|
| nf-core/rnafusion   | RNA-seq analysis pipeline for detection gene-fusions |
| nf-core/rnaseq      | RNA sequencing analysis pipeline using STAR or HISAT2, with gene counts and quality control |
| nf-core/hlatyping   | Precision HLA typing from next-generation sequencing data |
| nf-core/eager       | A fully reproducible and state of the art ancient DNA analysis pipeline. |
| nf-core/mhcquant    | Identify and quantify peptides from mass spectrometry raw data |
| nf-core/methylseq   | Methylation (Bisulfite-Sequencing) analysis pipeline using Bismark or bwa-meth + MethylDackel  |
| nf-core/ampliseq    | 16S rRNA amplicon sequencing analysis workflow using QIIME2 |
| nf-core/deepvariant | Google's DeepVariant variant caller as a Nextflow pipeline |
All of their pipelines with links to their individual docs can be found here: [https://nf-co.re/pipelines](https://nf-co.re/pipelines)

To run any of the nf-core pipelines, you first have to install nextflow in your home directory:
```
cd
curl -s https://get.nextflow.io | bash 
```
The first time you're running a specific pipeline (e.g. rnaseq), run it in an interactive session:
```
srun --cpus-per-task=4 --time=03:00:00
```
You won't have to do this again for future runs of the same pipeline, just for the first time.
Once the interactive session starts, you have to load two modules:
```
module load singularity/3.1 squashfs-tools/4.3
```
Perfect, now you can run your pipeline!
I recommend doing this from your scratch directory:
```
cd /gpfs/scratch/$USER
```
Let's run for example the rnaseq pipeline (for info on all the parameters, see the [documentation](https://github.com/nf-core/rnaseq/blob/master/docs/usage.md)):
```
nextflow run nf-core/rnaseq -profile bigpurple --reads 'path/to/data/sample_*_{1,2}.fastq' --fasta my_reference.fasta --gtf my_reference.gtf --outdir myresults
```
And off you go! The first time you run this, it will take some time to pull the docker container and convert it to a singularity container, which is what our cluster is using. Once you have run the pipeline once, the singularity container will be cached in the directory `/gpfs/scratch/$USER/singularity_images_nextflow` and nextflow will not download it again unless it gets deleted. As a result you can run the pipeline from the login node now :D, no need for an interactive session.
The results will be in a folder called `myresults`.
 
This pipeline also let's you use some often used genomes like for human, drosophila and yeast, without having to specify a fasta reference and a corresponding gtf, you can find more info in the [documentation](https://github.com/nf-core/rnaseq/blob/master/docs/usage.md)), including how to use HISAT2 instead of the default aligner which is STAR.

Notice also how we specified `-profile bigpurple` - that makes sure the pipeline runs smoothly on bigpurple. At the time of writing, this does not work for all of the pipelines, just try it out, the pipeline will throw you an error telling you the that the `profile bigpurple is not defined`. Not to worry though, you can always run nf-core pipelines with the config file for bigpurple. The rnaseq pipeline does this automatically when you specify `-profile bigpurple` and eventually all of the pipelines will follow suit. Until then, here is the work-around:
```
git clone https://github.com/nf-core/configs.git /gpfs/scratch/$USER/nf-core/configs
```
Now you can run the nf-core pipelines like this:
```
nextflow run nf-core/<name of the pipeline> -c /gpfs/scratch/$USER/nf-core/configs/conf/bigpurple.config <additional pipeline parameters>
```