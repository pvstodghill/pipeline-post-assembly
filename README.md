# pipeline-post-assembly

Post-assembly pipeline for genome QC and annotation.

## Cloning the repo

This pipeline uses Git submodules. The easiest way to clone this repo
(with a recent version of `git`) is

```
git clone --recurse-submodules https://github.com/pvstodghill/pipeline-post-assembly.git
```

## Installing prereqs

You will want to install the following:

- [Snakemake](https://snakemake.github.io/)
- [Conda](https://conda.io)
- [Perl](https://www.perl.org/)

## Configuring the pipeline

To run the pipeline on your own data,

1. Copy `config.template.yaml` to `config.yaml`.  Edit `config.yaml`
   according to your needs and local environment.

## Running the pipeline

`snakemake --use-conda`

