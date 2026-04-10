# Overview

Polars pipeline for Maf merging and deduplication from what to me looks like a potential datalake structure into parquet files for easy checkpointing and reading for later. The Polars streams straight from compressed Maf files after reading some metadata files to determine which ones
should be kept and syncs them back to disk. This keeps RAM usage down and was perfect for Nextflow in the other project. Since the original merge is
synced to Parquet subsiquent reads are faster and the streaming engine of Polars can take advantage of the columnar format and only read the columns
it needs for analysis/transformation. 

Originally the plan was to use Docker DIND rootless as a way to orchestrate this with Jenkins but other options like K3s and Nextflow and AWS were explored. Pixi was used as a package manager to deal with software drift but this was traded out for Nix in the later project for containerization rather than Docker Compose. I will probably revisit K3s and Docker DIND rootless but will focus on AWS S3, Lambda, and Batch for now.