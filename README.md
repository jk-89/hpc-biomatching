## Biomatching with CUDA

This repository contains an implementation of genome variant matching against a genome database, simulating the diagnostic process that matches many (thousands) of genome variants to the same database.

This project implements a solution for this problem using CUDA, with the primary goal of assessing whether GPU operations can speed up the matching process compared to a standard CPU-based approach.

## Results

Benchmark results, file descriptions, and a detailed report can be found in the `report.pdf` file.

## How to run
1. Specify the path to your CUDA Compiler in the second line of the `Makefile`.
2. Run the code to index the database.
3. Once indexing is complete, proceed to match the provided genomes.
```
make
./gpugenv -i path/to/database.tsv path/to/your/output/index
./gpugenv path/to/variant.vcf path/to/your/output/index path/to/matched/output.tsv
```
