bcbio_prepare_samples.py --out fastq --csv samples.csv
bcbio_nextgen.py -w template --only-metadata template.yaml samples-merged.csv fastq/*gz
cd samples-merged/work
# if you want to run in parallel or using a cluster read docs: https://bcbio-nextgen.readthedocs.io/en/latest/contents/parallel.html
bcbio_nextgen.py ../samples-merged.yaml
