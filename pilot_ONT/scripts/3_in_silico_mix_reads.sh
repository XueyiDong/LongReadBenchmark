cd /vast/scratch/users/dong.x/long_read_benchmark/pilot_ONT/

mkdir -p fq_mix
cat fq_subsample/barcode01.025.fq.gz fq_subsample/barcode05.075.fq.gz > fq_mix/rep1-075.fq.gz
cat fq_subsample/barcode01.050.fq.gz fq_subsample/barcode05.050.fq.gz > fq_mix/rep1-050.fq.gz
cat fq_subsample/barcode01.075.fq.gz fq_subsample/barcode05.025.fq.gz > fq_mix/rep1-025.fq.gz
cat fq_subsample/barcode06.025.fq.gz fq_subsample/barcode10.075.fq.gz > fq_mix/rep2-075.fq.gz
cat fq_subsample/barcode06.050.fq.gz fq_subsample/barcode10.050.fq.gz > fq_mix/rep2-050.fq.gz
cat fq_subsample/barcode06.075.fq.gz fq_subsample/barcode10.025.fq.gz > fq_mix/rep2-025.fq.gz
