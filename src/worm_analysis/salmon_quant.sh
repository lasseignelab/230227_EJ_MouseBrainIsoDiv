#!/bin/bash
cd /data/project/lasseigne_lab/EmmaJones/worm_data/rawData

for fn in SRR36572{24..33};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i /data/project/lasseigne_lab/GENOME_dir/C_elegans_WBPS18/c_elegans_index/ -l A \
         -r ${samp}.fastq \
         -p 8 --validateMappings \
         -o /data/user/efjones/230227_EJ_MouseBrainIsoDiv/data/salmon/${samp}_quant
done