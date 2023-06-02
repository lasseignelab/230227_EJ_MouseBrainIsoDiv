# open and run docker
docker run -d -it \
--name run_pfam \
--mount type=bind,source=/data/user/efjones/230227_EJ_MouseBrainIsoDiv/data/switchlist_fasta,target=/opt/PfamScan/fasta,readonly \
--mount type=bind,source=/data/user/efjones/230227_EJ_MouseBrainIsoDiv/data/Pfam,target=/opt/PfamScan/results \
pfam_test:1.0.36

# run actual script
cd /opt/PfamScan/

perl pfam_scan.pl \
-fasta /data/user/efjones/230227_EJ_MouseBrainIsoDiv/data/switchlist_fasta/region_region_AA.fasta \
-dir ~/PfamScan/PfamFiles \
| tee /data/user/efjones/230227_EJ_MouseBrainIsoDiv/results/PfamScan/region_region_pfamscan_out.txt