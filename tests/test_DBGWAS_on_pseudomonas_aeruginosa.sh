set -eux

#run the tested version on pseudomonas_aeruginosa_full_dataset
wget https://www.dropbox.com/s/0g1llvdbfv1jys6/pseudomonas_aeruginosa_full_dataset.zip?dl=1 -O pseudomonas_aeruginosa_full_dataset.zip
unzip pseudomonas_aeruginosa_full_dataset.zip
wget https://www.dropbox.com/s/mt3g4oh0bt5jwmr/Resistance_DB_for_DBGWAS.fasta?dl=1 -O Resistance_DB_for_DBGWAS.fasta
wget https://www.dropbox.com/s/9y1p0yw918ips6k/uniprot_sprot_bacteria_for_DBGWAS.fasta?dl=1 -O uniprot_sprot_bacteria_for_DBGWAS.fasta
./DBGWAS -strains pseudomonas_aeruginosa_full_dataset/strains -newick pseudomonas_aeruginosa_full_dataset/strains.newick -nc_db Resistance_DB_for_DBGWAS.fasta -pt_db uniprot_sprot_bacteria_for_DBGWAS.fasta

#get the correct output
wget https://www.dropbox.com/s/jp24184zony1977/correct_output.zip?dl=1 -O correct_output.zip
unzip correct_output.zip

#compare both output
diff -rq correct_output output | grep -v ".png" | grep -v ".log.txt" | grep -v "visualisations/index.html"

# If this previous command does not output anything, then it is fine
