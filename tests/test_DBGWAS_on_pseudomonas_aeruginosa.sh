set -eux

#run the tested version on pseudomonas_aeruginosa_full_dataset
wget https://www.dropbox.com/s/0g1llvdbfv1jys6/pseudomonas_aeruginosa_full_dataset.zip?dl=1 -O pseudomonas_aeruginosa_full_dataset.zip
unzip pseudomonas_aeruginosa_full_dataset.zip
wget https://www.dropbox.com/s/mt3g4oh0bt5jwmr/Resistance_DB_for_DBGWAS.fasta?dl=1 -O Resistance_DB_for_DBGWAS.fasta
wget https://www.dropbox.com/s/9y1p0yw918ips6k/uniprot_sprot_bacteria_for_DBGWAS.fasta?dl=1 -O uniprot_sprot_bacteria_for_DBGWAS.fasta
./DBGWAS -strains pseudomonas_aeruginosa_full_dataset/strains -newick pseudomonas_aeruginosa_full_dataset/strains.newick -nc-db Resistance_DB_for_DBGWAS.fasta -pt-db uniprot_sprot_bacteria_for_DBGWAS.fasta

#get the correct output

#to compare with v0.5.2:
#wget https://www.dropbox.com/s/pr5vn76xksdbtzj/correct_output_v0.5.2.zip?dl=1 -O correct_output.zip

#to compare with v0.5.4:
wget https://www.dropbox.com/s/akzz2jor9pmw6yo/correct_output_v0.5.4.zip?dl=1 -O correct_output.zip

unzip correct_output.zip

#compare both outputs
echo "Differences found (check them to see if they are really a problem):"
diff -rq correct_output output | grep -v ".png" | grep -v ".log.txt" | grep -v "nucl_db_fixed.nin" | grep -v "prot_db_fixed.pin"
