#This file contains all scripts run for the publication "Microbial iron metabolism as revealed by gene expression profiles in contrasted Southern Ocean regimes". Pavla Debeljak, Eve Toulza, Sara Beier, Stephane Blain, Ingrid Obernosterer. The sequence data have been submitted to the EMBL databases under accession number PRJEB30315. The authors thank the Roscoff Bioinformatics platform ABiMS (http://abims.sb-roscoff.fr) for providing computational resources. Only final databases are provided (NCBI&GOS&KEGG curated), previous versions of databases upon request via pavla.debeljak@gmail.com

# subset blasts on all reads for 1% of sequences
fasta-subsample HDC-1b_nonrrna_Format_Fasta.fasta 318724
fasta-subsample HDC-2b_nonrrna_Format_Fasta.fasta 315194
fasta-subsample HDC-3b_nonrrna_Format_Fasta.fasta 360706
fasta-subsample HDC-4b_nonrrna_Format_Fasta.fasta 320999
fasta-subsample HDC-5b_nonrrna_Format_Fasta.fasta 298242
fasta-subsample HDC-6b_nonrrna_Format_Fasta.fasta 346207
# diamond blast for subsets
#blastx align translated DNA query sequences against a protein reference database
diamond blastx -d nr.dmnd -q subset1b.fasta -a subset1b -t diamond/tmp -k 1 -e 10 -p 12 &&
#generate formatted output (xml=5) from DAA files
diamond view -a subset1b.daa -o ~/diamond/subset/subset1b.xml -outfmt 5


# blast p of each database against GOS protein database

makeblastdb -in GOS_short.fa -dbtype 'prot' -out GOS
blastp+ -query FE_genes.fasta -db GOS -out FE_genes_GOS -evalue 0.001 -max_target_seqs 25 -outfmt 6
blastp+ -query Iso_NCBIonly.fasta -db GOS -out Iso_NCBI_GOS -evalue 0.001 -max_target_seqs 25 -outfmt 6
blastp+ -query Aco_NCBIonly.fasta -db GOS -out Aco_NCBI_GOS -evalue 0.001 -max_target_seqs 25 -outfmt 6


# blast GOS sequences against KEGG
## for top e-value
awk '!x[$1]++' FS="\t" Iso_NCBI_GOS > Iso_tophit
awk -F "\t" '$12<0.00001' KEGGbluni.FE_GOS_final.tab > test
cd /data/sara/KEOPS

#blast searches
nohup blastp -db /data/KEGG/official_KEGG/19052016/pro_eu.with.ko.pep -query Aco_GOS.fa -outfmt 6 -num_alignments 1 -num_threads 10 -out blastout.Aco.tab &

nohup blastp -db /data/KEGG/official_KEGG/19052016/pro_eu.with.ko.pep -query Iso_GOS.fa -outfmt 6 -num_alignments 1 -num_threads 10 -out blastout.Iso.tab &

nohup blastp -db /data/KEGG/official_KEGG/19052016/pro_eu.with.ko.pep -query FE_GOS_final.fa -outfmt 6 -num_alignments 1 -num_threads 10 -out blastout.FE_GOS_final.tab &

#sort out unique hits
cat blastout.Aco.tab|sort -nk12,12gr -nk3,3gr -nk4,4gr -t $'\t' |awk '!x[$1]++' FS="\t" >bluni.Aco.tab
cat blastout.Iso.tab|sort -nk12,12gr -nk3,3gr -nk4,4gr -t $'\t' |awk '!x[$1]++' FS="\t" >bluni.Iso.tab
cat blastout.FE_GOS_final.tab|sort -nk12,12gr -nk3,3gr -nk4,4gr -t $'\t' |awk '!x[$1]++' FS="\t" >bluni.FE_GOS_final.tab

#match blast output with KEGG identifiers
perl -e ' $col1=0; $col2=1; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' /data/KEGG/official_KEGG/19052016/genes_ko2.list bluni.Aco.tab|cut -f 2,4,5,6,7,8,9,10,11,12,13,14,15 >KEGGbluni.Aco.tab


perl -e ' $col1=0; $col2=1; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' /data/KEGG/official_KEGG/19052016/genes_ko2.list bluni.Iso.tab|cut -f 2,4,5,6,7,8,9,10,11,12,13,14,15 >KEGGbluni.Iso.tab

perl -e ' $col1=0; $col2=1; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' /data/KEGG/official_KEGG/19052016/genes_ko2.list bluni.FE_GOS_final.tab|cut -f 2,4,5,6,7,8,9,10,11,12,13,14,15 >KEGGbluni.FE_GOS_final.tab

###After Kegg blast uniq on KEGG Id and grep on only these (look in supplementary table 5).
Then GOS Id taken out for perl on .fa file 

awk 'NR%2==0' Iso_ordered.fa | paste -d'\n' Iso_Affil_done.txt - > Iso_onlyGOS_final.fa

#Kegg Id in same order as file output from blast as text file > in beginning
#to order seqeunces:
python reorder_fasta.py Iso_GOS_KEGG.fa KEGG_Iso_GOSn
#sequences in right order
#replace fasta names with names from text file
awk 'NR%2==0' Iso_ordered.fa | paste -d'\n' Iso_Affil_done.txt - > Iso_onlyGOS_final.fa
#finally
cat Iso_NCBIonly.fasta Iso_onlyGOS_final.fa > Iso_DB_NCBI_GOS.fa
#
#Aconitase
python reorder_fasta.py Acoonly_GOS_KEGG.fa Acoonly_Aff.txt
awk 'NR%2==0' Aco_ordered.fa.txt | paste -d'\n' Aco_FinalAff.txt - > Aco_onlyGOS_final.fa
cat Aco_NCBIonly.fasta Aco_onlyGOS_final.fa > Aco_DB_NCBI_GOS.fa
#
#FE_Database
python reorder_fasta.py FE_GOS_KEGG_finalnew.fa FE_GOS_KEGG.txt
awk 'NR%2==0' New_order_FE_GOS.fa.txt | paste -d'\n' KEGG_AFF_fasta2.txt - > FE_GOS_KEGG_final_new.fa
cat FE_NCBI_new.fa FE_GOS_KEGG_final_new.fa > FE_DB_NCBI_GOS.fa
#
#for Aconitase delimiter was after second "," and then from 3 on
cut -f 2 HDC_1b_Aco.tab > 1b_2
cut -d',' -f 3 1b_2 > 1b_species
sort 1b_species > 1b_Spe_sorted
uniq -c 1b_Spe_sorted > 1b_Aco_count

#for Isocitrate Lyase delimiter was after "_"
cut -f 2 HDC_1b_Iso.tab > 1b_2
cut -d'_' -f 3-7 1b_2 > 1b_species
sort 1b_species > 1b_Spe_sorted
uniq -c 1b_Spe_sorted > 1b_Iso_count

#for FE I was much smarter and had exact Identifier so procesuder as following for each HDC file (diamond tab result file
cut -f 2 HDC_1b_FE.tab > 1b_2_FE
grep ":FL:" 1b_2_FE > 1b_FL
grep ":ST:" 1b_2_FE > 1b_ST
grep ":F2:" 1b_2_FE > 1b_F2
grep ":F3:" 1b_2_FE > 1b_F3
grep ":SU:" 1b_2_FE > 1b_SU

#then for each
cut -d ':' -f 3 1b_F3 | sort | uniq -c > 1b_F3_Phylum
cut -d ':' -f 4 1b_F3 | sort | uniq -c > 1b_F3_Class
cut -d ':' -f 5 1b_F3 | sort | uniq -c > 1b_F3_Order
cut -d ':' -f 6 1b_F3 | sort | uniq -c > 1b_F3_Family
cut -d ':' -f 7 1b_F3 | sort | uniq -c > 1b_F3_Genus
