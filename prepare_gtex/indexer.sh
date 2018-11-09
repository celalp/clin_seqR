#! /bin/bash


databases=$(ls | grep ".db$")

for database in $databases
   do
     sqlite3 $database "create index subject_index on subjects ('SUBJID')"
     sqlite3 $database "create index sample_index on samples ('SAMPID')"
     sqlite3 $database "create index gene_reads_index on gene_reads ('Description')"
     sqlite3 $database "create index gene_tpm_index on gene_tpm ('Description')"
     sqlite3 $database "create index transcript_read_index on transcript_reads('gene_id')"
     sqlite3 $database "create index transcript_tpm_index on transcript_tpm('gene_id')"
     sqlite3 $database "create index exon_reads_index on exon_reads ('gene_id')"
     sqlite3 $database "create index junctions_index on junctions ('gene_id')"

     echo "$databse done"
done

echo "all done"
