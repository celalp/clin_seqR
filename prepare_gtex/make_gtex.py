import argparse
import pandas as pd
import os
import sqlalchemy as sql
from multiprocessing import Pool
import sqlalchemy_utils as utils
#import ruamel.yaml as yam #unused for now because I need to set other things first


# TODO this needs some serious refactoring
# TODO add yaml for config options this will all go in to a giant shell script to set things up


def get_sample_names(tiss, sampleinfo, tissue_col, eng, filter, filter_col=None,
                     to_db=True, ret=False, schema=None):
    '''Using the path of the sample info file get the names of the samples'''
    all_samples = pd.read_csv(sampleinfo, header=0, sep="\t")
    tissue_samples=all_samples.loc[all_samples[tissue_col] == tiss]
    if filter_col is not None:
        tissue_samples=tissue_samples.loc[tissue_samples[filter_col] == filter]
    if to_db:
        if schema is not None:
            tissue_samples.to_sql("samples", eng, if_exists='fail',schema=schema, chunksize=5000)
        else:
            tissue_samples.to_sql("samples", eng, if_exists='fail', chunksize=5000)
    if ret:
        return(tissue_samples.iloc[:,0])

def get_subject_info(samplenames, subjectinfo, column, eng, to_db=True, ret=False, schema=None):
    '''Using the path of the subject info file get the subject properties'''
    all_subjects = pd.read_csv(subjectinfo, header=0, sep="\t")
    subjects=list()
    for sample in samplenames:
        subject = "GTEX-"+sample.split("-")[1] #hard coded for v7 but it's unlikely to change
        subjects.append(subject)
    subj_info=all_subjects[all_subjects[column].isin(subjects)]
    if to_db:
        if schema is not None:
            subj_info.to_sql("subjects", eng, if_exists='fail', chunksize=5000, schema=schema)
        else:
            subj_info.to_sql("subjects", eng, if_exists='fail', chunksize=5000)
    if ret:
        return(subj_info.iloc[:,0])


def get_gene_info(samplenames, gene_reads, gene_tpm, eng, to_db=True, schema=None):
    '''get gene info, tpm and read count from rnaseqc, uses file paths
    this function gets both reads and tpm for genes, I am not bothering with medians'''
    cols=samplenames.copy()
    cols.extend(["Name", "Description"])
    files={'gene_reads':gene_reads, 'gene_tpm':gene_tpm}
    for file in files.keys():
        for chunk in pd.read_csv(files[file], header=0, sep="\t", skiprows=2, usecols=cols, chunksize=5000):
            chunk.columns.values[[0,1]] = ["gene_id", "gene_name"]
            if to_db:
                if schema is not None:
                    chunk.to_sql(file, eng, if_exists='append', chunksize=5000, schema=schema)
                else:
                    chunk.to_sql(file, eng, if_exists='append', chunksize=5000)


def get_transcript_info(samplenames, transcript_reads, transcript_tpm, eng, to_db=True, schema=None):
    '''get transcript info tpm and read expected count from RSEM uses file paths
    this function gets both tpms and reads for transcripts'''
    cols = samplenames.copy()
    cols.extend(['transcript_id', 'gene_id'])
    files = {'transcript_reads': transcript_reads, 'transcript_tpm': transcript_tpm}
    for file in files.keys():
        for chunk in pd.read_csv(files[file], header=0, sep="\t", usecols=cols, chunksize=5000):
            if to_db:
                if schema is not None:
                    chunk.to_sql(file, eng, if_exists='append', chunksize=5000, schema=schema)
                else:
                    chunk.to_sql(file, eng, if_exists='append', chunksize=5000)

def get_junction_info(samplenames, junction_reads, eng, to_db=True, schema=None):
    '''get junctions uses paths'''
    cols = samplenames.copy()
    cols.extend(["junction_id", "Description"])
    for chunk in pd.read_csv(junction_reads, header=0, sep="\t", usecols=cols, skiprows=2, chunksize=5000):
        newcols = ["chrom", "start", "end"]
        for col in newcols:
            index = newcols.index(col)
            junctions = chunk.iloc[:, 0]
            chunk[col] = [var.split("_")[index] for var in junctions]
        chunk.columns.values[1]="gene_id"
        if to_db:
            if schema is not None:
                chunk.to_sql("junction_reads", eng, if_exists='append', chunksize=5000, schema=schema)
            else:
                chunk.to_sql("junction_reads", eng, if_exists='append', chunksize=5000)

def get_exon_info(samplenmes, exon_reads, eng, to_db=True, schema=None):
    '''get exon read counts uses paths write to db'''
    cols=samplenmes.copy()
    cols.append("exon_id")
    for chunk in pd.read_parquet(exon_reads,  header=0, sep="\t", usecols=cols, chunksize=5000):
        chunk["gene_id"]=[exon.split("_")[0] for exon in chunk.iloc[:,0]]
        if to_db:
            if schema is not None:
                chunk.to_sql("exon_reads", eng, if_exists='append', chunksize=5000, schema=schema)
            else:
                chunk.to_sql("exon_reads", eng, if_exists='append', chunksize=5000)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='split gtex data into tissue specific databases')
    parser.add_argument('-d', '--datadir', type=str, help='path of the dir that has all the UNCOMPRESSED DATA')
    parser.add_argument('-o', '--outdir', type=str, help='path of the output directory will be created if doesnt exist')
    parser.add_argument('-f','--files', type=str, help='path of the file that contains the paths of files to be processed')
    parser.add_argument('-t', '--threads', type=int, help='number of threads to spawn')
    args=parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # this is ugly

    files=pd.read_csv("files_file.csv", header=None, sep=',', index_col=0, skiprows=1, squeeze=True).to_dict()

    # tissues file is manually generated bc I'm lazy
    tissues=pd.read_csv(files['tissues'], squeeze=True, header=None)
    tissues=tissues.values.tolist()

    def make_db(tissue, files=files, d=args.outdir):
        with open("logs.txt", "a") as log:
            l=tissue+" started\n"
            log.write(l)
            try:
                # I am not creating any indexes because this is temporary storage it will move to postgres
                fname = tissue.replace("(", " ")
                fname = fname.replace(")", " ")
                fname = fname.replace("-", " ")
                fname = " ".join(fname.split())
                fname = fname.replace(" ", "_")
                fname = fname.lower()

                dbname = d + "/" + fname + "_gtex.db"
                conn = "sqlite:///" + dbname
                engine = sql.create_engine(conn, echo=False)
                samples = get_sample_names(tissue, files['samples'], "SMTSD", engine,
                                   filter="RNASEQ", filter_col="SMAFRZE", to_db=True, ret=True)
                samples = samples.tolist()
                get_subject_info(samples, files['subjects'], "SUBJID", engine, to_db=True, ret=False)
                get_gene_info(samples, files['gene_reads'], files['gene_tpm'], engine, to_db=True)
                get_transcript_info(samples, files['transcript_reads'], files['transcript_tpm'], engine,
                                    to_db=True)
                get_junction_info(samples, files['junction'], engine, to_db=True)
                get_exon_info(samples, files['exon_reads'], engine, to_db=True)
            except:
                log.exception(tissue+" there was an error\n")

    with Pool(args.threads) as p:
        p.map(make_db, tissues)

