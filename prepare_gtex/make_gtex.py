import argparse
import pandas as pd
import os
from sqlalchemy import create_engine
import HTSeq as ht
from multiprocessing import Pool




# exon reads, gene reads, transcript reads(rsem)
# gene tpm, transcript tpm, junctions
# eqtl
# skipping normalized expression


# I need to look into indexes and foreign keys for better db integration but that is not that necessary initially
# I need to get gene_id in all expression tables and I also need subject id in samples table, that will be done ad hoc
# in the samples table with the usual string processing

# TODO refactor reading and writing to SQL
# TODO multithread for later use


def get_sample_names(tiss, sampleinfo, tissue_col, eng, filter, filter_col=None, to_db=True, ret=False):
    '''Using the path of the sample info file get the names of the samples'''
    all_samples = pd.read_csv(sampleinfo, header=0, sep="\t")
    tissue_samples=all_samples.loc[all_samples[tissue_col] == tiss]
    if filter_col is not None:
        tissue_samples=tissue_samples.loc[tissue_samples[filter_col] == filter]
    if to_db:
        tissue_samples.to_sql("samples", eng, if_exists='fail', chunksize=5000)
    if ret:
        return(tissue_samples.iloc[:,0])

def get_subject_info(samplenames, subjectinfo, column, eng, to_db=True, ret=False):
    '''Using the path of the subject info file get the subject properties'''
    all_subjects = pd.read_csv(subjectinfo, header=0, sep="\t")
    subjects=list()
    for sample in samplenames:
        subject = sample[0:9] #hard coded for v7 but it's unlikely to change
        subjects.append(subject)
    subj_info=all_subjects[all_subjects[column].isin(subjects)]
    if to_db:
        subj_info.to_sql("subjects", eng, if_exists='fail', chunksize=5000)
    if ret:
        return(subj_info.iloc[:,0])


def get_gene_info(samplenames, gene_reads, gene_tpm, eng, to_db=True):
    '''get gene info, tpm and read count from rnaseqc, uses file paths
    this function gets both reads and tpm for genes, I am not bothering with medians'''
    cols=samplenames.copy()
    cols.extend(["Name", "Description"])
    files={'gene_reads':gene_reads, 'gene_tpm':gene_tpm}
    for file in files.keys():
        for chunk in pd.read_csv(files[file], header=0, sep="\t", skiprows=2, usecols=cols, chunksize=5000):
            chunk.columns.values[[0,1]] = ["gene_id", "gene_name"]
            if to_db:
                chunk.to_sql(file, eng, if_exists='append', chunksize=5000)


def get_transcript_info(samplenames, transcript_reads, transcript_tpm, eng, to_db=True):
    '''get transcript info tpm and read expected count from RSEM uses file paths
    this function gets both tpms and reads for transcripts'''
    cols = samplenames.copy()
    cols.extend(['transcript_id', 'gene_id'])
    files = {'transcript_reads': transcript_reads, 'transcript_tpm': transcript_tpm}
    for file in files.keys():
        for chunk in pd.read_csv(files[file], header=0, sep="\t", usecols=cols, chunksize=5000):
            if to_db:
                chunk.to_sql(file, eng, if_exists='append', chunksize=5000)

def get_junction_info(samplenames, junction_reads, eng, to_db=True):
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
            chunk.to_sql("junctions", eng, if_exists='append', chunksize=5000)

def get_exon_info(samplenmes, exon_reads, eng, to_db=True):
    '''get exon read counts uses paths write to db'''
    cols=samplenmes.copy()
    cols.append("exon_id")
    for chunk in pd.read_csv(exon_reads, header=0, sep="\t", usecols=cols, chunksize=5000):
        chunk["gene_id"]=[exon.split("_")[0] for exon in chunk.iloc[:,0]]
        if to_db:
            chunk.to_sql("exon_reads", eng, if_exists='append', chunksize=5000)

def get_eqtl(dir, tissue, eng, to_db=True):
    '''this is just an import with gene names as keys'''
    loc=dir+"/"+tissue+".v7.signif_variant_gene_pairs.txt"
    eqtl=pd.read_csv(loc, header=0, sep="\t")
    newcols=["chrom", "location", "reference", "variant"]
    for col in newcols:
        index=newcols.index(col)
        variants=eqtl.iloc[:,0]
        eqtl[col]=[var.split("_")[index] for var in variants]
    if to_db:
        eqtl.to_sql("eqtl", eng, if_exists='fail', chunksize=5000)

# this is too slow would need to be re-factored the obj cant be more than 1gb
def get_annot(gtf, path=None, eng=None, to_db=False, write=True, ret=False):
    """this is to parse the gft file into a more
    reasonable table to be included in the db, This prob is not
    the best way to do it since it will be shared among
    all dbs"""
    exons=[]
    file=ht.GFF_Reader(gtf)
    for feature in file:
        if feature.type == 'exon':
            exon_details=feature.attr
            exon_details['start']=feature.iv.start
            exon_details['end'] = feature.iv.end
            exon_details['strand'] = feature.iv.strand
            exon_details=pd.DataFrame(exon_details, index=[0])
            exon_details=exon_details[["gene_id", "gene_name", "transcript_id", "exon_id",
                                       "exon_number", "start", "end", "strand"]]
            if write and path is not None:
                with open(path , 'a') as file:
                    exon_details.to_csv(file, header=False, index=False)
            else:
                exons.append(exon_details)
    if len(exons) > 0:
        exon_df=pd.concat(exons)
        if ret:
            return exon_df
        if to_db and eng is not None:
            exon_df.to_sql("annotations", eng, if_exists='fail', chunksize=5000)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='split gtex data into tissue specific databases')
    parser.add_argument('-d', '--datadir', type=str, help='path of the dir that has all the UNCOMPRESSED DATA')
    parser.add_argument('-o', '--outdir', type=str, help='path of the output directory will be created if doesnt exist')
    parser.add_argument('-f','--files', type=str, help='path of the file that contains the paths of files to be processed')
    parser.add_argument('-t', '--threads', type=int, help='number of threads to spawn')
    args=parser.parse_args()

    print(args)
    os.makedirs(args.outdir, exist_ok=True)

    # this is ugly

    files=pd.read_csv("files_file.csv", header=None, sep=',', index_col=0, skiprows=1, squeeze=True).to_dict()

    # tissues file is manually generated bc I'm lazy
    tissues=pd.read_csv(files['tissues'], squeeze=True, header=None)
    tissues=tissues.values.tolist()

    # make annoation once keep in memory and add to all dbs later
    #annot=get_annot(annotation_file, ret=True, write=False)

    # wrap it around
    def make_db(tissue, files=files, d=args.outdir):
        with open("logs.txt", "a") as log:
            l=tissue+" started"
            log.write(l)
            try:
                fname = tissue.replace("(", " ")
                fname = fname.replace(")", " ")
                fname = fname.replace("-", " ")
                fname = " ".join(fname.split())
                fname = fname.replace(" ", "_")

                dbname = d + "/" + fname + "_gtex.db"
                conn = "sqlite:///" + dbname
                engine = create_engine(conn, echo=False)

        # this part needs refactoring

        # get samples and subjects
                try:
                    samples = get_sample_names(tissue, files['samples'], "SMTSD", engine,
                                   filter="RNASEQ", filter_col="SMAFRZE", to_db=True, ret=True)
                except:
                    log.write(tissue+" error in samples")

                try:
                    samples = samples.tolist()
                    get_subject_info(samples, files['subjects'], "SUBJID", engine, to_db=True, ret=False)
                except:
                    log.write(tissue+" error in subjects")

        # expression(gene, transcript)
                try:
                    get_gene_info(samples, files['gene_reads'], files['gene_tpm'], engine, to_db=True)
                except:
                    log.write(tissue+" error in genes")

                try:
                    get_transcript_info(samples, files['transcript_reads'], files['transcript_tpm'], engine, to_db=True)
                except:
                    log.write(tissue+ "error in transcripts")
        # exon reads, junctions

                try:
                    get_junction_info(samples, files['junction'], engine)
                except:
                    log.write(tissue+" error in junctions")

                try:
                    get_exon_info(samples, files['exon_reads'], engine)
                except:
                    log.write(tissue + " error in exons")

                try:
                    get_eqtl(files['eqtl'], fname, engine)
                except:
                    log.write(tissue + " error in eqtl")

                log.write(tissue + "done")
            except:
                log.write(tissue+" there was an error")


    with Pool(args.threads) as p:
        p.map(make_db, tissues)

