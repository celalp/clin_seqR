import pandas as pd
import os
import sqlalchemy_utils as utils
import sqlalchemy as sql
import argparse as arg
from collections import OrderedDict as od
import logging as log

# this is a little manual process but you can wrap this with other default values in the future
# this for the vcf by default this creates a sqlite database and only imports the variants and samples table
# we are operating under the assumption that at most you have a trio and may be other siblings but not the
# whole village and even if you do they should be considered other samples to begin with


# attr takes a file with single column with the following values:
# family_id, individual_id,paternal_id,
# Maternal_ID, Sex, Phenotype, Project name and any other keys

def make_ped(samplename, project_name, tissue):
    filename = samplename + ".ped"
    ped_dat = od([('#family_id', samplename),
                  ('individual_id', samplename),
                  ('paterntal_id', -9),
                  ('maternal_id', -9),
                  ('sex', -9),
                  ("phenotype", -9),
                  ('ethnicity', -9),
                  ('project_name', project_name),
                  ('tissue', tissue)])
    ped_df = pd.DataFrame(ped_dat, index=[0])
    ped_df.to_csv(filename, index=False, sep="\t")
    return (filename, ped_df)

# this is now process gemini
#def process_vcf(vcf, ped, samplename, python2, vcf2db, filepath="."):
#    dbpath=filepath+"/"+samplename+".db"
#    cmd = " ".join([python2, vcf2db, vcf, ped, dbpath])
#    os.system(cmd)
#    return dbpath

# alot of  hard coded things in here need to talk to sergey
def process_genes(file, samplename):
    genes=pd.read_csv(file, header=0, sep="\t")
    genes.columns.values=["gene_id", "counts_counts", "gene_name"]
    genes["sample"]=samplename
    cols=["gene_id", "gene_name", "genes_counts", "samplename"]
    genes=genes[cols]
    return(genes)

def process_transcripts(file, samplename):
    transcripts=pd.read_csv(file, header=0, sep="\t")
    transcripts["sample"]=samplename
    return transcripts

def process_junctions(file, samplename):
    junctions=pd.read_csv(file, header=None, sep="\t")
    junctions.columns.values=["chr", "start", "end", "count", "gene_name"] #this prob needs to be matched with gene_id
    junctions['sample']=samplename
    cols=['chr', 'start', 'end', 'gene_name', 'count', 'sample']
    junctions=junctions[cols]
    return junctions


if __name__=="__main__":
    parser = arg.ArgumentParser(description='make an new database or add to exisitng samples database')
    parser.add_argument('--gemini', type=str, help="where the geminidb is")
    parser.add_argument('--genes', type=str, help="genes_file")
    parser.add_argument('--transcripts', type=str, help="transcripts file")
    parser.add_argument('--junctions', type=str, help="junctions file")
    parser.add_argument('-s', '--samplename', type=str)
    parser.add_argument('-p', '--projectname', type=str)
    parser.add_argument('-t', '--tissue', type=str)
    parser.add_argument('--user', type=str, help="username", default=None)
    parser.add_argument('--pwd', type=str, help="password")
    parser.add_argument('-n', '--dbname', type=str, help="Name of the posgtres database", default="samples")
    args = parser.parse_args()

    #create the database
    pgres = 'postgresql://' + user + ":" + pwd + "@localhost/" + dbname
    pgres_engine = sql.create_engine(pgres)
    if not utils.database_exists(pgres_engine.url):
        utils.functions.create_database(pgres_engine.url)


    #put logger here
    ped=make_ped(samplename=args.samplename, project_name=args.projectname, tissue=args.tissue)
    ped[1].to_sql("samples", con=pgres_engine, if_exists='append')

    conn = "sqlite:///" + args.gemini
    slite_engine = sql.create_engine(conn, echo=False)

    #from one db to another in chunks
    for chunk in pd.read_sql_table('variants', slite_engine, chunksize=50000):
        chunk['project_name']=ped[1]['project_name']
        cols_to_keep=[col for col in chunk.columns.values if "gt" not in col or "vep" not in col]
        chunk.to_sql("variants", pgres_engine, if_exists='append', chunksize=50000)

    with pgres_engine.connect() as con:
    # TODO figure out which columns to index
        index = "create index on variants "+index_cols
        con.execute(index)


