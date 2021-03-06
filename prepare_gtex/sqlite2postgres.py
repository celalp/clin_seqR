import argparse
import pandas as pd
import logging as log
import sqlalchemy as sql
import sqlalchemy_utils as utils

def to_pgres(tissue, datadir, user, pwd, dbname, samples):
    fname = tissue.replace("(", " ")
    fname = fname.replace(")", " ")
    fname = fname.replace("-", " ")
    fname = " ".join(fname.split())
    fname = fname.replace(" ", "_")
    sname = fname.lower()


    pgres = 'postgresql://' + user + ":" + pwd + "@localhost/" + dbname
    pgres_engine = sql.create_engine(pgres)
    if not utils.database_exists(pgres_engine.url):
        utils.functions.create_database(pgres_engine.url)
    pgres_engine.execute(sql.schema.CreateSchema(sname))
    pgres_engine.execute(sql.schema.CreateSchema('samples'))

    try:
        samples.to_sql('samples', pgres_engine, schema='samples', if_exists='fail', index=False)
        with pgres_engine.connect() as con:
            con.execute("create index on samples.samples (sampid)")
    except:
        log.exception("samples table error")

    dbname = datadir + "/" + fname + "_gtex.db"
    conn = "sqlite:///" + dbname
    slite_engine = sql.create_engine(conn, echo=False)

    logger = log.getLogger()
    handler = log.StreamHandler()
    logger.addHandler(handler)
    for table in slite_engine.table_names():
        if table == "samples" or table == 'subjects':
            continue
        else:
            try:
                tabdf=pd.read_sql_table(table, slite_engine)
                tabdf.rename(columns=lambda x: x.lower(), inplace=True)
                tabdf.to_sql(table, pgres_engine, schema=sname, if_exists='fail', index=False)
                msg=tissue + " " + table+" export done \n"
                logger.info(msg)
                with pgres_engine.connect() as con:
                    index = "create index on " + sname + "." + table + " (gene_id)"
                    con.execute(index)
                    msg = tissue + " " + table + " indexing done \n"
                    logger.info(msg) #this also is not working
            except Exception:
                msg=tissue + table + "error"
                logger.exception(msg)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='split gtex data into tissue specific databases')
    parser.add_argument('-d', '--datadir', type=str, help='path of the dir that has all the sqlite dbs')
    parser.add_argument('-u', '--username', type=str, help="If using postgres, it's hard coded", default=None)
    parser.add_argument('-p', '--pwd', type=str, help="If using a database type that is not sqlite", default=None)
    parser.add_argument('-n', '--dbname', type=str, help="Name of the posgtres database")
    parser.add_argument('-t', '--tissue', type=str, help="Name of file with tissue names")
    parser.add_argument('--sample', type=str, help="samples file", default=None)
    parser.add_argument('--subject', type=str, help="subject file", default=None)
    args=parser.parse_args()

    tissues = pd.read_csv(args.tissue, squeeze=True, header=None)
    tissues = tissues.values.tolist()

    samples=pd.read_csv(args.sample, header=0, sep="\t")
    samples.rename(columns=lambda x: x.lower(), inplace=True)
    samples = samples.loc[samples.SMAFRZE =="RNASEQ"]
    subjects=[]
    for sample in samples["sampid"]:
        subjid="GTEX-"+sample.split("-")[1]
        subjects.append(subjid)
    samples["subjid"]=subjects

    subjects=pd.read_csv(args.subjects, header=0, sep="\t")
    subjects.rename(columns=lambda x: x.lower(), inplace=True)

    samples=samples.merge(subjects, on='subjid', how='left')

    for tissue in tissues:
        to_pgres(tissue, datadir=args.datadir, user=args.username, pwd=args.pwd, dbname=args.dbname, samples=samples)

