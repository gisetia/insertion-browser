import pandas as pd


def load_gene_annotations(data_path, assembly='hg38', known=True,
                          coding=True):

    filters = [("known", "==", known), ("coding", "==", coding)]
    refseq = pd.read_parquet(f'{data_path}/refseq/ncbi-genes-{assembly}.pq',
                             filters=filters)

    return refseq


def load_insertions(data_path, screen_name, chrom, start, end,
                    assembly='hg38'):

    filters = [("pos", ">=", start), ("pos", "<=", end)]
    insertions = pd.read_parquet(f'{data_path}/screen-insertions/'
                                 f'{screen_name}/{assembly}/insertions.pq/'
                                 f'{chrom}', filters=filters)

    return insertions
