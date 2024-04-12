import json
from io import BytesIO, StringIO
import requests
import pandas as pd
# import streamlit as st

session = requests.Session()

# @st.cache_data()


def request_gdrive_data(file_id, session=session):
    url = f'https://drive.google.com/uc?export=download&id={file_id}'
    return requests.get(url)


def load_gene_annotations__(annot_name):

    with open('file_ids/gene_annotations.json', 'r') as f:
        gene_annotations = json.load(f)

    file_id = gene_annotations[annot_name]
    resp = request_gdrive_data(file_id)

    refseq = pd.read_parquet(BytesIO(resp.content))

    return refseq


def load_insertions__(screen_name, assembly='hg38'):

    with open('file_ids/ins_file_ids.json', 'r') as f:
        ins_file_ids = json.load(f)

    channels = ['high', 'low']
    insertions_list = list()
    for channel in channels:

        file_id = ins_file_ids[screen_name][assembly][channel]
        resp = request_gdrive_data(file_id)

        insertions = pd.read_parquet(BytesIO(resp.content))
        insertions['chan'] = channel
        insertions_list.append(insertions)

        insertions = pd.concat(insertions_list)

    return insertions


def load_gene_annotations(data_path, assembly='hg38',
                          known=True, coding=True):

    filters = [("known", "==", known), ("coding", "==", coding)]
    refseq = pd.read_parquet(f'{data_path}/refseq/ncbi-genes-{assembly}.pq',
                             filters=filters)

    return refseq


def load_insertions(screen_name, chrom, start, end, assembly='hg38',):

    filters = [("pos", ">=", start), ("pos", "<=", end)]
    insertions = pd.read_parquet(f'{data_path}/screen-insertons/{screen_name}/'
                                 f'{assembly}/insertions.pq/{chrom}',
                                 filters=filters)

    return insertions
