import pandas as pd
from itertools import groupby
from operator import itemgetter
from typing import Optional


# def read_refseq(data_dir: str, assembly: str,
#                 coding_only: Optional[bool] = True,
#                 known_only: Optional[bool] = True) -> pd.DataFrame:
#     '''Returns data from refseq file as dataframe. Keeps only standard
#     chromosomes. With coding_only or known_only, non-coding and predicted
#     entries can be included or excluded. (Excluded by default).

#     nooo
#     If 'name_chrom' is True, adds column 'name_chr' with name2_chrom.
#     eg. C1orf141_chr1.
#     '''

#     print(f'Reading refseq file for {assembly}.')

#     filename = f'{data_dir}/ncbi-genes-{assembly}.parquet.snappy'
#     refseq = pd.read_parquet(filename, engine='pyarrow')

#     if coding_only:
#         refseq = refseq[refseq['coding']]

#     if known_only:
#         refseq = refseq[refseq['known']]

#     return refseq


def collapse_gene_refseq(gene_refseq, outer_cds=True, mode='first-last'):
    # Collapse multiple gene transcripts into a single one. Mode first-last
    # takes first txStart and last txEnd,

    # print(gene_refseq.name2.iloc[0])

    if mode == 'first-last':

        if outer_cds:
            cdsStart = gene_refseq.cdsStart.min()
            cdsEnd = gene_refseq.cdsEnd.max()
        else:
            cdsStart = gene_refseq.cdsStart.max()
            cdsEnd = gene_refseq.cdsEnd.min()

        gene_pos = pd.Series({'gene': gene_refseq.name2.iloc[0],
                              'chrom': gene_refseq.chrom.iloc[0],
                              'txStart': gene_refseq.txStart.min(),
                              'cdsStart': cdsStart,
                              'txEnd': gene_refseq.txEnd.max(),
                              'cdsEnd': cdsEnd,
                              'strand': gene_refseq.strand.iloc[0]})

    elif mode == 'longest-cds':
        exons = get_exon_regions(gene_refseq)

        # Get transcripts with the longest cds
        tx_name_cds = get_longest_cds(exons)

        # Dataframe with all transcripts with the longest cds
        longest_cds_tx = gene_refseq.query('name in @tx_name_cds')

        # Choose transcripts with the most upstream cds start
        if longest_cds_tx.iloc[0].strand == '+':
            tx_name = longest_cds_tx['name'].loc[longest_cds_tx['cdsStart']
                                                 == longest_cds_tx['cdsStart']
                                                 .min()].values
            first_start = gene_refseq['name'].loc[gene_refseq['cdsStart']
                                                  == gene_refseq['cdsStart']
                                                  .min()].values
        else:
            tx_name = longest_cds_tx['name'].loc[longest_cds_tx['cdsEnd']
                                                 == longest_cds_tx['cdsEnd']
                                                 .max()].values
            first_start = gene_refseq['name'].loc[gene_refseq['cdsEnd']
                                                  == gene_refseq['cdsEnd']
                                                  .max()].values

        # print(tx_name_cds, tx_name)
        # print(set(tx_name_cds).intersection(set(tx_name)))
        if len(set(tx_name_cds).intersection(set(first_start))) == 0:
            print(f'{gene_refseq.name2.iloc[0]}: different longest cds '
                  'and first cds start', tx_name_cds, first_start)

        # if len(tx_name) > 1:
        #     print(f'{gene_refseq.name2.iloc[0]}: {len(tx_name)} '
        #           f'tx - {tx_name}')

        cols = ['name2', 'chrom', 'txStart', 'cdsStart', 'txEnd', 'cdsEnd',
                'strand', 'name']
        gene_pos = (gene_refseq.query('name in @tx_name[0]')[cols]
                    .rename(columns={'name2': 'gene'}).set_index('gene'))

    return gene_pos


def get_exon_length(tx_exons: pd.DataFrame) -> pd.Series:
    length = tx_exons['reg_lims'].apply(
        lambda x: x[1] - x[0])
    exon_length = length.sum()

    return exon_length


def get_longest_exon(exons: pd.DataFrame) -> str:
    # Returns name of transcript(s) with longest exon
    lengths = exons.groupby('name').apply(get_exon_length)
    # print(lengths)
    tx_name = lengths.index[lengths == lengths.max()].values

    return tx_name


def get_cds_length(tx_exons: pd.DataFrame) -> pd.Series:
    # Returns transcript names as index and their cds length
    tx_cds_exons = tx_exons.query('reg_type == "exCds"')
    tx_cds_exons['length'] = tx_cds_exons['reg_lims'].apply(
        lambda x: x[1] - x[0])
    cds_length = tx_cds_exons['length'].sum()

    return cds_length


def get_longest_cds(exons: pd.DataFrame) -> str:
    # Returns name of transcript(s) with longest cds
    cds_lengths = exons.groupby('name').apply(get_cds_length)
    # print(cds_lengths)
    tx_name = cds_lengths.index[cds_lengths == cds_lengths.max()].values

    return tx_name


# def collapse_gene_refseq_old(gene_refseq):

#     gene_pos = pd.Series({'gene': gene_refseq.name2.iloc[0],
#                           'chrom': gene_refseq.chrom.iloc[0],
#                           'txStart': gene_refseq.txStart.min(),
#                           'cdsStart': gene_refseq.cdsStart.min(),
#                           'txEnd': gene_refseq.txEnd.max(),
#                           'cdsEnd': gene_refseq.cdsEnd.max(),
#                           'strand': gene_refseq.strand.iloc[0]})

#     return gene_pos


def contig_list_lims(lst):
    '''Divide a list into contiguous sublists and return the lower and
    upper limits of such sublists. (Right open???)
    '''
    lims = []
    for k, g in groupby(enumerate(lst), lambda x: x[0] - x[1]):
        x = list(map(itemgetter(1), g))
        lims.append([x[0], x[-1] + 1])
    return lims


def get_exon_regions(gene_pos: pd.DataFrame) -> pd.DataFrame:
    '''Returns dataframe with start and end positions of all gene exons and
    whether they are cds or utr. Eg.:
    name2	name	        exon_id	tx_id	reg_type	reg_lims
    JAK2	NM_001322195.1	2	    1	    exCds	    [5021987, 5022212]
    JAK2	NM_001322195.1	3	    1	    exCds	    [5029782, 5029905]
    JAK2	NM_004972.3	    126	    6	    exUtr	    [5021962, 5021986]
    '''

    exon = gene_pos.copy(deep=True)
    exon = exon.reset_index(drop=True)
    exon['tx_id'] = exon.index + 1

    cols = ['exonStarts', 'exonEnds']
    exon[cols] = exon[cols].applymap(lambda x: x.split(','))

    exon['cdsRange'] = exon.apply(lambda x: range(x.cdsStart, x.cdsEnd),
                                  axis=1)
    exon['exonRange'] = (exon
                         .apply(lambda t:
                                [range(int(x), int(y)) for i, x
                                 in enumerate(t['exonStarts']) for j, y
                                 in enumerate(t['exonEnds']) if i == j
                                 and x != '' and y != ''], axis=1))

    exon = exon.explode('exonRange')
    exon = exon.reset_index(drop=True)
    exon['exon_id'] = exon.index + 1
    exon = exon.drop(columns=['#bin', 'exonStarts', 'exonEnds', 'score',
                              'cdsStartStat', 'cdsEndStat', 'exonFrames',
                              'cdsStart', 'cdsEnd'])

    exon['exCds_lst'] = exon.apply(lambda t: [x for x in t.exonRange
                                              if x in t.cdsRange], axis=1)
    exon['exUtr_lst'] = exon.apply(lambda t: [x for x in t.exonRange
                                              if x not in t.cdsRange], axis=1)

    exon['exCds'] = exon.apply(lambda t: contig_list_lims(t.exCds_lst), axis=1)
    exon['exUtr'] = exon.apply(lambda t: contig_list_lims(t.exUtr_lst), axis=1)

    exon_regions = exon.melt(id_vars=['name2', 'name', 'exon_id', 'tx_id'],
                             value_vars=['exCds', 'exUtr'],
                             var_name='reg_type', value_name='reg_lims')
    exon_regions = exon_regions.explode('reg_lims').dropna()

    return exon_regions


def get_refseq_overlaps(gene_overlaps):
    overlap = pd.Series()
    gene = gene_overlaps.name
    overlap['gene'] = gene
    gene_overlaps['strands_dir'] = (gene_overlaps
                                    .apply(lambda x: 'same' if x.strand1 ==
                                           x.strand2 else 'opposite', axis=1))

    overlap['overlapping_genes'] = ', '.join(gene_overlaps['gene2'])
    overlap['overlapping_strands'] = ', '.join(
        gene_overlaps['strands_dir'])
    overlap['overlapping_lengths'] = ', '.join(gene_overlaps['overlap']
                                               .apply(lambda x: str(x)))
    return overlap
