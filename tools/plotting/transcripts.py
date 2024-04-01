import pandas as pd
import numpy as np
from math import pi
from typing import Optional
from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, HoverTool, Range1d,
                          NumeralTickFormatter, LabelSet, Arrow, NormalHead,
                          Title, CrosshairTool)

from tools.refseq import collapse_gene_refseq, get_exon_regions


class TranscriptPlot():

    def __init__(self, refseq: pd.DataFrame, assembly: str, chrom: str,
                 start: int, end: int,
                 load_padd: Optional[int] = 100000,
                 dashed_edges: Optional[bool] = True,
                 load_gene: Optional[str] = None,
                 x_axis: Optional[bool] = True) -> None:

        self.refseq = refseq
        self.assembly = assembly
        self.chrom = chrom
        self.start = start
        self.end = end

        self.load_padd = load_padd
        self.load_start = start - self.load_padd
        self.load_end = end + self.load_padd

        self.load_gene = load_gene

        self.transcripts = self.load_transcripts()
        self.exons = self.load_exons()

        if self.load_gene is None:
            plot_height = 300
            plot_title = ''
        else:
            tx_num = len(self.transcripts)
            plot_height = (tx_num)*20
            plot_title = self.load_gene

        plot_height = (len(self.transcripts))*20
        print('plot_height', len(self.transcripts), plot_height)

        # Set plot
        self.ylim = (0.5, len(self.transcripts) + 0.5)
        self.plt = figure(plot_width=1000, frame_height=plot_height,
                          x_range=Range1d(start-(end-start)/60,  # - 1000,
                                          end+(end-start)/60,  # + 1000,
                                          bounds=(self.load_start,
                                                  self.load_end),
                                          min_interval=5
                                          ),
                          y_range=self.ylim,
                          title=plot_title,
                          min_border_left=150, min_border_right=50,
                          #   min_border_top=100,
                          x_axis_location='below',
                          tools='reset, save, xwheel_pan, xwheel_zoom',
                          active_scroll='xwheel_pan')
        self.plt.xaxis.formatter = NumeralTickFormatter(format='0,0')
        self.plt.ygrid.grid_line_color = None
        self.plt.xgrid.grid_line_color = None
        self.plt.yaxis.visible = False
        self.plt.outline_line_color = None

        if not x_axis:
            self.plt.xaxis.visible = False

        # Plot dashed edges at start and end
        if dashed_edges:
            self.plt.line(x=(self.start, self.start), y=self.ylim,
                          color='#E2E2E4', line_dash=[5, 2])
            self.plt.line(x=(self.end, self.end), y=self.ylim, color='#E2E2E4',
                          line_dash=[5, 2])

        # Plot transcript lines
        self.transcripts.apply(self.plot_tx_line, axis=1)
        self.plt.add_tools(HoverTool(tooltips=[('Gene', '@name_chrom'),
                                               ('Transcript', '@name'),
                                               ('Strand', '@strand'),
                                               #    ('txStart', '@start'),
                                               #    ('txEnd', '@end'),
                                               ('Length', '@length{0,0}'),
                                               ('Coding', '@coding'),
                                               ('Known', '@known')],
                                     line_policy='interp',
                                     names=['tx_line']))

        # Plot exons
        self.exons.apply(self.plot_exon, axis=1)

        self.transcripts.apply(self.plot_tx_arrow, axis=1)

        # self.plot_scale_bar()

    def plot_scale_bar(self):
        bar_size = 10**int(np.log10(self.end - self.start) - 1)

        self.plt.rect(x=self.start, width=bar_size, y=self.ylim[1]-1,
                      height=0.2,
                      color='#2d3436')

        lbl_dict = {'size': [f'{int(bar_size/1000)}kb'],
                    'xpos': [self.start], 'ypos': [self.ylim[1]-.8]}
        lbl_source = ColumnDataSource(lbl_dict)
        labels = LabelSet(x='xpos', y='ypos', text='size',
                            x_offset=0, y_offset=0, source=lbl_source,
                            text_font_size='8pt')
        self.plt.add_layout(labels)

    def plot_exon_rect(self, ex):

        # Generate data source
        ex_df = pd.DataFrame()
        ex_df['ypos'] = (self.transcripts.tx_id[self.transcripts['name']
                                                == ex['name']])
        ex_df['width'] = ex.reg_lims[1] - ex.reg_lims[0]
        ex_df['xpos'] = ex.reg_lims[0] + ex_df['width']/2 - 0.5

        # ex_df['color'] = '#fdae61' if ex.reg_type == 'exCds' else '#3288bd'
        # ex_df['height'] = 0.8 if ex.reg_type == 'exCds' else 0.5
        ex_df['color'] = 'red' if ex.reg_type == 'exCds' else 'red'
        ex_df['height'] = 1 if ex.reg_type == 'exCds' else 0.8

        ex_source = ColumnDataSource(ex_df)

        self.plt.rect(x='xpos', y='ypos', width='width', height='height',
                      color='color', source=ex_source)

    def plot_exon(self, ex):
        ex_df = pd.DataFrame()
        y_center = (self.transcripts.tx_id[self.transcripts['name']
                                           == ex['name']])
        y_off = 0.3 if ex.reg_type == 'exCds' else 0.15
        ex_df['ypos'] = [y_center + y_off, y_center - y_off,
                         y_center - y_off, y_center + y_off]
        ex_df['xpos'] = [ex.reg_lims[0] - 0.5, ex.reg_lims[0] - 0.5,
                         ex.reg_lims[1] - 0.5, ex.reg_lims[1] - 0.5]
        # ex_color = '#fdae61' if ex.reg_type == 'exCds' else '#3288bd'

        ex_color = '#e58e26' if ex.reg_type == 'exCds' else '#0c2461'

        ex_source = ColumnDataSource(ex_df)

        self.plt.patch(x='xpos', y='ypos', color=ex_color, source=ex_source)

    def hide_tools(self) -> None:
        self.plt.toolbar_location = None

    def load_exons(self) -> pd.DataFrame:

        exons = (self.transcripts.groupby('name_chrom')
                 .apply(get_exon_regions).reset_index(drop=True))

        return exons

    def load_transcripts(self) -> pd.DataFrame:

        chrom = self.chrom
        chrom_refseq = self.refseq.query('chrom == @chrom')
        coll_refseq = chrom_refseq.groupby(
            'name_chrom').apply(collapse_gene_refseq)

        load_start = self.load_start
        load_end = self.load_end
        if self.load_gene is None:
            load_genes = coll_refseq.query('(txEnd >= @load_start '
                                           '& txEnd <= @load_end) |'
                                           '(txStart >= @load_start '
                                           '& txStart <= @load_end) |'
                                           '(txStart <= @load_start '
                                           '& txStart >= @load_end)').index
            q_refseq = chrom_refseq.query('name_chrom in @load_genes')
        else:
            q_refseq = chrom_refseq.query('name_chrom == @self.load_gene')

        q_refseq = q_refseq.sort_values(by='name_chrom')
        q_refseq['tx_id'] = q_refseq.reset_index().index + 1

        return q_refseq

    def plot_tx_line(self, tx) -> None:

        # Generate data source
        tx_df = pd.DataFrame()
        # tx_df = dict()

        # tx_df['xpos'] = [tx.txStart - 0.5, tx.txEnd - 0.5]
        if tx['strand'] == '-':
            tx_df['xpos'] = [tx.txStart - (self.end-self.start)/120,
                            tx.txEnd + 0.5]
        else:
            tx_df['xpos'] = [tx.txStart - 0.5,
                            tx.txEnd + (self.end-self.start)/120]

        tx_df['ypos'] = [tx.tx_id, tx.tx_id]
        tx_df['name'] = [tx['name'], tx['name']]
        tx_df['name_chrom'] = [tx.name_chrom, tx.name_chrom]
        tx_df['length'] = [tx.txEnd - tx.txStart, tx.txEnd - tx.txStart]
        tx_df['coding'] = [tx.coding, tx.coding]
        tx_df['known'] = [tx.known, tx.known]
        tx_df['strand'] = [tx.strand, tx.strand]
        tx_df['start'] = [tx.txStart, tx.txStart]
        tx_df['end'] = [tx.txEnd, tx.txEnd]

        dash = 'solid' if tx.coding else [3, 1]
        tx_line_color = '#BCBCBF' if tx.known else '#B2D1F0'

        tx_source = ColumnDataSource(tx_df)

        # Plot line
        self.plt.line(x='xpos', y='ypos', source=tx_source, line_dash=dash,
                      line_color=tx_line_color, line_width=1, name='tx_line')

    def plot_tx_arrow(self, tx):

        if tx['strand'] == '-':
            # x_trian = tx['txEnd'] - 0.5
            x_trian = tx['txStart'] - (self.end-self.start)/120
            angle = pi/2
        else:
            # x_trian = tx['txStart'] - 0.5

            x_trian = tx['txEnd'] + (self.end-self.start)/120
            angle = -pi/2

        self.plt.triangle(x=[x_trian], y=[tx['tx_id']], angle=angle, size=6,
                          line_color='#636e72', fill_color='#636e72')

        # Add gene name as label at end of line
        if self.load_gene is None:
            lbl_dict = {'name_chrom': [tx.name_chrom], 'xpos': [tx.txEnd - 0.5],
                        'ypos': [tx.tx_id]}
            lbl_source = ColumnDataSource(lbl_dict)
            labels = LabelSet(x='xpos', y='ypos', text='name_chrom',
                              x_offset=6, y_offset=-6, source=lbl_source,
                              text_font_size='8pt')
            self.plt.add_layout(labels)

    def link_ins(self, ins_plot) -> None:

        self.plt.x_range = ins_plot.x_range
        self.plt.xaxis.visible = False
        self.plt.margin = (0, 0, 40, 0)
        self.plt.add_tools(ins_plot.tools[5])

        # Remove title from insertion plot
        t = Title()
        t.text = ''
        # ins_plot.title = t

    def hover_position(self) -> None:
        # Add dummy line for hover tool with position
        self.plt.line(x=(self.load_start, self.load_end), y=[.1, .1],
                      color='white', line_width=2, name='needshover')

        self.plt.add_tools(HoverTool(tooltips=[(f'Position {self.chrom}',
                                                '$x{0,0}')],
                                     mode='vline',
                                     line_policy='interp',
                                     names=['needshover']),
                           CrosshairTool(dimensions='height',
                                         line_color='#363638',
                                         line_width=1))
