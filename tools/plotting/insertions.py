# from token import OP
import pandas as pd
import numpy as np
from math import pi
from typing import Optional
from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource, HoverTool, CrosshairTool,
                          RangeTool, Range1d, LinearAxis, NumeralTickFormatter,
                          Div, LabelSet)
from bokeh.palettes import PiYG8
from bokeh.transform import jitter

# from insertools.refseq import collapse_gene_refseq, get_exon_regions

pd.options.mode.chained_assignment = None


class InsertionPlot():

    def __init__(self, insertions: pd.DataFrame, screen_name: str,
                 assembly: str, chrom: str, start: int,
                 end: int, load_padd: Optional[int] = 300000,
                 screen_type: Optional[str] = 'ip',
                 jitter_ins: Optional[bool] = False,
                 plot_height: Optional[int] = None,
                 dashed_edges: Optional[bool] = True,
                 strand: Optional[str] = None) -> None:


        self.load_padd = load_padd
        self.insertions = insertions
        self.screen = screen_name
        self.assembly = assembly
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.screen_type = screen_type
        self.load_start = start - self.load_padd
        self.load_end = end + self.load_padd
        self.jitter = jitter_ins

        # Set plot
        if self.screen_type == 'ip' or self.screen_type == 'pa':
            if self.strand:
                self.ylim = (0.5, 3)
                height = plot_height or 100
            else:
                self.ylim = (0, 5.5)
                height = plot_height or 200

        elif self.screen_type == 'sl':
            self.ylim = (0, 12)
            height = plot_height or 300

        margins = (self.load_end-self.load_start)/60
        self.plt = figure(
            # title=(f'Insertions of screen {screen_name} - {assembly} at '
            #        f'{chrom}:{start + 1:,} - {end:,} ({end-start+1:,} bp)'),
            plot_width=1000, plot_height=height,
            x_range=Range1d(start - margins,
                            end + margins,
                            bounds=(start-self.load_padd, end+self.load_padd),
                            min_interval=5
                            ),
            y_range=self.ylim,
            x_axis_location='below',
            min_border_left=150, min_border_right=50,
            tools='reset, save, xpan, xwheel_zoom',
            active_scroll='xwheel_zoom', active_drag='xpan')

        self.plt.ygrid.grid_line_color = None
        self.plt.xgrid.grid_line_color = None
        self.plt.yaxis.major_tick_line_color = None
        self.plt.yaxis.minor_tick_line_color = None
        self.plt.yaxis.axis_line_color = None
        self.plt.xaxis.visible = False
        self.plt.outline_line_color = None
        self.plt.xaxis.formatter = NumeralTickFormatter(format='0,0')
        self.plt.toolbar.logo = None

        if self.screen_type == 'ip' or self.screen_type == 'pa':

            if self.strand:
                self.plt.yaxis.ticker = [1, 2]
                self.plt.yaxis.major_label_overrides = {1: 'Low',
                                                        2: 'High'}
            else:
                self.plt.yaxis.ticker = [1, 2, 3.5, 4.5]
                self.plt.yaxis.major_label_overrides = {1: 'Low - strand',
                                                        2: 'Low + strand',
                                                        3.5: 'High - strand',
                                                        4.5: 'High + strand'}
        elif self.screen_type == 'sl':
            self.plt.yaxis.ticker = [1, 2, 4, 5, 7, 8, 10, 11]
            self.plt.yaxis.major_label_overrides = {1: 'replicate 4 -',
                                                    2: 'replicate 4 +',
                                                    4: 'replicate 3 -',
                                                    5: 'replicate 3 +',
                                                    7: 'replicate 2 -',
                                                    8: 'replicate 2 +',
                                                    10: 'replicate 1 -',
                                                    11: 'replicate 1 +'}

        self.set_source()

        # Plot insertions
        ins_line_color = '#BCBCBF'
        x_line = (self.load_start, self.load_end)
        # print('aaa', self.plt.yaxis.ticker.ticks)

        # self.plt.line(x=x_line, y=[1, 1], color=ins_line_color, line_width=2,
        #               name='ins_line')
        # self.plt.line(x=x_line, y=[2, 2], color=ins_line_color, line_width=2,
        #               name='ins_line')
        # self.plt.line(x=x_line, y=[4, 4], color=ins_line_color, line_width=2,
        #               name='ins_line')
        # self.plt.line(x=x_line, y=[5, 5], color=ins_line_color, line_width=2,
        #               name='ins_line')

        if jitter_ins:

            for y in self.plt.yaxis.ticker.ticks:
                # self.plt.line(x=x_line, y=[y, y], color=ins_line_color,
                #               line_width=2, name='ins_line')
                self.plt.rect(x=(x_line[1]-x_line[0])/2 + x_line[0],
                              width=x_line[1]-x_line[0],
                              y=y, height=0.8, color=ins_line_color,
                              line_width=0, name='ins_line', alpha=0.2)
            self.plt.circle(x='xpos', y=jitter('ypos', width=0.75,
                                               range=self.plt.y_range),
                            color='color', source=self.source,
                            angle=pi/2, line_width=1, size=1.5,
                            name='insertions_dash')
        else:
            for y in self.plt.yaxis.ticker.ticks:
                self.plt.line(x=x_line, y=[y, y], color=ins_line_color,
                              line_width=2, name='ins_line')
            self.plt.dash(x='xpos', y='ypos', color='color',
                          source=self.source,
                          angle=pi/2, line_width=1, size=15,
                          name='insertions_dash')

        # Plot dashed edges
        if dashed_edges:
            self.plt.line(x=(self.start, self.start), y=self.ylim,
                          color='#E2E2E4', line_dash=[5, 2])
            self.plt.line(x=(self.end, self.end), y=self.ylim, color='#E2E2E4',
                          line_dash=[5, 2])

        # Change title when x range changes
        # self.plt.x_range.on_change('start', self.update_title)
        # self.plt.x_range.on_change('end', self.update_title)

        # Additional title as div
        self.div_title = Div(text='', margin=(0, 0, 0, 90),
                             width=self.plt.plot_width,
                             style={'color': '#494A4A'})

        self.update_div_title('dummy', 'dummy', 'dummy')
        self.plt.x_range.on_change('start', self.update_div_title)
        self.plt.x_range.on_change('end', self.update_div_title)

        self.plot_scale_bar()

    def plot_scale_bar(self):
        bar_size = max(1000, 10**int(np.log10(self.end - self.start) - 1))

        self.plt.rect(x=self.start + bar_size/2, width=bar_size, y=self.ylim[1]-.4,
                      height=0.05,
                      color='#2d3436', name='scale_bar')

        lbl_dict = {'size': [f'{int(bar_size/1000)}kb'],
                    'xpos': [self.start], 'ypos': [self.ylim[1]-.35]}
        lbl_source = ColumnDataSource(lbl_dict)
        labels = LabelSet(x='xpos', y='ypos', text='size',
                            x_offset=0, y_offset=0, source=lbl_source,
                            text_font_size='8pt', name='scale_bar_labels')
        self.plt.add_layout(labels)

    def update_title(self, attr, old, new):

        start = int(self.plt.x_range.start)
        end = int(self.plt.x_range.end)
        self.plt.title.text = (f'Insertions of screen {self.screen} - '
                               f'{self.assembly} at {self.chrom}:{start + 1:,}'
                               f' - {end:,} ({end-start+1:,} bp)')

    def set_source(self) -> ColumnDataSource:
        load_start = self.load_start
        load_end = self.load_end
        chrom = self.chrom
        strand = self.strand

        if self.strand:
            # q_ins = self.insertions.query('chr == @chrom '
            #                               '& pos >= @load_start '
            #                               '& pos <= @load_end '
            #                               '& strand == @strand')
            q_ins = self.insertions.query('strand == @strand')
        else:
            # q_ins = self.insertions.query('chr == @chrom '
            #                               '& pos >= @load_start '
            #                               '& pos <= @load_end')
            q_ins = self.insertions

        if self.screen_type == 'ip' or self.screen_type == 'pa':

            # ins_colors = {'h+': PiYG8[0], 'l+': PiYG8[-1],
            #               'h-': PiYG8[2], 'l-': PiYG8[-3]}

            if self.strand:
                ins_colors = {'h': '#f7b784', 'l': '#3381bd'}
                ins_pos = {'h': 2, 'l': 1}
                q_ins['color'] = q_ins.apply(
                    lambda x: ins_colors[x['chan'][0]],
                    axis=1)
                q_ins['ypos'] = q_ins.apply(
                    lambda x: ins_pos[x['chan'][0]], axis=1)
                q_ins['xpos'] = q_ins.apply(lambda x: x['pos'], axis=1)

            else:
                ins_colors = {'h+': '#f7b784', 'l+': '#3381bd',
                              'h-': '#f7b784', 'l-': '#3381bd'}
                ins_pos = {'h+': 4.5, 'l+': 2,
                           'h-': 3.5, 'l-': 1}
                q_ins['color'] = q_ins.apply(
                    lambda x: ins_colors[x['chan'][0] + x['strand'][0]],
                    axis=1)
                q_ins['ypos'] = q_ins.apply(
                    lambda x: ins_pos[x['chan'][0] + x['strand'][0]], axis=1)
                q_ins['xpos'] = q_ins.apply(lambda x: x['pos'], axis=1)

        elif self.screen_type == 'sl':

            ins_colors = {'+': PiYG8[-1], '-': PiYG8[0]}
            ins_pos = {'4-': 1, '4+': 2, '3-': 4, '3+': 5, '2-': 7, '2+': 8,
                       '1-': 10, '1+': 11}

            q_ins['color'] = q_ins.apply(
                lambda x: ins_colors[x['strand'][0]], axis=1)
            q_ins['ypos'] = q_ins.apply(
                lambda x: ins_pos[str(x['replicate'][0]) + x['strand'][0]],
                axis=1)
            q_ins['xpos'] = q_ins.apply(lambda x: x['pos'], axis=1)

        source_ins = ColumnDataSource(q_ins)

        self.source = source_ins

    def hide_tools(self) -> None:
        self.plt.toolbar_location = None

    def hover_position(self) -> None:
        # Add dummy line for hover tool with position
        self.plt.line(x=(self.load_start, self.load_end), y=[3, 3],
                      color='white', line_width=2, name='needshover')

        self.plt.add_tools(HoverTool(tooltips=[(f'Position {self.chrom}',
                                                '$x{0,0}')],
                                     mode='vline',
                                     line_policy='interp',
                                     names=['needshover']),
                           CrosshairTool(dimensions='height',
                                         line_color='#363638',
                                         line_width=1))

    def for_selection(self, other_plot: figure,
                      padd: Optional[int] = None) -> None:
        # Use self plot for selecting x ranges of other_plot

        padd = padd or self.load_padd

        margins = (self.load_end-self.load_start)/60

        self.plt.tools = []
        self.plt.margin = (40, 0, 0, 0)
        self.plt.plot_height = 150
        self.plt.title = None
        self.plt.xaxis.axis_label = 'Absolute position (0-based)'

        self.plt.select('scale_bar').visible = False
        self.plt.select('scale_bar_labels').visible = False

        self.plt.extra_x_ranges = {'relative_pos':
                                   Range1d(start=-self.load_padd,
                                           end=self.end
                                           - self.start
                                           + self.load_padd)}
        self.plt.add_layout(LinearAxis(x_range_name='relative_pos',
                                       axis_label='Relative position'),
                            'above')

        self.plt.x_range = Range1d(self.start-padd, self.end+padd,
                                   bounds=(self.load_start, self.load_end))
        self.plt.xaxis.formatter = NumeralTickFormatter(format='0,0')
        self.plt.yaxis.visible = False

        if self.jitter:
            ins_glyph = self.plt.select(name='insertions_dash')
            ins_glyph.glyph.size = 0.5

            ins_glyph = self.plt.select(name='ins_line')
            ins_glyph.glyph.line_width = 0         

        else:
            ins_glyph = self.plt.select(name='insertions_dash')
            ins_glyph.glyph.size = 5

            ins_glyph = self.plt.select(name='ins_line')
            ins_glyph.glyph.line_width = 1

        range_tool = RangeTool(x_range=other_plot.x_range)
        range_tool.overlay.fill_color = "navy"
        range_tool.overlay.fill_alpha = 0.1
        self.plt.add_tools(range_tool)
        self.plt.toolbar.active_multi = range_tool
        self.plt.toolbar_location = None

    def update_div_title(self, attr, old, new):

        # Add title of figure
        start = int(self.plt.x_range.start)
        end = int(self.plt.x_range.end)

        self.div_title.text = (f'<b>Insertions of screen {self.screen} - '
                               f'{self.assembly} at {self.chrom}:{start + 1:,}'
                               f'-{end:,} ({end - start:,} bp)</b>')
