# %%
import os
from importlib import reload
from bokeh.layouts import column, row
from bokeh.plotting import output_file, show, curdoc
from bokeh.models import (AutocompleteInput, Div, TextInput,
                          RadioButtonGroup, Dropdown, Select)

from tools.load_data import *
import tools.plotting.insertions as pltins
import tools.plotting.transcripts as plttx


# col1, col2 = st.columns([3, 2])

# x = [1, 2, 3, 4, 5]
# y = [6, 7, 2, 4, 5]

# p = figure(
#     title='simple line example',
#     x_axis_label='x',
#     y_axis_label='y', frame_height=200)

# p.line(x, y, legend_label='Trend', line_width=2)

# col1.bokeh_chart(p, use_container_width=True)
# col2.write('This is a sidebar')
# %%
screen_name = 'screenA'
gene = 'EPC2'

refseq = load_gene_annotations('refseq_hg38').query(
    'known == True and coding == True')

gene_pos = refseq.query('name_chrom == @gene')
chrom = gene_pos.chrom.head(1).values[0]
start = gene_pos.txStart.min()
end = gene_pos.txEnd.max()

# %%

# position = 'chr9:5,450,503-5,470,567'
# chrom = f'chr{position.split(":")[0][3:]}'

# # Convert to 0-based left-closed right-open
# start = int(position.split(':')[1].split('-')[0].replace(',', '')) - 1
# end = int(position.split(':')[1].split('-')[1].replace(',', ''))

insertions = load_insertions(screen_name)


# %%
reload(pltins)
reload(plttx)


def plot_ins(insertions, screen_name, chrom, start, end, refseq,
             screen_type='ip', assembly='hg38'):

    print(f'Plotting for {screen_name} - {assembly} - {chrom}:{start}-{end}')

    # print('Creating plots...')
    ins = pltins.InsertionPlot(insertions, screen_name, assembly, chrom,
                               start, end, screen_type=screen_type,
                               jitter_ins=True)
    ins.hover_position()

    select = pltins.InsertionPlot(insertions, screen_name, assembly, chrom,
                                  start, end, screen_type=screen_type)
    select.for_selection(ins.plt)

    transcript = plttx.TranscriptPlot(refseq, assembly, chrom, start, end)
    transcript.link_ins(ins.plt)

    plots = column(ins.div_title, transcript.plt, ins.plt, select.plt,
                   name='plots')

    # plots = column(transcript.plt, ins.plt)

    return plots


plots = plot_ins(insertions, screen_name, chrom, start, end, refseq)

# Menus
menu_margins = (20, 30, 0, 10)
menu_width = 150
pos_input = TextInput(title='1-based position', value='',
                      placeholder='eg. chr9:5,450,503-5,470,567',
                      width=menu_width, margin=menu_margins)

# screen_opts = ['Screen A', 'Screen B']
# screen_menu = AutocompleteInput(title='Screen', value=screen_name,
#                                 completions=screen_opts, width=menu_width,
#                                 min_characters=1, case_sensitive=False,
#                                 margin=menu_margins)

screen_opts = ['screenA', 'screenB']
screen_menu = Select(title='Select screen', value='screenA',
                     options=screen_opts, width=menu_width,
                     margin=menu_margins)


assembly_opts = ['hg19', 'hg38']
assembly_menu = RadioButtonGroup(labels=assembly_opts, active=1,
                                 width=menu_width, margin=menu_margins)
gene_opts = sorted(refseq.name2.unique())
gene_menu = AutocompleteInput(title='Gene (refseq symbol)', value=gene,
                              completions=gene_opts, width=menu_width,
                              min_characters=1, case_sensitive=False,
                              margin=menu_margins)

txt_out = Div(text='', margin=menu_margins, width=menu_width)


def load_gene(attr, old, new):
    txt_out.text = 'Loading gene...'
    curdoc().add_next_tick_callback(update_gene)


def update_gene():
    txt_out.text = 'Finished loading gene.'
    gene = gene_menu.value
    screen_name = screen_menu.value
    assembly = assembly_opts[assembly_menu.active]

    gene_pos = refseq.query('name2 == @gene')
    global chrom
    global start
    global end
    chrom = gene_pos.chrom.head(1).values[0]
    start = gene_pos.txStart.min()
    end = gene_pos.txEnd.max()

    plots = plot_ins(insertions, screen_name, assembly, chrom, start, end,
                     refseq)
    layout.children[1] = plots


menus = column(screen_menu, assembly_menu, gene_menu, pos_input, txt_out)

layout = row(menus, plots)

curdoc().add_root(layout)
# st.bokeh_chart(layout, use_container_width=False)

# title = st.text_input('Select gene', gene, on_change=update_gene)

# st.button("Rerun")
# %%
