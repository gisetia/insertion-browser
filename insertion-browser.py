# %%
# from bokeh.events import DocumentReady
import os
from importlib import reload
from bokeh.layouts import column, row
from bokeh.plotting import output_file, show, curdoc
from bokeh.models import (AutocompleteInput, Div, TextInput,
                          RadioButtonGroup, Dropdown, Select)

from tools.load_data import *
import tools.plotting.insertions as pltins
import tools.plotting.transcripts as plttx

data_path = 'gs://gisetia-insertion-browser/processed_data'
if os.uname().sysname == 'Darwin':
    data_path = 'processed_data'

# %%
menu_margins = (20, 0, 0, 0)
menu_width = 150
txt_out = Div(text='', margin=menu_margins, width=menu_width)

assembly = 'hg38'
screen_name = 'screenA'
gene = 'JAK2'

padd = 200000

txt_out.text = 'Loading gene annotations...'

print('Loading refseq...')

refseq = load_gene_annotations(data_path, assembly=assembly)

gene_pos = refseq.query('name_chrom == @gene')
chrom = gene_pos.chrom.head(1).values[0]
start = gene_pos.txStart.min()
end = gene_pos.txEnd.max()

txt_out.text = 'Loading screen insertions...'
print('Loading insertions...')
insertions = load_insertions(data_path, screen_name, chrom, start-padd,
                             end+padd, assembly=assembly)
txt_out.text = ''

# %%
reload(pltins)
reload(plttx)


def plot_ins(insertions, screen_name, chrom, start, end, refseq,
             screen_type='ip', assembly='hg38'):

    print(f'Plotting for {screen_name} - {assembly} - {chrom}:{start}-{end}')

    # print('Creating plots...')
    ins = pltins.InsertionPlot(insertions, screen_name, assembly, chrom,
                               start, end, screen_type=screen_type,
                               jitter_ins=True, load_padd=padd)
    ins.hover_position()

    select = pltins.InsertionPlot(insertions, screen_name, assembly, chrom,
                                  start, end, screen_type=screen_type,
                                  jitter_ins=True, load_padd=padd)
    select.for_selection(ins.plt)

    transcript = plttx.TranscriptPlot(refseq, assembly, chrom, start, end,
                                      load_padd=padd)
    transcript.link_ins(ins.plt)
    global transcr
    transcr = transcript.transcripts

    plots = column(ins.div_title, transcript.plt, ins.plt, select.plt,
                   name='plots')

    # print(transcr)
    # plots = column(transcript.plt, ins.plt)

    return plots


plots = plot_ins(insertions, screen_name, chrom, start, end, refseq)


# Menus
pos_input = TextInput(title='1-based position', value='',
                      placeholder='eg. chr9:5,450,503-5,470,567',
                      width=menu_width, margin=menu_margins)

screen_opts = ['screenA', 'screenB']
screen_menu = Select(title='Select screen', value='screenA',
                     options=screen_opts, width=menu_width,
                     margin=menu_margins)

assembly_opts = ['hg19', 'hg38']
assembly_menu = RadioButtonGroup(labels=assembly_opts, active=1,
                                 width=menu_width, margin=menu_margins)
gene_opts = sorted(refseq.name_chrom.unique())
gene_menu = AutocompleteInput(title='Gene (refseq symbol)', value=gene,
                              completions=gene_opts, width=menu_width,
                              min_characters=1, case_sensitive=False,
                              margin=menu_margins)


def load_gene(attr, old, new):
    txt_out.text = 'Loading gene...'
    curdoc().add_next_tick_callback(update_gene)


def update_gene():
    
    gene = gene_menu.value
    screen_name = screen_menu.value
    assembly = assembly_opts[assembly_menu.active]

    gene_pos = refseq.query('name_chrom == @gene')
    global chrom
    global start
    global end
    chrom = gene_pos.chrom.head(1).values[0]
    start = gene_pos.txStart.min()
    end = gene_pos.txEnd.max()

    insertions = load_insertions(data_path, screen_name, chrom, start-padd,
                                 end+padd, assembly=assembly)
    plots = plot_ins(insertions, screen_name, chrom, start, end, refseq,
                     assembly=assembly)
    layout.children[1] = plots
    txt_out.text = 'Finished loading gene.'


def load_position(attr, old, new):
    txt_out.text = 'Loading position...'
    curdoc().add_next_tick_callback(update_position)


def update_position():
    
    screen_name = screen_menu.value
    assembly = assembly_opts[assembly_menu.active]
    position = pos_input.value

    # Convert to 0-based left-closed right-open
    global chrom
    global start
    global end
    chrom = f'chr{position.split(":")[0][3:]}'
    start = int(position.split(':')[1].split('-')[0].replace(',', '')) - 1
    end = int(position.split(':')[1].split('-')[1].replace(',', ''))

    insertions = load_insertions(data_path, screen_name, chrom, start-padd,
                                 end+padd, assembly=assembly)
    plots = plot_ins(insertions, screen_name, chrom, start, end, refseq,
                     assembly=assembly)
    layout.children[1] = plots
    txt_out.text = 'Finished loading position.'


def load_assembly(attr, old, new):
    txt_out.text = 'Loading assembly...'
    curdoc().add_next_tick_callback(update_refseq)


def update_refseq():
    
    screen_name = screen_menu.value
    assembly = assembly_opts[assembly_menu.active]

    global refseq
    refseq = load_gene_annotations(data_path, assembly=assembly)
    try:
        plots = plot_ins(insertions, screen_name, chrom, start, end, refseq,
                        assembly=assembly)
        layout.children[1] = plots
        txt_out.text = 'Finished loading assembly.'
    except OSError:
        txt_out.text = 'No data found for these parameters.'


def load_screen(attr, old, new):
    txt_out.text = 'Loading screen...'
    curdoc().add_next_tick_callback(update_screen_ins)

def update_screen_ins():
    screen_name = screen_menu.value
    assembly = assembly_opts[assembly_menu.active]

    try:
        global insertions
        insertions = load_insertions(data_path, screen_name, chrom, start-padd,
                                     end+padd, assembly=assembly)
        plots = plot_ins(insertions, screen_name, chrom, start, end, refseq,
                     assembly=assembly)
        layout.children[1] = plots
        txt_out.text = 'Finished loading screen.'
    except OSError:
        txt_out.text = 'No data found for these parameters.'


screen_menu.on_change('value', load_screen)
assembly_menu.on_change('active', load_assembly)
gene_menu.on_change('value', load_gene)
pos_input.on_change('value', load_position)

menus = column(screen_menu, assembly_menu, gene_menu, pos_input, txt_out)

layout = row(menus, plots)

curdoc().add_root(layout)
# st.bokeh_chart(layout, use_container_width=False)

# title = st.text_input('Select gene', gene, on_change=update_gene)

# st.button("Rerun")
# %%
