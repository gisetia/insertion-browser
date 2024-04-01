from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.plotting import output_file, show, curdoc
from bokeh.models import AutocompleteInput, Div, TextInput, RadioButtonGroup


x = [1, 2, 3, 4, 5]
y = [6, 7, 2, 4, 5]

p = figure(
    title='simple line example',
    x_axis_label='x',
    y_axis_label='y')

p.line(x, y, legend_label='Trend', line_width=2)


menu = AutocompleteInput(title='menu', value='lala', min_characters=1,
                         completions=['lala', 'lolo', 'lulu'])
layout = row(menu, p)

curdoc().title = 'title blop'
curdoc().add_root(layout)
