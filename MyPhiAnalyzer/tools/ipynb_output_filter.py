#! /usr/bin/python

import sys
from IPython.nbformat.current import read, write

json_in = read(sys.stdin, 'json')

for sheet in json_in.worksheets:
    for cell in sheet.cells:
        if "outputs" in cell:
            cell.outputs = []
        if "prompt_number" in cell:
            cell.prompt_number = ''

write(json_in, sys.stdout, 'json')

# does not work yet on lxplus, because IPython.nbformat.current does not exist.
# if you want to filter output from ipython notebook using this file,
# work on the repository locally on your machine and run the following git commands there:

#git config filter.dropoutput_ipynb.clean /path_to_filter/ipynb_output_filter.py
#git config filter.dropoutput_ipynb.smudge cat
