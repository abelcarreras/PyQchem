
import numpy as np
import sys

from pyqchem.parsers.parser_rasci import parser_rasci

with open('../../../Desktop/my_programs/g-tensor/g-tensor_final_results/cucl4_2-_def2tzvp_11_6.out', encoding="utf8") as f:
    output = f.read()

parse_data = parser_rasci(output)
data = parse_data['interstate_properties']

# for i in range(1, len(data)):
#     for j in range(i+1, len(data)):
print(data)
# print(data[(1,2)]['total_soc_mat'])
print()
