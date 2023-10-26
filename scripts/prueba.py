
import numpy as np

from pyqchem.parsers.parser_rasci import parser_rasci

with open('../../../Desktop/my_programs/g-tensor/g-tensor_final_results/h2o/h2o_def2tzvp_5_5.out', encoding="utf8") as f:
    output = f.read()

output = parser_rasci(output)
data = output['interstate_properties']

# for i in range(1, 20):
#     for j in range(1, 20):
#         if i != j:
#             print(data[(i,j)]['total_soc_mat'])

print('\nsoc_tot (cm-1)')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
for i in range(1, 10):
    line = '{:3}'.format(i)
    for j in range(1, 10):
        try:
            line += '{:18.3f}'.format(np.array(output['interstate_properties'][(i, j)]['total_soc_mat'])[0, 0])
        except KeyError:
            line += '        -         '
    print(line)
