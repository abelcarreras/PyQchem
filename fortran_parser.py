import os
from graphviz import Digraph

root_dir = '/Users/abel/Programes/qchem/rasman2/'
dir_list = os.listdir(root_dir)
dir_list = [ file for file in dir_list if file.endswith('.F')]
print(dir_list)

dot = Digraph(comment='Test')

sub_names = {}
for file in dir_list:
    with open(root_dir+file, 'r') as f:
        for line in f.readlines():
            if line.upper().find('SUBROUTINE') > 0 and line[0].upper() != 'C':
                # print(line)
                n = line.upper().find('SUBROUTINE')
                sub_line = line[n:].upper().strip()
                sub_name = sub_line.split('(')[0].replace('SUBROUTINE', '').strip()
                # print(sub_names, file)
                sub_names.update({sub_name: file})
                dot.node(file, file)

pairs = []
files = []
for file in dir_list:
    with open(root_dir + file, 'r') as f:
        for line in f.readlines():
            if line.upper().find('CALL') > 0 and line[0].upper() != 'C':
                n = line.upper().find('CALL')
                call_lines = line[n:].upper().strip()
                call_name = call_lines.split('(')[0].replace('CALL', '').strip()
                print(call_name)
                try:
                    if not file+sub_names[call_name] in pairs:
                        print(call_name, file, sub_names[call_name])
                        dot.edge(file, sub_names[call_name], label=call_name)
                        pairs.append(file+sub_names[call_name])
                        files.append(sub_names[call_name])
                        files.append(file)

                except KeyError:
                    pass
                    # dot.edge(file, 'EXTERNAL', label=call_name)


for b in dot.body:
    if not '{}'.format(b).split('"')[1] in files:
        dot.body.remove(b)


#dot.engine = 'fdp'
dot.engine = 'dot'

dot.graph_attr['pad'] = "0.5"
dot.graph_attr['nodesep'] = "0.25"
dot.graph_attr['ranksep'] = "2"

print(dot.source)
dot.render('test.gv', view=True)
