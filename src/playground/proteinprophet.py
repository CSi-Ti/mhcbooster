from pyteomics import protxml

data = protxml.read('/mnt/d/test/interact.prot_normal_JY_new.xml')
data = list(data)
sequence_map = {}
max_proteins = 1
for entry in data:
    if len(entry['protein']) > max_proteins:
        max_proteins = len(entry['protein'])
    protein = entry['protein'][0]['protein_name']
    for peptide in entry['protein'][0]['peptide']:
        if peptide['peptide_sequence'] not in sequence_map.keys():
            sequence_map[peptide['peptide_sequence']] = {}
        if protein in sequence_map[peptide['peptide_sequence']].keys():
            sequence_map[peptide['peptide_sequence']][protein] = max(sequence_map[peptide['peptide_sequence']][protein], entry['probability'])
        else:
            sequence_map[peptide['peptide_sequence']][protein] = entry['probability']
print('debug')