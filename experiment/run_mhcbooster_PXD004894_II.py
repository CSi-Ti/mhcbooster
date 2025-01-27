
import sys
from pathlib import Path
from mhcnames import normalize_allele_name

sys.path.append('/home/rw762/workspace/mhc-booster')
from mhcvalidator.validator import MhcValidator

allele_map_I = {
    'MM5' :	    ['A*01', 'A*25', 'B*08', 'B*18'],
    'MM8' :	    ['A*0101', 'A*0301', 'B*0702', 'B*0801', 'C*0701', 'C*0702'],
    'MM12' :	['A*0101', 'B*0801', 'C*0701'],
    'MM15' :	['A*0301', 'A*6801', 'B*2705', 'B*3503', 'C*0202', 'C*0401'],
    'MM16' :	['A*0101', 'A*2402', 'B*0702', 'B*0801', 'C*0701', 'C*0702']
}

allele_map_II = { }

allele_map = allele_map_II

avail_alleles = [line.split(' ')[0] for line in open('/mnt/d/workspace/mhc-validator-2/third_party/netMHCpan-4.1/Linux_x86_64/data/MHC_pseudo.dat')]

pin_file_dir = '/mnt/d/workspace/mhc-validator-2/data/PXD004894/HLA-II'
mzml_folder = Path('/mnt/d/workspace/mhc-validator-2/data/PXD004894/HLA-II')
output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/PXD004894/HLA-II/test')
output_folder.mkdir(parents=True, exist_ok=True)

pin_files = list(Path(pin_file_dir).rglob('*.pin'))

rt_all_models = ['Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic', 'Prosit_2019_irt', 'Prosit_2024_irt_cit']
ms2_all_models = ['Prosit_2019_intensity', 'Prosit_2024_intensity_cit', 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_CID', 'Prosit_2020_intensity_HCD',
                  'ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_Immuno_HCD', 'ms2pip_timsTOF2024']
ccs_top_models = ['AlphaPeptDeep_ccs_generic']

for i in range(len(pin_files)):
    print(f'TracerID: {i}')
    pin = pin_files[i]
    file_name = pin.stem
    if 'edited' in file_name:
        continue
    alleles = []
    for donor in allele_map.keys():
        if donor in file_name:
            for allele in allele_map[donor]:
                allele = normalize_allele_name(allele).replace('*', '')
                if allele in avail_alleles:
                    alleles.append(allele)
                else:
                    allele = allele.split(':')[0]
                    for j in range(len(avail_alleles)):
                        if avail_alleles[j].startswith(allele):
                            alleles.append(avail_alleles[j])
                            break
            break

    use_mhc_scores = True
    if len(alleles) == 0 or 'CLL' in file_name:
        use_mhc_scores = False
        print('Warning: no allele specified for {}'.format(file_name))

    print(file_name)

    validator = MhcValidator(max_threads=20)  # Open a MHCvalidator instance, a new one has to be opened for each .pin file
    validator.set_mhc_params(alleles=alleles, mhc_class='II', min_pep_len=9)  # Load the alleles you specified above
    validator.load_data(pin, filetype='pin')  # Load the pin file
    validator.run(sequence_encoding=True, mhcflurry=use_mhc_scores, bigmhc=use_mhc_scores,
                  netmhcpan=use_mhc_scores, mixmhc2pred=use_mhc_scores,
                  autort=False, deeplc=False,
                  im2deep=False,
                  peptdeep=False,
                  # koina_predictors=rt_all_models + ms2_all_models + ccs_top_models,
                  koina_predictors=['Prosit_2019_irt', 'AlphaPeptDeep_rt_generic', 'ms2pip_Immuno_HCD', 'ms2pip_HCD2021'],
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}')
