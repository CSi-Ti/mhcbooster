import subprocess
from pathlib import Path


JAVA = '/vast/palmer/apps/avx2/software/Java/17.0.4/bin/java'
MSFragger = '/home/rw762/software/MSFragger-4.1/MSFragger-4.1.jar'
MSFragger_SPLIT = '/home/rw762/software/MSFragger-4.1/msfragger_pep_split.py'
PARAMS = '/home/rw762/palmer_scratch/PXD019643/fragger.params'

raw_folder = Path('/home/rw762/palmer_scratch/PXD019643/raw')

for raw_file in raw_folder.glob('*.raw'):
    command = f'python {MSFragger_SPLIT} 4 "{JAVA} -jar -Dfile.encoding=UTF-8 -Xmx30G" {MSFragger} {PARAMS} {raw_file}'
    print(command)
    subprocess.run(command, shell=True)
    break
