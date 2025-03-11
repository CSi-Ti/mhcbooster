
import psutil
import subprocess
import tempfile

from pathlib import Path

def get_msfragger_command(param_path, fasta_path, raw_path, num_threads):
    java_exe_path = Path(__file__).parent.parent / 'third_party' / 'jre-17.0.14' / 'bin' / 'java'
    msfragger_exe_path = Path(__file__).parent.parent / 'third_party' / 'MSFragger-4.1' / 'MSFragger-4.1.jar'
    msfragger_split_path = Path(__file__).parent.parent / 'third_party' / 'msfragger_pep_split.py'
    avail_mem = psutil.virtual_memory().available / (1024 ** 3)
    print(f'Avail Memory = {avail_mem:.2f} GB')
    fasta_size = Path(fasta_path).stat().st_size / (1024 ** 2)
    split = int(fasta_size / avail_mem * 8)
    print(f'Splitting fasta to {split} slices...')

    param_data = list(open(param_path).read().splitlines())
    param_data = [line for line in param_data if not line.startswith('database_name')
                  and not line.startswith('num_threads')]
    param_data.insert(0, f'database_name = {fasta_path}')
    param_data.insert(1, f'num_threads = {num_threads}')
    with tempfile.NamedTemporaryFile('w', delete=False) as tmp_param_file:
        tmp_param_file.write('\n'.join(param_data))
        tmp_param_file.flush()

    commands = []
    for raw_file in Path(raw_path).iterdir():
        if raw_file.suffix not in {'.d', '.raw', '.RAW', '.wiff'}:
            continue
        if raw_file.with_suffix('.pin').exists():
            continue
        if split > 1:
            command = (
                f'python {msfragger_split_path} {split} "{java_exe_path} -Xmx{int(avail_mem)}G -jar -Dfile.encoding=UTF-8"'
                f' {msfragger_exe_path} {tmp_param_file.name} {raw_file}')
        else:
            command = (f'{java_exe_path} -Xmx{int(avail_mem)}G -jar -Dfile.encoding=UTF-8 {msfragger_exe_path} '
                       f'{tmp_param_file.name} {raw_file}')
        commands.append(command)
    return commands

def run_msfragger(param_path, fasta_path, raw_path, num_threads):
    commands = get_msfragger_command(param_path, fasta_path, raw_path, num_threads)
    for command in commands:
        subprocess.run(command, shell=True)