from pathlib import Path
import requests
import zipfile
import io
import tempfile
from tqdm import tqdm

path = Path(__file__)

print(path)
print(path.parent.parent.parent)
autort_url = 'https://github.com/bzhanglab/AutoRT/archive/refs/heads/master.zip'
autort_root = Path(__file__).parent.parent.parent / 'third_party' / 'AutoRT-master'
third_party_root = autort_root.parent

print('Start to download AutoRT from GitHub. It will take several minutes...')
response = requests.get(autort_url, stream=True)
if response.status_code == 200:
    destination_path = third_party_root / 'autort.zip'
    with open(destination_path, 'wb') as file:
        # Initialize tqdm progress bar with total size and chunk size
        with tqdm(unit='B', unit_scale=True) as pbar:
            # Download the file in chunks
            for chunk in response.iter_content(chunk_size=1024):
                # Write the chunk to the file
                file.write(chunk)
                # Update the progress bar with the size of the chunk
                pbar.update(len(chunk))

        print("\nDownload completed!")
        zip_file = zipfile.ZipFile(destination_path)
        zip_file.extractall(third_party_root)
        print(f'Successfully downloaded AutoRT from GitHub and extracted to {third_party_root}')
else:
    print(f'Failed to download AutoRT from GitHub. Status code: {response.status_code}')