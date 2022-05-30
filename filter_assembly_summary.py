
import os, sys, subprocess, datetime
from core import CAPTURE

NCBI_dir = '/Users/matt/Desktop/NCBI'
os.makedirs(NCBI_dir, exist_ok=True)

current_date = datetime.datetime.now().strftime("%Y-%m-%d")

current_subdir = f'{NCBI_dir}/{current_date}'
os.makedirs(current_subdir, exist_ok=True)

# DOWNLOAD NCBI ASSEMBLY INFO
info_urls = [
    'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt',
    'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt'
    ]

info_downloads = []

os.chdir(current_subdir)

for url in info_urls:

    download = os.path.basename(url)
    info_downloads.append(download)

    if not os.path.exists(download):
    
        subprocess.call(f'curl -O {url}', shell=True)

info_prefix = os.path.commonprefix(info_downloads)



# EXTRACT RELEVANT GENOMES

genus_list = [
    'Fusarium',
    'Rhizoctonia',
    'Pythium',
    'Globisporangium',
    'Microdochium',
    'Oculimacula',
    'Parastagonospora',
    'Pyrenophora',
    'Zymoseptoria'
    ]

for genus in genus_list:
    
    genus_subdir = f'{current_subdir}/{genus}'
    os.makedirs(genus_subdir, exist_ok=True)

    outputs = [ 
        f'{genus_subdir}/{genus}_{description}.txt' 
        for description in [
            'all'
            ,'representative'
            ]
        ]

    all_genomes, representative_genomes = outputs

    headers = CAPTURE(f'head -n 2 {info_prefix}* | tail -n 2')

    for output in outputs:

        prep_cmd = f'echo "{headers}" > {output}'
        subprocess.call(prep_cmd, shell=True)

    filter_cmd = f'grep "{genus}" {info_prefix}* | grep -v "virus" >> {all_genomes}'
    subprocess.call(filter_cmd, shell=True)

    filter_cmd = f'grep "{genus}" {info_prefix}* | grep -v "virus" | grep "representative genome" >> {representative_genomes}'
    subprocess.call(filter_cmd, shell=True)

