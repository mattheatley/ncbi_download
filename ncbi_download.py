import os, sys, subprocess, argparse
from re import sub
from core import SBATCH, ezSub, CAPTURE

parser = argparse.ArgumentParser(
                                 description='Download taxon library'
                                ,prog=os.path.basename(__file__)
                                ,usage='%(prog)s [options]'
                                ,epilog='see readme for further details.'
                                )

taxons = [
          'bacteria'
         ,'fungi'
         ,'protozoa'
         ,'archaea'
         ,'viral'
         ,'plant'
         ,'invertebrate'
         ]

databases = [
             'genbank'
            ,'refseq'
            ]

parser.add_argument('-p',        metavar='</path/to/working_directory>', type=str, default='~/NCBI',                              help='specify path to working directory')
parser.add_argument('-t',        metavar='<taxon>',                      type=str, default='all',     choices=['all',*taxons],    help='specify taxon to download')
parser.add_argument('-d',        metavar='<database>',                   type=str, default='all',     choices=['all',*databases], help='specify database source')
parser.add_argument('-download', action='store_true',                                                                             help='download files')
parser.add_argument('-review',   action='store_true',                                                                             help='review file status')

PATH, SELECTED, DATABASE, downloading, reviewing = vars(parser.parse_args()).values() # define user inputs

taxons    = taxons    if SELECTED == 'all' else [SELECTED]
databases = databases if DATABASE == 'all' else [DATABASE]

WRK_DIR = f'/{PATH.strip("/")}'

os.makedirs(WRK_DIR, exist_ok=True)

for SOURCE in databases:

    # format path as per rsync_from_ncbi.pl

    RSYNC_PATH = 'https://ftp.ncbi.nlm.nih.gov/genomes/'
    SUFFIX='_genomic.fna.gz'

    BASE_URL = f'{RSYNC_PATH}{SOURCE}'

    SOURCE_DIR = f'{WRK_DIR}/{SOURCE}'

    os.makedirs(SOURCE_DIR, exist_ok=True)

    for TAXON in taxons:

        TAXON_DIR = f'{SOURCE_DIR}/{TAXON}'
        GEN_DIR = f'{TAXON_DIR}/all'

        ASSEMBLY_URL = f'{BASE_URL}/{TAXON}/assembly_summary.txt'

        ASSEMBLY_FILE = f'{TAXON_DIR}/assembly_summary.txt'
        MANIFEST_FILE = f'{TAXON_DIR}/manifest.txt'
                
        [ os.makedirs(subdir, exist_ok=True) for subdir in [TAXON_DIR, GEN_DIR] ]

        if downloading:

            # download assembly file
            subprocess.run(f'wget -O {ASSEMBLY_FILE} {ASSEMBLY_URL}', shell=True) # linux
            #subprocess.run(f'cURL -o {ASSEMBLY_FILE} {ASSEMBLY_URL}', shell=True) # unix
        
            with open(MANIFEST_FILE, 'w') as MANIFEST:

                for line in open(ASSEMBLY_FILE, 'r').readlines():
                    
                    if line.startswith('#'): continue
                    
                    line = line.strip().split('\t')

                    taxid, = line[5:6]
                    assembly_level, = line[11:12]
                    ftp_path, = line[19:20]

                    if assembly_level in ['Chromosome', 'Complete Genome']:

                        # path recorded in assembly file
                        if ftp_path != 'na':

                            # remove server path
                            ftp_path = ftp_path.replace(f'{RSYNC_PATH}','')

                            *_, basename = os.path.split(ftp_path)
                            ftp_file = f'{basename}{SUFFIX}'

                            manifest_path = f'{ftp_path}/{ftp_file}'

                            print(manifest_path, file=MANIFEST)

        scripts = []

        batch_dirs = [ f'{TAXON_DIR}/{description}' for description in ['id','sh','oe','completed'] ]

        id_dir, sh_dir, oe_dir, downloaded_dir = batch_dirs

        [ os.makedirs(path, exist_ok=True) for path in batch_dirs ]

        id_file = f'{id_dir}/{TAXON}.id'


        if reviewing:
            
            print('REVIEWING...')
            GEN_FILES = [ contents.name for contents in os.scandir(GEN_DIR) ]

        for line in open(MANIFEST_FILE, 'r').readlines():

            URL = line.strip()

            *_, TARGET_FILE = os.path.split(URL)
            
            SAMPLE = TARGET_FILE.replace(SUFFIX,'')

            if reviewing:

                if not TARGET_FILE in GEN_FILES:

                    print(f'{TARGET_FILE} missing!')

                    err_file = f'{oe_dir}/{SAMPLE}.err'   

                    status = CAPTURE(f'cat {err_file}')

                    print(status)

            # write slurm file
            if downloading:

                sh_file = f'{sh_dir}/{SAMPLE}.sh'

                with open(sh_file, 'w') as sh:

                    # specify hpcc directives (slurm)
                    hpcc_directives = SBATCH(
                        job_id=SAMPLE
                        ,partition='defq'
                        ,nodes=1
                        ,ntasks=1
                        ,memory='1gb'
                        ,walltime='01:00:00'
                        ,out_err=oe_dir
                        )
                    
                    sh.write(
                        hpcc_directives
                    +'echo STARTED `date`\n'
                    +f'wget --spider {RSYNC_PATH}{URL}\n'
                    +f'wget -O {GEN_DIR}/{TARGET_FILE} {RSYNC_PATH}{URL}\n'
                    #+f'mv {sh_file} {downloaded_dir}/\n'
                    +'echo FINISHED `date`\n'
                    )

                    scripts.append(sh_file)

        # submit tasks
        if downloading:

            print('\nSUBMITTING TASKS...')                    
            with open(id_file, 'a+') as id_output:
                for i,script in enumerate(scripts,1): 
                    ezSub(i=i, check=120, user='mbzmch', limit=500) # maintain tasks below parellel task limit
                    sub_id = CAPTURE(f'sbatch -d singleton {script}') # submit task
                    if os.path.basename(script).startswith('buddy'): continue
                    else: print(sub_id, script, sep='\t', file=id_output, flush=True) # record task job id & shell script
            print('\tALL TASKS SUBMITTED\n')

