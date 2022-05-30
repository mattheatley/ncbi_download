import os, sys, subprocess, argparse, datetime
from re import sub
from core import SBATCH, ezSub, CAPTURE

parser = argparse.ArgumentParser(
                                 description='Download taxon library'
                                ,prog=os.path.basename(__file__)
                                ,usage='%(prog)s [options]'
                                ,epilog='see readme for further details.'
                                )

taxons = [
          'bacteria', 'fungi', 'protozoa', 'archaea', 'viral', 'plant', 'invertebrate'
         ]

databases = [
             'genbank', 'refseq'
            ]

current_date = datetime.datetime.now().strftime("%Y-%m-%d")
parser.add_argument('-p', metavar='</path/to/directory>',   type=str, default=f'~/NCBI/{current_date}',                               help='specify path to working directory')
parser.add_argument('-t', metavar='<taxon>',                type=str, default='all',    nargs='+', choices=['all',*taxons, 'custom'], help='specify taxon to download')
parser.add_argument('-d', metavar='<database>',             type=str, default='all',    nargs='+', choices=['all',*databases],        help='specify database source')
parser.add_argument('-q', metavar='<quality>',              type=str, default='good',              choices=['good','best'],           help='select only Complete/Chromosome (good) or Representative (best) genomes')
parser.add_argument('-c', metavar='</path/to/custom_list>', type=str, required=False,                                                 help='specify path to custom list')

modes = parser.add_mutually_exclusive_group(required=True) # run modes
modes.add_argument( '-download', action='store_true', help='download files'     )
modes.add_argument( '-review',   action='store_true', help='review file status' )

PATH, SELECTED_TAXONS, SELECTED_DATABASES, QUALITY, CUSTOM_LIST, downloading, reviewing = vars(parser.parse_args()).values() # define user inputs

if 'custom' in SELECTED_TAXONS:
    if not CUSTOM_LIST: 
        print("ERROR!")
        sys.exit(0)
else:
    if CUSTOM_LIST:
        
        custom_accessions = [ ]

        for line in open(CUSTOM_LIST,'r').readlines():
        
            if line.startswith('#'): continue
        
            accession, *_ = info = line.strip().split('\t')

            *_, accession = accession.split(':')

            custom_accessions.append(accession)

taxons    =                                                taxons    if 'all' in SELECTED_TAXONS else SELECTED_TAXONS
databases = ['custom'] if 'custom' in SELECTED_TAXONS else databases if 'all' in SELECTED_DATABASES else SELECTED_DATABASES

WRK_DIR = f'/{PATH.strip("/")}'
# LOCAL TESTING
WRK_DIR = f'/Users/matt/Desktop{WRK_DIR}'
# LOCAL TESTING
os.makedirs(WRK_DIR, exist_ok=True)

for SOURCE in databases:

    # format path as per rsync_from_ncbi.pl

    RSYNC_PATH = 'https://ftp.ncbi.nlm.nih.gov/genomes'
    ASSEMBLY_SUFFIX='_genomic.fna.gz'

    BASE_URL = f'{RSYNC_PATH}/{SOURCE}'

    SOURCE_DIR = f'{WRK_DIR}/{SOURCE}/{QUALITY}'
    os.makedirs(SOURCE_DIR, exist_ok=True)


    for TAXON in taxons:

        TAXON_DIR = f'{SOURCE_DIR}/{TAXON}'
        GEN_DIR = f'{TAXON_DIR}/all'
        [ os.makedirs(subdir, exist_ok=True) for subdir in [TAXON_DIR, GEN_DIR] ]

        MANIFEST_FILE = f'{TAXON_DIR}/manifest.txt'
        SELECTED_FILE = f'{TAXON_DIR}/selected.txt'
        ASSEMBLY_FILE = CUSTOM_LIST if TAXON == 'custom' else f'{TAXON_DIR}/assembly_summary.txt'

        if downloading:

            # specify assembly summary to download
            if TAXON != 'custom':

                ASSEMBLY_URL = f'{BASE_URL}/{TAXON}/assembly_summary.txt'

                # download assembly file
                #subprocess.run(f'wget -O {ASSEMBLY_FILE} {ASSEMBLY_URL}', shell=True) # linux
                subprocess.run(f'cURL -o {ASSEMBLY_FILE} {ASSEMBLY_URL}', shell=True) # unix
            
            with open(MANIFEST_FILE, 'w') as MANIFEST, open(SELECTED_FILE, 'w') as SELECTED:

                for i, line in enumerate(open(ASSEMBLY_FILE, 'r').readlines()):

                    if i == 1: 
                        
                        line = line.strip('\n')
                        print(line, file=SELECTED)

                        # extract relevant column headers
                        assembly_header = line.strip('# ').split('\t')

                        header_idx = [ 
                            assembly_header.index(name)
                            for name in [
                                 'assembly_accession'
                                ,'taxid'
                                ,'organism_name'
                                ,'assembly_level'
                                ,'refseq_category'
                                ,'ftp_path'
                                ]
                            ]

                    # ignore other comments
                    if line.startswith('#'): continue

                    # extract relevant line info
                    info = line.strip().split('\t')

                    assembly_info = [
                                     info[idx]
                                     for idx in header_idx
                                    ]

                    assembly_accession, taxid, organism_name, assembly_level, refseq_category, ftp_path = assembly_info



                    # specify filter criteria; any listed assembly (custom) or only best quality (taxon) 
                    assembly_level_criteria  = ('') if TAXON == 'custom' or QUALITY == 'best' else ('chromosome', 'complete genome')
                    refseq_category_criteria = ('') if TAXON == 'custom' or QUALITY == 'good' else ('reference', 'representative')

                    if assembly_level.lower().startswith(assembly_level_criteria):

                        if refseq_category.lower().startswith(refseq_category_criteria):

                            # ignore assembly if also part of custom set
                            if TAXON != 'custom' and CUSTOM_LIST:

                                if assembly_accession in custom_accessions: 

                                    print(f'skipping {organism_name}; {assembly_accession} (included in custom list)')
                                    continue


                            # path recorded in assembly file
                            if ftp_path != 'na':

                                print(*info, file=SELECTED)

                                # remove server from path
                                ftp_path = ftp_path.replace(f'{RSYNC_PATH}/','')

                                *_, basename = os.path.split(ftp_path)
                                ftp_file = f'{basename}{ASSEMBLY_SUFFIX}'

                                manifest_path = f'{ftp_path}/{ftp_file}'

                                print(manifest_path, file=MANIFEST)



        scripts = []

        batch_dirs = [ f'{TAXON_DIR}/{description}' for description in ['id','sh','oe'] ]

        id_dir, sh_dir, oe_dir = batch_dirs

        [
         os.makedirs(path, exist_ok=True)
         for path in batch_dirs
        ]

        id_file = f'{id_dir}/{TAXON}.id'


        if reviewing:
            
            print('REVIEWING...')
            
            GEN_FILES = [
                         contents.name
                         for contents in os.scandir(GEN_DIR)
                        ]

        for line in open(MANIFEST_FILE, 'r').readlines():

            URL = line.strip()

            *_, TARGET_FILE = os.path.split(URL)
            
            SAMPLE = TARGET_FILE.replace(ASSEMBLY_SUFFIX,'')
            DESCRIPTION = f'{TAXON}_{SAMPLE}'

            if reviewing:

                if not TARGET_FILE in GEN_FILES:

                    print(f'{TARGET_FILE} missing!')

                    err_file = f'{oe_dir}/{DESCRIPTION}.err'   

                    status = CAPTURE(f'cat {err_file}')

                    print(status)

            # write slurm file
            if downloading:

                sh_file = f'{sh_dir}/{DESCRIPTION}.sh'

                with open(sh_file, 'w') as sh:

                    # specify hpcc directives (slurm)
                    hpcc_directives = SBATCH(
                        job_id=DESCRIPTION
                        ,partition='defq'
                        ,nodes=1
                        ,ntasks=1
                        ,memory='5gb'
                        ,walltime='02:00:00'
                        ,out_err=oe_dir
                        )
                    
                    sh.write(
                        hpcc_directives
                    +'echo STARTED `date`\n'
                    +f'wget --spider {RSYNC_PATH}/{URL}\n'
                    +f'wget -O {GEN_DIR}/{TARGET_FILE} {RSYNC_PATH}/{URL}\n'
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
