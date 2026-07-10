import urllib.request
import gzip
import os

os.makedirs('resources', exist_ok=True)
pfams = [
    'PF01338', 'PF03006', 'PF07968', 'PF00765', 'PF00196', 'PF08660', # Screening/QS
    'PF00109', 'PF02801', 'PF00668', 'PF00501', 'PF05147', 'PF04738', # PKS/NRPS/RiPPs
    'PF04183', 'PF06283', 'PF06316', 'PF00582', 'PF03936',            # Siderophores/Phenazines/Terpenes
    'PF00589', 'PF00872', 'PF00239'                                   # MGEs
]
combined = ''

for p in pfams:
    print(f'Downloading {p}...')
    url = f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{p}?annotation=hmm'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            raw = response.read()
            try: 
                combined += gzip.decompress(raw).decode('utf-8') + '\n'
            except: 
                combined += raw.decode('utf-8') + '\n'
    except Exception as e: 
        print(f'Failed {p}: {e}')

with open('resources/essential_bgc.hmm', 'w', encoding='utf-8') as f:
    f.write(combined)

print('Done downloading HMMs!')
