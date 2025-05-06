##@title Analyze your protein

import os
from google.colab import files
import datetime
import re

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

########## input
INPUT = "MSHRKFSAPRHGHLGFLPHKRSHRHRGKVKTWPRDDPSQPVHLTAFLGYKAGMTHTLREVHRPGLKISKREEVEAVTIVETPPLVVVGVVGYVATPRGLRSFKTIFAEHLSDECRRRFYKDWHKSKKKAFTKACKRWRDTDGKKQLQKDFAAMKKYCKVIRVIVHTQMKLLPFRQKKAHIMEIQLNGGTVAEKVAWAQARLEKQVPVHSVFSQSEVIDVIAVTKGRGVKGVTSRWHTKKLPRKTHKGLRKVACIGAWHPARVGCSIARAGQKGYHHRTELNKKIFRIGRGPHMEDGKLVKNNASTSYDVTAKSITPLGGFPHYGEVNNDFVMLKGCIAGTKKRVITLRKSLLVHHSRQAVENIELKFIDTTSKFGHGRFQTAQEKRAFMGPQKKHLEKETPETSGDL" #@param ["RPL3L", "MYC"] {allow-input: true}

#@markdown - Input format: one raw protein sequence; space allowed
#@markdown - Example: copy & paste a multi-line sequence from a FASTA file (without the header)
#@markdown - To run: click `Runtime` -> `Run all` in the menu bar, or click the triangle play/run button on the left

seq = INPUT

if seq == "RPL3L":
  seq = "MSHRKFSAPRHGHLGFLPHKRSHRHRGKVKTWPRDDPSQPVHLTAFLGYKAGMTHTLREVHRPGLKISKREEVEAVTIVETPPLVVVGVVGYVATPRGLRSFKTIFAEHLSDECRRRFYKDWHKSKKKAFTKACKRWRDTDGKKQLQKDFAAMKKYCKVIRVIVHTQMKLLPFRQKKAHIMEIQLNGGTVAEKVAWAQARLEKQVPVHSVFSQSEVIDVIAVTKGRGVKGVTSRWHTKKLPRKTHKGLRKVACIGAWHPARVGCSIARAGQKGYHHRTELNKKIFRIGRGPHMEDGKLVKNNASTSYDVTAKSITPLGGFPHYGEVNNDFVMLKGCIAGTKKRVITLRKSLLVHHSRQAVENIELKFIDTTSKFGHGRFQTAQEKRAFMGPQKKHLEKETPETSGDL"
elif seq == "MYC":
  seq = "MDFFRVVENQQPPATMPLNVSFTNRNYDLDYDSVQPYFYCDEEENFYQQQQQSELQPPAPSEDIWKKFELLPTPPLSPSRRSGLCSPSYVAVTPFSLRGDNDGGGGSFSTADQLEMVTELLGGDMVNQSFICDPDDETFIKNIIIQDCMWSGFSAAAKLVSEKLASYQAARKDSGSPNPARGHSVCSTSSLYLQDLSAAASECIDPSVVFPYPLNDSSSPKSCASQDSSAFSPSSDSLLSSTESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEEEIDVVSVEKRQAPGKRSESGSPSAGGHSKPPHSPLVLKRCHVSTHQHNYAAPPSTRKDYPAAKRVKLDSVRVLRQISNNRKCTSPRSSDTEENVKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILSVQAEEQKLISEEDLLRKRREQLKHKLEQLRNSCA"
else:
  seq = seq.upper().replace(' ','')
  if not all(char in 'ACDEFGHIKLMNPQRSTVWY' for char in seq):
    print('\n'+ bcolors.BOLD +bcolors.FAIL + "WARNING: Your sequence contains letters other than ACDEFGHIKLMNPQRSTVWY!"+bcolors.ENDC)
    L0  = len(seq)
    seq = re.sub('[^ACDEFGHIKLMNPQRSTVWY]+', '', seq)
    L1 = len(seq)
    print(L0-L1,'non-aa letters removed!'+bcolors.ENDC)
    exit()

########## Options

MODEL = "esm1b_t33_650M_UR50S" #@param ["esm1v_t33_650M_UR90S_1", "esm1b_t33_650M_UR50S"]

if os.path.exists("ESMScan-all-mutants.txt"):
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+': Removing files from a previous run')
  !rm ESMScan-* res.zip run.sh

if not os.path.exists("ESM-Scan"):
  print('\n\n'+ bcolors.BOLD +bcolors.OKBLUE + "Installing packages"  +bcolors.ENDC)
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
  !pip install biopython
  !pip install fair-esm
  !git clone https://github.com/xuebingwu/ESM-Scan.git
  !mv /content/ESM-Scan/esm1b_t33_650M_UR50S-contact-regression.pt /content/

# 🔧 Patch for PyTorch 2.6+ unpickling restriction
!sed -i 's/torch.load(str(model_location), map_location="cpu")/torch.load(str(model_location), map_location="cpu", weights_only=False)/' /usr/local/lib/python3.*/dist-packages/esm/pretrained.py

model_location="/content/"+MODEL+".pt"
if not os.path.exists(model_location):
  print('\n\n'+ bcolors.BOLD +bcolors.OKBLUE + "Downloading pre-trained ESM model"  +bcolors.ENDC)
  if MODEL == "esm1b_t33_650M_UR50S":
    !wget https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt
  else:
    !wget https://dl.fbaipublicfiles.com/fair-esm/models/esm1v_t33_650M_UR90S_1.pt

print('\n\n'+ bcolors.BOLD +bcolors.OKBLUE + "Running saturation mutagenesis"  +bcolors.ENDC)

cmd="python /content/ESM-Scan/esmscan.py --model-location "+model_location+" --sequence "+seq

print(cmd)

with open("run.sh",'w') as f:
  f.write(cmd+'\n')

!chmod +x /content/run.sh
!/content/run.sh

print('\n\n'+ bcolors.BOLD +bcolors.OKBLUE + "Downloading results"  +bcolors.ENDC)

if os.path.exists('ESMScan-res-in-matrix.csv'):
  os.system(f'zip res.zip *.pdf *.csv')
  files.download(f"res.zip")
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+': Done! Please see results in res.zip')
else:
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+': No output files generated')
