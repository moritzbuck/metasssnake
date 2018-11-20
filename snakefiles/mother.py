import os
from glob import glob
from os.path import join as pjoin
import Bio
from tqdm import tqdm

configfile : "params.json"
workdir : config['home']

include : "read_processing.py"
