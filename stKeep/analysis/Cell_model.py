# -*- coding: utf-8 -*-
"""

@author: chunman zuo
@modified_by: Pedro Castillo
"""
import anndata
import numpy
import torch
import os
import random
import time
import pandas as pd 
import numpy as np
import pyreadr


import torch
import random
import sys
from __main__ import *

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import stKeep

parser  =  stKeep.parameter_setting()
args    =  parser.parse_args()

## random seed 
numpy.random.seed( args.seed )
random.seed( args.seed )
torch.manual_seed( args.seed )
torch.cuda.manual_seed( args.seed )

start = time.time()
# args.outPath         = args.inputPath + 'stKeep/'
args.outPath         = outPath
args.use_cuda        = args.use_cuda and torch.cuda.is_available()

args.spotGene        = args.outPath + args.spotGene
args.spotGroup       = args.outPath + args.spotGroup
args.spotLatent      = args.outPath + args.spotLatent
args.visualFeature   = args.outPath + args.visualFeature
args.spatialLocation = args.outPath + args.spatialLocation
args.pos_pair        = args.outPath + args.pos_pair 

args.patience        = 90 #30
args.tau             = 0.3 #0.05 #0.3
args.feat_drop       = 0.2
args.sample_rate     = [30,1] #[18,9] #[30, 1] [Spot, gene]
args.attn_drop       = 0.1
args.lr              = 0.02 
args.lam             = 0.1

stKeep.Train_cell_model( args )

duration = time.time() - start
print('Finish training, total time is: ' + str(duration) + 's' )