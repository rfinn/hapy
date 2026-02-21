#!/usr/bin/env python

import os
import glob
# change files that are pointing_XX to pointing-XX

ffiles = glob.glob('pointing_*')
#print(ffiles)
for f in ffiles:
    #print(f)
    new_name = f.replace('pointing_','pointing-')
    os.rename(f,new_name)
