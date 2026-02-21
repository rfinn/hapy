#!/usr/bin/env python

"""
NOTE: 
this is superceded by ~/github/mucho-galfit/mksbatchfilemask.py

GOAL: 
run mask1galaxy.py on ALL galaxies!

USAGE:

python ~/github/halphagui/maskallgalaxies.py


this will move to /mnt/astrophysics/muchogalfit-output and start running

"""
import multiprocessing as mp
import os
import subprocess


def get_mask(dir):
    """ this is the worker function to call mask1galaxy.py  """
    homedir = os.getenv("HOME")
    cmd = f"python {homedir}/github/halphagui/mask1galaxy.py {dir}"
    os.system(cmd)
    #try:
    #    subprocess.call(cmd)
    #except:
    #    print(f"WARNING: trouble running mask for {dir}")


topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
try:
    os.chdir(topdir)
except FileNotFoundError: # assuming that we are running on virgo vms
    topdir = '/mnt/astrophysics/muchogalfit-output/'
    os.chdir(topdir)

# read in Dirs.txt to get list of directories
dlist = []
dfile = open('Dirs.txt','r')
for line in dfile:
    dlist.append(line.rstrip())

# set up multiprocessing
my_pool = mp.Pool(mp.cpu_count())
myresults = [my_pool.apply_async(get_mask(dir)) for dir in dlist]
    
my_pool.close()
my_pool.join()
