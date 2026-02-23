#!/usr/bin/env python

import os

'''
os.chdir('/home/rfinn/research/Virgo/gui-output-2017/')
# 2017 data
pointings = ['1','3','4','6','7','10','12','13','14','16','19','23','27','28']
prefix = 'v17'
obsyear = '2017'
for p in pointings:
    pprefix = prefix+'p%02d'%(int(p))
    os.system('python ~/github/halphagui/testing/halphamain.py --virgo --laptop --pointing '+p+' --auto --prefix '+pprefix+' --obsyear '+obsyear)

os.chdir('/home/rfinn/research/Virgo/gui-output-2018/')
pointings = ['04','19','54','56','57']
prefix = 'v18'
obsyear = '2018'
for p in pointings:
    pprefix = prefix+'p%02d'%(int(p))
    os.system('python ~/github/halphagui/testing/halphamain.py --virgo --laptop --pointing '+p+' --auto --prefix '+pprefix+' --obsyear '+obsyear)

pointings = ['LM-09','lowmass-4','002','0035','003','004','005','006','007','008','012','021','038','09','18','80']
'''

os.chdir('/home/rfinn/research/Virgo/gui-output-2020/')
### SKIPPING LOWMASS FOR NOW B/C THEY THROW A WRENCH IN NAMING CONVENTION
pointings = ['002','0035','003','004','005','006','007','008','012','021','038','09','18','80']
#pointings = ['1','3','4','6','7','10','12','13','14','16','19','23','27','28']
prefix = 'v20'
obsyear = '2020'
for p in pointings:
    try:
        pprefix = prefix+'p%02d'%(int(p))
    except:
        pprefix=prefix+p
    os.system('python ~/github/halphagui/testing/halphamain.py --virgo --laptop --pointing '+p+' --auto --prefix '+pprefix+' --obsyear '+obsyear)
