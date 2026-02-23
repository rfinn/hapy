    #####################################
    ## SETUP COMMAND-LINE PARAMETERS
    #####################################
    import argparse    
    parser = argparse.ArgumentParser(description ='Run gui for analyzing Halpha images')

    parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    
    parser.add_argument('--rimage',dest = 'rimage', default=None,help='r-band image')
    parser.add_argument('--haimage',dest = 'haimage', default=None,help='Halpha image')
    parser.add_argument('--csimage',dest = 'csimage', default=None,help='Continuum-subtracted Halpha image')    
    parser.add_argument('--filter',dest = 'filter', default=None,help='filter. options are 4, 8, 12, 16, inthalpha, or intha6657')
    parser.add_argument('--tabledir',dest = 'tabledir', default=None,help='table directory. something like /home/rfinn/research/Virgo/tables-north/v1/')
    parser.add_argument('--psfdir',dest = 'psfdir', default=None,help='set this to the directory containing PSF images')        
    parser.add_argument('--prefix',dest = 'prefix', default='v17p03',help='prefix associated with the coadded image.  Default is v17p03. required when running auto.')
    parser.add_argument('--auto',dest = 'auto', action='store_true',default=False,help='set this to process the images automatically, without the gui')
    
    parser.add_argument('--virgo',dest = 'virgo', action='store_true',default=False,help='set this if running on virgo data.  The virgo filaments catalog will be used as input.')
    parser.add_argument('--uat',dest = 'uat', action='store_true',default=False,help='set this if running on uat halpha groups.  The AGC (210720) will be used as the parent catalog.')     
    parser.add_argument('--draco',dest = 'draco', action='store_true',default=False,help='set this if running on draco.')   
    parser.add_argument('--nebula',dest = 'nebula', action='store_true',default=False,help='set this if running on open nebula virtual machine.  catalog paths will be set accordingly.')
    parser.add_argument('--laptop',dest = 'laptop', action='store_true',default=False,help="custom setting for running on Rose's laptop. catalog paths will be set accordingly.")
    
    parser.add_argument('--obsyear',dest = 'obsyear', default=None,help='year that data were taken.  this finds the right image directory if you are building the image name in pieces..  ')
    parser.add_argument('--pointing',dest = 'pointing', default=None,help='Pointing number that you want to load.  ONLY FOR VIRGO DATA, and only if you are buildling the image name in pieces.')
    
    parser.add_argument('--testing',dest = 'testing', action='store_true',default=False,help='set this if running on open nebula virtual machine')
    parser.add_argument('--onegal',dest = 'onegal', default=None, help='provide galaxy name to run halpha gui just on one galaxy')    
    parser.add_argument('--verbose',dest = 'verbose', action='store_true',default=False,help='set this for extra print statements')    
        
    args = parser.parse_args()
    
    logger = log.get_logger("example1", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)


    sepath = os.getenv('HOME')+'/github/halphagui/astromatic/'

