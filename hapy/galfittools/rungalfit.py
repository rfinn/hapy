#!/usr/bin/env python

'''


REQUIRED MODULES:
pyds9


'''

#import pyds9
import os
from astropy.io import fits

#import rungalfit

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, List, Tuple, Dict
from astropy.io import fits


@dataclass(frozen=True)
class FitParam:
    """ A fitted value with an uncertainty """
    value: float
    error: float = 0.0

    def as_tuple(self) -> Tuple[float, float]:
        return (self.value, self.error)

@dataclass(frozen=True)
class ComponentResult:
    """ Results for a single-component model (e.g. sersic) """
    xc: FitParam
    yc: FitParam
    mag: FitParam
    re: FitParam
    n: FitParam
    ar: FitParam
    pa: FitParam
    numerical_error_flag: int = 0

@dataclass(frozen=True)
class AsymmetryResult:
    """ Fourier mode 1 asymmetry, if enabled """
    f1: FitParam
    f1pa: FitParam

@dataclass(frozen=True)
class GalfitResult:
    """ Full parsed GALFIT results """
    components: List[ComponentResult]
    sky: Optional[FitParam] = None
    chi2nu: Optional[float] = None
    error: Optional[str] = None
    asymmetry: Optional[AsymmetryResult] = None

    # helpers
    @property
    def ncomp(self) -> int:
        retun len(self.components)

    def component(self, i: int = 0) -> ComponentResult:
        return self.components[i]


def _parse_fitparam_from_header_value(raw) -> Tuple[FitParam, int]:
    """
    Parse a GALFIT header value into FitParam.
    Returns (FitParam, numerical_error_flag_increment).
    """
    numerical_flag = 0
    s = str(raw).strip()

    # Fixed parameters are often written like [12.34]
    if "[" in s and "]" in s:
        s = s.replace("[", "").replace("]", "").strip()
        try:
            return FitParam(float(s), 0.0), 0
        except ValueError:
            # fall through
            pass

    # Typical case: "12.34 +/- 0.56"
    if "+/-" in s:
        left, right = [t.strip() for t in s.split("+/-", 1)]

        # Numerical issues flagged with "*"
        if "*" in left or "*" in right:
            numerical_flag = 1
            left = left.replace("*", "")
            right = right.replace("*", "")

        return FitParam(float(left), float(right)), numerical_flag

    # Sometimes CHI2NU is a raw float string, and ERROR can be a string.
    # We'll let callers special-case those.
    try:
        return FitParam(float(s), 0.0), 0
    except ValueError:
        return FitParam(float("nan"), 0.0), 0



def parse_galfit_results_dc(
    galfit_outimage: str,
    ncomp: int = 1,
    asymflag: bool = False,
) -> GalfitResult:
    """
    Parse GALFIT output FITS into a GalfitResult dataclass.

    Assumes GALFIT writes fitted params into extension 2 header with keys like:
      1_XC, 1_YC, 1_MAG, 1_RE, 1_N, 1_AR, 1_PA, (n+1)_SKY, CHI2NU, ERROR
    and if asymflag:
      1_F1, 1_F1PA (often in component 1)
    """
    hdr = fits.getheader(galfit_outimage, 2)

    components: List[ComponentResult] = []
    asymmetry: Optional[AsymmetryResult] = None

    for ci in range(1, ncomp + 1):
        numerical_error_flag = 0

        def g(key: str) -> FitParam:
            nonlocal numerical_error_flag
            fp, inc = _parse_fitparam_from_header_value(hdr[key])
            numerical_error_flag = max(numerical_error_flag, inc)
            return fp

        comp = ComponentResult(
            xc=g(f"{ci}_XC"),
            yc=g(f"{ci}_YC"),
            mag=g(f"{ci}_MAG"),
            re=g(f"{ci}_RE"),
            n=g(f"{ci}_N"),
            ar=g(f"{ci}_AR"),
            pa=g(f"{ci}_PA"),
            numerical_error_flag=numerical_error_flag,
        )
        components.append(comp)

    # Sky is usually component ncomp+1 in your legacy scheme, e.g. 2_SKY for 1comp, 3_SKY for 2comp
    sky_key = f"{ncomp + 1}_SKY"
    sky: Optional[FitParam] = None
    if sky_key in hdr:
        sky, _ = _parse_fitparam_from_header_value(hdr[sky_key])

    # Asymmetry Fourier mode is usually stored on component 1 keys 1_F1 and 1_F1PA
    if asymflag:
        if "1_F1" in hdr and "1_F1PA" in hdr:
            f1, _ = _parse_fitparam_from_header_value(hdr["1_F1"])
            f1pa, _ = _parse_fitparam_from_header_value(hdr["1_F1PA"])
            asymmetry = AsymmetryResult(f1=f1, f1pa=f1pa)

    # CHI2NU often appears as a float string
    chi2nu = None
    if "CHI2NU" in hdr:
        try:
            chi2nu = float(str(hdr["CHI2NU"]).split()[0])
        except Exception:
            chi2nu = None

    # ERROR can be stringy; keep raw
    err = None
    if "ERROR" in hdr:
        err = str(hdr["ERROR"]).strip()

    return GalfitResult(
        components=components,
        sky=sky,
        chi2nu=chi2nu,
        error=err,
        asymmetry=asymmetry,
    )



def parse_galfit_results(galfit_outimage, asymflag=0, ncomp=1, return_keywords=False):
    """
    Backwards-compatible wrapper around parse_galfit_results_dc.

    Returns:
      - fit_parameters: list where each component contributes 7 FitParam tuples,
        then NUMERICAL_ERROR_FLAG (int), then SKY FitParam tuple, maybe asym terms,
        ERROR (string), CHI2NU (float)
      - header_keywords: parallel list of labels (optional)
    """
    res = parse_galfit_results_dc(galfit_outimage, ncomp=ncomp, asymflag=bool(asymflag))

    fit_parameters = []
    header_keywords = []

    # components
    fields = ["XC", "YC", "MAG", "RE", "N", "AR", "PA"]
    for i, c in enumerate(res.components, start=1):
        for f in fields:
            fp = getattr(c, f.lower() if f != "AR" else "ar")
            fit_parameters.append(fp.as_tuple())
            header_keywords.append(f"{i}_{f}")
        fit_parameters.append(c.numerical_error_flag)
        header_keywords.append("NUMERICAL_ERROR_FLAG")

    # sky
    if res.sky is not None:
        fit_parameters.append(res.sky.as_tuple())
        header_keywords.append(f"{ncomp+1}_SKY")

    # asymmetry (if requested)
    if bool(asymflag) and res.asymmetry is not None:
        fit_parameters.append(res.asymmetry.f1.as_tuple())
        header_keywords.append("1_F1")
        fit_parameters.append(res.asymmetry.f1pa.as_tuple())
        header_keywords.append("1_F1PA")

    # ERROR and CHI2NU
    fit_parameters.append(res.error if res.error is not None else "")
    header_keywords.append("ERROR")
    fit_parameters.append(res.chi2nu if res.chi2nu is not None else float("nan"))
    header_keywords.append("CHI2NU")

    if return_keywords:
        return fit_parameters, header_keywords
    return fit_parameters




class galfit:
    def __init__(self,galname=None,image=None,sigma_image=None,psf_image=None,psf_oversampling=None,mask_image=None,xminfit=None,yminfit=None,xmaxfit=None,ymaxfit=None,convolution_size=None,magzp=None,pscale=None,convflag=1,constraintflag=1,fitallflag=False,ncomp=1,asym=False):
        self.galname=galname
        self.image=image

        self.sigma_image=sigma_image
        self.psf_image=psf_image
        self.psf_oversampling=psf_oversampling
        self.mask_image=mask_image
        self.xminfit=xminfit
        self.yminfit=yminfit
        self.xmaxfit=xmaxfit
        self.ymaxfit=ymaxfit
        self.convolution_size=convolution_size
        self.magzp=magzp
        self.pscale=pscale
        self.convflag=convflag
        self.constraintflag=constraintflag
        self.fitallflag=fitallflag
        self.ncomp=ncomp
        self.asymmetry=asym
        if self.sigma_image == None:
            self.sigma_image = 'none'

        #print(xminfit,xmaxfit,yminfit,ymaxfit,convolution_size)
        #print(self.xminfit,self.xmaxfit,self.yminfit,self.ymaxfit,self.convolution_size)
        #print('psf_image = ',psf_image)
        #print('self.fitall = ',self.fitallflag)
        #print('***%%%%%%%%%%%%%%%%%')

        

    def create_output_names(self):
        if self.asymmetry:
            output_image=str(self.galname)+'-'+ str(self.ncomp) +'Comp-galfit-out-asym.fits'
        else:
            output_image=str(self.galname)+'-'+ str(self.ncomp) +'Comp-galfit-out.fits'

        self.output_image=output_image
        # create galfit input file
        self.galfile=str(self.galname)+'-galfit.input.'+str(self.ncomp)+'Comp'


    def open_galfit_input(self):
        self.galfit_input=open(self.galfile,'w')


    def write_image_params(self):#,input_image,output_image,sigma_image,psf_image,psf_oversampling,mask_image,xminfit,xmaxfit,yminfit,ymaxfit,convolution_size,magzp,pscale,convflag=1,constraintflag=1,fitallflag=0):
        self.galfit_input.write('# IMAGE PARAMETERS\n')
        self.galfit_input.write('A) '+self.image+'              # Input data image (FITS file)\n')
        self.galfit_input.write('B) '+self.output_image+'       # Name for the output image\n')
        self.galfit_input.write('C) %s                # Sigma image name (made from data if blank or "none") \n'%(self.sigma_image))
        if self.convflag:
            self.galfit_input.write('D) '+self.psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
            self.galfit_input.write('E) %i                   # PSF oversampling factor relative to data\n'%(self.psf_oversampling))
        if self.fitallflag:
            self.galfit_input.write('F)            # Pixel mask (ASCII file or FITS file with non-0 values)\n')
        else:
            self.galfit_input.write('F) '+self.mask_image+'           # Pixel mask (ASCII file or FITS file with non-0 values)\n')

        self.galfit_input.write('H) '+str(int(round(self.xminfit)))+' '+str(int(round(self.xmaxfit)))+' '+str(int(round(self.yminfit)))+' '+str(int(round(self.ymaxfit)))+'     # Image region to fit (xmin xmax ymin ymax)\n')
        if self.convflag:
            self.galfit_input.write('I) '+str(int(round(self.convolution_size)))+' '+str(int(round(self.convolution_size)))+'             # Size of convolution box (x y)\n')
        self.galfit_input.write('J) %5.2f              # Magnitude photometric zeropoint \n'%(self.magzp))
        self.galfit_input.write('K) %6.5f   %6.5f         # Plate scale (dx dy)  [arcsec/pix]\n'%(self.pscale,self.pscale))
        self.galfit_input.write('O) regular                # Display type (regular, curses, both)\n')
        self.galfit_input.write('P) 0                   # Create output image only? (1=yes; 0=optimize) \n')
        self.galfit_input.write('S) 0                   # Modify/create objects interactively?\n')


    def set_sersic_params(self,xobj=None,yobj=None,mag=None,rad=None,nsersic=None,BA=None,PA=None,fitmag=1,fitcenter=1,fitrad=1,fitBA=1,fitPA=1,fitn=1,first_time=0):
        self.xobj=xobj
        self.yobj=yobj
        self.mag=mag
        self.rad=rad
        self.nsersic=nsersic
        self.BA=BA
        self.PA=PA
        self.fitmag=fitmag
        self.fitn=fitn
        self.fitcenter=fitcenter
        self.fitrad=fitrad
        self.fitBA=fitBA
        self.fitPA=fitPA
        print('inside rungalfit, fitBA = ',self.fitBA)

        if first_time:
            self.xobj0=xobj
            self.yobj0=yobj
            self.mag0=mag
            self.rad0=rad
            self.nsersic0=nsersic
            self.BA0=BA
            self.PA0=PA
            self.fitmag0=fitmag
            self.fitn0=fitn
            self.fitcenter0=fitcenter
            self.fitrad0=fitrad
            self.fitBA0=fitBA
            self.fitPA0=fitPA
            self.asymmetry0=self.asymmetry
            
    def set_sersic_params_comp2(self,xobj=None,yobj=None,mag=None,rad=None,nsersic=None,BA=None,PA=None,fitmag=1,fitcenter=1,fitrad=1,fitBA=1,fitPA=1,fitn=1,first_time=0):
        self.xobj2=xobj
        self.yobj2=yobj
        self.mag2=mag
        self.rad2=rad
        self.nsersic2=nsersic
        self.BA2=BA
        self.PA2=PA
        self.fitmag2=fitmag
        self.fitn2=fitn
        self.fitcenter2=fitcenter
        self.fitrad2=fitrad
        self.fitBA2=fitBA
        self.fitPA2=fitPA
        if first_time:
            self.xobj02=xobj
            self.yobj02=yobj
            self.mag02=mag
            self.rad02=rad
            self.nsersic02=nsersic
            self.BA02=BA
            self.PA02=PA
            self.fitmag02=fitmag
            self.fitn02=fitn
            self.fitcenter02=fitcenter
            self.fitrad02=fitrad
            self.fitBA02=fitBA
            self.fitPA02=fitPA
            self.asymmetry02=self.asymmetry

    def reset_sersic_params(self):
        self.xobj=self.xobj0
        self.yobj=self.yobj0
        self.mag=self.mag0
        self.rad=self.rad0
        self.nsersic=self.nsersic0
        self.BA=self.BA0
        self.PA=self.PA0
        self.fitmag=self.fitmag0
        self.fitn=self.fitn0
        self.fitcenter=self.fitcenter0
        self.fitrad=self.fitrad0
        self.fitBA=self.fitBA0
        self.fitPA=self.fitPA0
        self.asymmetry=self.asymmetry0
        
    def set_sky(self,sky):
        self.sky=sky

    def write_sersic(self,objnumber,profile, nsersic=None):
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) %s             # Object type \n'%(profile))
        self.galfit_input.write(' 1) %8.1f  %8.1f %i %i  # position x, y        [pixel] \n'%(self.xobj,self.yobj,int(self.fitcenter),int(self.fitcenter)))
        self.galfit_input.write(' 3) %5.2f      %i       # total magnitude     \n'%(self.mag,self.fitmag))
        self.galfit_input.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(self.rad,self.fitrad))
        print('sersic n, fitsersicn = ',self.nsersic,self.fitn)
        if nsersic == None:
            self.galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(self.nsersic,int(self.fitn)))
        else:
            self.galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic,int(self.fitn)))
        print('BA, fitBA = ',self.BA,self.fitBA)
        self.galfit_input.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(self.BA,int(self.fitBA)))
        self.galfit_input.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(self.PA,int(self.fitPA)))
        if self.asymmetry:
            self.galfit_input.write('F1) 0.0001 0.00   1  1     # azim. Fourier mode 1, amplitude & phase angle \n')
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")

    def write_sersic_BD(self):
        
        ########################
        # write disk profile
        ########################    
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: 1 \n')
        self.galfit_input.write(' 0) sersic             # Object type \n')
        self.galfit_input.write(' 1) %8.1f  %8.1f %i %i  # position x, y        [pixel] \n'%(self.xobj,self.yobj,int(self.fitcenter),int(self.fitcenter)))
        self.galfit_input.write(' 3) %5.2f      %i       # total magnitude     \n'%(self.mag,self.fitmag))
        self.galfit_input.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(self.rad,self.fitrad))
        self.galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(self.nsersic,0))
        self.galfit_input.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(self.BA,int(self.fitBA)))
        self.galfit_input.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(self.PA,int(self.fitPA)))
        if self.asymmetry:
            self.galfit_input.write('F1) 0.0001 0.00   1  1     # azim. Fourier mode 1, amplitude & phase angle \n')
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")

        ########################
        # write bulge profile
        ########################
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: 2 \n')
        self.galfit_input.write(' 0) sersic             # Object type \n')
        self.galfit_input.write(' 1) %8.1f  %8.1f %i %i  # position x, y        [pixel] \n'%(self.xobj2,self.yobj2,int(self.fitcenter),int(self.fitcenter)))
        self.galfit_input.write(' 3) %5.2f      %i       # total magnitude     \n'%(self.mag2,self.fitmag))
        self.galfit_input.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(self.rad2,self.fitrad))
        self.galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(self.nsersic2,1)) # allow n to vary
        self.galfit_input.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(self.BA2,int(self.fitBA)))
        self.galfit_input.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(self.PA2,int(self.fitPA)))
        if self.asymmetry:
            self.galfit_input.write('F1) 0.0001 0.00   1  1     # azim. Fourier mode 1, amplitude & phase angle \n')
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
        

    def write_sky(self,objnumber):    
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) sky             # Object type \n')
        self.galfit_input.write(' 1) %8.1f   1  # sky background at center of fitting region [ADUs] \n'%(self.sky))
        self.galfit_input.write(' 2) 0      0       # dsky/dx (sky gradient in x)    \n')
        self.galfit_input.write(' 3) 0      0       # dsky/dy (sky gradient in y) \n')
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
    
        

    def write_input_file(self):
        self.create_output_names()
        self.open_galfit_input()
        print('in rungalfit.run_galfit, self.psf_image = ',self.psf_image)
        
        self.write_image_params()
                
        if (self.ncomp == 1):
            self.write_sersic(1,'sersic')
            self.write_sky(2)
            
        elif (self.ncomp == 2):
            self.write_sersic_BD()
            self.write_sky(3)
            
        if (self.fitallflag):
            print('%%%%%%%%%%%%%% HEY %%%%%%%%%%%%%')
            print('I think fitall is true, just sayin...')
            self.fitall()
        self.close_input_file()
        
    def run_galfit(self,displayflag=False):
        self.write_input_file()
        #print 'self.fitall = ',self.fitall
        s = 'galfit '+self.galfile
        print('run the following: ',s)


        errno=os.system(s)
        self.galfit_flag=1

        image_id=str(self.galname)+'-'
        self.galfit_log=image_id+str(self.ncomp)+'Comp-fit.log'
        s='cp fit.log '+self.galfit_log
        os.system(s)
        self.galfit_out=image_id+str(self.ncomp)+'Comp'+'-galfit.01'
        s='mv galfit.01 '+self.galfit_out
        try:
            os.rename('galfit.01',self.galfit_out)
        except:
            print("appears like galfit did not complete")
            #galflag[j]=0
            self.galfit_flag=0
            return
        if displayflag:
            self.display_results()


    def fitall(self,mindistance=8):
        os.system('cp '+homedir+'research/LocalClusters/sextractor/default.param .')
        os.system('cp '+homedir+'research/LocalClusters/sextractor/default.nnw .')
        s='sex '+self.image+'[1] -c '+homedir+'research/LocalClusters/sextractor/default.sex.24um.galfitsource -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE '+self.sigma_image+' -CATALOG_NAME '+self.galname+'test.cat -CATALOG_TYPE ASCII_HEAD'
        os.system(s)
        # read in SE table to get x,y for sources
        #fname=self.galname+'test.fits'
        fname=self.galname+'test.cat'
        print('FITALL CATALOG NAME = ',fname)
        objnumber=2
        profile='sersic'
        try:
            se=atpy.Table(fname,type='ascii')
            print('found ',len(se.X_IMAGE),' sources on the field of ',self.galname)
            nearbyobjflag=sqrt((se.X_IMAGE-self.xobj)**2+(se.Y_IMAGE-self.yobj)**2) > mindistance
            for k in range(len(se.X_IMAGE)):
                if nearbyobjflag[k]:
                    objnumer=objnumber+1
                    self.add_simple_sersic_object(objnumber,profile,se.X_IMAGE[k],se.Y_IMAGE[k],se.MAG_BEST[k],se.FLUX_RADIUS[k,0],2,se.B_IMAGE[k]/se.A_IMAGE[k],se.THETA_IMAGE[k])
        except AttributeError:
            print('WARNING: no sources detected in image!')
        input=('hit any key to continue \n')

    def add_simple_sersic_object(self,objnumber,profile,x,y,mag,rad,nsersic,BA,PA):
        self.galfit_input.write(' \n')
        self.galfit_input.write('# Object number: %i \n'%(objnumber))
        self.galfit_input.write(' 0) %s             # Object type \n'%(profile))
        self.galfit_input.write(' 1) %8.1f  %8.1f 1 1  # position x, y        [pixel] \n'%(x,y))
        self.galfit_input.write(' 3) %5.2f      1       # total magnitude     \n'%(mag))
        self.galfit_input.write(' 4) %8.2f       1       #     R_e              [Pixels] \n'%(rad))
        self.galfit_input.write(' 5) %5.2f       1       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic))
        self.galfit_input.write(' 9) %5.2f       1       # axis ratio (b/a)    \n'%(BA))
        self.galfit_input.write('10) %5.2f       1       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA))
        self.galfit_input.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
                                

    def close_input_file(self):
        self.galfit_input.close()

    def print_params(self):
        print('CURRENT INPUTS: \n mag = %5.2f %i \n Re = %5.2f %i \n n = %5.2f %i\n B/A = %5.2f %i \n PA = %5.2f %i \n fitall = %i \n fitcenter = %i \n'%(self.mag,self.fitmag,self.rad,self.fitrad,self.nsersic,self.fitn,self.BA,self.fitBA,self.PA,self.fitPA,self.fitallflag,self.fitcenter))


    def print_galfit_results(self, image):
        res = parse_galfit_results_dc(image, ncomp=self.ncomp, asymflag=self.asymmetry)
        for i, c in enumerate(res.components, start=1):
            print(f"Component {i}:")
            print(f"  XC  {c.xc.value:8.2f} +/- {c.xc.error:6.2f}")
            print(f"  YC  {c.yc.value:8.2f} +/- {c.yc.error:6.2f}")
            print(f"  MAG {c.mag.value:8.2f} +/- {c.mag.error:6.2f}")
            print(f"  RE  {c.re.value:8.2f} +/- {c.re.error:6.2f}")
            print(f"  N   {c.n.value:8.2f} +/- {c.n.error:6.2f}")
            print(f"  AR  {c.ar.value:8.2f} +/- {c.ar.error:6.2f}")
            print(f"  PA  {c.pa.value:8.2f} +/- {c.pa.error:6.2f}")
            print(f"  NUMERR {c.numerical_error_flag}")

        if res.sky is not None:
            print(f"Sky: {res.sky.value:.4g} +/- {res.sky.error:.4g}")
        if res.asymmetry is not None:
            print(f"F1: {res.asymmetry.f1.value:.4g} +/- {res.asymmetry.f1.error:.4g}")
            print(f"F1PA: {res.asymmetry.f1pa.value:.4g} +/- {res.asymmetry.f1pa.error:.4g}")

        print(f"ERROR: {res.error}")
        print(f"CHI2NU: {res.chi2nu}")

        
    def print_galfit_results_old(self,image):
        t, header_keywords = parse_galfit_results(image, ncomp=self.ncomp, asymflag=self.asymmetry, return_keywords=True)

        #if self.asymmetry:
        #    header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','1_F1','1_F1PA','ERROR','CHI2NU']
        #else:
        #    header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        #    #header_keywords = ['1_XC', '1_YC', '1_MAG', '1_RE', '1_N', '1_AR', '1_PA', '2_SKY', 'CHI2NU']
        #if self.ncomp == 2:
        #    header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_XC','2_YC','2_MAG','2_RE','2_N','2_AR','2_PA','3_SKY','ERROR','CHI2NU']

        print(f"in print_galfit_results, len(t) = {len(t)}, len(header_keywords) = {len(header_keywords)}")
        for i in range(len(header_keywords)):
            try:
                print('%6s : %5.2f +/- %5.2f'%(header_keywords[i],t[i][0],t[i][1]))
            except:
                print(f'WARNING: Problem parsing {header_keywords[i]:6s} : {t[i]}')
    def edit_params_menu(self):
        flag=str(input('What is wrong?\n o = nearby object (toggle fitall)  \n c = recenter \n f = hold values fixed \n a = toggle asymmetry parameter \n R = reset to original values \n g = go (run galfit) \n x=quit \n '))
        return flag

    def toggle_fitall(self):
        self.fitallflag=not(self.fitallflag)

    def toggle_asymmetry(self):
        self.asymmetry=not(self.asymmetry)

    def print_fix_menu(self):
        self.print_params()
        flag3=str(input('What do you want to hold fixed/toggle?\n n = fix sersic index \n r = fix Re \n b = fix B/A \n p = PA \n c = center \n f = use constraint file \n R = reset to original values \n g = go (run galfit) \n x=quit \n '))
        return flag3

    def fix_n(self):
        n=float(input('sersic exponent = '))
        self.nsersic=n
        self.fitn=not(self.fitn)

    def fix_rad(self):
        self.fitrad=not(self.fitrad)

    def fix_BA(self):
        self.fitBA=not(self.fitBA)
        print(self.fitBA, self.BA)
    def fix_PA(self):
        self.fitPA=not(self.fitPA)

    def fix_center(self):
        self.fitcenter=not(self.fitcenter)

    def add_constraint_file(self):
        self.constraintflag=not(self.constraintflag)
        

def run_galfit_headless(
    image,
    mask_image=None,
    sigma_image=None,
    psf_image=None,
    psf_oversampling=1,
    fit_box=None,                # (xmin, xmax, ymin, ymax)
    convolution_size=None,
    magzp=25.0,
    pscale=None,
    ncomp=1,
    asym=False,
    init=None,                   # list of component guesses
    fit_flags=None,              # which params to fit
    cwd=None,
    verbose=False,
):
    """
    Returns dict: {"out_fits": ..., "out_txt": ..., "log": ..., "results": parsed_dict}
    """
