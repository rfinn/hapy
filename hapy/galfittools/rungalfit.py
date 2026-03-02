from __future__ import annotations
'''
python wrapper to run galfit and parse results

'''

#import pyds9
import os
from dataclasses import dataclass, field
from typing import Optional, List, Tuple, Dict

from pathlib import Path
import subprocess
import shutil

from astropy.io import fits

from hapy.galfittools.results import GalfitComponent, GalfitResult


def _parse_galfit_header_pair(raw) -> Tuple[float, float, int]:
    """
    Parse a GALFIT header value into (value, error, numerical_error_flag).

    Handles:
      - fixed params like "[12.34]" -> (12.34, 0.0)
      - normal like "12.34 +/- 0.56"
      - numerical issues like "*12.34 +/- *0.56" -> numerical_error_flag=1
      - bare floats like "1.234" -> (1.234, 0.0)
    """
    numerical_flag = 0
    s = str(raw).strip()
    temp = None

    # Fixed parameters: [value]
    if "[" in s and "]" in s:
        s2 = s.replace("[", "").replace("]", "").strip()
        try:
            return float(s2), 0.0, 0
        except ValueError:
            return float("nan"), 0.0, 0

    # Typical: "val +/- err"
    if "+/-" in s:
        left, right = [t.strip() for t in s.split("+/-", 1)]

        if "*" in left or "*" in right:
            numerical_flag = 1
            left = left.replace("*", "")
            right = right.replace("*", "")

        try:
            return float(left), float(right), numerical_flag
        except ValueError:
            return float("nan"), 0.0, numerical_flag

    # Bare numeric
    try:
        return float(s), 0.0, 0
    except ValueError:
        return float("nan"), 0.0, 0


def parse_galfit_results_dc(
    galfit_outimage: str,
    ncomp: int = 1,
    asymflag: bool = False,
) -> GalfitResult:
    """
    Parse GALFIT output FITS into a GalfitResult dataclass (canonical results.py).

    Expects extension 2 header contains keys like:
      1_XC, 1_YC, 1_MAG, 1_RE, 1_N, 1_AR, 1_PA
      2_SKY (for 1-comp) or 3_SKY (for 2-comp), etc.
      CHI2NU
      ERROR
    If asymflag:
      1_F1, 1_F1PA
    """
    hdr = fits.getheader(galfit_outimage, 2)

    def comp(ci: int) -> GalfitComponent:
        numflag = 0

        xc, xc_err, nf = _parse_galfit_header_pair(hdr.get(f"{ci}_XC", "nan"))
        numflag = max(numflag, nf)
        yc, yc_err, nf = _parse_galfit_header_pair(hdr.get(f"{ci}_YC", "nan"))
        numflag = max(numflag, nf)
        mag, mag_err, nf = _parse_galfit_header_pair(hdr.get(f"{ci}_MAG", "nan"))
        numflag = max(numflag, nf)
        re, re_err, nf = _parse_galfit_header_pair(hdr.get(f"{ci}_RE", "nan"))
        numflag = max(numflag, nf)
        n, n_err, nf = _parse_galfit_header_pair(hdr.get(f"{ci}_N", "nan"))
        numflag = max(numflag, nf)

        # GALFIT uses AR (axis ratio b/a) in the header
        ba, ba_err, nf = _parse_galfit_header_pair(hdr.get(f"{ci}_AR", "nan"))
        numflag = max(numflag, nf)

        pa, pa_err, nf = _parse_galfit_header_pair(hdr.get(f"{ci}_PA", "nan"))
        numflag = max(numflag, nf)

        return GalfitComponent(
            xc=xc, yc=yc, mag=mag, re=re, n=n, ba=ba, pa=pa,
            xc_err=xc_err, yc_err=yc_err, mag_err=mag_err, re_err=re_err,
            n_err=n_err, ba_err=ba_err, pa_err=pa_err,
            numerical_error_flag=numflag,
        )

    comp1 = comp(1)
    comp2: Optional[GalfitComponent] = comp(2) if ncomp >= 2 else None

    # Sky is usually (ncomp+1)_SKY in your convention
    sky_key = f"{ncomp + 1}_SKY"
    sky, sky_err, _ = _parse_galfit_header_pair(hdr.get(sky_key, "0.0+/-0.0"))

    # CHI2NU
    try:
        chi2nu = float(str(hdr.get("CHI2NU", 0.0)).split()[0])
    except Exception:
        chi2nu = 0.0

    # ERROR keyword (GALFIT sometimes writes numeric or string; your results.py stores float)
    try:
        error = float(str(hdr.get("ERROR", 0.0)).split()[0])
    except Exception:
        error = 0.0

    # Asymmetry terms
    f1 = f1_err = f1_pa = f1_pa_err = None
    if asymflag:
        if "1_F1" in hdr:
            f1, f1_err, _ = _parse_galfit_header_pair(hdr["1_F1"])
        if "1_F1PA" in hdr:
            f1_pa, f1_pa_err, _ = _parse_galfit_header_pair(hdr["1_F1PA"])

    return GalfitResult(
        ncomp=ncomp,
        comp1=comp1,
        comp2=comp2,
        sky=sky,
        sky_err=sky_err,
        chi2nu=chi2nu,
        error=error,
        f1=f1,
        f1_err=f1_err,
        f1_pa=f1_pa,
        f1_pa_err=f1_pa_err,
    )



def parse_galfit_results(galfit_outimage, asymflag=0, ncomp=1, return_keywords=False):
    """
    Backwards-compatible wrapper.

    Returns the old 'list of tuples' format that galfitwrapper currently expects.
    """
    res = parse_galfit_results_dc(galfit_outimage, ncomp=ncomp, asymflag=bool(asymflag))

    fit_parameters = []
    header_keywords = []

    # component 1 always
    fit_parameters += [
        (res.comp1.xc, res.comp1.xc_err),
        (res.comp1.yc, res.comp1.yc_err),
        (res.comp1.mag, res.comp1.mag_err),
        (res.comp1.re, res.comp1.re_err),
        (res.comp1.n, res.comp1.n_err),
        (res.comp1.ba, res.comp1.ba_err),
        (res.comp1.pa, res.comp1.pa_err),
    ]
    header_keywords += ["1_XC","1_YC","1_MAG","1_RE","1_N","1_AR","1_PA"]

    # component 2 if present
    if ncomp >= 2 and res.comp2 is not None:
        fit_parameters += [
            (res.comp2.xc, res.comp2.xc_err),
            (res.comp2.yc, res.comp2.yc_err),
            (res.comp2.mag, res.comp2.mag_err),
            (res.comp2.re, res.comp2.re_err),
            (res.comp2.n, res.comp2.n_err),
            (res.comp2.ba, res.comp2.ba_err),
            (res.comp2.pa, res.comp2.pa_err),
        ]
        header_keywords += ["2_XC","2_YC","2_MAG","2_RE","2_N","2_AR","2_PA"]

    # sky
    fit_parameters.append((res.sky, res.sky_err))
    header_keywords.append(f"{ncomp+1}_SKY")

    # asym terms
    if bool(asymflag):
        fit_parameters.append((res.f1 if res.f1 is not None else float("nan"),
                               res.f1_err if res.f1_err is not None else 0.0))
        header_keywords.append("1_F1")
        fit_parameters.append((res.f1_pa if res.f1_pa is not None else float("nan"),
                               res.f1_pa_err if res.f1_pa_err is not None else 0.0))
        header_keywords.append("1_F1PA")

    # ERROR and CHI2NU
    fit_parameters.append(res.error)
    header_keywords.append("ERROR")
    fit_parameters.append(res.chi2nu)
    header_keywords.append("CHI2NU")

    if return_keywords:
        return fit_parameters, header_keywords
    return fit_parameters




class RunGalfit:
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

        
    def disable_convolution(self):
        self.convflag = False
    def enable_convolution(self):
        self.convflag = True
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
        #print('inside rungalfit, fitBA = ',self.fitBA)

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
        #print('sersic n, fitsersicn = ',self.nsersic,self.fitn)
        if nsersic == None:
            self.galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(self.nsersic,int(self.fitn)))
        else:
            self.galfit_input.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic,int(self.fitn)))
        #print('BA, fitBA = ',self.BA,self.fitBA)
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
        #print('in rungalfit.run_galfit, self.psf_image = ',self.psf_image)
        
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



    def run_galfit(self, displayflag: bool = False):
        """
        Run GALFIT using subprocess, capture output, and rename outputs.

        Sets:
          self.galfit_flag = 1 on success, 0 on failure
          self.galfit_log, self.galfit_out
        """
        self.write_input_file()

        self.galfit_flag = 0  # pessimistic default

        # GALFIT command
        cmd = ["galfit", self.galfile]

        # Run GALFIT
        try:
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,  # we handle returncode ourselves for better logging
            )
        except FileNotFoundError as e:
            print("GALFIT executable not found (is it on PATH?)")
            print(e)
            return False
        except Exception as e:
            print("Unexpected error launching GALFIT")
            print(e)
            return False

        # Save stdout/stderr to a per-run log file (very helpful for CV crashes)
        #image_id = f"{self.galname}-"
        #runlog = Path(f"{image_id}{self.ncomp}Comp-galfit.stdout_stderr.txt")
        #try:
        #    runlog.write_text(
        #        "COMMAND:\n"
        #        + " ".join(cmd)
        #        + "\n\nSTDOUT:\n"
        #        + (proc.stdout or "")
        #        + "\n\nSTDERR:\n"
        #        + (proc.stderr or "")
        #    )
        #except Exception:
        #    pass

        if proc.returncode != 0:
            print(f"GALFIT failed (returncode={proc.returncode}). See {runlog}")
            return False

        # Verify outputs exist
        fitlog_src = Path("fit.log")
        out_src = Path("galfit.01")

        if not out_src.exists():
            print(f"GALFIT returned 0 but {out_src} is missing. See {runlog}")
            return False

        # Rename/copy outputs deterministically
        self.galfit_log = f"{image_id}{self.ncomp}Comp-fit.log"
        self.galfit_out = f"{image_id}{self.ncomp}Comp-galfit.01"

        try:
            if fitlog_src.exists():
                shutil.move(str(fitlog_src), self.galfit_log)
        except Exception as e:
            print(f"Warning: could not move fit.log -> {self.galfit_log}: {e}")

        try:
            shutil.move(str(out_src), self.galfit_out)
        except Exception as e:
            print(f"GALFIT produced {out_src} but could not move to {self.galfit_out}: {e}")
            return False

        # Success
        self.galfit_flag = 1

        if displayflag:
            self.display_results()
        return True

        
    def run_galfit_old(self,displayflag=False):
        self.write_input_file()
        #print 'self.fitall = ',self.fitall
        s = 'galfit '+self.galfile
        #print('run the following: ',s)


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

    def run_galfit_old(self,displayflag=False):
        self.write_input_file()
        #print 'self.fitall = ',self.fitall
        s = 'galfit '+self.galfile
        #print('run the following: ',s)


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

        c = res.comp1
        print("Component 1:")
        print(f"  XC  {c.xc:8.2f} +/- {c.xc_err:6.2f}")
        print(f"  YC  {c.yc:8.2f} +/- {c.yc_err:6.2f}")
        print(f"  MAG {c.mag:8.2f} +/- {c.mag_err:6.2f}")
        print(f"  RE  {c.re:8.2f} +/- {c.re_err:6.2f}")
        print(f"  N   {c.n:8.2f} +/- {c.n_err:6.2f}")
        print(f"  BA  {c.ba:8.2f} +/- {c.ba_err:6.2f}")
        print(f"  PA  {c.pa:8.2f} +/- {c.pa_err:6.2f}")
        print(f"  NUMERR {c.numerical_error_flag}")

        if res.comp2 is not None:
            c = res.comp2
            print("Component 2:")
            # same print block...

        print(f"Sky: {res.sky:.4g} +/- {res.sky_err:.4g}")
        print(f"ERROR: {res.error}")
        print(f"CHI2NU: {res.chi2nu}")
        if self.asymmetry and res.f1 is not None:
            print(f"F1: {res.f1:.4g} +/- {res.f1_err:.4g}")
            print(f"F1PA: {res.f1_pa:.4g} +/- {res.f1_pa_err:.4g}")


 

        
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
        
    def run_and_parse(self) -> GalfitResult:
        runok = self.run_galfit()


        if not runok:
            runlog = getattr(self, "galfit_runlog", None)
            msg = f"GALFIT failed for {getattr(self, 'galname', '')} (ncomp={self.ncomp})"
            if runlog:
                msg += f"; see {runlog}"
            raise RuntimeError(msg)

        if not Path(self.output_image).exists():
            raise RuntimeError(f"GALFIT reported success but output_image missing: {self.output_image}")
        return parse_galfit_results_dc(self.output_image,ncomp=self.ncomp,asymflag=self.asymmetry,)
