; NAME:
;     SYNTHETIC - Version 1.0
;
; PURPOSE: Parent class of the classes use to create synthetic images of different instruments
;
; CLASS ATTRIBUTES:       
;
;     + instrument:   Instrument that performed the observation
;     + filter:       Filter name
;     + obsmode:      Observation mode
;     + dimx, dimy:   Image array dimensions
;     + skyangle:     Angle of the detector in the sky
;     + image:        Pointer to an ideal image containing the objects without noise
;     + header:       Header of the fits file used to create the fits image
;     + sources:      Pointer to an array of structures defining the object
;   	    	    	--> Defined in child classes (Different for spectroscopy and imaging)
;                       --> error: defines if the source was defined correctly. Must be in the structure
;                                  that the child classes return.
;     + psf:          Pointer to an array containing the profile of the psf (in pixel units?)
;     + nlfactor:     Image non linearity factor
;     + snr:          Signal to noise
;     + fflat:        File name of the flat image
;     + fdark:        File name of the dark image
;     + fbadpix:      File name of the bad pixel mask
;     + fdistortion:  File name of the bad pixel mask
;     
;
; CLASS METHODS:
;     + INIT: Initialize object based on the initvals fits/text file or array. 
;     + CLEANUP: Clean heaps. Need to be implemented by child classes
;     + ImageSource: Must be implemented in child classes
;     + RefreshImage: Refresh the image attribute according to sources
;     + DefaultSource: Must be implemented in child classes
;     + AddSource:   Add a source to the array of sources
;     + CreateImage: Create a synthetic 2D image based on the values of the object attributes
;                    after convolution by psf and with noise
;     + CompareImage: Compare the input image to the theoretical image of the object (i.e. that does not contain noise)
;                     Returns the difference?
;     
; INHEREITS METHODS:
;     + From CLASSDEF
;
; MODIFICATION HISTORY:
;     Written by:  Miguel Charcos (mcharcos@sofia.usra.edu), USRA, May 16th 2012



;****************************************************************************
;     DISTCORR_MODEL - Create the distortion correction model using distcorr_model.pro
;                      used by drip
;****************************************************************************
FUNCTION synthetic::distcorr_model, image=image
  
  if self.fdistortion eq '' then begin
    self->Message,'Distortion file is not defined',priority='DEBUG',method='distcorr_model'
    if keyword_set(image) then return,image else return,-1
  endif
  if FILE_TEST(self.fdistortion) eq 0 then begin
    self->Message,'Distortion file does not exist',priority='DEBUG',method='distcorr_model'
    if keyword_set(image) then return,image else return,-1
  endif
  if keyword_set(image) eq 0 then begin
    if self.image eq ptr_new() then begin
      self->Message,'Image is empty',priority='DEBUG',method='distcorr_model'
      return,-1
    endif
    image = *self.image
  endif
  
  pin_npts = [12,12]
  spx=[3,3,3,2,5,5,6,6]
  spy=[1,2,3,3,5,6,5,6]
  
  pinpos = distcorr_model(self.fdistortion,NXIMG=self.dimx, NYIMG=self.dimy, $
                        NXPTS=pin_npts[0], NYPTS=pin_npts[1], SPX=spx, SPY=spy)  ;,/viewpin)
  
  xo = pinpos[1:pinpos[0]]
  yo = pinpos[pinpos[0]+1:2*pinpos[0]]
  xi = pinpos[2*pinpos[0]+1:3*pinpos[0]]
  yi = pinpos[3*pinpos[0]+1:4*pinpos[0]]
  
  if self.header eq ptr_new() then begin 
    self->Message,'Header is not completed',priority='WARNING',method='distcorr_model'
    self->Message,'Distortion has not been performed',priority='WARNING',method='distcorr_model'
    return,image
  endif 
  hdr = *self.header
  newimage = drip_undistort(image, hdr, PINPOS=[pinpos[0],xi,yi,xo,yo])
  
  nx = n_elements(newimage[*,0])
  ny = n_elements(newimage[0,*])
  
  nxmin = nx/2-self.dimx/2
  nymin = ny/2-self.dimy/2
  
  return, newimage[nxmin:nxmin+self.dimx-1,nymin:nymin+self.dimy-1]
  
END

;****************************************************************************
;     ImageSource - Return the image of the array of the input source
;****************************************************************************
FUNCTION synthetic::ImageSource,source, _EXTRA=extraProperties
  
  self->Message,'DefaultSource function should be defined in '+obj_class(self),priority='ERROR',method='SYNTHETIC::IMAGESORUCE'
    
  return,{error:1}
  
END

;****************************************************************************
;     AddBackground - Add the background to the image. By default it will
;                     be an homogeneous background but some intruments/modes
;                     may require differently. For example, high dispersion spc 
;                     or reduced field imaging.
;****************************************************************************
PRO synthetic::AddBackground

  if self.image eq ptr_new() then self.image = ptr_new(replicate(self.background, self.dimx,self.dimy)) $
  else *self.image = replicate(self.background, self.dimx,self.dimy)
  
END

;****************************************************************************
;     RefreshImage - Refresh the image attribute according to sources
;****************************************************************************
PRO synthetic::RefreshImage, nosources=nosources, _EXTRA=extraProperties
  
  self->Message,'Starting...',priority='DEBUG',method='RefreshImage'
  
  ; Verify that the dimension of the array is correct
  if self.dimx le 0 or self.dimy le 0. then begin
    self->Message,'Dimension must be greater than 0.',priority='ERROR',method='INIT'
    return
  endif
  
  self->AddBackground
  
  if self.sources eq ptr_new() then begin
    self->Message,'No source exist',priority='INFO',method='RefreshImage'
    return
  endif
  if n_elements(*self.sources) eq 0 then begin
    self->Message,'No source exist',priority='INFO',method='RefreshImage'
    return
  endif
  
  if keyword_set(nosources) then begin
    self->Message,'No sources were added to the image. Only background is returned',priority='INFO',method='RefreshImage'
    return
  endif
  
  self->Message,'Number of sources in object is '+strtrim(n_elements(*self.sources),2),priority='DEBUG',method='RefreshImage'
  for i=0,n_elements(*self.sources)-1 do begin
    self->Message,'Adding source '+strtrim(i+1,2)+' to image',priority='DEBUG',method='RefreshImage'
    auxim = self->ImageSource((*self.sources)[i], _EXTRA=extraProperties)
    if auxim.error eq 0 then *self.image = *self.image + auxim.image $
    else self->Message,'.... ERROR adding source',priority='DEBUG',method='RefreshImage'
  endfor
END


;****************************************************************************
;     DefaultSource  -  Returns a source with default values
;****************************************************************************
FUNCTION synthetic::DefaultSource, _EXTRA=extraProperties
  
  self->Message,'DefaultSource function should be defined in '+obj_class(self),priority='ERROR',method='SYNTHETIC::DEFAULTSORUCE'
  
  return, {error:1}
  
END


;****************************************************************************
;     AddSource  -  Add a source to the array of sources
;****************************************************************************
PRO synthetic::AddSource, _EXTRA=extraProperties
  
  source=self->DefaultSource(_EXTRA=extraProperties) 
  
  if self.sources eq ptr_new() then begin
    self->Message,'No source existed before (null pointer)',priority='INFO',method='AddSource'
    self->Message,'   adding current source.',priority='INFO',method='AddSource'
    self.sources = ptr_new([source])
    ;self->RefreshImage
    return
  endif
  if n_elements(*self.sources) eq 0 then begin
    self->Message,'No source exist',priority='INFO',method='AddSource'
    *self.sources = [source]
  endif else begin
    *self.sources = [*self.sources,source]
  endelse
  
  ;self->RefreshImage
  
END

;****************************************************************************
;     CreateImage - Create a synthetic 2D image based on the values of the 
;                   object attributes after convolution by psf and with noise
;****************************************************************************
FUNCTION synthetic::CreateImage, fname=fname, header=header, $
         dopsf=dopsf, nonoise=nonoise, bckplane=bckplane, _EXTRA=extraProperties
  
  ;if keyword_set(shiftsource) then begin
  ;  self->RefreshImage,shiftx=shiftsource[0],shifty=shiftsource[1]
  ;endif else begin
  ;  self->RefreshImage
  ;endelse
  
  if keyword_set(dopsf) then begin
    if dopsf ne 1 then self->RefreshImage, _EXTRA=extraProperties, dopsf=dopsf  $
    else self->RefreshImage, _EXTRA=extraProperties
  endif else begin
    self->RefreshImage, _EXTRA=extraProperties
  endelse
  
  ; First convolve the image with the psf
  imconv = *self.image 
  if keyword_set(dopsf) then begin
    self->Message,'DOPSF is set',priority='DEBUG',method='CreateImage'
    if dopsf eq 1 then begin
      self->Message,'     to 1',priority='DEBUG',method='CreateImage'
      if self.dispersion eq -1 then begin
	self->Message,'DISPERSION is set to -1',priority='DEBUG',method='CreateImage'
	if self.psf ne ptr_new() then begin
	  self->Message,'PSF array is initialized',priority='DEBUG',method='CreateImage'
	  imconv = convolve(*self.image,*self.psf)
	endif else begin
	  self->Message,'PSF array is empty',priority='DEBUG',method='CreateImage'
	  self->Message,'Creating PSF array of width '+strtrim(self.psffwhm,2),priority='DEBUG',method='CreateImage'
          self.psf = ptr_new(/allocate_heap)
	  *self.psf = psf_Gaussian( NPIXEL=[self.dimx,self.dimy], FWHM=self.psffwhm, /NORMALIZE, NDIMEN=2)
	  imconv = convolve(*self.image,*self.psf)
	endelse
      endif else begin
	imconv = *self.image
	if self.dispersion eq 0 or self.dispersion eq 1 then begin
          if self.psf ne ptr_new() then begin
	    cutpsf = (*self.psf)[self.dimx/2,*]
	  endif else begin
	    cutpsf = psf_Gaussian( NPIXEL=self.dimy, FWHM=self.psffwhm, /NORMALIZE, NDIMEN=1)
	  endelse
	  for i=0,n_elements(imconv[*,0])-1 do begin
	    imconv[i,*] = convol((*self.image)[i,*],cutpsf,1.,/edge_truncate)
	  endfor
	endif
	if self.dispersion eq 2 or self.dispersion eq 3 then begin
          if self.psf ne ptr_new() then begin
	    cutpsf = (*self.psf)[*,self.dimy/2]
	  endif else begin
	    cutpsf = psf_Gaussian( NPIXEL=self.dimx, FWHM=self.psffwhm, /NORMALIZE, NDIMEN=1)
	  endelse
          for i=0,n_elements(imconv[0,*])-1 do begin
	    imconv[*,i] = convol((*self.image)[*,i],cutpsf,1.,/edge_truncate)
	  endfor
	endif
      endelse
    endif  ; end if dopsf eq 1
  endif
  
  ; Scale by the conversion factor. In the future we should add something here about the exposure time
  self->Message,'Conversion factor is '+strtrim(self.conversion,2),priority='INFO',method='CreateImage'
  imconv = imconv/self.conversion
  sxaddpar,*self.header,'CONVFACT',self.conversion
  
  ; Add flat fields, darks and badpixels
  ; we may want to verify the size of the images when read from a file
  ; in order to match self.dimx and self.dimy
  imflat = replicate(1.,self.dimx,self.dimy)
  imdark = replicate(0.,self.dimx,self.dimy)
  imbpix = replicate(0.,self.dimx,self.dimy)
  
  if file_test(self.fbadpix) then begin
    imbpix = readfits(self.fbadpix,hbpix)
    s = size(imbpix)
    if s[0] eq 3 then begin
      auxbpix = total(imbpix,3)/s[3]
      imbpix = replicate(0.,self.dimx,self.dimy)
      k = where(auxbpix ne 0.)
      if k[0] ne -1 then imbpix[k] = 1.
    endif
  endif
  
  ; We may want to add noise to the images?!
  if file_test(self.fflat) then begin
    imflat = readfits(self.fflat,hflat)
    s = size(imflat)
    if s[0] eq 3 then begin
      case s[3] of
        4: begin
	   imflat = (imflat[*,*,0]+imflat[*,*,1])-(imflat[*,*,2]+imflat[*,*,3])
	 end
        2: begin
	   imflat = (imflat[*,*,0]-imflat[*,*,1])
	 end
        else: begin
	   imflat = imflat[*,*,0] 
	 end
      endcase
      k = where(imbpix eq 0)
      if k[0] ne -1 then imflat = imflat/median(imflat[k])
      k = where(imbpix eq 1)
      if k[0] ne -1 then imflat[k] = 1.
      
      ;imflat = total(imflat,3)/s[3]
    endif
  endif
  if file_test(self.fdark) then begin
    imdark = readfits(self.fdark,hdark)
    s = size(imdark)
    if s[0] eq 3 then begin
      imdark = total(imdark,3)/s[3]
    endif
  endif
  
  ;*self.image = (*self.image)*imflat+imdark
  ;imconv = imconv*imflat + imdark
  kbpix = where(imbpix eq 1) 
  if kbpix[0] ne -1 then imconv[k] = 0.
  
  ; Add random cosmic rays
  doCosmicRays = 0
  if doCosmicRays eq 1 then begin
    seed = 10.
    avgcrays = 2
    craysmin = 10000.
    craysmax = 10000.
    self->Message,'Averaged number of cosmic rays per image is '+strtrim(avgstrikes,2),priority='DEBUG',method='RefreshImage'
    self->Message,'Minimum value of cosmic rays is '+strtrim(craysmin,2),priority='DEBUG',method='RefreshImage'
    self->Message,'Maximum value of cosmic rays is '+strtrim(craysmax,2),priority='DEBUG',method='RefreshImage'
    ncrays = randomn(seed,POISSON=avgstrikes)
    craydns = long(randomu(seed,ncrays)*(maxcount-mincount)+mincount+0.5)
    x = (fix(randomu(seed,ncrays)*(self.dimx-1)+0.5) < (self.dimx-1)) > 0
    y = (fix(randomu(seed,ncrays)*(self.dimy-1)+0.5) < (self.dimy-1)) > 0

    imconv[x,y] = imconv[x,y] + craydns
  endif
  
  
  ; Add noise
  ;atv,imconv
  ;didscale = 0
  ;tmpscale = 1.e9
  ;if mean(imconv) lt 1./tmpscale then begin
  ;  didscale = 1
  ;  imconv = imconv*tmpscale
  ;endif
  
  ;k = where(imconv mod 1. ne 0.)
  ;if k[0] eq -1 then print,'LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL' else print,imconv[k]
  
  if not keyword_set(nonoise) then begin
    ;imres = POIDEV(imconv) 
    imnoise = sqrt(abs(imconv))
    kneg = where(imconv lt 0.)
    if kneg[0] ne -1 then imnoise[kneg] = 0.
    randf = 2.*randomu(seed,n_elements(imnoise))-1.
    auxidx = findgen(n_elements(randf))
    imnoise[auxidx] = randf[auxidx]*imnoise[auxidx]
    imres = imconv+imnoise
  endif else begin
    imres = imconv
  endelse
  ;if didscale eq 1 then begin help,
  ;  ;imconv = imconv/tmpscale
  ;  imres = imres/tmpscale
  ;endif
  
  
  ;atv,(imres-imconv)^2/imconv
  ;print,'MMMMMMMMMMMMMMMMMMMMMMMMMMMMM',mean((imres-imconv)^2/imconv),stddev((imres-imconv)^2/imconv)
  
  ; Now we add the non linearity effects due to the average counts of the image (=background)
  self->Message,'Linearity Correction Factor is '+strtrim(self.nlfactor,2),priority='INFO',method='CreateImage'
  imres = imres*self.nlfactor
  
  ; Apply distortion correction
  imres = self->distcorr_model(image=imres)
  
  
  ; Now add background image if set
  dobg = 0
  if self.fbackground ne '' then begin
    if file_test(self.fbackground) then begin
      imbackground = readfits(self.fbackground,/silent)
      sbg = size(imbackground)
      if sbg[0] ge 2 then begin
        if sbg[1] eq self.dimx and sbg[2] eq self.dimy then dobg = sbg[0]
      endif
    endif
  endif
  CASE dobg of
    0 : begin
        self->Message,'No background image was used as a reference',priority='INFO',method='RefreshImage'
      end
    2 : begin
        self->Message,'2D image was used as a background. File = '+self.fbackground,priority='INFO',method='RefreshImage'
	imres = imbackground + imres 
      end
    3 : begin
        self->Message,'First plane of 3D image was used as a background. File = '+self.fbackground,priority='INFO',method='RefreshImage'
	if not keyword_set(bckplane) then bckplane = 0
	if bckplane gt n_elements(imbackground[0,0,*])-1 then bckplane = 0
	imres = imbackground[*,*,bckplane] + imres 
      end
    else: begin
          self->Message,'This case is not handled by the code. Review synthetic.pro',priority='ERROR',method='RefreshImage'
        end
  endcase
  
  if keyword_set(fname) then begin
    self->Message,'Saving image to file'+fname,priority='INFO',method='CreateImage'
    ; Create basic header or update input header using values of object
    ; mkhdr,hdr,imconv + imnoise
    if self.header ne ptr_new() then begin
      namepos=strpos(fname,path_sep(),/reverse_search)
      if namepos ge 0  then sxaddpar,*self.header,'FILENAME',strmid(fname,namepos+1,strlen(fname)) $
      else sxaddpar,*self.header,'FILENAME',fname 
      writefits,fname,imres, *self.header 
    endif else begin
      writefits,fname,imres
    endelse
  endif
  
  ; We create the image with the default parameters: no shift and original flux
  ;if keyword_set(shiftsource) then begin
  ;  self->RefreshImage
  ;endif
  ;self->RefreshImage
  
  ;return,imnoise
  ;return,imconv
  return,imres
  
END

;****************************************************************************
;     CompareImage - Compare the input image to the theoretical image of the object 
;****************************************************************************
FUNCTION synthetic::CompareImage, cmpimg
    
  ; First convolve the image with the psf
  imconv = convolve( self.image, self.psf)
  
  ; We may want to return some values that help to analyze if they are similar within the error.
  
  return, cmpimg-imconv 
  
END

;****************************************************************************
;     CLEANUP - Call clean pointer heap variables. Requires implementation in child
;****************************************************************************
PRO synthetic::cleanup
  
  if self.image ne ptr_new() then ptr_free,self.image
  if self.header ne ptr_new() then ptr_free,self.header
  if self.sources ne ptr_new() then ptr_free,self.sources
  if self.psf ne ptr_new() then ptr_free,self.psf
  
END

;****************************************************************************
;     INIT - Initialize object 
;****************************************************************************
FUNCTION synthetic::init, dimx=dimx, dimy=dimy,snr=snr, psf=psf, background=background, $
                          instrument=instrument, obsmode=obsmode, filter=filter, skyangle=skyangle, $
			  ffits=ffits,fflat=fflat, fdark=fdark, fbadpix=fbadpix, fbackground=fbackground, $
			  _Extra=extraKeyword
  
  
  if self->classdef::init(_Extra=extraKeyword) eq 0 then begin
    self->Message,'Problem calling parent initition', priority='WARNING',method='INIT'
  endif
  
  ; We give a first initialization of the object attributes
  self.dimx = 256
  self.dimy = 256
  self.instrument = 'FORCAST'
  self.obsmode = 'NMC'
  self.filter = '24um'
  self.snr = 3.
  self.skyangle = 0.
  self.background = 0.
  self.conversion = 1.
  
  ; use ffits if set in the inputs
  if keyword_set(ffits) then begin
    if file_test(ffits) then begin
      im = readfits(ffits,header,/silent)
      self.filter = strtrim(sxpar(header,'WAVELNTH'),2)
      self.instrument = strtrim(sxpar(header,'INSTRUME'),2)
      self.obsmode = strtrim(sxpar(header,'INSTMODE'),2)
      self.skyangle = (180. -  float(sxpar(header,'SKY_ANGL')))/180.*!pi
      self.conversion = sxpar(header,'FRMRATE')*sxpar(header,'EPERADU')/1.e6
    endif else begin
      self->Message,'File does not exist: '+ffits,priority='ERROR',method='INIT'
      return,0
    endelse
  endif else begin
    mkhdr,header,replicate(0.,self.dimx,self.dimy)
  endelse
  sxaddpar,header,'HISTORY','--------- SYNTHETIC PARAMETERS -----------------'
  
  ; We add the file names of the instrument features if set
  if keyword_set(fflat) then self.fflat = fflat
  if keyword_set(fdark) then self.fdark = fdark
  if keyword_set(fbadpix) then self.fbadpix = fbadpix 
  if keyword_set(fdistortion) then self.fdistortion = fdistortion 
  if keyword_set(fbackground) then self.fbackground = fbackground 
  
  ; Now overwrite the values of the attributes with the inputs if set
  if keyword_set(dimx) then self.dimx = dimx
  if keyword_set(dimy) then self.dimy = dimy
  
  if keyword_set(instrument) then begin
    self.instrument = instrument
    sxaddpar,header,'INSTRUME',self.instrument
  endif
  if keyword_set(obsmode) then begin
    self.obsmode = obsmode
    sxaddpar,header,'INSTMODE',self.obsmode
  endif
  if keyword_set(filter) then begin
    self.filter = filter
    sxaddpar,header,'WAVELNTH',self.filter
  endif
  if n_elements(skyangle) gt 0 then begin
    self.skyangle = skyangle
    sxaddpar,header,'SKY_ANGL',180.-self.skyangle*180./!pi
  endif
  
  if keyword_set(snr) then self.snr = snr
  if n_elements(background) gt 0 then self.background = background
  
  ; Image non linearity factor is set to 1. because the calculation is specific 
  ; to each instrument
  self.nlfactor = 1.
  
  ; Calculate or read psf
  ;if keyword_set(psf) eq 0 then begin
  ;  ;psf = fltarr(dimx,dimy)
  ;  ; add gaussian
  ;  psf = psf_Gaussian( NPIXEL=[self.dimx,self.dimy], FWHM=4., /NORMALIZE, NDIMEN=2)
  ;  self.psf = ptr_new(psf)
  ;endif 
  if keyword_set(psf) then begin
    spsf = size(psf)
    self->Message,'PSF was initialized in the input parameters.',priority='INFO',method='INIT'
    if spsf[0] eq 2 then self.psf = ptr_new(psf)
  endif
  if keyword_set(psffwhm) then self.psffwhm = psffwhm else self.psffwhm = 6.
  
  
  
  if self.dimx le 0 or self.dimy le 0. then begin
    self->Message,'Dimension must be greater than 0.',priority='ERROR',method='INIT'
    return,0
  endif
  
  self.image = ptr_new(replicate(self.background, self.dimx,self.dimy))
  
  sxaddpar,header,'NAXIS1',self.dimx
  sxaddpar,header,'NAXIS2',self.dimy
  sxaddpar,header,'PROCSTAT','LEVEL_1'
  sxdelpar,header,'PRODTYPE'
  ;sxdelpar,header,'PIN_MOD'
  sxaddpar,header,'HISTORY','Synthetic Image created using '+OBJ_CLASS(self)
  self.header = ptr_new(header)
  
  self.dispersion = -1
  
  return, 1
  
END

;****************************************************************************
;     SYNTHETIC__DEFINE - Define the class structure for the class catalogue
;****************************************************************************

PRO synthetic__define

struct={synthetic, $
        instrument:'', $ 
        filter:'', $ 
        obsmode:'', $
	dimx: 0, $
	dimy: 0, $
	skyangle:0., $
	image: ptr_new(), $
	header: ptr_new(), $
        sources:ptr_new(), $   
        psf:ptr_new(), $    
	psffwhm: 4.,  $    ; used only if psf is not defined
	conversion:1. ,$   ; Conversion factor (=EPERADU*FRMRATE)
	background: 0., $  ; background in e/s
	nlfactor: 1., $    ; image non linearity factor
        snr:3., $
	fflat:'', $
	fdark: '', $
	fbadpix: '', $
	fdistortion: '', $
	fbackground: '', $
	dispersion:0, $       ; dispersion is set to -1 for imaging
        inherits classdef}

END

