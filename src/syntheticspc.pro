; NAME:
;     SYNTHETIC - Version 1.0
;
; PURPOSE: Parent class of the classes use to create synthetic images of different instruments
;
; CLASS ATTRIBUTES: 
;     + sources:      Pointer to an array of structures defining the object
;   	    	    	--> x,y: position of of the spectrum in the detector.
;                       --> w0: wavelength of the spectrum at x,y position
;                       --> dw: dispertion wavelength of the spectrum per pixel
;                       --> flux: pointer to array of the spectrum along the dispersion direction
;                       --> lambda: array of the wavelength for each flux array element (same size as flux)
;                       --> type: file (for user-given spectrum) or gaussian (for created spectrum with background and line)
;   	    	    	--> extension: properties of spectrum (Array A => A3 + A0*e{-[(X-A1)/A2]^2})
;                       --> fname: If the spectrum was made from a file given by user then here we store the path+file name of the file
;                       --> error: defines if the source was defined correctly.
;     
;     + dispersion:   Dispersion direction 
;                       --> 0 = x left to right, 
;                       --> 1 = x right to left 
;                       --> 2 = y left to right 
;                       --> 3 = y right to left
;
; CLASS METHODS:
;     + ImageSource: Return the image of the array of the input source
;     + DefaultSource:   Returns a source with default values
;     
; INHEREITS METHODS:
;     + From SYNTHETIC
;
; MODIFICATION HISTORY:
;     Written by:  Miguel Charcos (mcharcos@sofia.usra.edu), USRA, June 14th 2012


;****************************************************************************
;     IsInSlit - Return the index of the slit if the source is in or -1
;****************************************************************************
FUNCTION syntheticspc::IsInSlit, x, y
  
  if self.slits eq ptr_new() then begin
    self->Message,'Slit definition pointer is not defined',priority='ERROR',method='IsInSlit'
    return,-1
  endif
  
  for i=0,n_elements(*self.slits)-1 do begin
    slit = (*self.slits)[i]
    a = tan(slit.angslit)
    if a eq 0. then begin 
      dwidth = 0.
      dlength = abs(y-slit.y)
    endif else begin
      c = -(slit.y+a*slit.x)
      dwidth = abs((a*x+y+c)/sqrt(a^2+1.))
      
      a2 = -1./a
      c2 = slit.x/a-slit.y
      dlength = abs((a2*x+y+c2)/sqrt(a2^2+1.))
    endelse
    ;mkct,0
    ;window,0,retain=2
    ;absice = findgen(self.dimx)
    ;plot,absice,-a*absice-c,xrange=[0,256.],yrange=[0.,256.]
    ;oplot,[x],[y],psym=2,color=2
    self->Message,'Distances of source ['+strtrim(x,2)+','+strtrim(y,2)+'] to slit '+strtrim(i,2)+' is Dwidth='+strtrim(dwidth,2)+', DLength='+strtrim(dlength,2)+'.'
    if slit.width/2. gt dwidth and slit.length/2. gt dlength then begin
      self->Message,'Star found in slit of index '+strtrim(i,2),priority='INFO',method='IsInSlit'
      return,i
    endif
  endfor
  
  self->Message,'Star ['+strtrim(x,2)+','+strtrim(y,2)+'] is not in slit',priority='INFO',method='IsInSlit'
  return,-1
  
END


;****************************************************************************
;     UpdateMapDisp - Update the mapdisp array containing 1 where light 
;                     that is dispersed from the slits in the focal plane reach
;                     and 0 where no.
;****************************************************************************
PRO synthetic::UpdateMapDisp
  
  auxim = fltarr(self.dimx,self.dimy)
  
  if self.slits eq ptr_new() then begin
    self->Message,'No slit is input.',priority='INFO',method='UpdateMapDisp'
    return
  endif
  
  for i=0,n_elements(*self.slits)-1 do begin
    refim = fltarr(3.*self.dimx,3.*self.dimy)
    currentslit = (*self.slits)[i]
    refim[*,self.dimy*1.5-currentslit.length/2.:self.dimy*1.5+currentslit.length/2.] = 1.
    rotrefim = rot(refim,currentslit.angdisp*180./!pi,1.,1.5*self.dimx,1.5*self.dimy,/pivot)  ;, cubic=-0.5, missing=0)
    auxim = auxim + ArrayInsert(auxim,rotrefim,at=[currentslit.x,currentslit.y],pixel=[1.5*self.dimx,1.5*self.dimy])
  endfor
  
  if self.mapdisp eq ptr_new() then self.mapdisp = ptr_new(auxim) $
  else *self.mapdisp = auxim
  
END

;****************************************************************************
;     AddBackground - Add the background to the image. By default it will
;                     be an homogeneous background but some intruments/modes
;                     may require differently. For example, high dispersion spc 
;                     or reduced field imaging.
;****************************************************************************
PRO synthetic::AddBackground
  
  auxim = fltarr(self.dimx,self.dimy)
  
  if self.mapdisp ne ptr_new() then begin 
    k = where(*self.mapdisp eq 1)
    if k[0] ne -1 then auxim[k] = self.background
  endif
  
  if self.image eq ptr_new() then self.image = ptr_new(auxim) $
  else *self.image = auxim
  
END

;****************************************************************************
;     UpdateSource - Return the image of the array of the input source
;****************************************************************************
FUNCTION syntheticspc::UpdateSource,source, force=force
  
  if not keyword_set(force) then begin
    if lambda ne ptr_new() or flux ne ptr_new() then begin
      self->Message,'Sectrum is already filled',priority='INFO',method='UpdateSource'
      return,1
    endif
  endif
  
  CASE strlowcase(source.type) of
    'file': begin
             ; read file and update source flux and lambda
	     if source.fname eq '' then begin
	       self->Message,'File name is not set',priority='WARNING',method='UpdateSource'
	       return,0
	     endif
	     if file_test(source.fname) eq 0 then begin
	       self->Message,'File does not exist: '+source.fname,priority='WARNING',method='UpdateSource'
	       return,0
	     endif
	     
	     ; open file and fill source.lambda and source.flux
	     im = readfits(source.fname,h)
	     
	     s = size(im)
	     switch s[0] of
	       1: begin
		  fluxarr = im

		  ; Read the headers and figure out the initial wavelength and the step
		  readw0 = sxpar(h,'CRVAL1')
		  if readw0 eq 0 then begin
		    self->Message,'Problem reading CRVAL1',priority='WARNING',method='UpdateSource'
		    return,0
		  endif
		  readdw = sxpar(h,'CDELT1')
		  if readdw eq 0 then begin
		    self->Message,'Problem reading CDELT1',priority='WARNING',method='UpdateSource'
		    return,0
		  endif
		  lambda = readw0 + findgen(s[1])*readdw
		  break
	        end
	       2: 
	       3: 
	       4: begin
		  fluxarr = im[0,*]
	          fluxarr = fluxarr*(source.extension)[0] + (source.extension)[3]
		  lambda = im[1,*]
		  break
	        end
	      else: begin
	            self->Message,'This case is not handled. Problem calling method.',priority='ERROR',method='UpdateSource'
		    return,0
		  end
	    endswitch
	   end
    'gaussian': begin
		if source.w0 eq 0 then begin
		  self->Message,'W0 keyword is not set in source',priority='WARNING',method='UpdateSource'
		  return,0
		endif
		if source.dw eq 0 then begin
		  self->Message,'DW keyword is not set in source',priority='WARNING',method='UpdateSource'
		  return,0
		endif
               ; create a gaussian within the array
	       Nlambda = 2*max([self.dimx,self.dimy])
	       lambda = source.w0 + (findgen(Nlambda)-Nlambda/2.)*source.dw
	       fluxarr = gaussian(lambda,source.extension)
	   end
    'sav': begin
	   ; Example of Vega spectrum from Bill Vacca
	   ;FC2VIN          FLOAT     = Array[427977] - Don't remember
	   ;FCVIN           FLOAT     = Array[427977] - Continuum flux of Vega
	   ;FVIN            FLOAT     = Array[427977] - Flux of Vega (in ergs/cm2/s/A, I think...)
	   ;WVIN            FLOAT     = Array[427977] - wavelength (in microns, I believe)
             ; read file and update source flux and lambda
	     if source.fname eq '' then begin
	       self->Message,'File name is not set',priority='WARNING',method='UpdateSource'
	       return,0
	     endif
	     if file_test(source.fname) eq 0 then begin
	       self->Message,'File does not exist: '+source.fname,priority='WARNING',method='UpdateSource'
	       return,0
	     endif
	   restore,source.fname
	   lambda = wvin/1000.
	   fluxarr = fvin
	   fluxarr = fluxarr*(source.extension)[0] + (source.extension)[3]
	 end
    else: begin
          self->Message,'Source type '+source.type+' does not exist',priority='ERROR',method='UpdateSource'
	  return,{error:1,image:[0.,0.]}
	end
  endcase
  
  if n_elements(fluxarr) ne n_elements(lambda) then begin
    self->Message,'Dimension of FLUX array ('+strtrim(n_elements(fluxarr),2)+') does not match dimension of LAMBDA array ('+strtrim(n_elements(self.lambda),2)+')',priority='ERROR',method='UpdateSource'
    return,0
  endif
  
  if source.flux eq ptr_new() then begin
    source.flux = ptr_new(fluxarr)
    source.lambda = ptr_new(lambda)
  endif else begin
    if source.lambda eq ptr_new() then begin
      self->Message,'LAMBDA pointer was null while FLUX was filled',priority='ERROR',method='UpdateSource'
      return,0
    endif
    *source.flux = fluxarr
    *source.lambda = lambda
  endelse
  
  ; Fill the rest of the attributes of source
  ;source.w0 = (*source.lambda)[0]
  ;source.dw = (*source.lambda)[1]-(*source.lambda)[0]
  
  ; Update trace coefficients and calibration coefficients if it was not set before
  CalCoeff = [4.847, 0.011535156, 0.0, 0.0]
  if source.CalCoeff eq ptr_new() then source.CalCoeff = ptr_new(CalCoeff) ; else *source.CalCoeff = CalCoeff
  TraceCoeff = [0.] ;,0.1,0.001]
  if source.TraceCoeff eq ptr_new() then source.TraceCoeff = ptr_new(TraceCoeff) ; else *source.TraceCoeff = TraceCoeff
  
  return,1
  
END


;****************************************************************************
;     ImageSource - Return the image of the array of the input source
;****************************************************************************
FUNCTION syntheticspc::ImageSource,source, fsubimage=fsubimage, $
                        shiftsource=shiftsource, shiftx=shiftx, shifty=shifty, dopsf=dopsf
    
  ; Verify that the dimension of the array is correct
  if self.dimx le 0 or self.dimy le 0. then begin
    self->Message,'Dimension must be greater than 0.',priority='ERROR',method='ImageSource'
    return,{error:1,image:[0.,0.]}
  endif
  
  image = replicate(0., self.dimx, self.dimy)
  
  if self->UpdateSource(source, /force) eq 0 then begin
    self->Message,'Source was not correctly updated',priority='WARNING',method='ImageSource'
    return,{error:1,image:image}
  endif
  
  if keyword_set(shiftsource) then begin
    shiftx = shiftsource[0]
    shifty = shiftsource[1]
  endif
  if keyword_set(shiftx) eq 0 then shiftx = 0.
  if keyword_set(shifty) eq 0 then shifty = 0.
  
  xsource = round(source.x + shiftx)
  ysource = round(source.y + shifty)
  
  ; Check if source is in slit
  inslit = self->IsInSlit(xsource,ysource)
  if inslit eq -1 then  begin
    self->Message,'Star outside slit',priority='INFO',method='ImageSource'
    return,{error:0,image:image}
  endif
  currentslit = (*self.slits)[inslit]
  
  
  ; Create psf array for convolution of psf if necessary
  if keyword_set(dopsf) then begin
    if self.dispersion eq 0 or self.dispersion eq 1 then begin
      if self.psf ne ptr_new() then begin
	if dopsf eq 2 then cutpsf = (*self.psf)[self.dimx/2,*] $
	else cutpsf = *self.psf
      endif else begin
	if dopsf eq 2 then cutpsf = psf_Gaussian( NPIXEL=self.dimy, FWHM=self.psffwhm, /NORMALIZE, NDIMEN=1) $
	else cutpsf = psf_Gaussian( NPIXEL=self.dimy, FWHM=[currentslit.width,self.psffwhm], /NORMALIZE, NDIMEN=2)
      endelse
    endif
    if self.dispersion eq 2 or self.dispersion eq 3 then begin
      if self.psf ne ptr_new() then begin
	if dopsf eq 2 then cutpsf = (*self.psf)[*,self.dimy/2] $ 
	else cutpsf = (*self.psf)
      endif else begin
	if dopsf eq 2 then cutpsf = psf_Gaussian( NPIXEL=self.dimx, FWHM=self.psffwhm, /NORMALIZE, NDIMEN=1) $
	else cutpsf = psf_Gaussian( NPIXEL=self.dimx, FWHM=[self.psffwhm,currentslit.width], /NORMALIZE, NDIMEN=2)
      endelse
    endif
  endif
  
  ; Smooth spectra for the resolution of the slit/disperser
  dlambda = (*source.lambda)[1] - (*source.lambda)[0]
  fluxslit = *source.flux  ;smooth(*source.flux,round(currentslit.dw/dlambda))
  
  ; Create the sub image with the spectrum
  self->Message,'Adding '+source.type+' source ('+strtrim(xsource,2)+','+strtrim(ysource,2)+') to ('+strtrim(self.dimx,2)+','+strtrim(self.dimy,2)+') image',priority='DEBUG',method='ImageSource'
  CASE self.dispersion of 
    0: begin
       subimage = MkSpc(*source.lambda, flux=fluxslit, dimx=2.*self.dimx, dimy=2.*self.dimx, traceCoeff=(*source.traceCoeff), CalCoeff=(*source.CalCoeff), subpix=1000)
       imconv = subimage.image
       if keyword_set(dopsf) then begin
	 if dopsf eq 2 then begin
	   for i=0,n_elements(imconv[*,0])-1 do begin
	     imconv[i,*] = convol(reform((subimage.image)[i,*]),cutpsf,1.,/edge_truncate)
	   endfor
	 endif else begin
	   imconv = convol(subimage.image,cutpsf,1.,/edge_truncate)
	 endelse
       endif
       addarr = imconv
     end
    1: begin
       subimage = MkSpc(*source.lambda, flux=fluxslit, dimx=2.*self.dimx, dimy=2.*self.dimx, traceCoeff=(*source.traceCoeff), CalCoeff=(*source.CalCoeff), subpix=1000)
       imconv = subimage.image
       if keyword_set(dopsf) then begin
	 for i=0,n_elements(imconv[*,0])-1 do begin
	   imconv[i,*] = convol((subimage.image)[i,*],cutpsf,1.,/edge_truncate)
	 endfor
       endif
       addarr = reverse(subimage.image)
     end
    2: begin
       subimage = MkSpc(*source.lambda, flux=fluxslit, dimx=2.*self.dimy, dimy=2.*self.dimy, traceCoeff=(*source.traceCoeff), CalCoeff=(*source.CalCoeff), subpix=1000)
       imconv = subimage.image
       if keyword_set(dopsf) then begin
	 for i=0,n_elements(imconv[0,*])-1 do begin
	   imconv[*,i] = convol((*self.image)[*,i],cutpsf,1.,/edge_truncate)
	 endfor
       endif
              
       ;addarr = rot(subimage.image,-90.)
       ;addarr = transpose(reverse(subimage.image,2))
       addarr = transpose(subimage.image)
       ;fluxaux = (*source.flux)  ;gaussian(*source.lambda,[1.,7.,1.])
       ;window,1,retain=2
       ;plot,*source.lambda,fluxaux,xrange=[3,10]
       ;subimage = MkSpc(*source.lambda, flux=fluxaux, dimx=2.*self.dimy, dimy=source.dslit, traceCoeff=(*source.traceCoeff), CalCoeff=(*source.CalCoeff), subpix=1000)
     end
    3: begin
       subimage = MkSpc(*source.lambda, flux=fluxslit, dimx=2.*self.dimy, dimy=2.*self.dimy, traceCoeff=(*source.traceCoeff), CalCoeff=(*source.CalCoeff), subpix=1000)
       imconv = subimage.image
       if keyword_set(dopsf) then begin
	 for i=0,n_elements(imconv[0,*])-1 do begin
	   imconv[*,i] = convol((*self.image)[*,i],cutpsf,1.,/edge_truncate)
	 endfor
       endif
       addarr = reverse(subimage.image)
       addarr = rotate(addarr,1)
     end
    else: begin
          self->Message,'Dispersion '+strtrim(self.dispersion,2)+' is not available',priority='ERROR',method='ImageSource'
	  return,{error:1,image:[0.,0.]}
        end
  ENDCASE
  
  if subimage.error eq 1 then begin
    self->Message,'Problem inserting image in spectra',priority='ERROR',method='ImageSource'
    return,{error:1,image:image}
  endif
    
  ; I am going to set all the -1 values of subimage.wavemin to a high numbe
  ; so it does not get the indexes with -1 when looking for the pixel for w0
  wavemin = subimage.wavemin
  wavemax = subimage.wavemax
  kmin = where(wavemin eq -1)
  if kmin[0] ne -1 then (wavemin)[kmin] = 99999.
  
  ; find w0 in the array 
  k = where(wavemin le currentslit.w0 and wavemax ge currentslit.w0)
  if k[0] eq -1 then begin
    self->Message,'No pixel found cointaining wavelength '+strtrim(currentslit.w0,2),priority='WARNING',method='ImageSource'
    atarr = [0,n_elements(subimage.image[0.,*])/2]
  endif else begin
    nsubcols = n_elements(subimage.image[*,0.])
    atarr = [k[0] mod nsubcols,k[0]/nsubcols]
  endelse
  
  ;window,5,retain=2
  ;plot,subimage.image[0:255,atarr[1]]
  ;oplot,[atarr[0]],[subimage.image[atarr[0],atarr[1]]],psym=4
    
  CASE self.dispersion of 
    1: atarr[0] = subimage.image[*,0] - atarr[0]
    2: atarr = [atarr[1],atarr[0]]
    3: atarr = [n_elements(subimage.image[0,*])-atarr[1],n_elements(subimage.image[*,0])-atarr[0]]
    else: 
  ENDCASE
  self->Message,'Using pixel ['+strjoin(strtrim(atarr,2),',')+']',priority='WARNING',method='ImageSource'
  
  if atarr[0] lt 0. then pepe
  
  ;rotate around that pixel
  if currentslit.angdisp eq 0. then begin
    spcimage = addarr 
  endif else begin
    spcimage = rot(addarr,currentslit.angdisp*180./!pi,1.,atarr[0],atarr[1],/pivot, cubic=-0.5, missing=0, /interp)
  endelse
  ;atv,spcimage
  ;atvplot,atarr[0],atarr[1],psym=4
  
  ; Use the atarr index in the subimage based on w0 as well as the xsource and ysource indexes
  ; to calculate the position where the array should be inserted
  
  image = ArrayInsert(image,spcimage,at=[xsource,ysource],pixel=atarr)
  
  if keyword_set(fsubimage) then begin
    self->Message,'Saving intermediate image to '+fsubimage,priority='INFO',method='ImageSource'
    writefits,fsubimage,image,*self.header
  endif  
  
  return,{error:0,image:image}
  
END


;****************************************************************************
;     DefaultSource  -  Returns a source with default values
;****************************************************************************
FUNCTION syntheticspc::DefaultSource, _EXTRA=extraProperties
  
  ; w0 and dw are used in the gaussian case when the spectrum is made from scratch, no from reading a file
  
  defaultprop = {x:self.dimx/2, y:self.dimy/2, $
                 w0:0., dw:0., $  ;range:[5.0,8.0],  $   ; Central wavelength, wave resolution and range
		 CalCoeff:ptr_new(), TraceCoeff:ptr_new(), $
		 flux:ptr_new(), $    ; flux in e/s
                 lambda:ptr_new(), $    ; same unit as w0
		 type:'gaussian',extension:[1.,0.,1.,0.], fname:'', $
		 error:0}
  
  ;if self.dispersion eq 0 or self.dispersion eq 1 then defaultprop.dslit = self.dimy/2.
  
  return, fillstructure(defaultprop,_EXTRA=extraProperties)
  
END

;****************************************************************************
;     DefaultSlit  -  Returns a slit with default values
;****************************************************************************
FUNCTION syntheticspc::DefaultSlit, _EXTRA=extraProperties
  
  defaultprop = {x:self.dimx/2, y:self.dimy/2, angslit:0., $
                 length:self.dimx, width:2., $  ; length and width of the slit in pixels
		 w0:0., $          ; Wavelength (microns) for the slit at pixel x,y
		 dw:3./256., $     ; Wavelength resolution (microns)
		 angdisp:0.,  $    ; Anlge of the dispersion with respect to column or row (depending on dispersion)
		 error:0}
  
  if self.dispersion eq 0 or self.dispersion eq 1 then defaultprop.length = self.dimy
  
  return, fillstructure(defaultprop,_EXTRA=extraProperties)
  
END

;****************************************************************************
;     INIT - Initialize object 
;****************************************************************************
FUNCTION syntheticspc::init,  _Extra=extraKeyword
  
  res = self->synthetic::init(_Extra=extraKeyword)
  
  if res eq 0 then begin
    self->Message,'Problem calling parent initition', priority='WARNING',method='INIT'
    return,0
  endif
  
  ; For now, I only implement long slit
  slit1 = self->DefaultSlit(x=128,y=0,w0=6.0,width=4.,angslit=0.)
  self.slits = ptr_new([slit1])
  
  self->UpdateMapDisp
  
  self.dispersion = 0  ; This is by default but I guess each instrument has its own
  
  return,1
  
END

;****************************************************************************
;     SYNTHETICIMG__DEFINE - Define the class structure for the class catalogue
;****************************************************************************

PRO syntheticspc__define

struct={syntheticspc, $
        slits: ptr_new(), $ ; list of slits (in the detector). See DefaultSlit for definition
	mapdisp: ptr_new(), $
	inherits synthetic}

END

