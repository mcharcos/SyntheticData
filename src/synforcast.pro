; NAME:
;     SYNFORCAST - Version 1.0
;
; PURPOSE: Parent class of the classes use to create synthetic images of different instruments
;
; CLASS ATTRIBUTES:       
;
;     + instrument:   Instrument that performed the observation
;     + filter:       Filter name
;     + obsmode:      Observation mode
;     + dimx, dimy:   Image array dimensions
;     + image:        Pointer to an ideal image containing the objects without noise
;     + sources:      Pointer to an array of structures defining the object
;   	    	    	--> x, y: position of the object
;   	    	    	--> type: source or extended  
;   	    	    	--> extension: size of the object for extended sources
;     + psf:          Pointer to an array containing the profile of the psf (in pixel units?)
;     + snr:          Signal to noise
;     
;
; CLASS METHODS:
;     + INIT: Initialize object based on the initvals fits/text file or array. 
;     + CLEANUP: Clean heaps. Need to be implemented by child classes
;     + CreateImage: Create a synthetic 2D image based on the values of the object attributes
;                    after convolution by psf and with noise
;     
; INHEREITS METHODS:
;     + From CLASSDEF
;
; MODIFICATION HISTORY:
;     Written by:  Miguel Charcos (mcharcos@sofia.usra.edu), USRA, May 16th 2012


;****************************************************************************
;     CreateImage - Create a synthetic 2D image based on the values of the 
;                   object attributes after convolution by psf and with noise
;****************************************************************************
FUNCTION synforcast::CreateImage, fname=fname  ;, header=header
  
  xchop = -self.chopamp*sin(self.chopang)
  ychop = -self.chopamp*cos(self.chopang)
  xnod = -self.nodamp*sin(self.nodang)
  ynod = -self.nodamp*cos(self.nodang)
  self->Message,'Chop Amplitude: '+strtrim(self.chopamp,2)+' pixels, Nod Amplitude: '+strtrim(self.nodamp,2)+' pixels',priority='DEBUG',method='CreateImage'
  self->Message,'Chop Angle: '+strtrim(self.chopang*180./!pi,2 )+' degrees, Nod Angle: '+strtrim(self.nodang*180./!pi,2 )+' degrees',priority='DEBUG',method='CreateImage'
  
  ; Change x and y of the chop and nod shifts according to the coordinate system
  if (self.chopsys eq 2) then begin
    self->Message,'Chopping coordinate system is ERF',priority='DEBUG',method='CreateImage'
    erf_xchop = (xchop*cos(self.skyangle) + ychop*sin(self.skyangle))
    erf_ychop = (ychop*cos(self.skyangle) - xchop*sin(self.skyangle))
    xchop = erf_xchop
    ychop = erf_ychop
  endif
  if (self.nodsys eq 2) then begin
    self->Message,'Nodding coordinate system is ERF',priority='DEBUG',method='CreateImage'
    erf_xnod = (xnod*cos(self.skyangle) + ynod*sin(self.skyangle))
    erf_ynod = (ynod*cos(self.skyangle) - xnod*sin(self.skyangle))
    xnod = erf_xnod
    ynod = erf_ynod
  endif
  
  ; Create first plane 
  self->Message,'Creating plane 1, shift=[0,0]',priority='DEBUG',method='CreateImage'
  plane1 = self->synthetic::CreateImage(/dopsf)
  ; Create second plane (choped image)
  self->Message,'Creating plane 2, shift=['+strtrim(xchop,2)+','+strtrim(ychop,2)+']',priority='DEBUG',method='CreateImage'
  plane2 = self->synthetic::CreateImage(shiftsource=[xchop,ychop],/dopsf)
  ; Create second plane (nodded image)
  self->Message,'Creating plane 3, shift=['+strtrim(xnod,2)+','+strtrim(ynod,2)+']',priority='DEBUG',method='CreateImage'
  plane3 = self->synthetic::CreateImage(shiftsource=[xnod,ynod],/dopsf)
  ; Create second plane (choped+nodded image)
  self->Message,'Creating plane 4, shift=['+strtrim(xchop+xnod,2)+','+strtrim(ychop+ynod,2)+']',priority='DEBUG',method='CreateImage'
  plane4 = self->synthetic::CreateImage(shiftsource=[xchop+xnod,ychop+ynod],/dopsf)
  
  self->Message,'Observation mode is '+self.obsmode,priority='DEBUG',method='CreateImage'
  CASE strupcase(self.obsmode) of
    'C2': begin
	  self->Message,'Creating a 2 plane array',priority='DEBUG',method='CreateImage'
	  imres = [[[plane1]],[[plane2]]]
        end
    'C2N': begin
	  self->Message,'Creating a 4 plane array',priority='DEBUG',method='CreateImage'
          imres = [[[plane1]],[[plane2]],[[plane3]],[[plane4]]]
        end
    else : begin
	   self->Message,'Creating a single plane array with plane 1',priority='DEBUG',method='CreateImage'
	   imres = plane1
         end
  ENDCASE
  
  if keyword_set(fname) then begin
    self->Message,'Saving image to file'+fname,priority='INFO',method='CreateImage'
    ; Create basic header or update input header using values of object
    if self.header ne ptr_new() then begin
      namepos=strpos(fname,path_sep(),/reverse_search)
      if namepos ge 0  then sxaddpar,*self.header,'FILENAME',strmid(fname,namepos+1,strlen(fname)) $
      else sxaddpar,*self.header,'FILENAME',fname 
      writefits,fname,imres, *self.header 
      writefits,'plane1.fits',plane1, *self.header 
      writefits,'plane2.fits',plane2, *self.header 
      writefits,'plane3.fits',plane3, *self.header 
      writefits,'plane4.fits',plane4, *self.header 
    endif else begin
      writefits,fname,imres
    endelse
    ;writefits,fname,imres
  endif
  
  return,imres
  
END

;****************************************************************************
;     INIT - Initialize object 
;****************************************************************************
FUNCTION synforcast::init, nodamp=nodamp, chopamp=chopamp, nodang=nodang, chopang=chopang, $
                           ffits=ffits,  _Extra=extraKeyword
  
  self.pixscale = 0.77
  
  if keyword_set(ffits) then res = self->synthetic::init(_Extra=extraKeyword, ffits=ffits) $
  else res = self->synthetic::init(_Extra=extraKeyword)
  
  if res eq 0 then begin
    self->Message,'Problem calling parent initition', priority='WARNING',method='INIT'
    return,0
  endif
  ; At this stage header is initiated.
  
  if self.instrument ne 'FORCAST' and self.instrument ne '' then begin
    self->Message,'Instrument should be FORCAST ('+self.instrument+')',priority='WARN',method='INIT'
  endif
  self.instrument = 'FORCAST'
  
  
  ; Calculate image non linearity factor
  if self.background gt 0. then begin
    auxarr = drip_imgnonlin(replicate(1.,self.dimx,self.dimy),*self.header,siglev=self.background)
    self.nlfactor = 1./auxarr[0,0]
    self->Message,'Non linearity factor is '+strtrim(self.nlfactor,2)+' for background '+strtrim(self.background,2),priority='INFO',method='INIT'
  endif
  
  self.chopamp = 30.
  self.nodamp = 30.
  self.chopang = 0.
  self.nodang = 0.
  self.chopsys = 0
  self.nodsys = 0
  
  ; use ffits if set in the inputs
  if keyword_set(ffits) then begin
    self->Message,'Using file '+ffits,priority='DEBUG',method='INIT'
    if file_test(ffits) then begin
      im = readfits(ffits,header,/silent)
      self.chopamp = float(sxpar(header,'CHPAMP1'))*2./self.pixscale
      self.nodamp = float(sxpar(header,'NODAMP'))/self.pixscale
      self.chopang = float(sxpar(header,'CHPANGLE'))/180.*!pi
      self.nodang = float(sxpar(header,'NODANGLE'))/180.*!pi
      self.chopsys = float(sxpar(header,'CHPCOORD'))
      self.nodsys = float(sxpar(header,'NODCOORD'))
      
      self->Message,'Chop Amplitude: '+strtrim(self.chopamp,2)+' pixels, Nod Amplitude: '+strtrim(self.nodamp,2)+' pixels',priority='DEBUG',method='INIT'
      self->Message,'Chop Angle: '+strtrim(self.chopang*180./!pi,2 )+' degrees, Nod Angle: '+strtrim(self.nodang*180./!pi,2 )+' degrees',priority='DEBUG',method='INIT'
    endif else begin
      self->Message,'File does not exist: '+ffits,priority='ERROR',method='INIT'
      return,0
    endelse
  endif
  
  ;self.chopang = 0.
  ;self.nodang = 90./180.*!pi
  ;self.chopamp = 0.
  ;self.nodamp = 0.
  
  if n_elements(chopamp) gt 0 then begin
    self.chopamp = chopamp
    sxaddpar,*self.header,'CHPAMP1',self.chopamp*self.pixscale/2.
  endif 
  if n_elements(nodamp) gt 0 then begin
    self.nodamp = nodamp
    sxaddpar,*self.header,'NODAMP',self.nodamp*self.pixscale
  endif 
  if n_elements(chopang) gt 0 then begin
    self.chopang = chopang
    sxaddpar,*self.header,'CHPANGLE',self.chopang*180./!pi
  endif 
  if n_elements(nodang) gt 0 then begin
    self.nodang = nodang
    sxaddpar,*self.header,'NODANGLE',self.nodang*180./!pi
  endif 
  if n_elements(chopsys) gt 0 then begin
    self.chopsys = chopsys
    sxaddpar,*self.header,'CHPCOORD',self.chopsys
  endif 
  if n_elements(nodsys) gt 0 then begin
    self.nodsys = nodsys
    sxaddpar,*self.header,'NODCOORD',self.nodsys
  endif 
  
  self.dispersion = -1
  
  return, 1
  
END

;****************************************************************************
;     SYNFORCAST__DEFINE - Define the class structure for the class catalogue
;****************************************************************************

PRO synforcast__define

struct={synforcast, $
        chopamp: 3., $
	nodamp: 3., $
        chopang: 3., $
	nodang: 3., $
	chopsys: 0, $      ; Coordinage system of the chopping 0=SIRF 1=TARF 2=ERF
	nodsys: 0, $       ; Coordinage system of the nodding 0=SIRF 1=TARF 2=ERF
	pixscale: 0.77, $
        inherits syntheticimg}

END

