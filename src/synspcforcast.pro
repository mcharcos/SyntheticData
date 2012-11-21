; NAME:
;     SYNSPCFORCAST - Version 1.0
;
; PURPOSE: Used to create synthetic spectra for FORCAST
;
; CLASS ATTRIBUTES:       
;
;     + chopamp
;     
;
; CLASS METHODS:
;     + INIT: Initialize object based on the initvals fits/text file or array. 
;     + CLEANUP: Clean heaps. Need to be implemented by child classes
;     + CreateImage: Create a synthetic 2D image based on the values of the object attributes
;                    after convolution by psf and with noise
;     
; INHEREITS METHODS:
;     + From SYNTHETICSPC
;
; MODIFICATION HISTORY:
;     Written by:  Miguel Charcos (mcharcos@sofia.usra.edu), USRA, June 25th 2012


;****************************************************************************
;     CreateImage - Create a synthetic 2D image based on the values of the 
;                   object attributes after convolution by psf and with noise
;****************************************************************************
FUNCTION synspcforcast::CreateImage, fname=fname, _EXTRA=extraProperties  ;, header=header
  
  xchop = -self.chopamp*sin(self.chopang)
  ychop = -self.chopamp*cos(self.chopang)
  xnod = -self.nodamp*sin(self.nodang)
  ynod = -self.nodamp*cos(self.nodang)
  self->Message,'Chop Amplitude: '+strtrim(self.chopamp,2)+' pixels, Nod Amplitude: '+strtrim(self.nodamp,2)+' pixels',priority='DEBUG',method='CreateImage'
  self->Message,'Chop Angle: '+strtrim(self.chopang*180./!pi ,2)+' degrees, Nod Angle: '+strtrim(self.nodang*180./!pi ,2)+' degrees',priority='DEBUG',method='CreateImage'

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
  plane1 = self->syntheticspc::CreateImage(_EXTRA=extraProperties)
  ; Create second plane (choped image)
  self->Message,'Creating plane 2, shift=['+strtrim(xchop,2)+','+strtrim(ychop,2)+']',priority='DEBUG',method='CreateImage'
  plane2 = self->syntheticspc::CreateImage(shiftsource=[xchop,ychop],bckplane=1, _EXTRA=extraProperties)
  ; Create second plane (nodded image)
  self->Message,'Creating plane 3, shift=['+strtrim(xnod,2)+','+strtrim(ynod,2)+']',priority='DEBUG',method='CreateImage'
  plane3 = self->syntheticspc::CreateImage(shiftsource=[xnod,ynod],bckplane=2, _EXTRA=extraProperties)
  ; Create second plane (choped+nodded image)
  self->Message,'Creating plane 4, shift=['+strtrim(xchop+xnod,2)+','+strtrim(ychop+ynod,2)+']',priority='DEBUG',method='CreateImage'
  plane4 = self->syntheticspc::CreateImage(shiftsource=[xchop+xnod,ychop+ynod],bckplane=3, _EXTRA=extraProperties)
  
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
    'STARE': begin
	  self->Message,'Creating a 4 plane array',priority='DEBUG',method='CreateImage'
          imres = [[[plane1]],[[plane2]]]
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
FUNCTION synspcforcast::init, nodamp=nodamp, chopamp=chopamp, nodang=nodang, chopang=chopang, ffits=ffits, $
                           _Extra=extraKeyword
  
  self.pixscale = 0.77
  
  if keyword_set(ffits) then res = self->syntheticspc::init(_Extra=extraKeyword,ffits=ffits) $
  else res = self->syntheticspc::init(_Extra=extraKeyword)
  
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
  ; WARNING!!!!!!!!!!!!!!  
  ; We need to take a look to this correction and add the non-linearity effect.
  ;if self.background gt 0. then begin
  ;  auxarr = drip_imgnonlin(replicate(1.,self.dimx,self.dimy),*self.header,siglev=self.background)
  ;  self.nlfactor = 1./auxarr[0,0]
  ;  self->Message,'Non linearity factor is '+strtrim(self.nlfactor,2)+' for background '+strtrim(self.background,2)
  ;endif
  
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
      chopamp_read = float(sxpar(header,'CHPAMP1'))
      if chopamp_read ne -999 then begin
	self.chopamp = chopamp_read*2./self.pixscale
	self.nodamp = float(sxpar(header,'NODAMP'))/self.pixscale
	self.chopang = float(sxpar(header,'CHPANGLE'))/180.*!pi
	self.nodang = float(sxpar(header,'NODANGLE'))/180.*!pi
	self.chopsys = float(sxpar(header,'CHPCOORD'))
	self.nodsys = float(sxpar(header,'NODCOORD'))
      endif else begin
	self.chopamp = 0.
	self.nodamp = 0.
	self.chopang = 0.
	self.nodang = 0.
	self.chopsys = 0.
	self.nodsys = 0.
      endelse
      
      self->Message,'Chop Amplitude: '+strtrim(self.chopamp,2)+' pixels, Nod Amplitude: '+strtrim(self.nodamp,2)+' pixels',priority='DEBUG',method='INIT'
      self->Message,'Chop Angle: '+strtrim(self.chopang*180./!pi ,2)+' degrees, Nod Angle: '+strtrim(self.nodang*180./!pi ,2)+' degrees',priority='DEBUG',method='INIT'
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
  
  ; Now create slits according to the observation mode
  gmode = 2
  if(strtrim(sxpar(*self.header,'INSTCFGS'),2) eq '0' or $
     strtrim(sxpar(*self.header, 'SPECTEL1'),2) eq '0') then begin
    self->Message,'Problem with header keywords',priority='ERROR',method='INIT'
    return,0
  endif
  if(sxpar(*self.header,'INSTCFGS') eq 'GXD5-8' and $
     sxpar(*self.header, 'SPECTEL1') eq 'FOR_XG063') then gmode=0

  if(sxpar(*self.header,'INSTCFGS') eq 'GRISM_XD' and $
     sxpar(*self.header, 'SPECTEL1') eq 'FOR_XG111') then gmode=1

  if(sxpar(*self.header,'INSTCFGS') eq 'GRISM' and $
     sxpar(*self.header, 'SPECTEL1') eq 'FOR_G063') then gmode=2

  if(sxpar(*self.header,'INSTCFGS') eq 'GRISM' and $
     sxpar(*self.header, 'SPECTEL1') eq 'FOR_G111') then gmode=3

  if(sxpar(*self.header,'INSTCFGS') eq 'GRISM' and $
     sxpar(*self.header, 'SPECTEL1') eq 'FOR_G227') then gmode=4

  if(sxpar(*self.header,'INSTCFGS') eq 'GRISM' and $
     sxpar(*self.header, 'SPECTEL1') eq 'FOR_G329') then gmode=5

  ; Check if it is high dispersion or longslit spectroscopy
  if gmode eq 0 or gmode eq 1 then begin
    R = 1000.
    Lambda0 = 5.
    slit1 = self->DefaultSlit(x=0,y=202.+19./2.,w0=7.18061,dw=Lambda0/R,width=3.1, length=19.,angslit=0.,angdisp=0) ;-0.200848)
    slit2 = self->DefaultSlit(x=0,y=163.+19./2.,w0=6.72309,dw=Lambda0/R,width=3.1, length=19.,angslit=0.,angdisp=0) ;-0.183677)
    slit3 = self->DefaultSlit(x=0,y=128.+19./2.,w0=6.34376,dw=Lambda0/R,width=3.1, length=19.,angslit=0.,angdisp=0) ;-0.176025)
    slit4 = self->DefaultSlit(x=0,y=97.+19./2. ,w0=5.99245,dw=Lambda0/R,width=3.1, length=19.,angslit=0.,angdisp=0) ;-0.168352)
    slit5 = self->DefaultSlit(x=0,y=70.+19./2. ,w0=5.68018,dw=Lambda0/R,width=3.1, length=19.,angslit=0.,angdisp=0) ;-0.156805)
    slit6 = self->DefaultSlit(x=0,y=45.+19./2. ,w0=5.39719,dw=Lambda0/R,width=3.1, length=19.,angslit=0.,angdisp=0) ;-0.152946)
    slit7 = self->DefaultSlit(x=0,y=23.+19./2. ,w0=5.117  ,dw=Lambda0/R,width=3.1, length=19.,angslit=0.,angdisp=0) ;-0.145216)
    slit8 = self->DefaultSlit(x=0,y=0.+17./2.  ,w0=4.847  ,dw=Lambda0/R,width=3.1, length=19.,angslit=0.,angdisp=0) ;-0.145216)
    if self.slits eq ptr_new() then self.slits = ptr_new([slit1,slit2,slit3,slit4,slit5,slit6,slit7,slit8]) $
    else *self.slits = [slit1,slit2,slit3,slit4,slit5,slit6,slit7,slit8]
  endif else begin
    slit1 = self->DefaultSlit(x=128,y=0,w0=6.0,width=3.1, length=19.,angslit=0.)
    if self.slits eq ptr_new() then self.slits = ptr_new([slit1]) $
    else *self.slits = [slit1]
  endelse
  
  self->UpdateMapDisp
  
  return, 1
  
END

;****************************************************************************
;     SYNFORCAST__DEFINE - Define the class structure for the class catalogue
;****************************************************************************

PRO synspcforcast__define

struct={synspcforcast, $
        chopamp: 3., $
	nodamp: 3., $
        chopang: 3., $
	nodang: 3., $
	chopsys: 0, $      ; Coordinage system of the chopping 0=SIRF 1=TARF 2=ERF
	nodsys: 0, $       ; Coordinage system of the nodding 0=SIRF 1=TARF 2=ERF
	pixscale: 0.77, $
        inherits syntheticspc}

END

