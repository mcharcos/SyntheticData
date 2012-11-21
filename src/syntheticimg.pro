; NAME:
;     SYNTHETIC - Version 1.0
;
; PURPOSE: Parent class of the classes use to create synthetic images of different instruments
;
; CLASS ATTRIBUTES: 
;     + sources:      Pointer to an array of structures defining the object
;   	    	    	--> x, y: position of the object
;                       --> flux: theoretical observed flux of the source
;   	    	    	--> type: source or extended  
;   	    	    	--> extension: size of the object for extended sources
;                       --> error: defines if the source was defined correctly.
;     
;
; CLASS METHODS:
;     + ImageSource: Return the image of the array of the input source
;     + DefaultSource:   Returns a source with default values
;     
; INHEREITS METHODS:
;     + From SYNTHETIC
;
; MODIFICATION HISTORY:
;     Written by:  Miguel Charcos (mcharcos@sofia.usra.edu), USRA, May 16th 2012




;****************************************************************************
;     ImageSource - Return the image of the array of the input source
;****************************************************************************
FUNCTION syntheticimg::ImageSource,source, fluxsource=fluxsource, shiftsource=shiftsource, shiftx=shiftx, shifty=shifty
                                        
  
  ; Verify that the dimension of the array is correct
  if self.dimx le 0 or self.dimy le 0. then begin
    self->Message,'Dimension must be greater than 0.',priority='ERROR',method='INIT'
    return,{error:1,image:[0.,0.]}
  endif
  
  if keyword_set(shiftsource) then begin
    shiftx = shiftsource[0]
    shifty = shiftsource[1]
  endif
  if keyword_set(shiftx) eq 0 then shiftx = 0.
  if keyword_set(shifty) eq 0 then shifty = 0.
  if n_elements(fluxsource) eq 0 then fluxsource = source.flux
  
  image = replicate(0., self.dimx,self.dimy)
  
  xsource = source.x + shiftx
  ysource = source.y + shifty
  if xsource lt 0 or xsource ge self.dimx or ysource lt 0 or ysource ge self.dimy then begin
    self->Message,'Source is outside of image',priority='DEBUG',method='ImageSource'
    return,{error:0,image:image}
  endif
  
  
  self->Message,'Adding '+source.type+' source ('+strtrim(xsource,2)+','+strtrim(ysource,2)+') to ('+strtrim(self.dimx,2)+','+strtrim(self.dimy,2)+') image',priority='DEBUG',method='ImageSource'
  CASE strlowcase(source.type) of
    'point': begin
             self->Message,'    Flux = '+strtrim(fluxsource,2),priority='DEBUG',method='ImageSource'
	     image[xsource,ysource] = fluxsource 
	   end
    else: begin
          self->Message,'Source type '+source.type+' does not exist',priority='ERROR',method='INIT'
	  return,{error:1,image:[0.,0.]}
	end
  endcase
  
  return,{error:0,image:image}
  
END


;****************************************************************************
;     DefaultSource  -  Returns a source with default values
;****************************************************************************
FUNCTION syntheticimg::DefaultSource, _EXTRA=extraProperties
  
  defaultprop = {x:self.dimx/2, y:self.dimy/2, $
                 flux:500., $    ; flux in e/s
		 type:'point',extension:0., error:0}
  
  return, fillstructure(defaultprop,_EXTRA=extraProperties)
  
END

;****************************************************************************
;     SYNTHETICIMG__DEFINE - Define the class structure for the class catalogue
;****************************************************************************

PRO syntheticimg__define

struct={syntheticimg, $
        inherits synthetic}

END

