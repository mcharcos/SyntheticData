FUNCTION drip_getpar, header, keyword
  
  res = sxpar(header,keyword)
  
  if strtrim(res,2) eq '0' then return,'x'
  
  return,res
  
END
