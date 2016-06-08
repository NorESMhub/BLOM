c--------------------------------------------------------------------
c Arrays to keep a two time-level copy of sediment fields
c These array are copied back and forth in micom2hamocc.F
c and hamocc2micom.F in the same manner as the tracer field.
c Also, they written/read to and from restart files.
c There are probably more efficient/elegant solutions.
c 
      c o m m o n  /bgcs_hamocc/
     . sedlay2(idm,jdm,2*ks,nsedtra) ,powtra2(idm,jdm,2*ks,npowtra)
     .,burial2(idm,jdm,2,   nsedtra)
c
      real sedlay2,powtra2,burial2

