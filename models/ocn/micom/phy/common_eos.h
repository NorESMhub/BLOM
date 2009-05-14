c
c --- common blocks related to the equation of state
c
      c o m m o n /eos/
c
c --- coefficients for the functional fit of in situ density
c
     . a11,a12,a13,a14,a15,a16,b11,b12,b13
     .,a21,a22,a23,a24,a25,a26,b21,b22,b23
c
c --- reference pressure
     .,pref
c
c --- coefficients for potential density in sigma units with reference
c --- pressure at -pref-
     .,ap11,ap12,ap13,ap14,ap15,ap16
     .,ap21,ap22,ap23,ap24,ap25,ap26
c
c --- coefficients for potential density in sigma units with reference
c --- pressure at the surface
     .,ap110,ap120,ap130,ap140,ap150,ap160
     .,ap210,ap220,ap230,ap240,ap250,ap260
c
      real a11,a12,a13,a14,a15,a16,b11,b12,b13
     .    ,a21,a22,a23,a24,a25,a26,b21,b22,b23
     .    ,pref
     .    ,ap11,ap12,ap13,ap14,ap15,ap16
     .    ,ap21,ap22,ap23,ap24,ap25,ap26
     .    ,ap110,ap120,ap130,ap140,ap150,ap160
     .    ,ap210,ap220,ap230,ap240,ap250,ap260
