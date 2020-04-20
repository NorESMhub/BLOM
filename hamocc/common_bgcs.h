c Copyright (C) 2020  J. Schwinger
c
c This file is part of BLOM/iHAMOCC.
c
c BLOM is free software: you can redistribute it and/or modify it under the
c terms of the GNU Lesser General Public License as published by the Free 
c Software Foundation, either version 3 of the License, or (at your option) 
c any later version. 
c
c BLOM is distributed in the hope that it will be useful, but WITHOUT ANY 
c WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
c FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
c more details. 
c
c You should have received a copy of the GNU Lesser General Public License 
c along with BLOM. If not, see https://www.gnu.org/licenses/.


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

