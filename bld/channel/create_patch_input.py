# ------------------------------------------------------------------------------
# Copyright (C) 2021 Aleksi Nummelin
#
# This file is part of BLOM.
#
# BLOM is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
# more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BLOM. If not, see <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------------
#
# The purpose of this script is to generate patch.input file for idealized
# topography, which can then be used with the blom_dimensions script to create
# the dimensions.F BLOM module, containing all dimension parameters and grid
# partitioning information. Note that the domain size and chosen splitting will
# determine the processor count.
#
################################################## 
# MODIFY THE FOLLOWING PARAMETERS AS DESIRED
idm = 208
jdm = 512
ibig = 13 
jbig = 8
nreg = 1
##################################################
#
npe = idm//ibig
mpe = jdm//jbig
npes = npe*mpe
# OUTPUT PATCH FILENAME
fout = open('patch.input.'+str(npes)+'_python','w')
#
minsea = ibig*jbig-ibig
maxsea = ibig*jbig
avesea = int((minsea*2*idm/ibig+maxsea*(npes-2*idm/ibig))//npes)
#
ispt = list(range(1,idm+1,ibig))
iipe = [ibig]*len(ispt)
#
jspt=list(range(1,jdm+1,jbig))
jjpe=[jbig]*len(jspt)
#
ispt_str=''
iipe_str=''
for jj,j in enumerate(ispt):
   ispt_str=ispt_str+str(j).rjust(5)
   iipe_str=iipe_str+str(ibig).rjust(5)
   if (jj+1)%8==0 and (jj!=0 and jj<len(ispt)-1):
       ispt_str=ispt_str+'\n'+' '.ljust(12)
       iipe_str=iipe_str+'\n'+' '.ljust(12)

jspt_str=''
jjpe_str=''
for jj,j in enumerate(jspt):
   jspt_str=jspt_str+str(j).rjust(5)
   jjpe_str=jjpe_str+str(jbig).rjust(5)
   if (jj+1)%8==0 and (jj!=0 and jj<len(jspt)-1):
       jspt_str=jspt_str+'\n'+' '.ljust(12)
       jjpe_str=jjpe_str+'\n'+' '.ljust(12)

lineout0 = '  npes   npe   mpe   idm   jdm  ibig  jbig  nreg  minsea  maxsea  avesea\n'
lineout0 = lineout0+str(npes).rjust(6)+str(npe).rjust(6)+str(mpe).rjust(6)+str(idm).rjust(6)+str(jdm).rjust(6)+str(ibig).rjust(6)+str(jbig).rjust(6)+str(nreg).rjust(6)+str(minsea).rjust(8)+str(maxsea).rjust(8)+str(avesea).rjust(8)
fout.write(lineout0+'\n')
fout.write('\n')
#
for j in range(mpe):
   lineout1 = 'ispt('+str(j+1).rjust(3)+') = '+ispt_str
   lineout2 = 'iipe('+str(j+1).rjust(3)+') = '+iipe_str
   fout.write(lineout1+'\n')
   fout.write(lineout2+'\n')

fout.write('\n')
fout.write('\n')
lineout1 = 'jspt('+str(1).rjust(3)+') = '+jspt_str
lineout2 = 'jjpe('+str(1).rjust(3)+') = '+jjpe_str
fout.write(lineout1+'\n')
fout.write(lineout2+'\n')
