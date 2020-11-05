# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include Makefile.config

CPPDEFS = -DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE

.DEFAULT:
	-touch $@
all: blom
advect.o: ../phy/advect.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/advect.F
barotp.o: ../phy/barotp.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/barotp.F
bigrid.o: ../phy/bigrid.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/bigrid.F
blkdat.o: ../phy/blkdat.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/blkdat.F
blom.o: ../drivers/nocoupler/blom.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/nocoupler/blom.F
blom_init.o: ../phy/blom_init.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/param_bud.h ../phy/common_bud.h mod_xc.o mod_swabs.o mod_dia.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/blom_init.F
blom_step.o: ../phy/blom_step.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_bud.h ../phy/common_bud.h ../phy/timing.h mod_xc.o mod_swabs.o mod_dia.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/blom_step.F
chksum.o: ../phy/chksum.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/chksum.F
chksummsk.o: ../phy/chksummsk.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/chksummsk.F
chkvar.o: ../phy/chkvar.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/chkvar.F
clndr.o: ../phy/clndr.F ../phy/common_blocks.h ../phy/common_clndr.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/clndr.F
cmnfld.o: ../phy/cmnfld.F ../phy/common_blocks.h mod_xc.o mod_dia.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/cmnfld.F
cmnfld_bfsqf.o: ../phy/cmnfld_bfsqf.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/cmnfld_bfsqf.F
cmnfld_nslope.o: ../phy/cmnfld_nslope.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/cmnfld_nslope.F
convec.o: ../phy/convec.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o eosfun.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/convec.F
crc.o: ../phy/crc.c ../phy/crc.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/crc.c
data_mct.o: ../cesm/data_mct.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../cesm/data_mct.F
diapfl.o: ../phy/diapfl.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_tke.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/diapfl.F
difest.o: ../phy/difest.F ../phy/common_blocks.h ../phy/common_seaice.h ../phy/common_forc.h ../phy/common_geo.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_tke.h ../phy/common_tke.h mod_xc.o mod_nctools.o mod_dia.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/difest.F
diffus.o: ../phy/diffus.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/diffus.F
dimensions.o: ../phy/dimensions.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/dimensions.F
eddtra.o: ../phy/eddtra.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/eddtra.F
eddtra_gm.o: ../phy/eddtra_gm.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/eddtra_gm.F
eddtra_intdif.o: ../phy/eddtra_intdif.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/eddtra_intdif.F
eosfun.o: ../phy/eosfun.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/eosfun.F
external_abort.o: ../drivers/nocoupler/external_abort.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../drivers/nocoupler/external_abort.F
fill_global.o: ../phy/fill_global.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/fill_global.F
forfun_channel.o: ../phy/forfun_channel.F ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/forfun_channel.F
frcdat.o: ../phy/frcdat.F ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/frcdat.F
fwbbal.o: ../phy/fwbbal.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/fwbbal.F
geoenv_cesmextra.o: ../cesm/geoenv_cesmextra.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../cesm/geoenv_cesmextra.F
geoenv_channel.o: ../phy/geoenv_channel.F ../phy/common_blocks.h ../phy/common_geo.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/geoenv_channel.F
geoenv_file.o: ../phy/geoenv_file.F ../phy/common_blocks.h ../phy/common_geo.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/geoenv_file.F
geoenv_test.o: ../phy/geoenv_test.F ../phy/common_blocks.h ../phy/common_geo.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/geoenv_test.F
getfrc.o: ../phy/getfrc.F ../phy/common_blocks.h mod_xc.o mod_ben02.o mod_channel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/getfrc.F
getfrc_cesm.o: ../cesm/getfrc_cesm.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_seaice.h ../cesm/common_cesm.h ../phy/common_clndr.h mod_xc.o mod_nctools.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../cesm	../cesm/getfrc_cesm.F
getlpi.o: ../phy/getlpi.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/getlpi.F
gettime.o: ../phy/gettime.F ../phy/timing.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/gettime.F
hfharm.o: ../phy/hfharm.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/hfharm.F
icaux_ben02.o: ../ben02/icaux_ben02.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_seaice.h mod_xc.o mod_ben02.o mod_thdysi.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../ben02/icaux_ben02.F
icaux_cesm.o: ../cesm/icaux_cesm.F ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../cesm/icaux_cesm.F
icaux_channel.o: ../phy/icaux_channel.F ../phy/common_blocks.h ../phy/common_forc.h mod_xc.o mod_ben02.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/icaux_channel.F
ictsz_channel.o: ../phy/ictsz_channel.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/ictsz_channel.F
ictsz_file.o: ../phy/ictsz_file.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/ictsz_file.F
idarlx.o: ../phy/idarlx.F ../phy/common_blocks.h ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/idarlx.F
idlage_init.o: ../idlage/idlage_init.F ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../trc	../idlage/idlage_init.F
idlage_step.o: ../idlage/idlage_step.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h ../phy/common_clndr.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../idlage/idlage_step.F
inicon.o: ../phy/inicon.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_seaice.h mod_xc.o mod_swtfrz.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/inicon.F
inifrc.o: ../phy/inifrc.F ../phy/common_blocks.h mod_xc.o mod_ben02.o mod_channel.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/inifrc.F
inifrc_cesm.o: ../cesm/inifrc_cesm.F ../phy/common_forc.h ../cesm/common_cesm.h mod_xc.o mod_ben02.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../cesm	../cesm/inifrc_cesm.F
inigeo.o: ../phy/inigeo.F ../phy/common_blocks.h ../phy/common_geo.h mod_xc.o mod_nctools.o mod_dia.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/inigeo.F
iniphy.o: ../phy/iniphy.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/iniphy.F
initime.o: ../phy/initime.F ../phy/timing.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/initime.F
initke.o: ../phy/initke.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_tke.h ../phy/common_tke.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/initke.F
initrc.o: ../trc/initrc.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../trc/initrc.F
inivar.o: ../phy/inivar.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_seaice.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/inivar.F
intp1d.o: ../phy/intp1d.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/intp1d.F
isinf.o: ../phy/isinf.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/isinf.F
isnan.o: ../phy/isnan.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/isnan.F
mod_ben02.o: ../ben02/mod_ben02.F ../phy/common_blocks.h ../phy/common_geo.h ../phy/common_clndr.h ../phy/common_forc.h mod_xc.o mod_ben02func.o mod_bulktf.o mod_swtfrz.o mod_thdysi.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../ben02/mod_ben02.F
mod_ben02func.o: ../ben02/mod_ben02func.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../ben02/mod_ben02func.F
mod_bulktf.o: ../ben02/mod_bulktf.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../ben02/mod_bulktf.F
mod_channel.o: ../phy/mod_channel.F ../phy/common_forc.h ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/mod_channel.F
mod_dia.o: ../phy/mod_dia.F ../phy/common_clndr.h ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_seaice.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o mod_nctools.o mod_thdysi.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/mod_dia.F
mod_nctools.o: ../phy/mod_nctools.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/mod_nctools.F
mod_swabs.o: ../phy/mod_swabs.F ../phy/common_blocks.h ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/mod_swabs.F
mod_swtfrz.o: ../drivers/nocoupler/mod_swtfrz.F eosfun.o 
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/nocoupler/mod_swtfrz.F
mod_thdysi.o: ../ben02/mod_thdysi.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../ben02/mod_thdysi.F
mod_xc.o: ../phy/mod_xc.F ../phy/mod_xc_mp.inc ../phy/mod_xc_sm.inc dimensions.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/mod_xc.F
momtum.o: ../phy/momtum.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/momtum.F
mxlayr.o: ../phy/mxlayr.F ../phy/common_blocks.h  ../trc/param_trc.h ../trc/common_trc.h ../phy/param_tke.h mod_xc.o mod_swabs.o mod_nctools.o eosfun.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/mxlayr.F
ocntrc_init.o: ../trc/ocntrc_init.F ../trc/param_trc.h ../trc/common_trc.h ../phy/common_geo.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../trc -I../phy	../trc/ocntrc_init.F
pbtcor.o: ../phy/pbtcor.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/pbtcor.F
penint.o: ../phy/penint.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/penint.F
pgforc.o: ../phy/pgforc.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/pgforc.F
rdcsss.o: ../phy/rdcsss.F ../phy/common_blocks.h ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/rdcsss.F
rdlim.o: ../phy/rdlim.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h mod_xc.o mod_swabs.o mod_dia.o mod_ben02.o data_mct.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/rdlim.F
read_twedon.o: ../phy/read_twedon.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/read_twedon.F
remap_eitflx.o: ../phy/remap_eitflx.F ../trc/param_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../trc	../phy/remap_eitflx.F
remap_eitvel.o: ../phy/remap_eitvel.F ../trc/param_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../trc	../phy/remap_eitvel.F
restart_ocntrcrd.o: ../trc/restart_ocntrcrd.F ../phy/common_blocks.h ../phy/common_clndr.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o mod_dia.o mod_nctools.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../trc/restart_ocntrcrd.F
restart_ocntrcwt.o: ../trc/restart_ocntrcwt.F ../phy/common_blocks.h ../phy/common_clndr.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o mod_dia.o mod_nctools.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../trc/restart_ocntrcwt.F
restart_rd.o: ../phy/restart_rd.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h ../cesm/common_cesm.h ../phy/common_seaice.h ../trc/param_trc.h ../trc/common_trc.h ../phy/common_tke.h mod_xc.o mod_dia.o mod_ben02.o mod_thdysi.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../cesm -I../trc	../phy/restart_rd.F
restart_trcrd.o: ../trc/restart_trcrd.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../trc/restart_trcrd.F
restart_trcwt.o: ../trc/restart_trcwt.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../trc/restart_trcwt.F
restart_wt.o: ../phy/restart_wt.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h ../cesm/common_cesm.h ../phy/common_seaice.h ../trc/param_trc.h ../trc/common_trc.h ../phy/common_tke.h mod_xc.o mod_dia.o mod_ben02.o mod_thdysi.o        \
mod_nctools.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../cesm -I../trc	../phy/restart_wt.F
setlogunit.o: ../drivers/nocoupler/setlogunit.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../drivers/nocoupler/setlogunit.F
settemmin.o: ../phy/settemmin.F ../phy/common_blocks.h mod_xc.o eosfun.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/settemmin.F
sfcstr.o: ../phy/sfcstr.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/sfcstr.F
sfcstr_ben02.o: ../ben02/sfcstr_ben02.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_seaice.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../ben02/sfcstr_ben02.F
sfcstr_cesm.o: ../cesm/sfcstr_cesm.F ../phy/common_blocks.h ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../cesm/sfcstr_cesm.F
sfcstr_channel.o: ../phy/sfcstr_channel.F ../phy/common_blocks.h ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/sfcstr_channel.F
stdbud.o: ../phy/stdbud.F ../phy/common_blocks.h ../phy/param_bud.h ../phy/common_bud.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/stdbud.F
thermf.o: ../phy/thermf.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/thermf.F
thermf_ben02.o: ../ben02/thermf_ben02.F ../phy/common_blocks.h ../phy/common_seaice.h ../phy/common_forc.h ../phy/common_geo.h ../phy/common_clndr.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_tke.h mod_xc.o mod_swtfrz.o mod_swabs.o mod_ben02.o   \
mod_thdysi.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../ben02/thermf_ben02.F
thermf_cesm.o: ../cesm/thermf_cesm.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_geo.h ../phy/common_clndr.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_tke.h mod_xc.o mod_swtfrz.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../cesm/thermf_cesm.F
thermf_channel.o: ../phy/thermf_channel.F ../phy/common_blocks.h ../phy/common_seaice.h ../phy/common_forc.h ../phy/common_geo.h ../phy/common_clndr.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_tke.h mod_xc.o mod_ben02.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/thermf_channel.F
tmsmt1.o: ../phy/tmsmt1.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/tmsmt1.F
tmsmt2.o: ../phy/tmsmt2.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/tmsmt2.F
triint.o: ../phy/triint.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/triint.F
types.o: ../phy/types.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	../phy/types.F90
updtrc.o: ../trc/updtrc.F ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../trc	../trc/updtrc.F
wdiflx.o: ../phy/wdiflx.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h mod_xc.o mod_nctools.o mod_dia.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/wdiflx.F
wtime.o: ../phy/wtime.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/wtime.F
./param_tke.h: ../phy/param_tke.h
	cp ../phy/param_tke.h .
./mod_channel.F: ../phy/mod_channel.F
	cp ../phy/mod_channel.F .
./common_trc.h: ../trc/common_trc.h
	cp ../trc/common_trc.h .
./external_abort.F: ../drivers/nocoupler/external_abort.F
	cp ../drivers/nocoupler/external_abort.F .
./mod_dia.F: ../phy/mod_dia.F
	cp ../phy/mod_dia.F .
./difest.F: ../phy/difest.F
	cp ../phy/difest.F .
./isinf.F: ../phy/isinf.F
	cp ../phy/isinf.F .
./idlage_step.F: ../idlage/idlage_step.F
	cp ../idlage/idlage_step.F .
./advect.F: ../phy/advect.F
	cp ../phy/advect.F .
./common_forc.h: ../phy/common_forc.h
	cp ../phy/common_forc.h .
./blom_step.F: ../phy/blom_step.F
	cp ../phy/blom_step.F .
./geoenv_test.F: ../phy/geoenv_test.F
	cp ../phy/geoenv_test.F .
./getfrc.F: ../phy/getfrc.F
	cp ../phy/getfrc.F .
./geoenv_cesmextra.F: ../cesm/geoenv_cesmextra.F
	cp ../cesm/geoenv_cesmextra.F .
./eddtra_gm.F: ../phy/eddtra_gm.F
	cp ../phy/eddtra_gm.F .
./mod_swabs.F: ../phy/mod_swabs.F
	cp ../phy/mod_swabs.F .
./chksum.F: ../phy/chksum.F
	cp ../phy/chksum.F .
./restart_ocntrcrd.F: ../trc/restart_ocntrcrd.F
	cp ../trc/restart_ocntrcrd.F .
./thermf_cesm.F: ../cesm/thermf_cesm.F
	cp ../cesm/thermf_cesm.F .
./getfrc_cesm.F: ../cesm/getfrc_cesm.F
	cp ../cesm/getfrc_cesm.F .
./updtrc.F: ../trc/updtrc.F
	cp ../trc/updtrc.F .
./crc.c: ../phy/crc.c
	cp ../phy/crc.c .
./common_geo.h: ../phy/common_geo.h
	cp ../phy/common_geo.h .
./forfun_channel.F: ../phy/forfun_channel.F
	cp ../phy/forfun_channel.F .
./inifrc.F: ../phy/inifrc.F
	cp ../phy/inifrc.F .
./diapfl.F: ../phy/diapfl.F
	cp ../phy/diapfl.F .
./clndr.F: ../phy/clndr.F
	cp ../phy/clndr.F .
./ocntrc_init.F: ../trc/ocntrc_init.F
	cp ../trc/ocntrc_init.F .
./bigrid.F: ../phy/bigrid.F
	cp ../phy/bigrid.F .
./icaux_cesm.F: ../cesm/icaux_cesm.F
	cp ../cesm/icaux_cesm.F .
./wtime.F: ../phy/wtime.F
	cp ../phy/wtime.F .
./icaux_ben02.F: ../ben02/icaux_ben02.F
	cp ../ben02/icaux_ben02.F .
./restart_rd.F: ../phy/restart_rd.F
	cp ../phy/restart_rd.F .
./inicon.F: ../phy/inicon.F
	cp ../phy/inicon.F .
./initime.F: ../phy/initime.F
	cp ../phy/initime.F .
./isnan.F: ../phy/isnan.F
	cp ../phy/isnan.F .
./ictsz_file.F: ../phy/ictsz_file.F
	cp ../phy/ictsz_file.F .
./rdlim.F: ../phy/rdlim.F
	cp ../phy/rdlim.F .
./rdcsss.F: ../phy/rdcsss.F
	cp ../phy/rdcsss.F .
./diffus.F: ../phy/diffus.F
	cp ../phy/diffus.F .
./eosfun.F: ../phy/eosfun.F
	cp ../phy/eosfun.F .
./gettime.F: ../phy/gettime.F
	cp ../phy/gettime.F .
./cmnfld_bfsqf.F: ../phy/cmnfld_bfsqf.F
	cp ../phy/cmnfld_bfsqf.F .
./tmsmt2.F: ../phy/tmsmt2.F
	cp ../phy/tmsmt2.F .
./read_twedon.F: ../phy/read_twedon.F
	cp ../phy/read_twedon.F .
./mod_xc_sm.inc: ../phy/mod_xc_sm.inc
	cp ../phy/mod_xc_sm.inc .
./frcdat.F: ../phy/frcdat.F
	cp ../phy/frcdat.F .
./blom.F: ../drivers/nocoupler/blom.F
	cp ../drivers/nocoupler/blom.F .
./fill_global.F: ../phy/fill_global.F
	cp ../phy/fill_global.F .
./blom_init.F: ../phy/blom_init.F
	cp ../phy/blom_init.F .
./inifrc_cesm.F: ../cesm/inifrc_cesm.F
	cp ../cesm/inifrc_cesm.F .
./cmnfld_nslope.F: ../phy/cmnfld_nslope.F
	cp ../phy/cmnfld_nslope.F .
./inivar.F: ../phy/inivar.F
	cp ../phy/inivar.F .
./common_seaice.h: ../phy/common_seaice.h
	cp ../phy/common_seaice.h .
./thermf_ben02.F: ../ben02/thermf_ben02.F
	cp ../ben02/thermf_ben02.F .
./pgforc.F: ../phy/pgforc.F
	cp ../phy/pgforc.F .
./sfcstr_channel.F: ../phy/sfcstr_channel.F
	cp ../phy/sfcstr_channel.F .
./eddtra_intdif.F: ../phy/eddtra_intdif.F
	cp ../phy/eddtra_intdif.F .
./mod_swtfrz.F: ../drivers/nocoupler/mod_swtfrz.F
	cp ../drivers/nocoupler/mod_swtfrz.F .
./chkvar.F: ../phy/chkvar.F
	cp ../phy/chkvar.F .
./eddtra.F: ../phy/eddtra.F
	cp ../phy/eddtra.F .
./mod_ben02.F: ../ben02/mod_ben02.F
	cp ../ben02/mod_ben02.F .
./types.F90: ../phy/types.F90
	cp ../phy/types.F90 .
./mxlayr.F: ../phy/mxlayr.F
	cp ../phy/mxlayr.F .
./wdiflx.F: ../phy/wdiflx.F
	cp ../phy/wdiflx.F .
./triint.F: ../phy/triint.F
	cp ../phy/triint.F .
./getlpi.F: ../phy/getlpi.F
	cp ../phy/getlpi.F .
./sfcstr.F: ../phy/sfcstr.F
	cp ../phy/sfcstr.F .
./common_blocks.h: ../phy/common_blocks.h
	cp ../phy/common_blocks.h .
./blkdat.F: ../phy/blkdat.F
	cp ../phy/blkdat.F .
./settemmin.F: ../phy/settemmin.F
	cp ../phy/settemmin.F .
./icaux_channel.F: ../phy/icaux_channel.F
	cp ../phy/icaux_channel.F .
./mod_nctools.F: ../phy/mod_nctools.F
	cp ../phy/mod_nctools.F .
./param_trc.h: ../trc/param_trc.h
	cp ../trc/param_trc.h .
./restart_wt.F: ../phy/restart_wt.F
	cp ../phy/restart_wt.F .
./mod_xc_mp.inc: ../phy/mod_xc_mp.inc
	cp ../phy/mod_xc_mp.inc .
./thermf_channel.F: ../phy/thermf_channel.F
	cp ../phy/thermf_channel.F .
./iniphy.F: ../phy/iniphy.F
	cp ../phy/iniphy.F .
./penint.F: ../phy/penint.F
	cp ../phy/penint.F .
./tmsmt1.F: ../phy/tmsmt1.F
	cp ../phy/tmsmt1.F .
./thermf.F: ../phy/thermf.F
	cp ../phy/thermf.F .
./convec.F: ../phy/convec.F
	cp ../phy/convec.F .
./initke.F: ../phy/initke.F
	cp ../phy/initke.F .
./crc.h: ../phy/crc.h
	cp ../phy/crc.h .
./param_bud.h: ../phy/param_bud.h
	cp ../phy/param_bud.h .
./common_tke.h: ../phy/common_tke.h
	cp ../phy/common_tke.h .
./common_bud.h: ../phy/common_bud.h
	cp ../phy/common_bud.h .
./idlage_init.F: ../idlage/idlage_init.F
	cp ../idlage/idlage_init.F .
./mod_thdysi.F: ../ben02/mod_thdysi.F
	cp ../ben02/mod_thdysi.F .
./geoenv_channel.F: ../phy/geoenv_channel.F
	cp ../phy/geoenv_channel.F .
./timing.h: ../phy/timing.h
	cp ../phy/timing.h .
./remap_eitvel.F: ../phy/remap_eitvel.F
	cp ../phy/remap_eitvel.F .
./stdbud.F: ../phy/stdbud.F
	cp ../phy/stdbud.F .
./intp1d.F: ../phy/intp1d.F
	cp ../phy/intp1d.F .
./momtum.F: ../phy/momtum.F
	cp ../phy/momtum.F .
./data_mct.F: ../cesm/data_mct.F
	cp ../cesm/data_mct.F .
./geoenv_file.F: ../phy/geoenv_file.F
	cp ../phy/geoenv_file.F .
./pbtcor.F: ../phy/pbtcor.F
	cp ../phy/pbtcor.F .
./mod_bulktf.F: ../ben02/mod_bulktf.F
	cp ../ben02/mod_bulktf.F .
./cmnfld.F: ../phy/cmnfld.F
	cp ../phy/cmnfld.F .
./ictsz_channel.F: ../phy/ictsz_channel.F
	cp ../phy/ictsz_channel.F .
./setlogunit.F: ../drivers/nocoupler/setlogunit.F
	cp ../drivers/nocoupler/setlogunit.F .
./restart_trcwt.F: ../trc/restart_trcwt.F
	cp ../trc/restart_trcwt.F .
./initrc.F: ../trc/initrc.F
	cp ../trc/initrc.F .
./mod_xc.F: ../phy/mod_xc.F
	cp ../phy/mod_xc.F .
./inigeo.F: ../phy/inigeo.F
	cp ../phy/inigeo.F .
./remap_eitflx.F: ../phy/remap_eitflx.F
	cp ../phy/remap_eitflx.F .
./fwbbal.F: ../phy/fwbbal.F
	cp ../phy/fwbbal.F .
./hfharm.F: ../phy/hfharm.F
	cp ../phy/hfharm.F .
./restart_trcrd.F: ../trc/restart_trcrd.F
	cp ../trc/restart_trcrd.F .
./barotp.F: ../phy/barotp.F
	cp ../phy/barotp.F .
./sfcstr_ben02.F: ../ben02/sfcstr_ben02.F
	cp ../ben02/sfcstr_ben02.F .
./chksummsk.F: ../phy/chksummsk.F
	cp ../phy/chksummsk.F .
./dimensions.F: ../phy/dimensions.F
	cp ../phy/dimensions.F .
./idarlx.F: ../phy/idarlx.F
	cp ../phy/idarlx.F .
./common_cesm.h: ../cesm/common_cesm.h
	cp ../cesm/common_cesm.h .
./restart_ocntrcwt.F: ../trc/restart_ocntrcwt.F
	cp ../trc/restart_ocntrcwt.F .
./common_clndr.h: ../phy/common_clndr.h
	cp ../phy/common_clndr.h .
./mod_ben02func.F: ../ben02/mod_ben02func.F
	cp ../ben02/mod_ben02func.F .
./sfcstr_cesm.F: ../cesm/sfcstr_cesm.F
	cp ../cesm/sfcstr_cesm.F .
SRC = ../ben02/mod_bulktf.F ../phy/fwbbal.F ../cesm/sfcstr_cesm.F ../ben02/mod_thdysi.F ../phy/inicon.F ../phy/tmsmt1.F ../ben02/icaux_ben02.F ../drivers/nocoupler/setlogunit.F ../idlage/idlage_init.F ../phy/dimensions.F ../phy/penint.F                   \
../phy/mod_swabs.F ../phy/hfharm.F ../drivers/nocoupler/external_abort.F ../phy/wtime.F ../phy/momtum.F ../phy/intp1d.F ../phy/initime.F ../ben02/sfcstr_ben02.F ../phy/geoenv_channel.F ../phy/mod_channel.F ../phy/blom_init.F ../trc/restart_ocntrcwt.F     \
../phy/chksummsk.F ../phy/sfcstr.F ../phy/remap_eitflx.F ../cesm/geoenv_cesmextra.F ../phy/iniphy.F ../phy/mod_xc.F ../phy/gettime.F ../phy/thermf_channel.F ../phy/eddtra.F ../phy/bigrid.F ../phy/triint.F ../trc/restart_trcwt.F            \
../phy/initke.F ../phy/tmsmt2.F ../phy/isnan.F ../phy/getlpi.F ../phy/rdlim.F ../phy/rdcsss.F ../phy/mod_dia.F ../phy/remap_eitvel.F ../phy/ictsz_channel.F ../phy/diapfl.F ../phy/mxlayr.F ../phy/blkdat.F ../phy/stdbud.F ../phy/geoenv_file.F               \
../phy/clndr.F ../trc/restart_trcrd.F ../phy/getfrc.F ../phy/eddtra_intdif.F ../phy/barotp.F ../phy/cmnfld_nslope.F ../phy/frcdat.F ../phy/forfun_channel.F ../phy/eosfun.F ../phy/mod_nctools.F ../phy/icaux_channel.F ../phy/sfcstr_channel.F                \
../phy/convec.F ../phy/ictsz_file.F ../drivers/nocoupler/blom.F ../ben02/mod_ben02func.F ../phy/inigeo.F ../phy/pgforc.F ../cesm/getfrc_cesm.F ../drivers/nocoupler/mod_swtfrz.F ../cesm/data_mct.F ../cesm/inifrc_cesm.F ../phy/chkvar.F ../phy/pbtcor.F      \
../phy/fill_global.F ../phy/thermf.F ../phy/cmnfld_bfsqf.F ../phy/diffus.F ../trc/ocntrc_init.F ../phy/inivar.F ../phy/cmnfld.F ../trc/initrc.F ../phy/restart_wt.F ../phy/eddtra_gm.F ../phy/blom_step.F ../phy/read_twedon.F ../ben02/thermf_ben02.F         \
../trc/restart_ocntrcrd.F ../idlage/idlage_step.F ../phy/geoenv_test.F ../phy/idarlx.F ../phy/difest.F ../phy/settemmin.F ../phy/restart_rd.F ../phy/isinf.F ../cesm/thermf_cesm.F ../phy/crc.c ../trc/updtrc.F ../phy/wdiflx.F ../phy/advect.F                \
../phy/types.F90 ../phy/inifrc.F ../phy/chksum.F ../cesm/icaux_cesm.F ../ben02/mod_ben02.F ../phy/common_seaice.h ../phy/mod_xc_mp.inc ../phy/common_tke.h ../phy/common_clndr.h ../phy/common_bud.h ../phy/common_geo.h ../phy/timing.h                       \
../phy/common_blocks.h ../phy/param_tke.h ../phy/mod_xc_sm.inc ../trc/common_trc.h ../phy/crc.h ../phy/param_bud.h ../cesm/common_cesm.h ../trc/param_trc.h ../phy/common_forc.h
OBJ = mod_bulktf.o fwbbal.o sfcstr_cesm.o mod_thdysi.o inicon.o tmsmt1.o icaux_ben02.o setlogunit.o idlage_init.o dimensions.o penint.o mod_swabs.o hfharm.o external_abort.o wtime.o momtum.o intp1d.o initime.o sfcstr_ben02.o geoenv_channel.o              \
mod_channel.o blom_init.o restart_ocntrcwt.o chksummsk.o sfcstr.o remap_eitflx.o geoenv_cesmextra.o iniphy.o mod_xc.o gettime.o thermf_channel.o eddtra.o bigrid.o triint.o restart_trcwt.o initke.o tmsmt2.o isnan.o getlpi.o rdlim.o rdcsss.o       \
mod_dia.o remap_eitvel.o ictsz_channel.o diapfl.o mxlayr.o convec.o blkdat.o stdbud.o geoenv_file.o clndr.o restart_trcrd.o getfrc.o eddtra_intdif.o barotp.o cmnfld_nslope.o frcdat.o forfun_channel.o eosfun.o mod_nctools.o icaux_channel.o sfcstr_channel.o         \
ictsz_file.o blom.o mod_ben02func.o inigeo.o pgforc.o getfrc_cesm.o mod_swtfrz.o data_mct.o inifrc_cesm.o chkvar.o pbtcor.o fill_global.o thermf.o cmnfld_bfsqf.o diffus.o ocntrc_init.o inivar.o cmnfld.o initrc.o restart_wt.o eddtra_gm.o          \
blom_step.o read_twedon.o thermf_ben02.o restart_ocntrcrd.o idlage_step.o geoenv_test.o idarlx.o difest.o settemmin.o restart_rd.o isinf.o thermf_cesm.o crc.o updtrc.o wdiflx.o advect.o types.o inifrc.o chksum.o icaux_cesm.o mod_ben02.o
OFF = ../phy/param_tke.h ../phy/mod_channel.F ../trc/common_trc.h ../drivers/nocoupler/external_abort.F ../phy/mod_dia.F ../phy/difest.F ../phy/isinf.F ../idlage/idlage_step.F ../phy/advect.F ../phy/common_forc.h ../phy/blom_step.F ../phy/geoenv_test.F   \
../phy/getfrc.F ../cesm/geoenv_cesmextra.F ../phy/eddtra_gm.F ../phy/mod_swabs.F ../phy/chksum.F ../trc/restart_ocntrcrd.F ../cesm/thermf_cesm.F ../cesm/getfrc_cesm.F ../trc/updtrc.F ../phy/crc.c ../phy/common_geo.h ../phy/forfun_channel.F                \
../phy/inifrc.F ../phy/diapfl.F ../phy/clndr.F ../trc/ocntrc_init.F ../phy/bigrid.F ../cesm/icaux_cesm.F ../phy/wtime.F ../ben02/icaux_ben02.F ../phy/restart_rd.F ../phy/inicon.F ../phy/initime.F ../phy/isnan.F ../phy/ictsz_file.F ../phy/rdlim.F          \
../phy/rdcsss.F ../phy/diffus.F ../phy/gettime.F ../phy/cmnfld_bfsqf.F ../phy/tmsmt2.F ../phy/read_twedon.F ../phy/mod_xc_sm.inc ../phy/frcdat.F ../drivers/nocoupler/blom.F ../phy/fill_global.F ../phy/blom_init.F ../cesm/inifrc_cesm.F     \
../phy/cmnfld_nslope.F ../phy/inivar.F ../phy/common_seaice.h ../ben02/thermf_ben02.F ../phy/pgforc.F ../phy/sfcstr_channel.F ../phy/eddtra_intdif.F ../drivers/nocoupler/mod_swtfrz.F ../phy/chkvar.F ../phy/eddtra.F ../ben02/mod_ben02.F ../phy/types.F90   \
../phy/mxlayr.F ../phy/wdiflx.F ../phy/triint.F ../phy/getlpi.F ../phy/sfcstr.F ../phy/common_blocks.h ../phy/blkdat.F ../phy/settemmin.F ../phy/icaux_channel.F ../phy/mod_nctools.F ../trc/param_trc.h ../phy/restart_wt.F               \
../phy/mod_xc_mp.inc ../phy/thermf_channel.F ../phy/iniphy.F ../phy/penint.F ../phy/tmsmt1.F ../phy/thermf.F ../phy/convec.F ../phy/initke.F ../phy/crc.h ../phy/param_bud.h ../phy/common_tke.h ../phy/common_bud.h ../idlage/idlage_init.F ../phy/eosfun.F   \
../ben02/mod_thdysi.F ../phy/geoenv_channel.F ../phy/timing.h ../phy/remap_eitvel.F ../phy/stdbud.F ../phy/intp1d.F ../phy/momtum.F ../cesm/data_mct.F ../phy/geoenv_file.F ../phy/pbtcor.F ../ben02/mod_bulktf.F ../phy/cmnfld.F ../phy/ictsz_channel.F       \
../drivers/nocoupler/setlogunit.F ../trc/restart_trcwt.F ../trc/initrc.F ../phy/mod_xc.F ../phy/inigeo.F ../phy/remap_eitflx.F ../phy/fwbbal.F ../phy/hfharm.F ../trc/restart_trcrd.F ../phy/barotp.F ../ben02/sfcstr_ben02.F ../phy/chksummsk.F               \
../phy/dimensions.F ../phy/idarlx.F ../cesm/common_cesm.h ../trc/restart_ocntrcwt.F ../phy/common_clndr.h ../ben02/mod_ben02func.F ../cesm/sfcstr_cesm.F
clean: neat
	-rm -f .cppdefs $(OBJ) blom *.mod
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
blom: $(OBJ) 
	$(LD) $(OBJ) -o blom  $(LDFLAGS)
