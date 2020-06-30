# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include Makefile.config

CPPDEFS = -DMPI -DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE -DPNETCDF

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
convec.o: ../phy/convec.F ../phy/common_blocks.h ../trc/param_trc.h ../trc/common_trc.h mod_xc.o
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
eosdat.o: ../phy/eosdat.F ../phy/common_eos.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/eosdat.F
eosfun.o: ../phy/eosfun.F ../phy/common_blocks.h ../phy/common_eos.h mod_xc.o
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
micom.o: ../drivers/nocoupler/micom.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/nocoupler/micom.F
micom_init.o: ../phy/micom_init.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/param_bud.h ../phy/common_bud.h mod_xc.o mod_swabs.o mod_dia.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/micom_init.F
micom_step.o: ../phy/micom_step.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_bud.h ../phy/common_bud.h ../phy/timing.h mod_xc.o mod_swabs.o mod_dia.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I../trc	../phy/micom_step.F
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
mod_swtfrz.o: ../drivers/nocoupler/mod_swtfrz.F ../phy/common_eos.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/nocoupler/mod_swtfrz.F
mod_thdysi.o: ../ben02/mod_thdysi.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../ben02/mod_thdysi.F
mod_xc.o: ../phy/mod_xc.F ../phy/mod_xc_mp.inc ../phy/mod_xc_sm.inc dimensions.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/mod_xc.F
momtum.o: ../phy/momtum.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/momtum.F
mxlayr.o: ../phy/mxlayr.F ../phy/common_blocks.h ../phy/common_eos.h ../trc/param_trc.h ../trc/common_trc.h ../phy/param_tke.h mod_xc.o mod_swabs.o mod_nctools.o
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
settemmin.o: ../phy/settemmin.F ../phy/common_blocks.h ../phy/common_eos.h mod_xc.o
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
./mod_thdysi.F: ../ben02/mod_thdysi.F
	cp ../ben02/mod_thdysi.F .
./hfharm.F: ../phy/hfharm.F
	cp ../phy/hfharm.F .
./inigeo.F: ../phy/inigeo.F
	cp ../phy/inigeo.F .
./intp1d.F: ../phy/intp1d.F
	cp ../phy/intp1d.F .
./data_mct.F: ../cesm/data_mct.F
	cp ../cesm/data_mct.F .
./cmnfld_nslope.F: ../phy/cmnfld_nslope.F
	cp ../phy/cmnfld_nslope.F .
./common_seaice.h: ../phy/common_seaice.h
	cp ../phy/common_seaice.h .
./common_blocks.h: ../phy/common_blocks.h
	cp ../phy/common_blocks.h .
./momtum.F: ../phy/momtum.F
	cp ../phy/momtum.F .
./getfrc.F: ../phy/getfrc.F
	cp ../phy/getfrc.F .
./thermf.F: ../phy/thermf.F
	cp ../phy/thermf.F .
./settemmin.F: ../phy/settemmin.F
	cp ../phy/settemmin.F .
./initke.F: ../phy/initke.F
	cp ../phy/initke.F .
./diapfl.F: ../phy/diapfl.F
	cp ../phy/diapfl.F .
./ictsz_channel.F: ../phy/ictsz_channel.F
	cp ../phy/ictsz_channel.F .
./eddtra_intdif.F: ../phy/eddtra_intdif.F
	cp ../phy/eddtra_intdif.F .
./idlage_step.F: ../idlage/idlage_step.F
	cp ../idlage/idlage_step.F .
./idlage_init.F: ../idlage/idlage_init.F
	cp ../idlage/idlage_init.F .
./micom_step.F: ../phy/micom_step.F
	cp ../phy/micom_step.F .
./chkvar.F: ../phy/chkvar.F
	cp ../phy/chkvar.F .
./eddtra.F: ../phy/eddtra.F
	cp ../phy/eddtra.F .
./triint.F: ../phy/triint.F
	cp ../phy/triint.F .
./idarlx.F: ../phy/idarlx.F
	cp ../phy/idarlx.F .
./fill_global.F: ../phy/fill_global.F
	cp ../phy/fill_global.F .
./chksum.F: ../phy/chksum.F
	cp ../phy/chksum.F .
./mod_swtfrz.F: ../drivers/nocoupler/mod_swtfrz.F
	cp ../drivers/nocoupler/mod_swtfrz.F .
./geoenv_file.F: ../phy/geoenv_file.F
	cp ../phy/geoenv_file.F .
./wdiflx.F: ../phy/wdiflx.F
	cp ../phy/wdiflx.F .
./geoenv_test.F: ../phy/geoenv_test.F
	cp ../phy/geoenv_test.F .
./diffus.F: ../phy/diffus.F
	cp ../phy/diffus.F .
./mod_xc.F: ../phy/mod_xc.F
	cp ../phy/mod_xc.F .
./thermf_channel.F: ../phy/thermf_channel.F
	cp ../phy/thermf_channel.F .
./sfcstr.F: ../phy/sfcstr.F
	cp ../phy/sfcstr.F .
./mod_xc_mp.inc: ../phy/mod_xc_mp.inc
	cp ../phy/mod_xc_mp.inc .
./restart_rd.F: ../phy/restart_rd.F
	cp ../phy/restart_rd.F .
./initime.F: ../phy/initime.F
	cp ../phy/initime.F .
./blkdat.F: ../phy/blkdat.F
	cp ../phy/blkdat.F .
./convec.F: ../phy/convec.F
	cp ../phy/convec.F .
./eosfun.F: ../phy/eosfun.F
	cp ../phy/eosfun.F .
./common_eos.h: ../phy/common_eos.h
	cp ../phy/common_eos.h .
./restart_trcrd.F: ../trc/restart_trcrd.F
	cp ../trc/restart_trcrd.F .
./cmnfld_bfsqf.F: ../phy/cmnfld_bfsqf.F
	cp ../phy/cmnfld_bfsqf.F .
./mod_dia.F: ../phy/mod_dia.F
	cp ../phy/mod_dia.F .
./mod_ben02func.F: ../ben02/mod_ben02func.F
	cp ../ben02/mod_ben02func.F .
./geoenv_channel.F: ../phy/geoenv_channel.F
	cp ../phy/geoenv_channel.F .
./tmsmt1.F: ../phy/tmsmt1.F
	cp ../phy/tmsmt1.F .
./common_bud.h: ../phy/common_bud.h
	cp ../phy/common_bud.h .
./mod_swabs.F: ../phy/mod_swabs.F
	cp ../phy/mod_swabs.F .
./crc.h: ../phy/crc.h
	cp ../phy/crc.h .
./eosdat.F: ../phy/eosdat.F
	cp ../phy/eosdat.F .
./fwbbal.F: ../phy/fwbbal.F
	cp ../phy/fwbbal.F .
./barotp.F: ../phy/barotp.F
	cp ../phy/barotp.F .
./wtime.F: ../phy/wtime.F
	cp ../phy/wtime.F .
./iniphy.F: ../phy/iniphy.F
	cp ../phy/iniphy.F .
./inifrc.F: ../phy/inifrc.F
	cp ../phy/inifrc.F .
./inifrc_cesm.F: ../cesm/inifrc_cesm.F
	cp ../cesm/inifrc_cesm.F .
./ictsz_file.F: ../phy/ictsz_file.F
	cp ../phy/ictsz_file.F .
./icaux_cesm.F: ../cesm/icaux_cesm.F
	cp ../cesm/icaux_cesm.F .
./pbtcor.F: ../phy/pbtcor.F
	cp ../phy/pbtcor.F .
./isnan.F: ../phy/isnan.F
	cp ../phy/isnan.F .
./getfrc_cesm.F: ../cesm/getfrc_cesm.F
	cp ../cesm/getfrc_cesm.F .
./rdlim.F: ../phy/rdlim.F
	cp ../phy/rdlim.F .
./isinf.F: ../phy/isinf.F
	cp ../phy/isinf.F .
./micom.F: ../drivers/nocoupler/micom.F
	cp ../drivers/nocoupler/micom.F .
./advect.F: ../phy/advect.F
	cp ../phy/advect.F .
./timing.h: ../phy/timing.h
	cp ../phy/timing.h .
./restart_ocntrcwt.F: ../trc/restart_ocntrcwt.F
	cp ../trc/restart_ocntrcwt.F .
./param_trc.h: ../trc/param_trc.h
	cp ../trc/param_trc.h .
./ocntrc_init.F: ../trc/ocntrc_init.F
	cp ../trc/ocntrc_init.F .
./read_twedon.F: ../phy/read_twedon.F
	cp ../phy/read_twedon.F .
./initrc.F: ../trc/initrc.F
	cp ../trc/initrc.F .
./sfcstr_cesm.F: ../cesm/sfcstr_cesm.F
	cp ../cesm/sfcstr_cesm.F .
./updtrc.F: ../trc/updtrc.F
	cp ../trc/updtrc.F .
./micom_init.F: ../phy/micom_init.F
	cp ../phy/micom_init.F .
./setlogunit.F: ../drivers/nocoupler/setlogunit.F
	cp ../drivers/nocoupler/setlogunit.F .
./restart_wt.F: ../phy/restart_wt.F
	cp ../phy/restart_wt.F .
./mxlayr.F: ../phy/mxlayr.F
	cp ../phy/mxlayr.F .
./remap_eitvel.F: ../phy/remap_eitvel.F
	cp ../phy/remap_eitvel.F .
./param_bud.h: ../phy/param_bud.h
	cp ../phy/param_bud.h .
./pgforc.F: ../phy/pgforc.F
	cp ../phy/pgforc.F .
./eddtra_gm.F: ../phy/eddtra_gm.F
	cp ../phy/eddtra_gm.F .
./getlpi.F: ../phy/getlpi.F
	cp ../phy/getlpi.F .
./penint.F: ../phy/penint.F
	cp ../phy/penint.F .
./cmnfld.F: ../phy/cmnfld.F
	cp ../phy/cmnfld.F .
./restart_trcwt.F: ../trc/restart_trcwt.F
	cp ../trc/restart_trcwt.F .
./common_clndr.h: ../phy/common_clndr.h
	cp ../phy/common_clndr.h .
./mod_channel.F: ../phy/mod_channel.F
	cp ../phy/mod_channel.F .
./bigrid.F: ../phy/bigrid.F
	cp ../phy/bigrid.F .
./mod_ben02.F: ../ben02/mod_ben02.F
	cp ../ben02/mod_ben02.F .
./dimensions.F: ../phy/dimensions.F
	cp ../phy/dimensions.F .
./inicon.F: ../phy/inicon.F
	cp ../phy/inicon.F .
./sfcstr_ben02.F: ../ben02/sfcstr_ben02.F
	cp ../ben02/sfcstr_ben02.F .
./restart_ocntrcrd.F: ../trc/restart_ocntrcrd.F
	cp ../trc/restart_ocntrcrd.F .
./mod_bulktf.F: ../ben02/mod_bulktf.F
	cp ../ben02/mod_bulktf.F .
./mod_nctools.F: ../phy/mod_nctools.F
	cp ../phy/mod_nctools.F .
./rdcsss.F: ../phy/rdcsss.F
	cp ../phy/rdcsss.F .
./sfcstr_channel.F: ../phy/sfcstr_channel.F
	cp ../phy/sfcstr_channel.F .
./icaux_ben02.F: ../ben02/icaux_ben02.F
	cp ../ben02/icaux_ben02.F .
./types.F90: ../phy/types.F90
	cp ../phy/types.F90 .
./thermf_ben02.F: ../ben02/thermf_ben02.F
	cp ../ben02/thermf_ben02.F .
./external_abort.F: ../drivers/nocoupler/external_abort.F
	cp ../drivers/nocoupler/external_abort.F .
./tmsmt2.F: ../phy/tmsmt2.F
	cp ../phy/tmsmt2.F .
./remap_eitflx.F: ../phy/remap_eitflx.F
	cp ../phy/remap_eitflx.F .
./geoenv_cesmextra.F: ../cesm/geoenv_cesmextra.F
	cp ../cesm/geoenv_cesmextra.F .
./icaux_channel.F: ../phy/icaux_channel.F
	cp ../phy/icaux_channel.F .
./chksummsk.F: ../phy/chksummsk.F
	cp ../phy/chksummsk.F .
./thermf_cesm.F: ../cesm/thermf_cesm.F
	cp ../cesm/thermf_cesm.F .
./common_tke.h: ../phy/common_tke.h
	cp ../phy/common_tke.h .
./common_forc.h: ../phy/common_forc.h
	cp ../phy/common_forc.h .
./crc.c: ../phy/crc.c
	cp ../phy/crc.c .
./inivar.F: ../phy/inivar.F
	cp ../phy/inivar.F .
./common_trc.h: ../trc/common_trc.h
	cp ../trc/common_trc.h .
./difest.F: ../phy/difest.F
	cp ../phy/difest.F .
./mod_xc_sm.inc: ../phy/mod_xc_sm.inc
	cp ../phy/mod_xc_sm.inc .
./forfun_channel.F: ../phy/forfun_channel.F
	cp ../phy/forfun_channel.F .
./common_cesm.h: ../cesm/common_cesm.h
	cp ../cesm/common_cesm.h .
./clndr.F: ../phy/clndr.F
	cp ../phy/clndr.F .
./gettime.F: ../phy/gettime.F
	cp ../phy/gettime.F .
./common_geo.h: ../phy/common_geo.h
	cp ../phy/common_geo.h .
./frcdat.F: ../phy/frcdat.F
	cp ../phy/frcdat.F .
./param_tke.h: ../phy/param_tke.h
	cp ../phy/param_tke.h .
./stdbud.F: ../phy/stdbud.F
	cp ../phy/stdbud.F .
SRC = ../phy/mxlayr.F ../phy/diffus.F ../ben02/icaux_ben02.F ../phy/chksummsk.F ../phy/types.F90 ../phy/eddtra_intdif.F ../phy/inivar.F ../drivers/nocoupler/setlogunit.F ../trc/updtrc.F ../phy/restart_rd.F ../phy/settemmin.F ../phy/gettime.F              \
../phy/mod_xc.F ../phy/cmnfld.F ../phy/mod_swabs.F ../phy/blkdat.F ../trc/restart_trcwt.F ../drivers/nocoupler/micom.F ../phy/geoenv_file.F ../phy/mod_channel.F ../phy/rdlim.F ../idlage/idlage_init.F ../phy/getfrc.F ../phy/bigrid.F ../phy/isnan.F         \
../phy/pbtcor.F ../trc/restart_ocntrcwt.F ../drivers/nocoupler/external_abort.F ../phy/chksum.F ../ben02/thermf_ben02.F ../ben02/mod_ben02.F ../phy/fill_global.F ../phy/wdiflx.F ../trc/restart_trcrd.F ../phy/inigeo.F ../phy/advect.F                       \
../idlage/idlage_step.F ../phy/thermf.F ../phy/idarlx.F ../phy/initime.F ../trc/initrc.F ../phy/remap_eitflx.F ../phy/convec.F ../drivers/nocoupler/mod_swtfrz.F ../phy/frcdat.F ../phy/barotp.F ../phy/inifrc.F ../phy/wtime.F ../phy/getlpi.F                \
../phy/cmnfld_bfsqf.F ../cesm/sfcstr_cesm.F ../phy/icaux_channel.F ../phy/tmsmt2.F ../cesm/inifrc_cesm.F ../phy/mod_nctools.F ../phy/thermf_channel.F ../phy/crc.c ../phy/micom_init.F ../phy/forfun_channel.F ../phy/dimensions.F ../phy/stdbud.F             \
../trc/restart_ocntrcrd.F ../phy/ictsz_file.F ../phy/ictsz_channel.F ../phy/remap_eitvel.F ../ben02/mod_thdysi.F ../phy/restart_wt.F ../phy/eddtra_gm.F ../phy/cmnfld_nslope.F ../phy/iniphy.F ../phy/tmsmt1.F ../phy/pgforc.F ../phy/mod_dia.F                \
../phy/rdcsss.F ../phy/chkvar.F ../ben02/mod_bulktf.F ../trc/ocntrc_init.F ../phy/diapfl.F ../cesm/getfrc_cesm.F ../phy/geoenv_test.F ../phy/clndr.F ../phy/geoenv_channel.F ../cesm/geoenv_cesmextra.F ../phy/difest.F ../ben02/sfcstr_ben02.F                \
../phy/inicon.F ../cesm/data_mct.F ../phy/eosfun.F ../cesm/thermf_cesm.F ../phy/hfharm.F ../ben02/mod_ben02func.F ../phy/intp1d.F ../phy/penint.F ../phy/momtum.F ../cesm/icaux_cesm.F ../phy/eosdat.F ../phy/read_twedon.F ../phy/micom_step.F                \
../phy/initke.F ../phy/sfcstr_channel.F ../phy/sfcstr.F ../phy/triint.F ../phy/eddtra.F ../phy/isinf.F ../phy/fwbbal.F ../phy/common_seaice.h ../phy/timing.h ../phy/common_blocks.h ../phy/mod_xc_mp.inc ../trc/param_trc.h ../phy/common_eos.h               \
../phy/common_tke.h ../phy/common_forc.h ../phy/common_bud.h ../trc/common_trc.h ../phy/param_bud.h ../phy/crc.h ../phy/mod_xc_sm.inc ../cesm/common_cesm.h ../phy/common_geo.h ../phy/param_tke.h ../phy/common_clndr.h
OBJ = mxlayr.o diffus.o icaux_ben02.o chksummsk.o types.o eddtra_intdif.o inivar.o setlogunit.o updtrc.o restart_rd.o settemmin.o gettime.o mod_xc.o cmnfld.o mod_swabs.o blkdat.o restart_trcwt.o micom.o geoenv_file.o mod_channel.o rdlim.o idlage_init.o   \
getfrc.o bigrid.o isnan.o pbtcor.o restart_ocntrcwt.o external_abort.o chksum.o thermf_ben02.o mod_ben02.o fill_global.o wdiflx.o restart_trcrd.o inigeo.o advect.o idlage_step.o thermf.o idarlx.o initime.o initrc.o remap_eitflx.o convec.o mod_swtfrz.o    \
frcdat.o barotp.o inifrc.o wtime.o getlpi.o cmnfld_bfsqf.o sfcstr_cesm.o icaux_channel.o tmsmt2.o inifrc_cesm.o mod_nctools.o thermf_channel.o crc.o micom_init.o forfun_channel.o dimensions.o stdbud.o restart_ocntrcrd.o ictsz_file.o ictsz_channel.o       \
remap_eitvel.o mod_thdysi.o restart_wt.o eddtra_gm.o cmnfld_nslope.o iniphy.o tmsmt1.o pgforc.o mod_dia.o rdcsss.o chkvar.o mod_bulktf.o ocntrc_init.o diapfl.o getfrc_cesm.o geoenv_test.o clndr.o geoenv_channel.o geoenv_cesmextra.o difest.o               \
sfcstr_ben02.o inicon.o data_mct.o eosfun.o thermf_cesm.o hfharm.o mod_ben02func.o intp1d.o penint.o momtum.o icaux_cesm.o eosdat.o read_twedon.o micom_step.o initke.o sfcstr_channel.o sfcstr.o triint.o eddtra.o isinf.o fwbbal.o
OFF = ../ben02/mod_thdysi.F ../phy/hfharm.F ../phy/inigeo.F ../phy/intp1d.F ../cesm/data_mct.F ../phy/cmnfld_nslope.F ../phy/common_seaice.h ../phy/common_blocks.h ../phy/momtum.F ../phy/getfrc.F ../phy/thermf.F ../phy/settemmin.F ../phy/initke.F         \
../phy/diapfl.F ../phy/ictsz_channel.F ../phy/eddtra_intdif.F ../idlage/idlage_step.F ../idlage/idlage_init.F ../phy/micom_step.F ../phy/chkvar.F ../phy/eddtra.F ../phy/triint.F ../phy/idarlx.F ../phy/fill_global.F ../phy/chksum.F                         \
../drivers/nocoupler/mod_swtfrz.F ../phy/geoenv_file.F ../phy/wdiflx.F ../phy/geoenv_test.F ../phy/diffus.F ../phy/mod_xc.F ../phy/thermf_channel.F ../phy/sfcstr.F ../phy/mod_xc_mp.inc ../phy/restart_rd.F ../phy/initime.F ../phy/blkdat.F ../phy/convec.F  \
../phy/eosfun.F ../phy/common_eos.h ../trc/restart_trcrd.F ../phy/cmnfld_bfsqf.F ../phy/mod_dia.F ../ben02/mod_ben02func.F ../phy/geoenv_channel.F ../phy/tmsmt1.F ../phy/common_bud.h ../phy/mod_swabs.F ../phy/crc.h ../phy/eosdat.F ../phy/fwbbal.F         \
../phy/barotp.F ../phy/wtime.F ../phy/iniphy.F ../phy/inifrc.F ../cesm/inifrc_cesm.F ../phy/ictsz_file.F ../cesm/icaux_cesm.F ../phy/pbtcor.F ../phy/isnan.F ../cesm/getfrc_cesm.F ../phy/rdlim.F ../phy/isinf.F ../drivers/nocoupler/micom.F ../phy/advect.F  \
../phy/timing.h ../trc/restart_ocntrcwt.F ../trc/param_trc.h ../trc/ocntrc_init.F ../phy/read_twedon.F ../trc/initrc.F ../cesm/sfcstr_cesm.F ../trc/updtrc.F ../phy/micom_init.F ../drivers/nocoupler/setlogunit.F ../phy/restart_wt.F ../phy/mxlayr.F         \
../phy/remap_eitvel.F ../phy/param_bud.h ../phy/pgforc.F ../phy/eddtra_gm.F ../phy/getlpi.F ../phy/penint.F ../phy/cmnfld.F ../trc/restart_trcwt.F ../phy/common_clndr.h ../phy/mod_channel.F ../phy/bigrid.F ../ben02/mod_ben02.F ../phy/dimensions.F         \
../phy/inicon.F ../ben02/sfcstr_ben02.F ../trc/restart_ocntrcrd.F ../ben02/mod_bulktf.F ../phy/mod_nctools.F ../phy/rdcsss.F ../phy/sfcstr_channel.F ../ben02/icaux_ben02.F ../phy/types.F90 ../ben02/thermf_ben02.F ../drivers/nocoupler/external_abort.F     \
../phy/tmsmt2.F ../phy/remap_eitflx.F ../cesm/geoenv_cesmextra.F ../phy/icaux_channel.F ../phy/chksummsk.F ../cesm/thermf_cesm.F ../phy/common_tke.h ../phy/common_forc.h ../phy/crc.c ../phy/inivar.F ../trc/common_trc.h ../phy/difest.F                     \
../phy/mod_xc_sm.inc ../phy/forfun_channel.F ../cesm/common_cesm.h ../phy/clndr.F ../phy/gettime.F ../phy/common_geo.h ../phy/frcdat.F ../phy/param_tke.h ../phy/stdbud.F
clean: neat
	-rm -f .cppdefs $(OBJ) blom
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
