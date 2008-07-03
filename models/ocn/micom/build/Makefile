# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include Makefile.config.exe

CPPDEFS = -DCCSMCOUPLED -DMPI -DNCEP -UICEDYN -UTRC -UATRC -UTRCDIA -UNEST -DCOMPDIA -DCOMPRES

.DEFAULT:
	-touch $@
all: micom
advect.o: ../phy/advect.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/advect.F
albw_eval.o: ../phy/albw_eval.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/albw_eval.F
asfdat.o: ../phy/asfdat.F ../phy/common_asf.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/asfdat.F
asflux.o: ../phy/asflux.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_clim.h ../phy/common_syn.h ../phy/common_asf.h ../phy/common_geo.h ../phy/common_clndr.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/asflux.F
atmdat.o: ../phy/atmdat.F ../phy/param_atm.h ../phy/common_atm.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/atmdat.F
barotp.o: ../phy/barotp.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/barotp.F
bigrid.o: ../phy/bigrid.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/bigrid.F
blkdat.o: ../phy/blkdat.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/blkdat.F
bulktf.o: ../phy/bulktf.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/bulktf.F
chksum.o: ../phy/chksum.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/chksum.F
chksummsk.o: ../phy/chksummsk.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/chksummsk.F
chwbud.o: ../phy/chwbud.F ../phy/common_blocks.h ../phy/common_geo.h ../phy/common_clim.h ../phy/common_forc.h ../phy/common_ice.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/chwbud.F
clndr.o: ../phy/clndr.F ../phy/common_clndr.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/clndr.F
compweights.o: ../phy/compweights.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/compweights.F
convec.o: ../phy/convec.F ../phy/common_blocks.h ../phy/common_dia.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/convec.F
crc.o: ../phy/crc.c ../phy/crc.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/crc.c
data_mct.o: ../drivers/ccsm_cpl7/data_mct.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../drivers/ccsm_cpl7/data_mct.F
diaacc.o: ../phy/diaacc.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_ice.h ../phy/common_dia.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/diaacc.F
diaini2d.o: ../phy/diaini2d.F ../phy/common_dia.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/diaini2d.F
diaini3d.o: ../phy/diaini3d.F ../phy/common_blocks.h ../phy/common_dia.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/diaini3d.F
diapfl.o: ../phy/diapfl.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/diapfl.F
diffus.o: ../phy/diffus.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/diffus.F
domain_mct.o: ../drivers/ccsm_cpl7/domain_mct.F ../phy/common_geo.h ../phy/common_blocks.h types.o mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/ccsm_cpl7/domain_mct.F
eosdat.o: ../phy/eosdat.F ../phy/common_eos.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/eosdat.F
eosfun.o: ../phy/eosfun.F ../phy/common_eos.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/eosfun.F
export_mct.o: ../drivers/ccsm_cpl7/export_mct.F ../phy/common_blocks.h ../phy/common_geo.h types.o mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/ccsm_cpl7/export_mct.F
fill_global.o: ../phy/fill_global.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/fill_global.F
forfun.o: ../phy/forfun.F ../phy/common_blocks.h ../phy/param_atm.h ../phy/common_atm.h ../phy/common_forc.h ../phy/common_geo.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/forfun.F
forfun_clim.o: ../phy/forfun_clim.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/param_atm.h ../phy/common_atm.h ../phy/common_clim.h ../phy/common_geo.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/forfun_clim.F
forfun_syn.o: ../phy/forfun_syn.F ../phy/common_blocks.h ../phy/param_atm.h ../phy/common_atm.h ../phy/common_atmsyn.h ../phy/common_syn.h ../phy/common_clndr.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/forfun_syn.F
frcdat.o: ../phy/frcdat.F ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/frcdat.F
geopar.o: ../phy/geopar.F ../phy/common_blocks.h ../phy/common_geo.h ../phy/common_forc.h ../phy/common_ice.h ../phy/common_asf.h /work/apps/netcdf-3.6.2-cnl-pgi/include/netcdf.inc mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy -I/work/apps/netcdf-3.6.2-cnl-pgi/include	../phy/geopar.F
getflux_clim.o: ../phy/getflux_clim.F ../phy/common_forc.h ../phy/common_clim.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/getflux_clim.F
getflux_syn.o: ../phy/getflux_syn.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/param_atm.h ../phy/common_syn.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/getflux_syn.F
getlpi.o: ../phy/getlpi.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/getlpi.F
gettime.o: ../phy/gettime.F ../phy/timing.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/gettime.F
hfharm.o: ../phy/hfharm.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/hfharm.F
icedat.o: ../phy/icedat.F ../phy/common_ice.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/icedat.F
import_mct.o: ../drivers/ccsm_cpl7/import_mct.F ../phy/common_blocks.h ../phy/common_geo.h ../phy/common_forc.h types.o mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/ccsm_cpl7/import_mct.F
inicon.o: ../phy/inicon.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_asf.h ../phy/common_ice.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/inicon.F
initime.o: ../phy/initime.F ../phy/timing.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/initime.F
intpol.o: ../phy/intpol.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/intpol.F
intpol2d.o: ../phy/intpol2d.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/intpol2d.F
locwrite.o: ../phy/locwrite.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/locwrite.F
micom.o: ../drivers/ccsm_cpl7/testtrg/micom.F90 ocn_comp_mct.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	../drivers/ccsm_cpl7/testtrg/micom.F90
micom_init.o: ../phy/micom_init.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/param_bud.h ../phy/common_bud.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/micom_init.F
micom_step.o: ../phy/micom_step.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_dia.h ../phy/common_forc.h ../phy/param_bud.h ../phy/common_bud.h ../phy/timing.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/micom_step.F
mod_xc.o: ../phy/mod_xc.F ../phy/dimensions.h ../phy/mod_xc_mp.inc data_mct.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/mod_xc.F
modeltime.o: ../drivers/ccsm_cpl7/modeltime.F ../phy/common_blocks.h ../phy/common_clndr.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/ccsm_cpl7/modeltime.F
momtum.o: ../phy/momtum.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/momtum.F
mxlayr.o: ../phy/mxlayr.F ../phy/common_blocks.h ../phy/common_eos.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/mxlayr.F
ncdiaout.o: ../phy/ncdiaout.F ../phy/common_clndr.h ../phy/common_blocks.h ../phy/common_dia.h ../phy/common_forc.h ../phy/common_ice.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/ncdiaout.F
nctools.o: ../phy/nctools.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/nctools.F
ocn_comp_mct.o: ../drivers/ccsm_cpl7/ocn_comp_mct.F90 types.o data_mct.o mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	../drivers/ccsm_cpl7/ocn_comp_mct.F90
pbtcor.o: ../phy/pbtcor.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/pbtcor.F
penint.o: ../phy/penint.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/penint.F
pgforc.o: ../phy/pgforc.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/pgforc.F
rdatm_syn.o: ../phy/rdatm_syn.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/param_atm.h ../phy/common_atm.h ../phy/common_atmsyn.h ../phy/common_syn.h ../phy/common_geo.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/rdatm_syn.F
rdlim.o: ../phy/rdlim.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_dia.h ../phy/common_forc.h ../phy/param_atm.h ../phy/common_atmsyn.h mod_xc.o data_mct.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/rdlim.F
read_atm_field.o: ../phy/read_atm_field.F /work/apps/netcdf-3.6.2-cnl-pgi/include/netcdf.inc mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I/work/apps/netcdf-3.6.2-cnl-pgi/include	../phy/read_atm_field.F
read_atm_info.o: ../phy/read_atm_info.F /work/apps/netcdf-3.6.2-cnl-pgi/include/netcdf.inc
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I/work/apps/netcdf-3.6.2-cnl-pgi/include	../phy/read_atm_info.F
read_atm_topo.o: ../phy/read_atm_topo.F /work/apps/netcdf-3.6.2-cnl-pgi/include/netcdf.inc
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I/work/apps/netcdf-3.6.2-cnl-pgi/include	../phy/read_atm_topo.F
remap.o: ../phy/remap.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/remap.F
restart_rd.o: ../phy/restart_rd.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h ../phy/common_asf.h ../phy/common_ice.h ../phy/common_dia.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/restart_rd.F
restart_wt.o: ../phy/restart_wt.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h ../phy/common_asf.h ../phy/common_ice.h ../phy/common_dia.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/restart_wt.F
settemmin.o: ../phy/settemmin.F ../phy/common_blocks.h ../phy/common_eos.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/settemmin.F
sfcstr.o: ../phy/sfcstr.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_ice.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/sfcstr.F
sfcstr_ccsm.o: ../drivers/ccsm_cpl7/sfcstr_ccsm.F ../phy/common_blocks.h ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/ccsm_cpl7/sfcstr_ccsm.F
smooth.o: ../phy/smooth.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/smooth.F
spherdist.o: ../phy/spherdist.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/spherdist.F
stdbud.o: ../phy/stdbud.F ../phy/common_blocks.h ../phy/param_bud.h ../phy/common_bud.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/stdbud.F
sum_sbuff.o: ../drivers/ccsm_cpl7/sum_sbuff.F ../phy/common_blocks.h ../phy/common_forc.h types.o mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/ccsm_cpl7/sum_sbuff.F
thermf.o: ../phy/thermf.F ../phy/common_blocks.h ../phy/common_ice.h ../phy/common_forc.h ../phy/common_geo.h ../phy/common_clndr.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/thermf.F
thermf_ccsm.o: ../drivers/ccsm_cpl7/thermf_ccsm.F ../phy/common_blocks.h ../phy/common_forc.h ../phy/common_geo.h ../phy/common_ice.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../drivers/ccsm_cpl7/thermf_ccsm.F
thermf_func.o: ../phy/thermf_func.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/thermf_func.F
tmsmt1.o: ../phy/tmsmt1.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/tmsmt1.F
tmsmt2.o: ../phy/tmsmt2.F ../phy/common_blocks.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/tmsmt2.F
triint.o: ../phy/triint.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/triint.F
tsroot.o: ../phy/tsroot.F ../phy/common_eos.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/tsroot.F
types.o: ../drivers/ccsm_cpl7/types.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	../drivers/ccsm_cpl7/types.F90
uvrotr2g.o: ../phy/uvrotr2g.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/uvrotr2g.F
wdiflx.o: ../phy/wdiflx.F ../phy/common_blocks.h ../phy/common_clndr.h ../phy/common_forc.h mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I../phy	../phy/wdiflx.F
wtime.o: ../phy/wtime.F
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	../phy/wtime.F
./hfharm.F: ../phy/hfharm.F
	cp ../phy/hfharm.F .
./diaini3d.F: ../phy/diaini3d.F
	cp ../phy/diaini3d.F .
./uvrotr2g.F: ../phy/uvrotr2g.F
	cp ../phy/uvrotr2g.F .
./chwbud.F: ../phy/chwbud.F
	cp ../phy/chwbud.F .
./common_blocks.h: ../phy/common_blocks.h
	cp ../phy/common_blocks.h .
./netcdf.inc: /work/apps/netcdf-3.6.2-cnl-pgi/include/netcdf.inc
	cp /work/apps/netcdf-3.6.2-cnl-pgi/include/netcdf.inc .
./param_atm.h: ../phy/param_atm.h
	cp ../phy/param_atm.h .
./spherdist.F: ../phy/spherdist.F
	cp ../phy/spherdist.F .
./diaacc.F: ../phy/diaacc.F
	cp ../phy/diaacc.F .
./momtum.F: ../phy/momtum.F
	cp ../phy/momtum.F .
./micom.F90: ../drivers/ccsm_cpl7/testtrg/micom.F90
	cp ../drivers/ccsm_cpl7/testtrg/micom.F90 .
./thermf.F: ../phy/thermf.F
	cp ../phy/thermf.F .
./settemmin.F: ../phy/settemmin.F
	cp ../phy/settemmin.F .
./intpol2d.F: ../phy/intpol2d.F
	cp ../phy/intpol2d.F .
./diapfl.F: ../phy/diapfl.F
	cp ../phy/diapfl.F .
./asflux.F: ../phy/asflux.F
	cp ../phy/asflux.F .
./forfun.F: ../phy/forfun.F
	cp ../phy/forfun.F .
./forfun_syn.F: ../phy/forfun_syn.F
	cp ../phy/forfun_syn.F .
./remap.F: ../phy/remap.F
	cp ../phy/remap.F .
./dimensions.h: ../phy/dimensions.h
	cp ../phy/dimensions.h .
./common_syn.h: ../phy/common_syn.h
	cp ../phy/common_syn.h .
./micom_step.F: ../phy/micom_step.F
	cp ../phy/micom_step.F .
./import_mct.F: ../drivers/ccsm_cpl7/import_mct.F
	cp ../drivers/ccsm_cpl7/import_mct.F .
./diaini2d.F: ../phy/diaini2d.F
	cp ../phy/diaini2d.F .
./triint.F: ../phy/triint.F
	cp ../phy/triint.F .
./data_mct.F: ../drivers/ccsm_cpl7/data_mct.F
	cp ../drivers/ccsm_cpl7/data_mct.F .
./bulktf.F: ../phy/bulktf.F
	cp ../phy/bulktf.F .
./fill_global.F: ../phy/fill_global.F
	cp ../phy/fill_global.F .
./chksum.F: ../phy/chksum.F
	cp ../phy/chksum.F .
./forfun_clim.F: ../phy/forfun_clim.F
	cp ../phy/forfun_clim.F .
./wdiflx.F: ../phy/wdiflx.F
	cp ../phy/wdiflx.F .
./diffus.F: ../phy/diffus.F
	cp ../phy/diffus.F .
./mod_xc.F: ../phy/mod_xc.F
	cp ../phy/mod_xc.F .
./sfcstr.F: ../phy/sfcstr.F
	cp ../phy/sfcstr.F .
./mod_xc_mp.inc: ../phy/mod_xc_mp.inc
	cp ../phy/mod_xc_mp.inc .
./common_atm.h: ../phy/common_atm.h
	cp ../phy/common_atm.h .
./restart_rd.F: ../phy/restart_rd.F
	cp ../phy/restart_rd.F .
./getflux_clim.F: ../phy/getflux_clim.F
	cp ../phy/getflux_clim.F .
./geopar.F: ../phy/geopar.F
	cp ../phy/geopar.F .
./common_dia.h: ../phy/common_dia.h
	cp ../phy/common_dia.h .
./types.F90: ../drivers/ccsm_cpl7/types.F90
	cp ../drivers/ccsm_cpl7/types.F90 .
./initime.F: ../phy/initime.F
	cp ../phy/initime.F .
./blkdat.F: ../phy/blkdat.F
	cp ../phy/blkdat.F .
./convec.F: ../phy/convec.F
	cp ../phy/convec.F .
./tsroot.F: ../phy/tsroot.F
	cp ../phy/tsroot.F .
./eosfun.F: ../phy/eosfun.F
	cp ../phy/eosfun.F .
./common_eos.h: ../phy/common_eos.h
	cp ../phy/common_eos.h .
./nctools.F: ../phy/nctools.F
	cp ../phy/nctools.F .
./tmsmt1.F: ../phy/tmsmt1.F
	cp ../phy/tmsmt1.F .
./common_bud.h: ../phy/common_bud.h
	cp ../phy/common_bud.h .
./eosdat.F: ../phy/eosdat.F
	cp ../phy/eosdat.F .
./crc.h: ../phy/crc.h
	cp ../phy/crc.h .
./rdatm_syn.F: ../phy/rdatm_syn.F
	cp ../phy/rdatm_syn.F .
./barotp.F: ../phy/barotp.F
	cp ../phy/barotp.F .
./wtime.F: ../phy/wtime.F
	cp ../phy/wtime.F .
./locwrite.F: ../phy/locwrite.F
	cp ../phy/locwrite.F .
./pbtcor.F: ../phy/pbtcor.F
	cp ../phy/pbtcor.F .
./atmdat.F: ../phy/atmdat.F
	cp ../phy/atmdat.F .
./common_asf.h: ../phy/common_asf.h
	cp ../phy/common_asf.h .
./rdlim.F: ../phy/rdlim.F
	cp ../phy/rdlim.F .
./advect.F: ../phy/advect.F
	cp ../phy/advect.F .
./timing.h: ../phy/timing.h
	cp ../phy/timing.h .
./compweights.F: ../phy/compweights.F
	cp ../phy/compweights.F .
./domain_mct.F: ../drivers/ccsm_cpl7/domain_mct.F
	cp ../drivers/ccsm_cpl7/domain_mct.F .
./sfcstr_ccsm.F: ../drivers/ccsm_cpl7/sfcstr_ccsm.F
	cp ../drivers/ccsm_cpl7/sfcstr_ccsm.F .
./export_mct.F: ../drivers/ccsm_cpl7/export_mct.F
	cp ../drivers/ccsm_cpl7/export_mct.F .
./common_ice.h: ../phy/common_ice.h
	cp ../phy/common_ice.h .
./micom_init.F: ../phy/micom_init.F
	cp ../phy/micom_init.F .
./restart_wt.F: ../phy/restart_wt.F
	cp ../phy/restart_wt.F .
./mxlayr.F: ../phy/mxlayr.F
	cp ../phy/mxlayr.F .
./param_bud.h: ../phy/param_bud.h
	cp ../phy/param_bud.h .
./pgforc.F: ../phy/pgforc.F
	cp ../phy/pgforc.F .
./modeltime.F: ../drivers/ccsm_cpl7/modeltime.F
	cp ../drivers/ccsm_cpl7/modeltime.F .
./getlpi.F: ../phy/getlpi.F
	cp ../phy/getlpi.F .
./penint.F: ../phy/penint.F
	cp ../phy/penint.F .
./sum_sbuff.F: ../drivers/ccsm_cpl7/sum_sbuff.F
	cp ../drivers/ccsm_cpl7/sum_sbuff.F .
./common_clndr.h: ../phy/common_clndr.h
	cp ../phy/common_clndr.h .
./read_atm_info.F: ../phy/read_atm_info.F
	cp ../phy/read_atm_info.F .
./bigrid.F: ../phy/bigrid.F
	cp ../phy/bigrid.F .
./inicon.F: ../phy/inicon.F
	cp ../phy/inicon.F .
./intpol.F: ../phy/intpol.F
	cp ../phy/intpol.F .
./albw_eval.F: ../phy/albw_eval.F
	cp ../phy/albw_eval.F .
./read_atm_topo.F: ../phy/read_atm_topo.F
	cp ../phy/read_atm_topo.F .
./thermf_func.F: ../phy/thermf_func.F
	cp ../phy/thermf_func.F .
./tmsmt2.F: ../phy/tmsmt2.F
	cp ../phy/tmsmt2.F .
./thermf_ccsm.F: ../drivers/ccsm_cpl7/thermf_ccsm.F
	cp ../drivers/ccsm_cpl7/thermf_ccsm.F .
./chksummsk.F: ../phy/chksummsk.F
	cp ../phy/chksummsk.F .
./crc.c: ../phy/crc.c
	cp ../phy/crc.c .
./common_forc.h: ../phy/common_forc.h
	cp ../phy/common_forc.h .
./ncdiaout.F: ../phy/ncdiaout.F
	cp ../phy/ncdiaout.F .
./common_clim.h: ../phy/common_clim.h
	cp ../phy/common_clim.h .
./icedat.F: ../phy/icedat.F
	cp ../phy/icedat.F .
./ocn_comp_mct.F90: ../drivers/ccsm_cpl7/ocn_comp_mct.F90
	cp ../drivers/ccsm_cpl7/ocn_comp_mct.F90 .
./clndr.F: ../phy/clndr.F
	cp ../phy/clndr.F .
./gettime.F: ../phy/gettime.F
	cp ../phy/gettime.F .
./common_geo.h: ../phy/common_geo.h
	cp ../phy/common_geo.h .
./frcdat.F: ../phy/frcdat.F
	cp ../phy/frcdat.F .
./getflux_syn.F: ../phy/getflux_syn.F
	cp ../phy/getflux_syn.F .
./smooth.F: ../phy/smooth.F
	cp ../phy/smooth.F .
./asfdat.F: ../phy/asfdat.F
	cp ../phy/asfdat.F .
./common_atmsyn.h: ../phy/common_atmsyn.h
	cp ../phy/common_atmsyn.h .
./read_atm_field.F: ../phy/read_atm_field.F
	cp ../phy/read_atm_field.F .
./stdbud.F: ../phy/stdbud.F
	cp ../phy/stdbud.F .
SRC = ../phy/read_atm_topo.F ../phy/diffus.F ../phy/mxlayr.F ../drivers/ccsm_cpl7/thermf_ccsm.F ../phy/chksummsk.F ../drivers/ccsm_cpl7/types.F90 ../phy/ncdiaout.F ../phy/restart_rd.F ../phy/settemmin.F ../phy/getflux_syn.F ../phy/gettime.F               \
../drivers/ccsm_cpl7/import_mct.F ../phy/mod_xc.F ../phy/geopar.F ../phy/blkdat.F ../drivers/ccsm_cpl7/testtrg/micom.F90 ../phy/rdatm_syn.F ../phy/rdlim.F ../phy/bigrid.F ../phy/nctools.F ../phy/pbtcor.F ../phy/asfdat.F ../phy/chksum.F                    \
../phy/fill_global.F ../phy/wdiflx.F ../phy/advect.F ../phy/spherdist.F ../drivers/ccsm_cpl7/modeltime.F ../phy/thermf.F ../phy/initime.F ../phy/convec.F ../phy/frcdat.F ../phy/forfun_clim.F ../phy/barotp.F ../drivers/ccsm_cpl7/export_mct.F               \
../phy/forfun.F ../phy/wtime.F ../phy/getlpi.F ../phy/diaini2d.F ../phy/intpol2d.F ../phy/tmsmt2.F ../phy/uvrotr2g.F ../phy/albw_eval.F ../phy/diaini3d.F ../phy/remap.F ../phy/crc.c ../phy/micom_init.F ../phy/stdbud.F ../phy/forfun_syn.F ../phy/tsroot.F  \
../phy/read_atm_info.F ../phy/restart_wt.F ../phy/atmdat.F ../drivers/ccsm_cpl7/sfcstr_ccsm.F ../phy/diaacc.F ../phy/getflux_clim.F ../phy/tmsmt1.F ../phy/pgforc.F ../phy/read_atm_field.F ../phy/asflux.F ../phy/diapfl.F ../drivers/ccsm_cpl7/sum_sbuff.F   \
../phy/clndr.F ../phy/smooth.F ../drivers/ccsm_cpl7/ocn_comp_mct.F90 ../phy/bulktf.F ../phy/inicon.F ../phy/eosfun.F ../drivers/ccsm_cpl7/data_mct.F ../phy/hfharm.F ../phy/chwbud.F ../phy/thermf_func.F ../drivers/ccsm_cpl7/domain_mct.F ../phy/penint.F    \
../phy/momtum.F ../phy/icedat.F ../phy/locwrite.F ../phy/intpol.F ../phy/eosdat.F ../phy/micom_step.F ../phy/triint.F ../phy/sfcstr.F ../phy/compweights.F ../phy/timing.h ../phy/common_blocks.h ../phy/mod_xc_mp.inc                                         \
/work/apps/netcdf-3.6.2-cnl-pgi/include/netcdf.inc ../phy/common_atm.h ../phy/param_atm.h ../phy/common_dia.h ../phy/common_eos.h ../phy/common_ice.h ../phy/common_forc.h ../phy/common_bud.h ../phy/common_clim.h ../phy/dimensions.h ../phy/param_bud.h     \
../phy/common_syn.h ../phy/crc.h ../phy/common_geo.h ../phy/common_clndr.h ../phy/common_atmsyn.h ../phy/common_asf.h
OBJ = read_atm_topo.o diffus.o mxlayr.o thermf_ccsm.o chksummsk.o types.o ncdiaout.o restart_rd.o settemmin.o getflux_syn.o gettime.o import_mct.o mod_xc.o geopar.o blkdat.o micom.o rdatm_syn.o rdlim.o bigrid.o nctools.o pbtcor.o asfdat.o chksum.o        \
fill_global.o wdiflx.o advect.o spherdist.o modeltime.o thermf.o initime.o convec.o frcdat.o forfun_clim.o barotp.o export_mct.o forfun.o wtime.o getlpi.o diaini2d.o intpol2d.o tmsmt2.o uvrotr2g.o albw_eval.o diaini3d.o remap.o crc.o micom_init.o         \
stdbud.o forfun_syn.o tsroot.o read_atm_info.o restart_wt.o atmdat.o sfcstr_ccsm.o diaacc.o getflux_clim.o tmsmt1.o pgforc.o read_atm_field.o asflux.o diapfl.o sum_sbuff.o clndr.o smooth.o ocn_comp_mct.o bulktf.o inicon.o eosfun.o data_mct.o hfharm.o     \
chwbud.o thermf_func.o domain_mct.o penint.o momtum.o icedat.o locwrite.o intpol.o eosdat.o micom_step.o triint.o sfcstr.o compweights.o
OFF = ../phy/hfharm.F ../phy/diaini3d.F ../phy/uvrotr2g.F ../phy/chwbud.F ../phy/common_blocks.h /work/apps/netcdf-3.6.2-cnl-pgi/include/netcdf.inc ../phy/param_atm.h ../phy/spherdist.F ../phy/diaacc.F ../phy/momtum.F                                      \
../drivers/ccsm_cpl7/testtrg/micom.F90 ../phy/thermf.F ../phy/settemmin.F ../phy/intpol2d.F ../phy/diapfl.F ../phy/asflux.F ../phy/forfun.F ../phy/forfun_syn.F ../phy/remap.F ../phy/dimensions.h ../phy/common_syn.h ../phy/micom_step.F                     \
../drivers/ccsm_cpl7/import_mct.F ../phy/diaini2d.F ../phy/triint.F ../drivers/ccsm_cpl7/data_mct.F ../phy/bulktf.F ../phy/fill_global.F ../phy/chksum.F ../phy/forfun_clim.F ../phy/wdiflx.F ../phy/diffus.F ../phy/mod_xc.F ../phy/sfcstr.F                  \
../phy/mod_xc_mp.inc ../phy/common_atm.h ../phy/restart_rd.F ../phy/getflux_clim.F ../phy/geopar.F ../phy/common_dia.h ../drivers/ccsm_cpl7/types.F90 ../phy/initime.F ../phy/blkdat.F ../phy/convec.F ../phy/tsroot.F ../phy/eosfun.F ../phy/common_eos.h     \
../phy/nctools.F ../phy/tmsmt1.F ../phy/common_bud.h ../phy/eosdat.F ../phy/crc.h ../phy/rdatm_syn.F ../phy/barotp.F ../phy/wtime.F ../phy/locwrite.F ../phy/pbtcor.F ../phy/atmdat.F ../phy/common_asf.h ../phy/rdlim.F ../phy/advect.F ../phy/timing.h       \
../phy/compweights.F ../drivers/ccsm_cpl7/domain_mct.F ../drivers/ccsm_cpl7/sfcstr_ccsm.F ../drivers/ccsm_cpl7/export_mct.F ../phy/common_ice.h ../phy/micom_init.F ../phy/restart_wt.F ../phy/mxlayr.F ../phy/param_bud.h ../phy/pgforc.F                     \
../drivers/ccsm_cpl7/modeltime.F ../phy/getlpi.F ../phy/penint.F ../drivers/ccsm_cpl7/sum_sbuff.F ../phy/common_clndr.h ../phy/read_atm_info.F ../phy/bigrid.F ../phy/inicon.F ../phy/intpol.F ../phy/albw_eval.F ../phy/read_atm_topo.F ../phy/thermf_func.F  \
../phy/tmsmt2.F ../drivers/ccsm_cpl7/thermf_ccsm.F ../phy/chksummsk.F ../phy/crc.c ../phy/common_forc.h ../phy/ncdiaout.F ../phy/common_clim.h ../phy/icedat.F ../drivers/ccsm_cpl7/ocn_comp_mct.F90 ../phy/clndr.F ../phy/gettime.F ../phy/common_geo.h       \
../phy/frcdat.F ../phy/getflux_syn.F ../phy/smooth.F ../phy/asfdat.F ../phy/common_atmsyn.h ../phy/read_atm_field.F ../phy/stdbud.F
clean: neat
	-rm -f .cppdefs $(OBJ) micom
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
micom: $(OBJ) 
	$(LD) $(OBJ) -o micom  $(LDFLAGS)
