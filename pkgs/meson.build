fs = import('fs') # import fs system module: https://mesonbuild.com/Fs-module.html
sources += files('CVMix-src/src/shared/cvmix_background.F90',
                 'CVMix-src/src/shared/cvmix_convection.F90', 
                 'CVMix-src/src/shared/cvmix_ddiff.F90',
                 'CVMix-src/src/shared/cvmix_kinds_and_types.F90', 
                 'CVMix-src/src/shared/cvmix_kpp.F90',
                 'CVMix-src/src/shared/cvmix_math.F90',
                 'CVMix-src/src/shared/cvmix_put_get.F90', 
                 'CVMix-src/src/shared/cvmix_shear.F90', 
                 'CVMix-src/src/shared/cvmix_tidal.F90', 
                 'CVMix-src/src/shared/cvmix_utils.F90')

sources += files('GSW-Fortran/modules/gsw_mod_kinds.f90',
                 'GSW-Fortran/modules/gsw_mod_teos10_constants.f90',
                 'GSW-Fortran/modules/gsw_mod_specvol_coefficients.f90',
                 'GSW-Fortran/modules/gsw_mod_saar_data.f90',
                 'GSW-Fortran/modules/gsw_mod_baltic_data.f90',
                 'GSW-Fortran/modules/gsw_mod_toolbox.f90',
                 'GSW-Fortran/modules/gsw_mod_error_functions.f90',
                 'GSW-Fortran/toolbox/gsw_enthalpy_sso_0.f90',
                 'GSW-Fortran/toolbox/gsw_specvol_sso_0.f90',
                 'GSW-Fortran/toolbox/gsw_entropy_part.f90',
                 'GSW-Fortran/toolbox/gsw_entropy_part_zerop.f90',
                 'GSW-Fortran/toolbox/gsw_gibbs_pt0_pt0.f90',
                 'GSW-Fortran/toolbox/gsw_add_barrier.f90',
                 'GSW-Fortran/toolbox/gsw_add_mean.f90',
                 'GSW-Fortran/toolbox/gsw_util_indx.f90',
                 'GSW-Fortran/toolbox/gsw_util_xinterp1.f90',
                 'GSW-Fortran/toolbox/gsw_saar.f90',
                 'GSW-Fortran/toolbox/gsw_sa_from_sp_baltic.f90',
                 'GSW-Fortran/toolbox/gsw_p_from_z.f90',
                 'GSW-Fortran/toolbox/gsw_sa_from_sp.f90',
                 'GSW-Fortran/toolbox/gsw_pt0_from_t.f90')

if get_option('ecosys')
   sources += files('M4AGO-sinking-scheme/src/mo_m4ago_core.f90',
                    'M4AGO-sinking-scheme/src/mo_m4ago_physics.f90',
       		    'M4AGO-sinking-scheme/src/mo_ihamocc4m4ago.f90')
   if fs.exists('M4AGO-sinking-scheme/src/mo_m4ago_kind.F90')
     sources += files('M4AGO-sinking-scheme/src/mo_m4ago_kind.F90')
   endif
   if fs.exists('M4AGO-sinking-scheme/src/mo_m4ago_control.f90')
     sources += files('M4AGO-sinking-scheme/src/mo_m4ago_control.f90')
   endif
   if fs.exists('M4AGO-sinking-scheme/src/mo_m4ago_params.f90')
     sources += files('M4AGO-sinking-scheme/src/mo_m4ago_params.f90')
   endif
   if fs.exists('M4AGO-sinking-scheme/src/mo_m4ago_types.f90')
     sources += files('M4AGO-sinking-scheme/src/mo_m4ago_types.f90')
   endif
   if fs.exists('M4AGO-sinking-scheme/src/mo_m4ago_HAMOCCinit.F90')
     sources += files('M4AGO-sinking-scheme/src/mo_m4ago_HAMOCCinit.F90')
   endif
   if fs.exists('M4AGO-sinking-scheme/src/mo_m4ago_HAMOCCPrimPart.F90')
     sources += files('M4AGO-sinking-scheme/src/mo_m4ago_HAMOCCPrimPart.F90')
   endif
endif
