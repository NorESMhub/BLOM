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

if get_option('ecosys')
   sources += files('M4AGO-sinking-scheme/src/mo_m4ago_core.f90',
                    'M4AGO-sinking-scheme/src/mo_m4ago_physics.f90',
       		    'M4AGO-sinking-scheme/src/mo_ihamocc4m4ago.f90')
endif
