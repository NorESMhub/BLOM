# iHAMOCC restart remapping

The script `iHAMOCC_restart_remapping.py` allows to remap iHAMOCC restart 
files between

- different grids (hybrid versus isopycnic) 
- different number of pressure level layers (currently 53 in isopycnic versus 56 in hybrid)
- different grid resolutions (e.g. interpolation from lower to higher resolution) and
- supports the different unit systems (CGS versus MKS - needed to allow for remapping between old/new blom restart files)

and can perform an inventory adjustment, if needed.

## Usage
To enable remapping, fill the `yml_template.yml` file according to your needs. Information on settings is provided there.
 `coarse_blom` and `coarse_hamocc` refers to the source restart files, while `fine_blom` refers to the target BLOM 
restart file/grid,  on which the iHAMOCC restart files will be interpolated.

Due to interpolation and regridding, the inventories of tracers can change slighly, which can be adjusted for 
automatically, when setting `AdjustInventory` to `True`.

Call:
```
python iHAMOCC_resart_remapping.py my_yml.yml
```

to run the restart remapping.


**NOTE** The remapping is tested with some grids (`tnx1` and `tnx2`). 0.5 degree grid support can be easily added, but isn't supported right now. 
