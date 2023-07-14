#!/usr/bin/env python3

"""BLOM namelist creator
"""

# Typically ignore this.
# pylint: disable=invalid-name

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import,wrong-import-position

import os, shutil, sys, glob, filecmp, imp, re

CIMEROOT = os.environ.get("CIMEROOT")
if CIMEROOT is None:
    raise SystemExit("ERROR: must set CIMEROOT environment variable")
sys.path.append(os.path.join(CIMEROOT, "CIME", "Tools"))

from standard_script_setup import *
from CIME.case import Case
from CIME.nmlgen import NamelistGenerator
from CIME.utils import expect
from CIME.utils import run_cmd_no_fail, expect
from CIME.utils import run_cmd
from CIME.buildnml import create_namelist_infile, parse_input

import glob, shutil
logger = logging.getLogger(__name__)

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def _create_namelists(case, confdir, infile, nmlgen):
####################################################################################

    """Write out the namelist for this component.

    Most arguments are the same as those for `NamelistGenerator`.
    The `confdir` argument is used to specify the directory  in which output files will be placed.
    """

    CASEROOT                  = case.get_value("CASEROOT")
    OCN_GRID                  = case.get_value("OCN_GRID")
    BLOM_VCOORD               = case.get_value("BLOM_VCOORD")
    BLOM_UNIT                 = case.get_value("BLOM_UNIT")
    DIN_LOC_ROOT              = case.get_value("DIN_LOC_ROOT")
    RUN_TYPE                  = case.get_value("RUN_TYPE")
    CONTINUE_RUN              = case.get_value("CONTINUE_RUN")
    CASEBUILD                 = case.get_value("CASEBUILD")
    CCSM_CO2_PPMV             = case.get_value("CCSM_CO2_PPMV")
    OCN_NCPL                  = case.get_value("OCN_NCPL")
    BLOM_COUPLING             = case.get_value("BLOM_COUPLING")
    RUNDIR                    = case.get_value("RUNDIR")
    BLOM_TRACER_MODULES       = case.get_value("BLOM_TRACER_MODULES")
    BLOM_RIVER_NUTRIENTS      = case.get_value("BLOM_RIVER_NUTRIENTS")
    BLOM_N_DEPOSITION         = case.get_value("BLOM_N_DEPOSITION")
    BLOM_NDEP_SCENARIO        = case.get_value("BLOM_NDEP_SCENARIO")
    HAMOCC_VSLS               = case.get_value("HAMOCC_VSLS")
    HAMOCC_CISO               = case.get_value("HAMOCC_CISO")
    HAMOCC_SEDSPINUP          = case.get_value("HAMOCC_SEDSPINUP")
    HAMOCC_SEDSPINUP_YR_START = case.get_value("HAMOCC_SEDSPINUP_YR_START")
    HAMOCC_SEDSPINUP_YR_END   = case.get_value("HAMOCC_SEDSPINUP_YR_END")
    HAMOCC_SEDSPINUP_NCYCLE   = case.get_value("HAMOCC_SEDSPINUP_NCYCLE")
    RUN_STARTDATE             = case.get_value("RUN_STARTDATE")
    PIO_TYPENAME_OCN          = case.get_value("PIO_TYPENAME_OCN")
    PIO_NETCDF_FORMAT_OCN     = case.get_value("PIO_NETCDF_FORMAT_OCN")
    NINST_OCN                 = case.get_value("NINST_OCN")
    TEST                      = case.get_value("TEST")

    config['ocn_comp'] = case.get_value("COMP_OCN")

    YEAR0   = `echo $RUN_STARTDATE | cut -c1-4 `
    MONTH0  = `echo $RUN_STARTDATE | cut -c6-7 `
    DAY0    = `echo $RUN_STARTDATE | cut -c9-10`

    #----------------------------------------------------
    # Initialize namelist defaults
    #----------------------------------------------------
    nmlgen.init_defaults(infile, config)

    #----------------------------------------------------
    # Write out namelist groups
    #----------------------------------------------------
    groups=['limit']

    namelist_file = os.path.join(confdir, "ocn_in")
    nmlgen.write_output_file(namelist_file, data_list_path, groups=groups, sorted_groups=False)

    #logger.debug("blom: grid is %s" %(hgrid))
    #logger.debug("blom: physics is %s "%phys)

###############################################################################
def buildnml(case, caseroot, compname):
###############################################################################
    """Build the blom namelist """

    # Build the component namelist
    if compname != "blom":
        raise AttributeError
    comp_root_dir_ocn = case.get_value("COMP_ROOT_DIR_OCN")
    srcroot = case.get_value("SRCROOT")
    rundir = case.get_value("RUNDIR")
    ninst = case.get_value("NINST_OCN")

    #----------------------------------------------------
    # Construct the namelist generator
    #----------------------------------------------------
    # determine directory for user modified namelist_definitions.xml and namelist_defaults.xml
    user_xml_dir = os.path.join(caseroot, "SourceMods", "src.blom")
    expect (os.path.isdir(user_xml_dir),
            "user_xml_dir %s does not exist " %user_xml_dir)

    # user definition *replaces* existing definition.
    namelist_xml_dir = os.path.join(comp_root_dir_ocn, "cime_config")
    definition_file = [os.path.join(namelist_xml_dir, "namelist_definition_blom.xml")]
    user_definition = os.path.join(user_xml_dir, "namelist_definition_blom.xml")
    if os.path.isfile(user_definition):
        definition_file = [user_definition]
    for file_ in definition_file:
        expect(os.path.isfile(file_), "Namelist XML file %s not found!" % file_)

    # Create the namelist generator object - independent of instance
    nmlgen = NamelistGenerator(case, definition_file)

    #----------------------------------------------------
    # Loop over instances
    #----------------------------------------------------
    for inst_counter in range(1, ninst+1):

        # determine instance string
        inst_string = ""
        if ninst > 1:
            inst_string = '_' + '%04d' % inst_counter

        # If multi-instance case does not have restart file, use
        # single-case restart for each instance
        rpointer = "rpointer.ice"
        if (os.path.isfile(os.path.join(rundir,rpointer)) and
            (not os.path.isfile(os.path.join(rundir,rpointer + inst_string)))):
            shutil.copy(os.path.join(rundir, rpointer),
                        os.path.join(rundir, rpointer + inst_string))

        inst_string_label = inst_string
        if not inst_string_label:
            inst_string_label = "\"\""

        # create namelist_infile using user_nl_file as input
        user_nl_file = os.path.join(caseroot, "user_nl_blom" + inst_string)
        expect(os.path.isfile(user_nl_file),
               "Missing required user_nl_file %s " %(user_nl_file))
        infile = os.path.join(confdir, "namelist_infile")
        create_namelist_infile(case, user_nl_file, infile)
        namelist_infile = [infile]

        # create namelist
        _create_namelists(case, confdir, namelist_infile, nmlgen)

        # copy namelist files to rundir
        if os.path.isdir(rundir):
            file1  = os.path.join(confdir, "ice_in")
            file2 = os.path.join(rundir, "ice_in")
            if inst_string:
                file2 += inst_string
            logger.debug("BLOM namelist copy: file1 %s file2 %s " %(file1, file2))
            shutil.copy2(file1, file2)

def _strip_comments(fh, token="!"):
    ''' strip anything after token in each line of fh '''
    for line in fh:
        s = line.split(token, 1)[0].strip()
        if s:
            yield s


###############################################################################
def _main_func():

    caseroot = parse_input(sys.argv)
    with Case(caseroot, read_only=False) as case:
        buildnml(case, caseroot, "blom")

if __name__ == "__main__":
    _main_func()
