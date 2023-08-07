"""
Wrapper-class for the ParamGen
CIME tool, and associated methods
needed to generated the "ocn_in"
Fortran namelist file.

To run doctests on this file: python -m doctest ocn_in_paramgen.py
"""

#----------------------------------------
# Import generic python libraries/modules
#----------------------------------------

import os
import os.path
import sys
import re
from collections import OrderedDict

#----------------
# Import ParamGen
#----------------

_CIME_ROOT = os.environ.get("CIMEROOT")
if _CIME_ROOT is None:
    raise SystemExit("ERROR: must set CIMEROOT environment variable")
sys.path.append(os.path.join(_CIME_ROOT, "CIME", "Tools"))

#_PARAMGEN_ROOT = os.path.join(_CIME_ROOT, "CIME", "ParamGen")
# if not os.path.exists(_PARAMGEN_ROOT):
#     _EMSG = f"ERROR: Cannot find '{_PARAMGEN_ROOT}' directory.  Did you run checkout_externals?"
#     raise SystemExit(_EMSG)
#End if
#sys.path.append(_PARAMGEN_ROOT)
#pylint: disable=wrong-import-position
from paramgen import ParamGen
#pylint: enable=wrong-import-position

#Set of single and double quotes used by "_check_string_quotes" function:
_QUOTE_SET = {"'", '"'}

#Regular expression used by "remove_user_nl_comment" function:
_QUOTE_REGEX = re.compile(r"\".*?\"|'.*?'")

#Regular expression used by the "write" and "append_user_nl_file"
#methods to determine if the variable is an array, and what
#the dimensions of the array are:
_ARRAY_TYPE_REGEX = re.compile(r"[(][ ]*([0-9 ,]+)[ ]*[)]")

#Regular expression used to determine array indices in
#"find_arr_indices" function:
_ARR_INDEX_REGEX = re.compile(r"\((.+?)\)")

################################################################

class OcnInParamGenError(ValueError):
    """Class used to handle ocn_in ParamGen errors
    (e.g., log user errors without backtrace)"""

################################################################
#HELPER FUNCTIONS
################################################################

def _is_nml_logical_true(varname, var_val):

    """
    Checks if a "logical" XML namelist value is true or
    false.
    ----------
    varname -> The name of the variable being checked
    var_val -> The value of the variable being checked

    Returns a boolean that matches the value of the
    input logical.

    doctests:

    1. Check that a True value returns true:
    >>> _is_nml_logical_true("test", True)
    True

    2.  Check that a "true" value returns true:
    >>> _is_nml_logical_true("test", "true")
    True

    3.  Check that a ".true." value returns true:
    >>> _is_nml_logical_true("test", ".true.")
    True

    4.  Check that a "1" value returns true:
    >>> _is_nml_logical_true("test", "1")
    True

    5.  Check that a 1 (integer) value returns true:
    >>> _is_nml_logical_true("test", 1)
    True

    6.  Check that a False value returns false:
    >>> _is_nml_logical_true("test", False)
    False

    7.  Check that a "FALSE" value returns false:
    >>> _is_nml_logical_true("test", "FALSE")
    False

    8.  Check that a ".False." value returns false:
    >>> _is_nml_logical_true("test", ".False.")
    False

    9.  Check that a "0" value returns false:
    >>> _is_nml_logical_true("test", "0")
    False

    10.  Check that a 0 (integer) value returns false:
    >>> _is_nml_logical_true("test", 0)
    False

    11.  Check that a bad string value returns the correct error:
    >>> _is_nml_logical_true("test", "this_wont_work") # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:...
    XML namelist logical variable, 'test', must have a value of true, false, 1, or 0, not 'this_wont_work'

    12.  Check that a bad integer value returns the correct error:
    >>> _is_nml_logical_true("test", 3) # doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    ocn_in_paramgen.OcnInParamGenError:...
    XML namelist logical variable, 'test', must have a value of true, false, 1, or 0, not 3

    13.  Check that an unsupported type returns an error:
    >>> _is_nml_logical_true("test", 13.03) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:...
    XML namelist variable 'test' must have a value that is either a boolean, string, or integer, not float.

    """

    if isinstance(var_val, bool):
        return var_val
    #End if
    if isinstance(var_val, str):
        if var_val.lower() in {"true", ".true.", "1"}:
            return True
        #End if
        if var_val.lower() in {"false", ".false.", "0"}:
            return False
        #End if

        #Raise error if no match was found:
        emsg = f"\nXML namelist logical variable, '{varname}'"
        emsg += ", must have a value of true, false, 1, or 0, not"
        emsg += f" '{var_val}'"
        raise OcnInParamGenError(emsg)
    #End if

    if isinstance(var_val, int):
        if var_val == 1:
            return True
        #End if
        if var_val == 0:
            return False
        #End if

        #Raise error if no match was found:
        emsg = f"\nXML namelist logical variable, '{varname}'"
        emsg += ", must have a value of true, false, 1, or 0, not"
        emsg += f" {var_val}"
        raise OcnInParamGenError(emsg)
    #End if

    #Type is un-recognizeda, so raise an error:
    emsg = f"\nXML namelist variable '{varname}' must"
    emsg += " have a value that is either a boolean, string, or integer,"
    emsg += f" not {type(var_val).__name__}."
    raise OcnInParamGenError(emsg)

#####

def remove_user_nl_comment(user_string, comment_delim="!"):

    """
    Searches a one-line input string for a comment delimiter,
    and then returns the string with all text after the delimiter
    removed.
    ----------
    user_string   -> String that will be searched and processed for comments
    comment_delim -> Optional variable that sets the character type being used
                     as a comment delimiter. Defaults to the standard "!" fortran comment.

    Returns the input string, but with any commented text removed.

    doctests:

    1.  Check that a string with no comment delimiters returns full string:
    >>> remove_user_nl_comment("bananas")
    'bananas'

    2.  Check that a string with no comments outside quotes returns full string:
    >>> remove_user_nl_comment(" '!ban!anas!' ")
    " '!ban!anas!' "

    3.  Check that a string with no quotes but a comment returns string with no comment:
    >>> remove_user_nl_comment("bananas !But not apples")
    'bananas '

    4.  Check that a string with quotes and a comment returns string sans comment:
    >>> remove_user_nl_comment(" 'bananas' !But not apples")
    " 'bananas' "

    5.  Check that a string with a quoted comment and real comments returns proper string:
    >>> remove_user_nl_comment(" '!ba!na!nas!' !But not apples")
    " '!ba!na!nas!' "

    6.  Check that a string with a quoted comment and a real comment with multiple delimiters
        returns the proper string:
    >>> remove_user_nl_comment(" '!bananas!' !But not! apples!")
    " '!bananas!' "

    7.  Check that a string with a quoted comment and a commented quote returns the proper string:
    >>> remove_user_nl_comment(' "!bananas" !"But not apples" ')
    ' "!bananas" '

    8.  Check that a string with quotes inside quotes and multiple delimiters returns
        the proper string:
    >>> remove_user_nl_comment(''' "!bana'!'anas""other''fruit!" !But not '!Apples!' ''')
    ' "!bana\\'!\\'anas""other\\'\\'fruit!" '

    9.  Check that an array of strings returns the proper string:
    >>> remove_user_nl_comment(" 'bananas', 'apples', 'kiwis' ")
    " 'bananas', 'apples', 'kiwis' "

    10. Check that an array of strings with a comment returns the proper string:
    >>> remove_user_nl_comment(" 'bananas', 'apples', 'kiwis', !, and coconuts")
    " 'bananas', 'apples', 'kiwis', "

    11. Check that an array of of strings with comment delimiters and an actual comment
        returns the proper string:
    >>> remove_user_nl_comment(' , "!bananas", "app!les", "kiwis!", !And "Coconuts"!')
    ' , "!bananas", "app!les", "kiwis!", '

    12.  Check that a line with no comments or strings returns the proper string:
    >>> remove_user_nl_comment('5')
    '5'

    13.  Check that a line with a comment  but no internal strings returns the proper string:
    >>> remove_user_nl_comment(' .true. !And not .false.')
    ' .true. '

    14.  Check that an array of values with no comment returns the proper string:
    >>> remove_user_nl_comment('13.0d0, 15.0d0, 1100.35d0')
    '13.0d0, 15.0d0, 1100.35d0'

    15.  Check that an array of  values with a comment returns the proper string:
    >>> remove_user_nl_comment('13.0d0,! 15.0d0, 1100.35d0')
    '13.0d0,'

    16.  Check that a line that only contains a comment returns an empty string:
    >>> remove_user_nl_comment('! bananas and 13.0d0 5 .true. !@$#%*?')
    ''

    17.  Check that a line with an alternative comment delimiter returns the proper string:
    >>> remove_user_nl_comment('bananas #and 13.0d0 5 .true. !@$#%*?', comment_delim='#')
    'bananas '

    18. Check that some more unusual strings are handled correctly
    >>> remove_user_nl_comment("'Isn''t it a nice day'")
    "'Isn''t it a nice day'"
    >>> remove_user_nl_comment("'Isn''t it a nice day' !comment")
    "'Isn''t it a nice day' "
    >>> remove_user_nl_comment("'Isn!''!t it a nice! day'")
    "'Isn!''!t it a nice! day'"
    >>> remove_user_nl_comment("'Isn!''!t it a nice! day' ! comment")
    "'Isn!''!t it a nice! day' "
    >>> remove_user_nl_comment('''"This is 'one' string"''')
    '"This is \\'one\\' string"'
    >>> remove_user_nl_comment('''"This is 'one' string" !comment''')
    '"This is \\'one\\' string" '
    >>> remove_user_nl_comment("'This is \\"one\\" string'")
    '\\'This is "one" string\\''
    >>> remove_user_nl_comment("'This is \\"one\\" string'! comment")
    '\\'This is "one" string\\''
    >>> remove_user_nl_comment("'This! is \\"!one\\"! string'! comment")
    '\\'This! is "!one"! string\\''
    """

    #Create empty set for comment-delimiting indices:
    comment_delim_indices = set()

    #Search for all comment delimiters:
    for char_idx, char in enumerate(user_string):
        if char == comment_delim:
            #Add character index to set:
            comment_delim_indices.add(char_idx)
         #End if
    #End for

    #If no comments are present, then return string as-is:
    if not comment_delim_indices:
        return user_string
    #End if

    #Next, check if any single or double quotes are present:
    if not "'" in user_string and not '"' in user_string:
        #If no quotes, then cut-off string at first delimiter:
        return user_string[:sorted(comment_delim_indices)[0]]
    #End if

    #Create empty set for all character indices inside quotes:
    quote_text_indices = set()

    #Search for all text within quotes:
    quoted_text_matches = _QUOTE_REGEX.finditer(user_string)

    #Loop over all matches:
    for quote_match in quoted_text_matches:
        #Extract min/max indices of match:
        index_span = quote_match.span(0)
        #Add all indices to set:
        for index in range(index_span[0], index_span[1]):
            quote_text_indices.add(index)
        #End for
    #End for

    #Find all comment delimiters outside of quotes:
    non_quote_comment = comment_delim_indices.difference(quote_text_indices)

    if not non_quote_comment:
        #All comment delimiters are within quotes,
        #so return string as-is:
        return user_string
    #End if

    #Find first comment delimiter outside of quotes.
    #Everything to the right of it is part of the comment:
    return user_string[:sorted(non_quote_comment)[0]]

#####

def user_nl_str_to_int(string, var_name):

    """
    Checks if a string can be converted
    into an integer, and if not reports
    the relevant error.  This function
    is only used in the "get_user_nl_var_array_info"
    function below.
    ----------
    string   -> string to convert to integer.
    var_name -> name of the array variable
                associated with the string.

    doctests:

    1.  Check that a string with an integer can be
        converted properly:
    >>> user_nl_str_to_int("5", "banana")
    5

    2.  Check that a string with a float fails with
        the correct error:
    >>> user_nl_str_to_int("5.2", "banana") # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:...
    Invalid array index entry '5.2' used for variable 'banana' in 'user_nl_blom'.

    3.  Check that a string with a non-number fails with
        the correct error:
    >>> user_nl_str_to_int("a", "banana") # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:...
    Invalid array index entry 'a' used for variable 'banana' in 'user_nl_blom'.

    """

    #Attempt the conversion of the string to an integer:
    try:
        integer_val = int(string)
    except ValueError as verr:
        emsg = f"\nInvalid array index entry '{string}' "
        emsg += f"used for variable '{var_name}' in 'user_nl_blom'."
        raise OcnInParamGenError(emsg) from verr
    #End except

    #Return relevant integer value:
    return integer_val

#####

def check_dim_index(var_name, index_val, dim_size):

    """
    Checks that the user-specified index for the given
    variables is within the dimension size limit as
    specified in the namelist definition file.
    ----------
    var_name  -> Name of the array variable
                 associated with the string.
    index_val -> Index value provided by user
    dim_size  -> Maximum variable dimension size.

    doctests:

    1.  Check that an in-bounds index value
        returns nothing:
    >>> check_dim_index("banana", 5, 15)

    2.  Check that an index value that is
        too small returns the proper error:
    >>> check_dim_index("banana", 0, 15) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:...
    Variable 'banana' has index 0 in 'user_nl_blom', which is less than one (1), the minimal index value allowed.

    3.  Check that an index value that is
        too large returns the proper error:
    >>> check_dim_index("banana", 20, 15) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:...
    Variable 'banana' has index 20 in 'user_nl_blom', which is greater than the max dimension size of 15

    """

    #Make sure index is greater than zero:
    if index_val <= 0:
        emsg = f"\nVariable '{var_name}' has index {index_val}"
        emsg += " in 'user_nl_blom', which is less than one (1),"
        emsg += " the minimal index value allowed."
        raise OcnInParamGenError(emsg)
    #End if

    #Make sure index is not greater than max value:
    if index_val > dim_size:
        emsg = f"\nVariable '{var_name}' has index {index_val}"
        emsg += " in 'user_nl_blom', which is greater than the"
        emsg +=f" max dimension size of {dim_size}"
        raise OcnInParamGenError(emsg)
    #End if

#####

def parse_dim_spec(var_name, array_spec_text, dim_size):
    """
    Given the text of a single array dimension specification,
    return the range of values specified by the specification or
    raise an Exception if an error is detected.
    <var_name> is the variable name and is used for error messages
    <array_spec_text> is the text representation of the array spec
    <dim_size> is the size of that rank in <var_name>

    1. Check that a single, legal index returns the correct single value
    >>> parse_dim_spec('banana', '5', 10)
    [5]

    2. Check that a single, out-of-bounds index generates the proper error
    >>> parse_dim_spec('banana', '15', 10) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:
    Variable 'banana' has index 15 in 'user_nl_blom', which is greater than the max dimension size of 10
    >>> parse_dim_spec('banana', '0', 10) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:
    Variable 'banana' has index 0 in 'user_nl_blom', which is less than one (1), the minimal index value allowed.
    >>> parse_dim_spec('banana', '-2', 10) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:
    Variable 'banana' has index -2 in 'user_nl_blom', which is less than one (1), the minimal index value allowed.

    3. Check that a legal range returns the correct list of indices
    >>> parse_dim_spec('banana', '5:9', 10)
    [5, 6, 7, 8, 9]
    >>> parse_dim_spec('banana', ':9', 10)
    [1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> parse_dim_spec('banana', ':', 10)
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    >>> parse_dim_spec('banana', '6:', 10)
    [6, 7, 8, 9, 10]

    4. Check that an out-of-bounds range returns the correct list
    >>> parse_dim_spec('banana', '0:2', 10)
    [1, 2]
    >>> parse_dim_spec('banana', '7:11', 10)
    [7, 8, 9, 10]
    >>> parse_dim_spec('banana', '-1:11', 10)
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    5. Check that an empty range returns an empty list
    >>> parse_dim_spec('banana', '5:1', 10)
    []

    6. Check that a legal range with a stride returns the correct list
    >>> parse_dim_spec('banana', '5:9:2', 10)
    [5, 7, 9]
    >>> parse_dim_spec('banana', ':9:3', 10)
    [1, 4, 7]
    >>> parse_dim_spec('banana', '::3', 10)
    [1, 4, 7, 10]
    >>> parse_dim_spec('banana', '6:', 10)
    [6, 7, 8, 9, 10]
    >>> parse_dim_spec('banana', '9:1:-3', 10)
    [9, 6, 3]

    7. Check that a mismatched stride returns an empty list
    >>> parse_dim_spec('banana', '9:5:2', 10)
    []
    >>> parse_dim_spec('banana', '5:9:-2', 10)
    []
    >>> parse_dim_spec('banana', ':9:-3', 10)
    []
    >>> parse_dim_spec('banana', '::-2', 10)
    []
    >>> parse_dim_spec('banana', '6::-1', 10)
    []
    >>> parse_dim_spec('banana', '9:1:3', 10)
    []

    8. Check that a missing stride value generates an error
    >>> parse_dim_spec('banana', '2::', 10) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError:
    Two colons were provided for variable 'banana' in 'user_nl_blom', \
    but no stride value was provided.
    Please provide either a stride value, or remove the extra colon.
    """
    array_dims = [x.strip() for x in array_spec_text.split(':')]
    if len(array_dims) > 3:
        #Not sure what to do with three or more colons, so die here:
        emsg = f"Variable '{var_name}' has {len(array_dims) - 1} colons (:) "
        emsg += "listed in its dimension indexing in 'user_nl_blom'."
        emsg += " This is not a valid Fortran array section specification."
        raise OcnInParamGenError(emsg)
    #End if
    # Defaults
    arr_beg = 1
    arr_end = dim_size
    arr_stride = 1
    # Override start index?
    if array_dims[0]:
        arr_beg = user_nl_str_to_int(array_dims[0], var_name)
    # end if
    # Override end index?
    if len(array_dims) > 1:
        if array_dims[1].strip():
            arr_end = user_nl_str_to_int(array_dims[1], var_name)
        #End if (no else, blank means use default)
    else:
        # We only need to check this if it is only a single index
        check_dim_index(var_name, arr_beg, dim_size)
        # For a single index, the end is the same as the beginning
        arr_end = arr_beg
    #End if
    # Override stride?
    if len(array_dims) > 2:
        if array_dims[2]:
            arr_stride = user_nl_str_to_int(array_dims[2], var_name)
            if arr_stride == 0:
                emsg = f"Variable '{var_name}' has a stride of zero "
                emsg += "listed in its dimension indexing in 'user_nl_blom'."
                emsg += " This is not a valid Fortran stride."
                raise OcnInParamGenError(emsg)
            #End if
        else:
            emsg = f"Two colons were provided for variable '{var_name}'"
            emsg += " in 'user_nl_blom', but no stride value was provided."
            emsg += "\nPlease provide either a stride value, or remove the "
            emsg += "extra colon."
            raise OcnInParamGenError(emsg)
        #End if
    #End if (no else, just use default stride)
    # Now, create the set of entries
    # We need to modify the end to make the range function compatible with
    #    how Fortran uses it
    arr_end += int(arr_stride / abs(arr_stride))
    return [x for x in list(range(arr_beg, arr_end, arr_stride)) if
            ((x >= 1) and (x <= dim_size))]

#####

def _check_string_quotes(var_name, var_val):

    """
    Checks if a string is inside closed quotes,
    i.e. has both a starting and ending quote
    of the same type.  This function also
    raises an error if there are quotes but
    they aren't closed:

    doctests:

    1.  Check that a string with single quotes returns "True":
    >>> _check_string_quotes("Apple", "'Banana'")
    True

    2.  Check that a string with double quotes returns "True":
    >>> _check_string_quotes("Apple", '"Banana"')
    True

    3.  Check that a string without quotes returns "False":
    >>> _check_string_quotes("Apple", "Banana")
    False

    4.  Check that a string with mis-matching quote types raises
        the appropriate error:
    >>> _check_string_quotes("Apple", ''' "Banana' ''') # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError: Namelist entry 'Apple' is of type character but its input value:
    "Banana'
    has mis-matched quotes.  Please fix.

    5.  Check that a string with a missing ending quote type raises
        the appropriate error:
    >>> _check_string_quotes("Apple", "'Banana") # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError: Namelist entry 'Apple' is of type character but its input value:
    'Banana
    has mis-matched quotes.  Please fix.

    5.  Check that a string with a missing starting quote type raises
    the appropriate error:
    >>> _check_string_quotes("Apple", 'Banana"') # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError: Namelist entry 'Apple' is of type character but its input value:
    Banana"
    has mis-matched quotes.  Please fix.

    """

    #Make sure variable has been stripped:
    var_val_strip = var_val.strip()

    # If string is empty, just return
    if var_val_strip is '':
        return True

    #Set error message (just in case):
    emsg = f"Namelist entry '{var_name}' is of type character"
    emsg += " but its input value:"
    emsg += f"\n{var_val}\n"
    emsg += "has mis-matched quotes.  Please fix."

    #Check if starting and ending quotes exist and match:
    if var_val_strip[0] in _QUOTE_SET:
        if var_val_strip[0] == var_val_strip[-1]:
            #String is inside closed quotes:
            return True
        #End if

        #Starting and ending quotes don't match,
        #so raise an error:
        raise OcnInParamGenError(emsg)
    #End if

    #Check if there are ending quotes as well:
    if var_val_strip[-1] in _QUOTE_SET:
        #No starting quotes, raise an error:
        raise OcnInParamGenError(emsg)
    #End if

    #String is not inside quotes:
    return False

#####

def _get_nml_value_str(var_name, var_type, var_val):

    """
    Converts namelist variable inputs into their
    correct Fortran namelist value format
    ----------
    var_name -> Variable name (used for error message)
    var_type -> Variable type to convert to (logical, integer, real, character)
    var_val  -> Variable value to convert

    returns the fortran namelist-formatted variable
    value.

    doctests:

    1.  Check that a true logical variable outputs the correct value:
    >>> _get_nml_value_str("banana", "logical", "true")
    '.true.'

    2.  Check that a false logical variable outputs the correct value:
    >>> _get_nml_value_str("banana", "logical", "0")
    '.false.'

    3.  Check that an integer variable outputs the correct value:
    >>> _get_nml_value_str("banana", "integer", 5)
    '5'

    4.  Check that a real variable outputs the correct value:
    >>> _get_nml_value_str("banana", "real", "5d5")
    '5d5'

    5.  Check that a real variable with an integer value outputs
        the correct value:
    >>> _get_nml_value_str("banana", "real", 5)
    '5.d0'

    5.  Check that a character variable with no quotes outputs
        the correct value:
    >>> _get_nml_value_str("banana", "char*10", "apple")
    '"apple"'

    6.  Check that a character variable with quotes outputs
        the correct value:
    >>> _get_nml_value_str("banana", "char*250", " 'apple' ")
    "'apple'"

    7.  Check that a character variable with double quotes
        outputs the correct value:
    >>> _get_nml_value_str("banana", "char*N", ' "apple" ')
    '"apple"'

    8.  Check that a character variable with a quotation mark
        innternal to the string outputs the correct value:
    >>> _get_nml_value_str("banana", "char*31", ''' "app'le" ''')
    '"app\\'le"'

    9.  Check that a variable with an unknown type returns
        the proper error:
    >>> _get_nml_value_str("banana", "apple", "true") # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ocn_in_paramgen.OcnInParamGenError: Namelist type 'apple' for entry 'banana' is un-recognized.
    Acceptable namelist types are: logical, integer, real, or char*N.
    """

    #Create set for variable types
    #that don't need special handling:
    num_set = {"integer", "real"}

    #Check variable type:
    if var_type == 'logical':
        #If logical, then see if it is "True":
        if _is_nml_logical_true(var_name, var_val):
            return ".true."
        #End if
        #If not true, then must be false:
        return ".false."
    #End if

    if var_type in num_set:
        #Check if the variable value is an integer, but is being
        #used for a real-type variable:
        if var_type == "real" and isinstance(var_val, int):
            return f"{var_val}.d0"
        #End if

        #Otherwise, simply write value as-is:
        return f"{var_val}"
    #End if

    if "char" in var_type:
        #Remove extra white space:
        var_val_strip = var_val.strip()
        if '"' in var_val_strip:
            var_val_strip = var_val_strip.replace('"','')
            var_val_strip = "'" + var_val_strip + "'"

        #Check if string is wrapped in quotes:
        quoted_flag = _check_string_quotes(var_name, var_val_strip)

        #If not, then pass string with quotes:
        if not quoted_flag:
            return f'\'{var_val_strip}\''
        #End if

        #If so, then pass out original string as-is:
        return var_val_strip
    #End if

    #If one makes it here, then this is an un-recognized type option, so raise an error:
    emsg = f"Namelist type '{var_type}' for entry '{var_name}' is un-recognized.\n"
    emsg += "Acceptable namelist types are: logical, integer, real, or char*N."
    raise OcnInParamGenError(emsg)

################################################################
# MAIN "ocn_in" ParamGen class
################################################################

class OcnInParamGen(ParamGen):
    """
    Encapsulates data and read/write methods for
    the ocn_in Fortran namelist file and ParamGen
    object.
    """

    def __init__(self, pg_data_dict):

        """
        Initialize a ParamGen object directly
        using a ParamGen data dictionary, and
        create a new dictionary to match namelist
        variables to their associated groups
        when reading in "user_nl_blom".
        ----------
        pg_data_dict -> python dictionary with ParamGen keys/values

        """

        #Initialize ParamGen directly:
        super().__init__(pg_data_dict)

        #Create a namelist var/group dictionary,
        #which is used by the "append_user_nl_file"
        #method:
        self.__var_group_dict = {}

        #Create empty dictionaries that will contain
        #the namelist definition files and the set
        #of all namelist groups and variables:
        self.__nml_def_groups = {}
        self.__nml_def_vars   = {}

        #Set variables needed for ParamGen "reduction":
        self.__case = None
        self.__ocn_attr_dict = None

        #Initialize data structure for duplicate array
        #checking in user_nl_blom files.  This structure
        #is organized like so:
        # dict(var_name : list of sets)
        # list size = number of array dimensions specified
        # set contains the array indices specified for that
        # dimension:
        self.__set_index_vals = {}

    ####

    @classmethod
    def from_namelist_xml(cls, nml_xml_file):

        """
        Initialize ocn_in ParamGen object with XML file,
        ----------
        nml_xml_file -> path (str) to namelist definition XML file

        """

        #Create ParamGen object using base class:
        pg_xml = ParamGen.from_xml_nml(nml_xml_file, no_duplicates=True)

        #Initialize new "ocn_in" object:
        ocn_in_pg = OcnInParamGen(pg_xml.data)

        #Check if the new ParamGen object has all of the required
        #namelist elements:
        #----------------
        missing_elems = ocn_in_pg.check_nml_def_elems()

        if missing_elems:
            emsg = "The XML namelist definition file:\n"
            emsg += f"{nml_xml_file}\n"
            emsg += "has namelist entries that are missing required elements.\n"
            emsg += "Those entries and missing elements are:\n"
            for entry_id, missing_elems in missing_elems.items():
                emsg += f"{entry_id} : {', '.join(missing_elems)}\n"
            #End for
            raise OcnInParamGenError(emsg)
        #End if
        #----------------

        #Initialize file->group/var set dictionary:
        ocn_in_pg.__nml_def_groups[nml_xml_file] = set()
        ocn_in_pg.__nml_def_vars[nml_xml_file] = set()

        #Create namelist variable/group dictionary
        #and associated sets:
        #----------------
        for nml_group in ocn_in_pg._data:
            for var in ocn_in_pg._data[nml_group]:

                #Check if variable already exists in dictionary:
                if var in ocn_in_pg.__var_group_dict:
                    #No duplicate variables are allowed, even if
                    #in separate namelist groups, so raise an error.
                    #Please note that this error should always be
                    #caught earlier than this, so if it gets to this
                    #point something has gone seriously wrong:
                    emsg = f"Namelist entry id '{var}' exists"
                    emsg += f" in namelist group '{nml_group}'"
                    emsg += f" and '{ocn_in_pg.__var_group_dict[var]}'\n"
                    emsg += "Namelist variables can belong to only one group."
                    raise SystemError(emsg)
                #End if

                #If not, then add variable and group to dictionary:
                ocn_in_pg.__var_group_dict[var] = nml_group

                #Add namelist groups and variables to their
                #respective sets:
                ocn_in_pg.__nml_def_groups[nml_xml_file].add(nml_group)
                ocn_in_pg.__nml_def_vars[nml_xml_file].add(var)
            #End for
        #End for
        #----------------

        #Return object:
        return ocn_in_pg

    ####

    def check_nml_def_elems(self):

        """
        Function that checks if certain namelist definition
        file elements/tags that are optional for ParamGen
        but required are present
        """

        #Please note that "group" and "values" are automatically
        #required by the ParamGen schema.

        #Required namelist elements:
        req_elems = ["type", "desc", "category"]

        #Set missing  attributes dictionary:
        missing_elems = {}

        #Assume it is a ParamGen object, and loop over namelist groups:
        for nml_group in self._data:
            #Now loop over variables in group:
            for var in self._data[nml_group]:
                #Lastly loop over required namelist elements:
                for req_elem in req_elems:
                    #Check if required element is present:
                    if not req_elem in self._data[nml_group][var]:
                        #Add missing attribute to dictionary:
                        if var in missing_elems:
                            missing_elems[var].append(req_elem)
                        else:
                            missing_elems[var] = [req_elem]
                        #End if
                    #End if
                #End for
            #End for
        #End for

        #Return missing elements dictionary:
        return missing_elems

    ####

    def append_ocn_in_pg(self, ocn_pg_obj):

        """
        Append a new OcnInParamGen object
        to this one, ensuring that there are
        no duplicate namelist groups or variables.
        ----------
        ocn_pg_obj -> An OcnInParamGen object

        """
        #Loop over all XML files associated with input ocn_pg object:
        for input_file in ocn_pg_obj.__nml_def_groups:

            #Extract the group and variable sets from input PG object:
            input_groups = ocn_pg_obj.__nml_def_groups[input_file]
            input_vars   = ocn_pg_obj.__nml_def_vars[input_file]

            #Check that there are no matching namelist groups:
            #------------------------------------------------

            #Initialize error message string:
            emsg = ""

            #Loop over all namelist files and namelist group sets:
            for nml_file, nml_groups in self.__nml_def_groups.items():

                #Determine if any namelist groups are the same
                #between the two objects:
                same_groups = nml_groups.intersection(input_groups)

                #If so, then add to error message (as all namelist groups must be unique):
                if same_groups:
                    emsg += f"Cannot append:\n'{input_file}'\n"
                    emsg += " The following namelist groups conflict with those in"
                    emsg += f"\n'{nml_file} :'\n"
                    emsg += ", ".join(same_groups)
                #End if
            #End for

            #------------------------------------------------
            #Check that there are no matching namelist variables:
            #------------------------------------------------
            for nml_file, nml_vars in self.__nml_def_vars.items():

                #Determine if any namelist groups are the same
                #between the two objects:
                same_vars = nml_vars.intersection(input_vars)

                #If so, then add to error message (as all namelist variable ids must be unique):
                if same_vars:
                    emsg += f"Cannot append:\n'{input_file}'\n"
                    emsg += " The following namelist variablesconflict with those in"
                    emsg += f"\n'{nml_file} :'\n"
                    emsg += ", ".join(same_vars)
                #End if
            #End for
            #------------------------------------------------

        #End for (input files used to create input ocn_pb object)

        #Check if an error message was written.  If so then raise the
        #error(s) here:
        if emsg:
            raise OcnInParamGenError(emsg)
        #Endd if

        #Add input PG object dictionaries to this object's dicts:
        self.__nml_def_groups.update(ocn_pg_obj.__nml_def_groups)
        self.__nml_def_vars.update(ocn_pg_obj.__nml_def_vars)

        #Also combine PG object var-group dictionary needed for
        #appending "user_nl_blom":
        self.__var_group_dict.update(ocn_pg_obj.__var_group_dict)

        #Append input PG object to this object:
        self.append(ocn_pg_obj)

    ####

    def get_user_nl_var_array_info(self, var_str):

        """
        Checks whether the variable string
        is for a specific set of array
        indices.
        ----------
        var_str  -> variable name string.

        outputs:
        ----------
        is_array   -> Logical for whether variable
                      is an array.
        var_name   -> Name of variable
                      (with array indices stripped).
        arr_indxs  -> List of lists, with one list
                      for each array dimension. Each
                      dimension list contains all duplicated
                      indices for that dimension.
        data_group -> Namelist group for that particular
                      variable.

        """

        #Initialize variable name
        var_name = var_str

        #Initialize array index list:
        arr_indxs = []

        #Check for array syntax, i.e. parentheses:
        array_syntax_match = _ARR_INDEX_REGEX.search(var_str)

        #Extract variable name:
        if array_syntax_match:
            var_name = var_str[:array_syntax_match.start(0)]
        else:
            var_name = var_str
        #End if

        #Check that variable actually exists in ParamGen object:

        if var_name in self.__var_group_dict:
            #Extract namelist group list for variable:
            data_group = self.__var_group_dict[var_name]

        else:
            #Raise error that namelist variable isn't listed in
            #anywhere in a definition file:
            emsg = f"Variable '{var_name}' not found in any namelist definition files."
            emsg += " Please double-check 'user_nl_blom'."
            raise OcnInParamGenError(emsg)
        #End if

        #Extract variable type from ParamGen Object:
        var_type = self._data[data_group][var_name]["type"]

        #Search for array dimension specifications in type:
        array_type_dims = _ARRAY_TYPE_REGEX.search(var_type)

        #Determine if variable is actually an array or not:
        is_array = bool(array_type_dims)

        #Exit function here if no array indices were used in user_nl_blom file:
        if not array_syntax_match:

            #No parantheses used, so no indices need to be checked:
            return is_array, var_name, arr_indxs, data_group
        #End if

        #If variable is not an array, but array indices are being
        #used in user_nl_blom, then throw an error:
        if not is_array:
            emsg = f"Variable '{var_name}' is not an array, but array"
            emsg += " dimensions are being specified in 'user_nl_blom'."
            raise OcnInParamGenError(emsg)
        #End if

        #Extract array dimension information from variable type
        #as listed in the associated namelist definition file:
        #----------------------------------------------------

        #Pull out dimensions string:
        array_dim_text = array_type_dims.group(1)

        #Split text by number of commas (which should indicate dimensions):
        array_dims_list = array_dim_text.split(",")

        #Extract total number of dimensions:
        num_arr_dims = len(array_dims_list)

        #Create new list of max dim size:
        max_dim_sizes = []
        for dim_size in array_dims_list:
            max_dim_sizes.append(user_nl_str_to_int(dim_size, var_name))
        #End for

        #----------------------------------------------------

        #Now extract all text inside variable quotes:
        user_array_text = array_syntax_match.group(1)

        #Split text by number of commas (which should indicate dimensions):
        user_dim_text = user_array_text.split(",")

        #Check that the user hasn't listed the wrong number of dimensions:
        num_user_dims = len(user_dim_text)
        if num_user_dims != num_arr_dims:
            #Set proper grammar:
            if num_user_dims == 1:
                user_dim_str = "dimension"
            else:
                user_dim_str = "dimensions"
            #End if
            if num_arr_dims == 1:
                array_dim_str = "dimension."
            else:
                array_dim_str = "dimensions."
            #End if

            #Raise error with proper message:
            emsg = f"Variable '{var_name}' has {num_user_dims}"
            emsg += f" {user_dim_str} used in 'user_nl_blom', but is defined"
            emsg += f" to have {num_arr_dims} {array_dim_str}"
            raise OcnInParamGenError(emsg)
        #End if

        #Loop over dimensions:
        for dim_idx, array_index_text in enumerate(user_dim_text):
            #Create new array list entry:
            if array_index_text.strip() == ':':
                #Only a single colon provided.  In this case provide a special index
                #that indicates that specific indices can still be provided, but that the
                #whole array dimension cannot be written again:
                arr_indxs.append([-1])
            else:
                indices = parse_dim_spec(var_name, array_index_text,
                                         max_dim_sizes[dim_idx])
                if not indices:
                    ## Log a warning here if no values were returned?
                    pass
                #End if
                arr_indxs.append(indices)
            #End if
        #End for (dimensions)

        #Return relevant variables:
        return is_array, var_name, arr_indxs, data_group

    ####

    def check_array_indices(self, var_name, arr_index_list):

        """
        Checks whether the list of array indices has already
        been set for the given variable, and if so, raises
        an error.
        ----------
        var_name       -> Name of array variable being modified.

        arr_index_list -> A list of lists of array indices
                          with the first list representing
                          the dimensions, and the second list
                          containing the array indices being
                          set for that dimension.

        """

        #Initialize duplicated index flag,
        #We won't know the answer one way or the other
        #until either the end of the loop below, or
        #until certain conditions are met, so for now
        #initialize as "None":
        is_arr_dupl = None

        #Initialize "possible" array duplication flag:
        possible_dupl = False

        #Check if variable name exists in dictionary:
        if not var_name in self.__set_index_vals:
            #Create a new entry with an empty list,
            #it should then be filled out in the loop below:
            self.__set_index_vals[var_name] = []

            #Also set duplication to "False":
            is_arr_dupl = False
        #End if

        #Initialize dimension index for use in error-handling at
        #end of function:
        dim_indx = 0

        #Loop over each separate dimension list:
        for dim_indx, dim_arr_indxs in enumerate(arr_index_list):

            #Initialize duplicated index list (for last dimension checked):
            dup_indx_list = []

            #Check if dimension index exists for variable dictionary:
            if dim_indx == len(self.__set_index_vals[var_name]):
                #Create a new set of array indices for the new dimensions:
                self.__set_index_vals[var_name].append(set(dim_arr_indxs))

                #Since a new dimension is being specified, this is not a duplicate:
                is_arr_dupl = False
            else:
                #Loop over all array indices:
                for arr_indx in dim_arr_indxs:
                    #Check if array index has already been explicitly called:
                    if arr_indx in self.__set_index_vals[var_name][dim_indx]:
                        #Add array index to list of duplicated values:
                        dup_indx_list.append(arr_indx)

                        #This line is possibly a duplication,
                        #but will need to finish the loop to be sure:
                        possible_dupl = True
                    #End if

                    #Add index to "set index" set for variable:
                    self.__set_index_vals[var_name][dim_indx].add(arr_indx)
                #End for (array indices)

                #If there were no duplicates at this dimension, then this entry
                #is not a duplicate:
                if not possible_dupl:
                    is_arr_dupl = False
                #End if

            #End if (new dimension)
        #End for (dimensions)

        #If the duplication flag hasn't been set yet, then set it now:
        if is_arr_dupl is None:
            is_arr_dupl = possible_dupl
        #End if

        #Now raise an error if there is array duplication:
        if is_arr_dupl:
            if any(dup == -1 for dup in dup_indx_list):
                #This is a special case where a non-bounded
                #colon (:) was repeated twice, so write
                #the error message accordingly:
                emsg = f"Variable '{var_name}' has all values"
                emsg += " being set multiple times for"
                emsg += f" dimension {dim_indx+1}."
            else:
                emsg = f"Variable '{var_name}' has values"
                emsg += " at the following indices being"
                emsg += " set multiple times for dimension"
                emsg += f" ({dim_indx+1}) :\n"
                emsg += ", ".join(str(dup) for dup in dup_indx_list)
            #End if
            raise OcnInParamGenError(emsg)
        #End if

    ####

    def append_user_nl_file(self, user_nl_file):
        """
        Reads in user_nl_blom files and converts
        them to the proper ParamGen syntax.
        ----------
        user_nl_file -> Path (str) to user_nl_blom file.

        """

        #Create ordered dictionary to store namelist groups,
        #variables, and values from user_nl_blom file:
        _data = OrderedDict()

        #Initialize flag preventing duplicate namelist entries:
        no_duplicates = True

        #Initialize flag to mark whether a variable is an array:
        is_array = False

        #Initialize flag to mark whether the line is an array continuation line:
        is_continue_line = False

        #Open user_nl_blom file:
        with open(user_nl_file,'r', encoding='utf-8') as user_file:
            for line_num, line in enumerate(user_file):
                if len(line)>1:
                    #Split line into a list of words/characters:
                    line_s = line.split()

                    #If line is empty then go to next line:
                    if not line_s:
                        continue
                    #End if

                    #Check if a comment delimiter is somewhere in the string:
                    if '!' in line:
                        #Check if the entire line is a comment:
                        if line_s[0][0] == "!":
                            #Check if this comment is the duplicate keyword:
                            if "allow_duplicate_namelist_entries" in line_s:
                                #Next check if a user has set variable to True:
                                for word in line_s:
                                    if word.lower() == "true":
                                        #Allow duplicate namelist entries:
                                        no_duplicates = False
                                        break
                                    #End if
                                #End for
                            #End if
                            #Continue to next line in file:
                            continue
                        #End if
                        #Otherwise simply remove any part of the line that is commented out:
                        line = remove_user_nl_comment(line)
                    #End if

                    #Check if the first character on the line is a comma (,):
                    if line.strip()[0] == ",":
                        #Is this an array variable:
                        if is_array:
                            #Was a continuation line already provided:
                            if is_continue_line:
                                #Two commas were used in a row with a newline
                                #in-between. Technically this is allowed,
                                #but in practice it is VERY likely a mistake,
                                #so raise an error here:
                                emsg = f"Line number {line_num+1} in 'user_nl_blom'"
                                emsg += " starts with a comma (,) but the"
                                emsg += " previous line ended with a comma."
                                emsg += "\nPlease remove one of the commas."
                                raise OcnInParamGenError(emsg)
                            #End if

                            #If not, then set it to be a continuation line:
                            is_continue_line = True
                        else:
                            #The previous variable is not an array, so throw an error:
                            emsg  = f"Line number {line_num+1} in 'user_nl_blom'"
                            emsg += " starts with a comma (,) but the"
                            emsg += " associated namelist variable is not an array."
                            raise OcnInParamGenError(emsg)
                        #End if
                    #End if

                    #Now parse the line:
                    if "=" in line and (line.strip()[0] != "=") and not (is_array and is_continue_line):
                        line_ss  = line.split("=")       # Split line into before/after equals sign
                        var_str  = (line_ss[0]).strip()  # the first element is the variable name
                        var_str  = var_str.lower()

                        #Check if this variable is an array, and if so,
                        #then return what the variable name is, what indices (if any)
                        #are being specified, and what namelist (data) group it belongs to:
                        is_array, var_name, arr_indxs, data_group = \
                            self.get_user_nl_var_array_info(var_str)

                        #Check if variable can actually be set in user_nl_blom - or must be specified
                        #only via xmlchange
                        user_only_modify_via_xml = self._data[data_group][var_str]["modify_via_xml"] != 'unset'
                        if user_only_modify_via_xml:
                            xml_var = self._data[data_group][var_str]["modify_via_xml"]
                            emsg  = f"Cannot change {var_str} in user_nl_blom file,"
                            emsg += f" can only set {var_str} via xmlchange to {xml_var.upper()} \n"
                            raise OcnInParamGenError(emsg)

                        #Are there array indices specified:
                        if arr_indxs:

                            if no_duplicates:
                                #Check if any duplicate array indices are present:
                                self.check_array_indices(var_name, arr_indxs)
                            #End if

                            #ParamGen will think this is a "new" parameter variable, so we need
                            #to add a type in order for the "write" method to work properly.  The
                            #type can be copied directly from the original variable using "var_name":
                            var_type = self._data[data_group][var_name]["type"]
                        else:
                            #Variable doesn't need a type specified:
                            var_type = None
                        #End if (array indices)

                        #Extract value string:
                        val_str   = ' '.join(line_ss[1:]) # the rest is the value string

                        #Check if value string ends in array continuation:
                        if is_array:
                            #Check if the string ends in a comma (,):
                            if val_str.strip()[-1] == ",":
                                #If so, then make "array_continue_line" fully true:
                                is_continue_line = True
                            #End if
                        #End if

                        #Add the namelist group if not already in data dict:
                        if not data_group in _data:
                            _data[data_group] = {}
                        #End if

                        #Check if variable already exists in data dictionary:
                        if var_str in _data[data_group] and no_duplicates:
                            emsg = "Namelist variable '{}' set more than once in '{}'"
                            emsg += "\nPlease set each variable only once."
                            raise OcnInParamGenError(emsg.format(var_str, user_nl_file))
                        #End if

                        #Enter the parameter in the dictionary:
                        if var_type:
                            _data[data_group][var_str] = {'values':val_str, 'type':var_type}
                        else:
                            _data[data_group][var_str] = {'values':val_str}
                        #end if
                    elif (is_array and is_continue_line):
                        #See if there is an equals sign outside of quotes by treating it like
                        #a comment delimiter, and seeing if characters in the string are removed:
                        no_equals_line = remove_user_nl_comment(line, comment_delim='=')
                        if len(no_equals_line) != len(line):
                            #This looks like the start of a new namelist entry without the
                            #proper ending of the previous entry.  So raise an error here:
                            emsg = f"Line number {line_num+1} in 'user_nl_blom' appears"
                            emsg += " to be starting a new namelist entry,\nbut"
                            emsg += " the previous entry has a trailing comma (,).  Please fix."
                            raise OcnInParamGenError(emsg)
                        #End if

                        #This is an array continuation line, so append the line to previous
                        #variable's value as-is:
                        _data[data_group][var_str]['values'] += line

                        #Check if the line does NOT end in a comma (,):
                        if not line.strip()[-1] == ",":
                            #Notify loop to check the next line for a comma:
                            is_continue_line = False
                        #End if

                    else:
                        emsg = "Cannot parse the following line in '{}' :\n'{}'"
                        raise OcnInParamGenError(emsg.format(user_nl_file, line))
                    #End if ("=" sign check)
                #End if (len(line) > 1)
            #End for
        #End with

        #Check if any user_nl_blom data is present:
        if _data:
            #If so, then create new ParamGen object:
            pg_user = ParamGen(_data)

            #Append new user_nl_blom object to main ocn_in namelist object:
            self.append(pg_user)
        #End if
    ####

    def set_value(self, variable, value):
        """
        Reset value of namelist variable
        ----------
        """

        #Loop through namelist groups in alphabetical order:
        for nml_group in sorted(self._data):
            #Create function to properly sort variables with array indices:
            var_sort_key = lambda var : var[:var.index("(")] if "(" in var else var
            # Change value of variable
            for var in sorted(self._data[nml_group], key=var_sort_key):
                if var == variable:
                    self._data[nml_group][var]["values"] = value
                    break

    def get_value(self, variable):
        """
        Get value of namelist variable
        ----------
        """

        #Loop through namelist groups in alphabetical order:
        for nml_group in sorted(self._data):
            #Create function to properly sort variables with array indices:
            var_sort_key = lambda var : var[:var.index("(")] if "(" in var else var
            # Change value of variable
            for var in sorted(self._data[nml_group], key=var_sort_key):
                if var == variable:
                    value = self._data[nml_group][var]["values"]
                    break
        return value

    def write_nmlfile(self, output_path, groups):
        """
        Write data to Fortran namelist file.
        ----------
        output_path   -> path (str) to Fortran namelist (ocn_in) file
        """

        # Make sure ParamGen object has been reduced:
        if not self.reduced:
            emsg = "ParamGen object for ocn_in must be reduced before being "
            emsg += "written to file. Please check BLOM's buildnml script."
            raise OcnInParamGenError(emsg)
        #End if

        #Create function to properly sort variables with array indices:
        var_sort_key = lambda var : var[:var.index("(")] if "(" in var else var

        # Write Fortran namelist file:
        with open(os.path.join(output_path), 'w', encoding='utf-8') as ocn_in_fil:

            #Loop through namelist groups in alphabetical order:
            for nml_group in sorted(self._data):
                if nml_group in groups:
                    # Write namelist group:
                    ocn_in_fil.write("&"+nml_group+"\n")

                    # Write all variables within that group (sorted alphabetically):
                    for var in sorted(self._data[nml_group], key=var_sort_key):

                        #Extract variable value(s):
                        val = self._data[nml_group][var]["values"].strip()

                        #If no value is set then move to the next variable:
                        if val is None:
                            continue
                        #End if

                        #Extract variable type:
                        if "type" in self._data[nml_group][var]:
                            var_type = self._data[nml_group][var]["type"].strip()
                        else:
                            emsg = f"Namelist entry '{var}' is missing required 'type' element."
                            raise OcnInParamGenError(emsg)
                        #End if

                        # @ is used in a namelist to put the same namelist variable in multiple groups
                        # in the write phase, all characters in the namelist variable name after
                        # the @ and including the @ should be removed
                        if "@" in var:
                             var = re.sub("@.+$", "", var)

                        #Check if an array type:
                        array_type = _ARRAY_TYPE_REGEX.search(var_type)

                        if array_type:
                            #Grab all text before array regex match:
                            var_type = var_type[:array_type.start()].strip()

                            #Split the value into its array elements,
                            #this assumes that a comma (,) is the element
                            #delimiter:
                            array_elems = val.split(",")

                            #Write beginning of namelist entry:
                            nml_str = f"    {var} = "

                            #loop over array elements:
                            for elem in array_elems:

                                #Get properly-formatted variable value:
                                nml_val = _get_nml_value_str(var, var_type, elem)

                                #Add value string (with comma) to namelist string:
                                nml_str += (nml_val + ", ")
                            #End for

                            #There will always be a trailing comma and space (, ) so find it:
                            last_comma_idx = nml_str.rfind(", ")

                            #Write final string to file:
                            ocn_in_fil.write(nml_str[:last_comma_idx]+"\n")

                        else:  #Not an array

                            #Get properly-formatted variable value:
                            nml_val = _get_nml_value_str(var, var_type, val)

                            #Write variable to namelist file:
                            ocn_in_fil.write(f'    {var} = {nml_val}\n')

                        #End if (array type)
                    #End for (namelist variables)
                    # Add space for next namelist group:
                    ocn_in_fil.write('/\n\n')
                #End for (namelist groups)
            #End with (open ocn_in file)

    ####

    def write_inputdata(self, output_path, groups):

        """
        Write data to Fortran namelist file.
        ----------
        output_path   -> path (str) to blom input data list file

        """

        # Make sure ParamGen object has been reduced:
        if not self.reduced:
            emsg = "ParamGen object for ocn_in must be reduced before being "
            emsg += "written to file. Please check BLOM's buildnml script."
            raise OcnInParamGenError(emsg)
        #End if

        #Create function to properly sort variables with array indices:
        var_sort_key = lambda var : var[:var.index("(")] if "(" in var else var

        # Write Fortran namelist file:
        with open(os.path.join(output_path), 'w', encoding='utf-8') as ocn_in_fil:

            #Loop through namelist groups in alphabetical order:
            for nml_group in sorted(self._data):
                if nml_group in groups:

                    # Write all variables within that group (sorted alphabetically):
                    for var in sorted(self._data[nml_group], key=var_sort_key):

                        if self._data[nml_group][var]["is_inputdata"] == 'yes':
                            val = self._data[nml_group][var]["values"].strip()
                            val = val.replace('"',"'")
                            write_val = True
                            if val == "''":
                                write_val = False
                            elif 'unset' in val or 'UNSET' in val:
                                write_val = False
                            if write_val:
                                ocn_in_fil.write(f"{var} = {val} \n")


    def reduce_ocn_in(self, case, ocn_attr_dict):
        """
        Reduce XML namelist attributes
        (i.e. replace attribute/guard dictionary with value)
        ----------
        case          -> CIME case object
        ocn_attr_dict -> dictionary containing attribute values

        """

        # Set internal variables for use by "expand_func":
        self.__case = case
        self.__ocn_attr_dict = ocn_attr_dict

        # Reduce Param Data:
        self.reduce(self.__expand_func)

    ####

    def __expand_func(self, varname):

        """
        Function used to convert $XXX
        variables and XML attributes to
        their associated values.
        """

        #Check if varname matches a CIME case variable:
        val = self.__case.get_value(varname)

        #If not, then attempt to extract variable from
        #attribute dictionary:
        if val is None:
            if varname in self.__ocn_attr_dict:
                val = self.__ocn_attr_dict[varname]
            else:
                #Assume the XML attribute/guard is an empty string:
                val = ""
            #End if
        #End if

        #Return value if found:
        return val

############
#End of file
