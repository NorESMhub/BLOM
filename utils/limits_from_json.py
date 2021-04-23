import json
from distutils import util
from collections import OrderedDict
import os
#
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def str2value(value):
    """ """
    if value.lower() in ['true','false']:
        return util.strtobool(value.lower())
    elif value.isdigit():
        return int(value)
    elif isfloat(value):
        return float(value)
    else:
        return value
 

limits_in  = 'limits.json'
limits_out = 'limits_from_json'
user_in    = 'user_nl_blom'
skip_groups = ['MERDIA','SECDIA','CWMOD']
#
try:
    data = json.load(open(limits_in,'r'), object_pairs_hook=OrderedDict)
except FileNotFoundError:
    print(limits_in+' does not exist, aborting')
    exit
user_nl={}
if os.path.isfile(user_in):
    with open(user_in,'r') as f:
        user_data = f.readlines()
        for line in user_data:
            if line[0]!='#' and '=' in line:
                var = line.split('=')[0].strip()
                value = line.strip().split('=')[1].strip()
                print(var,value)
                if len(value.split(','))<=1:
                    user_nl[var] = [str2value(value)]
                    #value = distutils.util.strtobool(value.lower())
                elif len(value.split(','))>1:
                    dum = []
                    for v in value.split(','):
                        dum.append(str2value(v.strip()))
                    user_nl[var] = dum
else:
    print(user_in+' does not exist, continuing')

fout = open(limits_out,'w') 
linesize = 75 # limit line width
def_just = 12 # justify names
val_just = 12 # justify names
for group in data.keys():
    if group in skip_groups:
        continue
    # write definitions
    for var in data[group].keys():
        lineout = str(var)+' : '+data[group][var]['definition']
        c=0
        for w,word in enumerate(lineout.split()):
            c=c+len(word)+1
            if w==0:
                fout.write('! '+word.ljust(def_just)+' ')
            elif c>linesize: #linebreak if needed
                fout.write('\n')
                fout.write('! '.ljust(def_just+5)+word+' ')
                c=def_just+5
            # write the word
            else:
                fout.write(word+' ')
        fout.write('\n')
    #
    # write group tag and values
    fout.write('&'+group+'\n')
    for var in data[group].keys():
        if var in user_nl.keys():
            data[group][var]['value'] = user_nl[var]
        for v,val in enumerate(data[group][var]['value']):
            if v==0 and var not in ['MER_REGFLG']:
                lineout = '  '+var.ljust(val_just)+' = '
            elif var in ['MER_REGFLG']:
                dumstr = var+'('+str(v+1)+',:)'
                lineout = '  '+dumstr.ljust(val_just)+' ='
                for vv in val:
                    lineout = lineout+str(vv).rjust(2)
                    if vv!=val[-1]:
                        lineout = lineout+','
                fout.write(lineout+'\n')
                continue
            #
            if v>0:
                lineout = lineout+', '
            #
            if (val==True) and isinstance(val, bool):
                lineout = lineout+str('.true.').rjust(3)
            elif (val==False) and isinstance(val, bool):
                lineout = lineout+str('.false.').rjust(3)
            elif isinstance(val, float):
                if (abs(val)<1E-3 and val!=0.0) or abs(val)>1E3:
                    lineout = lineout+format(val,".3E").rjust(3)
                else:
                    lineout = lineout+str(val).rjust(3)
            else:
                lineout = lineout+str(val).rjust(3)
        if var not in ['MER_REGFLG']:
            fout.write(lineout+'\n')
    fout.write('/'+'\n')

fout.close()
