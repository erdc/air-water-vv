#!/usr/bin/env python
import sys
import re
ref_list = []
name_list = []
tol_list = ["no_tol_placeholder"]
name_list_no_tol = []
with open(sys.argv[1],'r') as r:
    for line in r.readlines():
        if "selfp_boomerAMG_" in line:
            fname=line[:-5]
            fnum=line[-3:-1]
            ftol=line[16:21]
            if not bool(re.search('preonly',fname)): # if not in the ignore list
                if fnum not in ref_list:
                    ref_list.append(fnum)
                if fname not in name_list and not bool(re.search('preonly',fname)):
                    name_list.append(fname)
                if ftol not in tol_list:
                    if bool(re.search('e-',ftol)):
                        tol_list.append(ftol)
                    elif fname not in name_list_no_tol:
                        name_list_no_tol.append(fname)
for take in [1,2]:
    print ""
    for tol in tol_list:
        name_list_tol = []
        if tol == "no_tol_placeholder":
            name_list_tol = name_list_no_tol
            for name in name_list_tol:
                print name
        else:
            for name in name_list:
                ftol=name[16:21]
                if bool(re.search(tol,ftol)):
                    print name
                    name_list_tol.append(name)
        trigger1 = False
        trigger2 = False
        for num in ref_list:
            first = True
            with open(sys.argv[1],'r') as r:
                for line in r.readlines():
                    if trigger2:
                        if first:
                            print str(num)+" "+line[:-1]
                            first = False
                        else:
                            print line[:-1]
                        trigger2 = False
                    elif trigger1:
                        if take == 1:
                            if first == True:
                                print str(num)+" "+line[:-1]
                                first = False
                            else:
                                print line[:-1]
                        elif take == 2:
                            trigger2 = True
                        trigger1 = False
                    elif "selfp_boomerAMG_" in line:
                        fname=line[:-5]
                        fnum=line[-3:-1]
                        if fnum == num:
                            if fname in name_list_tol:
                                trigger1 = True
