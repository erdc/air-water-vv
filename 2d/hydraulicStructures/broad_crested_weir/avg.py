#!/usr/bin/env python
import sys
newton_its=0.0
with open(sys.argv[1],'r') as f:
    for line in f.readlines():
        if "Iteration Statistics" in line:
            print "Average Iteration Statistics\n==================================="
        elif "_r" in line:
            if newton_its > 0.0:
                avg/=newton_its
                print "average ksp its ",avg," over ", newton_its," newton iterations"
                print "max ksp its ",maxits," over ", newton_its," newton iterations\n"
                print line[:-1]
            else:
                print line[:-1]
            avg=0.0
            maxits=0.0
            newton_its=0.0
        elif "reason = 4" in line:
            continue
        else:
            words = line.split()
            for i,w in enumerate(words):
                if w == "ksp.its=":
                    avg+=float(words[i+1])
                    maxits=max(float(words[i+1]),maxits)
                    newton_its+=1.0
    if newton_its > 0.0:
        avg/=newton_its
        print "average ksp its ",avg," over ", newton_its," newton iterations\n"


