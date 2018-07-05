from builtins import input
import os
import sys

#Listing files and storing them to list
files = os.listdir(".")


#Interactively asking for old and new name
oldname = input("Give old name:")
newname = input("Give new name: ")




for ff in files:
#Checking only python files
    if ".py" in ff:
#Reading python file and storing
        fid = open(ff,"r")
        pyfile = fid.readlines()
        fid.close()
#Openfinf python file and replacing with new name
        fid = open(ff,"w")
        for line in pyfile:
            fid.write(line.replace(oldname,newname))
        fid.close()
            
