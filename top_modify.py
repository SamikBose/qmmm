import numpy
import math
import os
fin = open("rhodopsin_charmm36.top","r")

lines = fin.readlines()

fout = open("rhodopsin_charmm36_EM_mod.top","w")
fout.write("Protein in water\n83703\n")

for i in xrange(4648): # This number is taken from top file(after this number of lines there is no info regarding the charges)
    toks = lines[i].split()
    if (len(toks) > 6 and toks[0] != ";" and toks[0][:1] != "#"):
        newline =(toks[4] + "   "+ toks[6]+"\n")
        fout.write(newline)
        
for i in xrange(122): # This number did not match with the top file. So we calculated the number of waters from the gro file
    newline = ("NA"+"   "+str(1.000)+"\n")
    fout.write(newline)


for i in xrange(120): # This number did not match with the top file. So we calculated the number of waters from the gro file
    newline = ("CL"+"   "+str(-1.000)+"\n")
    fout.write(newline)

for i in xrange(26371):# This number did not match with the top file. So we calculated the number of waters from the gro file
    newline = ("OW"+"   "+str(-0.834)+"\n"+"HW1"+"   "+str(0.417)+"\n"+"HW2"+"   "+str(0.417)+"\n")
    fout.write(newline)

fout.write("11.85500  11.80730   9.50000")

fout.close()
fin.close()

#os.system('paste md_withoutLIPID.gro rhodopsin_charmm36_mod.top > md_withoutLIPID_Charged-inc.txt')
