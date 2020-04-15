import math

#The cutoff distances of water and ions (should be modified as per the choice of the user).
#The distance is measured from a central atom of the QM region (not necessarily COM).
Rcutwater_sq = 400 #all water molecules within 20 angstrom of the QM region are incorporated in the MM charges of QMMM calculation.
Rcution_sq = 900   #all ions (NA and Cl) within 30 angstrom of the QM region are incorporated in the MM charges QMMM calculation.

# Open the input file
fin = open("em_topology_structure.txt","r")

# Store the lines of input file in a list
lines = fin.readlines()

# Length of the list containing the input
length =  len(lines)

#Close the input file
fin.close()

# Open the output file
fout = open("qmmm-setup_eomccsd_waterRcut225asq.in","w")

# Write the initial few lines in the input file
fout.write("$rem\njobtype sp\nmethod eom-ccsd\nee_states 2\nn_frozen_core fc\ncc_fno_thresh 9950\ncc_fno_usepop true\nmem_total 350000\nmem_static 80000\ncc_memory 150000\nbasis 3-21G\n$end\n\n$molecule\n1 1\n")

# Start the main loop over all the lines in the input file
for i in xrange(length):
    toks = lines[i].split()

    if (len(toks) ==8):
        if (toks[0] == "255LYR"):
            if (toks[1] == "CD"):
                coord_CDX = float(toks[3])*10
                coord_CDY = float(toks[4])*10
                coord_CDZ = float(toks[5])*10

            if (toks[1] != "CD" and toks[1] != "HD1" and toks[1] != "HD2" and (float(toks[2])> 3984 and float(toks[2]) < 4038)):

                newline= (toks[1][:1] + "   "+str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+"   "+str(float(toks[5])*10)+"\n")
                fout.write(newline)
                if (toks[1] == "C10"):
                    coord_C10X = float(toks[3])*10
                    coord_C10Y = float(toks[4])*10
                    coord_C10Z = float(toks[5])*10

                if (toks[1] == "CE"):
                    coord_CEX = float(toks[3])*10
                    coord_CEY = float(toks[4])*10
                    coord_CEZ = float(toks[5])*10

                    R = math.sqrt((coord_CEX - coord_CDX)*(coord_CEX - coord_CDX) + (coord_CEY - coord_CDY)*(coord_CEY - coord_CDY) + (coord_CEZ - coord_CDZ)*(coord_CEZ - coord_CDZ)) 

                    HX = coord_CEX + ((coord_CDX - coord_CEX)/R)*1.09  
                    HY = coord_CEY + ((coord_CDY - coord_CEY)/R)*1.09
                    HZ = coord_CEZ + ((coord_CDZ - coord_CEZ)/R)*1.09

                    newline= ("H" + "   "+str(HX)+"   "+str(HY)+"   "+str(HZ)+"\n")
                    fout.write(newline)
                if (toks[1] == "H173"):
                    fout.write("$end\n\n$external_charges\n") 
                    


###############################################################################################
# MM region construction: The external_charges part in the qchem input file

for i in xrange(length):
    toks = lines[i].split()

    ##### LOOP WHERE THE ATOMTYPE and ATOM NUMBER ARE NOT MERGED!!(PROTEIN and a few waters in here)(No CL/NA ions are here)
    if (len(toks) ==8):

        # Subsection for protein part of the MM region
        if (toks[0] != "255LYR" and i < 4350): # Number of lines to ensure that water and ions are left out of the loop
            newline= (str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+"   "+str(float(toks[5])*10)+ "   "+toks[7]+"\n")
            fout.write(newline)

        # Subsection for waters inside Rcut_sq of the chromophore
        if (toks[1] == "OW"):
            coord_OX = float(toks[3])*10
            coord_OY = float(toks[4])*10
            coord_OZ = float(toks[5])*10

            R_sq = (coord_OX - coord_C10X)*(coord_OX - coord_C10X) + (coord_OY - coord_C10Y)*(coord_OY - coord_C10Y) + (coord_OZ - coord_C10Z)*(coord_OZ - coord_C10Z)

            if (R_sq < Rcutwater_sq):         # Condition of the distance square
                toks2 = lines[i+1].split()  # Split the line that has the hydrogen (H1) of the accepted Oxygen
                toks3 = lines[i+2].split()  # Split the line that has the hydrogen (H2) of the accepted Oxygen

                newline=(str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+"   "+str(float(toks[5])*10)+ "   "+toks[7]+"\n")
                newline2=(str(float(toks2[3])*10)+"   "+str(float(toks2[4])*10)+"   "+str(float(toks2[5])*10)+ "   "+toks2[7]+"\n")
                newline3=(str(float(toks3[3])*10)+"   "+str(float(toks3[4])*10)+"   "+str(float(toks3[5])*10)+ "   "+toks3[7]+"\n")
                fout.write(newline)
                fout.write(newline2)
                fout.write(newline3)
            

        if (toks[0] == "255LYR"):

            if (toks[1] == "CG" or toks[1] == "HG1" or toks[1] == "HG2" or toks[1] == "CB" or toks[1] == "HB1" or toks[1] == "HB2"):
                newline= (str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+"   "+str(float(toks[5])*10)+ "   "+str(float(toks[7])+0.0115) +"\n")
                fout.write(newline)

            if (float(toks[2]) < 3976): 
                newline= (str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+"   "+str(float(toks[5])*10)+ "   "+toks[7]+"\n")
                fout.write(newline)

            if (float(toks[2]) >4037):
                newline= (str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+"   "+str(float(toks[5])*10)+ "   "+toks[7]+"\n")
                fout.write(newline)


    ##### LOOP WHERE THE ATOMTYPE and ATOM NUMBER ARE MERGED!! (No PROTEIN here) && (WATER and CL/NA ions are here)
    if (len(toks) == 7):
        if ("OW" in toks[1]):
            coord_OX = float(toks[2])*10
            coord_OY = float(toks[3])*10
            coord_OZ = float(toks[4])*10

            R_sq = (coord_OX - coord_C10X)*(coord_OX - coord_C10X) + (coord_OY - coord_C10Y)*(coord_OY - coord_C10Y) + (coord_OZ - coord_C10Z)*(coord_OZ - coord_C10Z)
            if (R_sq < Rcutwater_sq):         # Condition of the distance square
                toks2 = lines[i+1].split()  # Split the line that has the hydrogen (H1) of the accepted Oxygen
                toks3 = lines[i+2].split()  # Split the line that has the hydrogen (H2) of the accepted Oxygen

                newline=(str(float(toks[2])*10)+"   "+str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+ "   "+toks[6]+"\n")
                newline2=(str(float(toks2[2])*10)+"   "+str(float(toks2[3])*10)+"   "+str(float(toks2[4])*10)+ "   "+toks2[6]+"\n")
                newline3=(str(float(toks3[2])*10)+"   "+str(float(toks3[3])*10)+"   "+str(float(toks3[4])*10)+ "   "+toks3[6]+"\n")
                fout.write(newline)
                fout.write(newline2)
                fout.write(newline3)


        if ("NA" in toks[1]):
            coord_NAX = float(toks[2])*10
            coord_NAY = float(toks[3])*10
            coord_NAZ = float(toks[4])*10

            R_sq = (coord_NAX - coord_C10X)*(coord_NAX - coord_C10X) + (coord_NAY - coord_C10Y)*(coord_NAY - coord_C10Y) + (coord_NAZ - coord_C10Z)*(coord_NAZ - coord_C10Z)
#            print(R_sq)
            if (R_sq < Rcution_sq):         # Condition of the distance square
                newline=(str(float(toks[2])*10)+"   "+str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+ "   "+toks[6]+"\n")
                fout.write(newline)

        if ("CL" in toks[1]):
            coord_CLX = float(toks[2])*10
            coord_CLY = float(toks[3])*10
            coord_CLZ = float(toks[4])*10

            R_sq = (coord_CLX - coord_C10X)*(coord_CLX - coord_C10X) + (coord_CLY - coord_C10Y)*(coord_CLY - coord_C10Y) + (coord_CLZ - coord_C10Z)*(coord_CLZ - coord_C10Z)

            if (R_sq < Rcution_sq):         # Condition of the distance square
                newline=(str(float(toks[2])*10)+"   "+str(float(toks[3])*10)+"   "+str(float(toks[4])*10)+ "   "+toks[6]+"\n")
                fout.write(newline)

fout.write("$end")
fout.close()

