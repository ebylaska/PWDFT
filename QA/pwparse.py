#!/usr/bin/env python3

def evalnum(s):
   try:
      return int(s)
   except ValueError:
      return float(s)

def is_number(s):
   try:
      float(s)
      return True
   except ValueError:
      return False


def parse_energies(filename):
   estr = ""
   eoln = "\n"
   with open(filename) as file:
      a = file.read()
      for ln in a.split("\n"):
         if (("total     energy    :" in ln) or
             ("total orbital energy:" in ln) or
             ("hartree energy      :" in ln) or
             ("exc-corr energy     :" in ln) or
             ("ion-ion energy      :" in ln)):
            estr += ln + eoln
   return estr

def parse_optimization_energies(filename):
   ostr = ""
   eoln = "\n"
   with open(filename) as file:
      a = file.read()
      for ln in a.split("\n"):
         if "@ " in ln:
            ostr += ln + eoln
   return ostr

def parse_xyz_geometries(filename):
   xstr = ""
   eoln = "\n"
   with open(filename) as file:
      a = file.read()
      for xyz in a.split("XYZ format geometry\n")[1:]:
          try:
             nion =   evalnum(xyz.split("\n")[1])
             xstr += "\n".join(xyz.split("\n")[1:nion+3]) + eoln
          except:
             print("xyz=",xyz) 
   return xstr

def parse_ion_forces(filename):
   fstr = ""
   eoln = "\n"
   with open(filename) as file:
      a = file.read()
      for forces in a.split(" ion forces (au):\n")[1:]:
          try:
             fstr += forces.split("C.O.M.")[0] + eoln
          except:
             print("forces=",forces) 
   return fstr

def parse_apc_charges(filename):
   astr = ""
   eoln = "\n"
   with open(filename) as file:
      a = file.read()
      for apc in a.split("APC Potential:\n")[1:]:
         try:
            astr += "APC Potential: " + apc.split("\n")[0] + eoln
         except:
            print("apc potential=",apc) 
      for apc in a.split("APC Point Charges:\n")[1:]:
         try:
            astr += "APC Point Charges: " + apc.split("\n")[0] + eoln
         except:
            print("apc charges=",apc) 
   return astr



def compare_strings(str1,str2):
   result = "Failed"
   if (str1==str2): 
      result = "OK"
   else:
       list1 = str1.split()
       list2 = str2.split()
       if (len(list1)==len(list2)):
          result = "OK"
          for i in range(len(list1)):
              if (is_number(list1[i]) and is_number(list2[i])):
                  de = abs(evalnum(list1[i]) - evalnum(list2[i]))
                  if (de>1.0e-4): result = "Failed"
              else:
                  if (list1[i]!=list2[i]): result = "Failed"
   return result



def main():
   import sys,time,os
   import getopt,subprocess
   usage = \
   """
   This program compares the energies, optimization_energies, xyz_geometries, 
   and ion_forces of two pwdft output decks. If they are different the string "Failed" is output.
 
   Usage: pwparse pw1.out pw2.out
 


   Returned when checking output decks:
   Failed
   OK

   Example output when the files agree:
   Comparing Filenames: benzene/benzene.out benzene/check.out
   parse_energies: OK
   parse_optimization_energies: OK
   parse_xyz_geometries: OK
   parse_ion_forces: OK

   Example output when the files disagree:
   Comparing Filenames: benzene/benzene.out benzene/check.out
   parse_energies: OK
   parse_optimization_energies: OK
   parse_xyz_geometries: Failed
   parse_ion_forces: OK

   """

   opts, args = getopt.getopt(sys.argv[1:], "A:q")

   if len(args) < 2:
      print(usage)
      return

   filename1 = args[0]
   filename2 = args[1]
   print("Comparing Filenames: " + filename1 + " " + filename2)

   estr1 = parse_energies(filename1)
   estr2 = parse_energies(filename2)
   print("parse_energies: " + compare_strings(estr1,estr2))

   ostr1 = parse_optimization_energies(filename1)
   ostr2 = parse_optimization_energies(filename2)
   print("parse_optimization_energies: " + compare_strings(ostr1,ostr2))

   xstr1 = parse_xyz_geometries(filename1)
   xstr2 = parse_xyz_geometries(filename2)
   print("parse_xyz_geometries: " + compare_strings(xstr1,xstr2))

   fstr1 = parse_ion_forces(filename1)
   fstr2 = parse_ion_forces(filename2)
   print("parse_ion_forces: " + compare_strings(fstr1,fstr2))

   astr1 = parse_apc_charges(filename1)
   astr2 = parse_apc_charges(filename2)
   print("parse_apc_charges: " + compare_strings(astr1,astr2))



if __name__ == "__main__":
  main()

