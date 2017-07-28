# this is a temporary script to convert incorrect pdb files generated with a buggy version
# of ptools into hopefully correct pdb files.
"""pdb2pdb.py: converts incorrect pdb files generated with a buggy version of ptools into hopefully correct pdb files. 
Usage:  pdb2pdb.py incorrect.pdb > correct.pdb
"""

# usage:  pdb2pdb.py incorrect.pdb > correct.pdb


import sys

# try to see if the pdb has bugs...
# we check for Calpha atoms a see if they are declared as calcium would be
# however we cannot rely on ptools to do that (yet)

lines = open(sys.argv[1]).readlines()

atoms = [i for i in lines if i[:4] == "ATOM"]

buggy = None
for line in atoms:
    name = line[12:16]
    if name.rstrip() == "CA":
        buggy = True
        break
    elif name.rstrip() == " CA":
        buggy = False

if buggy is None:
    raise RuntimeError("could not decide if it is a buggy pdb file or not")

if not buggy:
    print "ok"
else:
    # we shift 1 to 3-letters atom types to the right:
    for line in atoms:
        name = line[12:16]
        lstline = list(line)

        if len(name.strip()) <= 3:
            to_insert = name.strip()
            lstline[12] = ' '
            lstline[13:13 + len(to_insert)] = list(to_insert)
            print "".join(lstline),
