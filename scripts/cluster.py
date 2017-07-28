#!/usr/bin/env python
#
"""cluster.py: clusters attract output predictions. 
Usage: cluster.py attract.out ligand.pdb [options]
type   cluster.py --help for more information
"""  

import re
import sys
from optparse import OptionParser

import extract
from ptools import Rigidbody, Rmsd, WritePDB


class Struct:
    pass


# ====================
# from extract.py

class StructureI:
    def __cmp__(self, other):
        if self.trans < other.trans:
            return -1
        if self.trans > other.trans:
            return 1
        return cmp(self.rot, other.rot)
    pass


# ====================


def cluster(lig, structures, nclusters, cluster_memory, rmsd_cutoff, energy_cutoff):
    cluster_memory_ = cluster_memory + 1
    thecluster = []

    structures.sort(key=lambda i: i.ener)

    for s in structures:
        if s.ener > 0:
            break

        new = True
        sc = extract.rigidXMat44(lig, s.matrix)

        for c in reversed(thecluster[-cluster_memory:]):
            if (c.ext.ener - s.ener) < energy_cutoff and Rmsd(sc, c.structure) < rmsd_cutoff:
                c.count += 1
                new = False
                break

        if new:
            if len(thecluster) == nclusters + cluster_memory:
                # don't create more than nclusters + cluster_memory clusters
                # so the first nclusters  clusters will be fully filled
                break

            c = Struct()
            c.structure = sc
            c.ext = s
            c.count = 1
            thecluster.append(c)
            if len(thecluster) > cluster_memory:
                del thecluster[-cluster_memory_].structure

    for i in thecluster:
        if hasattr(i, "structure"):
            del i.structure
    thecluster.sort(key=lambda i: i.ext.ener)
    return thecluster[:nclusters]


if __name__ == "__main__":
    parser = OptionParser()
    parser.usage = 'cluster.py <out_file> <lig_file> [options]'
    parser.add_option("-e", "--energy_cutoff", action="store", type="float", dest="energy_cutoff",
                      help="Energy cutoff value (default=1000.0)", default=1000.0)
    parser.add_option("-r", "--rmsd_cutoff", action="store", type="float", dest="rmsd_cutoff",
                      help="Rmsd cutoff value (default=1.0)", default=1.0)
    parser.add_option("--nclusters", action="store", type="int", dest="nclusters",
                      help="number of cluster to output (default=200)", default=200)
    parser.add_option("-m", "--memory", action="store", type="int", dest="cluster_memory", default=50,
                      help="only the latest m clusters are compared during the clustering process, an increase of this value will increase significantly the time processing  (default=50)")
    parser.add_option("--extract", action="store", dest="extractTo", default=None, help="Extract the [nclusters] to files named 'prefix_rank_trans_rot.pdb'")

    (options, args) = parser.parse_args()

    outputfile = sys.argv[1]
    ligandfile = sys.argv[2]

    lig = Rigidbody(ligandfile)

    e = extract.Extractor(outputfile)  # extracts output structures or reuse the generated database
    validkeys = []
    regexp = re.compile("[0-9]+:[0-9]+")  # filter keys of the form "23:356"
    for k in e.d.keys():
        if regexp.match(k):
            validkeys.append(k)

    structures = []
    for k in validkeys:
        structures.append(e.d[k])

    thecluster = cluster(lig, structures,
                         options.nclusters,
                         options.cluster_memory,
                         options.rmsd_cutoff,
                         options.energy_cutoff)

    print "%-4s %6s %6s %13s %13s %6s %8s" % (" ", "Trans", "Rot", "Ener", "RmsdCA_ref", "Rank", "Weight")
    for i in range(len(thecluster)):
        print "%-4s %6s %6s %13.7f %s %6i %8s" % ("==", str(thecluster[i].ext.trans), str(thecluster[i].ext.rot), float(thecluster[i].ext.ener), thecluster[i].ext.rmsd, i + 1, str(thecluster[i].count))

    if options.extractTo is not None:
        for i, s in enumerate(thecluster):
            structure = extract.rigidXMat44(lig, s.ext.matrix)
            WritePDB(structure, "%s_%i_%i_%i.pdb" % (options.extractTo, i, s.ext.trans, s.ext.rot))
