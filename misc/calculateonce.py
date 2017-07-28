# class to avoid recalculating if a database of precalculated results is available
"""Defines the class NoRecalc that avoids recalculating results that are stored in an available database.
"""

import os
import shelve
import sys
from exceptions import NotImplementedError
from stat import ST_MTIME


def devnull(string):
    # no output at all!
    return


class NoRecalc(object):

    def __init__(self, databasefile, dependencies, args=[], info_output=devnull, error_output=sys.stderr.write):
        """
        databasefile: file that contains precalculated results
        dependencies: list of files that are used to generate the databasefile
        args: arguments used to generate the output file
        generate: function that is called with (dependencies) to generate results
        load: function that loads pregenerated results from databasefile
        generate and load must both return the same values, ie data extracted from
        the database
        """

        self.dependencies = dependencies
        self.databasefile = databasefile
        self.info_output = info_output

        self.db = None
        self.norecalcCacheFile = ".%s.norecalc" % databasefile

        flag1 = os.path.exists(databasefile)
        flag2 = os.path.exists(self.norecalcCacheFile)

        if flag1 and flag2:
            self.info_output("Database file (%s) exists\n" % (databasefile))
            # check if up-to-date:
            uptodateflag = True  # if database is more recent than dependencies, no need to regenerate it
            statdb = os.stat(databasefile)
            for dep in dependencies:
                try:
                    statout = os.stat(dep)
                except OSError:
                    error_output("dependency file %s is missing !\n" % dep)
                    raise
                if statdb[ST_MTIME] < statout[ST_MTIME]:
                    uptodateflag = False
                    self.info_output("%s is more recent than database, going to regenerate database\n" % dep)

            # check if current arguments are compatible with stored database:
            if not self.checkArgs(databasefile, args):
                print "new arguments invalidate database"
                uptodateflag = False

            if uptodateflag:
                # we can use the database here :-)
                self.db = self.load(databasefile)
                return

        # we must delete the database if it exists and generate a new one
        self.info_output("Removing the old database files (if present)\n")
        if flag1:
            os.remove(databasefile)
        if flag2:
            os.remove(self.norecalcCacheFile)

        # we must now generate database and cache file
        self.db = self.generate(databasefile, dependencies, args)
        d = shelve.open(self.norecalcCacheFile)
        d["dependencies"] = dependencies
        d["databasefile"] = databasefile
        # d["args"]=args   args must be picklable but ptools is not (yet)
        d["status"] = "finished"
        return

    def generate(self, databasefile, dependencies, args):
        """generate database result file, must return data extracted
        from database (like load).

        To be implemented in children classes.
        """
        raise NotImplementedError()

    def checkArgs(self, databasefile, args):
        """check if the database must be depreciated because input parameters invalidate the previous simulation
           return True if the database is OK
           return False if the database need to be rebuilt"""
        self.info_output("Warning: database validity with respect to your arguments was not checked\n")
        return True

    def load(self, databasefile):
        """To be implemented in children classes."""
        raise NotImplementedError()
