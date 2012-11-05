import Bio.SeqIO as SeqIO
import sys
import os
import time
import logging
import csv
import collections

log = logging.getLogger("NRDatabaseFilter")

def remove_common_words(words):
    """ Remove common words for the list of words given"""
    common_words = ["the", "for", "of" ]
    return [w for w in words if w not in common_words ]

def contains_all_keywords(words, keywords):
    """
        Returns true in all the words are present in the keywords.
        @param words The words to test.
        @param keywords The keywords that the words must contain
        @return Returns true if the words contain all the keywords
    """
    for keyword in keywords:
        if not keyword in words:
            return False
    return True


class NRDatabaseFilter:
    """
        Filter a NR database by splitting it according to the annotations that
        match the description of COGs obtained from IMG/M.
    """

    def __init__(self, fn_nr_database):
        """
            file of the NR database
        """
        self.nr = fn_nr_database
        self.ids = []

    def set_ids(self, ids):
        """
            Sets a set of ids for the descriptions. The are used to build filenames.
            If they are not provided, the files with be called database.nrXXX, where
            XXX are numeric characters.
            @param ids
        """
        log.debug("Ids for the descriptions set: %s", ids)
        self.ids = ids

    def set_descriptions(self, descriptions):
        """
            @param set of descriptions used to filter sequences. The descriptions are
            split in keywords and common words are removed from the keywords

        """
        log.debug("Setting descriptions to filter")
        if len(descriptions) == 0:
            raise ValueError("The list of descriptions is empty")
        self.descriptions_keywords = []
        for d in descriptions:
            keywords = d.lower().split(" ")
            self.descriptions_keywords.append(keywords)

    def do_filtering(self,overwrite=True):
        """
            Does the filtering seleting the sequences with descriptions matching
            the descriptions stored in the class
            @param overwrite It true, all the previous NR files created for the COGs
            are overwritten.
        """
        log.info("Filtering database with the descriptions provided")
        if len(self.descriptions_keywords) == 0:
           raise ValueError("The database cannot be filtered. No descriptions " \
           "to extract are given")
        if len(self.ids) == 0 and len(self.ids) != len(self.descriptions_keywords):
            log.debug("self.ids: %s",self.ids)
            log.debug("self.descriptions_keywords: %s",self.descriptions_keywords)
            raise ValueError("The number of ids provided do not match the number " \
                "of descriptions. An id per description is needed. ")
        if overwrite:
            map(os.remove, self.get_database_files_created())
        for seq_record in SeqIO.parse(self.nr, "fasta"):
            found, ind = self.compare_to_descriptions(seq_record.description)

            # if the description for the sequence is found in the descriptions
            # to extract, write to its database file
            n_ids = len(self.ids)
            if found:
                log.debug("Found match for %s" % seq_record.description)
                if not n_ids:
                    fn = "nr.filtered.{0:06d}".format(ind)
                else:
                    fn = self.get_database_name(self.ids[ind])
                    log.debug("Appending %s (index %d) to %s",self.ids[ind],ind,fn)
                f = open(fn, "a")
                f.write(seq_record.format("fasta"))
                f.close()

    def compare_to_descriptions(self, description):
        """ Compare a description to all the descriptions that are going to
            be filtered. If the description matches any of the descriptions
                             stored in the class, returns True.
            @param description The description to compare
            @return A bool value and an index. The bool is true if the
            description is found within the descriptions stored in the class.
            The index is the index of the description that matches the
            argument
        """
#        log.debug("Checking description %s",description)
        words = description.lower().split(" ")
        for i, keywords in enumerate(self.descriptions_keywords):
#            log.debug("Checking keywords: %s",keywords)
            if contains_all_keywords(words, keywords):
                return True, i
        return False, 0

    def get_database_files_created(self):
        """
            returns the names of all the files created by the function
            do_filtering()
        """
        if len(self.ids) == 0:
            return glog.glob("nr.filtered*")
        else:
            fns = []
            for i in self.ids:
                fn = self.get_database_name(i)
                if os.path.exists(fn):
                    fns.append(fn)
            return fns

    def get_database_name(self,identifier):
        """
            @param identifier A string
            @returns A name for a database file based on the identifier
        """
        return "nr." + identifier