""" Prep a series of graphs for the algorithm diagram figure
"""
from collections import defaultdict
from itertools import combinations

class otuMap(dict):
    def __init__(self):
        """
        Maps sequences to OTU assignments. The dict should be of the form {seqID: otuIndex}, e.g. {'seqA': 1, 'seqB': 2, 'seqC': 3}
        """
        super().__init__()

    @classmethod
    def from_list(cls, otu_list):
        """
        :param otu_list: should be of the form [{seqA, seqB}, {seqC, secD}, {secE}] where the index in the list is the OTU ID.
        """
        return cls({seq: idx for idx, otu in enumerate(otu_list)
                             for seq in otu})

    @property
    def otu_to_seq(self):
        otu_dict_set = defaultdict(set)
        for seq, otu_id in self.items():
            otu_dict_set[otu_id].add(seq)
        return otu_dict_set

    def conf_mat(self, dist_mat): # TODO
        """
        :param dist_mat: dict of sets of the form {seqA: {seqB, seqC}, seqB: {seqA}}
        """
        for seq1, seq2 in combinations(self.keys(), 2)

    def mcc(self, dist_mat): # TODO
        """
        Calculate the Matthews Correlation coefficient given a distance matrix.
        """


class OptiFit: # TODO
    def __init__(self, ref_otus, query_seqs, dist_mat):
        self.ref_otus = ref_otus
        self.query_seqs = query_seqs
        self.dist_mat = dist_mat