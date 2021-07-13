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

    def conf_mat(self, dist_mat):
        """
        :param dist_mat: dict of sets, e.g. {seqA: {seqB, seqC}, seqB: {seqA}}
        """
        conf_mat = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}
        for seq1, seq2 in combinations(self.keys(), 2):
            # are they actually similar (true assignment)
            true_otu = (seq2 in dist_mat[seq1]) or (seq1 in dist_mat[seq2])
            # were they assigned to the same OTU (predicted assignment)
            pred_otu = self[seq1] == self[seq2]

            if true_otu and pred_otu:
                result = 'tp'
            elif true_otu and not pred_otu:
                result = 'fn'
            elif not true_otu and pred_otu:
                result = 'fp'
            elif not true_otu and not pred_otu:
                result = 'tn'
            else:
                raise ValueError("this should never happen")
            conf_mat[result] += 1
        return conf_mat

    def mcc(self, conf_mat):
        """
        Calculate the Matthews Correlation coefficient given a confusion matrix.
        """
        return round((tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)), 2)


class OptiFit: # TODO
    def __init__(self, ref_otus, query_seqs, dist_mat):
        self.ref_otus = ref_otus
        self.query_seqs = query_seqs
        self.dist_mat = dist_mat