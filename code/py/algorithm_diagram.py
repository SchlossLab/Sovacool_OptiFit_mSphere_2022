""" Prep a series of graphs for the algorithm diagram figure
"""
from collections import defaultdict, Counter
from itertools import combinations


def mcc(conf_mat):
    """
    Calculate the Matthews Correlation coefficient given a confusion matrix.
    """
    tp = conf_mat['tp']
    tn = conf_mat['tn']
    fp = conf_mat['fp']
    fn = conf_mat['fn']
    return round((tp * tn - fp * fn) /
                 sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)),
                 2)


def classify(true_otu, pred_otu):
    """
    Classify a prediction as a true positive (tp), true negative (tn),
    false positive (fp), or false negataive (fn).
    """
    if true_otu and pred_otu:
        result = 'tp'
    elif true_otu and not pred_otu:
        result = 'fn'
    elif not true_otu and pred_otu:
        result = 'fp'
    elif not true_otu and not pred_otu:
        result = 'tn'
    else:
        raise ValueError("this should never ever happen")
    return result


def dist_pairs_to_sets(dframe):
    """
    Convert a dict of 2 lists, with each sequence at the same index being within
    the distance threshold, to a dict of sets.
    """
    dist_set = defaultdict(set)
    for seq1, seq2 in zip(dframe['seq1'], dframe['seq2']):
        dist_set[seq1].add(seq2)
        dist_set[seq2].add(seq1)
    return dist_set

class otuMap(dict):
    def __init__(self, *args, **kwargs):
        """
        Maps sequences to OTU assignments. The dict should be of the form
          {seqID: otuIndex}, e.g. {'seqA': 1, 'seqB': 2, 'seqC': 3}
        """
        super().__init__(*args, **kwargs)

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__.items())

    @classmethod
    def from_list(cls, otu_list):
        """
        :param otu_list: should be of the form
          [{seqA, seqB}, {seqC, secD}, {secE}]
          where the index in the list is the OTU ID.
        """
        return cls({seq: idx for idx, otu in enumerate(otu_list)
                             for seq in otu})

    @property
    def otu_to_seq(self):
        """
        :return: dictionary of sets, with keys as OTU IDs and values as sets
          containing sequence IDs.
        """
        otu_dict_set = defaultdict(set)
        for seq, otu_id in self.items():
            otu_dict_set[otu_id].add(seq)
        return otu_dict_set

    def conf_mat(self, dist_mat, baseline_conf_mat = False):
        """
        :param dist_mat: dict of sets, e.g. {seqA: {seqB, seqC}, seqB: {seqA}}
        :param baseline_conf_mat: provide a confusion matrix (as a dictionary)
          with starting values other than zero.
        """
        # build list of tp, tn, fp, & fn
        classes = [classify((seq2 in dist_mat[seq1]) or (seq1 in dist_mat[seq2]), self[seq1] == self[seq2]) for seq1, seq2 in combinations(self.keys(), 2)
        # convert to a dictionary of counts
        conf_mat = Counter(classes)
        if baseline_conf_mat:
            conf_mat.update(baseline_conf_mat)
        return conf_mat


class OptiFit: # TODO
    def __init__(self, ref_otus, query_seqs, dist_mat):
        self.ref_otus = ref_otus
        self.query_seqs = query_seqs
        self.dist_mat = dist_mat

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__.items())


def main():
    conf_mat = {'tp': 12, 'tn': 1210, 'fp': 0, 'fn': 3}
    dist_frame = {'seq1': ["D", "F", "G", "H", "I", "I", "J", "J", "N", "O", "P",
                         "P", "P", "Q", "Q"],
                  'seq2': ["B", "E", "C", "A", "B", "D", "A", "H", "M", "L", "K",
                         "L", "O", "E", "F"]}
    dist_mat = dist_pairs_to_sets(dist_frame)


if __name__ == "__main__":
    main()
