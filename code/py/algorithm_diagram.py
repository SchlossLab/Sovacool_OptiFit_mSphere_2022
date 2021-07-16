""" Prep a series of graphs for the algorithm diagram figure.
"""
from collections import defaultdict, Counter
from itertools import combinations


def mcc(conf_mat):
    """
    Calculate the Matthews Correlation coefficient given a confusion matrix.
    """
    tp = conf_mat["tp"]
    tn = conf_mat["tn"]
    fp = conf_mat["fp"]
    fn = conf_mat["fn"]
    return round(
        (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)), 2
    )


def classify(true_otu, pred_otu):
    """
    Classify a prediction as a true positive (tp), true negative (tn),
        false positive (fp), or false negataive (fn).
    """
    if true_otu and pred_otu:
        result = "tp"
    elif true_otu and not pred_otu:
        result = "fn"
    elif not true_otu and pred_otu:
        result = "fp"
    elif not true_otu and not pred_otu:
        result = "tn"
    else:
        raise ValueError("this should never ever happen")
    return result


def dist_pairs_to_sets(dframe):
    """
    Convert a dict of 2 lists, with each sequence at the same index being within
        the distance threshold, to a dict of sets.
    :param dframe: {'seq1': ['A', 'A', 'B'],
                    'seq2': ['B', 'C', 'D']}
    :return: {'A': {'B','C'}, 'B': {'A','C','D'}, 'C': {'A','B'}, 'D': {'B'}}
    """
    dist_set = defaultdict(set)
    for seq1, seq2 in zip(dframe["seq1"], dframe["seq2"]):
        dist_set[seq1].add(seq2)
        dist_set[seq2].add(seq1)
    return dist_set


class otuMap():
    """
    Maps sequences to OTU assignments. The dict should be of the form
        {seqID: otuIndex}, e.g. {'seqA': 1, 'seqB': 2, 'seqC': 3}
    """

    def __init__(self, seqs_to_otus = None):
        self.seqs_to_otus = seqs_to_otus if seqs_to_otus else dict()

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.seqs_to_otus)

    @classmethod
    def from_list(cls, otu_list):
        """
        :param otu_list: should be of the form
            [{seqA, seqB}, {seqC, secD}, {secE}]
            where the index in the list is the OTU ID.
        :return: an otuMap
        """
        return cls(seqs_to_otus = {seq: idx for idx, otu in enumerate(otu_list) 
                                            for seq in otu})

    @property
    def otus_to_seqs(self):
        """
        :return: dictionary of sets, with keys as OTU IDs and values as sets
            containing sequence IDs.
        """
        otu_dict_set = defaultdict(set)
        for seq, otu_id in self.seqs_to_otus.items():
            otu_dict_set[otu_id].add(seq)
        return otu_dict_set

    def conf_mat(self, dist_mat, baseline_conf_mat=False):
        """
        :param dist_mat: dict of sets, e.g. {seqA: {seqB, seqC}, seqB: {seqA}}
        :param baseline_conf_mat: provide a confusion matrix (as a dictionary)
            with starting values other than zero.
        :return: a confusion matrix as a dictionary of counts containing
            the keys 'tp', 'tn', 'fp', and 'fn'
        """
        # build list of tp, tn, fp, & fn
        classes = [
            classify(
                (seq2 in dist_mat[seq1]) or (seq1 in dist_mat[seq2]),
                self.seqs_to_otus[seq1] == self.seqs_to_otus[seq2],
            )
            for seq1, seq2 in combinations(self.seqs_to_otus.keys(), 2)
        ]
        # convert to a dictionary of counts
        conf_mat = Counter(classes)
        if baseline_conf_mat:
            conf_mat.update(baseline_conf_mat)
        return conf_mat


class OptiFit:  # TODO
    def __init__(self, ref_otus, query_seqs, dist_mat):
        self.ref_otus = ref_otus
        self.query_seqs = query_seqs
        self.dist_mat = dist_mat

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__.items())


def main():
    conf_mat = {"tp": 12, "tn": 1210, "fp": 0, "fn": 3}
    dist_frame = {
        "seq1": [
            "D",
            "F",
            "G",
            "H",
            "I",
            "I",
            "J",
            "J",
            "N",
            "O",
            "P",
            "P",
            "P",
            "Q",
            "Q",
        ],
        "seq2": [
            "B",
            "E",
            "C",
            "A",
            "B",
            "D",
            "A",
            "H",
            "M",
            "L",
            "K",
            "L",
            "O",
            "E",
            "F",
        ],
    }
    dist_mat = dist_pairs_to_sets(dist_frame)
    otu_list = [{}, {'I', 'D', 'B'}, {'F', 'E', 'Q'}, {'C', 'G'}, {'H', 'J', 'A'}, {'M', 'N'}, {'P', 'L', 'O'}, {'K'}]
    otus = otuMap.from_list(otu_list)
    print(otus)


if __name__ == "__main__":
    main()
