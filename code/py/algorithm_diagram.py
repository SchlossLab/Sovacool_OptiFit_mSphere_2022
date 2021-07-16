""" Prep a series of graphs for the algorithm diagram figure.
"""
from collections import defaultdict, Counter
from copy import deepcopy
from itertools import combinations
from math import comb, sqrt
import pandas as pd


def opticlust_example():
    dist_frame = {
        "seq1": [ "D", "F", "G", "H", "I", "I", "J", "J", "N", "O", "P", "P", "P", "Q", "Q", ],
        "seq2": [ "B", "E", "C", "A", "B", "D", "A", "H", "M", "L", "K", "L", "O", "E", "F", ],
    }
    assert len(dist_frame['seq1']) == len(dist_frame['seq2'])
    dist_mat = dist_pairs_to_sets(dist_frame)
    otu_list = [{}, {'I', 'D', 'B'}, {'F', 'E', 'Q'}, {'C', 'G'},
                {'H', 'J', 'A'}, {'M', 'N'}, {'P', 'L', 'O'}, {'K'}]
    otus = otuMap.from_list(otu_list, dist_mat = dist_mat, n_seqs = 50)
    #print(otus)
    #conf_mat = otus.conf_mat(dist_mat)
    #print('mcc current:', mcc(conf_mat),
    #      '\nmcc from correct conf mat:',
    #      mcc({"tp": 14, "tn": 1210, "fp": 0, "fn": 1}),
    #      '\nmcc correct: 0.97')
    return otus


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


class otuMap:
    """
    Maps sequences to OTU assignments.
    """

    def __init__(self, seqs_to_otus = None, dist_mat = None, n_seqs = 0):
        """
        :param seqs_to_otus: dict of the form {seqID: otuIndex},
            e.g. {'seqA': 1, 'seqB': 2, 'seqC': 3}.
        :param n_seqs: total number of sequences in dataset.
        """
        self.seqs_to_otus = seqs_to_otus if seqs_to_otus else dict()
        self.dist_mat = dist_mat if dist_mat else dict()
        self.n_seqs = n_seqs
        assert n_seqs >= len(self.seqs_to_otus)

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.seqs_to_otus)

    @classmethod
    def from_list(cls, otu_list, **kwargs):
        """
        :param otu_list: should be of the form
            [{seqA, seqB}, {seqC, secD}, {secE}]
            where the index in the list is the OTU ID.
        :return: an otuMap
        """
        return cls(seqs_to_otus = {seq: idx for idx, otu in enumerate(otu_list)
                                            for seq in otu}, **kwargs)

    @property
    def seqs(self):
        return set(self.seqs_to_otus.keys())

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

    @property
    def ghost_pairs(self):
        """
        ghost_pairs: number of pairs from ghost sequences, calculated from the
            distance matrix and total number of seqs (n_seqs).
            Ghost sequences are not similar enough to any other seq
            to be included in the distance matrix, thus they form singleton OTUs
            and contribute to the number of true negatives.
            These are not shown in the otuMap in order to save space.
        """
        n_sim_seqs = len(self.dist_mat)
        n_unsim_seqs = self.n_seqs - n_sim_seqs
        # number of distances within the distance threshold, i.e. they're included in dist_mat
        #n_dists = sum([len(dist_mat[s1]) for s1 in dist_mat]) / 2
        return comb(n_unsim_seqs, 2) + n_unsim_seqs * n_sim_seqs

    @property
    def conf_mat(self):
        """
        :return: a confusion matrix as a dictionary of counts containing
            the keys 'tp', 'tn', 'fp', and 'fn'
        """
        # build list of tp, tn, fp, & fn
        classes = [
            classify(
                (seq2 in self.dist_mat[seq1]) or (seq1 in self.dist_mat[seq2]),
                self.seqs_to_otus[seq1] == self.seqs_to_otus[seq2],
            )
            for seq1, seq2 in combinations(self.seqs_to_otus.keys(), 2)
        ]
        # account for additional singleton sequences
        classes.extend('tn' for i in range(self.ghost_pairs))
        # convert to a dictionary of counts
        return Counter(classes)

    @property
    def mcc(self):
        return mcc(self.conf_mat)


class OptiFit:
    def __init__(self, ref_otus, query_seqs, query_dist_mat, n_seqs = 0):
        self.ref_otus = ref_otus
        self.query_seqs = query_seqs

        # create merged dist mat
        dist_mat = self.ref_otus.dist_mat.copy()
        for s1 in query_dist_mat:
            dist_mat[s1].update(query_dist_mat[s1])
            for s2 in query_dist_mat[s1]:
                dist_mat[s2].add(s1)
        # initialize OTUs from reference
        seqs_to_otus =  self.ref_otus.seqs_to_otus.copy()
        # seed each query sequence as a singelton OTU
        n_otus = len(self.ref_otus.otus_to_seqs)
        for seq in self.query_seqs:
            n_otus += 1
            seqs_to_otus[seq] = n_otus

        self.fitmap = otuMap(seqs_to_otus = seqs_to_otus,
                             dist_mat = dist_mat,
                             n_seqs = n_seqs)
        self.iterations = list()  # list of dataframes for ggraph

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.fitmap)

    @property
    def mcc(self):
        return self.fitmap.mcc

    def iterate(self):
        return [OptiIter(self.fitmap, seq) for seq in self.query_seqs]

"""
g1 <-
  tbl_graph(nodes = data.frame(name = c("A B C", "D E F", "G **H** I", "J K L"),
                               i = 1:4),
            edges = data.frame(from = c(2, 2, 2),
                               to = c(3, 4, 2),
                               mcc = c(0.97, 0.84, 0.80)) %>%
              mutate(is_loop = from == to))
"""


class OptiIter:
    def __init__(self, fitmap, curr_seq):
        """
        Calculate possible MCCs if the current seq is moved to different OTUs.
        Store the nodes and edges as dictionaries for tidygraph.
        """
        sim_seqs = fitmap.dist_mat[curr_seq]
        options = list()
        self.edges = {'from': [], 'to': [], 'mcc': []}
        for sim_seq in sim_seqs:
                option = OptiOption(fitmap, curr_seq, sim_seq)
                edges['from'].append(option.from)
                edges['to'].append(option.to)
                edges['mcc'].append(potential_map.mcc)
        self.nodes = {'name': [' '.join([s for s in otu
                                         if s != curr_seq else f"**{s}**"]
                                         )
                               for otu in optifit.fitmap.otus_to_seqs.values()],
                      'id': list(fitmap.otus_to_seqs.keys())}


class OptiOption:
    def __init__(self, curr_otu_map, curr_seq, sim_seq):
        """
        Calculate the hypothetical MCC if curr_seq were moved to the OTU containing sim_seq
        """
        self.curr_seq = curr_seq
        self.sim_seq = sim_seq
        self.from = curr_otu_map.seqs_to_otus[curr_seq]
        self.to = curr_otu_map.seqs_to_otus[sim_seq]
        otu_map = deepcopy(curr_otu_map)
        otu_map.seqs_to_otus[curr_seq] = to
        self.mcc = otu_map.mcc


def main():
    ref_otus = opticlust_example()
    query_seqs = {'X', 'Y', 'Z'}
    query_dist_mat = dist_pairs_to_sets({'seq1': ['X', 'X', 'X', 'X', 'Y'],
                                         'seq2': ['Y', 'C', 'G', 'K', 'C']})
    optifit = OptiFit(ref_otus, query_seqs, query_dist_mat, n_seqs = 53)
    print(optifit.mcc)


if __name__ == "__main__":
    main()
