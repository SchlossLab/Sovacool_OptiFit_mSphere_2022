""" Prep a series of graphs for the algorithm diagram figure.
"""
from collections import defaultdict, Counter
from copy import deepcopy
from itertools import combinations
from math import sqrt, factorial#, comb
import numpy
import pandas


def binom_coeff(n, k):
    """
    apparently `from math import comb` isn't working in reticulate right now,
    so we gotta do it this way
    """
    return 0 if k > n else int(factorial(n) / (factorial(k) * factorial(n - k)))


def create_optifit():
    ref_otus = opticlust_example()
    query_seqs = ["W", "X", "Y", "Z"]
    query_dist_mat = dist_pairs_to_sets(
        {
            "seq1": ["X", "X", "X", "X", "Y", "W", "W", "W"],
            "seq2": ["Y", "C", "G", "N", "C", "M", "N", "F"],
        }
    )
    # original opticlust supplement example had 50 sequences, with 15 pairs
    # being within the distance threshold.
    # we now want to fit 4 query sequences to those reference OTUs.
    # thus n_seqs = 50 + 4 = 54.
    optifit = OptiFit(ref_otus, query_seqs, query_dist_mat, n_seqs=54)
    return optifit

def write_iters(optifit_iters, outdir = "figures/algorithm_steps"):
    for i, iter_dir in enumerate(optifit_iters):
        iter_dir['nodes'].to_csv(f"{outdir}/{i}_nodes")
        iter_dir['edges'].to_csv(f"{outdir}/{i}_edges")

def opticlust_example():
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
    assert len(dist_frame["seq1"]) == len(dist_frame["seq2"])
    dist_mat = dist_pairs_to_sets(dist_frame)
    otu_list = [
        {},
        {"I", "D", "B"},
        {"F", "E", "Q"},
        {"C", "G"},
        {"H", "J", "A"},
        {"M", "N"},
        {"P", "L", "O"},
        {"K"},
    ]
    otus = otuMap.from_list(otu_list, dist_mat=dist_mat, n_seqs=50)
    # print(otus)
    # conf_mat = otus.conf_mat(dist_mat)
    # print('mcc current:', mcc(conf_mat),
    #      '\nmcc from correct conf mat:',
    #      mcc({"tp": 14, "tn": 1210, "fp": 0, "fn": 1}),
    #      '\nmcc correct: 0.97')
    return otus


def get_dists():
    return pandas.DataFrame.from_dict(
        {
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
                "X",
                "X",
                "X",
                "X",
                "Y",
                "W",
                "W",
                "W",
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
                "Y",
                "C",
                "G",
                "N",
                "C",
                "M",
                "N",
                "F",
            ],
        }
    )


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


def dist_array_to_dframe(dist_array):
    """
    :param dist_array: pd.DataFrame(dist_array,
                          columns=['A', 'B', 'C', 'D'],
                          index=['A', 'B', 'C', 'D'])
    :return: pd.DataFrame.from_dict({'seq1': ['A', 'A', 'B'],
                                     'seq2': ['B', 'C', 'D']})

    """
    return  # TODO


def otu_list_to_dict(otu_list):
    """
    :param otu_list: should be of the form
        [{seqA, seqB}, {seqC, secD}, {secE}]
        where the index in the list is the OTU ID.
    :return: an otuMap
    """
    return {seq: idx for idx, otu in enumerate(otu_list) for seq in otu}


def format_seq(
    s,
    curr_seq,
    query_seqs,
    color_curr_seq=False,
    do_color=True,
    base_color="#000000",
    ref_color="#D95F02",
    query_color="#1B9E77",
):
    """
    format sequences with bold for the current iteration seq and
    color-code reference and query sequences.
    """
    if s == curr_seq:
        s = f"**{s}**"
        color = query_color if color_curr_seq else base_color
    else:
        color = query_color if s in query_seqs else ref_color

    return f"<span style = 'color:{color};'>{s}</span>" if do_color else s


class otuMap:
    """
    Maps sequences to OTU assignments.
    """

    def __init__(self, seqs_to_otus=None, dist_mat=None, n_seqs=0):
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

    def renumber_otus(self):
        """
        re-number OTU ids so they're continuous
        """
        old_otus = self.otus_to_seqs
        if max(old_otus) != len(old_otus):  # then OTU IDs are not continuous
            self.seqs_to_otus = otu_list_to_dict(
                [{}] + [otu for idx, otu in old_otus.items() if otu]
            )

    @classmethod
    def from_list(cls, otu_list, **kwargs):
        """
        :param otu_list: should be of the form
            [{seqA, seqB}, {seqC, secD}, {secE}]
            where the index in the list is the OTU ID.
        :return: an otuMap
        """
        return cls(seqs_to_otus=otu_list_to_dict(otu_list), **kwargs)

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
        # n_dists = sum([len(dist_mat[s1]) for s1 in dist_mat]) / 2
        return binom_coeff(n_unsim_seqs, 2) + n_unsim_seqs * n_sim_seqs

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
        classes.extend("tn" for i in range(self.ghost_pairs))
        # convert to a dictionary of counts
        return Counter(classes)

    @property
    def mcc(self):
        return mcc(self.conf_mat)

    @property
    def dists_to_array(self):
        """
        :param dist_mat: {'A': {'B','C'}, 'B': {'A','C','D'}, 'C': {'A','B'}, 'D': {'B'}}
        :return: pd.DataFrame(dist_array,
                              columns=['A', 'B', 'C', 'D'],
                              index=['A', 'B', 'C', 'D'])
        """
        dist_sets = self.dist_mat
        seqs = {seq2 for seq1 in dist_sets for seq2 in dist_sets[seq1]}
        seqs.update({seq for seq in dist_sets})
        seqs = list(sorted(seqs))
        len_seqs = len(seqs)
        dist_array = numpy.zeros(len_seqs, len_seqs)
        for i, seqi in enumerate(seqs):
            for j, seqj in enumerate(seqs):
                if i == j or seqi in dist_sets[seqj] or seqj in dist_sets[seqi]:
                    dist_array[i][j] = 1
                    dist_array[j][i] = 1
        return pandas.DataFrame(dist_array, columns=seqs, index=seqs)


class OptiFit:
    def __init__(self, ref_otus, query_seqs, query_dist_mat, n_seqs=0):
        self.ref_otus = ref_otus
        self.query_seqs = query_seqs

        # create merged dist mat
        dist_mat = self.ref_otus.dist_mat.copy()
        for s1 in query_dist_mat:
            dist_mat[s1].update(query_dist_mat[s1])
            for s2 in query_dist_mat[s1]:
                dist_mat[s2].add(s1)
        # initialize OTUs from reference
        seqs_to_otus = self.ref_otus.seqs_to_otus.copy()
        # seed each query sequence as a singelton OTU
        n_otus = len(self.ref_otus.otus_to_seqs)
        for seq in self.query_seqs:
            n_otus += 1
            seqs_to_otus[seq] = n_otus

        self.fitmap = otuMap(
            seqs_to_otus=seqs_to_otus, dist_mat=dist_mat, n_seqs=n_seqs
        )

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.fitmap)

    @property
    def mcc(self):
        return self.fitmap.mcc

    @property
    def iterate(self):
        iterations = list()
        curr_fitmap = self.fitmap
        prev_mcc = 0
        while not numpy.isclose(prev_mcc, self.mcc):
            prev_mcc = self.mcc
            for seq in self.query_seqs:
                if seq in self.fitmap.dist_mat:
                    iteration = OptiIter(curr_fitmap, seq, self.query_seqs)
                    iterations.append(iteration.to_dict)
                    curr_fitmap = iteration.best_map
        return 

    @property   
    def iterate_obj(self):
        iterations = list()
        curr_fitmap = self.fitmap
        prev_mcc = 0
        while not numpy.isclose(prev_mcc, self.mcc):
            prev_mcc = self.mcc
            for seq in self.query_seqs:
                if seq in self.fitmap.dist_mat:
                    iteration = OptiIter(curr_fitmap, seq, self.query_seqs)
                    iterations.append(iteration)
                    curr_fitmap = iteration.best_map
        return iterations


class OptiIter:
    def __init__(self, curr_fitmap, curr_seq, query_seqs):
        """
        Calculate possible MCCs if the current seq is moved to different OTUs.
        Store the nodes and edges as dictionaries for tidygraph.
        """
        sim_seqs = curr_fitmap.dist_mat[curr_seq]
        sim_seqs.add(curr_seq)
        options = list()
        edges = {"from": [], "to": [], "mcc": []}
        for sim_seq in sim_seqs:
            option = OptiOption(curr_fitmap, curr_seq, sim_seq)
            options.append(option)
            edges["from"].append(option.from_otu)
            edges["to"].append(option.to_otu)
            edges["mcc"].append(option.mcc)

        self.edges = pandas.DataFrame.from_dict(edges)
        self.nodes = pandas.DataFrame.from_dict(
            {
                "name": [
                    " ".join(
                        [
                            format_seq(s, curr_seq, query_seqs, color_curr_seq=True)
                            for s in sorted(otu)
                        ]
                    )
                    for otu in curr_fitmap.otus_to_seqs.values()
                ],
                "id": list(curr_fitmap.otus_to_seqs.keys()),
            }
        )
        best_option = max(options)
        # make sure OTU numbers are continuous
        self.best_map = best_option.fitmap

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__.items())

    @property
    def to_dict(self):
        return {"nodes": self.nodes, "edges": self.edges}

class OptiOption:
    def __init__(self, curr_fitmap, curr_seq, sim_seq):
        """
        Calculate the hypothetical MCC if curr_seq were moved to the OTU containing sim_seq
        """
        self.curr_seq = curr_seq
        self.sim_seq = sim_seq
        self.from_otu = curr_fitmap.seqs_to_otus[curr_seq]
        self.to_otu = curr_fitmap.seqs_to_otus[sim_seq]

        self.fitmap = deepcopy(curr_fitmap)
        self.fitmap.seqs_to_otus[curr_seq] = self.to_otu
        self.fitmap.renumber_otus()
        self.mcc = self.fitmap.mcc

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__.items())

    def __eq__(self, other):
        return self.mcc == other.mcc

    def __ge__(self, other):
        return self.mcc >= other.mcc

    def __gt__(self, other):
        return self.mcc > other.mcc
