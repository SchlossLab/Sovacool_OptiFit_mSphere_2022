#!/usr/local/bin/python3
""" Select weighted subsets of sequences to be used as references and samples for OptiFit """
from collections import defaultdict
import numpy as np
import shutil


def main(
    seed,
    fasta_file,
    count_file,
    dist_file,
    dissim_thresh,
    ref_frac,
    query_frac,
    ref_weight,
    ref_accnos_file,
    query_accnos_file,
    all_accnos_file,
):

    np.random.seed(int(seed))

    all_seqs = SeqList.from_files(
        fasta_file,
        count_file,
        dist_file,
        threshold=dissim_thresh,
    )
    num_all_seqs = len(all_seqs)
    ref_frac = float(ref_frac)
    query_frac = float(query_frac)
    assert ref_frac + query_frac <= 1
    query_size = round_subset_size(query_frac, num_all_seqs)
    ref_size = round_subset_size(ref_frac, num_all_seqs)

    ref_list = all_seqs.get_sample(ref_size, ref_weight)
    assert check_subsample(ref_frac, len(ref_list), num_all_seqs)
    ref_list.write_ids(ref_accnos_file)

    remaining_seqs = SeqList.set_diff(all_seqs, ref_list)
    query_list = (
        remaining_seqs
        if ref_frac + query_frac == 1
        else remaining_seqs.get_sample(query_size, "simple")
    )
    assert check_subsample(query_frac, len(query_list), num_all_seqs)
    query_list.write_ids(query_accnos_file)

    all_seqs = SeqList(
        [seq for seqlist in [ref_list, query_list] for seq in seqlist.seqs]
    )
    all_seqs.write_ids(all_accnos_file)


def round_subset_size(fraction, total_size):
    return int(round(fraction * total_size, 0))


def check_subsample(fraction, subset_size, total_size):
    return np.isclose(subset_size, fraction * total_size, rtol=0.05)


class MetaSeq:
    def __init__(self, seq_id, abs_abun, sum_sim):
        self.seq_id = seq_id
        self.abs_abun = abs_abun
        self.sum_sim = sum_sim

    def __repr__(self):
        return f"{self.__class__}({self.__dict__})"

    def __hash__(self):
        return hash(self.__repr__())


class SeqList:
    def __init__(self, seqs):
        self.seqs = seqs

    def __len__(self):
        return len(self.seqs)

    def __repr__(self):
        return f"{self.__class}({self.__dict__})"

    @property
    def ids(self):
        return [seq.seq_id for seq in self.seqs]

    @property
    def rel_abuns(self):
        total_abun = sum(seq.abs_abun for seq in self.seqs)
        return [seq.abs_abun / total_abun for seq in self.seqs]

    @property
    def rel_sims(self):
        total_sim = sum(seq.sum_sim for seq in self.seqs)
        return [seq.sum_sim / total_sim for seq in self.seqs]

    @classmethod
    def from_files(cls, fasta_fn, count_fn, dist_fn, threshold=0.03):
        with open(dist_fn, "r") as dist_file:
            sum_sims = defaultdict(int)
            for line in dist_file:
                line = line.strip().split()
                seq_id1 = line[0]
                seq_id2 = line[1]
                is_similar = int(float(line[2]) < threshold)
                sum_sims[seq_id1] += is_similar
                sum_sims[seq_id2] += is_similar
        with open(count_fn, "r") as count_file:
            line = next(count_file)  # toss out the header
            seq_list = [
                MetaSeq(
                    seq_id=line.strip().split()[0],
                    abs_abun=int(line.strip().split()[1]),
                    sum_sim=sum_sims[line.strip().split()[0]],
                )
                for line in count_file
            ]
        return cls(seq_list)

    @classmethod
    def set_diff(cls, lhs, rhs):
        return cls(list(set(lhs.seqs) - set(rhs.seqs)))

    def get_sample(self, query_size, weight_method):
        random_weight_probs = {
            "simple": None,
            "abundance": self.rel_abuns,
            "distance": self.rel_sims,
        }
        sample_seqs = np.random.choice(
            self.seqs,
            replace=False,
            size=query_size,
            p=random_weight_probs[weight_method],
        )
        return SeqList(sample_seqs)

    def write_ids(self, output_fn):
        with open(output_fn, "w") as outfile:
            for seq_id in self.ids:
                outfile.write(f"{seq_id}\n")


if __name__ == "__main__":
    main(
        snakemake.wildcards.seed,
        snakemake.input.fasta,
        snakemake.input.count,
        snakemake.input.dist,
        snakemake.params.dissim_thresh,
        snakemake.wildcards.ref_frac,
        snakemake.wildcards.sample_frac,
        snakemake.wildcards.ref_weight,
        snakemake.output.ref_accnos,
        snakemake.output.query_accnos,
        snakemake.output.all_accnos,
    )
