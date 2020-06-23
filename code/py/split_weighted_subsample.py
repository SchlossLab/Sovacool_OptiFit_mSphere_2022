#!/usr/local/bin/python3
""" Select weighted subsets of sequences to be used as references and samples for OptiFit """
from collections import defaultdict
import numpy as np
import pandas as pd
import shutil


def main():
    np.random.seed(int(snakemake.wildcards.seed))

    all_seqs = SeqList.from_files(
        snakemake.input.fasta, snakemake.input.count, snakemake.input.dist
    )
    print("all_seqs read")

    sample_size = int(round(float(snakemake.wildcards.sample_frac) * len(all_seqs), 0))
    print(sample_size)
    ref_size = int(round(float(snakemake.wildcards.ref_frac) * len(all_seqs), 0))
    ref_list = all_seqs.get_sample(ref_size, snakemake.wildcards.ref_weight)
    print("ref_list made")
    ref_list.write_ids(snakemake.output.ref_accnos)
    print("ref_list written")

    remaining_seqs = SeqList.set_diff(all_seqs, ref_list)
    print("remaining seqs extracted")
    sample_list = remaining_seqs.get_sample(sample_size, "simple")
    print("sample seqs obtained")
    sample_list.write_ids(snakemake.output.sample_accnos)
    print("sample accnos written")


class MetaSeq:
    def __init__(self, seq_id, abs_abun, sum_dist):
        self.seq_id = seq_id
        self.abs_abun = abs_abun
        self.sum_dist = sum_dist

    def __repr__(self):
        return f"{self.__class}({self.__dict__})"

    @property
    def sum_sim(self):
        return 1 - sum_dist


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
    def rel_dists(self):
        total_dist = sum(seq.sum_dist for seq in self.seqs)
        return [seq.sum_dist / total_dist for seq in self.seqs]

    @property
    def rel_sims(self):
        return [dist - 1 for dist in self.rel_dists]

    @classmethod
    def from_files(cls, fasta_fn, count_fn, dist_fn):
        with open(dist_fn, "r") as dist_file:
            line = next(dist_file)
            sum_dists = defaultdict(int)
            for line in dist_file:
                line = line.strip().split()
                seq_id1 = line[0]
                seq_id2 = line[1]
                dist = float(line[2])
                sum_dists[seq_id1] += dist
                sum_dists[seq_id2] += dist
        with open(count_fn, "r") as count_file:
            line = next(count_file)  # toss out the header
            seq_dict = {
                line.strip().split()[0]: MetaSeq(
                    seq_id=line.strip().split()[0],
                    abs_abun=int(line.strip().split()[1]),
                    sum_dist=sum_dists[line.strip().split()[0]],
                )
                for line in count_file
            }
        return cls(list(sorted(seq_dict.values(), key=lambda seq: seq.seq_id)))

    @classmethod
    def set_diff(cls, lhs, rhs):
        return cls([seq for seq in lhs.seqs if seq.seq_id not in set(rhs.ids)])

    def get_sample(self, sample_size, weight_method):
        random_weight_probs = {
            "simple": None,
            "abundance": self.rel_abuns,
            "distance": self.rel_dists,
        }
        sample_seqs = np.random.choice(
            self.seqs,
            replace=False,
            size=sample_size,
            p=random_weight_probs[weight_method],
        )
        return SeqList(sample_seqs)

    def write_ids(self, output_fn):
        with open(output_fn, "w") as outfile:
            for seq_id in self.ids:
                outfile.write(f"{seq_id}\n")


main()
