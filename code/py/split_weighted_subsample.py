#!/usr/local/bin/python3
import Bio.SeqIO
import numpy as np

class MetaSeq:
    def __init__(self, seq_id, avg_abun, avg_dist):
        self.seq_id = seq_id
        self.avg_abun = avg_abun
        self.avg_dist = avg_dist

    @property
    def avg_sim:
        return 1 - avg_dist

def main():
    with open(input.fasta, 'r') as fasta_file:
        seq_dict = {seq_record.id: MetaSeq(seq_record.id, np.nan, np.nan)
                    for seq_record in Bio.SeqIO.read(fasta_file, 'fasta')}
    for src, dest in [[input.fasta, output.fasta],
                      [input.count, output.count],
                      [input.dist, output.dist]
                      ]:
        shutil.copyfile(src, dest)
    numpy.random.seed(wildcards.seed)
    with open(input.count, 'r') as count_file:
        line = next(count_file)
        for line in count_file:
            line = line.strip().split('\t')
            seq_id = line[0]
            seq_dict[seq_id].avg_abun = np.mean(float(count) for count in line[1:])
    with open(input.dist, 'r') as dist_file:
        line = next(count_file)
        

main()