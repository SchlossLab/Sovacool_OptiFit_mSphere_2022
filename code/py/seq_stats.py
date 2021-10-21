from split_weighted_subsample import MetaSeq, SeqList


def main(fasta_file, count_file, dist_file, output_file, dissim_thresh=0.03):

        all_seqs = SeqList.from_files(
            fasta_file,
            count_file,
            dist_file,
            threshold=dissim_thresh,
        )
        with open(output_file, 'w') as file:
            file.write('seq\tabundance\tsimilarity\n')
            for seq in all_seqs.seqs:
                file.write(f"{seq.seq_id}\t{seq.abs_abun}\t{seq.sum_sim}\n")

if __name__ == "__main__":
    main(
        snakemake.input.fasta,
        snakemake.input.count,
        snakemake.input.dist,
        snakemake.output.tsv,
        dissim_thresh = snakemake.params.dissim_thresh
    )
