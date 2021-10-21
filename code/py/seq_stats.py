from split_weighted_subsample import MetaSeq, SeqList


def main(fasta_file, count_file, dist_file, output_file, dissim_thresh=0.03):

    all_seqs = SeqList.from_files(
        fasta_file,
        count_file,
        dist_file,
        threshold=dissim_thresh,
    )
    with open(output_file, "w") as file:
        file.write("seq\trel_abundance\trel_similarity\n")
        for seq_id, rel_abun, rel_sim in zip(
            all_seqs.ids, all_seqs.rel_abuns, all_seqs.rel_sims
        ):
            file.write(f"{seq_id}\t{rel_abun}\t{rel_sim}\n")


if __name__ == "__main__":
    main(
        snakemake.input.fasta,
        snakemake.input.count,
        snakemake.input.dist,
        snakemake.output.tsv,
        dissim_thresh=snakemake.params.dissim_thresh,
    )
