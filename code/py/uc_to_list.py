import collections
import sys


def uc_to_list(clustered_filename, list_filename):
    clusters = collections.defaultdict(set)
    with open(clustered_filename, "r") as clustered_file:
        for line in clustered_file:
            line = line.split("\t")
            # https://drive5.com/usearch/manual/opt_uc.html
            record_type = line[0]
            cluster_num = line[1]
            query_label = line[8].split(";size=")[0]
            if record_type in {"H", "S"}:
                clusters[cluster_num].add(query_label)
    with open(list_filename, "w") as list_file:
        list_file.write("userLabel\t" + str(len(clusters)))
        for cluster_id, seqs in clusters.items():
            list_file.write("\t" + ",".join(sorted(seqs)))
        list_file.write("\n")


if __name__ == "__main__":
    if "snakemake" in globals() or "snakemake" in locals():
        uc_to_list(snakemake.input.uc, snakemake.output.list)
    else:
        uc_to_list(sys.argv[1], sys.argv[2])
