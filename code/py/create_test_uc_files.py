def main(uc1_filename, uc2_filename, list_filename):
    otus1 = {1: {'a', 'b'},
             2: {'c', 'd'},
             3: {'e'}
             }
    otus2 = {1: {'f', 'g'},
             2: {'h'}
             }
    for filename, otus in ((uc1_filename, otus1),
                           (uc2_filename, otus2)):
        with open(filename, 'w') as uc_file:
            for otu_id, seqs in otus1.items():
                for seq_id in seqs:
                    uc_file.write(f"H\t{otu_id}\t1\t100\t+\t-\t-\t=\t{seq_id}\t{otu_id}\n")
    combined = [','.join(seqs) for otus in [otus1, otus2] for seqs in otus.values()]
    with open(list_filename, 'w') as listfile:
        listfile.write(f"userLabel\t{str(len(combined))}")
        for otu in combined:
            listfile.write(f"\t{otu}")
        listfile.write('\n')

if __name__ == "__main__":
    main('code/tests/data/closed.uc',
         'code/tests/data/denovo.uc',
         'code/tests/data/open.list')