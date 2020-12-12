
def parse_accnos(filename):
    with open(filename, 'r') as file:
        seqs = {line.strip() for line in file}
    return seqs


def main():
    ref = parse_accnos(snakemake.input.ref)
    query = parse_accnos(snakemake.input.query)
    all = parse_accnos(snakemake.intput.all)

    assert not ref.intersection(query)
    assert all - ref.intersection(all)
    assert all - query.intersection(all)
    assert ref + query == all

    with open(snakemake.output.txt, 'w') as file:
        file.write('check_split_passed\n')
        file.write(True)

main()