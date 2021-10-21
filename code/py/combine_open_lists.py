""" Combine two mothur list files into one with renamed OTUs """
import sys


def main(infilename1, infilename2, outfilename):
    mlist = mothurList.from_list_file(infilename1)
    mlist.combine(mothurList.from_list_file(infilename2))
    mlist.write(outfilename)


class mothurList:
    def __init__(self, label, otu_assigns):
        self.label = label
        self.otu_assigns = otu_assigns

    def __eq__(self, other):
        return all([self.__class__ ==  other.__class__,
                    self.label == other.label,
                    self.otu_assigns == other.otu_assigns])

    @classmethod
    def from_list_file(cls, list_filename):
        with open(list_filename, "r") as listfile:
            lines = [line.strip().split('\t') for line in listfile]
            len_lines = len(lines)
            dat_line = []
            if len_lines == 1: # for when there's no header line
                dat_line = lines[0]
            elif len_lines == 2: # ditch the header line
                dat_line == lines[1]
            else: # this shouldn't happen
                raise ValueError(f"List file contains {len_lines} lines, but only 1 or 2 were expected.")
            list_label = dat_line[0]
            otu_assigns = dat_line[2:]
        return cls(list_label, otu_assigns)

    @property
    def num_otus(self):
        return len(self.otu_assigns)

    @property
    def get_new_otu_names(self):
        return ["OTU_" + str(i) for i in range(1, self.num_otus + 1)]

    def combine(self, other):
        assert other.__class__ == self.__class__
        assert self.label == other.label
        self.otu_assigns += other.otu_assigns

    def write(self, list_filename, use_header=False):
        with open(list_filename, "w") as listfile:
            if use_header:
                listfile.write(
                    "\t".join(["label", "numOTUs"] + self.get_new_otu_names) + "\n"
                )
            listfile.write(f"{self.label}\t{str(self.num_otus)}")
            for otu in self.otu_assigns:
                listfile.write(f"\t{otu}")
            listfile.write('\n')


if __name__ == "__main__":
    if "snakemake" in locals() or "snakemake" in globals():
        main(snakemake.input.list_closed, snakemake.input.list_denovo, snakemake.output.list)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
