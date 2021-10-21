"""
Usage:
    ./remove_hyphens_dist.py <input.dist> <output.ng.dist>
"""
import sys


def replace_chars(string, bad_chars=[".", "-"], replacement="_"):
    for bad_char in bad_chars:
        string = string.replace(bad_char, replacement)
    return string


def main(infilename, outfilename):
    with open(outfilename, "w") as outfile:
        with open(infilename, "r") as infile:
            for line in infile:
                # mothur dist files aredelmited by spaces
                line_split = line.split(" ")
                # first & second columns are seq names
                for i in range(2):
                    line_split[i] = replace_chars(line_split[i])
                outfile.write(" ".join(line_split))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
