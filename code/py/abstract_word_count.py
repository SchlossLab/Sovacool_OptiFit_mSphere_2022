def word_count(infilename, starter, stopper):
    with open(infilename, "r") as infile:
        words = []
        line = next(infile)
        while starter != line:
            line = next(infile)
        line = next(infile)  # make sure the first parsed line is not the starter
        while stopper != line:
            words += [word for word in line.strip().split()]
            line = next(infile)
    return len(words)


def check_wc(section_name, num_words, word_limit):
    if num_words > word_limit:
        raise ValueError(
            f"The {section_name} section is {num_words} words. You need to cut {num_words - word_limit} words."
        )


def main(src_filename, log_filename):
    with open(log_filename, "w") as outfile:
        outfile.write("section\tword_count\n")
        for section, word_limit, starter, stopper in zip(
            ["abstract", "importance"],
            [250, 150],
            ["## Abstract\n", "### Importance\n"],
            ["### Importance\n", "\\newpage\n"],
        ):
            wc = word_count(src_filename, starter, stopper)
            check_wc(section, wc, word_limit)
            outfile.write(f"{section}\t{wc}\n")


if __name__ == "__main__":
    if "snakemake" in locals() or "snakemake" in globals():
        main(snakemake.input.src, snakemake.output.txt)
    else:
        main("paper/paper.Rmd", "log/count_words_abstract.log")
