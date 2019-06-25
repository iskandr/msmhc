def write_fasta(filename, name_to_seq_dict, maxwidth=60):
    with open(filename, "w") as f:
        for (name, seq) in name_to_seq_dict.items():
            assert ">" not in name
            if isinstance(name, tuple):
                name = " ".join(name)
            f.write(">%s\n" % name)
            for i in range(len(seq) // maxwidth + 1):
                subseq = seq[i * maxwidth:(i + 1) * maxwidth]
                if len(subseq) > 0:
                    f.write(subseq)
            f.write("\n")