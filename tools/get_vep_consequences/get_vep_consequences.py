import sys
import json
import gzip

def get_vep_consequences(config_file, input_ld_file, output_file):
    
    with open(config_file) as f:
        config = json.load(f)["get_vep_consequences"]

    vars_to_keep = {}

    with open(input_ld_file) as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            chrom = data[0]
            pos = data[1]
            locus = data[2]
            ld = data[3]
            vars_to_keep[(chrom, pos)] = (locus, ld)

    i = 0
    with open(output_file, "w") as w:
        w.write("chrom\tpos\tlocus\tld\tref\talt\tconsequence\tgene\ttranscript\thgnc\n")
        with gzip.open(config["cadd_file"]) as f:
            f.readline()
            f.readline()
            for line in f:
                i += 1
                if i % 1000000 == 0:
                    print i
                data = line.strip().split()
                chrom = data[0]
                pos = data[1]
                ref = data[2]
                alt = data[3]
                consequence = data[7]
                gene = data[18]
                transcript = data[19]
                hgnc = data[20]
                if ((chrom, pos)) not in vars_to_keep:
                    continue

                locus = vars_to_keep[(chrom, pos)][0]
                ld = vars_to_keep[(chrom, pos)][1]
                w.write("\t".join([chrom, pos, locus, ld, ref, alt, consequence, gene, transcript, hgnc]) + "\n")

config_file = sys.argv[1]
input_ld_file = sys.argv[2]
output_file = sys.argv[3]

get_vep_consequences(config_file, input_ld_file, output_file)
