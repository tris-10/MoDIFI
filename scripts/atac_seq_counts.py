import pysam
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="The script counts up reads per peak")
    parser.add_argument("input_intervals")
    parser.add_argument("input_bam")
    parser.add_argument("output_file")
    parser.add_argument("sample_name")
    parser.add_argument("--minQ", type=int, default=0, help="Minimum mapping quality to include a read (default: 0)")
    args = parser.parse_args()

    with pysam.AlignmentFile(args.input_bam, 'rb') as in_sam, open(args.input_intervals, "r") as ipf, open(args.output_file, "w") as opf:
        opf.write("{0}\n".format(args.sample_name))
        for line in ipf:
            items = line.strip().split("\t")
            read_set = set()
            for read in in_sam.fetch(items[0], int(items[1]), int(items[2])):
                if read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < args.minQ:
                    continue
                read_set.add(read.query_name)
            #items.append(str(len(read_set)))
            #opf.write("\t".join(items) + "\n")
            opf.write("{0}\n".format(len(read_set)))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)
