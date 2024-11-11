import sys
from Bio import SeqIO
import pysam
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

sns.set_theme(style="whitegrid")
colors = [sns.color_palette("bright")[i] for i in [0, 1, 2, 6]]


def canonical(kmer):
    kmer = kmer
    kmer_rc = kmer.reverse_complement()
    kmer = str(kmer)
    kmer_rc = str(kmer_rc)
    return kmer if kmer < kmer_rc else kmer_rc


def load(fp):
    kmers = {}
    for line in open(fp):
        kmer, w = line.strip("\n").split("\t")
        kmers[kmer] = int(w)
    return kmers


def load_bg(fp):
    data = []
    for line in open(fp):
        idx, s, e, v = line.strip("\n").split("\t")
        s, e = int(s), int(e)
        v = float(v)
        for p in range(s, e):
            data.append(v)
    return data


def load_pileup(fp, chrom, s, e):
    print(chrom, s, e)
    data = [0] * (e - s + 1)
    bam = pysam.AlignmentFile(fp, "rb")
    i = 0
    for pileupcolumn in bam.pileup(
        chrom,
        s - 1,
        e,
        truncate=True,
        flag_filter=0,
        ignore_orphans=False,
        min_base_quality=0,
    ):
        skipcount = 0
        for pileupread in pileupcolumn.pileups:
            skipcount += pileupread.is_del or pileupread.is_refskip
        data[i] = len(pileupcolumn.pileups) - skipcount
        i += 1
    return data


def parse_region(r):
    chrom = r.split(":")[0]
    s, e = r.split(":")[1].split("-")
    return chrom, int(s), int(e)


def main():
    fa_fn = sys.argv[1]

    record = next(SeqIO.parse(fa_fn, "fasta"))
    idx = record.id
    chrom, s, e = parse_region(idx)
    seq = record.seq

    KDBS = [load(sys.argv[2]), load(sys.argv[3]), load(sys.argv[4]), load(sys.argv[5])]
    BGS = [load_bg(sys.argv[6]), load_bg(sys.argv[7])]
    BAMS = [
        load_pileup(sys.argv[8], chrom, s, e),
        load_pileup(sys.argv[9], chrom, s, e),
        load_pileup(sys.argv[10], chrom, s, e),
    ]

    n = len(KDBS)
    k = len(next(iter(KDBS[0].keys())))

    labels = ["Ref", "Ctrl", "Case1", "Case2"]
    data = []
    for p in range(0, len(seq) - k + 1):
        kmer = canonical(seq[p : p + k].upper())
        for i in range(n):
            w = KDBS[i][kmer] if kmer in KDBS[i] else 0
            for _ in range(w + 1):
                data.append([labels[i], p])
    df = pd.DataFrame(data, columns=["Type", "p"])

    labels2 = labels[2:]
    n2 = len(labels2)
    data2 = []
    for p in range(0, len(seq) - k + 1):
        for i in range(n2):
            w = BGS[i][p]
            data2.append([labels2[i], p, w])
    df2 = pd.DataFrame(data2, columns=["Type", "p", "ratio"])

    labels3 = labels[1:]
    n3 = len(labels3)
    data3 = []
    for p in range(0, len(seq) - k + 1):
        for i in range(n3):
            w = BAMS[i][p]
            data3.append([labels3[i], p, w])
    df3 = pd.DataFrame(data3, columns=["Type", "p", "pileup"])

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

    sns.histplot(
        df,
        x="p",
        hue="Type",
        bins=len(seq) - k,  #  range(0, len(seq) - k + 1),
        discrete=True,
        element="step",
        # kde=True,
        palette=colors,
        ax=ax1,
    )

    sns.barplot(df2, x="p", y="ratio", hue="Type", palette=colors[2:], ax=ax2)

    sns.barplot(df3, x="p", y="pileup", hue="Type", palette=colors[1:], ax=ax3)

    ax1.set_xticks(range(len(seq)), list(seq))
    ax2.xaxis.tick_top()
    plt.suptitle(idx)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.10)
    plt.show()


if __name__ == "__main__":
    main()
