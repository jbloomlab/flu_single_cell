"""Make FASTA and GTF files for vRNAs from Genbank plasmid maps."""


import re
import os
import glob
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import Bio.SeqFeature
import BCBio.GFF


def main():
    """Main body of script."""
    u12 = re.compile('AGC[GA]AAAGCAGG') # U12 with position 4 polymorphism
    u13 = 'CCTTGTTTCTACT' # reverse of U13
    segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

    for syntype in ['', '-syn']:
        vrnas = []
        for seg in segs:
            fplasmid = glob.glob('*_pHW18*-{0}{1}.gb'.format(seg, syntype))
            assert len(fplasmid) == 1, "Didn't find 1 map for {0}{1}".format(
                    seg, syntype)
            fplasmid = fplasmid[0]
            plasmid = str(Bio.SeqIO.read(fplasmid, 'genbank').seq)
            assert plasmid.count(u13) == 1, "Not exactly one U13 for {0}{1}".format(
                    seg, syntype)
            u12matches = list(u12.finditer(plasmid))
            assert len(u12matches) == 1, "Not exactly one U12 for {0}{1}".format(
                    seg, syntype)
            vrna = plasmid[u12matches[0].start(0) : plasmid.index(u13) + len(u13)]
            print("vRNA of length {0} for {1}{2}".format(len(vrna), seg, syntype))
            vrna = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(vrna), 
                    id='flu' + seg + syntype,
                    description='{0} vRNA from plasmid {1}'.format(seg, fplasmid),
                    features=[Bio.SeqFeature.SeqFeature(
                        Bio.SeqFeature.FeatureLocation(0, len(vrna)),
                        type='exon',
                        strand=1,
                        qualifiers={
                            "source":"JesseBloom",
                            "gene_id":'flu' + seg,
                            "transcript_id":'flu' + seg,
                            "gene_biotype":"vRNA",
                            },
                        )],
                    )
            vrnas.append(vrna)
        fastafile = 'flu-wsn{0}.fasta'.format(syntype)
        print("Writing the {0} vRNAs to {1}".format(len(vrnas), fastafile))
        Bio.SeqIO.write(vrnas, fastafile, 'fasta')
        gtffile = os.path.splitext(fastafile)[0] + '.gtf'
        print("Writing the vRNA annotations to {0}\n".format(gtffile))
        with open(gtffile, 'w') as f:
            BCBio.GFF.write(vrnas, f)
        # Hacky conversion of GFF3 file to GTF, changing 9th column to
        # delimit qualifiers with spaces / quotes rather than equals sign.
        # Conversion is probably not robust to all qualifiers, but works here.
        with open(gtffile) as f:
            lines = f.readlines()
        newlines = []
        for line in lines:
            if line[0] == '#':
                newlines.append(line)
            else:
                entries = line.strip().split('\t')
                assert len(entries) == 9, str(len(entries)) + '\n' + line
                newqualifiers = []
                for qualifier in entries[-1].split(';'):
                    (key, value) = qualifier.split('=')
                    newqualifiers.append('{0} "{1}"'.format(key, value))
                entries[-1] = '; '.join(newqualifiers)
                newlines.append('\t'.join(entries) + '\n')
        with open(gtffile, 'w') as f:
            f.write(''.join(newlines))


# run the script
if __name__ == '__main__':
    main()

