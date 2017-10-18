from __future__ import print_function
import re, argparse
import pandas as pd
import numpy as np

COLUMN_MODSEQ = "Modified sequence"
COLUMN_PROTEINS = "Proteins"
COLUMN_MODPROB = "Acetyl (K) Probabilities"
MODTYPE = "(ac)"

COLUMN_SITES = "Sites"
COLUMN_PROB = "Probability"

def run_sites(fn_fasta, fn_evidence, fn_output, conventional_counting):
    fa = Fasta()
    fa.set_file(fn_fasta)
    fa.parse_fasta()

    df = pd.read_csv(fn_evidence, sep='\t', low_memory=False)
    df.dropna(axis=0, how="all", inplace=True)

    df["pepseq"] = df[COLUMN_MODSEQ].apply(lambda aaseq: aaseq.replace("_", "").replace(MODTYPE, ""))
    df["pepseq"] = df["pepseq"].apply(remove_modification_in_parentheses)
    df["start_pos"] = df.apply(get_start_position_of_sequence_proteinGroups, args=(fa,), axis=1)

    if conventional_counting > 0:
        df["start_pos"] = df["start_pos"].apply(start_counting_from_num, args=(conventional_counting, ))

    COLUMN_MODSEQ_index = df.columns.tolist().index(COLUMN_MODSEQ) + 1  # due to index
    start_pos_index = df.columns.tolist().index("start_pos") + 1  # due to index
    Sites_pos_within_pep_list = []
    for row in df.itertuples():
        mod_seq = row[COLUMN_MODSEQ_index]
        start_pos_list = row[start_pos_index].split(";")  # string for every protein in proteinGroups the start position of the peptide
        sites_list = parse_sites_within_pepseq(mod_seq, [])

        sites_per_row = ""
        for protein_start_pos in start_pos_list:
            try:
                protein_start_pos = int(float(protein_start_pos))
            except ValueError:
                sites_per_row += "(nan)" + ";"
                continue
            sites_per_protein = "(" + "+".join([str(site + protein_start_pos) for site in sites_list]) + ")"
            sites_per_row += sites_per_protein + ";"
        Sites_pos_within_pep_list.append(sites_per_row[:-1])
    df["Positions within Proteins"] = Sites_pos_within_pep_list
    df.to_csv(fn_output, sep='\t', header=True, index=False)


class Fasta(object):
    """
    not instantiated with file, in order to add multiple files, by setting and parsing iteratively.
    search for an AccessionNumber, Organism, DataBase, sequence, or partial sequence
    get proteotypic peptides for Organism
    write fasta of subset (e.g. AccessionNumbers, Organism, etc.)
    """

    def __init__(self):
        """
        an2aaseq_dict: key=AccessionNumber, val=AAseq
        org2an_dict: key=Organism, val=ListOfAN
        db2an_dict: key=DataBase, val=ListOfAN
        an2header_dict: key=AccessionNumber, val=header # original header
        :return: None
        """
        organism_regexes = [r"^>([s|t][p|r])\|(.*)\|.*?OS=(.*?) (?:[A-Z][A-Z]=)?"]
        self.my_regex = re.compile('|'.join(organism_regexes))
        self.an2aaseq_dict = {}
        self.an2header_dict = {}

    def set_file(self, fn):
        self.fasta_file = fn

    def get_file(self):
        return self.fasta_file

    def parse_fasta(self, unique=True):
        """
        parse fasta-file by entering infos into dicts
        print summary information
        number of entries total, non-redundant aaseqs, entries per DataBase, entries per Organism
        :return: Str
        """
        for entry in self.yield_entry():
            header, aaseq = entry
            match = re.search(self.my_regex, header)
            if match:  # at least one regex found something
                match_groups = filter(None, match.groups())
                if len(match_groups) == 3:
                    database, an, TaxName = match_groups
                    an = an.strip()
                else:
                    print("Nothing could be extracted from the following fasta-entry: ", "\n", header, aaseq)
                    raise StopIteration
            else:
                print("Nothing could be extracted from the following fasta-entry: ", "\n", header, aaseq)
                raise StopIteration
            # now fill dicts
            if not self.an2aaseq_dict.has_key(an):
                self.an2aaseq_dict[an] = aaseq
            else:  # AccessionNumbers should be unique identifiers
                # self.an2aaseq_dict[an].append(aaseq)
                if unique:
                    print("AccessionNumbers should be unique identifiers!")
                    print(an, "\n", self.an2aaseq_dict[an], "\n", header)
                    raise StopIteration
                else:
                    print("AN double entry: {}".format(an))

            if not self.an2header_dict.has_key(an):
                self.an2header_dict[an] = header
            else:
                if unique:
                    print("AccessionNumbers should be unique identifiers!")
                    print(an, "\n", self.an2aaseq_dict[an], "\n", aaseq)
                    raise StopIteration
                else:
                    pass

    def yield_entry(self):
        """
        generator that yields one entry of a fasta-file at a time,
        as tuple (header, aaseq)
        :return: tuple(Str, Str)
        """
        with open(self.get_file(), "r") as fh:
            aaseq = ""
            header = ""
            did_first = False
            for line in fh:
                if line[0] == ">":
                    if did_first == True:
                        if len(aaseq) > 0:
                            yield (header, aaseq)
                    else:
                        did_first = True
                    header = line.rstrip()
                    aaseq = ""
                else:
                    aaseq += line.strip().upper()
            if len(aaseq) > 0:
                yield (header, aaseq)

    def get_aaseq_from_an(self, an):
        """
        return AAseq (protein sequence) of given AccessionNumber as String
        :param an: String
        :return: String
        """
        return self.an2aaseq_dict[an]


def remove_modification_in_parentheses(aaseq):
    try:
        start_index = aaseq.index("(")
    except ValueError:
        return aaseq
    stop_index = aaseq.index(")")
    aaseq = aaseq[:start_index] + aaseq[stop_index + 1:]
    return remove_modification_in_parentheses(aaseq)

def clean_an(an):
    return an.replace("CON__", "").split("-")[0]

def get_start_position_of_sequence_proteinGroups(df, fa):
    proteinGroup = df[COLUMN_PROTEINS]
    peptide = df["pepseq"]

    positions_list = []
    for an in proteinGroup.split(";"):
        try:
            aaseq = fa.get_aaseq_from_an(an)
        except KeyError:
            try:
                aaseq = fa.get_aaseq_from_an(clean_an(an))
            except KeyError:
                aaseq = False
        if aaseq:
            positions_list.append(get_start_position_from_aaseq(aaseq, peptide))
        else:
            positions_list.append(np.nan)
    return ";".join([str(ele) for ele in positions_list])

def get_start_position_from_aaseq(aaseq, peptide):
    try:
        start = aaseq.index(peptide)
    except ValueError:
        return np.nan
    return start

def start_counting_from_num(num_string, conventional_counting):
    try:
        num_string_split = num_string.split(";")
    except ValueError:
        return ""
    numbers_2_join = []
    for num in num_string_split:
        try:
            numbers_2_join.append(str(int(float(num)) + conventional_counting))
        except ValueError:
            numbers_2_join.append("nan")
    return ";".join(numbers_2_join)

def parse_sites_within_pepseq(aaseq, sites):
    try:
        start_index = aaseq.index(MODTYPE)
        start_index -= 1
    except ValueError:
        return sites
    stop_index = start_index + len(MODTYPE)
    sites_new = [start_index]
    sites += sites_new
    aaseq = aaseq[:start_index + 1] + aaseq[stop_index + 1:]
    return parse_sites_within_pepseq(aaseq, sites)

def parse_probabilities_grep_pos_2_prob(aaseq, pos_2_prob_dict):
    """
    the following didn't work as expected since list is not initialized as expected
    df[COLUMN_PROB] = df[COLUMN_MODPROB].apply(parse_probabilities, args=([], ))
    """
    try:
        start_index = aaseq.index("(")
    except ValueError:
        return pos_2_prob_dict
    stop_index = aaseq.index(")")
    probabilities_new = float(aaseq[start_index + 1 : stop_index])
    aaseq_position = start_index - 1
    pos_2_prob_dict[aaseq_position] = probabilities_new
    aaseq = aaseq[:start_index] + aaseq[stop_index + 1:]
    return parse_probabilities_grep_pos_2_prob(aaseq, pos_2_prob_dict)


if __name__ == "__main__":
    debug = True

    parser = argparse.ArgumentParser()

    parser.add_argument("-fa", "--fn_fasta", help="FASTA file absolute path", type=str, default=r"O:\Proteomics\USERS\BTW\FASTA\MOUSE20150706.fasta")

    parser.add_argument("-ev", "--fn_evidence", help="MaxQuant evidence file absolute path", type=str, default=r"O:\Proteomics\USERS\BTW\FASTA\evidence.txt")

    parser.add_argument("-o", "--fn_output", help="Output file name absolute path (defaults to 'EvidenceFileName_sites.txt')", type=str, default=None)

    parser.add_argument("-modtype", "--modification_type", help="name of the amino acid modification (default='(ac)'", type=str, default=None)

    parser.add_argument("-thr", "--probability_threshold", help="Filter threshold based on localisation probability (default=1.0 nothing removed), e.g. 0.9 for 90 percent", type=float, default=1.0)

    parser.add_argument("-cc", "--conventional_counting", help="Start counting from 0 or from 1 (default=1, meaning start from 1)", type=int, default=1)

    args = parser.parse_args()
    fn_fasta = args.fn_fasta
    fn_evidence = args.fn_evidence
    fn_output = args.fn_output
    probability_threshold = args.probability_threshold
    conventional_counting = args.conventional_counting

    if debug:
        fn_fasta = r"/Volumes/Speedy/FASTA/HUMAN20150706.fasta"
        fn_evidence = r"/Users/dblyon/CloudStation/CPR/BTW_sites/sites_positionsWithinProteins_input_v2.txt"
        fn_output = r"/Users/dblyon/CloudStation/CPR/BTW_sites/sites_positionsWithinProteins_output.txt"
        conventional_counting = 1

    if args.fn_output is None:
        fn_output = fn_evidence.replace(".txt", "_PosWithinProt.txt")
        args.fn_output = fn_output

    if args.modification_type is not None:
        MODTYPE = args.modification_type
    else:
        args.modification_type = MODTYPE

    print("#" * 80)
    for arg in sorted(vars(args)):
        print(arg, ": ", getattr(args, arg))

    run_sites(fn_fasta, fn_evidence, fn_output, conventional_counting)
    print("#" * 80)
    import this
    print("#" * 80)