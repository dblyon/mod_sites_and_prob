from __future__ import print_function
import re, sys, argparse
import pandas as pd
import numpy as np


COLUMN_MODSEQ = "Modified sequence"
COLUMN_LEADRAZPROT = "Leading razor protein"
COLUMN_MODPROB = "Acetyl (K) Probabilities"
MODTYPE = "(ac)"

##### new columns
COLUMN_SITES = "Sites"
COLUMN_PROB = "Probability"


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
                match_groups = filter(None, match.groups())  # #!!! start, since Fasta entries not consistent
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


def select_first_of_colon_sep(ans_sep_string):
    try:
        return ans_sep_string.split(";")[0]
    except AttributeError:
        return ans_sep_string

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

def get_start_position_of_sequence(df, fa):
    an = df[COLUMN_LEADRAZPROT]
    try:
        aaseq = fa.get_aaseq_from_an(an)
    except KeyError:
        try:
            aaseq = fa.get_aaseq_from_an(clean_an(an))
        except KeyError:
            return np.nan
    peptide = df["pepseq"]
    try:
        start = aaseq.index(peptide)
    except ValueError:
        return np.nan
    return start

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

def remove_modifications_not_MODTYPE(aaseq, my_regex):
    strings_2_replace = []
    aaseq = aaseq.replace("_", "")
    for match in my_regex.finditer(aaseq):
        modification = match.group()
        if modification != MODTYPE:
            strings_2_replace.append(modification)
    for group in strings_2_replace:
        aaseq = aaseq.replace(group, "")
    return aaseq

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

def add_COLUMN_SITES_and_PROB_2_df(df):
    COLUMN_SITES_list, COLUMN_PROB_list = [], []
    for index_, row in df.iterrows():
        start_pos = row["start_pos"]
        pepseq_mod = row["pepseq_mod"]
        pepseq_prob = row[COLUMN_MODPROB]

        sites_list = parse_sites_within_pepseq(pepseq_mod, [])
        COLUMN_SITES_list.append(";".join([str(ele + start_pos) for ele in sites_list]))

        pos_2_prob_dict = parse_probabilities_grep_pos_2_prob(pepseq_prob, {})
        COLUMN_PROB_list.append(";".join([str(pos_2_prob_dict[key]) for key in pos_2_prob_dict if key in sites_list]))
    df[COLUMN_SITES] = pd.Series(COLUMN_SITES_list)
    df[COLUMN_PROB] = pd.Series(COLUMN_PROB_list)
    return df

def is_any_above_threshold(numbers_as_string, probability_threshold):
    for num in numbers_as_string.split(";"):
        if float(num) >= probability_threshold:
            return True
    return False

def error_(parser):
    sys.stderr.write("The arguments passed are invalid.\nPlease check the input parameters.\n\n")
    parser.print_help()
    sys.exit(2)

def start_counting_from_num(num_string, conventional_counting):
    try:
        num_string_split = num_string.split(";")
    except ValueError:
        return ""
    return ";".join([str(float(num) + conventional_counting) for num in num_string_split])

def run_sites(fn_fasta, fn_evidence, fn_output, probability_threshold, conventional_counting):
    fa = Fasta()
    fa.set_file(fn_fasta)
    fa.parse_fasta()

    df = pd.read_csv(fn_evidence, sep='\t')
    cols_needed = [COLUMN_MODSEQ, COLUMN_MODPROB, COLUMN_LEADRAZPROT]
    df = df[cols_needed]
    df.dropna(axis=0, how="all", inplace=True)

    # choose first AN if mulitple ANs in COLUMN_LEADRAZPROT
    df[COLUMN_LEADRAZPROT] = df[COLUMN_LEADRAZPROT].apply(select_first_of_colon_sep)

    df["pepseq"] = df[COLUMN_MODSEQ].apply(lambda aaseq: aaseq.replace("_", "").replace(MODTYPE, ""))
    df["pepseq"] = df["pepseq"].apply(remove_modification_in_parentheses)
    df["start_pos"] = df.apply(get_start_position_of_sequence, args=(fa, ), axis=1)
    df = df[df[COLUMN_MODPROB].notnull()]

    # add sites and probabilities
    my_regex = re.compile(r"(\(\w+\))")
    df["pepseq_mod"] = df[COLUMN_MODSEQ].apply(remove_modifications_not_MODTYPE, args=(my_regex,))
    df = add_COLUMN_SITES_and_PROB_2_df(df)
    df = df[df[COLUMN_SITES].notnull()]

    if probability_threshold > 0:
        df = df[df[COLUMN_PROB].apply(is_any_above_threshold, args=(probability_threshold,))]

    if conventional_counting > 0:
        df[COLUMN_SITES] = df[COLUMN_SITES].apply(start_counting_from_num, args=(conventional_counting, )) #lambda num_string: ";".join([str(int(float(num))) + conventional_counting for num in num_string.split(";")])

    # keep only relevant columns and write to file
    df2write = df[[COLUMN_MODSEQ, COLUMN_MODPROB, COLUMN_LEADRAZPROT, COLUMN_SITES, COLUMN_PROB]]
    df2write.to_csv(fn_output, sep='\t', header=True, index=False)
    return df2write


if __name__ == "__main__":
    debug = False

    parser = argparse.ArgumentParser()

    parser.add_argument("-fa", "--fn_fasta", help="FASTA file absolute path", type=str, default=r"O:\Proteomics\USERS\BTW\FASTA\MOUSE20150706.fasta")

    parser.add_argument("-ev", "--fn_evidence", help="MaxQuant evidence file absolute path", type=str)

    parser.add_argument("-o", "--fn_output", help="Output file name absolute path (defaults to 'EvidenceFileName_sites.txt')", type=str, default=None)

    parser.add_argument("-modprob", "--modification_probability_column", help="Name of the modification probability column (default='Acetyl (K) Probabilities'", type=str, default=None)

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
        fn_fasta = r"/Users/dblyon/CloudStation/CPR/BTW_sites/MOUSE20150706.fasta"
        fn_evidence = r"/Users/dblyon/CloudStation/CPR/BTW_sites/evidence.txt"
        fn_output = r"/Users/dblyon/CloudStation/CPR/BTW_sites/evidence_output.txt"

    if args.fn_output is None:
        fn_output = fn_evidence.replace(".txt", "_sites.txt")
        args.fn_output = fn_output

    if args.modification_probability_column is not None:
        COLUMN_MODPROB = args.modification_probability_column
    else:
        args.modification_probability_column = COLUMN_MODPROB

    if args.modification_type is not None:
        MODTYPE = args.modification_type
    else:
        args.modification_type = MODTYPE

    print("#" * 80)
    print("\tawesome script is running like crazy\n".upper())
    for arg in sorted(vars(args)):
        print(arg, ": ", getattr(args, arg))

    run_sites(fn_fasta, fn_evidence, fn_output, probability_threshold, conventional_counting)
    print("\n...finished processing.\nClose this shell window and publish big.")
    print("#" * 80)
