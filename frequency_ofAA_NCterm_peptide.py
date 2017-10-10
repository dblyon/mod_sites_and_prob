from __future__ import print_function
import re
from collections import deque
import itertools
from tqdm import tqdm


def run(fn_fasta, fn_out, aa_2_count="K", enzyme="trypin", missed_cleavages=2, min_length=7, max_length=25):
    print("#" * 80)
    print("\tfantastic script is doing stuff\n".upper())
    print("\tplease stay where you are, the successful completion is not dependendant on your location, but as soon as you finished reading this it will most probably be finished anyway...\n")
    fa = Fasta()
    fa.set_file(fn_fasta)
    fa.parse_fasta()
    peptide_list = []
    for an, aaseq in tqdm(fa.an2aaseq_dict.items()):
        peptide_list += cleave(aaseq, expasy_rules[enzyme], missed_cleavages=missed_cleavages, min_length=min_length)

    if max_length is not None:
        peptide_list = [ele for ele in peptide_list if len(ele) <= max_length]

    N_0, N_1, N_2, N_3 = 0, 0, 0, 0
    C_0, C_1, C_2, C_3 = 0, 0, 0, 0
    for peptide in peptide_list:
        peptide_as_list = list(peptide)
        if peptide_as_list[0] == aa_2_count:
            N_0 += 1
        if peptide_as_list[1] == aa_2_count:
            N_1 += 1
        if peptide_as_list[2] == aa_2_count:
            N_2 += 1
        if peptide_as_list[3] == aa_2_count:
            N_3 += 1

        if peptide_as_list[-1] == aa_2_count:
            C_0 += 1
        if peptide_as_list[-2] == aa_2_count:
            C_1 += 1
        if peptide_as_list[-3] == aa_2_count:
            C_2 += 1
        if peptide_as_list[-4] == aa_2_count:
            C_3 += 1

    with open(fn_out, "w") as fh_out:
        fh_out.write("Position\tCount\n")
        fh_out.write("N1\t{}\n".format(str(N_0)))
        fh_out.write("N2\t{}\n".format(str(N_1)))
        fh_out.write("N3\t{}\n".format(str(N_2)))
        fh_out.write("N4\t{}\n".format(str(N_3)))
        fh_out.write("C4\t{}\n".format(str(C_3)))
        fh_out.write("C3\t{}\n".format(str(C_2)))
        fh_out.write("C2\t{}\n".format(str(C_1)))
        fh_out.write("C1\t{}\n".format(str(C_0)))

    print("\n...finished processing.\nClose the terminal and spend some time with your kids.\n")
    print("or entertain yourself: {}\n you deserve it.".format("https://xkcd.com/"))
    print("#" * 80)


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
        organism_regexes = [
            r"^>(dbj|emb|na|sp|tr|gb|ref)\|(.*?)\|.*?OS=(.*?)\s(?:(?:\w+=\w+\s?)+)(?:NCBI_TaxID=(\d+))?",
            # r"^>([s|t][p|r])\|([\w]+).*?OS=(.*?) [A-Z][A-Z]=",
            r"^>([s|t][p|r])\|(.*)\|.*?OS=(.*?) (?:[A-Z][A-Z]=)?",  # get Isoforms e.g. "A0AUZ9-2", make OS blabla optional
            r"^>gi\|\d+\|(\w+)\|(\w+\.\d).*\[(.*)\]",
            r"^>([^ ]+) ([^\[]+)"]
        self.my_regex = re.compile('|'.join(organism_regexes))
        self.an2aaseq_dict = {}
        self.org2an_dict = {}
        self.an2org_dict = {}
        self.db2an_dict = {}
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
            if match: # at least one regex found something
                match_groups = filter(None, match.groups())
                if len(match_groups) == 2:
                    an, organism = match_groups
                    database = "na"
                elif len(match_groups) == 3:
                    database, an, organism = match_groups
                    database = database.strip()
                else:
                    print("Nothing could be extracted from the following fasta-entry: ", "\n", header, aaseq)
                    raise StopIteration
                an = an.strip()
                organism = organism.strip() # #!!! stop
            else:
                print("Nothing could be extracted from the following fasta-entry: ", "\n", header, aaseq)
                raise StopIteration
            # now fill dicts
            if not self.an2aaseq_dict.has_key(an):
                self.an2aaseq_dict[an] = aaseq
            else: # AccessionNumbers should be unique identifiers
                if unique:
                    print("AccessionNumbers should be unique identifiers!")
                    print(an, "\n", self.an2aaseq_dict[an], "\n", header)
                    raise StopIteration
                else:
                    pass
                    # print("AN double entry: {}".format(an))
            if not self.org2an_dict.has_key(organism):
                self.org2an_dict[organism] = [an]
            else:
                self.org2an_dict[organism].append(an)

            if not self.an2org_dict.has_key(an):
                self.an2org_dict[an] = organism
            else:
                if unique:
                    print("AccessionNumbers should be unique identifiers!")
                    print(an, "\n", self.an2org_dict[an], "\n", aaseq)
                    raise StopIteration
                else:
                    pass
            if not self.db2an_dict.has_key(database):
                self.db2an_dict[database] = [an]
            else:
                self.db2an_dict[database].append(an)

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
        fh = open(self.get_file(), "r")
        aaseq = ""
        header = ""
        did_first = False
        for line in fh:
            if line[0] == ">":
                if did_first:
                    if len(aaseq) > 0:
                        yield(header, aaseq)
                else:
                    did_first = True
                header = line.rstrip()
                aaseq = ""
            else:
                aaseq += line.strip().upper()
        if len(aaseq) > 0:
            yield(header, aaseq)
        fh.close()

##### adapted from Pyteomics
##### https://bitbucket.org/levitsky/pyteomics
std_amino_acids = ['Q', 'W', 'E', 'R', 'T', 'Y', 'I', 'P', 'A', 'S',
                   'D', 'F', 'G', 'H', 'K', 'L', 'C', 'V', 'N', 'M']

expasy_rules = {
    'arg-c': r'R',
    'asp-n': r'\w(?=D)',
    'bnps-skatole': r'W',
    'caspase 1': r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2': r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3': r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4': r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5': r'(?<=[LW]EH)D',
    'caspase 6': r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7': r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8': r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9': r'(?<=LEH)D',
    'caspase 10': r'(?<=IEA)D',
    'chymotrypsin high specificity': r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain': r'R',
    'cnbr': r'M',
    'enterokinase': r'(?<=[DE]{3})K',
    'factor xa': r'(?<=[AFGILTVM][DE]G)R',
    'formic acid': r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b': r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc': r'K',
    'ntcb': r'\w(?=C)',
    'pepsin ph1.3': r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                    r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'pepsin ph2.0': r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                    r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k': r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin': r'[^DE](?=[AFILMV])',
    'thrombin': r'((?<=G)R(?=G))|'
                r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin': r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    'trypsin_no_exceptions': r'([KR])'
}
"""
This dict contains regular expressions for cleavage rules of the most
popular proteolytic enzymes. The rules were taken from the
`PeptideCutter tool
<http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_
at Expasy.
"""

"""modX labels for the 20 standard amino acids."""

std_nterm = 'H-'
"""modX label for the unmodified N-terminus."""

std_cterm = '-OH'
"""modX label for the unmodified C-terminus."""

std_labels = std_amino_acids + [std_nterm, std_cterm]
"""modX labels for the standard amino acids and unmodified termini."""

def is_term_mod(label):
    """Check if `label` corresponds to a terminal modification.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool
    """
    return label[0] == '-' or label[-1] == '-'

def length(sequence, **kwargs):
    """Calculate the number of amino acid residues in a polypeptide
    written in modX notation.

    Parameters
    ----------
    sequence : str or list or dict
        A string with a polypeptide sequence, a list with a parsed sequence or
        a dict of amino acid composition.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Examples
    --------
    >>> length('PEPTIDE')
    7
    >>> length('H-PEPTIDE-OH')
    7
    """
    if not sequence: return 0

    if isinstance(sequence, str) or isinstance(sequence, list):
        if isinstance(sequence, str):
            parsed_sequence = parse(sequence, **kwargs)
        else:
            parsed_sequence = sequence
        num_term_groups = 0
        if is_term_mod(parsed_sequence[0]):
            num_term_groups += 1
        if is_term_mod(parsed_sequence[-1]):
            num_term_groups += 1
        return len(parsed_sequence) - num_term_groups
    elif isinstance(sequence, dict):
        return sum(amount for aa, amount in sequence.items()
                   if not is_term_mod(aa))

    print('Unsupported type of sequence.')


_modX_sequence = re.compile(r'^([^-]+-)?((?:[a-z]*[A-Z])+)(-[^-]+)?$')
_modX_group = re.compile(r'[a-z]*[A-Z]')
_modX_split = re.compile(r'([a-z]*)([A-Z])')

def parse(sequence, show_unmodified_termini=False, split=False,
          allow_unknown_modifications=False,
          **kwargs):
    """Parse a sequence string written in modX notation into a list of
    labels or (if `split` argument is :py:const:`True`) into a list of
    tuples representing amino acid residues and their modifications.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    show_unmodified_termini : bool, optional
        If :py:const:`True` then the unmodified N- and C-termini are explicitly
        shown in the returned list. Default value is :py:const:`False`.
    split : bool, optional
        If :py:const:`True` then the result will be a list of tuples with 1 to 4
        elements: terminal modification, modification, residue. Default value is
        :py:const:`False`.
    allow_unknown_modifications : bool, optional
        If :py:const:`True` then do not raise an exception when an unknown
        modification of a known amino acid residue is found in the sequence.
        This also includes terminal groups.
        Default value is :py:const:`False`.

        .. note::
            Since version 2.5, this parameter has effect only if `labels`
            are provided.
    labels : container, optional
        A container of allowed labels for amino acids,
        modifications and terminal modifications.
        If not provided, no checks will be done.
        Separate labels for modifications (such as 'p' or 'ox')
        can be supplied, which means they are applicable to all residues.

        .. warning::
            If `show_unmodified_termini` is set to :py:const:`True`, standard
            terminal groups need to be present in `labels`.

        .. warning::
            Avoid using sequences with only one terminal group, as they are
            ambiguous. If you provide one, `labels` (or :py:const:`std_labels`)
            will be used to resolve the ambiguity.

    Returns
    -------
    out : list
        List of tuples with labels of modifications and amino acid residues.

    Examples
    --------
    >>> parse('PEPTIDE', split=True)
    [('P',), ('E',), ('P',), ('T',), ('I',), ('D',), ('E',)]
    >>> parse('H-PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse('PEPTIDE', show_unmodified_termini=True)
    ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH']
    >>> parse('TEpSToxM', labels=std_labels + ['pS', 'oxM'])
    ['T', 'E', 'pS', 'T', 'oxM']
    >>> parse('zPEPzTIDzE', True, True, labels=std_labels+['z'])
    [('H-', 'z', 'P'), ('E',), ('P',), ('z', 'T'), ('I',), ('D',), ('z', 'E', '-OH')]
    """
    sequence = str(sequence)

    try:
        n, body, c = re.match(_modX_sequence, sequence).groups()
    except AttributeError:
        print('Not a valid modX sequence: ' + sequence)

    # Check for allowed labels, if they were explicitly given
    labels = kwargs.get('labels')
    # labels help save the day when only one terminal group is given
    if c is None and n is not None:
        if labels is None:
            labels = std_labels
        # we can try to resolve the ambiguity
        if n != std_nterm and n not in labels:
            # n is the body then
            c = '-' + body
            body = n[:-1]
            n = None

    # Actual parsing
    if split:
        parsed_sequence = [g if g[0] else (g[1],) for g in re.findall(
            _modX_split, body)]
    else:
        parsed_sequence = re.findall(_modX_group, body)
    nterm, cterm = (n or std_nterm), (c or std_cterm)

    # Check against `labels` if given
    if labels is not None:
        labels = set(labels)
        for term, std_term in zip([n, c], [std_nterm, std_cterm]):
            if term and term not in labels and not allow_unknown_modifications:
                print(
                    'Unknown label: {}'.format(term))
        for group in parsed_sequence:
            if split:
                mod, X = group if len(group) == 2 else ('', group[0])
            else:
                mod, X = re.match(_modX_split, group).groups()
            if ((not mod) and X not in labels) or not ((mod + X in labels) or (
                            X in labels and (
                                    mod in labels or allow_unknown_modifications))):
                print(
                    'Unknown label: {}'.format(group))

    # Append terminal labels
    if show_unmodified_termini or nterm != std_nterm:
        if split:
            parsed_sequence[0] = (nterm,) + parsed_sequence[0]
        else:
            parsed_sequence.insert(0, nterm)
    if show_unmodified_termini or cterm != std_cterm:
        if split:
            parsed_sequence[-1] = parsed_sequence[-1] + (cterm,)
        else:
            parsed_sequence.append(cterm)

    return parsed_sequence

def cleave(sequence, rule, missed_cleavages=0, min_length=None, **kwargs):
    """Cleaves a polypeptide sequence using a given rule.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    rule : str or compiled regex
        A regular expression describing the site of cleavage. It is recommended
        to design the regex so that it matches only the residue whose C-terminal
        bond is to be cleaved. All additional requirements should be specified
        using `lookaround assertions
        <http://www.regular-expressions.info/lookaround.html>`_.
    missed_cleavages : int, optional
        Maximum number of allowed missed cleavages. Defaults to 0.
    min_length : int or None, optional
        Minimum peptide length. Defaults to :py:const:`None`.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Returns
    -------
    out : set
        A set of unique (!) peptides.

    Examples
    --------
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0, labels='ABK') == {'AK', 'BK'}
    True
    >>> cleave('GKGKYKCK', expasy_rules['trypsin'], 2) == \
    {'CK', 'GKYK', 'YKCK', 'GKGK', 'GKYKCK', 'GK', 'GKGKYK', 'YK'}
    True

    """
    sequence = sequence.replace("*", "")
    return sorted(set(_cleave(sequence, rule, missed_cleavages, min_length, **kwargs)))

def _cleave(sequence, rule, missed_cleavages=2, min_length=7, **kwargs):
    """Like :py:func:`cleave`, but the result is a list. Refer to
    :py:func:`cleave` for explanation of parameters.
    """
    peptides = []
    cleavage_sites = deque([0], maxlen=missed_cleavages + 2)
    for i in itertools.chain(map(lambda x: x.end(), re.finditer(rule, sequence)), [None]):
        cleavage_sites.append(i)
        for j in range(len(cleavage_sites) - 1):
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or length(seq, **kwargs) >= min_length:
                    peptides.append(seq)
    return peptides


if __name__ == "__main__":
    missed_cleavages = 2
    min_length = 7
    max_length = 25 # None
    enzyme = "trypsin" # "trypsin_no_exceptions"
    # Please note: "trypsin" adheres to expasy rules (http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#Tryps)
    # "trypsin_no_exceptions" simply cuts at "K" or "R" without any exceptions.
    aa_2_count = "K"
    # fn_fasta = r"/Volumes/cpr-1/Proteomics/USERS/BTW/FASTA/YEAST_s288c_20150706.fasta"
    fn_fasta = r"/Volumes/Speedy/FASTA/Escherichia_coli_RefProt_20170505.fasta"
    # fn_out = r"/Volumes/cpr-1/Proteomics/USERS/BTW/FASTA/YEAST_s288c_20150706_internal_Lysines.txt"
    fn_out = r"/Users/dblyon/Downloads/test_output.txt"
    run(fn_fasta, fn_out, aa_2_count, enzyme, missed_cleavages, min_length, max_length)
