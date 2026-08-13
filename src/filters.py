"""Filters that narrow a cancer gene list down to cancer-selective candidates."""
import re
import pandas as pd


def cancer_specific_genes(gene_map_cancer, gene_map_normal):
    """Genes that were hits in cancer screens but never showed up as a hit
    in any normal-cell screen."""
    normal_ids = set(gene_map_normal.keys())
    return [g for gid, g in gene_map_cancer.items() if gid not in normal_ids]

# why did we have to get the crisper essential list? Becuase we only had ~95 normal screens 
# and through that you cant extract all the essential genes
def load_depmap_common_essentials(csv_path):
    """DepMap's list of genes essential in nearly all cell lines (pan-essential
    housekeeping genes: ribosome, splicing, DNA replication, etc.)."""
    df = pd.read_csv(csv_path)
    entrez_ids = df.iloc[:, 0].str.extract(r"\((\d+)\)")[0]
    return set(entrez_ids.dropna().astype(int))


def exclude_ids(genes, exclude_id_set):
    return [g for g in genes if g.identifier_id not in exclude_id_set]


# Symbol patterns for common housekeeping gene families DepMap's list
# doesn't always fully cover (used as a supplementary cleanup pass).
HOUSEKEEPING_PATTERN = re.compile(
    r"^("
    r"RPL\d|RPS\d|MRPL\d|MRPS\d"        # ribosomal proteins
    r"|EIF\d|EEF\d"                      # translation factors
    r"|POLR\d|POLE\d?|POLD\d?"           # RNA/DNA polymerase subunits
    r"|PSMA\d|PSMB\d|PSMC\d|PSMD\d"      # proteasome subunits
    r"|SNRNP\d|SNRPA|SNRPB|SNRPD|SF3[AB]" # spliceosome components
    r"|[A-Z]{1,3}ARS\d?$"                # aminoacyl-tRNA synthetases
    r"|COX\d|NDUF|ATP5"                  # mitochondrial respiratory chain
    r"|HIST\d|H2A|H2B|H3F|H4C"           # core histones
    r")"
)
KNOWN_PAN_ESSENTIAL_SYMBOLS = {"PRPF19", "PCNA", "KIF11", "DCAF1"}


def is_housekeeping_symbol(symbol):
    if not symbol:
        return False
    if symbol in KNOWN_PAN_ESSENTIAL_SYMBOLS:
        return True
    return bool(HOUSEKEEPING_PATTERN.match(symbol))


def exclude_housekeeping(genes):
    result_genes = []
    for gene in genes:
        if is_housekeeping_symbol(gene.symbol) != True:
            result_genes.append(gene)
    return result_genes
