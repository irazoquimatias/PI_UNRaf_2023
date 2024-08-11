from typing import Dict, List, Tuple
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
import json

path_ref_fasta: str = r'referencia.fasta'
path_ref_gbk: str = r'referencia.gbk'
path_muestras_fasta: str = r'muestras.fasta'

FORMAT_MAPPING: Dict[str, str] = {'gbk': 'genbank',
                                  'fasta': 'fasta'}

def get_codon_mutation_type(codon1: str, codon2: str) -> bool:
    """
    Determine if two codons are synonymous or nonsynonymous mutations.

    Args:
        codon1 (str): First codon sequence (3-letter string).
        codon2 (str): Second codon sequence (3-letter string).

    Returns:
        tuple: A tuple containing a string indicating the type of mutation ('synonymous' or 'nonsynonymous'),
               and the amino acids corresponding to each codon.
               Example: ('synonymous', 'A', 'A') if the codons are synonymous and code for the same amino acid,
                        or ('nonsynonymous', 'A', 'T') if they code for different amino acids.
    """
    seq1 = Seq.Seq(codon1)
    seq2 = Seq.Seq(codon2)

    aa1 = seq1.translate()
    aa2 = seq2.translate()
    
    if aa1 == aa2:
        return ('synonymous', aa1, aa2)
    else:
        return ('nonsynonymous', aa1, aa2)


def annotate_codon_mutations(seq1: str, seq2: str) -> Dict[str, List[str]]:
    """
    Compare two DNA sequences codon by codon and classify differences as synonymous or nonsynonymous.

    Args:
        seq1 (str): First DNA sequence.
        seq2 (str): Second DNA sequence.

    Returns:
        Dict[str, List[str]]: A dictionary with keys 'synonymous' and 'nonsynonymous', each containing
                              a list of strings describing the differences found.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")

    differences = {"synonymous": [], "nonsynonymous": []}

    for i in range(0, len(seq1), 3):
        codon1 = seq1[i:i+3]
        codon2 = seq2[i:i+3]

        if codon1 != codon2:
            mutation_type, amino1, amino2 = get_codon_mutation_type(codon1, codon2)
            differences[mutation_type].append(f"{amino1}{i+1}{amino2}")
            
    return differences

def extract_gene_locations(seq_record: SeqRecord) -> List[Tuple[str, FeatureLocation]]:
    """
    Extract gene names and their locations from a SeqRecord object.

    Args:
        seq_record (SeqRecord): A BioPython SeqRecord object containing features.

    Returns:
        List[Tuple[str, FeatureLocation]]: A list of tuples, where each tuple contains
        the gene name as the first element and the corresponding FeatureLocation as the second.
    """
    gene_location_pairs: List[Tuple[str, FeatureLocation]] = []
    for feature in seq_record.features:
        if feature.type == 'gene':
            gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
            gene_location_pairs.append((gene_name, feature.location))

    return gene_location_pairs

def get_format_for_extension(file_path: str) -> str:
    """
    Determine the format type from the file extension.

    Args:
        file_path (str): Path to the file.

    Returns:
        str: Format type ('genbank' or 'fasta').
    """
    extension: str = file_path.split('.')[-1]
    return FORMAT_MAPPING[extension]

def get_genes_locations(sequence: SeqRecord) -> List[FeatureLocation]:
    """
    Extract locations of gene features from a sequence record.

    Args:
        sequence (SeqRecord): BioPython SeqRecord object.

    Returns:
        List[FeatureLocation]: List of FeatureLocation objects representing gene locations.
    """
    gene_locations: List[FeatureLocation] = []
    for feature in sequence.features:
        if feature.type == 'gene':
            gene_locations.append(feature.location)
    return gene_locations

def compare_genes_between_reference_and_samples(reference_file_path: str, sample_file_path: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Compare gene sequences between a reference genome and multiple sample genomes to identify synonymous and nonsynonymous mutations.

    Args:
        reference_file_path (str): Path to the reference genome file in GenBank or FASTA format.
        sample_file_path (str): Path to the file containing sample genomes in GenBank or FASTA format.

    Returns:
        Dict[str, Dict[str, List[str]]]: A dictionary where:
            - Keys are sample IDs from the sample genome files.
            - Values are dictionaries where:
                - Keys are gene names from the reference genome.
                - Values are dictionaries with two keys: 'synonymous' and 'nonsynonymous'.
                    - 'synonymous' contains a list of descriptions of synonymous mutations.
                    - 'nonsynonymous' contains a list of descriptions of nonsynonymous mutations.

    The descriptions of mutations are in the format of "A1N1A2", where:
        - A1 is the amino acid from the reference sequence.
        - N1 is the codon position in the gene.
        - A2 is the amino acid from the sample sequence.

    Example:
        {
            "sample1": {
                "geneA": {
                    "synonymous": ["A15A"],
                    "nonsynonymous": ["A45T"]
                },
                "geneB": {
                    "synonymous": [],
                    "nonsynonymous": ["M30V"]
                }
            },
            "sample2": {
                "geneA": {
                    "synonymous": ["A15A"],
                    "nonsynonymous": ["A45G"]
                }
            }
        }
    """
    mutation_summary_by_sample = {}
    
    # Load the reference genome
    reference_seq_record: SeqRecord = SeqIO.read(reference_file_path, get_format_for_extension(reference_file_path))
    reference_genes_with_locations: List[Tuple[str, FeatureLocation]] = extract_gene_locations(reference_seq_record)

    # Load the sample genomes
    with open(sample_file_path) as sample_file:
        sample_seq_records = list(SeqIO.parse(sample_file, get_format_for_extension(sample_file_path)))

    # Iterate over each gene in the reference genome
    for gene_name, reference_feature_location in reference_genes_with_locations:
        reference_gene_sequence = reference_seq_record.seq[reference_feature_location.start:reference_feature_location.end]

        # Compare the gene sequences for each sample
        for sample_seq_record in sample_seq_records:
            sample_id = sample_seq_record.id
            if sample_id not in mutation_summary_by_sample:
                mutation_summary_by_sample[sample_id] = {}

            sample_gene_sequence = sample_seq_record.seq[reference_feature_location.start:reference_feature_location.end]

            if reference_gene_sequence != sample_gene_sequence:
                if gene_name not in mutation_summary_by_sample:
                    mutation_summary_by_sample[sample_id][gene_name] = {'synonymous': [], 'nonsynonymous': []}
    
                # Classify differences in codons between reference and sample
                codon_diff = annotate_codon_mutations(str(reference_gene_sequence), str(sample_gene_sequence))
                mutation_summary_by_sample[sample_id][gene_name]['synonymous'].extend(codon_diff['synonymous'])
                mutation_summary_by_sample[sample_id][gene_name]['nonsynonymous'].extend(codon_diff['nonsynonymous'])

    return mutation_summary_by_sample

print(
    json.dumps(
        compare_genes_between_reference_and_samples(path_ref_gbk, path_muestras_fasta),
        indent=4
    )
)
