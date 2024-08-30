# This module provides tools for analyzing genetic mutations by comparing reference and sample DNA sequences.
# It includes classes for representing codons, genes, and mutations, as well as a service for annotating mutations across multiple sample sequences.
from typing import Dict, List, Optional, Tuple
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from Bio import SeqIO, Align
import re
import json
import sys
import argparse
import os

__author__ = "Hernán Bernal"


class Codon:
    """
    Represents a codon, a sequence of three nucleotides in a DNA or RNA sequence.

    Attributes:
        sequence (str): The nucleotide sequence of the codon.
        position (int): The position of the codon in the sequence.
    """

    def __init__(self, sequence: str, position: int):
        """
        Initializes a Codon object with a given codon sequence and the codon position in the parent sequence.

        Args:
            sequence (str): The nucleotide sequence of the codon.
            position (int): The starting position of the codon in the parent sequence.
        """
        self.sequence = sequence
        self.position = position
        self.amino_acid = str(Seq(sequence).translate())

    def compare(self, other: "Codon") -> Tuple[str]:
        """
        Compares this codon with another codon to identify nucleotide differences.

        The position of the mutations is always in reference to the parent sequence.

        For example, 'a23g' indicates a nucleotide change from adenine to guanine at position 23.

        If the comparison involves a codon containing a gap ('-'), the method will still
        return the differences.

        Args:
            other (Codon): The other codon to compare with.

        Returns:
            Tuple[str]: A tuple containing strings that represent the differences between
                        the two codons.
                        (e.g., 'a12g' for a nucleotide change from 'a' to 'g' at position 12).

                        If there are no differences, an empty tuple is returned.

        Raises:
            ValueError: If the codons are not of the same length.
        """
        if len(self.sequence) != len(other.sequence):
            raise ValueError("Codons must be of the same length.")

        differences = [
            f"{ref.lower()}{self.position + i + 1}{sample.lower()}"
            for i, (ref, sample) in enumerate(zip(self.sequence, other.sequence))
            if ref != sample
        ]

        return tuple(differences)

    def translate(self) -> str:
        """
        Translates the codon sequence into its corresponding amino acid.

        Returns:
            str: The single-letter code of the amino acid translated from the codon.
        """
        return self.amino_acid


class Gene:
    """
    Represents a gene in a DNA sequence.

    Attributes:
        name (str): The name of the gene.
        sequence (Seq): The nucleotide sequence of the gene.
    """

    def __init__(self, name: str, sequence: Seq, location: FeatureLocation):
        """
        Initializes a Gene object with a given name and sequence.

        Args:
            name (str): The name of the gene.
            sequence (Seq): The nucleotide sequence of the gene.
        """
        self.name = name
        self.sequence = sequence
        self.location = location


class MutationService:
    """
    Service for analyzing mutations between a reference sequence and sample sequences.

    Attributes:
        ref_seq (Seq): The nucleotide sequence of the reference genome.
        genes (List[Gene]): A list of Gene objects representing the genes in the reference genome.
    """

    def __init__(self, seq_reference: SeqRecord, print_progress: bool = False):
        """
        Initializes a MutationService object with a reference sequence.

        Args:
            seq_reference (SeqRecord): A SeqRecord object containing the reference sequence and its features.
        """
        self.ref_seq = seq_reference.seq
        self.genes = self.__populate_genes(seq_reference)
        self._print_progress = print_progress
        self._aligner = Align.PairwiseAligner()
        self._aligner.mismatch_score = -0.2
        self._aligner.open_gap_score = -2

    def __populate_genes(self, sequence: SeqRecord) -> List[Gene]:
        """
        Populates a list of Gene objects from the given SeqRecord.

        Args:
            sequence (SeqRecord): A SeqRecord object containing the nucleotide sequence and features.

        Returns:
            List[Gene]: A list of Gene objects extracted from the SeqRecord.
        """
        genes: List[Gene] = []
        for feature in sequence.features:
            if feature.type == "gene":
                gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
                gene_seq = sequence.seq[
                    feature.location.start : feature.location.end
                ]
                genes.append(Gene(gene_name, gene_seq, feature.location))
        return genes

    def __annotate_gene_mutations(
        self, seq_sample: Seq, print_progress: bool = False
    ) -> Dict[str, List[str]]:
        """
        Annotates mutations between the reference sequence and a sample sequence for each gene.

        Args:
            seq_sample (Seq): The nucleotide sequence of the sample.
            print_progress (bool, optional): Whether to print progress updates. Defaults to False.

        Returns:
            Dict[str, List[str]]: A dictionary where keys are gene names and values are lists of mutation annotations
                                  categorized as 'synonymous', 'nonsynonymous', 'deletions' and 'insertions'.

        Raises:
            ValueError: If the length of the reference sequence and the sample sequence do not match.
        """
        # Length of reference sequence and sample sequence must be equal.
        if len(self.ref_seq) != len(seq_sample):
            raise ValueError(
                f"Sequences must be of equal length.\nReference sequence len: {len(self.ref_seq)}\nSample sequence len: {len(seq_sample)}"
            )

        annotations = {}
        for gene in self.genes:
            if gene.name not in annotations:
                annotations[gene.name] = {
                    "synonymous": [],
                    "nonsynonymous": [],
                    "insertions": [],
                    "deletions": [],
                }

            # Print progress
            if self._print_progress:
                print(
                    f"Processing gene: {gene.name}, Numbers of codons: {len(gene.sequence) / 3}"
                )  # Print gene name and total codons for this gene.

            num_codons = len(gene.sequence) // 3

            for pos in range(0, len(gene.sequence), 3): # The position is relative to the gene
                # Print progress
                if self._print_progress:
                    progress = (pos // 3) + 1
                    percent_complete = (progress / num_codons) * 100
                    print(
                        f"Processing codon {progress}/{num_codons} ({percent_complete:.2f}%)"
                    )

                # Getting sample gene sequence
                sample_gene_sequence = seq_sample[
                    gene.location.start : gene.location.end
                ]

                # Getting codons for mutation calculation.
                ref_codon = Codon(gene.sequence[pos : pos + 3], pos)
                ref_codon_next = Codon(gene.sequence[pos + 3 : pos + 6], pos + 3)

                sample_codon = Codon(sample_gene_sequence[pos : pos + 3], pos)
                sample_codon_next = Codon(sample_gene_sequence[pos + 3 : pos + 6], pos + 3)

                # Start calculating mutations.
                ## Insertion case
                if ref_codon_next.sequence == '---':
                    annotated_insertion = f"{ref_codon.amino_acid}{pos}{sample_codon.amino_acid}{sample_codon_next.amino_acid}"
                    annotations[gene.name]['insertions'].append(annotated_insertion)
                    pos = pos + 3

                ## Deletion case
                elif sample_codon.sequence == '---':
                    annotated_deletion = f"{pos}{ref_codon.amino_acid}"
                    annotations[gene.name]['deletions'].append(annotated_deletion)

                ## Synonymous or Nonsynonymous case
                elif ref_codon.sequence != sample_codon.sequence:
                    ## Synonymous case
                    if ref_codon.amino_acid == sample_codon.amino_acid:
                        nucleotid_differ = ref_codon.compare(sample_codon)
                        annotations[gene.name]['synonymous'].extend(nucleotid_differ)
                    ## Nonsynonymous case
                    else:
                        aminoacid_differ = (
                            f"{ref_codon.amino_acid}{ref_codon.position}{sample_codon.amino_acid}",
                        )
                        annotations[gene.name]['nonsynonymous'].extend(aminoacid_differ)

        return annotations

    def __align_sequences(self, seq_reference: Seq, seq_sample: Seq) -> Seq:
        """
        Aligns a sample sequence to the reference sequence using pairwise alignment.

        This method uses the `Align.PairwiseAligner` to align the provided sample sequence to the reference sequence. 
        It adjusts the reference sequence to match the alignment and returns the aligned sample sequence.

        Args:
            seq_reference (Seq): The nucleotide sequence of the reference genome.
            seq_sample (Seq): The nucleotide sequence of the sample to be aligned.

        Returns:
            Seq: The aligned nucleotide sequence of the sample.

        Raises:
            ValueError: If the reference and sample sequences are of different lengths.
        """
        # Perform the alignment
        alignments = self._aligner.align(seq_reference, seq_sample)
        
        # Extract the first (best) alignment
        aligned_reference = alignments[0][0]
        aligned_sample = alignments[0][1]
        
        # Update the reference sequence if it has been modified during alignment
        if self.ref_seq != aligned_reference:
            self.ref_seq = aligned_reference

        return aligned_sample

    def analyze_mutations(
        self, seq_samples: List[SeqRecord]
    ) -> Dict[str, Dict[str, List[str]]]:
        """
        Analyzes mutations between the reference sequence and multiple sample sequences.

        Args:
            seq_samples (List[SeqRecord]): A list of SeqRecord objects, each containing a sample sequence and its ID.

        Returns:
            Dict[str, Dict[str, List[str]]]: A dictionary where keys are sample IDs and values are dictionaries with
                                             gene names as keys and lists of mutation annotations as values.
        """
        annotations = {}
        for seq_record_sample in seq_samples:
            if seq_record_sample.id not in annotations:
                annotations[seq_record_sample.id] = {}

            seq_sample = seq_record_sample.seq
            aligned_seq_sample = self.__align_sequences(self.ref_seq, seq_sample)
            # Print progress
            if self._print_progress:
                print(f"Processing sample: {seq_record_sample.id}")
            gene_mutation_details = self.__annotate_gene_mutations(aligned_seq_sample)
            annotations[seq_record_sample.id] = gene_mutation_details

        return annotations

FORMAT_MAPPING: Dict[str, str] = {'gbk': 'genbank',
                                  'fasta': 'fasta'}

def get_format_for_extension(file_path: str) -> str:
    """
    Determine the format type from the file extension.

    Args:
        file_path (str): Path to the file.

    Returns:
        str: Format type ('genbank' or 'fasta').
    """
    extension: str = file_path.rsplit('.', 1)[-1].lower()
    if extension in FORMAT_MAPPING:
        return FORMAT_MAPPING[extension]
    else:
        raise ValueError(f"Unsupported file format: {extension}")

def parse_options():
    parser = argparse.ArgumentParser(description='Genetic Mutation Analyzer')
    parser.add_argument('--reference-file-path', required=True, help='Path to the reference genome file.')
    parser.add_argument('--samples-file-path', required=True, help='Path to the sample sequences file.')
    parser.add_argument('--print-progress', action='store_true', help='Print progress updates during analysis.')
    return parser.parse_args()

# Main code
if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_options()

    # Attempt to read the reference sequence file
    try:
        # Load the reference sequence from the specified file path
        seq_record_reference = SeqIO.read(
            args.reference_file_path,
            get_format_for_extension(args.reference_file_path)
        )
    except FileNotFoundError:
        # Handle the case where the reference file is not found
        print(f"ERROR: Reference file '{args.reference_file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        # Handle other potential errors when reading the reference file
        print(f"ERROR: Failed to read reference file '{args.reference_file_path}': {e}")
        sys.exit(1)

    # Attempt to read the sample sequences file
    try:
        # Open and parse the sample sequences file
        with open(args.samples_file_path) as sample_file:
            sample_seq_records = list(
                SeqIO.parse(sample_file, get_format_for_extension(args.samples_file_path))
            )
    except FileNotFoundError:
        # Handle the case where the sample file is not found
        print(f"ERROR: Samples file '{args.samples_file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        # Handle other potential errors when reading the sample file
        print(f"ERROR: Failed to read samples file '{args.samples_file_path}': {e}")
        sys.exit(1)

    samples_with_gaps = []
    for seq in sample_seq_records:
        seq.seq = seq.seq.replace('-', '')
        samples_with_gaps.append(seq)

    # Create a MutationService instance with the loaded reference sequence
    mutation_service = MutationService(seq_record_reference, args.print_progress)

    # Analyze mutations across all sample sequences
    mutation_report = mutation_service.analyze_mutations(samples_with_gaps)

    # Print the resulting mutation report in a formatted JSON
    print(json.dumps(mutation_report, indent=4))
