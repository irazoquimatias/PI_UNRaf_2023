# This module provides tools for analyzing genetic mutations by comparing reference and sample DNA sequences.
# It includes classes for representing codons, genes, and mutations, as well as a service for annotating mutations across multiple sample sequences.
from typing import Dict, List, Optional, Tuple
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from Bio import SeqIO
import re
import json
import sys

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

    def compare(self, other: "Codon") -> Tuple[str]:
        """
        Compares this codon with another codon to identify nucleotide differences.

        The position of the mutations is always in reference to the parent sequence.

        Nucleotide changes are reported in lowercase letters, while amino acid changes
        are reported in uppercase letters.

        For example, 'a23g' indicates a nucleotide change from adenine to guanine at position 23,
        and 'A23G' indicates an amino acid change from alanine to glycine at the corresponding position.

        If the comparison involves a codon containing a gap ('-'), the method will still
        return the differences, but such cases are typically classified as 'unknown'
        because gaps can indicate complex events like insertions or deletions that do not
        directly correspond to a single nucleotide substitution.

        For example, 'a23-' indicates a nucleotide change from adenine to a gap
        at position 23, which could suggest a deletion. Such a case would be
        classified as 'unknown'.

        Args:
            other (Codon): The other codon to compare with.

        Returns:
            Tuple[str]: A tuple containing strings that represent the differences between
                        the two codons. Each string is formatted as 'nucleotide_position_change'
                        (e.g., 'a12g' for a nucleotide change from 'a' to 'g' at position 12).

                        If there are no differences, an empty tuple is returned.

        Raises:
            ValueError: If the codons are not of the same length.
        """
        if len(self.sequence) != len(other.sequence):
            raise ValueError("Codons must be of the same length.")

        differences = []
        for i in range(len(self.sequence)):
            if self.sequence[i] != other.sequence[i]:
                nucleotid_ref = f"{self.sequence[i].lower()}"
                nucleotid_sample = f"{other.sequence[i].lower()}"
                nucleotid_position = f"{self.position + i + 1}"
                difference = (
                    f"{nucleotid_ref}{nucleotid_position}{nucleotid_sample}"
                )
                differences.append(difference)

        return tuple(differences)

    def translate(self) -> str:
        """
        Translates the codon sequence into its corresponding amino acid.

        Returns:
            str: The single-letter code of the amino acid translated from the codon.
        """
        return str(Seq(self.sequence).translate())


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


class Mutation:
    """
    Represents a mutation between a reference codon and a sample codon.

    Attributes:
        ref_codon (Codon): The codon from the reference sequence.
        sample_codon (Codon): The codon from the sample sequence.
    """

    def __init__(self, ref_codon: Codon, sample_codon: Codon):
        """
        Initializes a Mutation object with reference and sample codons.

        Args:
            ref_codon (Codon): The codon from the reference sequence.
            sample_codon (Codon): The codon from the sample sequence.
        """
        self.ref_codon = ref_codon
        self.sample_codon = sample_codon

    def annotate(self) -> Tuple[str, Tuple[str]]:
        """
        Annotates the type of mutation and identifies the differences between the reference and sample codons.

        The method classifies mutations as 'synonymous', 'nonsynonymous', or 'unknown'.
        - Synonymous mutations do not change the amino acid produced by the codon.
        - Nonsynonymous mutations result in a different amino acid.
        - 'Unknown' is returned when the sample codon contains a gap ('-'), indicating a potential insertion or deletion event.

        Returns:
            Tuple[str, Tuple[str]]: A tuple containing:
                - str: The type of mutation ('synonymous', 'nonsynonymous', or 'unknown').
                - Tuple[str]: A tuple of strings that describe the differences between the reference and sample codons.
                            (e.g., 'a12g' for a nucleotide change from 'a' to 'g' at position 12).
                            (e.g., 'A12G' for an amino acid change from 'A' to 'G' at position 12).
                            (e.g., 'a12-' for a nucleotide change from 'a' to -' at position 12).
        """
        if "-" in self.sample_codon.sequence:
            nucleotid_differ = self.ref_codon.compare(self.sample_codon)
            return "unknown", nucleotid_differ

        ref_aminoacid = self.ref_codon.translate()
        sample_aminoacid = self.sample_codon.translate()

        if ref_aminoacid == sample_aminoacid:
            nucleotid_differ = self.ref_codon.compare(self.sample_codon)
            return "synonymous", nucleotid_differ
        else:
            aminoacid_differ = (
                f"{ref_aminoacid}{self.ref_codon.position}{sample_aminoacid}",
            )
            return "nonsynonymous", aminoacid_differ


class MutationService:
    """
    Service for analyzing mutations between a reference sequence and sample sequences.

    Attributes:
        ref_seq (Seq): The nucleotide sequence of the reference genome.
        genes (List[Gene]): A list of Gene objects representing the genes in the reference genome.
    """

    def __init__(self, seq_reference: SeqRecord):
        """
        Initializes a MutationService object with a reference sequence.

        Args:
            seq_reference (SeqRecord): A SeqRecord object containing the reference sequence and its features.
        """
        self.ref_seq = seq_reference.seq
        self.genes = self.__populate_genes(seq_reference)

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
                # print(gene_name, len(gene_seq) / 3) # How many codons in gene.
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
                                  categorized as 'synonymous', 'nonsynonymous', or 'unknown'.

        Raises:
            ValueError: If the length of the reference sequence and the sample sequence do not match.
        """
        if len(self.ref_seq) != len(seq_sample):
            raise ValueError(
                f"Sequences must be of equal length.\nReference sequence: {self.ref_seq}"
            )

        annotations = {}
        for gene in self.genes:
            if gene.name not in annotations:
                annotations[gene.name] = {
                    "synonymous": [],
                    "nonsynonymous": [],
                    "unknown": [],
                }

            print(
                f"Processing gene: {gene.name}, Numbers of codons: {len(gene.sequence) / 3}"
            )  # Print gene name an total codons for this gene.

            num_codons = len(gene.sequence) // 3
            for i in range(0, len(gene.sequence), 3):
                if print_progress:
                    progress = (i // 3) + 1
                    percent_complete = (progress / num_codons) * 100
                    print(
                        f"Processing codon {progress}/{num_codons} ({percent_complete:.2f}%)"
                    )

                ref_codon = Codon(gene.sequence[i : i + 3], i)
                sample_gene_secuence = seq_sample[
                    gene.location.start : gene.location.end
                ]
                sample_codon = Codon(sample_gene_secuence[i : i + 3], i)

                if ref_codon.sequence != sample_codon.sequence:
                    mutation = Mutation(ref_codon, sample_codon)
                    mutation_type, differences = mutation.annotate()
                    annotations[gene.name][mutation_type].extend(differences)
        print("")
        return annotations

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
            print(f"Processing sample: {seq_record_sample.id}")
            gene_mutation_details = self.__annotate_gene_mutations(seq_sample)
            annotations[seq_record_sample.id] = gene_mutation_details

        return annotations


path_ref_gbk: str = r"referencia.gbk"
path_muestras_fasta: str = r"muestras.fasta"

# Load reference genome
seq_record_reference: SeqRecord = SeqIO.read(
    path_ref_gbk, "genbank"
)  # Unlike FASTA, GenBank provides the specific position for each gene, which is necessary to locate the genes in the sequence.

# Load the sample genomes
with open(path_muestras_fasta) as sample_file:
    seq_record_samples = list(SeqIO.parse(sample_file, "fasta"))

# Get report
mutation_service = MutationService(seq_record_reference)
mutation_report = mutation_service.analyze_mutations(seq_record_samples)
print(json.dumps(mutation_report, indent=4))
