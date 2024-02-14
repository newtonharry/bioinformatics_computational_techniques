import pyfastx
from random import randint, seed
from numpy import random as nprandom
import math


protein_alphabet = "ACDEFGHIKLMNPQRSTVWY"

# Motif length
M = 12

seed(42)

# List of protein sequences
# Choose a random starting position for each sequence
sequences: list[tuple[str, int]] = [
    (seq, randint(0, len(seq) - M)) for _, seq in pyfastx.Fastx("./mbl_seqs.fa")
]
num_of_sequences = len(sequences)
previous_sequence_indexes: list[int] = []
left_out_sequence_index = 0
iterations = 5000
pseudocount = 1e-6


# This document contains the proper way to calculate the PPM
# https://www.bioconductor.org/packages/devel/bioc/vignettes/universalmotif/inst/doc/IntroductionToSequenceMotifs.pdf
def build_ppm(
    sequences: list[tuple[str, int]],
    left_out_sequence_index: int,
) -> dict[str, list[float]]:
    ppm = {base: [0.0] * M for base in protein_alphabet}

    for j, (seq, starting_pos) in enumerate(sequences):
        if j != left_out_sequence_index:
            for i, seq_pos in enumerate(range(starting_pos, starting_pos + M)):
                ppm[seq[seq_pos]][i] += 1

    # Normalize each column in the PPM to convert it from a frequency matrix to a probability matrix
    total_at_each_pos = [
        sum([ppm[base][i] for base in protein_alphabet]) for i in range(M)
    ]

    for i in range(M):
        for base in protein_alphabet:
            # The document specified says to add a pseudocount divided by the number of residues in the alphabet
            ppm[base][i] += pseudocount / len(protein_alphabet)
            # The document specified says to add a pseudocount to the denominator of the normalization
            ppm[base][i] /= (
                total_at_each_pos[i] + pseudocount
            )  # Normalize the residue counts at each position

    return ppm


# This document shows how to calculate the PWM
# https://www.bioconductor.org/packages/devel/bioc/vignettes/universalmotif/inst/doc/IntroductionToSequenceMotifs.pdf
# This function is to be used after using the PPM to assign the new starting motif position  in the left out sequence
# as the numpy random choice function can't take negative values (as produced by the log likelihood during PWM creation)
def ppm_to_pwm(
    ppm: dict[str, list[float]], background_probabilities: dict[str, float]
) -> dict[str, list[float]]:
    pwm = {base: [0.0] * M for base in protein_alphabet}
    for base in protein_alphabet:
        for i in range(M):
            # Take into consideration background residue probabilities in order to weigh
            # residue positions that that are less likely to occur than  background residue positions
            # and vice versa
            pwm[base][i] = math.log(ppm[base][i] / background_probabilities[base], 2)

    return pwm


def calculate_background_probabilities(
    sequences: list[tuple[str, int]], M: int
) -> dict[str, float]:
    background_frequencies = {base: 0.0 for base in protein_alphabet}
    for seq, start in sequences:
        for base in protein_alphabet:
            background_frequencies[base] += seq[:start].count(base) + seq[
                start + M :
            ].count(base)

    # Normalize the background frequencies
    total_background_probs = sum(background_frequencies.values())
    for base in protein_alphabet:
        background_frequencies[base] /= total_background_probs

    return background_frequencies


# This function is to be used after using the PPM to assign the new starting motif position  in the left out sequence
def generate_distribution_for_left_out_sequence(
    PWM: dict[str, list[float]],
    background_probabilities: dict[str, float],
    left_out_sequence: str,
) -> list[float]:
    distribution: list[float] = [1e-6] * (len(left_out_sequence) - M)
    for i in range(len(left_out_sequence) - M):
        for j, base in enumerate(left_out_sequence[i : i + M]):
            # Weigh the probability of the starting position of the motif using the background residue probabilities
            # We don't want to perform a log likelihood here as we are only interested in the relative probabilities
            # Negative values will not work with numpy's random choice function
            distribution[i] *= PWM[base][j] / background_probabilities[base]

    # Normalize the distribution
    distribution_sum = sum(distribution)
    distribution = list(map(lambda x: x / distribution_sum, distribution))

    return distribution


def choose_left_out_sequence(
    sequences: list[tuple[str, int]], previous_sequence_indexes: list[int]
) -> tuple[int, str]:
    left_out_sequence_index = randint(0, num_of_sequences - 1)

    # Wait until we get a sequence that hasn't been left out before
    while left_out_sequence_index in previous_sequence_indexes:
        left_out_sequence_index = randint(0, num_of_sequences - 1)
    return left_out_sequence_index, sequences[left_out_sequence_index][0]


def choose_new_position(distribution: list[float], left_out_sequence: str) -> int:
    # Choose a new position for the left out sequence weighted by the probability of each starting position of the motif
    draw = nprandom.choice(len(left_out_sequence) - M, 1, p=distribution)[0]
    return draw


# NOTE: START HERE to begin understanding the code
for iteration in range(iterations):
    # Choose a sequence to be left out
    left_out_sequence_index, left_out_sequence_string = choose_left_out_sequence(
        sequences, previous_sequence_indexes
    )

    # Calculate background probabilities for each base across all sequences,
    # only counting residues that are outside the motif for each sequence
    background_probabilities = calculate_background_probabilities(sequences, M)

    # Build a PPM from the n - 1 sequences (not including left out sequence)
    PPM = build_ppm(sequences, left_out_sequence_index)

    # Generate a distribution for the left out sequence
    distribution = generate_distribution_for_left_out_sequence(
        PPM, background_probabilities, left_out_sequence_string
    )

    # Choose a new position for the left out sequence given a weighted distribution of starting positions
    new_position = choose_new_position(distribution, left_out_sequence_string)

    # Update the left out sequence with its motif starting position
    sequences[left_out_sequence_index] = (
        sequences[left_out_sequence_index][0],
        new_position,
    )

    # Generate a PWM from the PPM
    PWM = ppm_to_pwm(PPM, background_probabilities)

    # Use the PWM to find the most likely sequence in the profile
    # as the PPM doesn't weigh profile positions with background residue probabilities
    most_likely_motif = "".join(
        [max(PWM, key=lambda base: PWM[base][i]) for i in range(M)]
    )

    LL = sum([max([PWM[base][i] for base in protein_alphabet]) for i in range(M)])

    print(
        f"Most likely motif at iteration {iteration} is {most_likely_motif} with log-likelihood of {LL} "
    )
