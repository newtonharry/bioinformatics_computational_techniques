from typing import List, Dict
from pprint import pprint
from random import choice


def initialize_pfm(
    sequences: List[str], motif_length: int
) -> (Dict[str, List[float]], Dict[str, float]):
    """
    Initialize the Position Frequency Matrix (PFM) and background probabilities.

    The PFM represents the frequency of each base at each position within a motif.
    The background probabilities represent the frequency of each base in the entire set of sequences.
    """

    # Initialize PFM and background probabilities with zeros
    pfm = {base: [0] * motif_length for base in "ATGC"}
    background = {base: 0 for base in "ATGC"}

    # Calculate the total number of bases across all sequences
    total_bases = sum(len(seq) for seq in sequences)

    # Calculate the weight for each potential motif site
    weight = 1 / (len(sequences[0]) - motif_length + 1)

    # Calculate initial PFM based on a uniform distribution of potential motif sites
    for sequence in sequences:
        for start in range(len(sequence) - motif_length + 1):
            site = sequence[start : start + motif_length]
            for i, base in enumerate(site):
                pfm[base][i] += weight

    # Normalize the PFM such that the probabilities at each position sum to 1
    for i in range(motif_length):
        total_at_position = sum(pfm[base][i] for base in "ATGC")
        for base in "ATGC":
            pfm[base][i] /= total_at_position

    # Calculate the background probabilities for each base across all sequences
    for sequence in sequences:
        for base in "ATGC":
            background[base] += sequence.count(base)
    for base in "ATGC":
        background[base] /= total_bases

    return pfm, background


def motif_probability(site: str, pfm: Dict[str, List[float]]) -> float:
    """
    Calculate the probability of the binding site using the PFM.

    This is essentially the likelihood of observing a given sequence based on the current PFM.
    """
    probability = 1e-6
    for i, base in enumerate(site):
        probability *= pfm[base][i]
    return probability


def background_probability(site: str, background: Dict[str, float]) -> float:
    """
    Calculate the background probability for a given site.

    This represents the likelihood of observing a sequence based solely on the background frequencies,
    without any specific motif pattern.
    """
    probability = 1e-6
    for base in site:
        probability *= background[base]
    return probability


def calculate_expectation(
    sequences: List[str],
    pfm: Dict[str, List[float]],
    background: Dict[str, float],
    motif_length: int,
) -> List[List[float]]:
    """
    Calculate the expected motif locations using the current PFM and background probabilities.

    For each potential motif site in each sequence, this function calculates:
    1. The probability that the site is a motif (binding probability).
    2. The probability that the site is just background (background probability).

    The final output is a weighted probability that each site is a motif.
    """
    total_binding_probs = []
    total_background_probs = []

    for sequence in sequences:
        binding_probs, background_probs = [], []
        for start in range(len(sequence) - motif_length + 1):
            site = sequence[start : start + motif_length]
            binding_probs.append(motif_probability(site, pfm))
            background_probs.append(background_probability(site, background))

        total_binding_probs.append(binding_probs)
        total_background_probs.append(background_probs)

    weighted_binding_probs = []
    for binding_probs, background_probs in zip(
        total_binding_probs, total_background_probs
    ):
        weighted_probs = [
            bp / (bp + bg) for bp, bg in zip(binding_probs, background_probs)
        ]
        total = sum(weighted_probs)
        normalized_probs = [wp / total for wp in weighted_probs]
        weighted_binding_probs.append(normalized_probs)

    return weighted_binding_probs


def maximize_expectation(
    sequences: List[str],
    background: Dict[str, float],
    weighted_probs: List[List[float]],
    motif_length: int,
) -> (Dict[str, List[float]], Dict[str, float]):
    """
    Adjust the PFM and background probabilities to maximize the likelihood of the observed data.

    This step updates the PFM and background based on the weighted motif probabilities
    calculated in the expectation step.
    """
    pfm = {base: [0] * motif_length for base in "ATGC"}
    new_background = {base: 0 for base in "ATGC"}

    for seq_id, sequence in enumerate(sequences):
        for start in range(len(sequence) - motif_length + 1):
            site = sequence[start : start + motif_length]
            weights = weighted_probs[seq_id]
            for i, base in enumerate(site):
                pfm[base][i] += weights[start]
        for base in "ATGC":
            new_background[base] += sequence.count(base)

    # Normalize the PFM and background probabilities
    for base in "ATGC":
        for i in range(motif_length):
            pfm[base][i] /= len(sequences)
    total = sum(new_background.values())
    for base in "ATGC":
        new_background[base] /= total

    return pfm, new_background


def has_converged(
    old_pfm: Dict[str, List[float]], new_pfm: Dict[str, List[float]], motif_length: int
) -> bool:
    """
    Check if the EM algorithm has converged.

    Convergence is determined by the Euclidean distance between the old and new PFMs.
    If the change is smaller than a set threshold, the algorithm is considered to have converged.
    """
    return (
        sum(
            (new_pfm[base][i] - old_pfm[base][i]) ** 2
            for base in old_pfm
            for i in range(motif_length)
        )
        < 1e-6
    )


def expectation_maximization(sequences: List[str], motif_length: int):
    """
    Run the Expectation-Maximization (EM) algorithm to find motifs within sequences.

    The EM algorithm iteratively estimates the most likely motifs and their positions
    within a set of sequences. It does this by:
    1. Estimating the expected positions of motifs based on the current model (Expectation).
    2. Adjusting the model to fit the data as closely as possible (Maximization).
    """
    pfm, background = initialize_pfm(sequences, motif_length)
    while True:
        old_pfm = pfm.copy()
        weighted_probs = calculate_expectation(sequences, pfm, background, motif_length)
        pfm, background = maximize_expectation(
            sequences, background, weighted_probs, motif_length
        )
        if has_converged(old_pfm, pfm, motif_length):
            break
    pprint(background)


# Generate random sequences for testing
sequences = ["".join(choice("ATGC") for _ in range(3000)) for _ in range(10)]

# Running the corrected algorithm
expectation_maximization(sequences, 45)
