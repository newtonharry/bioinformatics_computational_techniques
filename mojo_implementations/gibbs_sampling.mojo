from python import Python
import random
from collections import Dict, KeyElement


@value
struct Sequence(CollectionElement, Stringable):
    var sequence: String
    var motif_starting_index: UInt64

    fn __init__(inout self, sequence: String, index: UInt64):
        self.sequence = sequence
        self.motif_starting_index = index

    fn __str__(self) -> String:
        return self.sequence

    fn __getitem__(self, i: Int) -> String:
        return self.sequence[i]


struct Sequences:
    var data: DynamicVector[Sequence]

    fn __init__(inout self):
        self.data = DynamicVector[Sequence]()

    fn add_element(inout self, a: Sequence) -> None:
        self.data.append(a)

    fn clear(inout self):
        self.data.clear()

    fn length(self) -> UInt64:
        return self.data.__len__()

    fn contains(self, a: Sequence) -> Int:
        for i in range(self.length()):
            if self.data[i].sequence == a.sequence:
                return 1
        return 0

    fn __getitem__(self, i: Int) -> Sequence:
        return self.data[i]


@value
struct StringKey(KeyElement):
    var s: String

    fn __init__(inout self, owned s: String):
        self.s = s ^

    fn __init__(inout self, s: StringLiteral):
        self.s = String(s)

    fn __hash__(self) -> Int:
        let ptr = self.s._buffer.data.value
        return hash(DTypePointer[DType.int8](ptr), len(self.s))

    fn __eq__(self, other: Self) -> Bool:
        return self.s == other.s


struct Profile:
    var matrix: Dict[StringKey, Tensor[DType.int64]]
    var alphabet: String

    fn __init__(inout self, alphabet: StringLiteral, length: Int):
        self.alphabet = String(alphabet)
        self.matrix = Dict[StringKey, Tensor[DType.int64]]()
        for char in range(len(alphabet)):
            self.matrix[self.alphabet[char]] = Tensor[DType.int64](
                self.alphabet.__len__(), length
            )

    fn __getitem__(self, key: StringKey) raises -> Reference[Tensor[DType.int64]]:
        return Reference(self.matrix.find(key).value())

    fn __moveinit__(inout self, owned other: Self):
        self.alphabet = other.alphabet
        self.matrix = Dict[StringKey, Tensor[DType.int64]]()
        for char in range(len(other.alphabet)):
            self.matrix[self.alphabet[char]] = (
                other.matrix.find(other.alphabet[char]).or_else(0).value()
            )

    fn __copyinit__(inout self, other: Self):
        self.alphabet = other.alphabet
        self.matrix = Dict[StringKey, Tensor[DType.int64]]()
        for char in range(len(other.alphabet)):
            self.matrix[self.alphabet[char]] = (
                other.matrix.find(other.alphabet[char]).or_else(0).value()
            )

    fn __str__(self):
        return


fn calculate_background_probabilities(
    borrowed sequences: Sequences, motif_length: Int
) -> Profile:
    let PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
    var background_probabilities = Profile(PROTEIN_ALPHABET, 1)
    return background_probabilities ^


fn build_ppm(
    borrowed sequences: Sequences, left_out_sequence_index: UInt64, motif_length: Int
) raises -> Profile:
    let PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
    var ppm = Profile(PROTEIN_ALPHABET, motif_length)
    for i in range(len(sequences.length())):
        if i != left_out_sequence_index.to_int():
            let sequence = sequences[i]
            let starting_index = sequences[i].motif_starting_index
            let offset = 0
            for j in range(starting_index, starting_index + motif_length):
                ppm[sequence[j]][offset] += 1


# Import pyfastx to process fasta file
fn main() raises:
    let pyfastx = Python.import_module("pyfastx")
    var protein_sequences = pyfastx.Fastx("./mbl_seqs.fa")
    var sequences = Sequences()
    for sequence in protein_sequences:
        let py_str = str(sequence[1])
        let random_index = random.random_ui64(0, py_str.__len__())
        let seq = Sequence(py_str, random_index)
        sequences.add_element(seq)

    let iterations = 1000
    let motif_length = 12
    let seqs_len = sequences.length()

    let left_out_sequence_index = random.random_ui64(0, sequences.length() - 1)
    let background_probabilities = calculate_background_probabilities(sequences, 1)
    let PPM = build_ppm(sequences, left_out_sequence_index, motif_length)
