from ..graph import (
    overlapping_patterns, de_bruijn_path, de_bruijn_graph, parse_graph,
    eulerian_cycle, path_to_graph, euler_path, build_di_graph,
    reconstruct_string, universal_circular_string, is_universal
)
from ._lib import find_datasets


def sorted_graph(graph: dict) -> dict:
    return {x: sorted(y) for x, y in graph.items()}


def test_overlapping_patterns():
    dataset = [
        ("ATGCG GCATG CATGC AGGCA GGCAT",
         "AGGCA -> GGCAT CATGC -> ATGCG GCATG -> CATGC GGCAT -> GCATG")
    ] + find_datasets("overlapping_patterns")

    for input_sample, output_sample in dataset:
        patterns = input_sample.split()
        output_sample = parse_graph(output_sample)
        output = overlapping_patterns(patterns)
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def test_de_bruijn_path():
    dataset = [
        ((4, "AAGATTCTCTAAGA"),
         "AAG -> AGA,AGA AGA -> GAT ATT -> TTC CTA -> TAA CTC -> TCT GAT -> "
         "ATT TAA -> AAG TCT -> CTA,CTC TTC -> TCT")
    ] + find_datasets("de_bruijn_path", lambda x: x.splitlines())

    for input_sample, output_sample in dataset:
        k, text = input_sample
        k = int(k)
        output_sample = parse_graph(output_sample)
        output = sorted_graph(de_bruijn_path(k, text))
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def test_de_bruijn_graph():
    dataset = [
        (("GAGG CAGG GGGG GGGA CAGG AGGG GGAG"),
         "AGG -> GGG CAG -> AGG,AGG GAG -> AGG GGA -> GAG GGG -> GGA,GGG")
    ] + find_datasets("de_bruijn_graph")

    for input_sample, output_sample in dataset:
        kmers = input_sample.split()
        output_sample = parse_graph(output_sample)
        output = sorted_graph(de_bruijn_graph(kmers))
        assert output == output_sample, "{0} != {1}".format(output, output_sample)


def test_eulerian_cycle():
    dataset = [
        (("0 -> 3 1 -> 0 2 -> 1,6 3 -> 2 4 -> 2 5 -> 4 6 -> 5,8 7 -> 9 "
          "8 -> 7 9 -> 6"), "6->8->7->9->6->5->4->2->1->0->3->2->6")
    ] + find_datasets("eulerian_cycle")

    for input_sample, output_sample in dataset:
        graph = parse_graph(input_sample)
        output = "->".join(eulerian_cycle(graph))
        assert path_to_graph(output_sample) == path_to_graph(output), \
            "{0} != {1}".format(output, output_sample)


def test_euler_path():
    dataset = [
        ("0 -> 2 1 -> 3 2 -> 1 3 -> 0,4 6 -> 3,7 7 -> 8 8 -> 9 9 -> 6",
         "6->7->8->9->6->3->0->2->1->3->4")
    ] + find_datasets("euler_path")

    for input_sample, output_sample in dataset:
        graph = build_di_graph(parse_graph(input_sample))
        output = euler_path(graph)
        output_sample = output_sample.split("->")
        assert path_to_graph(output_sample) == path_to_graph(output), \
            "{0} != {1}".format(output, output_sample)


def test_reconstruct_string():
    def parse_input(text):
        lines = text.splitlines()
        return lines[0], "\n".join(lines[1:])

    dataset = [
        ((4, "CTTA ACCA TACC GGCT GCTT TTAC"), "GGCTTACCA")
    ] + find_datasets("reconstruct_string", parse_input)

    for input_sample, output_sample in dataset:
        k, kmers = input_sample
        k, kmers = int(k), kmers.split()
        output = reconstruct_string(kmers, k)
        assert output_sample == output, \
            "{0} != {1}".format(output, output_sample)


def test_universal_circular_string():
    dataset = [
        (4, "0000110010111101")
    ] + find_datasets("universal_circular_string")

    for input_sample, output_sample in dataset:
        k = int(input_sample)
        output = universal_circular_string(k)
        assert is_universal(output, k) and is_universal(output_sample, k), \
            "{0} != {1}".format(output, output_sample)


def main():
    #test_overlapping_patterns()
    test_de_bruijn_path()
    test_de_bruijn_graph()
    test_eulerian_cycle()
    test_euler_path()
    test_reconstruct_string()
    test_universal_circular_string()
    print("Success!")

if __name__ == '__main__':
    main()
