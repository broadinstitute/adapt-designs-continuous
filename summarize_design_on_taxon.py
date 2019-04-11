"""Functions for computing statistics from one or more designs.
"""

import argparse
import glob
import math
import statistics
import os

__author__ = 'Hayden Metsky <hayden@mit.edu>'


# Set the number of targets to read from each design
NUM_TARGETS = 10


class DesignTarget:
    """Store information on a design of a single target.
    """

    def __init__(self, target_start, target_end, guide_seqs,
            left_primer_seqs, right_primer_seqs, cost):
        self.target_start = target_start
        self.target_end = target_end
        self.target_length = target_end - target_start
        self.guide_seqs = tuple(sorted(guide_seqs))
        self.left_primer_seqs = tuple(sorted(left_primer_seqs))
        self.right_primer_seqs = tuple(sorted(right_primer_seqs))
        self.cost = cost

    def has_similar_endpoints(self, other, nearby_nt=20):
        """
        Determine if this target is similar to another based on coordinates.

        This compares the start/end coordinates of this target to
        another. Note that because coordinates depend on the alignment,
        this should only be compared against designs from the same
        or highly similar alignment (so coordinates are in the same
        space).

        An alternative way to compute this would be based on the left/right
        primer sequences: for example, by seeing if the longest common
        substring up to k mismatches is at least some length long.

        Args:
            other: DesignTarget object
            nearby_nt: say that other is similar to self if the start/end
                coordinates are within nearby_nt nt of each other

        Returns:
            True iff the start/end of other are both near the start/end
            of self
        """
        start_diff = math.fabs(self.target_start - other.target_start)
        end_diff = math.fabs(self.target_end - other.target_end)
        return start_diff <= nearby_nt and end_diff <= nearby_nt

    def __eq__(self, other):
        return (self.guide_seqs == other.guide_seqs and
                self.left_primer_seqs == other.left_primer_seqs and
                self.right_primer_seqs == other.right_primer_seqs)

    def __hash__(self):
        return (hash(self.guide_seqs) +
                hash(self.left_primer_seqs) +
                hash(self.right_primer_seqs))


class Design:
    """Store information on a design encompassing multiple possible targets.

    This stores the targets as a set (unordered). As a result, two designs
    can be equal if all their targets are equal, even if they are ordered
    differently (e.g., different costs for each).
    """

    def __init__(self, targets):
        self.best_target = targets[0]
        self.targets = frozenset(set(targets))

    def __eq__(self, other):
        return self.targets == other.targets

    def __hash__(self):
        return hash(self.targets)

    def jaccard_similarity(self, other):
        """Compute Jaccard similarity between this design and another.

        Args:
            other: Design object

        Returns:
            Jaccard similarity, by comparing targets, between self and other
        """
        intersection = self.targets & other.targets
        union = self.targets | other.targets
        return float(len(intersection)) / float(len(union))

    def jaccard_similarity_with_loose_equality(self, other):
        """Compute Jaccard similarity between this design and another
        by counting two designs as the same if their endpoint coordinates
        are almost identical.

        Note that this is currently highly inefficient, by comparing
        all pairs of DesignTargets. It might be more efficient to sort
        by coordinate and iterate through the combination of DesignTargets
        from both self and other.

        Also, note that designs should be produced from the same alignments
        (e.g., exactly the same inputs); otherwise, the coordinates
        between them may not be comparable.

        Args:
            other: design object

        Returns:
            Jaccard similarity, by comparing targets, between self and other
        """
        def find_unique_targets(targets):
            unique = set()
            for target in targets:
                # Check if an 'equivalent' target already exists
                already_exists = False
                for ut in unique:
                    if target.has_similar_endpoints(ut):
                        already_exists = True
                        break
                if not already_exists:
                    unique.add(target)
            return unique
        unique_self = find_unique_targets(self.targets)
        unique_other = find_unique_targets(other.targets)
        union = find_unique_targets(self.targets.union(other.targets))
        union_size = len(union)
        intersection_size = len(unique_self) + len(unique_other) - union_size
        return float(intersection_size) / float(union_size)

    @staticmethod
    def from_file(fn, num_targets=None):
        """Read a collection of targets from a file.

        Args:
            fn: path to a TSV file giving targets
            num_targets: only construct a Design from the top num_targets
                targets, as ordered by cost (if None, use all)

        Returns:
            Design object
        """
        rows = []
        with open(fn) as f:
            col_names = {}
            for i, line in enumerate(f):
                line = line.rstrip()
                ls = line.split('\t')
                if i == 0:
                    # Parse header
                    for j in range(len(ls)):
                        col_names[j] = ls[j]
                else:
                    # Read each column as a variable
                    cols = {}
                    for j in range(len(ls)):
                        cols[col_names[j]] = ls[j]
                    rows += [(cols['cost'], cols['target-start'],
                             cols['target-end'], cols)]

        # Sort rows by cost (first in the tuple); in case of ties, sort
        # by target start and target end positions (second and third in
        # the tuple)
        # Pull out the best N targets
        rows = sorted(rows)
        if num_targets != None:
            if len(rows) < num_targets:
                raise Exception(("The number of rows in a design (%d) is fewer "
                    "than the number of targets to read (%d)") %
                    (len(rows), num_targets))
            rows = rows[:num_targets]

        targets = []
        for row in rows:
            _, _, _, cols = row
            targets += [DesignTarget(
                int(cols['target-start']),
                int(cols['target-end']),
                cols['guide-target-sequences'].split(' '),
                cols['left-primer-target-sequences'].split(' '),
                cols['right-primer-target-sequences'].split(' '),
                float(cols['cost'])
            )]

        return Design(targets)


def read_designs(design_dir, timestamp):
    """Read designs from a timestamp for a taxonomy.

    Args:
        design_dir: path to directory containing designs for a taxonomy
        timestamp: timestamp of design to read

    Returns:
        list of Design objects (one per cluster); or None if no design
        exists
    """
    def fp(clust_num):
        return os.path.join(design_dir, str(timestamp),
                'design.tsv.' + str(clust_num))

    if not os.path.isfile(fp(0)):
        return None

    # Read all designs in this directory
    designs = []
    clust_num = 0
    while os.path.isfile(fp(clust_num)):
        designs += [Design.from_file(fp(clust_num), NUM_TARGETS)]
        clust_num += 1
    return designs


def read_most_recent_designs(design_dir):
    """Read the most recent design from a taxonomy directory.

    Args:
        design_dir: path to directory containing designs for a taxonomy

    Returns:
        tuple (t, d) where t is the timestamp of when the designs were
        created and d is list of Design objects (one per cluster)
    """
    # design_dir contains directories named by timestamp
    timestamps = [int(t) for t in os.listdir(design_dir) if
            os.path.isdir(os.path.join(design_dir, t))]
    timestamps = sorted(timestamps)

    # Step backward from the most recent timestamps to find one that
    # contains a design
    for t in timestamps[::-1]:
        designs = read_designs(design_dir, t)
        if designs is not None:
            return (t, designs)
    raise Exception("Could not find design inside %s" % design_dir)


def find_most_recent_time_designs_changed(design_dir, jaccard_thres):
    """Find the most recent time designs changed to produce the latest designs
    for a taxonomy.

    In particular, this finds the oldest timestamp (t) of designs that
    match the most recent designs. t indicates when the designs last
    changed.

    Args:
        design_dir: path to directory containing designs for a taxonomy
        jaccard_thres: value of Jaccard similarity, such that any designs
            with less than this Jaccard similarity are deemed to not
            be equal

    Returns:
        oldest timestamp of designs that 'equals' the most recent designs
    """
    _, most_recent = read_most_recent_designs(design_dir)

    # design_dir contains directories named by timestamp
    timestamps = [int(t) for t in os.listdir(design_dir) if
            os.path.isdir(os.path.join(design_dir, t))]
    timestamps = sorted(timestamps)
    timestamp_of_last_equal_design = timestamps[-1]

    # Step backward from the most recent timestamps to find the oldest
    # that 'equals' the most recent (but stop when finding one that differs;
    # thus, this may not really find the oldest, but rather than oldest
    # where there was a change and the designs matched the most recent
    # from there on)
    for t in timestamps[::-1]:
        designs = read_designs(design_dir, t)
        if designs is None:
            # Could not find designs for this timestamp
            continue
        if len(designs) != len(most_recent):
            # The number of clusters differs
            return timestamp_of_last_equal_design
        # Compare each cluster, assuming they are ordered the
        # same
        for i in range(len(designs)):
            js = designs[i].jaccard_similarity(most_recent[i])
            if js < jaccard_thres:
                # Cluster i differs
                return timestamp_of_last_equal_design
        # All clusters are the same between designs and most_recent
        timestamp_of_last_equal_design = t

    # All designs are equal to the most recent, so return the timestamp
    # of the oldest design
    return timestamp_of_last_equal_design


def main(args):
    # Read the most recent designs
    timestamp, designs = read_most_recent_designs(args.design_dir)

    # Determine when the designs last changed
    last_changed_timestamp = find_most_recent_time_designs_changed(
            args.design_dir, args.jaccard_thres)

    # Print a row for each cluster
    for i, design in enumerate(designs):
        mean_cost = statistics.mean(target.cost for target in design.targets)
        mean_num_primers5 = statistics.mean(len(target.left_primer_seqs)
                for target in design.targets)
        mean_num_primers3 = statistics.mean(len(target.right_primer_seqs)
                for target in design.targets)
        mean_num_guides = statistics.mean(len(target.guide_seqs)
                for target in design.targets)

        row = [i, timestamp, last_changed_timestamp,
                design.best_target.cost,
                design.best_target.target_length,
                len(design.best_target.left_primer_seqs),
                len(design.best_target.right_primer_seqs),
                len(design.best_target.guide_seqs)]
        print('\t'.join(str(x) for x in row))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('design_dir',
            help=("Path to directory containing designs (subdirectories "
                  "should be timestamps"))
    parser.add_argument('jaccard_thres',
            type=float,
            help=("Consider two designs whose Jaccard similarity to be "
                  "less than this to be different"))

    args = parser.parse_args()
    main(args)
