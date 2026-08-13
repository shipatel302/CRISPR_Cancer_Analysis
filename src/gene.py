class Gene:
    """Tracks how essential one gene looks across a set of CRISPR screens."""

    def __init__(self, identifier_id, symbol=None, total_screens=1):
        self.identifier_id = identifier_id
        self.symbol = symbol
        self.scores = []          # Bayes Factor from each screen where the gene was a hit
        self.repetitions = []     # screen IDs where the gene was a hit
        self.total_screens = total_screens

    def add_score(self, score):
        self.scores.append(score)

    def add_repetition(self, screen_id):
        self.repetitions.append(screen_id)

    def num_repetitions(self):
        return len(self.repetitions)

    def avg_occurance(self):
        """Fraction of all screens where this gene was flagged essential."""
        if self.repetitions:
            return len(self.repetitions) / self.total_screens
        return 0

    def average_score(self):
        """Mean Bayes Factor across the screens where it was a hit."""
        if self.scores:
            return sum(self.scores) / len(self.scores)
        return 0

    def combined_score(self):
        """Average score weighted by how consistently the gene shows up.
        Rewards genes that are both strongly and repeatedly essential."""
        return self.avg_occurance() * self.average_score()
