"""Bar-chart plots of the top candidate genes."""
import matplotlib
matplotlib.use("Agg")  # no display needed, just save PNG files
import matplotlib.pyplot as plt


def plot_top_genes(genes, score_fn, title, ylabel, out_path, top_n=10):
    top = sorted(genes, key=score_fn, reverse=True)[:top_n]
    symbols = [g.symbol for g in top]
    scores = [score_fn(g) for g in top]

    plt.figure(figsize=(10, 6))
    plt.bar(symbols, scores, color="#4C72B0")
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel("Gene")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
