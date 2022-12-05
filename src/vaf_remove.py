import pandas as pd
import numpy as np
from scipy import stats
import sys

mutation_map_input_file = sys.argv[1]
vaf = float(sys.argv[2])
significance = float(sys.argv[3])
flag = sys.argv[4]

try:

    mutation_map = pd.read_csv(mutation_map_input_file, header=None, sep="\t")

    mutation_map = mutation_map.rename(columns=dict(zip(list(range(7)), ["chr", "pos", "ref", "alt", "context", "coverage", "alt_count"])))
    if mutation_map.shape[1] < 8:
        mutation_map["filter"] = "PASS"
    else:
        mutation_map = mutation_map.rename(columns={7: "filter"})

    try:
        p_value = stats.binom.cdf(k=mutation_map["alt_count"], n=mutation_map["coverage"], p=vaf)
        mutation_map["filter"] = np.where(p_value > significance, np.where(mutation_map["filter"] == "PASS", flag, mutation_map["filter"]+";"+flag), mutation_map["filter"])
    except TypeError:
        pass

    for m in range(mutation_map.shape[0]):
        print("\t".join(map(str, mutation_map.loc[m])))

except pd.errors.EmptyDataError:
    pass
