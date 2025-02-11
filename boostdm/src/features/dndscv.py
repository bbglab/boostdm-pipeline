
import pandas as pd


def excess(w):
    if w <= 1:
        return 0
    else:
        return (w - 1) / w


def filter(df, dndscv, xs_thresh):
    dnds = pd.read_csv(dndscv, sep="\t")
    # TODO do not make an outer merge
    merged = df.merge(
        dnds[['gene_name', 'wmis_cv', 'wnon_cv', 'wspl_cv', 'pallsubs_cv', 'qallsubs_cv']],
        left_on='gene', right_on='gene_name', how='outer').dropna()

    # select all mutations in genes with high excess (from dNdS) for the respective csqn type
    merged['xs_mis'] = merged['wmis_cv'].apply(excess)
    merged['xs_non'] = merged['wnon_cv'].apply(excess)
    merged['xs_spl'] = merged['wspl_cv'].apply(excess)
    mask_mis = (merged['xs_mis'] > xs_thresh) & (merged['impact'] == 'Missense')
    mask_non = (merged['xs_non'] > xs_thresh) & (merged['impact'] == 'Nonsense')
    mask_spl = (merged['xs_spl'] > xs_thresh) & (merged['impact'] == 'Essential_Splice')

    merged = merged[mask_mis | mask_non | mask_spl]

    return merged
