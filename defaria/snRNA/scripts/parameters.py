filters = {
    'min_counts': 200,
    'max_counts': 6000,
    'max_total_counts': 40000,
    'max_pct_counts_mt': 15
}

pc_parameters = {
    'n_pcs': 30,
    'n_neighbors': 20,
    'leiden_resolution': 0.8
    }


mnn_parameters = {
    'k': 20,
    "sigma": 1,
    "batch_key": "dataset",
    "do_concatenate": True
}