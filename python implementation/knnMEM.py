import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from joblib import Parallel, delayed


def zero_ref_func(exp_data, num_pops, num_markers):
    MAGref = np.zeros((num_pops, num_markers))
    # Compute IQR for each column
    IQRs = np.subtract(*np.percentile(exp_data, [75, 25], axis=0))
    medIQRref = np.median(IQRs)
    IQRref = np.full((num_pops, num_markers), medIQRref)
    return MAGref, IQRref

def IQR_thresh_func(MAGpop, MAGref, IQRpop, IQRref, num_markers):
    """
    Efficiently calculate an IQR threshold using universal IQR thresholding.
    """
    # Exclude the last column (assumed to be 'cluster' or non-marker)
    cols = slice(0, num_markers - 1)

    # Compute 25th percentiles for each marker (vectorized)
    magpop_q2 = np.quantile(MAGpop[:, cols], 0.25, axis=0)
    magref_q2 = np.quantile(MAGref[:, cols], 0.25, axis=0)

    # Boolean masks for values below or equal to 25th percentile (vectorized)
    MAGpop_belowThresh = MAGpop[:, cols] <= magpop_q2
    MAGref_belowThresh = MAGref[:, cols] <= magref_q2

    # For each marker, get the minimum IQR for those below threshold, else min of all
    IQR_thresh_pop = np.where(
        MAGpop_belowThresh.any(axis=0),
        np.where(
            MAGpop_belowThresh.any(axis=0),
            np.min(np.where(MAGpop_belowThresh, IQRpop[:, cols], np.inf), axis=0),
            np.min(IQRpop[:, cols], axis=0)
        ),
        np.min(IQRpop[:, cols], axis=0)
    )
    IQR_thresh_ref = np.where(
        MAGref_belowThresh.any(axis=0),
        np.where(
            MAGref_belowThresh.any(axis=0),
            np.min(np.where(MAGref_belowThresh, IQRref[:, cols], np.inf), axis=0),
            np.min(IQRref[:, cols], axis=0)
        ),
        np.min(IQRref[:, cols], axis=0)
    )

    # Combine and average
    IQR_thresh1 = np.mean(np.concatenate([IQR_thresh_pop, IQR_thresh_ref]))
    return IQR_thresh1

def minimal_MEM(
    exp_data,
    transform=False,
    cofactor=1,
    choose_ref=False,
    zero_ref=False,
    input_IQR_ref='',
    IQR_thresh=None,
    choose_markers=False,
    markers="all",
    rename_markers=False,
    new_marker_names="none",
    file_is_clust=False,
    add_fileID=False,
    output_prescaled_MEM=True,
    scale_matrix="linear",
    scale_factor=0
):
    
    #file_order = exp_data if file_is_clust else 0

    # Get markers to include in analysis

    # markerList = list(range(exp_data.shape[1]))
    # exp_data = exp_data.iloc[:, markerList]
    marker_names = list(exp_data.columns)



    # Initialize variables
    #marker_names = list(marker_names)
    num_markers = exp_data.shape[1]
    #exp_data.columns = marker_names
    num_cells = exp_data.shape[0]
    num_pops = exp_data['cluster'].nunique()
    pop_names = exp_data['cluster'].unique()

    MAGpop = np.zeros((num_pops, num_markers))
    MAGref = np.zeros((num_pops, num_markers))
    IQRpop = np.zeros((num_pops, num_markers))
    IQRref = np.zeros((num_pops, num_markers))



    # Get population medians and IQRs
    marker_cols = [col for col in exp_data.columns if col != 'cluster']
    grouped = exp_data.groupby('cluster')
    medians = grouped[marker_cols].median(numeric_only=True).abs()
    iqr = grouped[marker_cols].quantile(0.75) - grouped[marker_cols].quantile(0.25)

    MAGpop[:, :-1] = medians.loc[pop_names].values
    IQRpop[:, :-1] = iqr.loc[pop_names].values
    
    # for i, pop in enumerate(pop_names):
    #     subset = exp_data[exp_data['cluster'] == pop]
    #     MAGpop[i, :] = np.abs(subset.median(skipna=True))
    #     IQRpop[i, :] = subset.quantile(0.75) - subset.quantile(0.25)

    # Get reference population medians and IQRs

    if zero_ref:
        MAGref, IQRref = zero_ref_func(exp_data.iloc[:, :-1], num_pops, num_markers)
    # else:
    #     # Efficient leave-one-cluster-out medians and IQRs
    #     for i, pop in enumerate(pop_names):
    #         subset = exp_data[exp_data['cluster'] != pop]
    #         MAGref[i, :] = subset.median(numeric_only=True).abs().values
    #         IQRref[i, :] = (subset.quantile(0.75, numeric_only=True) - subset.quantile(0.25, numeric_only=True)).values
    marker_cols = [col for col in exp_data.columns if col != 'cluster']
    
    data = exp_data[marker_cols].values
    clusters = exp_data['cluster'].values
    def leave_one_out_stats(pop):
        mask = clusters != pop
        subset = data[mask]
        med = np.abs(np.median(subset, axis=0))
        iqr = np.percentile(subset, 75, axis=0) - np.percentile(subset, 25, axis=0)
        return med, iqr

    results = Parallel(n_jobs=-1)(
        delayed(leave_one_out_stats)(pop) for pop in pop_names
    )
    MAGref[:, :-1], IQRref[:, :-1] = map(np.vstack, zip(*results))

    # User supplied IQR reference value
    if input_IQR_ref != '':
        if not zero_ref:
            print("Warning: You can only specify an IQR reference value when using zero_ref MEM.")
        if not isinstance(input_IQR_ref, (int, float)) or np.size(input_IQR_ref) != 1:
            print("Warning: Please supply a single number to be used as the IQR reference value.")
        IQRref = np.full((num_pops, num_markers), input_IQR_ref)
    else:
        if num_pops == 1:
            print("Warning: You are attempting to run zero_ref MEM on one cluster. Consider specifying input_IQR_ref.")

    # Set and apply IQR threshold
    if IQR_thresh is None or num_pops < 4:
        IQR_thresh = 0.5
    if IQR_thresh == "auto":
        IQR_thresh = IQR_thresh_func(MAGpop, MAGref, IQRpop, IQRref, num_markers)

    # Vectorized thresholding for IQRpop and IQRref 
    IQRpop[:, :-1] = np.maximum(IQRpop[:, :-1], IQR_thresh)
    IQRref[:, :-1] = np.maximum(IQRref[:, :-1], IQR_thresh)

    # Calculate the second term in the MEM equation
    IQRcomp = np.divide(IQRref, IQRpop, out=np.zeros_like(IQRref), where=IQRpop != 0) - 1
    IQRcomp[IQRcomp < 0] = 0

    if zero_ref:
        IQRcomp[MAGpop < 1.44] = 0

    # Calculate MEM scores
    MAG_diff = MAGpop - MAGref
    MEM_matrix = np.abs(MAGpop - MAGref) + IQRcomp
    MEM_matrix[~(MAG_diff >= 0)] = -MEM_matrix[~(MAG_diff >= 0)]
    if zero_ref:
        MEM_matrix[MEM_matrix < 0] = 0

    # Put MEM values on -10 to +10 scale
    prescaled_MEM_matrix = MEM_matrix.copy()
    # if scale_matrix == "linear":
    #     scale_max = np.max(np.abs(MEM_matrix[:, :-1]))
    #     if num_pops == 1 and zero_ref:
    #         MEM_matrix = np.c_[((MEM_matrix[:, :-1].T / scale_max) * 10).T, MEM_matrix[:, -1]]
    #     else:
    #         MEM_matrix = np.c_[(MEM_matrix[:, :-1] / scale_max) * 10, MEM_matrix[:, -1]]
    # elif scale_matrix == "log":
    #     scaled_matrix = np.log(MEM_matrix) / np.exp(scale_factor)
    #     scale_max = np.max(np.abs(scaled_matrix[:, :-1]))
    #     scale_min = np.min(np.abs(scaled_matrix[:, :-1]))
    #     if num_pops == 1 and zero_ref:
    #         MEM_matrix = np.c_[((scaled_matrix[:, :-1].T - scale_min) / scale_max * 10).T, scaled_matrix[:, -1]]
    #     else:
    #         MEM_matrix = np.c_[(scaled_matrix[:, :-1] - scale_min) / scale_max * 10, scaled_matrix[:, -1]]
    # elif scale_matrix == "arcsinh":
    #     scaled_matrix = np.arcsinh(MEM_matrix / scale_factor)
    #     scale_max = np.max(np.abs(scaled_matrix[:, :-1]))
    #     scale_min = np.min(np.abs(scaled_matrix[:, :-1]))
    #     if num_pops == 1 and zero_ref:
    #         MEM_matrix = np.c_[((scaled_matrix[:, :-1].T - scale_min) / scale_max * 10).T, scaled_matrix[:, -1]]
    #     else:
    #         MEM_matrix = np.c_[(scaled_matrix[:, :-1] - scale_min) / scale_max * 10, scaled_matrix[:, -1]]

    # Rename rows and columns of all matrices
    def rename_table(x):
        df = pd.DataFrame(x)
        df.columns = marker_names[:-1]
        df.index = pop_names
        return df

    # if num_pops == 1 and zero_ref:
    #     object_list_labeled = [
    #         rename_table(MAGpop[:, :-1].T),
    #         rename_table(MAGref[:, :-1].T),
    #         rename_table(IQRpop[:, :-1].T),
    #         rename_table(IQRref[:, :-1].T),
    #         rename_table(MEM_matrix[:, :-1].T),
    #         rename_table(prescaled_MEM_matrix[:, :-1].T)
    #     ]
    #else:
    object_list_labeled = [
            # rename_table(MAGpop[:, :-1]),
            # rename_table(MAGref[:, :-1]),
            # rename_table(IQRpop[:, :-1]),
            # rename_table(IQRref[:, :-1]),
            # rename_table(MEM_matrix[:, :-1]),
            rename_table(prescaled_MEM_matrix[:, :-1])
        ]
    #object_list_labeled.append(file_order)

    all_values = {
        # "MAGpop": object_list_labeled[0],
        # "MAGref": object_list_labeled[1],
        # "IQRpop": object_list_labeled[2],
        # "IQRref": object_list_labeled[3],
        # "MEM_matrix": object_list_labeled[4],
        "prescaled_MEM_matrix": object_list_labeled[0],
        # "File Order": object_list_labeled[6]
    }


    # Export pre-scaled MEM values
    # if output_prescaled_MEM:
    #     prescaled_df = pd.DataFrame(prescaled_MEM_matrix, columns=marker_names)
    #     os.makedirs("output files", exist_ok=True)
    #     prescaled_df.iloc[:, :-1].to_csv(
    #         f"./output files/{pd.Timestamp.now().strftime('%Y-%m-%d_%H%M%S')} Pre-scaled MEM matrix.txt",
    #         sep="\t"
    #     )

    return all_values



def minimal_process_row_joblib(r, neighbor_index, transformed_data, ref_IQR):
    cell_indices = neighbor_index[r, :]
    MEM_data = transformed_data.iloc[cell_indices, :].copy()
    MEM_data['cluster'] = r
    MEM_values = minimal_MEM(MEM_data, zero_ref=True, input_IQR_ref=ref_IQR)
    return MEM_values["prescaled_MEM_matrix"].iloc[0].values

def minimal_knn_MEM_joblib(tsne_data, transformed_data, kvalue=60, n_jobs=-1):
    """
    Parallel Python version of knn_MEM using joblib.
    """
    ref_IQR = np.median(np.subtract(*np.percentile(transformed_data.values, [75, 25], axis=0)))
    nbrs = NearestNeighbors(n_neighbors=kvalue, algorithm='ball_tree').fit(tsne_data)
    neighbor_index = nbrs.kneighbors(tsne_data, return_distance=False)
    n_rows = neighbor_index.shape[0]
    
    results = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(minimal_process_row_joblib)(r, neighbor_index, transformed_data, ref_IQR)
        for r in range(n_rows)
    )
    
    all_labels = np.vstack(results)
    scale_max = np.max(np.abs(all_labels))
    MEM_scores = (all_labels / scale_max) * 10
    # Use columns from the last MEM_values (safe if all are the same)
    MEM_scores_df = pd.DataFrame(MEM_scores, columns=minimal_MEM(
        transformed_data.iloc[:kvalue, :].assign(cluster=0), zero_ref=True, input_IQR_ref=ref_IQR
    )["prescaled_MEM_matrix"].columns)
    return MEM_scores_df