import numpy as np
import pandas as pd
import os


def choose_markers_func(exp_data):
    print("Numbered column names, in order they appear in file:")
    marker_cols = list(exp_data.columns[:-1])
    for i, col in enumerate(marker_cols):
        print(f"{i+1}: {col}")
    markers1 = input("Enter column numbers to include (e.g. 1:5,6,8:10):\n")
    sep_vals = [val.strip() for val in markers1.split(",")]
    list_vals = []
    for val in sep_vals:
        if ":" in val:
            start, end = map(int, val.split(":"))
            new_vals = list(range(start, end + 1))
            list_vals.extend(new_vals)
        elif val:
            list_vals.append(int(val))
    # Convert to zero-based indices and always include last column (cluster)
    markerList = [v - 1 for v in list_vals] + [exp_data.shape[1] - 1]
    print("Selected marker indices (zero-based):", markerList)

    return markerList  # Retrieve output.markerList after selection

def rename_markers_func(exp_data, marker_names):
    # Show current marker names (excluding 'cluster')
    marker_cols = list(exp_data.columns[:-1])
    print("Current marker names (excluding 'cluster'):")
    print(", ".join(marker_cols))
    user_input_names = input(
        "Enter new marker names, in same order as above, separated by commas.\nNo spaces allowed in name.\n"
    )
    new_marker_names = [name.strip() for name in user_input_names.split(",")]

    if len(new_marker_names) != (len(marker_names) - 1):
        print("Warning: Number of new marker names does not match number of markers.")
        # Recursively prompt again
        return rename_markers_func(exp_data, marker_names)

    # Add 'cluster' as last column
    new_marker_names.append("cluster")
    return new_marker_names

def zero_ref_func(exp_data, num_pops, num_markers):
    MAGref = np.zeros((num_pops, num_markers))
    # Compute IQR for each column
    IQRs = np.subtract(*np.percentile(exp_data, [75, 25], axis=0))
    medIQRref = np.median(IQRs)
    IQRref = np.full((num_pops, num_markers), medIQRref)
    return MAGref, IQRref

def IQR_thresh_func(MAGpop, MAGref, IQRpop, IQRref, num_markers):
    """
    Automatically calculate an IQR threshold using universal IQR thresholding.
    """
    IQR_thresh_pop = []
    IQR_thresh_ref = []
    for i in range(num_markers - 1):
        # 2nd quantile (25th percentile) for column i
        magpop_q2 = np.quantile(MAGpop[:, i], 0.25)
        magref_q2 = np.quantile(MAGref[:, i], 0.25)
        # Boolean mask for values below or equal to 25th percentile
        MAGpop_belowThresh = MAGpop[:, i] <= magpop_q2
        MAGref_belowThresh = MAGref[:, i] <= magref_q2
        # Minimum IQR for those below threshold
        if np.any(MAGpop_belowThresh):
            IQR_thresh_pop.append(np.min(IQRpop[:, i][MAGpop_belowThresh]))
        else:
            IQR_thresh_pop.append(np.min(IQRpop[:, i]))
        if np.any(MAGref_belowThresh):
            IQR_thresh_ref.append(np.min(IQRref[:, i][MAGref_belowThresh]))
        else:
            IQR_thresh_ref.append(np.min(IQRref[:, i]))
    IQR_thresh1 = np.mean(IQR_thresh_pop + IQR_thresh_ref)
    return IQR_thresh1

def pythonMEM(
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
    
    file_order = exp_data if file_is_clust else 0

     # Get markers to include in analysis
    marker_names = list(exp_data.columns[:-1]) + ["cluster"]
    if choose_markers:
        markerList = choose_markers_func(exp_data)
    elif not choose_markers and markers == "all":
        markerList = list(range(exp_data.shape[1]))
    else:
        sep_vals = markers.split(",")
        list_vals = []
        for val in sep_vals:
            if ":" in val:
                start, end = map(int, val.split(":"))
                list_vals.extend(range(start, end + 1))
            else:
                list_vals.append(int(val))
        markerList = list_vals + [exp_data.shape[1] - 1]

    exp_data = exp_data.iloc[:, markerList]
    marker_names = list(exp_data.columns)

    # Rename markers
    if rename_markers:
        new_marker_names = rename_markers_func(exp_data, marker_names)
    elif not rename_markers and new_marker_names == "none":
        new_marker_names = marker_names

    marker_names = new_marker_names

    # Initialize variables
    marker_names = list(marker_names)
    num_markers = exp_data.shape[1]
    exp_data.columns = marker_names
    num_cells = exp_data.shape[0]
    num_pops = exp_data['cluster'].nunique()
    pop_names = exp_data['cluster'].unique()

    MAGpop = np.zeros((num_pops, num_markers))
    MAGref = np.zeros((num_pops, num_markers))
    IQRpop = np.zeros((num_pops, num_markers))
    IQRref = np.zeros((num_pops, num_markers))

    # Transform values if specified
    if transform:
        exp_data.iloc[:, :-1] = np.arcsinh(exp_data.iloc[:, :-1] / cofactor)

    # Get population medians and IQRs
    for i, pop in enumerate(pop_names):
        subset = exp_data[exp_data['cluster'] == pop]
        MAGpop[i, :] = np.abs(subset.median(skipna=True))
        IQRpop[i, :] = subset.quantile(0.75) - subset.quantile(0.25)

    # Get reference population medians and IQRs
    if choose_ref:
        raise NotImplementedError("choose_ref logic not implemented")
    elif zero_ref:
        MAGref, IQRref = zero_ref_func(exp_data.iloc[:, :-1], num_pops, num_markers)
    else:
        for i, pop in enumerate(pop_names):
            subset = exp_data[exp_data['cluster'] != pop]
            MAGref[i, :] = np.abs(subset.median(skipna=True))
            IQRref[i, :] = subset.quantile(0.75) - subset.quantile(0.25)

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

    for i in range(num_markers - 1):
        IQRpop[:, i] = np.maximum(IQRpop[:, i], IQR_thresh)
        IQRref[:, i] = np.maximum(IQRref[:, i], IQR_thresh)

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
    if scale_matrix == "linear":
        scale_max = np.max(np.abs(MEM_matrix[:, :-1]))
        if num_pops == 1 and zero_ref:
            MEM_matrix = np.c_[((MEM_matrix[:, :-1].T / scale_max) * 10).T, MEM_matrix[:, -1]]
        else:
            MEM_matrix = np.c_[(MEM_matrix[:, :-1] / scale_max) * 10, MEM_matrix[:, -1]]
    elif scale_matrix == "log":
        scaled_matrix = np.log(MEM_matrix) / np.exp(scale_factor)
        scale_max = np.max(np.abs(scaled_matrix[:, :-1]))
        scale_min = np.min(np.abs(scaled_matrix[:, :-1]))
        if num_pops == 1 and zero_ref:
            MEM_matrix = np.c_[((scaled_matrix[:, :-1].T - scale_min) / scale_max * 10).T, scaled_matrix[:, -1]]
        else:
            MEM_matrix = np.c_[(scaled_matrix[:, :-1] - scale_min) / scale_max * 10, scaled_matrix[:, -1]]
    elif scale_matrix == "arcsinh":
        scaled_matrix = np.arcsinh(MEM_matrix / scale_factor)
        scale_max = np.max(np.abs(scaled_matrix[:, :-1]))
        scale_min = np.min(np.abs(scaled_matrix[:, :-1]))
        if num_pops == 1 and zero_ref:
            MEM_matrix = np.c_[((scaled_matrix[:, :-1].T - scale_min) / scale_max * 10).T, scaled_matrix[:, -1]]
        else:
            MEM_matrix = np.c_[(scaled_matrix[:, :-1] - scale_min) / scale_max * 10, scaled_matrix[:, -1]]

    # Rename rows and columns of all matrices
    def rename_table(x):
        df = pd.DataFrame(x)
        df.columns = marker_names[:-1]
        df.index = pop_names
        return df


    object_list_labeled = [
            rename_table(MAGpop[:, :-1]),
            rename_table(MAGref[:, :-1]),
            rename_table(IQRpop[:, :-1]),
            rename_table(IQRref[:, :-1]),
            rename_table(MEM_matrix[:, :-1]),
            rename_table(prescaled_MEM_matrix[:, :-1])
        ]
    object_list_labeled.append(file_order)

    all_values = {
        "MAGpop": object_list_labeled[0],
        "MAGref": object_list_labeled[1],
        "IQRpop": object_list_labeled[2],
        "IQRref": object_list_labeled[3],
        "MEM_matrix": object_list_labeled[4],
        "prescaled_MEM_matrix": object_list_labeled[5],
        "File Order": object_list_labeled[6]
    }


    # Export pre-scaled MEM values
    if output_prescaled_MEM:
        prescaled_df = pd.DataFrame(prescaled_MEM_matrix, columns=marker_names)
        os.makedirs("output files", exist_ok=True)
        prescaled_df.iloc[:, :-1].to_csv(
            f"./output files/{pd.Timestamp.now().strftime('%Y-%m-%d_%H%M%S')} Pre-scaled MEM matrix.txt",
            sep="\t"
        )

    return all_values


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
    
    file_order = exp_data if file_is_clust else 0

     # Get markers to include in analysis

    markerList = list(range(exp_data.shape[1]))
    exp_data = exp_data.iloc[:, markerList]
    marker_names = list(exp_data.columns)



    # Initialize variables
    marker_names = list(marker_names)
    num_markers = exp_data.shape[1]
    exp_data.columns = marker_names
    num_cells = exp_data.shape[0]
    num_pops = exp_data['cluster'].nunique()
    pop_names = exp_data['cluster'].unique()

    MAGpop = np.zeros((num_pops, num_markers))
    MAGref = np.zeros((num_pops, num_markers))
    IQRpop = np.zeros((num_pops, num_markers))
    IQRref = np.zeros((num_pops, num_markers))



    # Get population medians and IQRs
    for i, pop in enumerate(pop_names):
        subset = exp_data[exp_data['cluster'] == pop]
        MAGpop[i, :] = np.abs(subset.median(skipna=True))
        IQRpop[i, :] = subset.quantile(0.75) - subset.quantile(0.25)

    # Get reference population medians and IQRs
    if choose_ref:
        raise NotImplementedError("choose_ref logic not implemented")
    elif zero_ref:
        MAGref, IQRref = zero_ref_func(exp_data.iloc[:, :-1], num_pops, num_markers)
    else:
        for i, pop in enumerate(pop_names):
            subset = exp_data[exp_data['cluster'] != pop]
            MAGref[i, :] = np.abs(subset.median(skipna=True))
            IQRref[i, :] = subset.quantile(0.75) - subset.quantile(0.25)

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

    for i in range(num_markers - 1):
        IQRpop[:, i] = np.maximum(IQRpop[:, i], IQR_thresh)
        IQRref[:, i] = np.maximum(IQRref[:, i], IQR_thresh)

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
    object_list_labeled.append(file_order)

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


