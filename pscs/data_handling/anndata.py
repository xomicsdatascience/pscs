# This file contains code for converting user-supplied files (.csv, .tsv) into a
# file containing AnnData objects.

import anndata
import pandas as pd


def load_multitext(quant_data_filepath: str,
                   obs_filepath: str = None,
                   var_filepath: str = None,
                   delimiter: str = ",",
                   first_col_is_index: bool = True) -> anndata.AnnData:
    """
    Loads the specified files and returns the combined data as an AnnData object.
    specified.
    Parameters
    ----------
    quant_data_filepath : str
        Path to the quantitative data
    obs_filepath : str
        Path to the observation-level metadata.
    var_filepath : str
        Path to the variable-level metadata.
    delimiter : str
        Character separator to use for loading the text files. .csv use ",", .tsv use "\t". Default: ","
    first_col_is_index : bool
        Whether the first column in the datafiles is the index.

    Returns
    -------
    anndata.AnnData
        AnnData object with the combined data specified by the input.
    """
    # Read quant data
    data = anndata.read_text(quant_data_filepath, delimiter=delimiter, first_column_names=first_col_is_index)

    if first_col_is_index:
        idx = 0
    else:
        idx = False

    # Read obs data, if any
    if obs_filepath is not None:
        data.obs = pd.read_csv(obs_filepath, delimiter=delimiter, index_col=idx)
    # Read var data, if any
    if var_filepath is not None:
        data.var = pd.read_csv(var_filepath, delimiter=delimiter, index_col=idx)
    return data


