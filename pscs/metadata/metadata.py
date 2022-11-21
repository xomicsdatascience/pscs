from .meta_extractor import TableDims, TableHash


def get_metadata(filepath: str) -> dict:
    """
    Examines the file to extract metadata (sample count, gene count, etc.).
    Parameters
    ----------
    filepath : str
        Path to the file to load

    Returns
    -------
    dict
        Dict containing meta data

    Raises
    ------
    ValueError
        If filepath does not correspond to a .csv or .tsv file.
    """
    if filepath.endswith('.csv'):
        sep = ','
    elif filepath.endswith('.tsv'):
        sep = '\t'
    else:
        raise ValueError('Input filetype not recognized')

    meta_data = MetaTable(filepath=filepath, sep=sep)
    meta_data.process()
    return meta_data.meta_dict


class MetaTable:
    """
    Wrapper for the different types of metadata extractors operating on tables.
    """
    def __init__(self,
                 filepath: str,
                 sep: str = ',',
                 chunk_size: int = 2**12,
                 meta_extractors: list = None):
        """
        Initializes the MetaTable and stores relevant parameters.
        Parameters
        ----------
        filepath : str
            Path to the file to process.
        sep : str
            Separator used to denote fields in the table.
        chunk_size : int
            Number of bytes to load per chunk.
        meta_extractors : list
            List of meta data extractors to use.
        """
        self.filepath = filepath
        self.meta_dict = {}  # used to store meta data using extractor names
        self.chunk_size = chunk_size
        if meta_extractors is None:
            self.meta_extractors = [TableDims(sep=sep), TableHash()]
        else:
            self.meta_extractors = meta_extractors
        return

    def process(self):
        """
        Process the data and return the metadata.
        """
        f = open(self.filepath, 'rb')
        chunk = f.read(self.chunk_size)
        while chunk:
            for extractor in self.meta_extractors:
                extractor.process_chunk(chunk)
            chunk = f.read(self.chunk_size)
        f.close()
        for extractor in self.meta_extractors:
            self.meta_dict[extractor.name] = extractor.finalize()
        return