from abc import ABC, abstractmethod
import hashlib

str_encoding = 'utf-8'
# NOTE: The principle here is to operate on chunks of data. This avoids having to load the same file multiple times.

class MetaExtractor(ABC):
    @abstractmethod
    def __init__(self,
                 name: str = 'extractor_name'):
        raise NotImplementedError('__init__ method not implemented')
        return

    @abstractmethod
    def process_chunk(self,
                      chunk: bytes):
        """
        Process a chunk of binary data of arbitrary length.
        Parameters
        ----------
        chunk : bytes
            Chunk of data to process, as obtained from an IO stream (e.g. open('example', 'rb').read(1024))

        Returns
        -------
        None
        """
        raise NotImplementedError('process_chunk method not implemented')
        return

    @abstractmethod
    def finalize(self):
        """
        Signals the end of data. Meta data should be returned here.
        Returns
        -------
        object
            Data structure containing the metadata (int, float, str, tuple, etc.)
        """
        raise NotImplementedError('finalize method not implemented')
        return


class TableDims(MetaExtractor):
    """Object to examine chunks of a table and determine overall size (number of rows, columns)."""
    def __init__(self, sep: str = ',', name='table_dimensions'):
        """
        Initializes object.
        Parameters
        ----------
        sep : str
            Separator used for separating columns.
        name : str
            Name to associate with this instance.
        """
        self.name = name
        self.num_cols = -1  # number of columns in the table
        self.num_rows = -1  # number of rows in the table
        self.sep = sep
        self._num_seps = 0  # internal; number of separator characters in the first row
        self._num_newline = 0  # internal; number of newlines in the file
        return

    def process_chunk(self,
                      chunk: bytes):
        """
        Process a data chunk and store relevant data.
        Parameters
        ----------
        chunk : bytes
            Bytes from the file to examine.
        """
        str_chunk = chunk.decode(str_encoding)
        self._num_newline += str_chunk.count('\n')
        if self.num_cols == -1:  # First row is still being processed
            if '\n' not in str_chunk:
                self._num_seps += str_chunk.count(self.sep)
            else:  # First row ends in this chunk
                newline_loc = str_chunk.find('\n')
                self._num_seps += str_chunk[:newline_loc].count(self.sep)  # number of seps before newline
                self.num_cols = self._num_seps + 1
        return

    def finalize(self) -> tuple:
        """
        Signal end of chunks; compute relevant values and return
        Returns
        -------
        tuple
            Tuple containing (num_rows, num_cols). num_rows excludes the header
        """
        self.num_rows = self._num_newline - 1  # exclude header; num_cols was already computed
        return self.num_rows, self.num_cols


class TableHash(MetaExtractor):
    """Wrapper for the SHA3 256 object from hashlib"""
    def __init__(self,
                 hash_func: callable = hashlib.sha3_256,
                 name: str = 'table_hash'):
        """
        hash_func : callable
            Hash function to use; expects structure similar to hashlib hash functions.
        name : str
            Name to associate with this instance.
        Parameters
        ----------
        hash_func
        name
        """
        self.hash = hash_func()
        self.name = name
        return

    def process_chunk(self,
                      chunk: bytes):
        """
        Digest the chunk into the running hash.
        Parameters
        ----------
        chunk : bytes
            Chunk of data from the file to process.
        """
        self.hash.update(chunk)
        return

    def finalize(self):
        return self.hash.hexdigest()
