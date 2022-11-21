import unittest
import os
import hashlib
from pscs.metadata.meta_extractor import TableHash, TableDims
from pscs.metadata.metadata import MetaTable, get_metadata

test_topdir = os.path.dirname(__file__)
test_data_dir = os.path.join(test_topdir, 'test_data', 'metadata')

# Test data was generated as follows:
# import random
# import string
# import numpy as np
# import pandas as pd
# columns = [''.join(random.choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for _ in range(20)) for _ in range(50)]
# data_list = []
# num_samples = 100
# for _ in range(num_samples):
#     col_subset = random.sample(columns, k=random.randint(20, 40))
#     data_dict = {}
#     for col in col_subset:
#         data_dict[col] = np.random.random()
#     data_list.append(data_dict)
# df = pd.DataFrame(data_list, columns=columns)
# df.to_csv('sample_csv.csv')
## df.to_csv('sample_tsv.tsv', sep='\t')  # change values above to get different shape for tsv
## Hash values were then obtained using openssl

csv_shape = (100, 50)
tsv_shape = (104, 46)
csv_path = os.path.join(test_data_dir, 'sample_csv.csv')
csv_sha3_256 = open(os.path.join(test_data_dir, 'sample_csv.sha3-256'), 'r').read().strip().split(' ')[-1]
csv_sha3_512 = open(os.path.join(test_data_dir, 'sample_csv.sha3-512'), 'r').read().strip().split(' ')[-1]
tsv_path = os.path.join(test_data_dir, 'sample_tsv.tsv')
tsv_sha3_512 = open(os.path.join(test_data_dir, 'sample_tsv.sha3-512'), 'r').read().strip().split(' ')[-1]
tsv_sha3_256 = open(os.path.join(test_data_dir, 'sample_tsv.sha3-256'), 'r').read().strip().split(' ')[-1]

### Extractor tests
class TestHashExtractor(unittest.TestCase):
    def test_init(self):
        TableHash()
        hash_table = TableHash(name='test_name')
        self.assertEqual('test_name', hash_table.name)
        hash_table = TableHash(hash_func=hashlib.sha3_512)
        self.assertEqual(hash_table.hash.name, hashlib.sha3_512().name)
        hash_table = TableHash(name='test_name', hash_func=hashlib.sha3_512)
        self.assertEqual(hash_table.name, 'test_name')
        self.assertEqual(hash_table.hash.name, hashlib.sha3_512().name)
        return

    def test_file(self):
        # Default hash func test
        hash_table = TableHash(name='test_hash')
        chunk_size = 2**10
        f = open(csv_path, 'rb')
        chunk = f.read(chunk_size)
        while chunk:
            hash_table.process_chunk(chunk=chunk)
            chunk = f.read(chunk_size)
        f.close()
        hash_val = hash_table.finalize()
        self.assertEqual(hash_val, csv_sha3_256)

        # Different hash func test
        hash_table = TableHash(hash_func=hashlib.sha3_512, name='test_hash_sha3512')
        chunk_size = 2**10
        f = open(csv_path, 'rb')
        chunk = f.read(chunk_size)
        while chunk:
            hash_table.process_chunk(chunk=chunk)
            chunk = f.read(chunk_size)
        f.close()
        hash_val = hash_table.finalize()
        self.assertEqual(hash_val, csv_sha3_512)
        return

class TestTableDimsExtractor(unittest.TestCase):
    def test_init(self):
        TableDims()
        dim_table = TableDims(name='test_name')
        self.assertEqual(dim_table.name, 'test_name')
        dim_table = TableDims(sep='\t')
        self.assertEqual(dim_table.sep, '\t')
        dim_table = TableDims(sep=',', name='test_name_dims')
        self.assertEqual(dim_table.sep, ',')
        self.assertEqual(dim_table.name, 'test_name_dims')
        return

    def test_csv(self):
        dim_table = TableDims()
        f = open(csv_path, 'rb')
        chunk_size = 2**10
        chunk = f.read(chunk_size)
        while chunk:
            dim_table.process_chunk(chunk)
            chunk = f.read(chunk_size)
        f.close()
        dim_table_shape = dim_table.finalize()
        self.assertEqual(dim_table_shape, csv_shape)
        return

    def test_tsv(self):
        dim_table = TableDims(sep='\t')
        f = open(tsv_path, 'rb')
        chunk_size = 2**10
        chunk = f.read(chunk_size)
        while chunk:
            dim_table.process_chunk(chunk)
            chunk = f.read(chunk_size)
        f.close()
        dim_table_shape = dim_table.finalize()
        self.assertEqual(dim_table_shape, tsv_shape)
        return

class TestMetaTable(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, MetaTable)
        meta_table = MetaTable(filepath='',
                               sep=',',
                               chunk_size=2**10,
                               meta_extractors=[TableDims()])
        return

    def test_dims_and_hash(self):
        meta_table = MetaTable(filepath=csv_path,
                               sep=',',
                               meta_extractors=[TableDims(sep=',', name='dims'),
                                                TableHash(hash_func=hashlib.sha3_256, name='hash')])
        meta_table.process()
        self.assertEqual(meta_table.meta_dict['dims'], csv_shape)
        self.assertEqual(meta_table.meta_dict['hash'], csv_sha3_256)

        # Test different parameters
        meta_table = MetaTable(filepath=tsv_path,
                               sep='\t',
                               meta_extractors=[TableDims(sep='\t', name='other_dims'),
                                                TableHash(hash_func=hashlib.sha3_512, name='other_hash')])
        meta_table.process()
        self.assertEqual(meta_table.meta_dict['other_dims'], tsv_shape)
        self.assertEqual(meta_table.meta_dict['other_hash'], tsv_sha3_512)
        return

class TestGetMeta(unittest.TestCase):
    def test_csv(self):
        # Note: this test will fail if the default parameter values of the extractors are changed.
        meta_dict = get_metadata(csv_path)
        self.assertEqual(meta_dict['table_dimensions'], csv_shape)
        self.assertEqual(meta_dict['table_hash'], csv_sha3_256)
        return

    def test_tsv(self):
        # Note: this test will fail if the default parameter values of the extractors are changed.
        meta_dict = get_metadata(tsv_path)
        self.assertEqual(meta_dict['table_dimensions'], tsv_shape)
        self.assertEqual(meta_dict['table_hash'], tsv_sha3_256)
        return
