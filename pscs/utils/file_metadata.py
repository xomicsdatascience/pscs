import os
from collections import defaultdict as dd
from pathlib import Path
from typing import Union

image_exts = {"apng", "jpg", "jpeg", "png", "avif", "gif", "svg", "webp", "bmp", "ico", "tiff"}
text_exts = {"txt", "doc", "docx"}
table_exts = {"csv", "tsv"}
data_exts = {"h5ad"}

ext_mapper = dd(lambda: "unknown")
for ext in image_exts:
    ext_mapper[ext] = "image"
for ext in text_exts:
    ext_mapper[ext] = "text"
for ext in table_exts:
    ext_mapper[ext] = "table"
for ext in data_exts:
    ext_mapper[ext] = "data"

def get_file_type(filename: Union[Path, str]) -> str:
    """Examines the extension of the given filename and returns the type of file (e.g., "image", "table", etc.). If
    the type is not recognized, the function returns "unknown"."""
    ext = Path(filename).suffix[1:]
    return ext_mapper[ext]
