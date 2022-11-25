import numpy as np

GT_TO_BED_BITS = {
    0: 3,
    1: 2,
    np.nan: 1,
    2: 0,
}


def rand_bed_file(n: int, p: int, filepath: str, magic_numbers=True):
    """
    Generates a random genotype matrix of dimensions n x p
    where n: number of individuals, p: number of markers
    and stores it directly on disk in .bed format

    This uses less memory than generating a full genotype matrix
    in memory and writing it to disk in one go, but it might take
    more time.
    """
    barr = bytearray()
    # add magic numbers for the plink col-major format
    # (sample-major has a 0 as the third byte)
    if magic_numbers:
        barr.extend([108, 27, 1])

    rng = np.random.default_rng()

    with open(filepath, "wb") as out:
        for _col_ix in range(p):
            col = _rand_gt_col(n, rng)
            _add_col_to_bytearray_with_padding(
                col, barr)
            out.write(barr)
            barr = bytearray()


def rand_gt_mat(n: int, p: int) -> np.ndarray:
    """
    Generates a random genotype matrix of dimensions n x p
    where n: number of individuals, p: number of markers
    """
    rng = np.random.default_rng()
    maf = rng.uniform(0.0, 0.5, p)
    return np.column_stack([rng.binomial(2, f, n).astype(np.byte) for f in maf])


def write_m_to_bed(m: np.ndarray, filepath: str, magic_numbers=True):
    """
    Writes a genotype matrix stored in a np.ndarray to disk in .bed format
    """
    barr = bytearray()
    # add magic numbers for the plink col-major format
    # (sample-major has a 0 as the third byte)
    if magic_numbers:
        barr.extend([108, 27, 1])

    with open(filepath, "wb") as out:
        for col_ix in range(m.shape[1]):
            _add_mat_col_to_bytearray_with_padding(
                m, col_ix, barr)
            out.write(barr)
            barr = bytearray()


def _rand_gt_col(n: int, rng) -> np.ndarray:
    maf = rng.uniform(0.0, 0.5)
    return rng.binomial(2, maf, n).astype(np.byte)


def _to_bytes(m: np.ndarray):
    barr = bytearray()
    for col_ix in range(m.shape[1]):
        _add_mat_col_to_bytearray_with_padding(m, col_ix, barr)
    return barr


def _add_mat_col_to_bytearray(m: np.ndarray, col_ix: int, barr: bytearray):
    assert m.shape[0] % 4 == 0
    row_ix = 0
    while row_ix < m.shape[0]:
        b = 0
        row_ix += 4
        for _offset in range(4):
            row_ix -= 1
            gt = m[row_ix, col_ix]
            b = b << 2
            if np.isnan(gt):
                b += 1
            else:
                b += (GT_TO_BED_BITS[gt])
        row_ix += 4
        barr.append(b)


def _add_mat_col_to_bytearray_with_padding(m: np.ndarray, col_ix: int, barr: bytearray):
    row_ix = 0
    while row_ix < m.shape[0]:
        b = 0
        row_ix += 4
        for _offset in range(4):
            row_ix -= 1
            gt = 0 if row_ix >= m.shape[0] else m[row_ix, col_ix]
            b = b << 2
            if np.isnan(gt):
                b += 1
            else:
                b += (GT_TO_BED_BITS[gt])
        row_ix += 4
        barr.append(b)


def _add_col_to_bytearray_with_padding(col: np.ndarray, barr: bytearray):
    row_ix = 0
    n = col.shape[0]
    while row_ix < n:
        b = 0
        row_ix += 4
        for _offset in range(4):
            row_ix -= 1
            gt = 0 if row_ix >= n else col[row_ix]
            b = b << 2
            if np.isnan(gt):
                b += 1
            else:
                b += (GT_TO_BED_BITS[gt])
        row_ix += 4
        barr.append(b)
