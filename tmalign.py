import numpy as np
import gemmi
import _tmalign


def tmalign(
    chain1: gemmi.Chain,
    chain2: gemmi.Chain,
    norm: str = "shorter",
) -> dict:
    """
    Align two protein chains and return TM-score and related metrics.

    Parameters
    ----------
    chain1 : gemmi.Chain
    chain2 : gemmi.Chain
    norm : str
        How to normalise the TM-score. One of:
          "shorter"  – normalise by the length of the shorter chain (default)
          "longer"   – normalise by the length of the longer chain
          "average"  – normalise by the average length of the chains
          "chain1"   – normalise by the length of chain1
          "chain2"   – normalise by the length of chain2

    Returns
    -------
    dict with keys:
        tm_score  : TM-score under the requested normalisation
        rmsd      : RMSD over aligned residue pairs (Å)
        seq_id    : sequence identity fraction over aligned pairs
        n_aligned : number of aligned residue pairs within d8 cutoff
    """
    VALID_NORMS = ("shorter", "longer", "average", "chain1", "chain2")
    if norm not in VALID_NORMS:
        raise ValueError(f"norm must be one of {VALID_NORMS}, got {norm!r}")

    xa, seqx = _extract_ca(chain1)
    ya, seqy = _extract_ca(chain2)

    xlen, ylen = len(seqx), len(seqy)
    norm_length = {
        "shorter": min(xlen, ylen),
        "longer":  max(xlen, ylen),
        "average": 0.5 * (xlen + ylen),
        "chain1":  xlen,
        "chain2":  ylen,
    }[norm]

    raw = _tmalign.tmalign(xa, ya, seqx, seqy, norm_length)  # type: ignore[attr-defined]

    return {
        "tm_score":  raw["tm_score"],
        "rmsd":      raw["rmsd"],
        "seq_id":    raw["seq_id"],
        "n_aligned": raw["n_aligned"],
    }


def _extract_ca(chain: gemmi.Chain) -> tuple[np.ndarray, str]:
    """
    Extract Cα coordinates and single-letter sequence from a gemmi Chain.
    Uses first conformer only; skips residues with no CA atom.

    Returns
    -------
    coords : np.ndarray, shape (N, 3), dtype float64
    seq    : str, length N
    """
    coords, seq = [], []
    for residue in chain.first_conformer():
        ca = residue.find_atom("CA", "*")
        if ca is None:
            continue
        coords.append([ca.pos.x, ca.pos.y, ca.pos.z])
        seq.append(gemmi.one_letter_code([residue.name]) or "X")
    return np.array(coords, dtype=np.float64), "".join(seq)
