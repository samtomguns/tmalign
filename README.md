# tmalign

Python bindings for the [TM-align](https://zhanggroup.org/TM-align/) protein structure alignment algorithm. Stripped of all features except for TM-score and Cα RMSD, and written entirely by Claude Code.

## Install

```bash
uv add git+https://github.com/...
# or
pip install .
```

## Usage

```python
import tmalign
import gemmi

st1 = gemmi.read_structure("a.pdb")
st2 = gemmi.read_structure("b.pdb")

result = tmalign.tmalign(st1[0][0], st2[0][0], norm="shorter")
# result: {"tm_score": float, "rmsd": float, "seq_id": float, "n_aligned": int}
```

`norm` options: `"shorter"` (default), `"longer"`, `"average"`, `"chain1"`, `"chain2"`.

Low-level interface (takes Nx3 and Mx3 Cα coordinate arrays directly):

```python
import _tmalign
result = _tmalign.tmalign(xa, ya, seqx, seqy, norm_length)
```

## Credits 
If you use the package for research, you should cite the original TM-score papers:

Y. Zhang, J. Skolnick, Scoring function for automated assessment of protein structure template quality, Proteins, 57: 702-710 (2004).
J. Xu, Y. Zhang, How significant is a protein structure similarity with TM-score=0.5? Bioinformatics, 26, 889-895 (2010).
