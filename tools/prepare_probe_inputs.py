"""Reassemble lossless layer shards into the original probe input arrays."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifest = json.loads((ROOT / "data/esm_probe/embedding_manifest.json").read_text())
    pieces = {key: [] for key in manifest["arrays"]}
    cursor = 0
    for row in manifest["shards"]:
        path = ROOT / row["path"]
        assert row["start_layer"] == cursor
        with path.open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == row["sha256"]
        with np.load(path, allow_pickle=False) as z:
            for key in pieces:
                pieces[key].append(z[key])
        cursor = row["end_layer_exclusive"]
    arrays = {key: np.concatenate(values, axis=1) for key, values in pieces.items()}
    for key, arr in arrays.items():
        expected = manifest["arrays"][key]
        assert list(arr.shape) == expected["shape"] and str(arr.dtype) == expected["dtype"]
        assert hashlib.sha256(arr.tobytes()).hexdigest() == expected["sha256_c_order"]
    output = ROOT / "outputs/esm-probe/esm_probe_embeddings.npz"
    if output.exists():
        with np.load(output, allow_pickle=False) as z:
            assert set(z.files) == set(arrays)
            assert all(np.array_equal(z[k], arrays[k]) for k in arrays), "Existing probe input differs"
    else:
        np.savez_compressed(output, **arrays)
    with output.open("rb") as stream:
        exact = hashlib.file_digest(stream, "sha256").hexdigest() == manifest["original_sha256"]
    print(json.dumps({"arrays_exact": True, "original_npz_bytes_exact": exact, "layers": cursor}))


if __name__ == "__main__":
    main()

