#!/usr/bin/env python3
"""
benchmark_harness.py

Benchmark harness for the Phase 1 + Phase 2 refactor:
  - Reader throughput in dict mode vs array mode
  - Core per-qid accumulation kernel in:
      (a) pure Python (array-based)
      (b) Numba JIT (array-based), if numba is installed

This script is designed to be practical on HPC:
  - Lets you cap work (max_qids, max_pairs_per_qid, max_chunks)
  - Supports warmup + repeats
  - Prints throughput (pairs/sec) and wall time

Usage examples
--------------
# Benchmark reading/parsing only (dict vs array) on real file:
python benchmark_harness.py --dist-file distances.tsv --mode reader

# Benchmark kernel on synthetic data (default sizes are safe):
python benchmark_harness.py --mode synthetic

# Benchmark kernel + reader arrays on real file (if you have updated reader):
python benchmark_harness.py --dist-file distances.tsv --mode real-arrays

Notes
-----
- This harness does NOT require your full assign() class. It focuses on the
  hot-loop primitives that dominate runtime (qid×rid accumulation).
- For 200 qids × 100,000 rids = 20M pairs, you may need to reduce repeats or
  limit chunks depending on node memory/time policies.
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np


try:
    import numba as nb  # type: ignore

    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False


# ---------------------------- Data structures --------------------------------


@dataclass(frozen=True)
class DistChunkArrays:
    """
    Array representation of a distance chunk.

    qids:
        int64 array of query IDs (encoded)
    rid_ids:
        object array of arrays, each element is int64 rid indices for that qid
    dists:
        object array of arrays, each element is float32/float64 distances for that qid
    rid_cluster:
        object array of arrays, each element is int32 cluster IDs (finest-rank) per rid
        aligned with rid_ids and dists
    """
    qids: np.ndarray
    rid_ids: np.ndarray
    dists: np.ndarray
    rid_cluster: np.ndarray


# ---------------------------- Synthetic generator ----------------------------


def generate_synthetic_chunk(
    n_qids: int,
    n_pairs_per_qid: int,
    n_clusters: int,
    seed: int = 0,
    dist_dtype: np.dtype = np.float32,
) -> DistChunkArrays:
    """
    Generate synthetic per-qid arrays approximating your scale.

    Parameters
    ----------
    n_qids
        Number of queries in chunk.
    n_pairs_per_qid
        Number of rid distances per query.
    n_clusters
        Number of cluster IDs (finest rank) that rid samples map to.
        Higher means more fragmented bucket maps.
    seed
        RNG seed.
    dist_dtype
        dtype for distances.

    Returns
    -------
    DistChunkArrays
        A synthetic chunk.
    """
    rng = np.random.default_rng(seed)

    qids = np.arange(n_qids, dtype=np.int64)

    rid_ids = np.empty(n_qids, dtype=object)
    dists = np.empty(n_qids, dtype=object)
    rid_cluster = np.empty(n_qids, dtype=object)

    # Synthetic rids: we only need stable integer IDs for benchmarking
    for i in range(n_qids):
        rid_ids[i] = rng.integers(0, 10_000_000, size=n_pairs_per_qid, dtype=np.int64)
        # Distances: roughly small integers / floats, typical for allele Hamming
        d = rng.integers(0, 50, size=n_pairs_per_qid).astype(dist_dtype)
        dists[i] = d
        rid_cluster[i] = rng.integers(0, n_clusters, size=n_pairs_per_qid, dtype=np.int32)

    return DistChunkArrays(qids=qids, rid_ids=rid_ids, dists=dists, rid_cluster=rid_cluster)


# ---------------------------- Core accumulation ------------------------------


def accumulate_stats_python(
    cluster_ids: np.ndarray,
    dists: np.ndarray,
    thresholds: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Pure Python / NumPy accumulation (array-based) for ONE qid.

    Computes per-cluster:
      - count
      - sum
      - min
      - max
      - matched_counts per rank (n_ranks)

    Implementation strategy:
      - Remap sparse cluster_ids to dense [0..k-1] for this qid
      - Use bincount for count/sum and per-rank matched_counts
      - Use np.minimum.at / np.maximum.at for min/max (requires dense indices)

    Parameters
    ----------
    cluster_ids
        int32 array of cluster IDs (sparse/global).
    dists
        float array of distances.
    thresholds
        float array length n_ranks.

    Returns
    -------
    dense_cluster_ids
        int32 array length k of original cluster IDs for each dense index.
    count
        int32 array length k.
    sum_dist
        float64 array length k.
    min_dist
        float array length k.
    max_dist
        float array length k.
    matched
        int32 array shape (n_ranks, k).
    """
    # Remap to dense indices per qid
    uniq, inv = np.unique(cluster_ids, return_inverse=True)
    k = uniq.size

    count = np.bincount(inv, minlength=k).astype(np.int32)
    sum_dist = np.bincount(inv, weights=dists.astype(np.float64), minlength=k)

    # min/max via "at" reductions
    min_dist = np.full(k, np.inf, dtype=dists.dtype)
    max_dist = np.full(k, -np.inf, dtype=dists.dtype)
    np.minimum.at(min_dist, inv, dists)
    np.maximum.at(max_dist, inv, dists)

    # matched counts per rank
    n_ranks = thresholds.size
    matched = np.zeros((n_ranks, k), dtype=np.int32)
    for r in range(n_ranks):
        mask = (dists <= thresholds[r]).astype(np.int32)
        matched[r, :] = np.bincount(inv, weights=mask.astype(np.int64), minlength=k).astype(np.int32)

    return uniq.astype(np.int32), count, sum_dist, min_dist, max_dist, matched


if NUMBA_AVAILABLE:
    @nb.njit(cache=True)
    def _reduce_sorted_numba(
        cluster_ids_sorted: np.ndarray,
        dists_sorted: np.ndarray,
        thresholds: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Numba kernel: reduce over sorted cluster IDs for ONE qid.

        Parameters
        ----------
        cluster_ids_sorted
            int32 array sorted by cluster id.
        dists_sorted
            float array aligned with cluster_ids_sorted.
        thresholds
            float array length n_ranks.

        Returns
        -------
        out_cluster_ids
            int32 array of unique cluster ids.
        count
            int32 array.
        sum_dist
            float64 array.
        min_dist
            float array.
        max_dist
            float array.
        matched
            int32 array shape (n_ranks, n_clusters_in_qid).
        """
        n = cluster_ids_sorted.shape[0]
        n_ranks = thresholds.shape[0]

        # First pass: count unique clusters
        n_unique = 0
        last = -2147483648
        for i in range(n):
            cid = cluster_ids_sorted[i]
            if cid != last:
                n_unique += 1
                last = cid

        out_cluster_ids = np.empty(n_unique, dtype=np.int32)
        count = np.zeros(n_unique, dtype=np.int32)
        sum_dist = np.zeros(n_unique, dtype=np.float64)
        min_dist = np.empty(n_unique, dtype=dists_sorted.dtype)
        max_dist = np.empty(n_unique, dtype=dists_sorted.dtype)
        matched = np.zeros((n_ranks, n_unique), dtype=np.int32)

        # Second pass: reduce
        idx = -1
        last = -2147483648
        for i in range(n):
            cid = cluster_ids_sorted[i]
            d = dists_sorted[i]
            if cid != last:
                idx += 1
                out_cluster_ids[idx] = cid
                count[idx] = 0
                sum_dist[idx] = 0.0
                min_dist[idx] = d
                max_dist[idx] = d
                last = cid

            count[idx] += 1
            sum_dist[idx] += float(d)
            if d < min_dist[idx]:
                min_dist[idx] = d
            if d > max_dist[idx]:
                max_dist[idx] = d

            for r in range(n_ranks):
                if d <= thresholds[r]:
                    matched[r, idx] += 1

        return out_cluster_ids, count, sum_dist, min_dist, max_dist, matched


def accumulate_stats_numba(
    cluster_ids: np.ndarray,
    dists: np.ndarray,
    thresholds: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Wrapper for Numba kernel. Sorts by cluster id, then reduces.

    Parameters
    ----------
    cluster_ids
        int32 array of cluster IDs (global).
    dists
        float array.
    thresholds
        float array length n_ranks.

    Returns
    -------
    Same as kernel.
    """
    order = np.argsort(cluster_ids, kind="mergesort")
    return _reduce_sorted_numba(cluster_ids[order], dists[order], thresholds)


# ---------------------------- Reader benchmarks ------------------------------


def benchmark_reader(
    dist_file: str,
    max_chunks: int,
    n_records: int,
) -> None:
    """
    Benchmark reader throughput (dict mode vs array mode), if available.

    This expects your updated reader to live at:
      genomic_address_service.classes.reader.dist_reader

    It times:
      - iterating chunks in dict mode (baseline)
      - iterating chunks in array mode (id_mapper required)
    """
    try:
        from genomic_address_service.classes.reader import dist_reader  # type: ignore
    except Exception as e:
        print(f"[ERROR] Cannot import dist_reader: {e}", file=sys.stderr)
        sys.exit(2)

    # Dummy ID mapper: maps string IDs to stable ints on the fly
    # For benchmarking, we keep this minimal.
    id_map: Dict[str, int] = {}
    next_id = 0

    def id_mapper(x: str) -> int:
        nonlocal next_id
        if x in id_map:
            return id_map[x]
        id_map[x] = next_id
        next_id += 1
        return id_map[x]

    # Dict mode
    t0 = time.perf_counter()
    pairs = 0
    chunks = 0
    for chunk in dist_reader(f=dist_file, n_records=n_records, delim="\t").read_data():
        chunks += 1
        for qid, rd in chunk.items():
            pairs += len(rd)
        if chunks >= max_chunks:
            break
    t1 = time.perf_counter()
    print(f"[reader:dict] chunks={chunks} pairs={pairs:,} time={t1-t0:.3f}s rate={pairs/(t1-t0):,.0f} pairs/s")

    # Array mode (if implemented)
    t0 = time.perf_counter()
    pairs = 0
    chunks = 0
    for chunk in dist_reader(f=dist_file, n_records=n_records, delim="\t").read_data(id_mapper=id_mapper):
        chunks += 1
        # Updated reader should yield a DistChunk-like object or tuple; handle both
        if hasattr(chunk, "qids") and hasattr(chunk, "dists"):
            # chunk.dists is object array of per-qid arrays
            for arr in chunk.dists:
                pairs += arr.shape[0]
        else:
            # Unknown shape; can't count, but still time
            pass
        if chunks >= max_chunks:
            break
    t1 = time.perf_counter()
    print(f"[reader:array] chunks={chunks} pairs={pairs:,} time={t1-t0:.3f}s rate={pairs/(t1-t0):,.0f} pairs/s")


# ---------------------------- Kernel benchmarks ------------------------------


def benchmark_kernel_synthetic(
    n_qids: int,
    n_pairs_per_qid: int,
    n_clusters: int,
    n_ranks: int,
    repeats: int,
    warmup: int,
    use_numba: bool,
) -> None:
    """
    Benchmark core accumulation kernel on synthetic arrays.

    Prints:
      - wall time
      - pairs/sec throughput
      - rough per-qid cost
    """
    thresholds = np.linspace(5, 20, num=n_ranks).astype(np.float32)

    # Prepare synthetic chunk once (so we benchmark compute, not RNG)
    chunk = generate_synthetic_chunk(
        n_qids=n_qids,
        n_pairs_per_qid=n_pairs_per_qid,
        n_clusters=n_clusters,
        seed=0,
        dist_dtype=np.float32,
    )

    total_pairs = n_qids * n_pairs_per_qid

    def run_python() -> float:
        t0 = time.perf_counter()
        for i in range(n_qids):
            accumulate_stats_python(
                chunk.rid_cluster[i].astype(np.int32),
                chunk.dists[i].astype(np.float32),
                thresholds,
            )
        return time.perf_counter() - t0

    def run_numba() -> float:
        if not NUMBA_AVAILABLE:
            raise RuntimeError("Numba not available")
        t0 = time.perf_counter()
        for i in range(n_qids):
            accumulate_stats_numba(
                chunk.rid_cluster[i].astype(np.int32),
                chunk.dists[i].astype(np.float32),
                thresholds,
            )
        return time.perf_counter() - t0

    # Warmup (important for numba compile and cache)
    for _ in range(warmup):
        _ = run_python()
        if use_numba and NUMBA_AVAILABLE:
            _ = run_numba()

    # Timed repeats
    best_py = float("inf")
    for _ in range(repeats):
        dt = run_python()
        best_py = min(best_py, dt)

    print(
        f"[kernel:python] qids={n_qids} pairs/qid={n_pairs_per_qid} "
        f"pairs={total_pairs:,} best={best_py:.3f}s rate={total_pairs/best_py:,.0f} pairs/s"
    )

    if use_numba:
        if not NUMBA_AVAILABLE:
            print("[kernel:numba] numba not installed; skipping")
            return

        best_nb = float("inf")
        for _ in range(repeats):
            dt = run_numba()
            best_nb = min(best_nb, dt)

        print(
            f"[kernel:numba]  qids={n_qids} pairs/qid={n_pairs_per_qid} "
            f"pairs={total_pairs:,} best={best_nb:.3f}s rate={total_pairs/best_nb:,.0f} pairs/s "
            f"speedup={best_py/best_nb:.2f}x"
        )


# ---------------------------- Main -------------------------------------------


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(description="Benchmark harness for array+Numba refactor.")
    p.add_argument(
        "--mode",
        choices=["synthetic", "reader"],
        default="synthetic",
        help="Benchmark mode: synthetic kernel or reader parsing.",
    )
    p.add_argument("--dist-file", default="", help="Distance file for reader mode.")
    p.add_argument("--max-chunks", type=int, default=3, help="Max chunks to read in reader mode.")
    p.add_argument("--n-records", type=int, default=1000, help="n_records passed to dist_reader.")

    # Synthetic sizing (safe defaults)
    p.add_argument("--n-qids", type=int, default=200, help="Synthetic number of qids.")
    p.add_argument("--pairs-per-qid", type=int, default=50_000, help="Synthetic pairs per qid.")
    p.add_argument("--n-clusters", type=int, default=50_000, help="Synthetic number of clusters.")
    p.add_argument("--n-ranks", type=int, default=4, help="Number of ranks/thresholds.")
    p.add_argument("--repeats", type=int, default=3, help="Timed repeats (best-of).")
    p.add_argument("--warmup", type=int, default=1, help="Warmup runs.")
    p.add_argument("--no-numba", action="store_true", help="Disable numba benchmark.")
    return p.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()

    if args.mode == "reader":
        if not args.dist_file:
            print("[ERROR] --dist-file is required for reader mode", file=sys.stderr)
            sys.exit(2)
        benchmark_reader(dist_file=args.dist_file, max_chunks=args.max_chunks, n_records=args.n_records)
        return

    # synthetic kernel mode
    use_numba = not args.no_numba
    benchmark_kernel_synthetic(
        n_qids=args.n_qids,
        n_pairs_per_qid=args.pairs_per_qid,
        n_clusters=args.n_clusters,
        n_ranks=args.n_ranks,
        repeats=args.repeats,
        warmup=args.warmup,
        use_numba=use_numba,
    )


if __name__ == "__main__":
    main()
