from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from genomic_address_service.constants import EXTENSIONS
from genomic_address_service.utils import is_file_ok
from genomic_address_service.classes.reader import DistChunk, dist_reader  # uses the updated reader.py


try:
    from numba import njit
    _HAVE_NUMBA = True
except Exception:
    _HAVE_NUMBA = False


# ----------------------------- Numba hot loop --------------------------------

if _HAVE_NUMBA:
    @njit(cache=True)
    def reduce_clusters_sorted(
        cluster_ids_sorted: np.ndarray,
        dists_sorted: np.ndarray,
        thresholds: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
        """
        Reduce per-cluster stats for sorted cluster ids.

        Parameters
        ----------
        cluster_ids_sorted
            int64 cluster ids, sorted non-decreasing.
        dists_sorted
            float32 distances aligned to cluster_ids_sorted.
        thresholds
            float32 thresholds of length n_ranks.

        Returns
        -------
        unique_ids
            Unique cluster ids in order of appearance.
        counts
            int32 counts per cluster.
        sums
            float64 sums per cluster.
        mins
            float32 min per cluster.
        maxs
            float32 max per cluster.
        matched
            int32 matched counts per cluster per rank, shape (n_clusters, n_ranks).
        min_dist
            float32 minimum distance across all inputs.
        """
        n = cluster_ids_sorted.shape[0]
        n_ranks = thresholds.shape[0]
        min_dist = np.float32(np.inf)
        if n == 0:
            return (
                np.empty(0, np.int64),
                np.empty(0, np.int32),
                np.empty(0, np.float64),
                np.empty(0, np.float32),
                np.empty(0, np.float32),
                np.empty((0, n_ranks), np.int32),
                min_dist,
            )

        # First pass to count clusters (runs)
        runs = 1
        for i in range(1, n):
            if cluster_ids_sorted[i] != cluster_ids_sorted[i - 1]:
                runs += 1

        unique_ids = np.empty(runs, np.int64)
        counts = np.zeros(runs, np.int32)
        sums = np.zeros(runs, np.float64)
        mins = np.full(runs, np.float32(np.inf), np.float32)
        maxs = np.full(runs, np.float32(-np.inf), np.float32)
        matched = np.zeros((runs, n_ranks), np.int32)

        run_idx = 0
        cur_id = cluster_ids_sorted[0]
        unique_ids[0] = cur_id

        for i in range(n):
            cid = cluster_ids_sorted[i]
            d = dists_sorted[i]
            if d < min_dist:
                min_dist = d

            if cid != cur_id:
                run_idx += 1
                cur_id = cid
                unique_ids[run_idx] = cid

            counts[run_idx] += 1
            sums[run_idx] += d
            if d < mins[run_idx]:
                mins[run_idx] = d
            if d > maxs[run_idx]:
                maxs[run_idx] = d

            # n_ranks is small (e.g., 4) so loop is fine
            for r in range(n_ranks):
                if d <= thresholds[r]:
                    matched[run_idx, r] += 1

        return unique_ids, counts, sums, mins, maxs, matched, min_dist


def _reduce_clusters_sorted_fallback(
    cluster_ids_sorted: np.ndarray,
    dists_sorted: np.ndarray,
    thresholds: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Non-numba fallback reducer. Same outputs as reduce_clusters_sorted().
    """
    n = cluster_ids_sorted.shape[0]
    n_ranks = int(thresholds.shape[0])
    if n == 0:
        return (
            np.empty(0, np.int64),
            np.empty(0, np.int32),
            np.empty(0, np.float64),
            np.empty(0, np.float32),
            np.empty(0, np.float32),
            np.empty((0, n_ranks), np.int32),
            float("inf"),
        )

    # Find run boundaries
    change = np.empty(n, dtype=bool)
    change[0] = True
    change[1:] = cluster_ids_sorted[1:] != cluster_ids_sorted[:-1]
    run_starts = np.where(change)[0]
    run_ends = np.append(run_starts[1:], n)

    unique_ids = cluster_ids_sorted[run_starts].astype(np.int64, copy=False)
    runs = unique_ids.shape[0]

    counts = (run_ends - run_starts).astype(np.int32, copy=False)
    sums = np.empty(runs, dtype=np.float64)
    mins = np.empty(runs, dtype=np.float32)
    maxs = np.empty(runs, dtype=np.float32)
    matched = np.zeros((runs, n_ranks), dtype=np.int32)

    min_dist = float(np.min(dists_sorted))

    for i, (s, e) in enumerate(zip(run_starts, run_ends)):
        seg = dists_sorted[s:e]
        sums[i] = float(np.sum(seg, dtype=np.float64))
        mins[i] = float(np.min(seg))
        maxs[i] = float(np.max(seg))
        for r in range(n_ranks):
            matched[i, r] = int(np.sum(seg <= thresholds[r]))

    return unique_ids, counts, sums, mins, maxs, matched, min_dist


# --------------------------- assignment structures ----------------------------

@dataclass
class ClusterStatsArray:
    """
    Array-backed per-cluster stats for a single rank bucket.

    Attributes
    ----------
    ids
        Cluster IDs (int64).
    counts
        Cluster counts (int32).
    sums
        Sum of distances (float64).
    mins
        Min distance (float32).
    maxs
        Max distance (float32).
    matched
        Matched counts per rank (int32), shape (n_clusters, n_ranks).
    """
    ids: np.ndarray
    counts: np.ndarray
    sums: np.ndarray
    mins: np.ndarray
    maxs: np.ndarray
    matched: np.ndarray

    @property
    def means(self) -> np.ndarray:
        """Mean distance per cluster (float64)."""
        return self.sums / np.maximum(self.counts.astype(np.float64), 1.0)


class UnionFind:
    """Disjoint-set union-find for connected components."""

    def __init__(self, items: Iterable[str]) -> None:
        self.parent: Dict[str, str] = {}
        self.rank: Dict[str, int] = {}
        for x in items:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x: str) -> str:
        """Find representative."""
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        """Union sets."""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1

    def components(self) -> List[List[str]]:
        """Return components as lists."""
        groups: Dict[str, List[str]] = {}
        for x in self.parent:
            r = self.find(x)
            groups.setdefault(r, []).append(x)
        return list(groups.values())


class assign:
    """
    Two-pass incremental assignment with Numba-accelerated pass-1 bucketing.

    Pass 1 (attach-to-existing):
      - Uses only pre-existing references (snapshot at chunk start).
      - Computes cluster stats at finest rank via array+Numba reduction.

    Pass 2 (de novo clustering):
      - Clusters remaining qids among themselves, matching requested linkage method,
        except: majority uses average linkage.
    """

    AVAILABLE_METHODS = ["average", "complete", "single", "majority"]

    ERROR_MISSING_DELIMITER = "delimiter was not found"
    ERROR_LENGTH = "genomic address length is incorrect"
    ERROR_NON_INTEGER = "address could not be converted to an integer"

    def __init__(
        self,
        dist_file: str,
        membership_file: str,
        threshold_map: Mapping[str, float],
        linkage_method: str,
        address_col: str,
        sample_col: str,
        batch_size: int,
        delimiter: str,
        marjority_fraction: float = 0.6,
    ) -> None:
        self.dist_file = dist_file
        self.batch_size = batch_size
        self.status = True
        self.threshold_map = dict(threshold_map)
        self.rank_ids: List[str] = list(self.threshold_map.keys())
        self.thresholds: np.ndarray = np.array(
            [self.threshold_map[r] for r in self.rank_ids],
            dtype=np.float32,
        )
        self.n_ranks = int(self.thresholds.shape[0])

        self.linkage_method = linkage_method
        if linkage_method not in self.AVAILABLE_METHODS:
            raise ValueError(f"linkage_method must be one of {self.AVAILABLE_METHODS}")

        self.delimiter = delimiter
        self.majority_fraction = marjority_fraction

        # Membership representations
        self.memberships_dict: Dict[str, str] = {}              # sample -> full address string
        self.sample_prefix_ids: Dict[str, np.ndarray] = {}      # sample -> int prefix ids per rank
        self.prefix_to_id: List[Dict[str, int]] = [dict() for _ in range(self.n_ranks)]
        self.id_to_prefix: List[List[str]] = [[] for _ in range(self.n_ranks)]
        self.parent_pid: List[Dict[int, int]] = [dict() for _ in range(self.n_ranks)]  # rank -> child->parent

        # New id counters per rank for *prefix IDs* (not the numeric address tokens)
        self._next_prefix_id: List[int] = [0] * self.n_ranks

        # Nomenclature tracker for numeric address tokens level_i
        self.nomenclature_cluster_tracker: Dict[str, int] = {}

        # Query tracking
        self.query_labels: set[str] = set()
        self.query_ids: set[str] = set()

        # Membership parse errors
        self.error_samples = {
            self.ERROR_MISSING_DELIMITER: [],
            self.ERROR_LENGTH: [],
            self.ERROR_NON_INTEGER: [],
        }

        if not is_file_ok(dist_file):
            self.status = False
            raise FileNotFoundError(f"Distance file missing/empty: {dist_file}")
        if not is_file_ok(membership_file):
            self.status = False
            raise FileNotFoundError(f"Membership file missing/empty: {membership_file}")

        df = self._read_membership(membership_file)
        cols = df.columns.values.tolist()
        if sample_col not in cols or address_col not in cols:
            self.status = False
            raise ValueError(f"Membership must include {sample_col} and {address_col}; got {cols}")

        addr_map = df[[sample_col, address_col]].set_index(sample_col).to_dict()[address_col]
        mem_df = self._format_membership_df(addr_map, self.delimiter)
        self._raise_membership_errors()

        self._process_memberships(mem_df)
        self._init_nomenclature_tracker(mem_df)

        # Initial run
        self.assign(n_records=batch_size)

    # ----------------------- membership parsing/loading -----------------------

    def _check_file_type(self, f: str) -> None:
        """Validate file extension is known to genomic_address_service."""
        extension = os.path.splitext(f)[1]
        valid_extensions = list(EXTENSIONS.keys())
        if extension not in valid_extensions:
            self.status = False
            raise ValueError(f"{f} has invalid extension {extension}; expected one of {valid_extensions}")

    def _read_membership(self, f: str) -> pd.DataFrame:
        """Read membership TSV."""
        self._check_file_type(f)
        return pd.read_csv(f, header=0, sep="\t", low_memory=False)

    def _format_membership_df(self, data: Mapping[str, str], delim: str) -> pd.DataFrame:
        """
        Convert membership address strings into a DataFrame of integer levels.

        Returns
        -------
        pd.DataFrame
            Index: sample IDs
            Columns: level_0 .. level_{n-1} (int)
        """
        membership: Dict[str, Dict[str, int]] = {}

        for sample_id, address_str in data.items():
            parts = str(address_str).split(delim)
            if len(parts) != self.n_ranks:
                if (delim not in str(address_str)) and self.n_ranks > 1:
                    self.error_samples[self.ERROR_MISSING_DELIMITER].append(sample_id)
                else:
                    self.error_samples[self.ERROR_LENGTH].append(sample_id)
                continue

            row: Dict[str, int] = {}
            ok = True
            for idx, value in enumerate(parts):
                try:
                    row[f"level_{idx}"] = int(value)
                except Exception:
                    self.error_samples[self.ERROR_NON_INTEGER].append(sample_id)
                    ok = False
                    break

            if ok:
                membership[sample_id] = row

        return pd.DataFrame.from_dict(membership, orient="index")

    def _raise_membership_errors(self) -> None:
        """Raise if membership addresses could not be parsed."""
        if self.error_samples[self.ERROR_MISSING_DELIMITER]:
            self.status = False
            raise ValueError(
                f"{self.ERROR_MISSING_DELIMITER}: {self.error_samples[self.ERROR_MISSING_DELIMITER]}"
            )
        if self.error_samples[self.ERROR_LENGTH]:
            self.status = False
            raise ValueError(
                f"{self.ERROR_LENGTH}: {self.error_samples[self.ERROR_LENGTH]} "
                f"(expected {self.n_ranks} ranks)"
            )
        if self.error_samples[self.ERROR_NON_INTEGER]:
            self.status = False
            raise ValueError(
                f"{self.ERROR_NON_INTEGER}: {self.error_samples[self.ERROR_NON_INTEGER]}"
            )

    # ----------------------- prefix interning (strings->ints) -----------------

    def _intern_prefix(self, rank: int, prefix: str) -> int:
        """
        Intern a prefix string into an integer ID for a specific rank.

        Parameters
        ----------
        rank
            Rank index [0..n_ranks-1].
        prefix
            Prefix string (e.g., '1.2.3').

        Returns
        -------
        int
            Integer prefix ID.
        """
        d = self.prefix_to_id[rank]
        pid = d.get(prefix)
        if pid is not None:
            return pid
        pid = self._next_prefix_id[rank]
        self._next_prefix_id[rank] += 1
        d[prefix] = pid
        self.id_to_prefix[rank].append(prefix)
        return pid

    def _build_prefixes(self, parts: List[str]) -> List[str]:
        """Build prefix strings for each rank from stringified parts."""
        delim = self.delimiter
        prefixes: List[str] = []
        cur = parts[0]
        prefixes.append(cur)
        for p in parts[1:]:
            cur = f"{cur}{delim}{p}"
            prefixes.append(cur)
        return prefixes

    def _process_memberships(self, mem_df: pd.DataFrame) -> None:
        """
        Build:
          - memberships_dict: sample -> full address string
          - sample_prefix_ids: sample -> prefix IDs per rank
          - parent_pid: rank -> child pid -> parent pid
        """
        delim = self.delimiter

        for sample_id, row in mem_df.iterrows():
            parts = [str(int(row[f"level_{i}"])) for i in range(self.n_ranks)]
            full_addr = delim.join(parts)
            self.memberships_dict[sample_id] = full_addr

            prefixes = self._build_prefixes(parts)
            pids = np.empty(self.n_ranks, dtype=np.int64)
            for r in range(self.n_ranks):
                pids[r] = self._intern_prefix(r, prefixes[r])

            # Parent map for prefix IDs
            # rank 0 parent is itself
            self.parent_pid[0].setdefault(int(pids[0]), int(pids[0]))
            for r in range(1, self.n_ranks):
                child = int(pids[r])
                parent = int(pids[r - 1])
                self.parent_pid[r].setdefault(child, parent)

            self.sample_prefix_ids[sample_id] = pids

    def _init_nomenclature_tracker(self, mem_df: pd.DataFrame) -> None:
        """Initialize next numeric cluster IDs per rank from membership max values."""
        tracker = mem_df.max().to_frame().T.to_dict()
        self.nomenclature_cluster_tracker = {}
        for col in tracker:
            self.nomenclature_cluster_tracker[col] = int(tracker[col][0]) + 1

    def add_memberships(self, sample_id: str, address: List[str]) -> None:
        """
        Persist a new membership and update interned prefix structures.
        """
        parts = [str(x) for x in address]
        full_addr = self.delimiter.join(parts)
        self.memberships_dict[sample_id] = full_addr

        prefixes = self._build_prefixes(parts)
        pids = np.empty(self.n_ranks, dtype=np.int64)
        for r in range(self.n_ranks):
            pids[r] = self._intern_prefix(r, prefixes[r])

        self.parent_pid[0].setdefault(int(pids[0]), int(pids[0]))
        for r in range(1, self.n_ranks):
            self.parent_pid[r].setdefault(int(pids[r]), int(pids[r - 1]))

        self.sample_prefix_ids[sample_id] = pids

    # ---------------------------- pass 1 attach -------------------------------

    def _reduce_bucket_pass1(
        self,
        rid_ids: np.ndarray,
        dists: np.ndarray,
        base_ref_set: set[str],
    ) -> Tuple[Optional[ClusterStatsArray], float]:
        """
        Build finest-rank bucket stats for a single qid using only base references.

        Parameters
        ----------
        rid_ids
            object array of rid strings.
        dists
            float32 array of distances aligned to rid_ids.
        base_ref_set
            Snapshot of pre-existing reference sample IDs.

        Returns
        -------
        stats, min_dist
            stats is None if no usable base refs exist for this qid.
        """
        # Filter to base refs, map rid -> finest rank prefix ID
        # This is the remaining Python-object overhead, but it is now *linear* and
        # the heavy reduction is numeric+Numba.
        fine_rank = self.n_ranks - 1

        # Pre-allocate worst-case arrays
        n = rid_ids.shape[0]
        cids = np.empty(n, dtype=np.int64)
        cdists = np.empty(n, dtype=np.float32)

        k = 0
        min_dist = float("inf")
        spids = self.sample_prefix_ids

        for i in range(n):
            rid = rid_ids[i]
            if rid not in base_ref_set:
                continue
            if rid not in spids:
                continue

            d = float(dists[i])
            if d < min_dist:
                min_dist = d

            cids[k] = int(spids[rid][fine_rank])
            cdists[k] = np.float32(d)
            k += 1

        if k == 0:
            return None, float("inf")

        cids = cids[:k]
        cdists = cdists[:k]

        # Sort by cluster id so we can reduce by runs
        order = np.argsort(cids)
        cids_s = cids[order]
        d_s = cdists[order]

        if _HAVE_NUMBA:
            ids, counts, sums, mins, maxs, matched, min_d = reduce_clusters_sorted(
                cids_s, d_s, self.thresholds
            )
        else:
            ids, counts, sums, mins, maxs, matched, min_d = _reduce_clusters_sorted_fallback(
                cids_s, d_s, self.thresholds
            )

        stats = ClusterStatsArray(
            ids=ids,
            counts=counts,
            sums=sums,
            mins=mins,
            maxs=maxs,
            matched=matched,
        )
        return stats, float(min_d)

    def _eligible_array(
        self,
        method: str,
        stats: ClusterStatsArray,
        rank_idx: int,
        threshold: float,
        majority_fraction: float,
    ) -> np.ndarray:
        """
        Vectorized eligibility mask over clusters for a given rank.

        Returns
        -------
        np.ndarray
            Boolean mask of length n_clusters.
        """
        if method == "complete":
            return stats.maxs <= threshold
        if method == "single":
            return stats.mins <= threshold
        if method == "average":
            return stats.means <= threshold

        # majority
        frac = stats.matched[:, rank_idx] / np.maximum(stats.counts.astype(np.float32), 1.0)
        return frac >= np.float32(majority_fraction)

    def _pick_optimal_cluster_array(
        self,
        eligible_ids: np.ndarray,
        stats: ClusterStatsArray,
        *,
        w_mean: float = 1.0,
        w_min: float = 0.5,
        w_size: float = 0.25,
        w_address: float = 0.001,
        eps: float = 1e-12,
    ) -> Optional[int]:
        """
        Pick optimal cluster among eligible prefix IDs.

        Tie-breaking uses last numeric token of the prefix string at the evaluated rank.
        (We decode prefix string only for eligible clusters; eligible set is usually small.)
        """
        if eligible_ids.size == 0:
            return None


        id_to_pos = {int(cid): i for i, cid in enumerate(stats.ids.tolist())}

        # Score and tie-break
        best: Optional[Tuple[float, int, int]] = None  # (score, last_token, cid)

        for cid in eligible_ids.tolist():
            cid_i = int(cid)
            pos = id_to_pos.get(cid_i)
            if pos is None:
                continue

            count = float(stats.counts[pos])
            if count <= 0:
                continue

            mean_d = float(stats.sums[pos] / max(count, 1.0))
            min_d = float(stats.mins[pos])

            last_tok = 0

            size_term = -math.log(max(count, 1.0) + eps)
            score = (w_mean * mean_d) + (w_min * min_d) + (w_size * size_term) + (w_address * last_tok)
            cand = (score, last_tok, cid_i)

            if best is None or cand < best:
                best = cand

        return None if best is None else best[2]

    def _aggregate_bucket_up_array(
        self,
        stats: ClusterStatsArray,
        rank_idx: int,
    ) -> ClusterStatsArray:
        """
        Aggregate ClusterStatsArray from rank_idx to rank_idx-1 using parent_pid mapping.

        Parameters
        ----------
        stats
            Stats at rank_idx.
        rank_idx
            Current rank index (>0).

        Returns
        -------
        ClusterStatsArray
            Aggregated stats at rank_idx-1.
        """
        if rank_idx <= 0:
            return stats

        parent_map = self.parent_pid[rank_idx]
        # Map each child id to parent id
        parent_ids = np.empty_like(stats.ids, dtype=np.int64)
        for i, cid in enumerate(stats.ids.tolist()):
            parent_ids[i] = int(parent_map.get(int(cid), int(cid)))

        # Reduce by parent id: sort + run reduce (cluster count is small relative to rid count)
        order = np.argsort(parent_ids)
        p_sorted = parent_ids[order]
        counts_sorted = stats.counts[order]
        sums_sorted = stats.sums[order]
        mins_sorted = stats.mins[order]
        maxs_sorted = stats.maxs[order]
        matched_sorted = stats.matched[order, :]

        # Manual run reduce (cluster-level; low volume; keep simple)
        out_ids: List[int] = []
        out_counts: List[int] = []
        out_sums: List[float] = []
        out_mins: List[float] = []
        out_maxs: List[float] = []
        out_matched: List[np.ndarray] = []

        cur = int(p_sorted[0])
        c_count = int(counts_sorted[0])
        c_sum = float(sums_sorted[0])
        c_min = float(mins_sorted[0])
        c_max = float(maxs_sorted[0])
        c_mat = matched_sorted[0, :].astype(np.int32, copy=True)

        for i in range(1, p_sorted.shape[0]):
            pid = int(p_sorted[i])
            if pid != cur:
                out_ids.append(cur)
                out_counts.append(c_count)
                out_sums.append(c_sum)
                out_mins.append(c_min)
                out_maxs.append(c_max)
                out_matched.append(c_mat)

                cur = pid
                c_count = int(counts_sorted[i])
                c_sum = float(sums_sorted[i])
                c_min = float(mins_sorted[i])
                c_max = float(maxs_sorted[i])
                c_mat = matched_sorted[i, :].astype(np.int32, copy=True)
            else:
                c_count += int(counts_sorted[i])
                c_sum += float(sums_sorted[i])
                if float(mins_sorted[i]) < c_min:
                    c_min = float(mins_sorted[i])
                if float(maxs_sorted[i]) > c_max:
                    c_max = float(maxs_sorted[i])
                c_mat += matched_sorted[i, :].astype(np.int32, copy=False)

        out_ids.append(cur)
        out_counts.append(c_count)
        out_sums.append(c_sum)
        out_mins.append(c_min)
        out_maxs.append(c_max)
        out_matched.append(c_mat)

        return ClusterStatsArray(
            ids=np.array(out_ids, dtype=np.int64),
            counts=np.array(out_counts, dtype=np.int32),
            sums=np.array(out_sums, dtype=np.float64),
            mins=np.array(out_mins, dtype=np.float32),
            maxs=np.array(out_maxs, dtype=np.float32),
            matched=np.vstack(out_matched).astype(np.int32, copy=False),
        )

    def _attach_to_existing_arrays(
        self,
        qid: str,
        rid_ids: np.ndarray,
        dists: np.ndarray,
        base_ref_set: set[str],
    ) -> Optional[List[str]]:
        """
        Pass 1 assignment for a single qid, using array-based reducer.

        Returns
        -------
        Optional[List[str]]
            Address tokens per rank, or None if not attachable.
        """
        stats, min_dist = self._reduce_bucket_pass1(rid_ids, dists, base_ref_set)
        if stats is None:
            return None

        # If no possible rank can match (min_dist too large), fail fast
        if min_dist > float(np.max(self.thresholds)):
            return None

        method = self.linkage_method
        maj_frac = self.majority_fraction

        # We need rank-specific decoding for last-token tie-break.
        def last_token_for(rank: int, pid: int) -> int:
            prefix = self.id_to_prefix[rank][pid]
            try:
                return int(prefix.split(self.delimiter)[-1])
            except Exception:
                return 10**18

        # Iterate from finest to coarsest rank
        current = stats
        for rank_idx in reversed(range(self.n_ranks)):
            thr = float(self.thresholds[rank_idx])

            if min_dist > thr:
                if rank_idx > 0:
                    current = self._aggregate_bucket_up_array(current, rank_idx)
                continue

            mask = self._eligible_array(method, current, rank_idx, thr, maj_frac)
            eligible_ids = current.ids[mask]
            if eligible_ids.size > 0:
                # Score with address tie-break using last token
                id_to_pos = {int(cid): i for i, cid in enumerate(current.ids.tolist())}
                best: Optional[Tuple[float, int, int]] = None  # (score,last,cid)

                for cid in eligible_ids.tolist():
                    cid_i = int(cid)
                    pos = id_to_pos[cid_i]
                    count = float(current.counts[pos])
                    mean_d = float(current.sums[pos] / max(count, 1.0))
                    min_d = float(current.mins[pos])
                    size_term = -math.log(max(count, 1.0) + 1e-12)
                    last_tok = last_token_for(rank_idx, cid_i)
                    score = (1.0 * mean_d) + (0.5 * min_d) + (0.25 * size_term) + (0.001 * last_tok)
                    cand = (score, last_tok, cid_i)
                    if best is None or cand < best:
                        best = cand

                if best is None:
                    return None

                chosen_pid = best[2]
                chosen_prefix = self.id_to_prefix[rank_idx][chosen_pid]
                ref_parts = chosen_prefix.split(self.delimiter)

                addr: List[Optional[str]] = [None] * self.n_ranks
                for i, v in enumerate(ref_parts):
                    addr[i] = v

                # Fill deeper ranks with new numeric IDs
                for i in range(self.n_ranks):
                    if addr[i] is None:
                        lvl = f"level_{i}"
                        addr[i] = str(self.nomenclature_cluster_tracker[lvl])
                        self.nomenclature_cluster_tracker[lvl] += 1

                return [str(x) for x in addr]

            if rank_idx > 0:
                current = self._aggregate_bucket_up_array(current, rank_idx)

        return None

    # ----------------------- pass 2: de novo clustering -----------------------

    def _get_pairwise(
        self,
        chunk: DistChunk,
        a: str,
        b: str,
        qid_to_idx: Dict[str, int],
    ) -> float:
        """
        Get symmetric distance between qids a and b from the chunk.

        Returns inf if missing.
        """
        if a == b:
            return 0.0
        ia = qid_to_idx.get(a)
        ib = qid_to_idx.get(b)
        if ia is None or ib is None:
            return float("inf")

        # Search within a's rid list for b (linear scan; n is small: ~200)
        rids_a = chunk.rids[ia]
        dists_a = chunk.dists[ia]
        for i in range(rids_a.shape[0]):
            if rids_a[i] == b:
                return float(dists_a[i])

        # Search within b's rid list for a
        rids_b = chunk.rids[ib]
        dists_b = chunk.dists[ib]
        for i in range(rids_b.shape[0]):
            if rids_b[i] == a:
                return float(dists_b[i])

        return float("inf")

    def _cluster_qids_threshold(
        self,
        qids: List[str],
        chunk: DistChunk,
        qid_to_idx: Dict[str, int],
        method: str,
        threshold: float,
    ) -> List[List[str]]:
        """
        Cluster qids at a fixed threshold using linkage method.

        majority is not used here; caller maps majority -> average.
        """
        if len(qids) <= 1:
            return [qids[:]]

        if method == "single":
            uf = UnionFind(qids)
            for i in range(len(qids)):
                a = qids[i]
                for j in range(i + 1, len(qids)):
                    b = qids[j]
                    if self._get_pairwise(chunk, a, b, qid_to_idx) <= threshold:
                        uf.union(a, b)
            return uf.components()

        # Agglomerative average/complete; n~200 so OK.
        n = len(qids)
        active = list(range(n))
        members: Dict[int, List[int]] = {i: [i] for i in range(n)}
        size: Dict[int, int] = {i: 1 for i in range(n)}

        sum_dist: Dict[Tuple[int, int], float] = {}
        max_dist: Dict[Tuple[int, int], float] = {}

        for i in range(n):
            for j in range(i + 1, n):
                d = self._get_pairwise(chunk, qids[i], qids[j], qid_to_idx)
                sum_dist[(i, j)] = d
                max_dist[(i, j)] = d

        def key_pair(x: int, y: int) -> Tuple[int, int]:
            return (x, y) if x < y else (y, x)

        def linkage(ci: int, cj: int) -> float:
            if method == "complete":
                return max_dist[key_pair(ci, cj)]
            return sum_dist[key_pair(ci, cj)] / (size[ci] * size[cj])

        next_id = n
        while True:
            best_pair = None
            best_val = float("inf")
            for i in range(len(active)):
                ci = active[i]
                for j in range(i + 1, len(active)):
                    cj = active[j]
                    v = linkage(ci, cj)
                    if v < best_val:
                        best_val = v
                        best_pair = (ci, cj)

            if best_pair is None or best_val > threshold:
                break

            a_id, b_id = best_pair
            new_id = next_id
            next_id += 1

            members[new_id] = members[a_id] + members[b_id]
            size[new_id] = size[a_id] + size[b_id]

            for other in active:
                if other in (a_id, b_id):
                    continue
                sa = sum_dist[key_pair(a_id, other)]
                sb = sum_dist[key_pair(b_id, other)]
                sum_dist[key_pair(new_id, other)] = sa + sb

                ma = max_dist[key_pair(a_id, other)]
                mb = max_dist[key_pair(b_id, other)]
                max_dist[key_pair(new_id, other)] = ma if ma >= mb else mb

            active = [x for x in active if x not in (a_id, b_id)]
            active.append(new_id)

        clusters: List[List[str]] = []
        for cid in active:
            clusters.append([qids[i] for i in members[cid]])
        return clusters

    def _pass2_denovo_assign(
        self,
        unassigned: List[str],
        chunk: DistChunk,
        qid_to_idx: Dict[str, int],
    ) -> Dict[str, List[str]]:
        """
        De novo assignment for unassigned qids.

        Uses linkage_method except:
          - majority => average
        """
        method = self.linkage_method
        if method == "majority":
            method = "average"

        remaining = sorted(unassigned)
        out: Dict[str, List[str]] = {}

        for rank_idx in reversed(range(self.n_ranks)):
            if len(remaining) <= 1:
                break

            thr = float(self.thresholds[rank_idx])
            clusters = self._cluster_qids_threshold(remaining, chunk, qid_to_idx, method, thr)
            multi = [c for c in clusters if len(c) >= 2]
            if not multi:
                continue

            multi.sort(key=lambda c: min(c))
            for cluster in multi:
                # shared prefix levels 0..rank_idx
                prefix_vals: List[str] = []
                for i in range(rank_idx + 1):
                    lvl = f"level_{i}"
                    prefix_vals.append(str(self.nomenclature_cluster_tracker[lvl]))
                    self.nomenclature_cluster_tracker[lvl] += 1

                for qid in sorted(cluster):
                    addr = prefix_vals[:]
                    for i in range(rank_idx + 1, self.n_ranks):
                        lvl = f"level_{i}"
                        addr.append(str(self.nomenclature_cluster_tracker[lvl]))
                        self.nomenclature_cluster_tracker[lvl] += 1
                    out[qid] = addr

            clustered = {q for c in multi for q in c}
            remaining = [q for q in remaining if q not in clustered]

        for qid in remaining:
            addr: List[str] = []
            for i in range(self.n_ranks):
                lvl = f"level_{i}"
                addr.append(str(self.nomenclature_cluster_tracker[lvl]))
                self.nomenclature_cluster_tracker[lvl] += 1
            out[qid] = addr

        return out

    # ----------------------------- main driver --------------------------------

    def cluster_voting(self, n_records: int = 1000, delim: str = "\t") -> None:
        """
        Two-pass assignment per DistChunk.

        Pass 1: attach to base references (snapshot at chunk start) using Numba reducer.
        Pass 2: de novo among remaining qids (majority treated as average).
        """
        reader_obj = dist_reader(f=self.dist_file, n_records=n_records, delim=delim)

        for chunk in reader_obj.read_data():
            qids = sorted(chunk.qids)
            self.query_ids.update(qids)
            for q in qids:
                self.query_labels.add(q)

            # Index mapping for pass 2 pair lookup
            qid_to_idx = {qid: i for i, qid in enumerate(chunk.qids)}

            # Snapshot base refs (deterministic)
            base_ref_set = set(self.memberships_dict.keys())

            pass1: Dict[str, List[str]] = {}
            unassigned: List[str] = []
            print(base_ref_set)
            for qid in qids:
                if qid in base_ref_set:
                    pass1[qid] = self.memberships_dict[qid].split(self.delimiter)
                    continue
                idx = qid_to_idx[qid]
                addr = self._attach_to_existing_arrays(
                    qid=qid,
                    rid_ids=chunk.rids[idx],
                    dists=chunk.dists[idx],
                    base_ref_set=base_ref_set,
                )
                if addr is None:
                    unassigned.append(qid)
                else:
                    pass1[qid] = addr

            pass2: Dict[str, List[str]] = {}
            if unassigned:
                pass2 = self._pass2_denovo_assign(unassigned, chunk, qid_to_idx)

            for qid in sorted(pass1.keys()):
                self.add_memberships(qid, pass1[qid])
            for qid in sorted(pass2.keys()):
                self.add_memberships(qid, pass2[qid])

    def assign(self, n_records: int = 1000, delim: str = "\t") -> None:
        """Entry point to execute assignment."""
        self.cluster_voting(n_records=n_records, delim=delim)
