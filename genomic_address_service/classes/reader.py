from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, List, Optional, Tuple

import numpy as np


@dataclass
class DistChunk:
    """
    Container for a distance chunk.

    Attributes
    ----------
    qids
        List of query IDs in this chunk.
    rids
        List of NumPy arrays. rids[i] are the reference IDs for qids[i].
    dists
        List of NumPy arrays. dists[i] are the distances aligned to rids[i].
    """
    qids: List[str]
    rids: List[np.ndarray]
    dists: List[np.ndarray]


class dist_reader:
    """
    Streaming reader for pairwise distance files.

    Expected input format (tab-delimited by default)
    -----------------------------------------------
    Header line (ignored), then rows:
        qid <delim> rid <delim> distance

    Notes
    -----
    - This reader groups rows by qid and yields up to `n_records` qids per chunk.
    - Output is array-based for efficient downstream numeric processing.
    """

    def __init__(self, f: str, n_records: int = 1000, delim: str = "\t") -> None:
        self.fpath = f
        self.delim = delim
        self.n_records = n_records

        self._fh = None
        self._row_number = 0

    def _open(self) -> None:
        self._fh = open(self.fpath, "r", encoding="utf-8")

    def _close(self) -> None:
        if self._fh is not None:
            self._fh.close()
        self._fh = None

    def read_pairs(self) -> Iterator[DistChunk]:
        """
        Read a 'pairs' distance file and yield chunks.

        Yields
        ------
        DistChunk
            Chunk containing qids and aligned arrays of (rids, dists) per qid.
        """
        self._open()
        try:
            # consume header
            _ = next(self._fh)

            qids: List[str] = []
            rids_list: List[List[str]] = []
            dists_list: List[List[float]] = []

            qid_to_idx = {}

            for line in self._fh:
                self._row_number += 1
                parts = line.rstrip("\n").split(self.delim)
                if len(parts) < 3:
                    continue

                qid = parts[0]
                rid = parts[1]
                try:
                    dist = float(parts[2])
                except ValueError:
                    continue

                idx = qid_to_idx.get(qid)
                if idx is None:
                    # Start a new qid group; if chunk full, flush.
                    if len(qids) >= self.n_records:
                        yield self._flush(qids, rids_list, dists_list)
                        qids, rids_list, dists_list = [], [], []
                        qid_to_idx = {}

                    qid_to_idx[qid] = len(qids)
                    qids.append(qid)
                    rids_list.append([rid])
                    dists_list.append([dist])
                else:
                    rids_list[idx].append(rid)
                    dists_list[idx].append(dist)

            # flush final
            if qids:
                yield self._flush(qids, rids_list, dists_list)
        finally:
            self._close()

    @staticmethod
    def _flush(
        qids: List[str],
        rids_list: List[List[str]],
        dists_list: List[List[float]],
    ) -> DistChunk:
        """
        Convert buffered python lists into NumPy arrays.

        Distances are converted to float32 for speed/memory.
        """
        rids_arr: List[np.ndarray] = []
        dists_arr: List[np.ndarray] = []
        for rids, dists in zip(rids_list, dists_list):
            rids_arr.append(np.array(rids, dtype=object))
            dists_arr.append(np.array(dists, dtype=np.float32))
        return DistChunk(qids=qids, rids=rids_arr, dists=dists_arr)

    def read_data(self) -> Iterator[DistChunk]:
        """
        Backwards-compatible entrypoint.

        """
        yield from self.read_pairs()
