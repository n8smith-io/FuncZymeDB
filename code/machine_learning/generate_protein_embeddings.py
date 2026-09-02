#!/usr/bin/env python3
"""Generate mean-pooled ESMC embeddings for a protein FASTA file."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from workflow_utils import load_fasta, save_embeddings


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate FuncPred-AI ESMC embeddings.")
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", choices=("esmc_300m", "esmc_600m"), default="esmc_600m")
    parser.add_argument("--device", default="cpu", help="Torch device, such as cpu or cuda.")
    return parser.parse_args()


def embed_sequences(records: dict[str, str], model_name: str, device: str) -> np.ndarray:
    import torch
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig

    model = ESMC.from_pretrained(model_name).to(device).eval()
    config = LogitsConfig(sequence=True, return_embeddings=True)
    vectors: list[np.ndarray] = []
    for index, (sequence_id, sequence) in enumerate(records.items(), start=1):
        with torch.no_grad():
            encoded = model.encode(ESMProtein(sequence=sequence))
            output = model.logits(encoded, config)
        tensor = output.embeddings[0].detach().float().cpu()
        residues = tensor[1:-1] if tensor.shape[0] > 2 else tensor
        vectors.append(residues.mean(dim=0).numpy())
        print(f"embedded {index}/{len(records)}: {sequence_id}")
    return np.vstack(vectors) if vectors else np.empty((0, 0), dtype=np.float32)


def main() -> None:
    args = parse_args()
    records = load_fasta(args.fasta)
    if not records:
        raise SystemExit("The input FASTA contains no sequences.")
    matrix = embed_sequences(records, args.model, args.device)
    save_embeddings(args.output, list(records), matrix, args.model)


if __name__ == "__main__":
    main()
