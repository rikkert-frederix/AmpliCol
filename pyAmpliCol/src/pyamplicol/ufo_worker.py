from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Isolated UFO model conversion worker")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--restriction", default="default")
    parser.add_argument(
        "--simplify",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    from importlib.metadata import version

    from ufo_model_loader.commands import load_model
    from ufo_model_loader.common import JSONLook

    restriction = None if args.restriction == "default" else args.restriction
    if restriction == "none":
        restriction = "full"
    model, parameter_card = load_model(
        str(args.source.resolve()),
        restriction,
        bool(args.simplify),
        wrap_indices_in_lorentz_structures=True,
    )
    payload = {
        "loader_version": version("ufo-model-loader"),
        "model": json.loads(model.to_json(JSONLook.COMPACT)),
        "parameter_card": {
            name: [value.real, value.imag]
            for name, value in sorted(parameter_card.items())
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
