"""CLI contract snapshot (plan C1): tests/data/cli_contract.json from ebscore.cli.build_parser().

The kumo-variant-calling-nf modules (design section 3 process table) call the ebscore subcommands with fixed flags.
The snapshot records, per (sub)command: the sorted option strings; per option (keyed by its first long flag) its
flags, required, the action kind (store, store_true, append, extend, help, version, ...), nargs, type, choices and the
default when it is a constant (None, bool, int, finite float, str, or a list of them); the positionals in order; the
mutually exclusive groups; and the subcommands. Help texts and metavars are not part of it: the wording is free, and
argparse formats --help differently between Python versions. tests/test_cli_contract.py rebuilds the structure and
compares it with the JSON.

After a deliberate CLI change, rerun this script, review the printed changes and the modules that call the changed
flags, and commit the JSON with the change:

    python tests/make_cli_contract.py     (env somvar; snapshots the CLI of this checkout's src/)
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).parent
CONTRACT = HERE / "data" / "cli_contract.json"
SCHEMA = "ebscore-cli-contract/1"

# Which kumo-variant-calling-nf module calls which subcommand (design section 3 process table, 2026-11-23).
NF_MODULES = {
    "pileup": "EBSCORE_PILEUP (modules/local/ebscore/pileup)",
    "pon-matrix": "EBSCORE_PONMATRIX (modules/local/ebscore/ponmatrix)",
    "pon-subset": "EBSCORE_PONSUBSET (modules/local/ebscore/ponsubset)",
    "score": "EBSCORE_SCORE (modules/local/ebscore/score)",
    "vcf2table": "EBSCORE_UNION (modules/local/ebscore/union)",
}
TOP_LEVEL_MODULES = "every ebscore module (versions.yml, the image checks `ebscore --version` / `--help`)"


def _is_constant(value) -> bool:
    if value is None or isinstance(value, (bool, int, str)):
        return True
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, (list, tuple)):
        return all(_is_constant(v) for v in value)
    return False


def _key(action: argparse.Action) -> str:
    """An option's key: its first long flag (`--threads` for `-t/--threads`), else its first flag; a positional's dest."""
    longs = [s for s in action.option_strings if s.startswith("--")]
    return (longs or action.option_strings or [action.dest])[0]


def _describe_action(action: argparse.Action, kinds: dict) -> dict:
    kind = kinds.get(type(action), type(action).__name__)
    entry = {
        "flags": sorted(action.option_strings),
        "required": bool(action.required),
        "action": kind,
        "nargs": action.nargs,
        "type": None if action.type is None else getattr(action.type, "__name__", type(action.type).__name__),
        "choices": None if action.choices is None else list(action.choices),
    }
    if kind in ("store_const", "append_const") and _is_constant(action.const):
        entry["const"] = action.const
    if _is_constant(action.default):  # a computed default (e.g. from the environment) is left out
        entry["default"] = action.default
    return entry


def describe_parser(parser: argparse.ArgumentParser) -> dict:
    """The contract of one parser (and, recursively, of its subcommands)."""
    # argparse's own action registry names the kinds: {_StoreAction: 'store', _AppendAction: 'append', ...}
    kinds = {cls: name for name, cls in parser._registries["action"].items() if isinstance(name, str)}
    out = {
        "option_strings": sorted(s for a in parser._actions for s in a.option_strings),
        "options": {},
        "positionals": [],
        "mutually_exclusive": sorted(
            ({"required": bool(g.required), "options": sorted(_key(a) for a in g._group_actions)}
             for g in parser._mutually_exclusive_groups),
            key=lambda g: g["options"],
        ),
    }
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            out["command_required"] = bool(action.required)
            out["commands"] = {name: describe_parser(sub) for name, sub in sorted(action.choices.items())}
        elif action.option_strings:
            out["options"][_key(action)] = _describe_action(action, kinds)
        else:
            out["positionals"].append({"dest": action.dest, **_describe_action(action, kinds)})
    return out


def build_contract(parser: argparse.ArgumentParser | None = None) -> dict:
    """The contract of `parser` (default: ebscore.cli.build_parser()), as plain JSON types."""
    if parser is None:
        from ebscore.cli import build_parser

        parser = build_parser()
    contract = {"schema": SCHEMA, "prog": parser.prog, **describe_parser(parser)}
    return json.loads(json.dumps(contract))  # tuples -> lists: equal to what json.load gives back


def _flat(obj) -> bool:
    def scalar(v):
        return v is None or isinstance(v, (bool, int, float, str))

    def leaf(v):
        return scalar(v) or (isinstance(v, list) and all(map(scalar, v)))

    if isinstance(obj, dict):
        return all(map(leaf, obj.values()))
    return leaf(obj)


def render(contract: dict) -> str:
    """JSON text with one line per option (and per group or flag list), so a changed flag is a one-line diff."""

    def dump(obj, level):
        if _flat(obj):
            return json.dumps(obj, sort_keys=True)
        pad, inner = "  " * level, "  " * (level + 1)
        if isinstance(obj, dict):
            items = [f"{inner}{json.dumps(k)}: {dump(obj[k], level + 1)}" for k in sorted(obj)]
            return "{\n" + ",\n".join(items) + "\n" + pad + "}"
        return "[\n" + ",\n".join(f"{inner}{dump(v, level + 1)}" for v in obj) + "\n" + pad + "]"

    return dump(contract, 0) + "\n"


def _label(path: tuple) -> str:
    # ("commands", "score", "options", "--fisher", "default") -> "ebscore score --fisher default": the container keys
    # "commands"/"options" sit at even positions, the command names and option keys after them
    return " ".join(["ebscore", *(p for i, p in enumerate(path) if not (i % 2 == 0 and p in ("commands", "options")))])


def contract_diff(old, new, path: tuple = ()) -> list[str]:
    """One line per added, removed or changed entry, e.g. `score --fisher default: "sf" -> "legacy"`."""
    if isinstance(old, dict) and isinstance(new, dict):
        lines = []
        for k in sorted(set(old) | set(new)):
            if k not in new:
                lines.append(f"{_label(path + (k,))}: removed")
            elif k not in old:
                lines.append(f"{_label(path + (k,))}: added {json.dumps(new[k], sort_keys=True)}")
            else:
                lines += contract_diff(old[k], new[k], path + (k,))
        return lines
    if old == new:
        return []
    if isinstance(old, list) and isinstance(new, list) and all(isinstance(v, str) for v in old + new):
        gone, added = [v for v in old if v not in new], [v for v in new if v not in old]
        parts = [f"{what} {', '.join(vs)}" for what, vs in (("added", added), ("removed", gone)) if vs]
        if parts:
            return [f"{_label(path)}: {'; '.join(parts)}"]
    return [f"{_label(path)}: {json.dumps(old, sort_keys=True)} -> {json.dumps(new, sort_keys=True)}"]


def affected_modules(old: dict, new: dict) -> list[str]:
    """The kumo-variant-calling-nf modules whose subcommand changed (design section 3)."""
    old_cmds, new_cmds = old.get("commands", {}), new.get("commands", {})
    out = [NF_MODULES.get(c, f"`{c}`: no module calls it yet")
           for c in sorted(set(old_cmds) | set(new_cmds)) if old_cmds.get(c) != new_cmds.get(c)]
    old_top = {k: v for k, v in old.items() if k != "commands"}
    new_top = {k: v for k, v in new.items() if k != "commands"}
    if old_top != new_top:
        out.insert(0, TOP_LEVEL_MODULES)
    return out


def main() -> int:
    src = HERE.parent / "src"
    if (src / "ebscore").is_dir():
        sys.path.insert(0, str(src))  # the CLI of this checkout, not whatever ebscore is installed
    import ebscore

    new = build_contract()
    old = json.loads(CONTRACT.read_text()) if CONTRACT.exists() else None
    CONTRACT.write_text(render(new))
    print(f"wrote {CONTRACT} (ebscore {ebscore.__version__} from {Path(ebscore.__file__).parent})")
    if old is None:
        return 0
    changes = contract_diff(old, new)
    if not changes:
        print("no contract change")
        return 0
    print("changes:", *(f"  {c}" for c in changes), sep="\n")
    print("review the kumo-variant-calling-nf modules that call the changed flags (design section 3 process table):",
          *(f"  {m}" for m in affected_modules(old, new)), sep="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
