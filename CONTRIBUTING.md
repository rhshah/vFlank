# Contributing to vflank

Thanks for your interest in improving vflank! This is the short version; the
full [Developer Guide](docs/DEVELOPER.md) covers setup, layout, and how to
extend the package.

## Quick start

```bash
git clone https://github.com/rhshah/vFlank.git
cd vFlank
pip install -e ".[dev]"
python -m pytest          # tests resolve `vflank` from src/ automatically
```

## The quality gate (run before every PR)

```bash
python -m ruff check src tests
python -m mypy src/vflank/core src/vflank/io
python -m pytest
```

All three must pass; CI runs them on Linux + macOS across Python 3.10–3.12.

## Conventions

- Work on a feature branch; keep `core/` pure and I/O-free.
- No silent failures — surface every error path (raise / log / report).
- Add or update tests in the matching `tests/` subtree.
- Match the surrounding style; the coordinate conventions in
  [CLAUDE.md](CLAUDE.md) are the #1 thing to get right when touching flanks.
- End commit messages and PRs with a clear, present-tense summary.

## Reporting issues

Open an issue at https://github.com/rhshah/vFlank/issues with a minimal
reproduction (a few-line MAF/TSV + the command) where possible.
