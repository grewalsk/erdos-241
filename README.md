# erdos

Attempts at open (and some solved) **Erdős problems**, one per git branch.

## Layout

- **`main`** (this branch) — hub only: this README and [`PROBLEMS.md`](PROBLEMS.md) (the index).
  No problem-specific code lives here.
- **`erdos-<number>`** — one branch per problem, each a self-contained project
  (verifier-first code, data, literature notes, writeups). Branch off `main` to start a new one.

Problem numbers follow [erdosproblems.com](https://www.erdosproblems.com/). Where a Lean
target exists, conventions are matched to
[google-deepmind/formal-conjectures](https://github.com/google-deepmind/formal-conjectures).

## Working conventions (per problem branch)

Each problem branch aims to hold:
- a **standalone brute-force verifier** first — every construction is checked against it;
- code with claims labeled **THEOREM** / **COMPUTATION** / **HEURISTIC**;
- a `literature_memo.md`, a `results.md`, data tables, and a final `assessment.md`.

## Start a new problem

```sh
git checkout main
git checkout -b erdos-<number>      # fresh branch off the clean hub
# ... work the problem ...
```

See [`PROBLEMS.md`](PROBLEMS.md) for the index and status of each attempt.
