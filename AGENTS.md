# Agent Notes for shield

This file contains a short tactical summary based on repository evidence. `../llm-wiki/scripts/refresh_wiki.py` rewrites only the machine-managed block.

<!-- BEGIN MACHINE:summary -->
## Quick start
- Repository-local guidance is sufficient: start with `AGENTS.md`, `README.md`, `docs/`, build/test/config files, and the source tree.
- Shield is a C program for performing radiation shielding analyses around a high-energy accelerator. It is a C port of the SHIELD11 program originally developed at SLAC.
- Primary work areas: `doc`, `examples`, `rpm`, `src`.

## Read first
- `README.md`: Primary project overview and workflow notes
- `src/Makefile`: Build system entry point or dependency manifest
- `LICENSE`: Repository configuration that affects local work
- `src/shield.c`: Source file named after the repository

## Build and test
- Documented setup/build commands: `make # builds all dependencies and the shield executable`.
- Detected build systems: GNU Make.
- Unknown: no test workflow evidence was found in the inspected files.
- Likely run commands or operator entry points: `./shield`.

## Operational warnings
- Local checkout layout appears significant; avoid casual changes to sibling-repo assumptions or relative paths.
- Platform-specific dependency setup matters; do not assume one platform's build recipe carries over unchanged.

## Compatibility constraints
- Local builds appear to assume sibling checkouts or fixed relative paths.
- Cross-platform support exists, but platform-specific dependency setup matters.
- Build and runtime behavior likely depends on neighboring core toolkit checkouts.

## Related knowledge
- Repository-local documentation should be treated as authoritative.
- If a shared `llm-wiki/` directory is present in this workspace or parent folder, consult [the matching repo page](../llm-wiki/repos/shield.md) for additional architectural context.
- If no shared wiki is present, continue using repository-local evidence only.
- If available, [the SDDS concept page](../llm-wiki/concepts/sdds.md) adds broader cross-repo context.
- If present in this workspace, [the cross-repo map](../llm-wiki/insights/cross-repo-map.md) helps explain related repositories.
<!-- END MACHINE:summary -->

## Human notes
Add durable repo-specific instructions here.
