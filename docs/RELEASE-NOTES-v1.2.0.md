# AutoZT v1.2.0 — release notes (ready to paste into the GitHub Release form)

**Title**: v1.2.0 — renamed to AutoZT (commands autozt / az, compatibility kept)

**Tag**: v1.2.0

## Changed

- The software is now called **AutoZT**. Package `autozt`, commands `autozt` and the short
  form `az`, environment variables `AUTOZT_*`, configuration in `~/.config/autozt`.
- The previous command names (`pa`, `phonoagent`) remain as thin wrappers that forward to
  the same program, so existing scripts and habits keep working.
- Version banner and metadata report AutoZT 1.2.0.

## Unchanged in behaviour

- Same 18 skills, same engine, same MCP interface (14 tools with an optional read-only
  profile exposing 6, 37 read-only resources, 2 prompts), same approval gate and audit.
- Same measured numbers: interception rate of hazardous actions 100% (8/8), hazardous
  actions executed 0, provenance verification 17/17 byte-identical on a cluster, local
  output determinism 100%.
- Tests: 33 assertions across five pytest suites pass, plus ten end-to-end suites.

## Upgrading

    pip install -e .        # installs autozt, az and pa entry points
    az --version            # AutoZT (autozt) version 1.2.0
    export AUTOZT_LANG=en   # English help and messages

Nothing else to change: configuration files keep their format, and a repository that used
`pa` keeps working through the compatibility wrapper.
