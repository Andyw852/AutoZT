# Releasing

1. Update CHANGELOG.md and the version in pyproject.toml and autozt/__init__.py.
2. Run the test suite: pytest -q (and pytest -q -m cluster with a configured cluster).
3. Commit, then tag: git tag -a vX.Y.Z -m "AutoZT X.Y.Z" && git push --tags.
4. Create the GitHub release from the tag; Zenodo picks it up through .zenodo.json and
   mints a DOI. Put the DOI badge and the concept DOI into README.md, README.en.md and
   CITATION.cff (identifiers section).
5. Build and optionally upload the distribution: python -m build.
