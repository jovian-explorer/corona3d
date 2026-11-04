# Releasing corona3d

corona3d is published on PyPI (https://pypi.org/project/corona3d/). A pushed version tag runs `.github/workflows/release.yml`, which builds the source distribution and the wheel, attaches them to a GitHub release, and uploads them to PyPI through a trusted publisher, so that no token is stored anywhere. PyPI knows the publisher as the repository `jovian-explorer/corona3d`, the workflow `release.yml`, and the environment `pypi`.

## For each release

1. Update `__version__` in `src/corona3d/__init__.py`, `version` and `date-released` in `CITATION.cff`, and the header of `CHANGELOG.md` (rename "Unreleased" to the version and date).
2. Run the checks on a clean tree:
   ```
   python -m pytest                              # offline suite in the tree
   CORONA3D_ONLINE=1 python -m pytest tests/test_online.py
   rm -rf build dist
   python -m build                               # writes dist/*.tar.gz and dist/*.whl
   python -m twine check dist/*
   ```
3. Inspect the wheel (`python -m zipfile -l dist/*.whl`): only `corona3d/` with `default_config.json` and `events/*.json`, and no tests, data, figures, or paper files.
4. Commit as `jovian_explorer <keshavaggarwal2310@gmail.com>`, tag (`git tag vX.Y.Z`), and push the commit and the tag. The workflow publishes the GitHub release and the PyPI version.
5. Verify `pip install "corona3d[all]"` in a new environment, and keep the repository, README, `CITATION.cff`, web site, and Zenodo consistent. A version uploaded to PyPI cannot be replaced, so that a mistake needs a new version number.

## Notes

- Python 3.12 is the lowest version because PyThea 1.3 and SunPy 7 and later require it.
- PyThea pulls in large packages (Streamlit, aiapy, and others), which is why it is limited to the extras `render` and `gcs`.
- The default working directory is the platform cache directory (platformdirs); `--home` and `CORONA3D_HOME` override it.
