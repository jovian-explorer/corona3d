# Releasing corona3d

The steps below publish a version on PyPI. Nothing in the repository uploads automatically; the tag workflow `.github/workflows/release.yml` only attaches the build to a GitHub release.

## Before the first upload

1. Make the repository public, because the project URLs and the images of the README point to it. The README images use relative paths (`docs/*.png`), which PyPI does not resolve; replace them by `https://raw.githubusercontent.com/jovian-explorer/corona3d/main/docs/...` links in the README (or in a copy used for the build).
2. Create accounts on https://test.pypi.org and https://pypi.org with two-factor authentication, and create an API token on each (a token for the whole account is needed only for the first upload, which creates the project; afterwards use a token limited to the project).
3. Check that the name `corona3d` is free on PyPI (`pip download corona3d --no-deps` must fail with "no matching distribution"); the first upload reserves it.
4. Enable the Zenodo GitHub integration for the repository so that the release is archived, and add the DOI to `CITATION.cff` and the README.

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
3. Inspect the wheel (`python -m zipfile -l dist/*.whl`): only `corona3d/` with `default_config.json` and `events/*.json`, and no tests, data, figures, or paper files. The sdist may hold the documentation text and the tests, and no images.
4. Install the wheel into fresh Python 3.12 and 3.13 virtual environments with `pip install "dist/corona3d-X.Y.Z-py3-none-any.whl[all,test]"`, and from a folder outside the repository run `corona3d --version`, `corona3d doctor`, and the test suite (copy `tests/` there). Delete the environments afterwards. Also check that the core install without extras runs `corona3d doctor` and that `corona3d render` names the extra `render`.
5. Run the regression of Fig. 7 (`fig7_regression/check.py`) with the installed package; it must be pixel-identical.
6. Commit as `jovian_explorer <keshavaggarwal2310@gmail.com>`, tag (`git tag vX.Y.Z`), and push the commit and the tag. The tag workflow builds the release on GitHub.
7. Upload to TestPyPI first and test the installation from there:
   ```
   python -m twine upload --repository testpypi dist/*
   python -m venv /tmp/t
   /tmp/t/bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ "corona3d[all]"
   ```
8. Upload to PyPI (`python -m twine upload dist/*`, with `__token__` as the user name and the token as the password), and verify `pip install "corona3d[all]"` in a new environment.
9. Afterwards, change the install text of the README and the user guide from "until then, from GitHub" to the PyPI command, and keep the repository, README, `CITATION.cff`, the web site, and Zenodo consistent. A version uploaded to PyPI cannot be replaced; a mistake needs a new version number.

## Notes

- Python 3.12 is the lowest version because PyThea 1.3 and SunPy 7 and later require it.
- PyThea pulls in large packages (Streamlit, aiapy, and others), which is why it is limited to the extras `render` and `gcs`.
- The default working directory is the platform cache directory (platformdirs); `--home` and `CORONA3D_HOME` override it.
