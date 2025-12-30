# Checklist for uploading to PyPI

> **Note**: This checklist should be used when preparing a release. All paths are relative to the project root directory.

## Pre-Release Preparation

### Version and Package Info
- [ ] Update version number in `pyproject.toml` (e.g., X.Y.Z)
- [ ] Check version number in `src/gwaslab/__init__.py` or `g_Sumstats.py` if present
- [ ] Verify all dependencies in `pyproject.toml` are correct and up-to-date
- [ ] Review Python version requirement (`requires-python`)

### Documentation Updates
- [ ] Update `docs/UpdateLogs.md`:
  - [ ] Remove "(Coming soon)" from version header if present
  - [ ] Add release date
  - [ ] Verify all major changes are documented
  - [ ] Review and update changelog entries
- [ ] Update `docs/index.md`:
  - [ ] Verify PyPI badge shows correct version
  - [ ] Check install command is correct (`pip install gwaslab`)
  - [ ] Verify Python version support information is accurate
- [ ] Verify all documentation links are working
- [ ] Review and update any version-specific examples or tutorials

### Code Quality Checks
- [ ] Run test suite from project root: `python test/run_all_tests.py`
- [ ] Verify all tests pass
- [ ] Check for any linting errors
- [ ] Clean up debug messages (remove all `[DEBUG ...]` log statements)
- [ ] Make sure log messages are consistent with code behavior
- [ ] Review breaking changes and ensure backward compatibility notes are documented
- [ ] Verify new features are properly documented
- [ ] Test new functionality manually
- [ ] Verify deprecated features have appropriate warnings

### Release-Specific Checks
- [ ] Test all new features and major changes
- [ ] Verify backward compatibility for existing functionality
- [ ] Check that deprecation warnings are in place for renamed/removed features
- [ ] Verify test modules are included in package distribution if needed
- [ ] Review API changes and ensure they are documented

## Build and Upload

### Build Package
- [ ] Clean previous builds: `rm -rf dist/ build/ *.egg-info/`
- [ ] Build source distribution: `python -m build --sdist`
- [ ] Build wheel: `python -m build --wheel`
- [ ] Verify build artifacts in `dist/` directory
- [ ] Test install from local wheel: `pip install dist/gwaslab-X.Y.Z-py3-none-any.whl --force-reinstall` (replace X.Y.Z with actual version)

### Upload to PyPI
- [ ] Upload to TestPyPI first (recommended): `python -m twine upload --repository testpypi dist/*`
- [ ] Test installation from TestPyPI: `pip install --index-url https://test.pypi.org/simple/ gwaslab==X.Y.Z` (replace X.Y.Z with actual version)
- [ ] Upload to PyPI: `python -m twine upload dist/*`
- [ ] Verify package appears on PyPI: https://pypi.org/project/gwaslab/
- [ ] Test installation from PyPI: `pip install gwaslab==X.Y.Z` (replace X.Y.Z with actual version)

## Post-Release

### Verification
- [ ] Verify PyPI page displays correct version and metadata
- [ ] Test installation in a clean environment
- [ ] Verify documentation website is updated (if applicable)
- [ ] Check that GitHub release notes are created (if using GitHub releases)

### Communication
- [ ] Update any relevant issue trackers
- [ ] Notify users of major changes (if applicable)
- [ ] Update any external documentation or tutorials