#/bin/bash
git add docs/UpdateLogs.md
git add docs/index.md
git add src/gwaslab/version.py
git add pyproject.toml
git commit -m "publish 3.4.0"
git push
python -m build
python -m twine upload gwaslab-3.4.0*
