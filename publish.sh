#/bin/bash
git add docs/UpdateLogs.md
git add docs/index.md
git add src/gwaslab/version.py
git add pyproject.toml
git add publish.sh
git commit -m "publish 3.4.11"
git push
python -m build
python -m twine upload dist/gwaslab-3.4.11*
