#/bin/bash

git add docs/UpdateLogs.md
git add docs/index.md
git add src/gwaslab/version.py
git add pyproject.toml
git commit -m "publish"
git push
python -m build
