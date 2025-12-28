#!/bin/bash
git pull
cp ./docs/index.md ./README.md
git add ./README.md
git commit -m "updated README"
git push
pip install zensical
zensical build --clean
# Note: For GitHub Pages deployment, use GitHub Actions workflow instead
# or manually push the site/ directory to gh-pages branch

