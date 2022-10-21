#!/bin/bash
git pull
cp ./docs/index.md ./README.md
git add ./README.md
git commit -m "updated README"
git push
mkdocs gh-deploy
