#! /bin/bash

# packages for building docs
pip install sphinx 

# Work from the docs directory
cd doc

# Put a working copy of the gh-pages where they are expected
mkdir ../doc_build/
PAGES_DIR=../doc_build/html
git clone -b gh-pages ${PAGES_URL} ${PAGES_DIR} > /dev/null 2>&1

# set up git so it can push
git config --global user.name "Travis-CI"
git config --global user.email "build@travis.org"

# html version is fine, push it later
make html

cd $PAGES_DIR
git add .
git commit -m "update docs"
git push origin gh-pages > /dev/null 2>&1

