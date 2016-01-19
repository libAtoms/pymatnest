#! /bin/bash

# packages for building docs
pip install sphinx nbconvert numpydoc

# Install in the virtualenv
python setup.py install
cd ../../../

# Work from the docs directory
cd doc

# Put a working copy of the gh-pages where they are expected
PAGES_DIR=../../pymatnest-pages
git clone -b gh-pages ${PAGES_URL} ${PAGES_DIR} > /dev/null 2>&1

# set up git so it can push
git config --global user.name "Travis-CI"
git config --global user.email "build@travis.org"

# For some reason, it won't import from the current directory;
export PYTHONPATH=`pwd`:$PYTHONPATH

# html version is fine, push it later
make docpush

