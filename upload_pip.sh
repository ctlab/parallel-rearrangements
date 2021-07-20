#!/bin/bash
rm dist/*
python setup.py sdist bdist_wheel
python -m twine upload dist/*