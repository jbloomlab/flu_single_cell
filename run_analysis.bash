#!/bin/bash

jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 align_and_annotate.ipynb
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 monocle_analysis.ipynb
