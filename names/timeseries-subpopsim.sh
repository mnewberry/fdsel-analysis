#!/bin/bash

# This component requires the binary parlin.native compiled from parlin.ml.
# Copy parlin.ml into the fdsel source repository and run:
#   ocamlbuild -pkgs pcre,gsl,batteries parlin.native
# Then copy parlin.native here.

./parlin.native
