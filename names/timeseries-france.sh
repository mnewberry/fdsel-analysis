#!/bin/bash

if [[ ! -e out/france ]] ; then mkdir out/france ; fi
unzip -o -q -d out/france inp/nat2018_csv.zip
