#!/bin/bash
g++-14 Usachev_md.cpp -fopenmp -O3 -pg -I./include -o Usachev_md

if [ -f "log.out" ]; then
    rm "log.out"
fi

