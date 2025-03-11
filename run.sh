#!/bin/bash

if [ -f "log.out" ]; then
    rm "log.out"
fi
    ./Usachev_md $1 >> log.out
