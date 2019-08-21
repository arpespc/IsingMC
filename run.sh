#!/bin/zsh

for T in {10..90..10}
do
    ./ising -T $T > out.$T  &
done
