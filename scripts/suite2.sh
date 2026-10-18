#!/usr/bin/env bash
# Расширенный набор: suite.sh + gen3 (sawtooth/bimodal/fifo/lifo/desc, seed 1..2). Аргументы — как у sim (без путей).
cd "$(dirname "$0")/.."
SETS="materials/tests"
for s in 1 2 3; do
  SETS="$SETS build/gen2_s$s/twophase build/gen2_s$s/sorted2 build/gen2_s$s/invert build/gen2_s$s/tiny build/gen_s$s/alternating build/gen_s$s/periodic build/gen_s$s/phases"
done
for s in 1 2; do for c in sawtooth bimodal fifo lifo desc; do SETS="$SETS build/gen3_s$s/$c"; done; done
./build/bin/sim.exe "$@" $SETS
