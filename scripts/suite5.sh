#!/usr/bin/env bash
# Набор №5: suite4 + concat (build/gen3_s{1,2,3}/concat). Бинарник: $SIM.
cd "$(dirname "$0")/.."
SETS="materials/tests"
for s in 1 2 3; do
  SETS="$SETS build/gen2_s$s/twophase build/gen2_s$s/sorted2 build/gen2_s$s/invert build/gen2_s$s/tiny build/gen_s$s/alternating build/gen_s$s/periodic build/gen_s$s/phases build/gen3_s$s/concat"
done
for s in 1 2; do for c in sawtooth bimodal fifo lifo desc sortedbig; do SETS="$SETS build/gen3_s$s/$c"; done; done
for m in inv scr shift half dbl; do SETS="$SETS build/var/$m"; done
${SIM:-./build/bin/sim.exe} "$@" $SETS
