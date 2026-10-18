#!/usr/bin/env python3
# cmp.py binA binB files... : сравнить ответы двух решений на одних и тех же командах (без интерактора)
import sys, subprocess, hashlib
a,b=sys.argv[1],sys.argv[2]; bad=0; n=0
for f in sys.argv[3:]:
    t=open(f).read().split(); k=int(t[1]); inp=("\n".join(t[2:2+k])+"\n0\n").encode()
    oa=subprocess.run([a],input=inp,capture_output=True).stdout; ob=subprocess.run([b],input=inp,capture_output=True).stdout
    n+=1
    if oa!=ob or not oa: bad+=1; print("DIFF",f)
print(f"tests {n} mismatches {bad}")
