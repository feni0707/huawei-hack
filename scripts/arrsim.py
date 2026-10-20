import random
random.seed(5)
N=3000  # shorts per length in window
def dead_frac(seq,p=0.5):
    d=tot=0;run=0
    for L in seq:
        if random.random()<p: run+=L
        else:
            if 0<run<8: d+=run
            elif run>=8: tot+=run
            run=0
    if 0<run<8: d+=run
    elif run>=8: tot+=run
    return d/(d+tot)
def zones():
    s=[]
    for L in range(1,8): s+=[L]*N
    return s
def pairs(a,b):
    s=[]
    for _ in range(N): s+=[a,b]
    return s
def pairzones():
    s=[]
    for a,b in ((1,7),(2,6),(3,5)): s+=pairs(a,b)
    s+=[4]*N
    return s
def tri():  # triples summing to 8-ish: (1,3,4)? 
    s=[]
    for _ in range(N): s+=[1,7,2,6]
    for _ in range(N): s+=[3,5,4,4]
    return s
def blocks(order):
    s=[]
    for _ in range(N): s+=order
    return s
def rnd():
    s=[random.randint(1,7) for _ in range(7*N)]
    return s
for name,seq in [('zones',zones()),('pairzones',pairzones()),('tri',tri()),('rep 1..7',blocks([1,2,3,4,5,6,7])),('rep 7..1',blocks([7,6,5,4,3,2,1])),('rep 1,7,2,6,3,5,4',blocks([1,7,2,6,3,5,4])),('random',rnd())]:
    print(name.ljust(20),'dead frac of victim cells: %.3f'%dead_frac(seq))
