import sys,collections
rows=[tuple(map(int,l.split())) for l in open(sys.argv[1])]
tot=sum(r[0] for r in rows);print('holes',len(rows),'total free inside',tot)
b=collections.Counter();cells=collections.Counter()
for sz,age,ll,ever in rows:
    cat=('short' if ll<=7 else 'long')
    sb=('1-3' if sz<=3 else '4-7' if sz<=7 else '8-15' if sz<=15 else '16-39' if sz<40 else '40+')
    b[(cat,sb)]+=1;cells[(cat,sb)]+=sz
for k in sorted(b): print(k,b[k],cells[k])
