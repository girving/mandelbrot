import sys
src=open(sys.argv[1]).read().split("def g(x)")[0]
sys.argv=[sys.argv[0],'8']
exec(src)
Fl=load('large.out')
for name,mk in [("[n,2]",lambda n:[n,2]),("[2,n]",lambda n:[2,n]),("[n]",lambda n:[n]),("[1,n,2]",None)]:
    if mk is None: continue
    out=[]
    for n in (8,12,14,16,20,32,64,128,256,512,1024):
        x=val(mk(n))
        if x in Fl: out.append("%d:%+.1e"%(n,Fmodel(mk(n))-Fl[x]))
    print("%-6s pred - actual: %s"%(name," ".join(out)))
