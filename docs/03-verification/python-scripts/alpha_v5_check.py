from fractions import Fraction as F
from itertools import product, combinations_with_replacement as C
import math
v=F(137)+F(35,972)-F(1,110592)
print("V5 =",v,float(v))
c22,s22=F(137035999177,10**9),F(21,10**9); c18=F(137035999084,10**9)
print("vs 2022: %.2f σ   vs 2018: %.2f σ"%(float((v-c22)/s22),float((v-c18)/s22)))
# 検証B 再現
V=[1,2,3,4,5,6,7,8,9,12,20,27,32,36,108,125,128,137]
hits=[t for t in product(V,repeat=4) if abs(F(137)+F(t[0]*t[1],108*t[2])-F(1,108*128*t[3])-c22)<s22]
print("検証B 再現:",len(hits),hits)
# 式の形を変えた場合（形の選び方の自由度）
forms={
 "137 + a·b/(c·d) − 1/(e·f·g)": None,
}
Vs=[2,4,5,6,7,8,9,12,20,108,125,128]
need=c22-137
cnt=0; ex=[]
fr=set()
for a,b in C(Vs,2):
  for c,d in C(Vs,2):
    x=F(a*b,c*d)
    if abs(x-need)>F(1,1000): continue
    for e,f,g in C(Vs,3):
      y=x-F(1,e*f*g)
      if abs(y-need)<s22: fr.add(y); ex.append((a,b,c,d,e,f,g))
print("形を少し一般化（Core=108 を固定しない）: 1σ内の組合せ",len(ex),"／異なる値",len(fr))
print(ex[:8])
