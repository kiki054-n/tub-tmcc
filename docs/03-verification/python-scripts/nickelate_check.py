# La3Ni2O7：面内の格子定数 a_p（擬正方、Ni–Ni 間）と Tc の照合（文献値）
import numpy as np
rows=[ # 名前, a(Å), b(Å), c(Å) または a_p 直接, Tc_onset(K), 角 Ni-O-Ni, 出典
 ("バルク 常圧（Amam）",5.439,5.376,20.403,None,"0（超電導なし）","168°","arXiv:2503.07929"),
 ("バルク 高圧（Fmmm）",5.289,5.218,19.734,None,"約80（14〜43 GPa）","180°","arXiv:2503.07929 / 2309.15078"),
 ("薄膜 La3Ni2O7 / SLAO",None,None,None,3.77,"約42（常圧）","約180°","arXiv:2507.13694, Ko et al. Nature 2025"),
 ("薄膜 (La,Pr)3Ni2O7 / SLAO",None,None,None,3.75,">40（常圧）","約180°","arXiv:2507.13694, Zhou et al."),
]
print("試料                         a_p(Å)   c/a_p   Tc onset(K)        Ni-O-Ni")
for n,a,b,c,ap,tc,ang,src in rows:
    if ap is None: ap=(a+b)/2/np.sqrt(2); ca=c/ap; cas=f"{ca:.3f}"
    else: cas="  -  "
    print(f"{n:26s} {ap:6.3f}  {cas}   {tc:18s} {ang}")
print("\nNi–頂点O 距離：2.297 Å（常圧）→ 2.122 Å（32.5 GPa）  変化 %.1f%%"%((2.122/2.297-1)*100))
print("a_p の変化（常圧→高圧）： %.1f%%"%(((5.289+5.218)/(5.439+5.376)-1)*100))
