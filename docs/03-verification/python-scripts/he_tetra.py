import numpy as np
Ry=13.605693; a0=0.529177
def model(theta):
    c=np.cos(theta); d=np.sqrt(2-2*c)          # 電子間距離/r
    rad=(1-c)/d/d**2                            # 電子反発の半径方向成分（ke²/r² 単位）
    tan=np.sin(theta)/d/d**2                    # 接線方向成分
    Z=2-rad; return d,rad,tan,Z
for name,th in [("向かい合う 180°",np.pi),("テトラ 109.47°",np.arccos(-1/3))]:
    d,rad,tan,Z=model(th)
    print(f"{name}: 電子間 d={d:.4f}r  反発(半径方向)={rad:.5f}  反発(接線方向)={tan:.5f}  Zeff={Z:.5f}  r={a0/Z:.4f}Å  E={-2*Z**2*Ry:.3f} eV")
print("2-√6/8 =",2-np.sqrt(6)/8, "  √(2/3)=",np.sqrt(2/3), " d=2√2/√3=",2*np.sqrt(2)/np.sqrt(3))
print("実測 He 全エネルギー -79.005 eV")
