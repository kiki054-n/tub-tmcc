"""
TUB-TMCC: CH4, NH3, H2O の結合角を F_R で再現する検証スクリプト

VSEPR理論の斥力モデルをTTTの関係項 F_R として解釈し、
CH4, NH3, H2O の結合角を数値的に再現する。

Usage:
    python 09-verify-chemical-bonding-angles.py

Date: 2026-11-07
"""

import numpy as np
from scipy.optimize import minimize_scalar, minimize

# ============================================================
# 実測値
# ============================================================
ANGLE_CH4 = 109.47   # °
ANGLE_NH3 = 107.0    # °
ANGLE_H2O = 104.5    # °


def chord_distance(cos_alpha):
    """単位球上の2点間の弦距離"""
    return np.sqrt(2 - 2 * np.clip(cos_alpha, -1, 1))


# ============================================================
# CH4: 正四面体（4結合対、等強度）
# ============================================================
def verify_CH4():
    """CH4 の結合角を計算（対称性から厳密解）"""
    cos_angle = -1 / 3
    angle = np.degrees(np.arccos(cos_angle))
    return angle


# ============================================================
# NH3: 3結合対 + 1孤電子対（C3v対称性）
# ============================================================
def NH3_geometry(theta_deg):
    """
    NH3 の幾何パラメータを計算
    theta: 孤電子対と結合対の角度（中心から見た角度）
    """
    theta = np.radians(theta_deg)
    cos_bb = 1.5 * np.cos(theta)**2 - 0.5
    alpha_bb = np.degrees(np.arccos(np.clip(cos_bb, -1, 1)))
    return cos_bb, alpha_bb


def NH3_potential(theta_deg, k, n):
    """
    NH3 のポテンシャルエネルギー
    k: 孤電子対/結合対の斥力比
    n: 斥力の距離依存性の指数
    """
    theta = np.radians(theta_deg)
    cos_bb = 1.5 * np.cos(theta)**2 - 0.5
    cos_bl = np.cos(theta)

    r_bb = chord_distance(cos_bb)
    r_bl = chord_distance(cos_bl)

    # 3組の結合対-結合対 + 3組の結合対-孤電子対
    U = 3 * 1 * 1 / (r_bb**n) + 3 * 1 * k / (r_bl**n)
    return U


def verify_NH3(k=1.24, n=1.5):
    """NH3 の平衡結合角を計算"""
    result = minimize_scalar(
        lambda theta: NH3_potential(theta, k, n),
        bounds=(30, 120),
        method='bounded'
    )
    theta_eq = result.x
    _, alpha_bb = NH3_geometry(theta_eq)
    return theta_eq, alpha_bb


# ============================================================
# H2O: 2結合対 + 2孤電子対（C2v対称性）
# ============================================================
def H2O_geometry(theta_deg, phi_deg):
    """
    H2O の各角度を計算
    theta: H-O-H 角の半角
    phi: L-O-L 角の半角
    """
    theta = np.radians(theta_deg)
    phi = np.radians(phi_deg)

    # H-O-H
    alpha_HH = 2 * theta_deg

    # L-O-L
    alpha_LL = 2 * phi_deg

    # H-O-L
    cos_HL = -np.cos(theta) * np.cos(phi)
    alpha_HL = np.degrees(np.arccos(np.clip(cos_HL, -1, 1)))

    return alpha_HH, alpha_LL, alpha_HL


def H2O_potential(params, k, n):
    """
    H2O のポテンシャルエネルギー
    params = (theta_deg, phi_deg)
    """
    theta_deg, phi_deg = params
    theta = np.radians(theta_deg)
    phi = np.radians(phi_deg)

    # H-H
    cos_HH = np.cos(2 * theta)
    r_HH = chord_distance(cos_HH)

    # L-L
    cos_LL = np.cos(2 * phi)
    r_LL = chord_distance(cos_LL)

    # H-L
    cos_HL = -np.cos(theta) * np.cos(phi)
    r_HL = chord_distance(cos_HL)

    # 1組の H-H + 1組の L-L + 4組の H-L
    U = 1 / (r_HH**n) + k**2 / (r_LL**n) + 4 * k / (r_HL**n)
    return U


def verify_H2O(k=1.28, n=2.0):
    """H2O の平衡結合角を計算"""
    result = minimize(
        lambda params: H2O_potential(params, k, n),
        x0=[50, 50],
        bounds=[(30, 70), (30, 70)],
        method='L-BFGS-B'
    )
    theta_eq, phi_eq = result.x
    alpha_HH, alpha_LL, alpha_HL = H2O_geometry(theta_eq, phi_eq)
    return theta_eq, phi_eq, alpha_HH, alpha_LL, alpha_HL


# ============================================================
# k の最適値を探索
# ============================================================
def find_optimal_k_NH3(n=1.5, k_range=None):
    """NH3 で最適な k を探索"""
    if k_range is None:
        k_range = np.linspace(1.0, 3.0, 401)

    best_k = 1.0
    best_err = 1e10

    for k in k_range:
        _, alpha_bb = verify_NH3(k, n)
        err = abs(alpha_bb - ANGLE_NH3)
        if err < best_err:
            best_err = err
            best_k = k

    return best_k, best_err


def find_optimal_k_H2O(n=2.0, k_range=None):
    """H2O で最適な k を探索"""
    if k_range is None:
        k_range = np.linspace(0.5, 3.0, 251)

    best_k = 1.0
    best_err = 1e10

    for k in k_range:
        _, _, alpha_HH, _, _ = verify_H2O(k, n)
        err = abs(alpha_HH - ANGLE_H2O)
        if err < best_err:
            best_err = err
            best_k = k

    return best_k, best_err


# ============================================================
# メイン実行
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("TUB-TMCC: 化学結合角の F_R 検")
    print("=" * 60)

    # CH4
    angle_ch4 = verify_CH4()
    print(f"\n[CH4] 計算値: {angle_ch4:.2f}° | 実測値: {ANGLE_CH4}° | 誤差: {abs(angle_ch4 - ANGLE_CH4):.2f}°")

    # NH3
    best_k_nh3, err_nh3 = find_optimal_k_NH3(n=1.5)
    theta_nh3, angle_nh3 = verify_NH3(k=best_k_nh3, n=1.5)
    print(f"[NH3] 計算値: {angle_nh3:.2f}° | 実測値: {ANGLE_NH3}° | 誤差: {err_nh3:.2f}° | k={best_k_nh3:.3f}")

    # H2O
    best_k_h2o, err_h2o = find_optimal_k_H2O(n=2.0)
    theta_h2o, phi_h2o, angle_hh, angle_ll, angle_hl = verify_H2O(k=best_k_h2o, n=2.0)
    print(f"[H2O] 計算値: {angle_hh:.2f}° | 実測値: {ANGLE_H2O}° | 誤差: {err_h2o:.2f}° | k={best_k_h2o:.3f}")
    print(f"      L-O-L = {angle_ll:.2f}° | H-O-L = {angle_hl:.2f}°")

    print(f"\n[k の一貫性] NH3: {best_k_nh3:.3f} | H2O: {best_k_h2o:.3f} | 差: {abs(best_k_nh3 - best_k_h2o):.3f}")
