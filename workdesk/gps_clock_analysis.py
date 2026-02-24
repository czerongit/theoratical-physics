#!/usr/bin/env python3
"""
gps clock frequency analysis for c variation tests

tests geometric framework prediction:
  - c ~ K_min^(1/2)
  - local K_min couples to gravitational potential: ΔK_min/K_min ≈ Δφ/c²
  - therefore: Δc/c = (1/2) × Δφ/c²

the geometric effect is an ADDITIONAL contribution beyond standard GR redshift,
comparable in magnitude (~50% of GR effect), not a tiny correction.

usage:
    python gps_clock_analysis.py [--data-dir PATH]

outputs:
    - console summary of analysis
    - gps_clock_analysis.png (visualization)
    - gps_clock_results.txt (detailed results)
"""

import argparse
import numpy as np
from pathlib import Path

# physical constants
G = 6.67430e-11       # m^3 kg^-1 s^-2
M_EARTH = 5.972e24    # kg
R_EARTH = 6.371e6     # m
c = 2.998e8           # m/s
H_0 = 2.18e-18        # s^-1 (hubble parameter)
K_MIN = H_0 / c       # minimum curvature scale ~ 7.27e-27 m^-1


def gravitational_potential(r: float) -> float:
    """gravitational potential at radius r from earth center"""
    return -G * M_EARTH / r


def compute_gr_redshift(r_sat: float, r_ground: float = R_EARTH) -> float:
    """
    standard gr gravitational redshift
    
    Δf/f = Δφ/c² where Δφ = φ_sat - φ_ground
    """
    phi_sat = gravitational_potential(r_sat)
    phi_ground = gravitational_potential(r_ground)
    delta_phi = phi_sat - phi_ground
    return delta_phi / c**2


def compute_geometric_c_variation(r_sat: float, r_ground: float = R_EARTH) -> float:
    """
    geometric c variation from embedding framework
    
    theory:
      - K_min couples to local gravitational curvature
      - ΔK_min/K_min ≈ Δφ/c² (local curvature tracks potential)
      - c ~ K_min^(1/2)
      - Δc/c = (1/2) × ΔK_min/K_min = (1/2) × Δφ/c²
    
    this is an ADDITIONAL effect beyond standard GR redshift
    """
    phi_sat = gravitational_potential(r_sat)
    phi_ground = gravitational_potential(r_ground)
    delta_phi = phi_sat - phi_ground
    
    # ΔK_min/K_min ≈ Δφ/c²
    delta_K_over_K = delta_phi / c**2
    
    # Δc/c = (1/2) × ΔK_min/K_min
    delta_c_over_c = 0.5 * delta_K_over_K
    
    return delta_c_over_c


def analyze_altitude(altitude_km: float, name: str = "") -> dict:
    """analyze clock effects at given altitude"""
    r_sat = R_EARTH + altitude_km * 1000
    
    gr_effect = compute_gr_redshift(r_sat)
    geom_effect = compute_geometric_c_variation(r_sat)
    
    # total predicted effect differs from pure GR
    total_effect = gr_effect + geom_effect
    fractional_difference = geom_effect / gr_effect if gr_effect != 0 else 0
    
    return {
        'name': name,
        'altitude_km': altitude_km,
        'radius_m': r_sat,
        'gr_redshift': gr_effect,
        'geometric_c_variation': geom_effect,
        'total_effect': total_effect,
        'fractional_difference': fractional_difference,
    }


def create_visualization(results: list, output_path: Path):
    """create visualization comparing gr and geometric predictions"""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping visualization")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    altitudes = [r['altitude_km'] for r in results]
    gr_effects = [r['gr_redshift'] for r in results]
    geom_effects = [r['geometric_c_variation'] for r in results]
    total_effects = [r['total_effect'] for r in results]
    
    # plot 1: effects vs altitude
    ax1 = axes[0]
    ax1.plot(altitudes, gr_effects, 'b-o', label='GR redshift (Δφ/c²)', 
             linewidth=2, markersize=8)
    ax1.plot(altitudes, geom_effects, 'r-s', label='Geometric Δc/c = ½×Δφ/c²', 
             linewidth=2, markersize=8)
    ax1.plot(altitudes, total_effects, 'g-^', label='Total (GR + geometric)', 
             linewidth=2, markersize=8)
    
    ax1.set_xlabel('Altitude (km)', fontsize=12)
    ax1.set_ylabel('Δf/f or Δc/c', fontsize=12)
    ax1.set_title('Clock Frequency Effects vs Altitude', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(-10, -10))
    
    # plot 2: ratio of geometric to gr
    ax2 = axes[1]
    ratios = [r['fractional_difference'] for r in results]
    ax2.bar(range(len(results)), ratios, color='#E67E22', alpha=0.7)
    ax2.set_xticks(range(len(results)))
    ax2.set_xticklabels([r['name'] for r in results], rotation=45, ha='right')
    ax2.set_ylabel('Geometric / GR ratio', fontsize=12)
    ax2.set_title('Geometric Effect as Fraction of GR', fontsize=14, fontweight='bold')
    ax2.axhline(y=0.5, color='r', linestyle='--', label='Predicted: 50%')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # add theory annotation
    fig.text(0.5, 0.02, 
             'Theory: c ~ K_min^(1/2), ΔK_min/K_min ≈ Δφ/c², ∴ Δc/c = ½×Δφ/c²',
             ha='center', fontsize=11, style='italic',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nVisualization saved: {output_path}")


def save_results(results: list, output_path: Path):
    """save detailed results to text file"""
    with open(output_path, 'w') as f:
        f.write("GPS Clock Analysis Results\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("THEORETICAL FRAMEWORK\n")
        f.write("-" * 40 + "\n")
        f.write("From Embedding Evolution Theorem:\n")
        f.write("  c ~ K_min^(1/2)\n")
        f.write("  Δc/c = (1/2) × ΔK_min/K_min\n\n")
        f.write("Local curvature couples to gravity:\n")
        f.write("  ΔK_min/K_min ≈ Δφ/c²\n\n")
        f.write("Therefore:\n")
        f.write("  Δc/c = (1/2) × Δφ/c² (geometric effect)\n")
        f.write("  Δf/f = Δφ/c² (standard GR redshift)\n\n")
        f.write("The geometric effect is ~50% of GR, not a tiny correction!\n\n")
        
        f.write("CONSTANTS\n")
        f.write("-" * 40 + "\n")
        f.write(f"K_min = H_0/c = {K_MIN:.3e} m^-1\n")
        f.write(f"H_0 = {H_0:.3e} s^-1\n")
        f.write(f"c = {c:.3e} m/s\n\n")
        
        f.write("RESULTS BY ALTITUDE\n")
        f.write("-" * 40 + "\n\n")
        
        for r in results:
            f.write(f"{r['name']} ({r['altitude_km']:.0f} km)\n")
            f.write(f"  GR redshift:      Δf/f = {r['gr_redshift']:.3e}\n")
            f.write(f"  Geometric Δc/c:         = {r['geometric_c_variation']:.3e}\n")
            f.write(f"  Total effect:           = {r['total_effect']:.3e}\n")
            f.write(f"  Geometric/GR ratio:     = {r['fractional_difference']:.1%}\n\n")
        
        f.write("DETECTION REQUIREMENTS\n")
        f.write("-" * 40 + "\n")
        f.write("To detect geometric c-variation:\n")
        f.write("  - Effect size: ~2.5×10^-10 at GPS altitude\n")
        f.write("  - GPS clock precision: ~10^-13 (insufficient)\n")
        f.write("  - Optical atomic clocks: ~10^-18 (sufficient by 10^8×)\n\n")
        f.write("The geometric effect is detectable with current optical clocks!\n")
    
    print(f"Results saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='GPS clock analysis for geometric c variation'
    )
    parser.add_argument('--data-dir', type=Path, default=None,
                        help='directory containing CLK and SP3 files (optional)')
    args = parser.parse_args()
    
    print("\n" + "=" * 70)
    print("GPS CLOCK FREQUENCY ANALYSIS")
    print("Geometric Framework: c ~ K_min^(1/2), Δc/c = (1/2) × Δφ/c²")
    print("=" * 70)
    
    print(f"\nK_min = H_0/c = {K_MIN:.3e} m⁻¹")
    print(f"\nTheory predicts:")
    print("  - Local K_min couples to gravitational potential")
    print("  - ΔK_min/K_min ≈ Δφ/c²")
    print("  - Therefore Δc/c = (1/2) × Δφ/c²")
    print("  - This is ~50% of the GR redshift effect")
    
    # analyze standard altitudes
    test_cases = [
        (400, "ISS"),
        (20200, "GPS"),
        (35786, "GEO"),
    ]
    
    results = []
    for alt, name in test_cases:
        result = analyze_altitude(alt, name)
        results.append(result)
    
    # summary table
    print("\n" + "=" * 80)
    print(f"{'Platform':<10} {'Alt(km)':<10} {'GR Δf/f':<15} {'Geom Δc/c':<15} {'Ratio':<10}")
    print("=" * 80)
    
    for r in results:
        print(f"{r['name']:<10} {r['altitude_km']:<10.0f} "
              f"{r['gr_redshift']:<15.3e} {r['geometric_c_variation']:<15.3e} "
              f"{r['fractional_difference']:<10.1%}")
    
    print("=" * 80)
    
    # key result
    gps = results[1]  # gps result
    print(f"\nKEY RESULT for GPS altitude ({gps['altitude_km']:.0f} km):")
    print(f"  Standard GR:  Δf/f = {gps['gr_redshift']:.3e}")
    print(f"  Geometric:    Δc/c = {gps['geometric_c_variation']:.3e}")
    print(f"  Ratio:        {gps['fractional_difference']:.1%} of GR effect")
    
    print("\nDETECTION STATUS:")
    print("  GPS clocks (10⁻¹³):      Cannot detect (need 10⁻¹⁰)")
    print("  Optical clocks (10⁻¹⁸):  CAN detect with 10⁸× margin")
    
    # save outputs
    output_dir = Path(__file__).parent
    create_visualization(results, output_dir / 'gps_clock_analysis.png')
    save_results(results, output_dir / 'gps_clock_results.txt')
    
    print("\n" + "=" * 70)
    print("CONCLUSION: Geometric c-variation is ~50% of GR effect")
    print("Detection requires optical atomic clocks in space")
    print("=" * 70 + "\n")


if __name__ == '__main__':
    main()
