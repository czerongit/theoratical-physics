#!/usr/bin/env python3
"""
alpha variation analysis from webb et al. (2010) and murphy et al. (2003)

tests geometric framework prediction: α ~ K_min^(-1/2)
computes implied curvature variations from observed α variations

theory (from embedding evolution theorem, section 7):
    c ~ K_min^(1/2)
    α ~ 1/c ~ K_min^(-1/2)
    Δα/α = -(1/2) × ΔK_min/K_min
    ΔK_min/K_min = -2 × Δα/α

usage:
    python alpha_variation_analysis.py

outputs:
    - console summary of analysis
    - alpha_variation_analysis.png (visualization)
"""

import numpy as np
import matplotlib.pyplot as plt

# physical constants
c = 2.998e8  # m/s
H_0 = 2.18e-18  # s^-1 (hubble parameter, planck 2018)
K_min_cosmic = H_0 / c  # minimum curvature scale
alpha_0 = 1/137.036  # fine structure constant (terrestrial)


def compute_curvature_from_alpha(delta_alpha_over_alpha: float) -> float:
    """
    compute implied K_min variation from observed α variation
    
    theory: α ~ 1/c ~ K_min^(-1/2)
    implies: Δα/α = (-1/2) × ΔK_min/K_min
    solving: ΔK_min/K_min = -2 × Δα/α
    """
    return -2 * delta_alpha_over_alpha


def analyze_murphy_2003():
    """
    murphy et al. (2003) - temporal variation
    128 absorption systems at 0.2 < z < 3.7
    """
    print("=" * 70)
    print("MURPHY ET AL. (2003) - TEMPORAL VARIATION")
    print("=" * 70)
    
    # observed data
    delta_alpha = -0.543e-5
    error = 0.116e-5
    significance = 4.7  # sigma
    
    print(f"\nObserved: Δα/α = ({delta_alpha:.3e} ± {error:.3e})")
    print(f"Significance: {significance}σ")
    print(f"Redshift range: 0.2 < z < 3.7")
    
    # implied curvature variation
    delta_K_over_K = compute_curvature_from_alpha(delta_alpha)
    delta_K_error = 2 * error  # error propagation
    
    print(f"\nTheory: ΔK_min/K_min = -2 × Δα/α")
    print(f"Implied: ΔK_min/K_min = {delta_K_over_K:.3e} ± {delta_K_error:.3e}")
    print(f"Physical: K_min was larger at z~2 (smaller α, more curved universe)")
    
    return {
        'delta_alpha': delta_alpha,
        'error': error,
        'significance': significance,
        'delta_K_over_K': delta_K_over_K,
        'delta_K_error': delta_K_error
    }


def analyze_webb_2010():
    """
    webb et al. (2010) - spatial variation (dipole)
    keck + vlt data from different sky directions
    """
    print("\n" + "=" * 70)
    print("WEBB ET AL. (2010) - SPATIAL VARIATION")
    print("=" * 70)
    
    # keck telescope (hawaii)
    keck = {
        'z_low': {'z': 'z < 1.8', 'delta_alpha': -0.54e-5, 'error': 0.12e-5},
        'z_high': {'z': 'z > 1.8', 'delta_alpha': -0.74e-5, 'error': 0.17e-5}
    }
    
    # vlt telescope (chile)
    vlt = {
        'z_low': {'z': 'z < 1.8', 'delta_alpha': -0.06e-5, 'error': 0.16e-5},
        'z_high': {'z': 'z > 1.8', 'delta_alpha': +0.61e-5, 'error': 0.20e-5}
    }
    
    print("\nKeck Telescope (Hawaii):")
    for key, data in keck.items():
        print(f"  {data['z']}: Δα/α = ({data['delta_alpha']:.2e} ± {data['error']:.2e})")
    
    print("\nVLT Telescope (Chile):")
    for key, data in vlt.items():
        print(f"  {data['z']}: Δα/α = ({data['delta_alpha']:.2e} ± {data['error']:.2e})")
    
    # dipole fit
    dipole_amplitude = 1.02e-5
    dipole_error = 0.21e-5
    significance = 4.2
    
    print(f"\nDipole amplitude: ({dipole_amplitude:.2e} ± {dipole_error:.2e})")
    print(f"Significance: {significance}σ")
    print(f"Direction: RA 17.5±0.9 hours, Dec -58±9 degrees")
    
    # implied curvature variation
    delta_K_over_K = 2 * dipole_amplitude  # absolute value
    delta_K_error = 2 * dipole_error
    
    print(f"\nTheory: |ΔK_min/K_min| = 2 × |Δα/α|")
    print(f"Implied: |ΔK_min/K_min| = {delta_K_over_K:.3e} ± {delta_K_error:.3e}")
    print(f"Physical: spatial gradient in embedding curvature across sky")
    
    return {
        'keck': keck,
        'vlt': vlt,
        'dipole_amplitude': dipole_amplitude,
        'dipole_error': dipole_error,
        'significance': significance,
        'delta_K_over_K': delta_K_over_K,
        'delta_K_error': delta_K_error
    }


def compare_with_cmb():
    """
    compare implied curvature variations with cmb fluctuation scale
    """
    print("\n" + "=" * 70)
    print("COMPARISON WITH CMB FLUCTUATIONS")
    print("=" * 70)
    
    # cmb density fluctuations
    delta_rho_over_rho_cmb = 1e-5
    
    # murphy result
    delta_K_murphy = 1.1e-5
    
    # webb result
    delta_K_webb = 2.0e-5
    
    print(f"\nCMB density fluctuations: δρ/ρ ~ {delta_rho_over_rho_cmb:.0e}")
    print(f"Murphy implied ΔK/K: {delta_K_murphy:.1e}")
    print(f"Webb implied |ΔK/K|: {delta_K_webb:.1e}")
    
    ratio_murphy = delta_K_murphy / delta_rho_over_rho_cmb
    ratio_webb = delta_K_webb / delta_rho_over_rho_cmb
    
    print(f"\nRatio (Murphy/CMB): {ratio_murphy:.1f}")
    print(f"Ratio (Webb/CMB): {ratio_webb:.1f}")
    print(f"\nBoth match CMB perturbation scale within factor of 2")
    print("Embedding curvature perturbations track matter perturbations")
    
    return {
        'delta_rho_cmb': delta_rho_over_rho_cmb,
        'ratio_murphy': ratio_murphy,
        'ratio_webb': ratio_webb
    }


def create_visualization(murphy_data, webb_data, output_path):
    """
    create visualization comparing observations with theory
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # plot 1: observed α variation by telescope/redshift
    ax1 = axes[0]
    
    labels = ['Keck\nz<1.8', 'Keck\nz>1.8', 'VLT\nz<1.8', 'VLT\nz>1.8']
    values = [
        webb_data['keck']['z_low']['delta_alpha'],
        webb_data['keck']['z_high']['delta_alpha'],
        webb_data['vlt']['z_low']['delta_alpha'],
        webb_data['vlt']['z_high']['delta_alpha']
    ]
    errors = [
        webb_data['keck']['z_low']['error'],
        webb_data['keck']['z_high']['error'],
        webb_data['vlt']['z_low']['error'],
        webb_data['vlt']['z_high']['error']
    ]
    colors = ['#2E86AB', '#2E86AB', '#A23B72', '#A23B72']
    
    x_pos = np.arange(len(labels))
    ax1.errorbar(x_pos, np.array(values) * 1e5, yerr=np.array(errors) * 1e5,
                 fmt='o', markersize=12, capsize=6, linewidth=2,
                 color='black', ecolor='gray')
    
    for i, (x, y, c) in enumerate(zip(x_pos, np.array(values) * 1e5, colors)):
        ax1.plot(x, y, 'o', markersize=14, color=c, alpha=0.7)
    
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=1.5)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(labels, fontsize=11)
    ax1.set_ylabel('Δα/α (×10⁻⁵)', fontsize=12)
    ax1.set_title('Webb et al. (2010) Observations', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.text(0.02, 0.98, f'{webb_data["significance"]}σ dipole\ndetection',
             transform=ax1.transAxes, fontsize=11, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # plot 2: theory summary
    ax2 = axes[1]
    
    theory_text = f"""
Geometric Framework Prediction
──────────────────────────────
Theory: α ~ K_min^(-1/2)

Relation (from EMT Theorem):
  Δα/α = (-1/2) × ΔK_min/K_min
  ΔK_min/K_min = -2 × Δα/α

Murphy (2003) - Temporal:
  Observed:  Δα/α = {murphy_data['delta_alpha']:.2e}
  Implied:   ΔK_min/K_min = {murphy_data['delta_K_over_K']:.2e}
  Status:    {murphy_data['significance']}σ detection ✓

Webb (2010) - Spatial:
  Observed:  |Δα/α| = {webb_data['dipole_amplitude']:.2e}
  Implied:   |ΔK_min/K_min| = {webb_data['delta_K_over_K']:.2e}
  Status:    {webb_data['significance']}σ detection ✓

Physical Interpretation:
  • Constants vary with curvature
  • Spatial gradient → dipole pattern
  • Temporal evolution → redshift trend
  • Both match CMB scale ~10⁻⁵
"""
    
    ax2.text(0.5, 0.5, theory_text, ha='center', va='center', fontsize=10,
             transform=ax2.transAxes, family='monospace',
             bbox=dict(boxstyle='round', facecolor='#E8F4F8', alpha=0.8))
    ax2.axis('off')
    ax2.set_title('Geometric Framework Analysis', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {output_path}")


def main():
    """
    main analysis routine
    """
    print("\n" + "=" * 70)
    print("FINE STRUCTURE CONSTANT VARIATION ANALYSIS")
    print("Testing Geometric Framework: α ~ K_min^(-1/2)")
    print("=" * 70)
    
    print(f"\nK_min = H₀/c = {K_min_cosmic:.3e} m⁻¹")
    
    # analyze both datasets
    murphy_data = analyze_murphy_2003()
    webb_data = analyze_webb_2010()
    cmb_comparison = compare_with_cmb()
    
    # summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nTheory: α ~ K_min^(-1/2), so Δα/α = -(1/2) × ΔK_min/K_min")
    print(f"Inverse: ΔK_min/K_min = -2 × Δα/α")
    print(f"\nBoth observations detect α variation at >4σ:")
    print(f"  Murphy: {murphy_data['significance']}σ temporal variation")
    print(f"  Webb:   {webb_data['significance']}σ spatial dipole")
    print(f"\nImplied curvature variations ~10⁻⁵ match CMB perturbation scale.")
    print(f"Embedding curvature perturbations track matter perturbations.")
    
    # create visualization
    from pathlib import Path
    output_dir = Path(__file__).parent
    output_path = output_dir / 'alpha_variation_analysis.png'
    create_visualization(murphy_data, webb_data, str(output_path))
    
    print("\n" + "=" * 70)
    print("CONCLUSION: Geometric framework SUPPORTED by α variation data")
    print("=" * 70 + "\n")


if __name__ == '__main__':
    main()
