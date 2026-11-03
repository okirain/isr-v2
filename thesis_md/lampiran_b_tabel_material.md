# LAMPIRAN B TABEL MATERIAL

Lampiran ini merangkum sifat material yang digunakan dalam simulasi dan analisis desain. Nilai-nilai tersebut merupakan nilai tipikal pada rentang suhu operasi dan digunakan sebagai input model. Untuk analisis struktur atau keselamatan yang membutuhkan akurasi tinggi, data material yang tersertifikasi harus digunakan.

## B.1 Fluida Kerja

**Tabel B.1** Sifat natrium cair pada berbagai suhu

| Suhu (K) | ρ_l (kg/m³) | c_p (J/kg·K) | k (W/m·K) | μ_l (Pa·s) | σ (N/m) | h_fg (MJ/kg) | P_sat (kPa) |
|----------|-------------|--------------|-----------|------------|---------|--------------|-------------|
| 800 | 826 | 1260 | 69,0 | 2,77 × 10⁻⁴ | 0,171 | 4,027 | 1,70 |
| 1000 | 779 | 1280 | 65,0 | 2,10 × 10⁻⁴ | 0,155 | 3,787 | 25,47 |
| 1073 | 762 | 1287 | 63,5 | 1,99 × 10⁻⁴ | 0,149 | 3,700 | 51,52 |
| 1200 | 732 | 1300 | 61,0 | 1,75 × 10⁻⁴ | 0,139 | 3,547 | 128,6 |

**Tabel B.2** Sifat kalium cair pada berbagai suhu

| Suhu (K) | ρ_l (kg/m³) | c_p (J/kg·K) | k (W/m·K) | μ_l (Pa·s) | σ (N/m) | h_fg (MJ/kg) | P_sat (kPa) |
|----------|-------------|--------------|-----------|------------|---------|--------------|-------------|
| 800 | 724 | 800 | 45,4 | 2,24 × 10⁻⁴ | 0,103 | 2,159 | 3,39 |
| 1000 | 680 | 810 | 43,0 | 1,66 × 10⁻⁴ | 0,095 | 1,979 | 38,6 |
| 1032 | 673 | 812 | 42,6 | 1,59 × 10⁻⁴ | 0,094 | 1,950 | 101,3 |
| 1200 | 636 | 820 | 40,6 | 1,36 × 10⁻⁴ | 0,087 | 1,799 | 186,6 |

**Tabel B.3** Sifat NaK-78 (22% Na, 78% K)

| Sifat | Nilai | Satuan |
|-------|-------|--------|
| Titik leleh | 260,5 | K |
| Titik didih (perkiraan) | ~1000 | K |
| Densitas pada 850 K | 733 | kg/m³ |
| Kalor jenis pada 850 K | 904 | J/(kg·K) |
| Konduktivitas pada 850 K | 49,6 | W/(m·K) |
| Viskositas pada 850 K | 2,3 × 10⁻⁴ | Pa·s |

## B.2 Material Struktur dan Bahan Bakar

**Tabel B.4** Sifat material struktur dan bahan bakar

| Material | ρ (g/cm³) | c_p (J/kg·K) | k (W/m·K) | T_leleh (K) | α (10⁻⁶/K) | Keterangan |
|----------|-----------|--------------|-----------|-------------|------------|------------|
| UN (55% U-235) | 14,264 (teoretis) | 200 | 21 pada 800 K; 20 pada 1200 K | 2523 (batas desain) | 8,0 | Bahan bakar |
| UN (90,6% TD) | 12,93 | 200 | 21 | 2523 | 8,0 | Massa pin 2,55 kg |
| Mo-30Re | 11,963 | 250 | 80 | 3063 | 5,5 | Cladding, blok inti, SSA |
| Haynes 230 | 8,97 | 500 | 22 pada 1000 K | 1573 | 13,0 | Shell heat pipe |
| Nikel | 8,90 | 460 | 60 | 1728 | 13,4 | Wick |
| BeO | 2,855 | 1000 | 30 pada 1000 K | 2800 | 8,0 | Reflektor |
| B4C (80% B-10) | 2,510 | 950 | 30 | 2723 | 5,0 | Drum kendali |
| SS316 | 8,00 | 500 | 16 pada 1000 K | 1673 | 16,0 | Bejana vakum |
| Rhenium | 21,02 | 137 | 48 | 3453 | 6,8 | Kanal pendingin, SSA |

## B.3 Data Tambahan untuk Model

**Tabel B.5** Parameter tambahan yang digunakan dalam simulasi

| Parameter | Nilai | Satuan | Sumber/asumsi |
|-----------|-------|--------|----------------|
| Emisivitas efektif MLI, ε_eff | 0,05–0,10 | — | asumsi desain |
| Area MLI efektif KRUSTY | 5,36 × 10⁻³ | m² | kalibrasi |
| Konduktivitas efektif blok Mo-30Re | 80 | W/(m·K) | literatur |
| Konduktivitas efektif wick Ni-Na | 15–25 | W/(m·K) | model |
| Panjang karakteristik spreading | 0,5 × pitch | m | model |
| Faktor keamanan transport | 1,5 | — | praktik desain |
| Kontingensi massa loop | 15 | % | praktik desain |
| Kontingensi massa heat pipe | 10 | % | praktik desain |

## B.4 Catatan Ketidakpastian

Sifat termofisika natrium dan kalium memiliki ketidakpastian sekitar ±5% pada rentang 800–1300 K, terutama untuk viskositas dan tegangan permukaan. Konduktivitas UN dan Mo-30Re memiliki ketidakpastian sekitar ±10% karena bergantung pada densitas, stoikiometri, dan riwayat fabrikasi. Permeabilitas wick arterial merupakan parameter yang paling tidak pasti, dengan rentang 5 × 10⁻¹⁰ hingga 5 × 10⁻⁹ m². Analisis sensitivitas pada Bab V.8 menunjukkan bahwa ambang permeabilitas minimum untuk desain acuan adalah sekitar 0,43 × 10⁻⁹ m², sehingga rentang tersebut masih dapat diterima tetapi harus diverifikasi secara eksperimental.
