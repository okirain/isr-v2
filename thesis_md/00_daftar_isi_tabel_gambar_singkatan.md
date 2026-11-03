# DAFTAR ISI

| Bagian | Halaman |
|--------|---------|
| HALAMAN JUDUL | i |
| PERNYATAAN BEBAS PLAGIASI | ii |
| HALAMAN PENGESAHAN | iii |
| KATA PENGANTAR | iv |
| DAFTAR ISI | v |
| DAFTAR TABEL | vii |
| DAFTAR GAMBAR | viii |
| DAFTAR SINGKATAN | ix |
| INTISARI | x |
| ABSTRACT | xi |
| **BAB I PENDAHULUAN** | **1** |
| I.1 Latar Belakang | 1 |
| I.2 Perumusan Masalah | 4 |
| I.2.1 Batasan Masalah | 5 |
| I.3 Tujuan Penelitian | 6 |
| I.4 Manfaat Penelitian | 7 |
| **BAB II TINJAUAN PUSTAKA** | **8** |
| II.1 Perkembangan Reaktor Ruang Angkasa dan Teknologi Pendinginnya | 8 |
| II.2 KRUSTY sebagai Benchmark Heat Pipe Logam Cair | 12 |
| II.3 Indonesian Space Reactor (ISR) sebagai Kasus Dasar | 16 |
| II.4 Dasar Teori Heat Pipe dan Batas Operasinya | 20 |
| II.5 Model Reduced-Order dan Validasi Termal | 24 |
| II.6 Celah Penelitian dan Kebaruan | 26 |
| **BAB III DASAR TEORI** | **28** |
| III.1 Termohidraulika Dua Fase pada Heat Pipe Logam Cair | 28 |
| III.2 Mekanisme Kapiler pada Kondisi Mikrogravitasi | 31 |
| III.3 Sifat Termofisika Natrium dan Kalium | 33 |
| III.4 Batas Operasi Heat Pipe | 36 |
| III.5 Neraca Energi Nodal dan Jaringan Termal | 40 |
| III.6 Konduksi Radial Bahan Bakar UN | 44 |
| III.7 Metode Numerik dan Kriteria Validasi | 46 |
| **BAB IV PELAKSANAAN PENELITIAN** | **48** |
| IV.1 Desain Penelitian | 48 |
| IV.2 Arsitektur Perangkat Lunak | 49 |
| IV.3 Diskretisasi dan Algoritma Numerik | 52 |
| IV.4 Protokol Benchmark KRUSTY | 54 |
| IV.5 Prosedur Desain ISR-V2 | 57 |
| IV.6 Asumsi Input dan Ketidakpastian | 60 |
| IV.7 Verifikasi dan Validasi | 62 |
| **BAB V HASIL DAN PEMBAHASAN** | **64** |
| V.1 Verifikasi Model | 64 |
| V.2 Validasi Sifat Termofisika | 66 |
| V.3 Kurva Batas Operasi | 68 |
| V.4 Validasi terhadap KRUSTY | 71 |
| V.5 Simulasi Transien Startup | 76 |
| V.6 Desain dan Kinerja ISR-V2 | 79 |
| V.7 Perbandingan dengan Loop NaK-78 | 85 |
| V.8 Analisis Sensitivitas | 89 |
| V.9 Keterbatasan dan Implikasi | 92 |
| **BAB VI KESIMPULAN DAN SARAN** | **96** |
| VI.1 Kesimpulan | 96 |
| VI.2 Saran | 98 |
| DAFTAR PUSTAKA | 100 |
| LAMPIRAN A KODE PROGRAM | 108 |
| LAMPIRAN B TABEL MATERIAL | 114 |

# DAFTAR TABEL

| Nomor | Judul | Halaman |
|-------|-------|---------|
| Tabel 2.1 | Perbandingan teknologi reaktor ruang angkasa dan sistem pendinginnya | 10 |
| Tabel 2.2 | Parameter desain utama KRUSTY | 14 |
| Tabel 2.3 | Parameter dasar Indonesian Space Reactor (ISR) | 18 |
| Tabel 3.1 | Koefisien korelasi sifat termofisika natrium dan kalium | 34 |
| Tabel 3.2 | Parameter geometri heat pipe acuan | 35 |
| Tabel 4.1 | Struktur modul perangkat lunak ISR-V2 | 50 |
| Tabel 4.2 | Titik benchmark KRUSTY yang digunakan | 56 |
| Tabel 4.3 | Asumsi input utama simulasi | 60 |
| Tabel 5.1 | Validasi sifat termofisika pada 1073 K | 67 |
| Tabel 5.2 | Batas operasi heat pipe pada 1050 K | 69 |
| Tabel 5.3 | Perbandingan suhu model dan data KRUSTY | 72 |
| Tabel 5.4 | Metrik validasi KRUSTY | 74 |
| Tabel 5.5 | Hasil pemindaian desain array heat pipe ISR-V2 | 80 |
| Tabel 5.6 | Perbandingan konfigurasi pendinginan ISR | 83 |
| Tabel 5.7 | Respons kegagalan heat pipe | 86 |
| Tabel 5.8 | Neraca massa sistem pendingin | 87 |
| Tabel 5.9 | Sensitivitas parameter desain ISR-V2 | 90 |
| Tabel 6.1 | Ringkasan jawaban tujuan penelitian | 97 |

# DAFTAR GAMBAR

| Nomor | Judul | Halaman |
|-------|-------|---------|
| Gambar 3.1 | Skema heat pipe silindris dan pembagian region aksial | 29 |
| Gambar 3.2 | Jaringan termal nodal model reduced-order | 41 |
| Gambar 4.1 | Arsitektur perangkat lunak simulator ISR-V2 | 50 |
| Gambar 4.2 | Diagram alir penelitian | 51 |
| Gambar 5.1 | Kurva batas operasi heat pipe natrium 3,0 cm | 70 |
| Gambar 5.2 | Validasi model terhadap data KRUSTY | 73 |
| Gambar 5.3 | Kurva startup satu heat pipe ISR-V2 | 78 |
| Gambar 5.4 | Ruang desain array heat pipe ISR-V2 | 82 |
| Gambar 5.5 | Profil suhu radial bahan bakar dan toleransi kegagalan | 84 |
| Gambar 5.6 | Perbandingan neraca massa sistem pendingin | 88 |
| Gambar 5.7 | Sensitivitas suhu puncak bahan bakar | 91 |

# DAFTAR SINGKATAN

| Singkatan | Kepanjangan |
|-----------|-------------|
| B4C | Boron karbida |
| BDF | Backward Differentiation Formula |
| BeO | Berilium oksida |
| BOL | Beginning of Life |
| CFD | Computational Fluid Dynamics |
| CPL | Capillary Pumped Loop |
| DTNTF | Departemen Teknik Nuklir dan Teknik Fisika |
| EM | Electromagnetic |
| EOL | End of Life |
| HP | Heat Pipe |
| ISR | Indonesian Space Reactor |
| ISR-V2 | Konfigurasi ISR dengan pendingin heat pipe natrium pasif |
| KRUSTY | Kilowatt Reactor Using Stirling TechnologY |
| LHP | Loop Heat Pipe |
| MAPE | Mean Absolute Percentage Error |
| MCNP | Monte Carlo N-Particle |
| MLI | Multilayer Insulation |
| NaK-78 | Paduan natrium-kalium eutektik (22 % Na, 78 % K) |
| ODE | Ordinary Differential Equation |
| PCS | Power Conversion System |
| RMSE | Root Mean Square Error |
| ROM | Reduced-Order Model |
| SSA | Spectral Shift Absorber |
| SS | Stainless Steel |
| TRL | Technology Readiness Level |
| UGM | Universitas Gadjah Mada |
| UN | Uranium mononitrida |
