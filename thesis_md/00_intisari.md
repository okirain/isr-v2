# INTISARI

Indonesian Space Reactor (ISR) merupakan reaktor cepat 500 kWth berpendingin NaK-78 dengan 61 pin bahan bakar uranium mononitrida. Sistem pendingin aktif berbasis pompa elektromagnetik menimbulkan massa, konsumsi daya, getaran, dan titik kegagalan tunggal. Heat pipe natrium menawarkan alternatif pasif, tetapi perilaku transiennya pada mikrogravitasi belum banyak dimodelkan.

Penelitian ini mengembangkan model reduced-order satu dimensi/lumped-parameter dalam Python untuk dinamika termal heat pipe natrium dan kalium pada 800–1300 K. Model mencakup neraca energi nodal bahan bakar, dinding, uap, dan sink, serta batas kapiler, sonik, entrainment, viskos, dan didih. Model divalidasi terhadap data uji 28 jam KRUSTY, lalu diterapkan pada desain ulang ISR dengan 19 heat pipe natrium 3,0 cm berwick arterial.

Validasi menghasilkan MAPE 3,13% dan galat maksimum 7,20% pada enam titik quasi-steady; suhu nominal 8 jam 1081,8 K dibandingkan 1073,2 K. Konfigurasi ISR-V2 memindahkan 26,3 kW per heat pipe dengan margin 1,85, menekan suhu puncak bahan bakar menjadi 1370 K (margin 45,7% terhadap 2523 K), dan tetap di bawah batas saat tiga heat pipe gagal (1795 K). Massa sistem hampir netral (536,4 kg berbanding 538,0 kg), tetapi pompa 6 kW(e) dan titik kegagalan tunggalnya dieliminasi.

**Kata kunci:** heat pipe, model reduced-order, mikrogravitasi, reaktor ruang angkasa, KRUSTY
