"""Literature metadata for the atomic airglow / aurora / meteor / geocoronal line inventory.

Every number below was read from a downloaded source (see REFS; local copies under
data/raw/airglow_atomic/lit/).  Wording conventions:
  typical_R : zenith (or stated viewing geometry) column emission rate in rayleighs, with the
              observing condition and the reference key in brackets.  Blank = no value found in
              the consulted sources (NOT zero).
  altitude_km : emitting-layer altitude as stated by the source.
Keys in `refs` point to REFS.
"""

REFS = {
 'ASD': 'Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024). NIST Atomic Spectra Database (ver. 5.12, released 7 Nov 2024), NIST, Gaithersburg, MD. https://physics.nist.gov/asd ; doi:10.18434/T4W30F (queried 2026-11-23)',
 'PR72': 'Peck, E. R. & Reeder, K. (1972), Dispersion of air, J. Opt. Soc. Am. 62, 958. doi:10.1364/JOSA.62.000958 (air refractive-index formula used by NIST ASD; ASD Lines help page)',
 'PALACE': 'Noll, S., Schmidt, C., Hannawald, P., Kausch, W., Kimeswenger, S. (2025), PALACE v1.0: Paranal Airglow Line And Continuum Emission model, Geosci. Model Dev. 18, 4353-4398. doi:10.5194/gmd-18-4353-2025 (reference intensities = Cerro Paranal climatology, zenith, 10.7-cm solar flux 100 sfu)',
 'Meier91': 'Meier, R. R. (1991), Ultraviolet spectroscopy and remote sensing of the upper atmosphere, Space Sci. Rev. 58, 1-185. doi:10.1007/BF01206000 ; ADS 1991SSRv...58....1M (Table II: UV dayglow column emission rate, nadir view from 600 km, moderately active Sun; Table III: UV nightglow, nadir from 600 km, equatorial)',
 'Noxon68': 'Noxon, J. F. (1968), Day airglow, Space Sci. Rev. 8, 92-134. ADS 1968SSRv....8...92N (Table I: dayglow zenith intensities and altitudes)',
 'Hunten67': 'Hunten, D. M. (1967), Spectroscopic studies of the twilight airglow, Space Sci. Rev. 6, 493-573. ADS 1967SSRv....6..493H',
 'Borovicka93': 'Borovicka, J. (1993), A fireball spectrum analysis, Astron. Astrophys. 279, 627-645. ADS 1993A&A...279..627B',
 'Borovicka94': 'Borovicka, J. (1994), Line identifications in a fireball spectrum, Astron. Astrophys. Suppl. 103, 83-96. ADS 1994A&AS..103...83B',
 'IMO': 'International Meteor Organization, IMO Photographic Handbook, Part 3: Meteor Spectra, Table 3-7 (lines frequently found in meteor spectra, relative intensities; after Ceplecha 1966, 1971). https://www.imo.net/docs/03spectra.pdf',
 'NASA-TM-2019-220142': 'Moorhead, A., Cooke, B., Blaauw, R., Moser, D., Ehlert, S. (2019), A Meteoroid Handbook for Aerospace Engineers and Managers, NASA/TM-2019-220142 (NTRS 20200000049) ("Meteors occur at altitudes of about 100 km")',
 'Whiter23': 'Whiter, D. K., Partamies, N., Gustavsson, B., Kauristie, K. (2023), The altitude of green OI 557.7 nm and blue N2+ 427.8 nm aurora, Ann. Geophys. 41, 1-12. doi:10.5194/angeo-41-1-2023',
 'Ivanescu21': 'Ivanescu, L., Baibakov, K., O\'Neill, N. T., Blanchet, J.-P., Schulz, K.-H. (2021), Accuracy in starphotometry, Atmos. Meas. Tech. 14, 6561-6599. doi:10.5194/amt-14-6561-2021 (quotes IBC I-IV = 1, 10, 100, 1000 kR of OI 557.7 nm after Chamberlain 1995)',
 'Barthelemy25': 'Barthelemy, M., Robert, E., Lamy, H. (2025), Synthetic spectra of the aurora: N2, N2+, N, N+, O2+ and O emissions, J. Space Weather Space Clim. 15, 19. doi:10.1051/swsc/2025013 (Table 2: modelled column intensities for 1 erg cm-2 s-1 electron precipitation)',
 'FUSE-IHB': 'FUSE Archive Instrument Handbook, Chapter 9: Airglow Appendix, MAST/STScI. https://archive.stsci.edu/fuse/instrumenthandbook/chapter9.html (Tables 9.2-1, 9.2-2)',
 'LK77': 'Lofthus, A. & Krupenie, P. H. (1977), The spectrum of molecular nitrogen, J. Phys. Chem. Ref. Data 6, 113. doi:10.1063/1.555546',
 'Krupenie72': 'Krupenie, P. H. (1972), The spectrum of molecular oxygen, J. Phys. Chem. Ref. Data 1, 423. doi:10.1063/1.3253101',
 'HITRAN2020': 'Gordon, I. E., et al. (2022), The HITRAN2020 molecular spectroscopic database, J. Quant. Spectrosc. Radiat. Transf. 277, 107949. doi:10.1016/j.jqsrt.2021.107949 (16OH line list, global iso id 48, retrieved via https://hitran.org/lbl/api on 2026-11-23)',
 'HH79-WebBook': 'Huber, K. P. & Herzberg, G. (1979), Constants of Diatomic Molecules, as served by NIST Chemistry WebBook SRD 69 (OH: https://webbook.nist.gov/cgi/cbook.cgi?ID=C3352576&Mask=1000; X2Pi G(v) polynomial of Maillard, Chauville et al. 1976, valid v<=5)',
 'Cosby06': 'Cosby, P. C., Sharpee, B. D., Slanger, T. G., Huestis, D. L., Hanuschik, R. W. (2006), High-resolution terrestrial nightglow emission line atlas from UVES/VLT: positions, intensities, and identifications for 2808 lines at 314-1043 nm, J. Geophys. Res. 111, A12307. doi:10.1029/2006JA012023 (cited via PALACE; not downloaded here)',
 'Sharpee05': 'Sharpee, B. D., Slanger, T. G., Cosby, P. C., Huestis, D. L. (2005), The N(2D-4S) 520 nm forbidden doublet in the nightglow: an experimental test of the theoretical intensity ratio, Geophys. Res. Lett. 32, L12106. doi:10.1029/2005GL023044 (cited via PALACE; not downloaded here)',
 'Gardner17': 'Gardner, D. D., Mierkiewicz, E. J., Roesler, F. L., Nossal, S. M., Haffner, L. M. (2017), Constraining Balmer alpha fine structure excitation measured in geocoronal hydrogen observations, J. Geophys. Res. Space Phys. 122, 10727. doi:10.1002/2017JA024055 (cited via PALACE; not downloaded here)',
 'Evans10': 'Evans, W. F. J., Gattinger, R. L., Slanger, T. G., Saran, D. V., Degenstein, D. A., Llewellyn, E. J. (2010), Discovery of the FeO orange bands in the terrestrial night airglow spectrum obtained with OSIRIS on the Odin spacecraft, Geophys. Res. Lett. 37, L22105 (cited via PALACE; not downloaded here)',
 'Gattinger10': 'Gattinger, R. L., et al. (2010), NO2 air afterglow and O and NO densities from Odin-OSIRIS night and ACE-FTS sunset observations in the Antarctic MLT region, J. Geophys. Res. 115, D12301. doi:10.1029/2009JD013205 (cited via PALACE; not downloaded here)',
 'AndersonBarth71': 'Anderson, J. G. & Barth, C. A. (1971), Rocket investigation of the Mg I and Mg II dayglow, J. Geophys. Res. 76, 3723 (cited via Meier 1991; not downloaded here)',
}

# group -> metadata.  'R' may be a dict {nominal_wavelength: text} for component-specific values.
M = {}
def g(key, category, mechanism, altitude_km='', typical_R='', note='', refs=''):
    M[key] = dict(category=category, mechanism=mechanism, altitude_km=altitude_km,
                  typical_R=typical_R, note=note, refs=refs)

g('OI_5577', 'nightglow|dayglow|twilight|aurora|meteor',
  'night: Barth two-step O+O+M->O2*, O2*+O->O(1S); aurora: electron impact + N2(A) energy transfer; meteor trains',
  'night ~97 (PALACE ref. height); aurora peak ~115 (114.84 mean, Whiter23); day 80-250 (Noxon68)',
  'night 163 R (Paranal, 100 sfu) [PALACE]; day ~2 kR [Noxon68]; aurora IBC I-IV = 1, 10, 100, 1000 kR [Ivanescu21]',
  'Forbidden E2 transition 1D2-1S0; same upper level as 2972/2958.  Short-duration meteor trains often dominated by 5577 [Borovicka93].',
  'PALACE;Noxon68;Whiter23;Ivanescu21;Borovicka93;ASD')
g('OI_6300', 'nightglow|dayglow|twilight|aurora',
  'F-region dissociative recombination O+ + O2 -> O2+ + O; O2+ + e -> O(1D) (night); photoelectron impact + photodissociation (day)',
  'peak ~250, ref. 300 (PALACE); day >125 (Noxon68)',
  {6300.304: 'night 123 R (Paranal) [PALACE]; day 2-60 kR [Noxon68]',
   6363.776: 'night 41 R (Paranal) [PALACE]; day: ~1/3 of 6300 (same upper level)',
   6391.733: 'night 0.019 R (PALACE, scaled from 6300 by NIST A ratio)'},
  '[O I] red doublet/triplet (M1+E2), lifetime ~2 min -> quenched below ~200 km.  6392 is J=0 E2-only component.',
  'PALACE;Noxon68;ASD')
for k in ('OI_6364', 'OI_6392'):
    M[k] = M['OI_6300']
g('OI_2972', 'nightglow|dayglow|aurora',
  'shares O(1S) upper level with 5577 (Barth mechanism at night; photoelectron/dissociative processes by day)',
  'as 5577', 'night 9 R (zenith; Hennes 1966 rocket, quoted in Meier91 Table III)',
  'Trans-auroral line 1S0-3P1; intensity ratio to 5577 fixed by A values (Meier91 Sect.2).', 'Meier91;ASD')
g('OI_2958', 'nightglow|dayglow|aurora', 'O(1S) -> 3P2 E2 branch (very weak)', 'as 5577', '',
  'A(2958)/A(2972) = 2.42e-4/7.54e-2 from ASD -> ~0.3% of 2972 (derived)', 'ASD')
g('OI_7774', 'nightglow|dayglow|aurora',
  'radiative recombination O+ + e -> O* (night, F region; quintet cascade); photoelectron impact (day)',
  'F-region (PALACE ref. height 300)', 'night 16.1 R (777 nm feature, Paranal) [PALACE]',
  'Quintet 3s 5S - 3p 5P triplet.  Predicted in dayglow (Noxon68 Sect. 4.3).', 'PALACE;Noxon68;ASD')
g('OI_8446', 'nightglow|dayglow|aurora',
  'night: radiative recombination; day: Bowen fluorescence (solar Ly-beta pumps O I 1027 -> 8446 -> 1304) + photoelectrons; aurora: e-impact on O and O2',
  'F-region (night)', 'night 6.7 R (845 nm feature, Paranal) [PALACE]; aurora model 0.173 kR per erg cm-2 s-1 at 1 keV [Barthelemy25]',
  'Triplet 3s 3S - 3p 3P.  Lower level is upper level of 1304.', 'PALACE;Noxon68;Barthelemy25;ASD')
g('OI_1304', 'dayglow|nightglow|aurora',
  'day: resonant scattering of solar O I 1304 + photoelectron impact; night: O+ + e radiative recombination',
  'day >100 (Noxon68)', 'day 10498 R (nadir from 600 km) [Meier91 T.II]; day 2-6 kR zenith [Noxon68]; night 54 R (nadir, equatorial) [Meier91 T.III]',
  'Optically thick resonance triplet.', 'Meier91;Noxon68;ASD')
g('OI_1356', 'dayglow|nightglow|aurora',
  'day: photoelectron impact on O (spin-forbidden, optically thin); night: O+ + e radiative recombination',
  'day >100 (Noxon68)', 'day 1068 R (nadir 600 km) [Meier91 T.II]; day 0.4 kR zenith [Noxon68]; night 52 R [Meier91 T.III]',
  '5S2 intercombination doublet; used with 1304 for O / F-region remote sensing (Meier91).', 'Meier91;Noxon68;ASD')
g('OI_1027', 'dayglow|nightglow', 'resonance scattering / Bowen fluorescence by solar Ly-beta (nearly coincident with H I 1025.72); photoelectrons',
  'thermosphere', 'day 340 R [Meier91 T.II]; night 2.3 R (blend with H Ly-beta) [Meier91 T.III]', '', 'Meier91;Noxon68;FUSE-IHB;ASD')
g('OI_989', 'dayglow|nightglow', 'resonance scattering + photoelectron impact (day); recombination (night)', 'thermosphere',
  'day 850 R [Meier91 T.II]; night 2.6 R [Meier91 T.III]', '', 'Meier91;FUSE-IHB;ASD')
g('OI_1152', 'dayglow', 'photoelectron impact / scattering from metastable 1D', 'thermosphere', 'day 150 R [Meier91 T.II]', '', 'Meier91;ASD')
g('OI_1040', 'dayglow', 'resonance scattering + photoelectron impact', 'thermosphere', 'day 60 R [Meier91 T.II]', 'Also listed in FUSE dayglow tables.', 'Meier91;FUSE-IHB;ASD')
g('OI_1641', 'dayglow', 'weak branch 3s 3S -> 1D of the 1304 upper level (photoelectron source)', 'thermosphere', 'day 50 R [Meier91 T.II]',
  '1641/1304 ratio used as O-density diagnostic (Meier91).', 'Meier91;ASD')
g('OI_1173', 'dayglow', 'excitation of 3s\' 3D from metastable O(1D)', 'thermosphere', '', 'Listed (no rate) in Meier91 Table II; FUSE Table 9.2-1.', 'Meier91;FUSE-IHB;ASD')
for k in ('OI_6157', 'OI_6455', 'OI_7002', 'OI_9264'):
    g(k, 'nightglow', 'O+ + e radiative recombination (F region)', 'F-region (PALACE ref. 300)', '',
      'Named in PALACE as recombination features (quintets 777 > 926 > 645 > 616 nm; triplets 845, 700 nm); PALACE models them from Escalante & Victor (1992) coefficients, no individual intensities quoted.',
      'PALACE;ASD')
g('OII_3727', 'aurora|dayglow|twilight', 'electron impact / photoionization-excitation of O to O+(2D)', 'F region', '',
  'Forbidden 4S-2D doublet; no intensity found in consulted sources.', 'ASD')
g('OII_7320', 'dayglow|twilight|aurora', 'photoionization-excitation of O producing O+(2P)', 'F region',
  'day (predicted) 220 R (7319/7320 pair) and 300 R (7330/7331 pair) [Hunten67]', 'Unresolved pairs at 7319-7320 and 7330-7331 A.', 'Hunten67;ASD')
g('OII_2470', 'dayglow|aurora', 'O+(2P) -> 4S branch (same upper level as 7320)', 'F region', '', 'No intensity found in consulted sources.', 'ASD')
g('OII_834', 'dayglow', 'photoionization-excitation of O + resonant scattering by O+', 'F region (O+ layer)', 'day 580 R [Meier91 T.II]',
  'Diagnostic of daytime O+ (Meier91 Sect. 5.2).', 'Meier91;ASD')
g('OII_617', 'dayglow', 'photoionization-excitation of O', 'F region', 'day 60 R [Meier91 T.II]', '', 'Meier91;ASD')
g('OII_539', 'dayglow', 'photoionization-excitation of O', 'F region', 'day 115 R [Meier91 T.II]', '', 'Meier91;ASD')
g('OII_718', 'dayglow', 'photoionization-excitation of O', 'F region', 'day 17 R [Meier91 T.II]', '', 'Meier91;ASD')
g('NI_5200', 'nightglow|twilight|dayglow|aurora', 'N(2D): dissociative recombination NO+ + e and O+/O2+ chemistry (night); fluorescence (twilight)',
  '~250 (PALACE ref.); ~300 (Hunten67)',
  {5197.902: 'night 1.91 R (Paranal) [PALACE]; twilight ~20 R and day ~100 R (doublet) [Hunten67]',
   5200.257: 'night 1.08 R (Paranal) [PALACE]; twilight ~20 R and day ~100 R (doublet) [Hunten67]'},
  'Intensity ratio 5198/5200 = 1.759+-0.014 (Sharpee05, via PALACE).  Upper level 2D lifetime ~10 h.', 'PALACE;Hunten67;Sharpee05;ASD')
g('NI_3466', 'dayglow|twilight|aurora', '4S-2P forbidden; fluorescence/dissociative excitation of N2', 'thermosphere', '',
  'Possible weak dayglow line (Meier91 Sect.2); absorption at 3466 + emission at 10400 proposed for twilight 5200 (Hunten67 Table I note 8).', 'Meier91;Hunten67;ASD')
g('NI_10400', 'twilight|aurora|dayglow', '2P -> 2D cascade (see 3466)', 'thermosphere', '', 'No intensity found in consulted sources.', 'Hunten67;ASD')
g('NI_1200', 'dayglow|aurora', 'resonance scattering + photodissociative excitation of N2 + photoelectron impact', 'thermosphere',
  'day 1500 R [Meier91 T.II]', '', 'Meier91;ASD')
g('NI_1134', 'dayglow', 'resonance scattering / dissociative excitation', 'thermosphere', 'day 400 R [Meier91 T.II]', 'In FUSE dayglow tables.', 'Meier91;FUSE-IHB;ASD')
g('NI_1493', 'dayglow|aurora', 'electron impact on N + dissociative excitation of N2', 'thermosphere', 'day 91 R (e + N2) [Meier91 T.II]', '', 'Meier91;ASD')
g('NI_1243', 'dayglow', 'dissociative excitation / electron impact', 'thermosphere', 'day 150 R [Meier91 T.II]', '', 'Meier91;ASD')
g('NI_1744', 'dayglow|aurora', 'dissociative excitation of N2', 'thermosphere', '', 'Not quantified in consulted sources.', 'ASD')
g('NII_6584', 'aurora|dayglow', 'N+(1D) from photo-/electron-impact dissociative ionization of N2', 'F region', '', 'No intensity found in consulted sources.', 'ASD')
g('NII_5755', 'aurora|dayglow', 'N+(1S) -> 1D', 'F region', '', 'No intensity found in consulted sources.', 'ASD')
g('NII_2143', 'dayglow|aurora', 'photodissociative ionization-excitation of N2 (the "auroral mystery feature")', 'thermosphere (>100)',
  'day 300 R (zenith above 100 km, equator, solar max) [Meier91 T.II]', 'Intercombination 3P - 5S2.', 'Meier91;ASD')
g('NII_1085', 'dayglow', 'photoionization-excitation of N2/N', 'thermosphere', 'day 500 R [Meier91 T.II]', 'In FUSE dayglow tables.', 'Meier91;FUSE-IHB;ASD')
g('NII_916', 'dayglow', 'photoionization-excitation of N2/N', 'thermosphere', 'day 185 R [Meier91 T.II]', '', 'Meier91;ASD')
g('HeI_10830', 'twilight|aurora|dayglow', 'resonance scattering of sunlight by metastable He(2 3S) (produced by photoelectrons, He+ recombination)',
  'screening height ~300 (Hunten67)', 'twilight ~1 kR (0.7-2.8 kR) [Hunten67]',
  'Blended with OH(5-2) Q branch at low resolution (Hunten67).', 'Hunten67;ASD')
g('HeI_5876', 'candidate', 'singlet-triplet cascade line of He I', '', '', 'Included at request; no airglow detection found in consulted sources.', 'ASD')
g('HeI_3889', 'twilight|nightglow', 'resonance scattering by He(2 3S)', 'as 10830',
  'twilight ~1 R [Hunten67]; night 0.54 R (UVES, Cosby06 via PALACE)', 'Superposed on N2+ 1N (0,0) 3914 band wing at low resolution.', 'Hunten67;PALACE;Cosby06;ASD')
g('HeI_584', 'dayglow|geocorona', 'resonance scattering of solar He I 584', 'exosphere/thermosphere', 'day 330 R [Meier91 T.II]', '', 'Meier91;FUSE-IHB;ASD')
g('HeI_537', 'dayglow', 'resonance scattering of solar He I 537', 'exosphere/thermosphere', '', 'Listed (no rate) in Meier91 Table II; FUSE Table 9.2-2.', 'Meier91;FUSE-IHB;ASD')
g('HeI_522', 'dayglow', 'resonance scattering of solar He I 522', 'exosphere/thermosphere', 'day 8.5 R [Meier91 T.II]', '', 'Meier91;ASD')
g('HeII_304', 'dayglow|nightglow|geocorona', 'resonance scattering of solar He II 304 by plasmaspheric He+; ring-current precipitation at night',
  'plasmasphere (several Earth radii)', 'day 5 R (zenith from 740 km) [Meier91 T.II]; night 0.1-5 R [Meier91 T.III]',
  'Plasmasphere imaging line (Meier91 Sect. 5.3.2).', 'Meier91;ASD')
g('HI_Lya', 'geocorona|dayglow|nightglow', 'resonant scattering of solar Ly-alpha by geocoronal H (multiply scattered into night side)',
  '>100 to >1000 (geocorona)', 'day 35000 R (nadir 600 km) [Meier91 T.II]; day 5-12 kR zenith [Noxon68]; night 1712 R [Meier91 T.III]',
  'ASD lists fine-structure components and the 1s-2 multiplet line.', 'Meier91;Noxon68;ASD')
g('HI_Lyb', 'geocorona|dayglow|nightglow', 'resonant scattering of solar Ly-beta', 'geocorona',
  'night 2.3 R (blend with O I 1027) [Meier91 T.III]', 'Present even in orbital night (FUSE-IHB).', 'Meier91;FUSE-IHB;ASD')
g('HI_Lyg', 'geocorona|dayglow', 'resonant scattering of solar Ly-gamma', 'geocorona', 'day 14 R [Meier91 T.II]', '', 'Meier91;FUSE-IHB;ASD')
g('HI_Ha', 'geocorona|nightglow|meteor', 'fluorescence of solar Ly-beta (n=3 excitation) in the sunlit geocorona; meteors',
  'geocorona (depends on shadow height)', 'night 5.29 R (Paranal mix of lines of sight) [PALACE]',
  'The two 3p 2P -> 2s 2S components carry ~95% of geocoronal Halpha (Gardner17, via PALACE). Meteor: IMO relative intensity 2*.', 'PALACE;Gardner17;IMO;ASD')
g('HI_Hb', 'geocorona|nightglow', 'cascade/fluorescence in geocorona', 'geocorona', 'night 0.51 R [PALACE]', '', 'PALACE;ASD')
g('NaI_D', 'nightglow|twilight|dayglow|meteor',
  'night: Chapman mechanism Na+O3->NaO+O2, NaO+O->Na(2P)+O2; twilight/day: resonance scattering; meteor ablation',
  '~92 (PALACE); twilight layer ~90-91 (Hunten67); day 85-95 (Noxon68)',
  {5889.951: 'night D2 23.1 R (Paranal) [PALACE]; twilight ~3 kR (D1+D2) [Hunten67]; day 5-40 kR [Noxon68]; meteor IMO rel. 40',
   5895.924: 'night D1 13.4 R (Paranal) [PALACE]; twilight ~3 kR (D1+D2) [Hunten67]; day 5-40 kR [Noxon68]; meteor IMO rel. 35'},
  'D2/D1 nightglow ratio variable (~1.7). Meteoric Na layer.', 'PALACE;Hunten67;Noxon68;IMO;ASD')
g('KI_7699', 'nightglow|twilight', 'night: K analogue of Chapman mechanism; twilight: resonance scattering', '~89 (PALACE)',
  {7664.899: 'night D2 1.58 R (from D1 x 1.67) [PALACE]; twilight 30-37.5 R zenith [Hunten67]',
   7698.965: 'night D1 0.95 R [PALACE]; twilight 30-37.5 R zenith [Hunten67]'},
  '7665 line sits on a strong O2 A-band absorption line (unobservable from ground; Hunten67).', 'PALACE;Hunten67;ASD')
g('LiI_6708', 'twilight', 'resonance scattering by meteoric (and artificially injected) Li', '~80-85 (Saskatoon; Hunten67)',
  'twilight ~20 R natural, much higher after high-altitude nuclear tests [Hunten67]',
  'ASD gives 0.06 A uncertainty (6Li/7Li isotope structure not resolved in ASD entry).', 'Hunten67;ASD')
g('CaI_4227', 'meteor', 'thermal emission in meteor plasma', '~100 (typical meteors)', 'meteor IMO rel. 11',
  'Not seen in twilight (Weill 1966, quoted by Hunten67).', 'IMO;Hunten67;NASA-TM-2019-220142;ASD')
g('CaII_HK', 'twilight|meteor', 'resonance scattering by meteoric Ca+ (sporadic clouds); meteor high-temperature component',
  'twilight emission at shadow heights 100-250 (Hunten67); meteors ~100', 'twilight few to ~100 R [Hunten67]; meteor IMO rel. K 40*, H 35*',
  '* = bright in fast meteors.', 'Hunten67;IMO;NASA-TM-2019-220142;ASD')
g('CaII_IRT', 'candidate', '', '', '', 'Included at request; no airglow/meteor detection documented in consulted sources.', 'ASD')
g('MgI_2852', 'dayglow|meteor', 'resonance scattering by thermospheric Mg', 'E/F region', '', 'Mg detected in dayglow at 2852 (Meier91 table of remote-sensing methods; AndersonBarth71).', 'Meier91;AndersonBarth71;ASD')
g('MgII_2800', 'dayglow|meteor', 'resonance scattering by thermospheric Mg+', 'E and F region', '', 'Mg+ 2798 doublet dayglow (Meier91; AndersonBarth71).', 'Meier91;AndersonBarth71;ASD')
g('MgI_b', 'meteor', 'thermal emission in meteor plasma (main ~4000-4700 K component)', '~100 (typical)', {5167.3213: 'meteor IMO rel. 17', 5172.684: 'meteor IMO rel. 25', 5183.604: 'meteor IMO rel. 28'}, '', 'IMO;Borovicka93;NASA-TM-2019-220142;ASD')
g('MgII_4481', 'meteor', 'high-temperature (~10000 K) meteor component', '~100 (typical)', 'meteor IMO rel. 15*', '', 'IMO;Borovicka93;ASD')
g('FeI_meteor_core', 'meteor', 'thermal emission in meteor plasma (main component, 3500-4700 K in EN151068)', '~100 typical; fireball EN151068 spectrum 57-30 km',
  {3719.93444: 'meteor IMO rel. 10', 3859.91106: 'meteor IMO rel. 11', 4045.81193: 'meteor IMO rel. 10', 5269.5366: 'meteor IMO rel. 14'},
  'Fe I dominates meteor spectra line count (Borovicka94: 305 lines identified).', 'IMO;Borovicka93;Borovicka94;NASA-TM-2019-220142;ASD')
g('FeI_meteor_extra', 'meteor', 'thermal emission in meteor plasma', '~100 typical',
  {3581.1928: '', 3734.86354: 'meteor IMO rel. 8', 3737.13125: 'meteor IMO rel. 9', 3820.42491: 'meteor IMO rel. 9', 4063.59365: 'meteor IMO rel. 9',
   4271.7599: 'meteor IMO rel. 10', 4307.9017: 'meteor IMO rel. 10', 4383.5444: 'meteor IMO rel. 14', 5328.0381: 'meteor IMO rel. 12', 5371.4891: 'meteor IMO rel. 9'},
  '', 'IMO;Borovicka94;ASD')
g('FeII_meteor', 'meteor', 'high-temperature (~10000 K) meteor component', '~100 typical',
  {4923.92123: 'meteor IMO rel. 2*', 5018.4358: 'meteor IMO rel. 3*', 5169.0282: ''}, 'Fe II multiplet a6S-z6P.', 'IMO;Borovicka93;ASD')
g('SiII_meteor', 'meteor', 'high-temperature (~10000 K) meteor component', '~100 typical', {6347.11: 'meteor IMO rel. 6*', 6371.37: 'meteor IMO rel. 3*'}, '', 'IMO;Borovicka93;ASD')
g('SiII_4130', 'meteor', 'high-temperature meteor component', '~100 typical', 'meteor IMO rel. 1* (4131.0 blend)', '', 'IMO;ASD')
g('CrI_meteor', 'meteor', 'thermal emission in meteor plasma', '~100 typical', {4254.33602: 'meteor IMO rel. 9', 4274.7971: 'meteor IMO rel. 8', 4289.7173: 'meteor IMO rel. 7'}, '', 'IMO;Borovicka94;ASD')
g('MnI_meteor', 'meteor', 'thermal emission in meteor plasma', '~100 typical', 'meteor IMO rel. 10', '', 'IMO;Borovicka94;ASD')
g('ArI_1048', 'candidate', 'resonance lines of Ar I', '', '', 'Included at request ("Si, Ar?"); no terrestrial-airglow measurement found in consulted sources.', 'ASD')
