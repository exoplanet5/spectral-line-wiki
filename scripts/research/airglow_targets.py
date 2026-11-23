"""Target list of atomic airglow / aurora / meteor / geocoronal emission lines.

Each target gives the ASD spectrum, a nominal wavelength in ASD's default display medium
(vacuum < 2000 A, standard air 2000-20000 A), a search tolerance, and (optionally) the expected
lower/upper term labels so that the correct ASD transition is selected unambiguously.
The nominal wavelengths here are ONLY search keys; every catalogued value is taken from the
downloaded ASD output (see fetch_airglow_asd.py / build_airglow_atomic.py).

group = key into the literature-metadata table (airglow_meta.py).
"""

# (spectrum, nominal_A, tol_A, term_i, term_k, group)
T = []
def add(spec, nominal, group, tol=0.15, ti=None, tk=None, all=False):
    """all=True: keep every transition (ti->tk) within +-tol of nominal (whole multiplet)."""
    T.append(dict(spec=spec, nominal=nominal, tol=tol, ti=ti, tk=tk, group=group, all=all))

# ---------------- O I ----------------
add('O I', 5577.339, 'OI_5577', ti='1D', tk='1S')
add('O I', 6300.304, 'OI_6300', ti='3P', tk='1D')
add('O I', 6363.776, 'OI_6364', ti='3P', tk='1D')
add('O I', 6391.733, 'OI_6392', ti='3P', tk='1D')
add('O I', 2972.29, 'OI_2972', ti='3P', tk='1S')
add('O I', 2958.37, 'OI_2958', ti='3P', tk='1S')
for w in (7771.944, 7774.166, 7775.388):
    add('O I', w, 'OI_7774', ti='5S*', tk='5P')
for w in (8446.247, 8446.359, 8446.758):
    add('O I', w, 'OI_8446', tol=0.05, ti='3S*', tk='3P')
for w in (1302.168, 1304.858, 1306.029):
    add('O I', w, 'OI_1304', ti='3P', tk='3S*')
for w in (1355.598, 1358.512):
    add('O I', w, 'OI_1356', ti='3P', tk='5S*')
add('O I', 1026.96, 'OI_1027', tol=1.3, ti='3P', tk='3D*', all=True)   # 6 components 1025.76-1028.16
add('O I', 989.69, 'OI_989', tol=1.2, ti='3P', tk='3D*', all=True)     # 6 components 988.58-990.80
add('O I', 1152.15, 'OI_1152', ti='1D', tk='1D*')
# ---------------- O II ----------------
for w in (3726.03, 3728.82):
    add('O II', w, 'OII_3727', ti='4S*', tk='2D*')
for w in (7318.92, 7319.99, 7329.67, 7330.73):
    add('O II', w, 'OII_7320', tol=0.1, ti='2D*', tk='2P*')
for w in (2470.22, 2470.34):
    add('O II', w, 'OII_2470', tol=0.05, ti='4S*', tk='2P*')
for w in (832.76, 833.33, 834.47):
    add('O II', w, 'OII_834', tol=0.05, ti='4S*', tk='4P')
# ---------------- N I ----------------
for w in (5197.90, 5200.26):
    add('N I', w, 'NI_5200', ti='4S*', tk='2D*')
for w in (3466.50, 3466.54):
    add('N I', w, 'NI_3466', tol=0.03, ti='4S*', tk='2P*')
for w in (10397.7, 10398.2, 10407.2, 10407.6):
    add('N I', w, 'NI_10400', tol=0.3, ti='2D*', tk='2P*')
add('N I', 1200.13, 'NI_1200', tol=0.7, ti='4S*', tk='4P', all=True)
add('N I', 1493.65, 'NI_1493', tol=1.1, ti='2D*', tk='2P', all=True)
add('N I', 1134.57, 'NI_1134', tol=0.5, ti='4S*', tk='4P', all=True)
add('N I', 1743.99, 'NI_1744', tol=1.4, ti='2P*', tk='2P', all=True)
add('N I', 1243.25, 'NI_1243', tol=0.2, ti='2D*', tk='2D', all=True)
# ---------------- N II ----------------
for w in (6548.05, 6583.45):
    add('N II', w, 'NII_6584', tol=0.1, ti='3P', tk='1D')
add('N II', 5754.59, 'NII_5755', tol=0.1, ti='1D', tk='1S')
for w in (2139.01, 2142.78):   # air (ASD shows air > 2000 A)
    add('N II', w, 'NII_2143', tol=0.3, ti='3P', tk='5S*')
add('N II', 1084.85, 'NII_1085', tol=0.9, ti='3P', tk='3D*', all=True)
# ---------------- He I / He II ----------------
for w in (10829.09, 10830.25, 10830.34):
    add('He I', w, 'HeI_10830', tol=0.03, ti='3S', tk='3P*')
add('He I', 5875.6, 'HeI_5876', tol=0.5, ti='3P*', tk='3D')
add('He I', 3888.6, 'HeI_3889', tol=0.2, ti='3S', tk='3P*')
add('He I', 584.334, 'HeI_584', tol=0.05, ti='1S', tk='1P*')
add('He I', 537.030, 'HeI_537', tol=0.05, ti='1S', tk='1P*')
add('He II', 303.78, 'HeII_304', tol=0.05)
# ---------------- H I ----------------
add('H I', 1215.67, 'HI_Lya', tol=0.02)
add('H I', 1025.72, 'HI_Lyb', tol=0.02)
add('H I', 972.537, 'HI_Lyg', tol=0.02)
add('H I', 6562.8, 'HI_Ha', tol=0.3)
add('H I', 4861.3, 'HI_Hb', tol=0.3)
# ---------------- metals (meteoric layers / meteors) ----------------
for w in (5889.951, 5895.924):
    add('Na I', w, 'NaI_D')
for w in (7664.90, 7698.96):
    add('K I', w, 'KI_7699', tol=0.1)
add('Li I', 6707.84, 'LiI_6708', tol=0.2, all=True)
add('Ca I', 4226.73, 'CaI_4227', tol=0.1)
for w in (3933.66, 3968.47):
    add('Ca II', w, 'CaII_HK', tol=0.1)
for w in (8498.02, 8542.09, 8662.14):
    add('Ca II', w, 'CaII_IRT', tol=0.1)
add('Mg I', 2852.13, 'MgI_2852', tol=0.1)
for w in (2795.53, 2802.70):
    add('Mg II', w, 'MgII_2800', tol=0.1)
for w in (3719.93, 3859.91, 4045.81, 5269.54):
    add('Fe I', w, 'FeI_meteor_core', tol=0.03)
for w in (3581.19, 3734.86, 3737.13, 3820.43, 4063.59, 4271.76, 4307.90, 4383.54, 5328.04, 5371.49):
    add('Fe I', w, 'FeI_meteor_extra', tol=0.03)
for w in (4923.92, 5018.44, 5169.03):
    add('Fe II', w, 'FeII_meteor', tol=0.03)
for w in (6347.11, 6371.37):
    add('Si II', w, 'SiII_meteor', tol=0.05)
for w in (1048.22, 1066.66):
    add('Ar I', w, 'ArI_1048', tol=0.05)

# ---------------- additions from Meier (1991) Table II / FUSE IHB / PALACE / IMO Table 3-7 ----------------
add('O I', 1040.46, 'OI_1040', tol=1.3, ti='3P', tk='3S*', all=True)      # 2p4 3P - 2p3 4s 3S
add('O I', 1641.305, 'OI_1641', tol=0.1, ti='1D', tk='3S*')               # 2p4 1D2 - 3s 3S1
add('O I', 1172.64, 'OI_1173', tol=0.3, ti='1D', tk='3D*', all=True)      # 2p4 1D2 - 3s' 3D
add('O I', 6156.8, 'OI_6157', tol=1.2, ti='5P', tk='5D*', all=True)       # 3p 5P - 4d 5D
add('O I', 6454.9, 'OI_6455', tol=1.6, ti='5P', tk='5S*', all=True)       # 3p 5P - 5s 5S
add('O I', 7002.1, 'OI_7002', tol=0.5, ti='3P', tk='3D*', all=True)       # 3p 3P - 4d 3D
add('O I', 9263.5, 'OI_9264', tol=3.0, ti='5P', tk='5D*', all=True)       # 3p 5P - 3d 5D
add('N II', 916.2, 'NII_916', tol=0.7, ti='3P', tk='3P*', all=True)
add('O II', 616.7, 'OII_617', tol=0.7, ti='2D*', tk='2P', all=True)
add('O II', 539.5, 'OII_539', tol=0.5, ti='4S*', tk='4P', all=True)
add('O II', 718.5, 'OII_718', tol=0.2, ti='2D*', tk='2D', all=True)
add('He I', 522.213, 'HeI_522', tol=0.05, ti='1S', tk='1P*')
for w in (5167.32, 5172.68, 5183.60):
    add('Mg I', w, 'MgI_b', tol=0.05, ti='3P*', tk='3S')
add('Mg II', 4481.2, 'MgII_4481', tol=0.2, ti='2D', tk='2F*', all=True)
for w in (4254.33, 4274.80, 4289.72):
    add('Cr I', w, 'CrI_meteor', tol=0.05, ti='7S', tk='7P*')
add('Mn I', 4030.75, 'MnI_meteor', tol=0.05, ti='6S', tk='6P*')
add('Si II', 4129.5, 'SiII_4130', tol=1.6, ti='2D', tk='2F*', all=True)
