class profil:
    def __init__(self,name, dl, caoptP, camaxP, cm0P, alpha0P, beta):
        self.dl=dl
        self.caoptP=caoptP
        self.camaxP=camaxP
        self.cm0P=cm0P
        self.alpha0P=alpha0P
        self.name=name
        self.beta=beta



ag25=profil('ag25', 0.0758, 0.84, 1.14, -0.0679, -2.75, 0.5)
ag35=profil("ag35", 0.0872, 0.90, 1.17, -0.0455, -2.13, 1.5 )
clarky10=profil("clarky10", 0.100, 0.950, 1.27, -0.0887, -2.23, 1.2 )
goe795=profil("goe795",0.0800, 0.735, 1.02, -0.0561, -2.31, 1.4)
mh32=profil("mh32", 0.866, 0.8, 1.14, -0.0500, -2.21, 1.5 )


'''
AG 25 2.41	7.58	-0.0679	-2.75	F3K	Thermikflug, Floater
AG 35	2.38	8.72	-0.0455	-2.13	RES	Thermikflug, offene Rippenbauweise
Aquila	4.00	9.38	-0.0638	-3.43	RES	Thermikflug, F3B WC 1977
Clark Y	3.55	11.72	-0.0838	-3.53		Thermikflug, Trainer, 1922
DU 86-084/18	1.12	8.45	-0.0203	-1.05	F5B	Speedflug, Laminarprofil, Wölbklappen
E-193	3.54	10.22	-0.0778	-3.39	F3B	Ehemaliges F3B Profil, F3B WC 1979
E-203	2.66	13.65	-0.0864	-3.12		Ehemaliges Großsegler Profil
E-205	3.00	10.48	-0.0462	-2.37	F3B	Ehemaliges F3B Profil, 1981
E-211	2.55	10.96	-0.1142	-4.18		Allround, Großsegler
E-214	4.03	11.11	-0.1536	-5.81		Thermikflug, Ente/Canards
E-374	2.26	10.92	-0.0353	-1.78		Segelkunstflug, Akrobatik
E-387	3.11	9.07	-0.0805	-3.54	F3E	Thermikflug, Elektroflug, 1986
FX 60-100	3.48	10.03	-0.1227	-4.65		Thermikflug, Großsegler, Oldtimer
FX 60-126	3.50	12.59	-0.1247	-4.76		Allround, Großsegler
FX 63-137	5.88	13.53	-0.2334	-9.05		Thermikflug, Großsegler, Oldtimer
Gö 602	3.41	9.95	-0.0984	-4.00	RES	Thermikflug
Gö 795	2.38	8.00	-0.0561	-2.31	RES	Thermikflug
HD-48 A	2.50	8.50	-0.0581	-2.50	F3B	Allround, Thermikflug
HN-1033 A	2.14	7.55	-0.0702	-2.74	F3K	Allround, Thermikflug, SAL
HN-354	1.92	7.88	-0.0699	-2.69	F3B	Allround, F3B WC 2003
HQ 1.0/9	1.00	9.00	-0.0314	-1.28	F3F	Speedflug, Hangflug, Dynamic Soaring
HQ 1.5/9	1.50	9.00	-0.0453	-1.86	F3F	Hangflug, F3B 1979-86
HQ 2.0/9	1.99	9.00	-0.0614	-2.51	F3B	Allround, F3B WC 1987
HQ 2.5/9	2.48	9.00	-0.0764	-3.11	F3B	Allround, F3B WC 1983 & 1985
HQ 3.5/12	3.49	12.00	-0.1065	-4.36	F3J	Thermikflug, ohne Wölbklappen
HQ/W-2.5/11	2.50	11.03	-0.0903	-3.55		Allround, Großsegler
MH 30	1.75	7.85	-0.0418	-1.78	F5B	Speedflug, Allround, F5D, F3D & F3F
MH 32	2.23	8.66	-0.0500	-2.21	F3B	Allround, Thermikflug, F3J
MH 43	1.59	8.47	-0.0148	-0.96	F5D	Elektroflug, Hotliner
NACA 2412	2.00	12.00	-0.0527	-2.11		Motorflug, Trainer
NACA 6409	6.00	9.00	-0.1568	-6.28	F2A	Thermikflug, Freiflug, Klasse Nordic (A2)
RG 14	1.58	8.48	-0.0460	-1.94	F3F	Hangflug
RG 14A 1.4/7.0	1.36	7.00	-0.0388	-1.70	F3F	Speedflug, Hangflug, F3B WC 1989
RG 15	1.76	8.92	-0.0665	-2.60	F3B	Allround, F3B WC 1991-1997
RG 8	2.22	10.80	-0.1020	-3.73		Allround, Großsegler
Ritz 2-30-10	2.00	9.91	-0.0303	-1.56		Hangflug
S 3021	2.96	9.47	-0.0546	-2.62	F3J	Thermikflug
S 4083	3.47	8.01	-0.1001	-4.17	F3K	Thermikflug, HLG
SA 7035	2.58	9.20	-0.0626	-2.73	F3J	Thermikflug, vgl. SD7080
SD 6060	1.87	10.37	-0.0292	-1.52		Segelkunstflug, Akrobatik
SD 7003	1.46	8.51	-0.0408	-1.76	F3B	Allround, Pfeilnurflügel 30°, 1993-95
SD 7037	2.99	9.20	-0.0784	-3.31	F3J	Thermikflug
SD 7080	2.45	9.16	-0.0644	-2.74	F3J	Thermikflug
'''