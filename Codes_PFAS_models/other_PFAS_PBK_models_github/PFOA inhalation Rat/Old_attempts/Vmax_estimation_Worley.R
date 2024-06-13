MKC = 0.0084	#!fraction mass of kidney (percent of BW); Brown 1997
protein = 2.0e-6	#!amount of protein in proximal tubule cells (mg protein/proximal tubule cell)
MK = MKC*BW*1000#	!mass of the kidney (g)
PTC = MK*6e7	#!number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
Vmax_apical_invitro = 9.3 #!Vmax of apical transporter (pmol/mg protein/min); invitro value for Oatp1a1 from Weaver, 2010
Vmax_apicalC = (Vmax_apical_invitro*PTC*protein*60*(MW/1e9)*1000) #!Vmax of basolateral transporters (mg/h/kg BW)
Vmax_apical = Vmax_apicalC * 0.3 #mg/h