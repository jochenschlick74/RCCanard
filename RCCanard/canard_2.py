#
# This file is part of the RCCanard distribution (https://github.com/jochenschlick74/RCCanard).
# Copyright (c) 2021 Jochen Schlick.
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

# Basic definitions and formulas following the book of Dieter Schall


from math import pi, asin, sqrt, cos, isclose
import numpy as np
from airfoils import profil, clarky10, ag25, ag35, goe795, mh32



class fluegel:
    def __init__(self,profil, b, li, lm, la, si, sa, x1, x2):
        self.profil=profil
        self.b, self.li, self.lm, self.la, self.si, self.sa, self.x1, self.x2= b, li, lm, la, si, sa, x1, x2
        self.doppeltrapez()
        self.foiltowing()
    
    def doppeltrapez(self):
        b, li, lm, la, si, sa, x1, x2= self.b, self.li, self.lm, self.la, self.si, self.sa, self.x1, self.x2
        # Fläche
        Fi=si*(li+lm)
        Fa=sa*(lm+la)
        F=Fi+Fa
        self.F=F

        # Streckung L
        self.L=b*b/F

        # Zuspitzung 
        lambda_i=lm/li
        lambda_a=la/lm
        self.lambdaF=(lm+la)/2/li

        # mittlere aerodynamische Flügeltiefe
        lmu_i=2/3*li*(1+lambda_i+lambda_i*lambda_i)/(1+lambda_i)
        lmu_a=2/3*lm*(1+lambda_a+lambda_a*lambda_a)/(1+lambda_a)
        self.lmu_eff=Fi/F*lmu_i + Fa/F*lmu_a

        #Neutralpunkt
        V1=1/3*si*x1*(li+2*lm)
        V2=sa*x1*(lm+la)
        V3=1/3*sa*(x2-x1)*(lm-2*la)
        self.xac=(V1+V2+V3)/F+0.25*self.lmu_eff

        # Zuspitzungswinkel 25% Linie
        self.phi25=np.rad2deg(asin((x2-li/4+la/4)/b/2))
        
    # liefert die Teilfläche eines Doppeltrapezes bei der Spannweite b1 aus
    def dp_teilflaeche(self,b1):
        b, li, lm, la, si, sa, x1, x2=self.b, self.li, self.lm, self.la, self.si, self.sa, self.x1, self.x2
        if b1>=b:
            Fi=si*(li+lm)
            Fa=sa*(lm+la)
            return Fi+Fa
        elif b1/2>si:
            lx=lm-((lm-la)/sa*(b1/2-si))
            Fi=si*(li+lm)
            Fa=(b1/2-si)*(lm+lx)
            return Fi+Fa
        else:
            lx=li-((li-lm)/si)*b1/2
            Fa=b1/2*(li+lx)
            return Fa

    # Block I Feld 3
    # camaxP: Max. Ca Wert Profil
    # cm0P: Nullmomentenbeiwert Profil
    # alpha0F: Nullauftriebswinkel Flügelprofil
    # alph0H: Nullauftriabswinkel HLW-Profil
    def foiltowing(self):
        caoptP, camaxP, cm0P, alpha0F, LF, phi25=self.profil.caoptP, self.profil.camaxP, self.profil.cm0P, self.profil.alpha0P, self.L, self.phi25 
        caalphaF=(pi*LF)/(1+sqrt(1+LF*LF/4))
        self.caalphaF=caalphaF
        self.camaxF=camaxP
        # caF=caalphaF*(alphaF-alpha0F)*0.0175
        cos_phi25=cos(np.deg2rad(phi25))

        # Nullmomentenbeiwert Flügel
        self.cm0F=(LF*cos_phi25*cos_phi25)/(LF+2*cos_phi25*cos_phi25)*cm0P

        # Einstellwinkel Flügel
        caF0=caoptP  # Betriebspunkt Flügel ist der optimale ca Wert des Profils
        self.caF0=caF0
        i= alpha0F + caF0/caalphaF*180.0/pi
        self.i=i

    def __str__(self):
        return f"F={self.F:.2f}, b={self.b:.2f}, lmu_eff={self.lmu_eff:.2f}, L={self.L:.2f}, phi25={self.phi25:.2f}, lambdaF={self.lambdaF:.3f}, caalphaF={self.caalphaF:.2f}, caF0={self.caF0}, xac={self.xac:.3f}, iF={self.i:.2f}"

def trapez(b, li, la, x1):
    return doppeltrapez(b=b, li=li, lm=li, la=la, si=0, sa=b/2, x1=0, x2=x1)



# class canard: definiert den Entenflügel und mache die dzu notwendigen Berechnungen
class canard:
    # Block II Feld 2
    # Fh: Flügelfläche Canard-flügel
    # lambda_h Zuspitzung Canard-flügel
    # epsilon: Abreiß-Winkel-Differenz in Grad
    # dlH: Profildicke (länge/höhe des Canard-Profils)
    # alphaH0: Nullauftriebswinkel Canard
    # alpha_max_P: Maximaler Anstellwinkel Canardprofil
    # camaxPMaximaler Ca des Canardprofile
    # GF: Flächenbelastung des gesamten Flugzeugs in N/m2
    # brH: Breite des Rumpfs
    def __init__(self,F,lambdaH,profil, GF, brH=0, hH=0):
        self.F=F
        FH=self.F
        self.profil=profil
        self.GF=GF
        self.brH=brH
        self.hH=hH
        self.lambdaH=lambdaH
        
        dlH=profil.dl
        self.la=8*dlH*sqrt(profil.camaxP/self.GF) # Achtung: Ergebnis in meter!!!
        self.li=self.la/self.lambdaH
        laH=self.la
        liH=self.li
        self.b=brH+2*self.F*lambdaH/self.la/(1+self.lambdaH)
        bH=self.b
        self.phi25=np.rad2deg(asin(((liH-laH)-liH/4+laH/4)/bH/2))
        self.L=bH*bH/FH
        self.lmu=2/3*liH*(1+lambdaH+lambdaH*lambdaH)/(1+lambdaH)
        V=1/3 *(bH-brH)/2*(liH-laH)*(liH+2*laH)
        self.xac=V/FH+0.25*self.lmu
        LH=self.L
        self.caalphaH=pi*LH/(1+sqrt(1+LH*LH/4))
        self.camaxH=self.profil.camaxP

    
    def __str__(self):
        return f"F={self.F:.2f}, b={self.b:.2f}, li={self.li:.2f}, la={self.la:.2f}, lmu={self.lmu:.3f}, L={self.L:.2f}, camaxH={self.camaxH:.2f}, xac={self.xac:.3f}, caalpha={self.caalphaH:.3f}"


# class entenflugzeug
# epsilon= Abreißdifferenz. Um wieviel grad früher reißt die Strömung am canard ab als am flügel
class entenflugzeug:
    def __init__(self,fluegel, canard, epsilon):
        self.fluegel=fluegel
        self.canard=canard
        self.epsilon=epsilon
        self.EWD=self.calc_EWD()
        self.betriebspunkthlw()
        
        kaw_new=0.7
        kaw_initial=0.0
        while abs(kaw_new-kaw_initial)/kaw_new>=0.05:
            kaw_initial=kaw_new
            self.lac=self.hebelarm(kaw=kaw_initial, SM=0.2) # erste grobe Abschätzung
            self.awg=self.abwindgradient()
            self.kaw=self.abwindfaktor()
            kaw_new=self.kaw

    # Block II Feld 4: EWD
    def calc_EWD(self):
        camaxF=self.fluegel.camaxF
        camaxH=self.canard.camaxH
        caalphaF=self.fluegel.caalphaF
        caalphaH=self.canard.caalphaH
        betaF=self.fluegel.profil.beta
        betaH=self.canard.profil.beta
        epsilon=self.epsilon
        return (camaxH/caalphaH - camaxF/caalphaF)*57.3 + betaH-betaF + epsilon
        
    # BlockII Feld 4
    def betriebspunkthlw(self):
        ewd=self.EWD
        self.iH=self.fluegel.i+ewd
        self.caH0=self.canard.caalphaH*(self.iH-self.canard.profil.alpha0P)*0.0175
    
    #todo: test Abwindfaktor!!
    # BlockIII Feld 1-4
    def abwindgradient(self):
        #liF=self.fluegel.li
        #bF=self.fluegel.b
        
        bH=self.canard.b
        hH=self.canard.hH
        lambdaH=self.canard.lambdaH
        LH=self.canard.L
        lac=self.lac
        phi25H=self.canard.phi25

        k1=1/LH-1/(1+pow(LH,1.7))
        k2=(10-3*lambdaH)/7
        k3=(1-abs(hH/bH))/pow(2*lac/bH,1/3)
        awg=4.44*pow(k1*k2*k3*sqrt(cos(phi25H)),1.19)
        return awg
    
    def abwindfaktor(self):
        caH0=self.caH0
        caalphaH=self.canard.caalphaH
        caF0=self.fluegel.caF0
        caalphaF=self.fluegel.caalphaF
        F=self.fluegel.F
        lambdaF=self.fluegel.lambdaF
        bH=self.canard.b
        liF=self.fluegel.li
        bF=self.fluegel.b
        awg=self.awg

        Fw=bH*liF*(1-bH/2/bF*(1-lambdaF))
        return 1-(((caH0/caF0)/(caalphaH/caalphaF))*(Fw/F)*awg)

    def hebelarm(self,kaw, SM=0.2):
        caH0=self.caH0
        caalphaH=self.canard.caalphaH
        FH=self.canard.F

        caF0=self.fluegel.caF0
        caalphaF=self.fluegel.caalphaF
        FF=self.fluegel.F
        
        cm0F=self.fluegel.cm0F
        lmu=self.fluegel.lmu_eff

        c=caH0/caF0
        g=caalphaH/caalphaF
        f=FH/FF
       
        lac=lmu*(kaw+f*g)/(c-g)/f*(SM*(1+c*f/kaw)-cm0F/caF0/kaw)
        return lac
    
    def __str__(self):
        return f"Flügel:\n   {self.fluegel.__str__()}\n" + \
            f"Canard:\n   {self.canard.__str__()}\n" + \
            f"EWD: {self.EWD:.2f}   " + \
            f"iH: {self.iH:.2f}   " +\
            f"caH0: {self.caH0:.2f}   " +\
            f"kaw: {self.kaw:.2f}   " +\
            f"lac: {self.lac:.3f}   "
            #            f"awg: {self.awg:.2f}   " + \
        





# End ------------------------------------------




def test_all():
    # Flügelform (einfaches Rechteck)
    b=2.0
    li=0.20
    lm=0.20
    la=0.20
    si=0.30
    sa=b/2-si
    x1=0
    x2=0
    # Profil Beispiel HQ2.5/9 mit 0 grad Klappe
    caoptP=0.6
    camaxP=1.1
    cm0P=-0.09
    alpha0F=-0.7
    caF0=caoptP # Betriebspunkt des Flügels ist der optimale ca des Profils

    F, L, phi25, lambda_eff, lmu_eff, xac = doppeltrapez(b=b, li=li, lm=lm, la=la, si=si, sa=sa, x1=x1, x2=x2)
    print(f"F={F:.2f}, L={L:.2f}, phi25={phi25:.2f}, lambda_eff={lambda_eff:.2f}, lmu_eff={lmu_eff:.2f}, xac={xac:.2f}")
    # Beispiel HQ2.5/9 mit 0 grad Klappe
    caaplhaF, cm0F, iF = foiltowing(caoptP=caoptP, camaxP=camaxP, cm0P=cm0P, alpha0F=alpha0F, LF=L, phi25=phi25)
    print(f"caalphaf={caaplhaF:.2f}, cmoF={cm0F:.2f}, iF={iF:.2f}")
    FH=F*0.3 # Blinde Definition Canard hat 30% der Flügelfläche
    #print(f"F={F}")
    lHmin,laH, liH, bH, LH, LmuH, xacH, phi25H= canard(FH=FH, lambda_H=1.0, epsilon=3, dlH=0.09, alpha0H=-0.7,alpha_max_P=10.0, camaxP=1.2, GF=10.0/F, brH=0)

    print(f"lHmin={lHmin:.2f},laH={laH:.2f}, liH={liH:.2f}, bH={bH:.2f}, LH={LH:.2f}, LmuH={LmuH:.2f}, xacH={xacH:.2f}")

    caalphaH, camaxH=foiltocanard(LH=LH, camaxF=1.1)
    alpha0H=-0.7
    EWD, iH, caH0=betriebspunkthlw(camaxF=camaxP, camaxH=camaxH, caalphaF=caaplhaF, caalphaH=caalphaH, betaF=3, betaH=3, epsilon=3, alpha0H=alpha0H, iF=iF)
    print(f"EWD={EWD:.2f}, iH={iH:.2f}, camaxH={camaxH:.2f}, caH0={caH0:.2f}")
    
    lac=hebelarm(caH0=caH0, caalphaH=caalphaH, FH=FH, caF0=caF0, caalphaF=caaplhaF, FF=F,kaw=1.0, cm0F=cm0F, caF=caF0, lmu=lmu_eff, SM=0.2)
    print(f"lac={lac}")

    # echter Abwindfaktor
    awg, Fw=abwindgradient(liF=li,bF=b,lambdaF=lambda_eff,bH=bH, hH=0.0, lambdaH=1.0,LH=LH,lac=lac,phi25H=phi25H)
    print(f"awg={awg}")
    #Fw=dp_teilflaeche(bH, b, li, lm, la, si, sa, x1, x2)
    kaw=abwindfaktor(caH0=caH0, caalphaH=caalphaH, caF0=caF0, caalphaF=caaplhaF, F=F, Fw=Fw, awg=awg)
    print(f"kaw={kaw}")


def test_classes():
    m=1.2 #kg
    f=fluegel(profil=mh32, b=1.8, li=0.22, lm=0.20, la=0.17, si=0.3, sa=0.9-0.3, x1=0.02, x2=0.05)
    print(f)
    c=canard(F=0.3*f.F, lambdaH=1.0, profil=clarky10, GF=m*9.81/1.3/f.F, brH=0.04, hH=-0.02)
    print(c)
    ef=entenflugzeug(fluegel=f, canard=c, epsilon=2)
    print(ef)




#test_doppeltrapez()
#test_foiltowing()
#test_canard()

#test_dp_teilflaeche();
#test_all()

test_classes()