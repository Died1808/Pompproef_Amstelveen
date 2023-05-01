"""
Onderstaande modules zijn nodig voor deze berekeningen en kunnen via pip install
aan Python worden toegevoegd
""" 

import math
import numpy as np
import pandas as pd
from ttim import *
import matplotlib.pyplot as plt
import fiona
import os
from cycler import cycler
pd.set_option('display.max_columns', None) 
pd.options.display.width=None
pd.set_option('display.max_rows',10) 
pd.set_option('display.expand_frame_repr', True)
pd.options.display.float_format = '{:,.2f}'.format

"""
Onderstaande def definieert de gebruikte kleuren en layout van de figuren
"""

def kleuren(n):
    cm = plt.get_cmap('gist_rainbow')
    color=[cm(1.*i/n) for i in range(n)]
    params = {'font.family': 'sans-serif',
              'font.sans-serif': 'arial',
              'axes.labelsize': 10,
              'axes.facecolor': '#ffffff', 
              'axes.labelcolor': 'black',
              'axes.prop_cycle': cycler('color', color),
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'lines.linewidth': 1,
              'grid.color': 'grey',
              'grid.linestyle': 'dashed',
              'grid.linewidth': 0.5,
              'text.usetex': False,
              'font.style': 'normal',
              'font.variant':'normal',
              'figure.facecolor': 'white',
              'font.size':8,
              'figure.autolayout': True,
              'figure.figsize': (8,6),
              'figure.dpi': 100,
              }
    plt.rcParams.update(params)

"""
Onderstaand blok is deels gebaseerd op een publicatie van Rob van Putten op Linkedin.
Het is zodanig aangepast en uitegbreid dat ze een willekeurige gef file (die in dezelfde
directory staat als dit script) kan lezen en omvormen tot een invoer dataframe voor TTIM.
Hierin worden de doorlatendheid en de anisotropie bepaald aan de hand van zelf ontwikelde
rekenregels. 
Vanzelfsprekend kan hier ook een eigen bodemindeling worden opgenomen
voor de aansturing van TTIM, het gegeven dat met de automatisch gegenereerde parameters
een goed fit wordt berekend wil niet zeggen dat ze ook een weergave van de werkelijkheid is.
Aan de andere kant is een geautomatiseerde invoer wel prettig om persoonlijke voorkueren een
beetje te mijden.
"""
class GEF:
    def __init__(self):
        self._data_seperator = ' '
        self._columns = {}
        self.x = 0.
        self.y = 0.
        self.z = 0.
        self.dz = []
        self.qc = []
        self.pw = []
        self.wg = []
        self.c  = []
        self.kb = []
        self.dist =[]
        
    def readFile(self, filename):
        lines = open(filename, 'r').readlines()
        for line in lines:
            reading_header = True
        for line in lines:   
            if reading_header:
                self._parseHeaderLine(line)
            else:
                self._parseDataLine(line)
            if line.find('#EOH') > -1:
                if self._check_header():
                    reading_header = False
                else:
                    print(filename,'bestaat al')
                    return
            
    def _check_header(self):
        if not 1 in self._columns:
            return False
        if not 2 in self._columns:
            return False
        return True

    def _parseHeaderLine(self, line):
        for xe in ['#COMMENT', 'Peil=', 'uitvoerder', 'materieel','WATERSTAND',
                    'opmerkingen#MEASUREMENTTEXT','==','= NAP','antropogeen']:
            if xe in line:
                return          
        if len(line.split()) == 0:
            return
        
        keyword, argline = line.split('=')         
        keyword = keyword.strip()
        argline = argline.strip()
        args = argline.split(',')
      
        if '#XYID' in line:
            argline = argline.replace('.','')        
        args = argline.split(',')

        if keyword=='#XYID':
            if float(args[1]) < 1e5:
                args[1] = args[1]
            else:
                args[1]=args[1].replace('.','')
            self.x = float(args[1])
            args[2]=args[2].replace('.','')
            self.y = float(args[2])
            if (len(str(int(self.x))))>5:
                 self.x=int(self.x/pow(10,len(str(int(self.x)))-6))
            if (len(str(int(self.y))))>5:
                 self.y=int(self.y/pow(10,len(str(int(self.y)))-6))
            if self.x > 3e5:
                self.x=self.x/10

        elif keyword=='#ZID':
            self.z = round(float(args[1]),3)
           
        elif keyword=='#COLUMNINFO':
            column = int(args[0])
            dtype = int(args[-1])
            if dtype==11:
                dtype = 10
            self._columns[dtype] = column - 1    
       
    def _parseDataLine(self, line):
        line=line.strip()
        line = line.replace('|',' ')
        line = line.replace(';',' ')
        line = line.replace('!',' ')
        line = line.replace(':',' ')
        args=line.split()
        for n, i in enumerate(args):
            if i in ['9.9990e+003','-9999.9900','-1.0000e+007','-99999','-99999.0',
                  '-99999.00','-99999.000','-9.9990e+003','999.000', '-9999.99', '999.999',
                  '99','9.999','99.999', '999.9']:
                args[n] = '0.1'       
        if len(line.split()) == 0:
            return
 
        zz  = round(abs(float(args[self._columns[1]])),4)
        dz = round(self.z - zz,4)
        qc = round(float(args[self._columns[2]]),4)  
        pw = float(args[self._columns[3]]) 
        if pw<-10:
            pw=0.1

        self.dz.append(dz)
        self.qc.append(qc)
        self.pw.append(pw)
        if qc<=0.001:
            qc=0.1
            self.wg.append(10.)
        else:
            wg = abs((pw / qc) * 100.)
        wg = abs((pw / qc) * 100.)

##############################################K-waarde        
        if wg>=0.0:
            if wg >5: wg=15
            ke=math.exp(wg)
        if ke <=0:  ke=1
        else:
            kb  = (qc / ke)*2.6
            self.kb.append(kb)
        
    """
    Onderstaande blokken verzamelen de data series en vormen deze om tot een dataframe
    met diepte (tov NAP) en doorlatendheid (m/dag)
    """
            
    def asNumpy(self):
        return np.transpose(np.array([self.dz, self.kb]))

    def asDataFrame(self):
        a = self.asNumpy()
        return pd.DataFrame(data=a, columns=['depth', 'k'])
    
    """
    Onderstaand blok bereidt het dataframe voor op invoer in TTIM. 
    Omdat TTIM bij zeer veel lagen te veel geheugen vraagt wordt de laagopbouw 
    uit de sondering (per 2 cm) omgevormd tot een laagopbouw per 0.2 meter).
    In het dan ontstane dataframe wordt de anisotropie op basis van een rekenregel
    bijgevoegd. De onderzijde van de aquifer wordt op basis van omgevingsinformatie
    (o.a. DINOLOKet) als harde waarde ingegeven (dzend). 
    De dikte van de aquifer en de Ss worden hard ingegeven (bepaald aan de hand van trial
    and error op basis van een verlaging/tijd grafiek in de waarnemingspeilbuis. Uit de dikte
    en de Ss volgt de totale berging (S) van de aquifer. 
    Het onttrokken debiet {m3/dag} wordt hard ingegeven

    """                                                     
    
    def df_TTIM(self, filename):
        df = self.asDataFrame()
        df = df.sort_values('depth', ascending=False)

        if df.empty:
            return df
        
        df = df.rolling(10).mean() 
        df = df.iloc[:: 10]
        df=df.dropna()
        df = df.reset_index(drop=True)
        
        df['kzkh']= 2/(33.653*np.exp(-0.066*(df['k'])))
        df['kzkh']= np.where(df['kzkh']>1, 1, df['kzkh'])

        print(df)

        dzend=-57
        print('Aangenomen diepte aq = ', dzend, '[m NAP]')
        dfd=df['depth']
        dfd.loc['ld']= dzend
        dfn = df.iloc[: , [1]].copy()   
        dfn = pd.concat([df, pd.DataFrame.from_records([{'depth':dzend},])], ignore_index=True)

        Dikte_aq = 46
        Ss  =  7e-5
        Saq =  Ss*Dikte_aq  
        print('Ss  = ' , '{:.2e}'.format(Ss), ' [-]')
        print('Saq = ' , '{:.2e}'.format(Saq), ' [-]')
        Qm     = 1195 # Debiet per dag
        
        """
        Hieronder wordt het TTIM model gemaakt De pompput met filterstelling
        wordt opgenomne, het script bepaalt dan zelf in welek lagen de filters in het
        model komen. 
        NB: Dit is een voorbeeld van hoe er met TTIM gerekend kan worden, nogmaals,
        de correlatie met de werkelijkheid is vooralsnog toevallig.
        Tenslotte wordt het model doorgerekend (duurt circa 5 minuten)
        """

        ml=Model3D(kaq=df['k'], z=dfn[ 'depth'], Saq= Ss, kzoverkh=df['kzkh'], tmin=0.0001, 
                    tmax=10, M=10)
        #pompput
        tf = -14 # top Filter
        bf = -29 # onderzijde filter
        tfind = df['depth'].sub(tf).abs().values.argmin()
        bfind = df['depth'].sub(bf).abs().values.argmin() 
        well1 = np.arange(tfind, bfind, 1)
        w1 = Well(ml, xw=0, yw=0, rw=0.25, tsandQ=[(0,Qm)], layers=well1)

        ml.solve()
        
        """
        Hieronder worden de uitkomsten getekend voor de verschillende lagen.
        Verder wordt de resulterende KD voor deze aquifer berekend.
        De bovenste figuur geeft een indicatie van de druklijn in de ondergrond, de
        onderste het verloop in de tijd van de verlagingen in de peilbuis op
        23 meter. Deze onderste figuur wordt sterk beinvloedt door de S.
        De meetgevens zijn opgenomen in een externe Excel (Pb-PP.xlsx)
      
        """
        laag = np.arange(61,62,1) #Zie vera de jong pp 118/158
        kleuren(len(laag))
        
        ##############  Gemeten peilbuizen
        mat=pd.read_excel('Pb_PP.xlsx', index_col=None, skiprows=0, engine='openpyxl')  
        x=mat['tijd']/1440
        y1=mat['P5']
        y2=mat['P10']
        y3=mat['P20']
        y4=mat['P40']
        y5=mat['P80']

        for la in laag:
            ml.xsection(x1=0, x2=1000, y1=0, y2=0, npoints=1000, t=5.25, layers=[la], lw=3, newfig = False)

        leg=plt.legend(round(df.loc[df.index[laag], 'depth'],1), loc=(0.85,0.05), fontsize = 8, title='Filter op [m NAP]')
        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)

        plt.grid(axis='both')
        plt.semilogx()
        plt.xlim(1,1000)
        plt.ylim(-3,0)
        plt.ylabel('Verlaging [m]')
        plt.xlabel('Afstand tot bron [m]')
        if dzend < df.iloc[-1,0]:
            extra = 0
        else:
            extra = abs(dzend)-abs(df.iloc[-1,0])
        KD = round((df['k'].sum()/5+extra*df.iloc[-1,1]),0)       
        print(KD)
        title = (filename +' PP_Amstelveen'+ '   KD = '+ str(int(KD))+ ' [m2/dag]')
        plt.title(title, fontsize=12, pad=10)
        
        a=w1.headinside(5.25)
 
        plt.plot(4.5,-1.89,'k+', markersize=20)
        plt.plot(9.5,-1.44,'k+', markersize=20)
        plt.plot(19.5,-1.11,'k+', markersize=20)
        plt.plot(34,-0.87,'k+', markersize=20)
        plt.plot(50,-0.76,'k+', markersize=20)
        plt.plot(71,-0.64,'k+', markersize=20)
        plt.plot(138,-0.42,'k+', markersize=20)
        plt.plot(253,-0.31,'k+', markersize=20)
        plt.plot(409.5,-0.24,'k+', markersize=20)
        plt.plot(915,-0.1,'k+', markersize=20)
        plt.savefig('Afstand.png', bbox_inches='tight')
        plt.show()
        
        pb= [4,9,18,33,50,71,139,253,409, 915]
        kleuren(len(pb))
        
    
        for p_b in pb:
            t=np.logspace(-4,1,1000)
            s1=ml.head(p_b,0,t)
            plt.semilogx(t,s1[61]*-1, lw=3, alpha=0.5)
           

        plt.grid(axis='both')
        plt.plot(x,y1, 'o', color='darkred',    markersize=10)
        plt.plot(x,y2, 'o', color='darkgoldenrod', markersize=10)
        plt.plot(x,y3, 'o', color='gold', markersize=10)
        plt.plot(x,y4, 'o', color='lightgreen', markersize=10)
        plt.plot(x,y5, 'o', color='darkgreen', markersize=10)
       
        
        
        # plt.plot(x,y2, '+', color='black', markersize=3)
        plt.ylim(2,-0.5)
        plt.xlabel('tijd [dagen]')
        plt.ylabel('verlaging [m]')
        plt.savefig('Tijd_dh.png', bbox_inches='tight')
        plt.show()
        plt.close()
        
for filename in os.listdir(os.getcwd()):
    if filename.endswith ('.GEF') or filename.endswith ('.gef'):
        if __name__=="__main__":
            g=GEF()
            g.readFile(filename)
            g.df_TTIM(filename)
