
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

"""

Test data for the ncmatparser.py test. Kept here in a separate file as it is
very messy.

"""

__all__ = ['testdata']

import numpy as np

#Put test contents here. Encapsulate in tuple to indicate that we expect no errors, otherwise we expect a BadInput error.
validsimple = """NCMAT v1
@DEBYETEMPERATURE
10.0
@CELL
lengths 1 1 1.01
angles 90 90 90
@SPACEGROUP
1
@ATOMPOSITIONS
Al 0 0 0
"""

validsimplev1 = validsimple
validsimplev2 = validsimple.replace("NCMAT v1","NCMAT v2")
validsimplev3 = validsimple.replace("NCMAT v1","NCMAT v3")
validsimplev4 = validsimple.replace("NCMAT v1","NCMAT v4").replace("@DEBYETEMPERATURE\n10.0","@DEBYETEMPERATURE\nAl 10.0")
validsimplev5 = validsimplev4.replace("NCMAT v4","NCMAT v5")
validsimplev6 = validsimplev4.replace("NCMAT v4","NCMAT v6")
validsimplev7 = validsimplev4.replace("NCMAT v4","NCMAT v7")
danishchars = '\u00e6\u00f8\u00e5'
danish_aring = '\u00e5'


def gen_fake_dyninfo_scatknl_fields( *, ename, frac=1.0, temp, badgridval=False ):
    return f'''element {ename}
fraction {frac}
type     scatknl
temperature {temp}
alphagrid { " ".join(str(e) for e in np.linspace(1e-6,100.0,10)) if not badgridval else "bla blu" }
betagrid { " ".join(str(e) for e in np.linspace(1e-6,100.0,10)) }
sab_scaled { "1.0 "*100 }
'''

def simplecubic( vx, sg, length ):
    return f"""NCMAT {vx}
@DEBYETEMPERATURE
Al 10.0
@CELL
cubic {length}
@SPACEGROUP
{sg}
@ATOMPOSITIONS
Al 0 0 0
"""

def alheader( vx ):
    return f"""NCMAT {vx}
@CELL
 lengths 4.04958 4.04958 4.04958
 angles 90 90 90
@ATOMPOSITIONS
  Al 0   1/2 1/2
  Al 0     0   0
  Al 1/2 1/2   0
  Al 1/2   0 1/2
"""

def beoheader(vx):
    return f"""NCMAT {vx}
@CELL
  lengths 2.698 2.698 4.359
  angles 90 90 120
@SPACEGROUP
  186
@ATOMPOSITIONS
  Be 0   0   0
  Be 1/3 2/3 1/2
  O  0   0   0.378
  O  1/3 2/3 0.878
"""

some_vdos = """
  type vdos
  vdos_egrid .0048501456072871 .039822248144042
  vdos_density .013232 .0141193 .0150065 .0158938 .0167811 .0176684 .0185557
    .0195536 .020725 .0218964 .0230678 .0242392 .0254106 .0265808 .027748
    .509862 .484214 .458806 .434668 .410583 .386357 .362131 .337836 .313452
    .289068 .261524 .233749 .20274 .16795 .125961 .0770796 .0135243 0
"""

def otherphasescfgstr(cfgstr):
    return validsimplev6 +f"""
@OTHERPHASES
0.001 {cfgstr}
"""

def atomdbsec(marker):
    return f'@ATOMDB\n{marker} 10.0u 1fm 2b 5b'

testdata = [

    (validsimplev1,),
    (validsimplev2,),
    (validsimplev3,),
    (validsimplev4,),
    (validsimplev5,),
    (validsimplev6,),
    (validsimplev7,),
    (validsimplev3+atomdbsec('Al'),),
    validsimplev3+atomdbsec('X'),
    (validsimplev4+atomdbsec('Al'),),
    validsimplev4+atomdbsec('X'),
    validsimplev1.replace('Al','Bg'),
    validsimplev2.replace('Al','Bg'),
    validsimplev3.replace('Al','Bg'),
    validsimplev1.replace('Al','D'),
    (validsimplev2.replace('Al','D'),),
    (validsimplev3.replace('Al','D'),),
    validsimplev1.replace('Al','T'),
    validsimplev2.replace('Al','T'),
    (validsimplev3.replace('Al','T'),),
    validsimplev1.replace('Al','He3'),
    validsimplev2.replace('Al','He3'),
    (validsimplev3.replace('Al','He3'),),
    validsimplev1.replace('Al','He03'),
    validsimplev2.replace('Al','He03'),
    validsimplev3.replace('Al','He03'),
    validsimplev1.replace('Al','Bg3'),
    validsimplev2.replace('Al','Bg3'),
    validsimplev3.replace('Al','Bg3'),
    validsimplev1.replace('Al','X'),
    validsimplev2.replace('Al','X'),
    validsimplev3.replace('Al','X'),
    validsimplev1.replace('Al','X')+atomdbsec('X'),
    validsimplev2.replace('Al','X')+atomdbsec('X'),
    validsimplev3.replace('Al','X')+atomdbsec('X'),
    (validsimplev3.replace('Al','Be')+atomdbsec('Be'),),
    validsimplev1.replace('Al','X17'),
    validsimplev2.replace('Al','X17'),
    validsimplev3.replace('Al','X17'),
    validsimplev3.replace('Al','X17')+atomdbsec('X17'),
    (validsimplev3.replace('Al','Be17')+atomdbsec('Be17'),),
    validsimplev1.replace('Al','X99'),
    validsimplev2.replace('Al','X99'),
    validsimplev3.replace('Al','X99'),
    validsimplev3.replace('Al','X99')+atomdbsec('X99'),
    validsimplev3.replace('Al','Po45')+atomdbsec('Po45'),#Po has more than 45 protons...
    (validsimplev3.replace('Al','Xe99')+atomdbsec('Xe99'),),
    validsimplev1.replace('Al','X100'),
    validsimplev2.replace('Al','X100'),
    validsimplev3.replace('Al','X100'),
    validsimplev1.replace('Al','X117'),
    validsimplev2.replace('Al','X117'),
    validsimplev3.replace('Al','X117'),
    validsimplev1.replace('Al','X01'),
    validsimplev2.replace('Al','X01'),
    validsimplev3.replace('Al','X01'),
    validsimplev1.replace('Al','D')+atomdbsec('D'),
    validsimplev1.replace('Al','T')+atomdbsec('T'),
    validsimplev2.replace('Al','D')+atomdbsec('D'),
    validsimplev2.replace('Al','T')+atomdbsec('T'),
    (validsimplev3.replace('Al','D')+atomdbsec('D'),),
    (validsimplev3.replace('Al','T')+atomdbsec('T'),),
    (validsimplev3+atomdbsec('D'),),
    (validsimplev3+atomdbsec('T'),),
    validsimplev3+atomdbsec('J'),
"",

"""
""",

"""NCMAT v1
#blabla
""",

"""bla""",

"NCMAT v1#not ok first line",

"NCMAT v2#ok first line",

" NCMAT v1",

"NCMAT v1 bla",

"""NCMAT v17
""",

"""NCMAT v1
#some comment
bla bla
@CELL
bla""",

"""NCMAT v1
@CELL
lala 0.1 0.1 0.1
""",

"""NCMAT v2
@CELL
  #dumb stuff
  lengths 0.0 0.0 0.0
""",

"""NCMAT v1
@CELL
  lengths 0.1 0.2 0.2
""",

"""NCMAT v1
@CELL
  angles 0.1 -0.2 0.2
""",

"""NCMAT v2
@CELL#bla bla
#bla bla
  lengths 0 0 90#bla
  angles 0.1 -0.2 0.2
""",

validsimple.replace('NCMAT v1','NCMAT v1 #comment'),#comments not allowed except as full-lines in in header in v1
validsimple.replace('angles 90 90 90','angles 90 90 90 #comment'),#comments not allowed except in header in v1
validsimple.replace('angles 90 90 90','#comment\nangles 90 90 90'),#comments not allowed except in header in v1
validsimple.replace('NCMAT v1','NCMAT v1\n #comment'),#comments not allowed except as full-lines in in header in v1
validsimple.replace('@DEBYETEMPERATURE','@DEBYETEMPERATURE #comment'),#comments not allowed except as full-lines in in header in v1

"""NCMAT v1
@SomeName
bla
""",

"""NCMAT v1
@ATOMPOSITIONS
bla
""",

"""NCMAT v1
@ATOMPOSITIONS
He
""",

"""NCMAT v1
@ATOMPOSITIONS
0.0 0.0 1.0
""",

"""NCMAT v1
@ATOMPOSITIONS
Si 0.2 0.6 0.9
He 0.0 zero 1.0
""",

"""NCMAT v1
@ATOMPOSITIONS
  Si 0.0 0.5 1.0
""",

"""NCMAT v1
@ATOMPOSITIONS
@CELL
  lengths 1 2 3
  angles 4 5 6
""",

"""NCMAT v1
@ATOMPOSITIONS
  Si 0.1 0.6 0.9
  O 0 0 1.0
@SPACEGROUP
  230
@CELL
  lengths 1 2 3
  angles 4 5 6
""",

("""NCMAT v1
@ATOMPOSITIONS
  Si 0.1 0.6 0.9
  O 0 0 1.0
@SPACEGROUP
  230
@CELL
  lengths 1 2 3
  angles 40 50 60
@DEBYETEMPERATURE
  300.0
""",),

("""NCMAT v2
@ATOMPOSITIONS
  Si 0.1 0.6 0.9
  O 0 0 1.0
@SPACEGROUP
  230
@CELL
  lengths 1 2 3
  angles 40 50 60
@DEBYETEMPERATURE
  300.0
""",),

"""NCMAT v2
@SPACEGROUP
0
""",

"""NCMAT v1
@SPACEGROUP
-1
""",

"""NCMAT v1
@SPACEGROUP
pb + 3mbla
""",

"""NCMAT v1
@SPACEGROUP
pb+3mbla
""",

"""NCMAT v1
@SPACEGROUP
""",

"""NCMAT v1
@SPACEGROUP
@CELL
""",

"""NCMAT v1
@SPACEGROUP
1
2
""",

"""NCMAT v1
@SPACEGROUP
1
""",

"""NCMAT v1
@SPACEGROUP
225
""",

"""NCMAT v1
@SPACEGROUP
225
@SPACEGROUP
225
""",

"""NCMAT v1
@DEBYETEMPERATURE
200.0
200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
200.0
Si 200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
Si 200.0
O 200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
SI 200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
Sii 200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
c 200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
C 200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
Si 200.0
""",

"""NCMAT v1
@DEBYETEMPERATURE
Si 200.0
50.0
""",

f"""NCMAT v1
hello {danishchars}
""",

f"""NCMAT  v1#bla bla444
#blaf sd sdf {danishchars}
@DEBYETEMPERATURE
Si 200.0
O 100.0
@CELL
  \tlengths\t2.0 2.0\t\t\t 2.0   \t#bla\r
angles 90 90 90
@SPACEGROUP
225
@ATOMPOSITIONS\r\nSi 0.0 0.0 0.0
Si 0.5 0.5 0.5
O 0.5 0 0.5
""",

(f"""NCMAT  v2#bla bla444
#blaf sd sdf {danishchars}
@DEBYETEMPERATURE
Si 200.0
O 100.0
@CELL
  \tlengths\t2.0 2.0\t\t\t 2.0   \t#bla\r
angles 90 90 90
@SPACEGROUP
225
@ATOMPOSITIONS\r\nSi 0.0 0.0 0.0
Si 0.5 0.5 0.5
O 0.5 0 0.5
""",),

validsimplev2+"\n  \r  #bla",# bad because of the loose \r
validsimplev2+"#bla\r#bla",# here the loose \r is trying to hide inside a comment
(validsimplev2+"#bla\r\n#bla",),# here the \r is part of a DOS line ending
validsimplev2+"#bla\n\r#bla",# not a DOS line ending
(validsimple,),
validsimple+"""
@DYNINFO
element Hej
fraction 0.5
type sterile
""",#@DYNINFO not allowed in v1

validsimplev2+"""
@DYNINFO
element = Hej
fraction = 0.5
type = sterile
""",

validsimplev2+"""
@DYNINFO
element Hej
fraction 0.5
type sterile
""",

validsimplev2+"""
@DYNINFO
element Si
fraction 0.5
type sterile
""",

validsimplev2+"""
@DYNINFO
element Al
fraction 0.5
type sterile
""",

(validsimplev2+"""
@DYNINFO
element Al
fraction 1.0
type sterile
""",),

validsimplev2+"""
@DYNINFO
element Al
fraction 0.5
type sterile
@DYNINFO
element Al
fraction 0.5
type sterile
""",

validsimplev2+"""
@DYNINFO
""",

"""NCMAT v2
@DEBYETEMPERATURE
  10.0
@CELL
  lengths 1 1 1
  angles 90 90 90
@ATOMPOSITIONS
  Al 0 0 0
  Al 0 0.5 0
  O 0 0 0.5
@DYNINFO
  type sterile
  element Al
  fraction 0.6666667
@DYNINFO
  type freegas
  element O
  fraction 0.3333333
""",

("""NCMAT v2
@DEBYETEMPERATURE
  10.0
@CELL
  lengths 1 1 1
  angles 90 90 90
@ATOMPOSITIONS
  Al 0 0 0
  Al 0 0.5 0
  O 0 0 0.5
@DYNINFO
  type sterile
  element Al
  fraction 0.66666666666666666666666667
@DYNINFO
  type freegas
  element O
  fraction 0.33333333333333333333333333
""",),

("""NCMAT v2
@DEBYETEMPERATURE
  10.0
@CELL
  lengths 1 1 1
  angles 90 90 90
@ATOMPOSITIONS
  Al 0 0 0
  Al 0 0.5 0
  O 0 0 0.5
@DYNINFO
  type sterile
  element Al
  fraction 0.66666666666666666666666667
@DYNINFO
  type freegas
  element O
  fraction 0.33333333333333333333333333
""",),

##"""NCMAT v2
##@DENSITY
##  0.4 atoms/Aa^3
##""",
##
##"""NCMAT v2
##@DENSITY
##  0.4 atoms/ Aa^3
##""",
##
##"""NCMAT v2
##@DENSITY
##  0.4 kg/m^3
##""",
##
##"""NCMAT v2
##@DENSITY
##  0.4 g/cm^3
##""",
##


"""NCMAT v2
@DENSITY
  0.4 atoms_per_aa3""",

"""NCMAT v2
@DENSITY
  0.4 kg_per_m3""",

"""NCMAT v2
@DENSITY
  0.4 g_per_cm3""",

"""NCMAT v2
@DENSITY
  0.4 atoms_per_gallon""",

"""NCMAT v2
@DENSITY
""",

"""NCMAT v2
@DENSITY
  3.3455 kg per m3
""",

(validsimplev2.replace('Al ','D '),),#Deuterium marked with D ok in NCMAT v2
validsimple.replace('Al ','D '),#Deuterium not ok in NCMAT v1
validsimple.replace('Al ','H2 '),#Deuterium marked with H2 not ok in NCMAT v2
(validsimplev3.replace('Al ','D '),),#both D and H2 ok in NCMAT v3
(validsimplev3.replace('Al ','H2 '),),#both D and H2 ok in NCMAT v3
(validsimplev4.replace('Al ','D '),),#both D and H2 ok in NCMAT v4
(validsimplev4.replace('Al ','H2 '),),#both D and H2 ok in NCMAT v4
(validsimplev4.replace('Al ','T '),),#both T and H3 ok in NCMAT v4
(validsimplev4.replace('Al ','H3 '),),#both T and H3 ok in NCMAT v4

"""
NCMAT v3
@ATOMDB
  Cr mix 0.99 Cr 0.01 Ti
  Al mix 0.99 Al 0.01 Cr
  B10 10.0u 5fm 1b 100b
  X mix 1.0 B10
  X2 mix 1.0 D
  X3 mix 0.5 He 0.2 T 0.3 Cr
@CELL
 lengths 4.04958 4.04958 4.04958
 angles 90 90 90
@SPACEGROUP
  225
@ATOMPOSITIONS
  Al 0   1/2 1/2
  X 0     0   0
  X2 1/2 1/2   0
  X1 1/2   0 1/2
@DEBYETEMPERATURE
  Al   410.4
@DYNINFO
  element  Al
  fraction 1
  type     vdos
  vdos_egrid .0048501456072871 .039822248144042
  vdos_density .013232 .0141193 .0150065 .0158938 .0167811 .0176684 .0185557
    .0195536 .020725 .0218964 .0230678 .0242392 .0254106 .0265808 .027748
    .0289152 .0300825 .0312497 .032417 .0338129 .0352559 .0366988 .0381418
    .0395847 .0412955 .043043 .0447906 .0465382 .0482858 .0500848 .0519404
    .0537959 .0556515 .057507 .0593626 .0613858 .0638226 .0662593 .0686961
    .0711328 .0735696 .0760063 .0785759 .081366 .0841561 .0869462 .0897363
    .0925265 .0954792 .0987494 .10202 .10529 .10856 .11183 .115538 .119309
    .12308 .126851 .130875 .135268 .139661 .144055 .148448 .152842 .157573
    .162305 .167036 .171767 .176499 .18201 .187579 .193147 .198833 .20494
    .211048 .217155 .223262 .229369 .236467 .243837 .251208 .258578 .266432
    .275333 .284233 .293133 .303592 .314275 .324959 .335642 .349394 .363213
    .377789 .396045 .414301 .439971 .471089 .514605 .610025 .707931 .696698
    .687436 .679569 .672862 .667612 .662494 .65946 .656426 .653391 .651511
    .649818 .648548 .648291 .648033 .648114 .648959 .649803 .650648 .652271
    .654251 .656231 .658385 .660883 .66338 .666831 .670734 .674636 .678939
    .68327 .687962 .693084 .698207 .703766 .710271 .716775 .72328 .73082 .738646
    .746472 .755309 .764205 .773061 .781132 .789068 .794645 .797985 .798822
    .798529 .79508 .774768 .71273 .67965 .649939 .627278 .607527 .590861 .57517
    .561495 .54884 .537956 .527173 .516685 .506197 .496343 .487181 .478019
    .46898 .460042 .451104 .442408 .434216 .426025 .415941 .403984 .39053
    .375124 .372167 .371561 .37199 .372419 .372847 .373079 .373154 .373228
    .373214 .372782 .372351 .37192 .371412 .370855 .370298 .369741 .368973
    .367973 .366974 .365974 .364975 .363901 .362208 .360515 .358822 .356884
    .354458 .352031 .349604 .347178 .34379 .340173 .336555 .332841 .327735
    .322629 .317522 .311459 .304041 .296623 .287643 .278031 .260712 .247663
    .243349 .241651 .241493 .245625 .252347 .260515 .272178 .288183 .308977
    .337849 .366583 .394604 .428828 .466863 .508577 .555451 .617322 .679027
    .79657 .954686 .930006 .908781 .916995 .932774 .950737 .968525 .988039 1
    .979429 .895711 .774741 .718589 .672154 .631992 .596623 .565164 .537389
    .509862 .484214 .458806 .434668 .410583 .386357 .362131 .337836 .313452
    .289068 .261524 .233749 .20274 .16795 .125961 .0770796 .0135243 0
""",

    #Leave this here (the missing 'f' in front of the f-string exposed some parser bug previously):
(alheader('v3')+"""
@DEBYETEMPERATURE
  Al   410.4
@DYNINFO
  element  Al
  fraction 1
  {some_vdos}
"""),

#V3: must always specify debye temps even when have vdos:

(alheader('v3')+f"""
@DEBYETEMPERATURE
  Al   410.4
@DYNINFO
  element  Al
  fraction 1
  {some_vdos}
""",),


alheader('v3')+f"""
@DYNINFO
  element  Al
  fraction 1
  {some_vdos}
""",

#V4: debye temps optional when have vdos, and will give warning if provided (but not here where we just parse):

(alheader('v4')+f"""
@DEBYETEMPERATURE
  Al   400.1234567
@DYNINFO
  element  Al
  fraction 1
  {some_vdos}
""",),

(alheader('v4')+f"""
@DYNINFO
  element  Al
  fraction 1
  {some_vdos}
""",),

#empty @DEBYETEMPERATURE is not allowed:

alheader('v4')+f"""
@DEBYETEMPERATURE
@DYNINFO
  element  Al
  fraction 1
  {some_vdos}
""",
#global @DEBYETEMPERATURE is not allowed:

alheader('v4')+f"""
@DEBYETEMPERATURE
300
@DYNINFO
  element  Al
  fraction 1
  {some_vdos}
""",

#Ok:
(alheader('v4')+"""
@DEBYETEMPERATURE
  Al 300
@DYNINFO
  element  Al
  fraction 1
  type freegas
""",),

#ok, one has vdos and one has debye temp:
(beoheader('v4')+f"""
@DEBYETEMPERATURE
  Be 300
@DYNINFO
  element  O
  fraction 1/2
  {some_vdos}
@DYNINFO
  element  Be
  fraction 1/2
  type sterile
""",),

#ok, both have vdos:
(beoheader('v4')+f"""
@DYNINFO
  element  O
  fraction 1/2
  {some_vdos}
@DYNINFO
  element  Be
  fraction 1/2
  {some_vdos}
""",),

#ok, both have vdos and one has debye temp (should emit warning if loaded not just parsed):
(beoheader('v4')+f"""
@DYNINFO
  element  O
  fraction 1/2
  {some_vdos}
@DYNINFO
  element  Be
  fraction 1/2
  {some_vdos}
@DEBYETEMPERATURE
  Be 150.0
""",),

#not ok, one has vdos and debye temp, the other has nothing:
beoheader('v4')+f"""
@DEBYETEMPERATURE
  Be 300
@DYNINFO
  element  O
  fraction 1/2
  type sterile
@DYNINFO
  element  Be
  fraction 1/2
  {some_vdos}
""",

#But ok in v5 with vdosdebye:
(beoheader('v5')+f"""
@DYNINFO
  element  O
  fraction 1/2
  type vdosdebye
  debye_temp 300
@DYNINFO
  element  Be
  fraction 1/2
  {some_vdos}
""",),

#Not ok in v5 with sterile/freegas:
beoheader('v5')+f"""
@DYNINFO
  element  O
  fraction 1/2
  type freegas
@DYNINFO
  element  Be
  fraction 1/2
  {some_vdos}
""",

#cubic keyword ok in v4:
(simplecubic('v4',225,4.01),),
(simplecubic('v4',195,4.01),),
(simplecubic('v4',230,4.01),),
#Not before
simplecubic('v1',225,4.01),
simplecubic('v2',225,4.01),
simplecubic('v3',225,4.01),
#Not if we mis-spell cubic:
simplecubic('v4',225,4.01).replace('cubic ','Cubic '),
#v4, but wrong somehow:
simplecubic('v4',194,4.01),
simplecubic('v4',1,4.01),
simplecubic('v4',225,4.01e10),
simplecubic('v4',225,0.0),
simplecubic('v4',225,4.01).replace('@CELL\n','@CELL\nlengths 4.0 4.0 4.0\n'),

("""NCMAT v4
@CELL
lengths 2.2 !! 3.3
angles 90 90 120
@ATOMPOSITIONS
Al 0 0 0
@DEBYETEMPERATURE
Al 10.0
""",),

("""NCMAT v4
@CELL
lengths 2.2 !! !!
angles 90 90 120
@ATOMPOSITIONS
Al 0 0 0
@DEBYETEMPERATURE
Al 10.0
""",),

"""NCMAT v4
@CELL
lengths !! 2.2 3.3
angles 90 90 120
@ATOMPOSITIONS
Al 0 0 0
@DEBYETEMPERATURE
Al 10.0
""",

"""NCMAT v4
@CELL
lengths 2.2 2.2 3.3
angles 90 !! 120
@ATOMPOSITIONS
Al 0 0 0
@DEBYETEMPERATURE
Al 10.0
""",

"""NCMAT v3
@CELL
lengths 2.2 !! 3.3
angles 90 90 120
@ATOMPOSITIONS
Al 0 0 0
@DEBYETEMPERATURE
Al 10.0
""",

]

vdtest="""NCMAT v5
@CELL
  lengths 4 4 4
  angles 90 90 90
@ATOMPOSITIONS
  H 0 0 0
@DYNINFO
  element H
  fraction 1
  type vdosdebye
  debye_temp 200.0
"""

#amorphus material, not all elements have MSD available:
vdtest_nonxtal="""NCMAT v5
@DENSITY
  1.2 g_per_cm3
@DYNINFO
  element H
  fraction 2/3
  type vdosdebye
  debye_temp 200.0
@DYNINFO
  element C
  fraction 1/3
  type freegas
"""

testdata += [
    (vdtest,),
    vdtest.replace('NCMAT v5','NCMAT v4'),
    vdtest.replace('NCMAT v5','NCMAT v3'),
    vdtest.replace('NCMAT v5','NCMAT v2'),
    vdtest.replace('NCMAT v5','NCMAT v1'),
    (vdtest+"  debye_temp 300.0"),
    vdtest.replace('debye_temp 200.0','debye_temp sdfsf'),
    vdtest.replace('debye_temp 200.0','debye_temp 200.0 K'),
    vdtest.replace('debye_temp 200.0','debye_temp 200.0K'),
    vdtest.replace('debye_temp 200.0','debye_temp 200.0  300'),
    (vdtest+"\n@DEBYETEMPERATURE\n  H 200.0\n"),
    (vdtest_nonxtal,),
    vdtest_nonxtal.replace('debye_temp 200.0',''),
]

testdata += [
#error on 0 density:
"""NCMAT v5
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  0 g_per_cm3
""",
    #errors on repeated lines:
"""NCMAT v5
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
  1 g_per_cm3
""",
"""NCMAT v5
@DEBYETEMPERATURE
 Al 10.0
@CELL
 cubic 4
@SPACEGROUP
 225
 225
@ATOMPOSITIONS
  Al 0 0 0
""",
"""NCMAT v5
@DYNINFO
  fraction 1
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@DYNINFO
  fraction 1
  element H
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@DYNINFO
  fraction 1
  element H
  type freegas
  type freegas
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@DYNINFO
  fraction 1
  element H
  type vdosdebye
  debye_temp 200
  debye_temp 200
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@DYNINFO
  fraction 1
  element H
  type freegas
  debye_temp 200
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@DYNINFO
  fraction 1
  element H
  type vdos
  vdos_egrid .005 .04
  vdos_density .013232 .0141193 .0150065 .0158938 .0167811 .0176684 .0185557
  vdos_density .013232 .0141193 .0150065 .0158938 .0167811 .0176684 .0185557
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@DYNINFO
  fraction 1
  element H
  type vdos
  vdos_egrid .005 .04
  vdos_egrid .005 .04
  vdos_density .013232 .0141193 .0150065 .0158938 .0167811 .0176684 .0185557
@DENSITY
  1 g_per_cm3
""",
#state of matter
("""NCMAT v5
@STATEOFMATTER
  solid
@DYNINFO
  fraction 1
  element H
  type vdos
  vdos_egrid .005 .04
  vdos_density .013232 .0141193 .0150065 .0158938 .0167811 .0176684 .0185557
@DENSITY
  1 g_per_cm3
""",),
"""NCMAT v5
@STATEOFMATTER
  gas
@DYNINFO
  fraction 1
  element H
  type vdos
  vdos_egrid .005 .04
  vdos_density .013232 .0141193 .0150065 .0158938 .0167811 .0176684 .0185557
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@STATEOFMATTER
  liquid
@DYNINFO
  fraction 1
  element H
  type vdos
  vdos_egrid .005 .04
  vdos_density .013232 .0141193 .0150065 .0158938 .0167811 .0176684 .0185557
@DENSITY
  1 g_per_cm3
""",
("""NCMAT v5
@STATEOFMATTER
  liquid
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",),
("""NCMAT v5
@STATEOFMATTER
  gas
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",),
("""NCMAT v5
@STATEOFMATTER
  solid
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",),
"""NCMAT v5
@STATEOFMATTER
  blabla
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v4
@STATEOFMATTER
  solid
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@STATEOFMATTER
  solid
  solid
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@STATEOFMATTER
  solid
  gas
@DYNINFO
  fraction 1
  element H
  type freegas
@DENSITY
  1 g_per_cm3
""",
"""NCMAT v5
@STATEOFMATTER
  solid
@DYNINFO
  fraction 1
  element H
  type freegas
@STATEOFMATTER
  solid
@DENSITY
  1 g_per_cm3
""",

(validsimplev6 +"""
@OTHERPHASES
0.001 void.ncmat
""",),
validsimplev6 +"""
@OTHERPHASES
0.001 void.ncmat
@OTHERPHASES
0.001 void.ncmat
""",

#UTF8 chars never allowed in ncmat data, even if normally allowed in filenames in cfg strings:
otherphasescfgstr(f'utf8_filen{danish_aring}me.ncmat'),
otherphasescfgstr('utf8_filename.ncmat ; \n dcutoff = 0.4'),#no newlines:
(otherphasescfgstr('f.ncmat;dcutoff=0.4'),),
(otherphasescfgstr('f.ncmat ; dcutoff =   0.4'),),
(otherphasescfgstr('f.ncmat dcutoff=0.8'),),#will not fail during parsing but later
(otherphasescfgstr('f.ncmat;dcutoff=0.4 dcutoffup=2.0'),),#will not fail during parsing but later
(otherphasescfgstr('f.ncmat;atomdb=H   is   D'),),
(otherphasescfgstr('  f.ncmat \t \t ; \t\t\tatomdb\t\t = H:: \t\t:  is:  :: D  '),),
(otherphasescfgstr('phases<0.4*f.ncmat;dcutoffup=2.0&0.6*g.ncmat> ;dcutoff=0.4'),),
(otherphasescfgstr('phases<0.4*f.ncmat;dcutoffup=2.0&0.6*g.ncmat> dcutoff=0.4'),),#will not fail during parsing but later
(otherphasescfgstr('phases<0.4*f.ncmat&0.6*g.ncmat>'),),
(otherphasescfgstr(' phases   <\t\t0.4 *   f.ncmat\t\t&0.6*g.ncmat >'),),
otherphasescfgstr(''),

(validsimplev7 +"\n@TEMPERATURE\n500\n",),
(validsimplev7 +"\n@TEMPERATURE\ndefault 500\n",),
(validsimplev7 +"\n@TEMPERATURE\n1e-99\n",),
(validsimplev7 +"\n@TEMPERATURE\n1e6\n",),
validsimplev7 +"\n@TEMPERATURE\n1.0000000001e6\n",
validsimplev7 +"\n@TEMPERATURE\n0.0\n",
validsimplev7 +"\n@TEMPERATURE\n-400\n",
(validsimplev7+'\n@DYNINFO\n'+gen_fake_dyninfo_scatknl_fields(ename='Al',temp=300)+"\n@TEMPERATURE\n300\n",),
validsimplev7+'\n@DYNINFO\n'+gen_fake_dyninfo_scatknl_fields(ename='Al',temp=300,badgridval=True)+"\n@TEMPERATURE\n300\n",
validsimplev7+'\n@DYNINFO\n'+gen_fake_dyninfo_scatknl_fields(ename='Al',temp=300)+"\n@TEMPERATURE\n400\n",

]

#should support only \r\n and \n, disallow \r unless it is immediately followed by \n (in which case it is ignore)
