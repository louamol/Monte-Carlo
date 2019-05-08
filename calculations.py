import numpy as np
from sympy import *

init_printing()

p11, p12, p22, p13, p23, p33, p14, p24, p34, p44 = symbols('φ_11 φ_12 φ_22 φ_13 φ_23 φ_33 φ_14 φ_24 φ_34 φ_44')
t12, t13, t23, t14, t24, t34 = symbols('θ_12 θ_13 θ_23 θ_14 θ_24 θ_34')

phi1 =exp(I*p11)

m1 =  Matrix([[exp(I*p12)*cos(t12),exp(I*p12)*sin(t12)],[-exp(I*p22)*sin(t12),exp(I*p22)*cos(t12)]])

phi2 = m1*Matrix([[phi1,0],[0,1]])

m2 = Matrix([[exp(I*p13)*cos(t13),exp(I*p13)*sin(t13)*cos(t23),exp(I*p13)*sin(t13)*sin(t23)],[-exp(I*p23)*sin(t13),exp(I*p23)*cos(t13)*cos(t23),exp(I*p23)*cos(t13)*sin(t23)],[0,-exp(I*p33)*sin(t23),exp(I*p33)*cos(t23)]])

phi3 = m2*Matrix([[phi2.row(0)[0],phi2.row(0)[1],0],[phi2.row(1)[0],phi2.row(1)[1],0],[0,0,1]])

m3 = Matrix([[exp(I*p14)*cos(t14),exp(I*p14)*sin(t14)*cos(t24),exp(I*p14)*sin(t14)*sin(t24)*cos(t34),exp(I*p14)*sin(t14)*sin(t24)*sin(t34)],[-exp(I*p24)*sin(t14),exp(I*p24)*cos(t14)*cos(t24),exp(I*p24)*cos(t14)*sin(t24)*cos(t34),exp(I*p24)*cos(t14)*sin(t24)*sin(t34)],[0,-exp(I*p34)*sin(t24),exp(I*p34)*cos(t24)*cos(t34),exp(I*p34)*cos(t24)*sin(t34)],[0,0,exp(I*p44)*cos(t34),exp(I*p44)*sin(t34)]])

phi4 = m3*Matrix([[phi3.row(0)[0],phi3.row(0)[1],phi3.row(0)[2],0],[phi3.row(1)[0],phi3.row(1)[1],phi3.row(1)[2],0],[phi3.row(2)[0],phi3.row(2)[1],phi3.row(2)[2],0],[0,0,0,1]])
pprint(phi4)

def com(A,B):
    return A*B-B*A

def acom(A,B):
    return A*B+B*A

def proba(U):
    U00 = Matrix([[U.row(0)[0],U.row(0)[1]],[U.row(1)[0],U.row(1)[1]]])
    U01 = Matrix([[U.row(0)[2],U.row(0)[3]],[U.row(1)[2],U.row(1)[3]]])
    U10 = Matrix([[U.row(2)[0],U.row(2)[1]],[U.row(3)[0],U.row(3)[1]]])
    U11 = Matrix([[U.row(3)[2],U.row(3)[3]],[U.row(3)[2],U.row(3)[3]]])
    print(U00,U01,U10,U11)

    a = com(U01,U00)**2
    b = 1/2*acom(com(U11,U00),com(U01,U00)) + 1/2*acom(com(U10,U01),com(U00,U01))
    c = 1/4*com(U10,U01)**2+1/4*com(U11,U00)**2 + 1/4*acom(com(U10,U01),com(U00,U11))
    d = 1/sqrt(2)*acom(com(U01,U00),com(U11,U10))
    e = 1/2*acom(com(U11,U10),com(U11,U00))+1/2*acom(com(U11,U10),com(U01,U10))
    f = com(U11,U10)**2

    norm = 1

    return 1/4*(a**2+b**2+c**2+d**2+e**2+f**2)/norm

U = 1/sqrt(2)*Matrix([[0,0,1,1],[0,0,1,-1],[1,-1,0,0],[-1,-1,0,0]])

pprint(proba(U))

pprint(proba(phi4))
