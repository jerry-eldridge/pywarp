from sympy import Matrix, Symbol,zeros,sin,cos, simplify

def Christoffel1_s(g,c,a,b):
    return simplify(0.5*(g[c,a].diff(Symbol(x[b])) +
                g[c,b].diff(Symbol(x[a])) +
                g[a,b].diff(Symbol(x[c]))).expand())
def Christoffel1(g,c,a,b):
    global christoffel1
    try:
        return christoffel1[c,a,b]
    except:
        christoffel1[c,a,b] = Christoffel1_s(g,c,a,b)
        return christoffel1[c,a,b]
def Christoffel2_s(g,i,k,l):
    ginv = g.inv()
    s = 0
    n = g.shape[0]
    for m in range(n):
        s = s + 0.5*ginv[i,m]*Christoffel1(g,m,k,l)
    return simplify(s.expand())
def Christoffel2(g,i,k,l):
    global christoffel2
    try:
        return christoffel2[i,k,l]
    except:
        christoffel2[i,k,l] = Christoffel2_s(g,i,k,l)
        return christoffel2[i,k,l]
def RiemannCurvature_s(g,rho,sig,mu,nu):
    val1 = Christoffel2(g,rho,nu,sig).diff(Symbol(x[mu])) - \
           Christoffel2(g,rho,mu,sig).diff(Symbol(x[nu]))
    val2 = 0
    val3 = 0
    n = g.shape[0]
    for lam in range(n):
        val2 = val2 + Christoffel2(g,rho,mu,lam)*Christoffel2(
            g,lam,nu,sig)
        val3 = val3 + Christoffel2(g,rho,nu,lam)*Christoffel2(
            g,lam,mu,sig)
    return simplify((val1 + val2 - val3).expand())
def RiemannCurvature(g,rho,sig,mu,nu):
    global riemanncurvature
    try:
        return riemanncurvature[rho,sig,mu,nu]
    except:
        riemanncurvature[rho,sig,mu,nu] = RiemannCurvature_s(g,rho,sig,mu,nu)
        return riemanncurvature[rho,sig,mu,nu]
        
def RicciTensor_s(g,i,j):
    s = 0
    n = g.shape[0]
    for k in range(n):
        s = s + RiemannCurvature(g,k,i,k,j)
    return simplify(s.expand())
def RicciTensor(g,i,j):
    global riccitensor
    try:
        return riccitensor[i,j]
    except:
        riccitensor[i,j] = RicciTensor_s(g,i,j)
        return riccitensor[i,j]

def ScalarCurvature(g):
    ginv = g.inv()
    s = 0
    n = g.shape[0]
    for i in range(n):
        for j in range(n):
            s = s + ginv[i,j]*RicciTensor(g,i,j)
    return s

def Base(i,base, bits):
    L = []
    j = 0
    n = 0
    ii = i
    for j in range(bits):
        a = ii%base
        ii = int(ii/base)
        n = n + a*base**j
        L.append(a)
    L.reverse()
    return L

def ShowPolar(reset=True):
    global x,X
    global christoffel1, christoffel2, riemanncurvature, riccitensor
    if reset:
        christoffel1 = {}
        christoffel2 = {}
        riemanncurvature = {}
        riccitensor = {}
    x = ["r","theta"]
    n = len(x)
    [r,theta] = map(Symbol,x)
    X = [r*cos(theta), r*sin(theta)]
    N = len(X)

    print "\nx = ", x
    print "\nX = ", X

    J = zeros( (len(x), len(X)))
    for i in range(N):
        for j in range(n):
            J[i,j] = X[i].diff(Symbol(x[j]))
    J = Matrix(J)
    print "\n(Jacobian) J=\n",J
    g = (J.T)*(J)
    for i in range(n):
        for j in range(n):
            g[i,j] = simplify(g[i,j].expand())
    print "\n(Metric Tensor) g_ij=\n",g
    print "\n(Inverse Metric Tensor) g_IJ=\n",g.inv()
    print "\n(Christoffel Symbols 1st) \n\nGamma1_ijk = ",
    for i in range(n**3):
        L = Base(i,n,3)
        c = L[0]
        a = L[1]
        b = L[2]
        print str(c)+str(a)+str(b)+":",Christoffel1(g,c,a,b),",",

    print "\n\n(Christoffel Symbols 2nd) Gamma2_Ijk = ",
    for i in range(n**3):
        L = Base(i,n,3)
        c = L[0]
        a = L[1]
        b = L[2]
        print str(c)+str(a)+str(b)+":",Christoffel2(g,c,a,b),",",
    print "\n\n(Riemann Curvature Tensor) R_Ijkl = ",
    for i in range(n**4):
        L = Base(i,n,4)
        rho = L[0]
        sig = L[1]
        mu  = L[2]
        nu  = L[3]
        print str(rho)+str(sig)+str(mu)+str(nu)+":",RiemannCurvature(g,rho,sig,mu,nu),",",
    print "\n\n(Ricci Tensor) R_ij = ",
    for i in range(n**2):
        L = Base(i,n,2)
        a = L[0]
        b = L[1]
        print str(a)+str(b)+":",RicciTensor(g,a,b),",",
    print "\n\n(Scalar Curvature) R = ", ScalarCurvature(g)
    return

def ShowPlanet(M, reset=True):
    global x, christoffel1, christoffel2, riemanncurvature, riccitensor
    if reset:
        christoffel1 = {}
        christoffel2 = {}
        riemanncurvature = {}
        riccitensor = {}
    x = ["t","r","theta","phi"]
    t,r,theta,phi = map(Symbol,x)
    G = Symbol('G') #6.673e-11 # gravitational constant
    c = Symbol('c') #299792458.0
    g = Matrix([
        [(1-(2*G*M)/(r*c**2)), 0, 0, 0],
         [0, -(r*c**2)/(r*c**2-2*G*M), 0, 0],
          [0,0,-r**2, 0],
          [0,0,0,-r**2*(sin(theta))**2]
          ])
    n = 4
    print "\n(Inverse Metric Tensor) g_IJ=\n",g.inv()
    print "\n(Christoffel Symbols 1st) \n\nGamma1_ijk = ",
    for i in range(n**3):
        L = Base(i,n,3)
        c = L[0]
        a = L[1]
        b = L[2]
        print str(c)+str(a)+str(b)+":",Christoffel1(g,c,a,b),",",

    print "\n\n(Christoffel Symbols 2nd) Gamma2_Ijk = ",
    for i in range(n**3):
        L = Base(i,n,3)
        c = L[0]
        a = L[1]
        b = L[2]
        print str(c)+str(a)+str(b)+":",Christoffel2(g,c,a,b),",",
    print "\n\n(Riemann Curvature Tensor) R_Ijkl = ",
    for i in range(n**4):
        L = Base(i,n,4)
        rho = L[0]
        sig = L[1]
        mu  = L[2]
        nu  = L[3]
        print str(rho)+str(sig)+str(mu)+str(nu)+":",RiemannCurvature(g,rho,sig,mu,nu),",",
    print "\n\n(Ricci Tensor) R_ij = ",
    for i in range(n**2):
        L = Base(i,n,2)
        a = L[0]
        b = L[1]
        print str(a)+str(b)+":",RicciTensor(g,a,b),",",
    print "\n\nCalculating with Schwartzschild metric g = \n", g
    R = ScalarCurvature(g)
    print "\n\nScalar Curvature R = ", R 
    return R


#ShowPolar()
M_earth = Symbol('M') #5.972e24 # mass of earth
R = ShowPlanet(M_earth, True)

