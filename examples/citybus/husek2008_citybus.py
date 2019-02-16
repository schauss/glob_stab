import sympy as sy

# parameters for citybus from Husek 2008
m = [9950,16000]
v = [1,20]

s, z = sy.symbols('s z')
m, v = sy.symbols('m v', real=True, positive=True)
z=1

# plant (from steering angle to measured displacement from the guiding wire)
n = 6.076e5*m*v**2*s**2 + 3.886e11*v*s + 4.803e10*v**2
d = s**2*(m**2*v**2*s**2+1.075e6*m*v*s+1.663e4*m*v**2+2.690e11)

# plant transfer function for actuator transfer function 1/s
G_p = n/(d*s)
G_p1 = sy.simplify(G_p)
print('G_p:')
print(sy.pretty(G_p))
print('G_p:')
print(sy.pretty(G_p1))

# controller
cn = 2.348e3*s**2+1.094e4*s+9.375e3
cd = s**3+50*s**2+1.25e3*s+1.563e4

# controller transfer function
G_c = cn/cd
print('G_c:')
print(sy.pretty(G_c))

# overall transfer function
G = G_p/(1+G_p*G_c*z)
print('G:')
print(sy.pretty(G))

G1 = sy.simplify(G)
print('G:')
print(sy.pretty(G1))

Gn, Gd = sy.fraction(G1)

p = sy.Poly(Gd, s)
print('G den:')
print(p.all_coeffs())
#for idx, a in p.all_coeffs():
#    print('a[{}][0] = {};'.format(idx,a))

