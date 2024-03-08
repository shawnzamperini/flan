# distutils: language = c++

from Impurity cimport Impurity

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.
cdef class PyImpurity:
	cdef Impurity c_imp  # Hold a C++ instance which we're wrapping
	
	def __init__(self, int imp_atom_num, float mass, float x, float y, 
                float z, float charge, int fbirth):
		self.c_imp = Impurity(imp_atom_num, mass, x, y, z, charge, fbirth)
		self.c_imp.vx = 0.0
		self.c_imp.vy = 0.0
		self.c_imp.vz = 0.0
	
	# Attribute access
	@property
	def x(self):
		return self.c_imp.x
	@x.setter
	def x(self, x):
		self.c_imp.x = x
	
	# Attribute access
	@property
	def y(self):
		return self.c_imp.y
	@y.setter
	def y(self, y):
		self.c_imp.y = y
	
	# Attribute access
	@property
	def z(self):
		return self.c_imp.z
	@z.setter
	def z(self, z):
		self.c_imp.z = z
    
    # Attribute access
	@property
	def charge(self):
		return self.c_imp.charge
	@charge.setter
	def charge(self, charge):
		self.c_imp.charge = charge            
    
	# Attribute access
	@property
	def imp_atom_num(self):
		return self.c_imp.imp_atom_num
	@imp_atom_num.setter
	def imp_atom_num(self, imp_atom_num):
		self.c_imp.imp_atom_num = imp_atom_num 
	
	# Attribute access
	@property
	def vx(self):
		return self.c_imp.vx
	@vx.setter
	def vx(self, vx):
		self.c_imp.vx = vx
	
	# Attribute access
	@property
	def vy(self):
		return self.c_imp.vy
	@vy.setter
	def vy(self, vy):
		self.c_imp.vy = vy
	
	# Attribute access
	@property
	def vz(self):
		return self.c_imp.vz
	@vz.setter
	def vz(self, vz):
		self.c_imp.vz = vz
	
	# Attribute access
	@property
	def mass(self):
		return self.c_imp.mass
	@mass.setter
	def mass(self, mass):
		self.c_imp.mass = mass
