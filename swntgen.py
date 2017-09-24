#!/usr/bin/env python
#********************************************************************
# Description: Generate a singel wall carbon nanotube from
#              helical and rotational symmetries as described in 
#              White, Robertson and Mintmire, 
#              PRB 47, 5485-5488 (1993).
# *******************************************************************


import math
import copy 
#import subprocess
import argparse

class Atom:

  def __init__(self,x,y,z,atom):
    self.x = x
    self.y = y
    self.z = z
    self.atom = atom

  def move(self,dx,dy,dz):
    self.x = self.x + dx
    self.y = self.y + dy
    self.z = self.z + dz

  def rotate(self,angle):

    new_x = self.x * math.cos(angle) - self.y * math.sin(angle)
    new_y = self.x * math.sin(angle) + self.y * math.cos(angle) 
    new_z = self.z

    self.x = new_x
    self.y = new_y
    self.z = new_z

  def rotate_move(self,angle,dx,dy,dz):
    self.rotate(angle)
    self.move(dx,dy,dz)

  def to_string(self):
    return "%s %3.15f %3.15f %3.15f" % (self.atom,self.x,self.y,self.z)


class Motif:
  def __init__(self,atomlist):
    self.atomlist = atomlist

  def rotate(self,angle):
    for atom in self.atomlist:
      atom.rotate(angle)
 
  def rotate_move(self,angle,dx,dy,dz):
    for atom in self.atomlist:
      atom.rotate_move(angle,dx,dy,dz)

  def to_string(self):
    result = []
    for atom in self.atomlist:
      result.append(atom.to_string())
    return result

  def __add__(self,next):
    return Motif(self.atomlist + next.atomlist)

  def size(self):
    return len(self.atomlist)


# Greatest common divisor
def GCD(n,m):
  if m != 0:
    return GCD(m,n % n)
  return n 

#Determine coefficients for helical vector (p1,p2)
def p_factors(n,m):
  sol = GCD(n,m)
  for p1 in range (0,n+2):
    for p2 in range (0,m+2):
      if p2*n - p1*m == sol:
       return (p1,p2)
  return (0,0)
  
# Default parameters
cc_length = 1.44
n = 5
m = 5
layers = 3

#command line inputs
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--chirality", help="Set chirality vector(n,m)", metavar="\"n m\"")
parser.add_argument("-l", "--layers", type = int, help="Number of layers", metavar="number") 
args = parser.parse_args()

if args.chirality:
 data = args.chirality.split()
 n = int(data[0].strip())
 m = int(data[1].strip())

if args.layers:
 layers = args.layers

# Cell parameter
a = cc_length * math.sqrt(3.0)

#Primitive lattice vectors
R1_x = a
R1_y = 0
R1_z = 0

R2_x = a*math.cos(math.pi/3.0)
R2_y = a*math.sin(math.pi/3.0)
R2_z = 0
   
#Atom positions
d_x = (R1_x + R2_x)/3.0
d_y = (R1_y + R2_y)/3.0
d_z = (R1_z + R2_z)/3.0


#Chiral vector
R_x = n*R1_x + m*R2_x
R_y = n*R1_y + m*R2_y
R_z = 0.0

R_squared = R_x*R_x + R_y*R_y
R_length = math.sqrt(R_squared)

#Rotation and translation of the second atom in the unit cell
rotation = 2.0*math.pi*(d_x*R_x+d_y*R_y)/R_squared
#print (rotation)
a_x = d_y*R_z - d_z*R_y
a_y = d_z*R_x - d_x*R_z
a_z = d_y*R_x - d_x*R_y
a_lenght = math.sqrt(a_x*a_x + a_y*a_y + a_z*a_z)
translation = a_lenght/R_length
#print (translation)

#Makeing of primitve cell motif
tube_radius = R_length/(2.0*math.pi)
motif = []
first_atom = Atom(tube_radius,0,0,"C")
second_atom = Atom(tube_radius,0,0,"C")
second_atom.rotate_move(rotation,0,0,translation)
motif.append(first_atom)
motif.append(second_atom)
basic_motif = Motif(motif)

#Make the helical motif
N = GCD(n,m)
helical_motif = copy.deepcopy(basic_motif)
for i in range(1,N):
    new_motif = copy.deepcopy(basic_motif)
    new_motif.rotate(i*2*math.pi/N)
    helical_motif = helical_motif + new_motif

#Translation and rotation of the  helical motif
(p1,p2)  = p_factors(n,m)
H_x = p1*R1_x + p2*R2_x 
H_y = p1*R1_y + p2*R2_y
H_z = p1*R1_z + p2*R2_z

b_x = H_y*R_z - H_z*R_y
b_y = H_z*R_x - H_x*R_z
b_z = H_y*R_x - H_x*R_y
b_lenght = math.sqrt(b_x*b_x+b_y*b_y+b_z*b_z)

h = b_lenght/R_length
a = 2*math.pi*(H_x*R_x+H_y*R_y+H_z*R_z)/R_squared

#Construction of the whole tube
tube = copy.deepcopy(helical_motif)
i = 1 
nc = 2*(n+m)*layers      # amount of atoms in the tube 
while tube.size() < nc:
   new_motif = copy.deepcopy(helical_motif)
   new_motif.rotate_move(i*a,0,0,i*h)
   tube = tube + new_motif
   i = i + 1
   
   
#Output tube in xyz coordinates
print (tube.size())
print ("(%d %d)" % (n, m))
for line in tube.to_string():
  print (line)

