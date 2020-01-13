

def rshift(b):
   br = (b>>1) + 4*(b&1)
   return br

def lshift(b):
   return rshift(rshift(b))

def parent(i):
   return (i/2)

def corner(i,j,k):
   return (4*(k%2) + 2*(j%2) + (i%2))


def hilbert3d(i,j,k,level,high):
   length = 1
   for c in range(high-level):
      length *= 8

   if (level==0):
      s = 0
      e = 4
      dist = 0
   else:
      (sp,ep,dist) = hilbert3d(parent(i),parent(j),parent(k),level-1,high)
      crnr = corner(i,j,k)
      pp = [0]*8
      esp = ep^sp
      pp[0] = sp
      pp[1] = pp[0]^(rshift(esp))
      pp[2] = pp[1]^(lshift(esp))
      pp[3] = pp[2]^(rshift(esp))
      pp[4] = pp[3]^(esp)
      pp[5] = pp[4]^(rshift(esp))
      pp[6] = pp[5]^(lshift(esp))
      pp[7] = pp[6]^(rshift(esp))
      if   (crnr==pp[0]):
         dist += 0*length; s = pp[0]; e = pp[1]
      elif (crnr==pp[1]):
         dist += 1*length; s = pp[0]; e = pp[3]
      elif (crnr==pp[2]):
         dist += 2*length; s = pp[0]; e = pp[3]
      elif (crnr==pp[3]):
         dist += 3*length; s = pp[2]; e = pp[5]
      elif (crnr==pp[4]):
         dist += 4*length; s = pp[2]; e = pp[5]
      elif (crnr==pp[5]):
         dist += 5*length; s = pp[4]; e = pp[7]
      elif (crnr==pp[6]):
         dist += 6*length; s = pp[4]; e = pp[7]
      elif (crnr==pp[7]):
         dist += 7*length; s = pp[6]; e = pp[7]
      else:
         print "CRAP=",crnr,pp
         s = pp[6]; e = pp[7]

   return (s,e,dist)


ppp = []
level = 7
N = 2**level
for k in range(N):
   for j in range(N):
      for i in range(N):
         dist = hilbert3d(i,j,k,level,level)
         ppp.append(dist[2])

dist = hilbert3d(0,0,N-1,level,level)
#print "ppp=",ppp
ppp.sort()
sss = set(ppp)
print "ppp2=",len(sss),N*N*N,dist


