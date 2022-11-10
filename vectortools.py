# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 01:00:19 2021

@author: junob
"""
import math as m
import numpy as np

class Tensor:
    def __init__(self, xx, xy, yx, yy):
        self.xx = xx
        self.xy = xy
        self.yx = yx
        self.yy = yy
        
    def __pos__(self):
        return self

    def __neg__(self):
        return Tensor(-self.xx, -self.xy, -self.yx, -self.yy)
        
    def __abs__(self):
        return self.xx*self.yy - self.xy*self.yx
    
    def __add__(self, other_Tensor):
        return Tensor(self.xx + other_Tensor.xx, self.xy + other_Tensor.xy,
                      self.yx + other_Tensor.yx, self.yy + other_Tensor.yy)
    
    def __sub__(self, other_Tensor):
        return Tensor(self.xx - other_Tensor.xx, self.xy - other_Tensor.xy,
                      self.yx - other_Tensor.yx, self.yy - other_Tensor.yy)
    
    def __mul__(self, num):
        return Tensor(self.xx*num, self.xy*num, self.yx*num, self.yy*num)
    
    def __rmul__(self, num):
        return Tensor(self.xx*num, self.xy*num, self.yx*num, self.yy*num)
    
    def __truediv__(self, num):
        return Tensor(self.xx/num, self.xy/num, self.yx/num, self.yy/num)
    
    def __pow__(self, num):
        result = Tensor(1, 0, 0, 1)
        if num < 0:
            pow_tensor = self.inverse()
        else:
            pow_tensor = self
        for i in range(abs(num)):
            result = result.dot(pow_tensor)
        return result
    
    def __eq__(self, other):
        if str(self) == str(other):
            return True
        else:
            return False
        
    def __ne__(self, other):
        if str(self) != str(other):
            return True
        else:
            return False
        
    def __str__(self):
        return 'Tensor(' + str(self.xx) + ', ' + str(self.xy) + ', ' + str(self.yx) + ', ' + str(self.yy) + ')'
    
    def inverse(self):
        if abs(self) == 0:
            return None
        else:
            return Tensor(self.yy, -self.xy, -self.yx, self.xx)/abs(self)
        
    def T(self):
        return Tensor(self.xx, self.yx, self.xy, self.yy)
        
    def dot(self, other_object):
        if isinstance(other_object, Tensor):
            return Tensor(self.xx*other_object.xx + self.xy*other_object.yx, self.xx*other_object.xy + self.xy*other_object.yy, 
                          self.yx*other_object.xx + self.yy*other_object.yx, self.yx*other_object.xy + self.yy*other_object.yy)
            
        elif isinstance(other_object, Vector):
            return Vector(self.xx*other_object.x + self.xy*other_object.y, self.yx*other_object.x + self.yy*other_object.y)
        
    def tuple(self):
        return ((self.xx, self.xy),
                (self.yx, self.yy))
    
    def list(self):
        return [[self.xx, self.xy],
                [self.yx, self.yy]]
    
    def array(self):
        return np.array([[self.xx, self.xy],
                         [self.yx, self.yy]])

class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y
     
    def __pos__(self):
        return self

    def __neg__(self):
        return Vector(-self.x, -self.y)
        
    def __abs__(self):
        return m.sqrt(self.x**2 + self.y**2)
    
    def __add__(self, other_Vector):
        return Vector(self.x + other_Vector.x, self.y + other_Vector.y)
    
    def __sub__(self, other_Vector):
        return Vector(self.x - other_Vector.x, self.y - other_Vector.y)
    
    def __mul__(self, num):
        return Vector(self.x*num, self.y*num)
    
    def __rmul__(self, num):
        return Vector(self.x*num, self.y*num)
    
    def __truediv__(self, num):
        return Vector(self.x/num, self.y/num)
    
    def __eq__(self, other):
        if str(self) == str(other):
            return True
        else:
            return False
        
    def __ne__(self, other):
        if str(self) != str(other):
            return True
        else:
            return False
    
    def __str__(self):
        return 'Vector(' + str(self.x) + ', ' + str(self.y) + ')'
    
    def dot(self, other_Vector):
        return self.x*other_Vector.x + self.y*other_Vector.y
    
    def cross(self, other_Vector):
        return self.x*other_Vector.y - self.y*other_Vector.x
    
    def tuple(self):
        return (self.x, self.y)
    
    def list(self):
        return [self.x, self.y]
    
    def array(self):
        return np.array([self.x, self.y])
    
def SO2(theta):
    return Tensor(m.cos(theta), -m.sin(theta), m.sin(theta), m.cos(theta))

if __name__ == '__main__':
    t1 = Tensor(1, 2,
                3, 4)
    t2 = SO2(m.pi/4)
    t3 = Tensor(1, 2,
                3, 4)
    print(abs(t1))
    print(t1 + t2)
    print(t1 - t2)
    print(t1*2)
    print(2*t1)
    print(t1/2)
    print(t1**3)
    print(t1**(-1))
    print(t1 == t3)
    print(t1 != t3)
    print(t1)
    print(t1.inverse())
    print(t1.T())
    print(t1.dot(t2))
    print(t1.tuple())
    print(t1.list())
    print(t1.array())
    
    v1 = Vector(3, 4)
    v2 = Vector(4, 3)
    print(t1.dot(v1))
    print(abs(v1))
    print(v1 + v2)
    print(v1 - v2)
    print(2*v1)
    print(v1*2)
    print(v1/2)
    print(v1 == v2)
    print(v1 != v2)
    print(v1)
    print(v1.dot(v2))
    print(v1.cross(v2))
    print(v1.tuple())
    print(v1.list())
    print(v1.array())
